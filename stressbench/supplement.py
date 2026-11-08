"""Sequential follow-ups after the main campaign; preserve all original failures."""
import json, os, subprocess, time
from pathlib import Path
import numpy as np
import psutil
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from common import ROOT, OUT, atomic, utc, sha

PY = str(ROOT / '.tools/enzyme-case/venv/bin/python')
VPY = str(ROOT / '.env/bin/python')

def read(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default

def validate(folder, query, count=1):
    q = read(query)
    expected = sum(len(c['sequence']) * len(c['chain_ids'])
                   for item in q['queries'].values() for c in item['chains'])
    files = sorted((folder / 'prediction').rglob('*_model.cif'))
    assert len(files) == count, (len(files), count)
    checks = []
    for file in files:
        cif = MMCIF2Dict(str(file))
        xyz = np.array([cif[f'_atom_site.Cartn_{c}'] for c in 'xyz'], dtype=float)
        assert np.isfinite(xyz).all()
        ca = {(ch, res) for atom, ch, res in zip(cif['_atom_site.label_atom_id'],
              cif['_atom_site.label_asym_id'], cif['_atom_site.label_seq_id']) if atom == 'CA'}
        assert len(ca) == expected, (len(ca), expected)
        checks.append({'file': str(file), 'sha256': sha(file),
                       'CA_residues': len(ca), 'finite_coordinates': True})
    config = read(folder / 'prediction/model_config.json')
    return {'checks': checks, 'expected_residues': expected,
            'model_recycles': config['architecture']['shared']['num_recycles'],
            'diffusion_steps': config['architecture']['shared']['diffusion']['no_full_rollout_steps'],
            'actual_msa_model_max_rows': config['architecture']['msa']['msa_module_embedder']['max_subsampled_all_msa']}

def main():
    original = read(OUT / 'completion.json')
    assert original is not None, 'Main campaign has not ended'
    assert not (OUT / 'controller.lock').exists(), 'Main controller still owns lock'
    pid = read(OUT / 'controller-start.json')['pid']
    assert not psutil.pid_exists(pid), 'Wait until main controller exits'
    assert not (OUT / 'STOP').exists(), 'Operator stop takes priority'
    os.environ['STRESSBENCH_OUTPUT'] = str(OUT)
    os.environ['STRESSBENCH_OWNER_PID'] = str(os.getpid())
    start = time.monotonic()
    caffeine = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(os.getpid())])
    stages = []

    def stage(name, timeout, command):
        folder = OUT / name
        assert not (OUT / 'STOP').exists(), 'Operator stop requested'
        # Completed evidence is never overwritten by a script restart.
        existing = read(folder / 'resources.json')
        if existing:
            stages.append({'name': name, **existing, 'already_recorded': True})
            return existing
        folder.mkdir(parents=True, exist_ok=True)
        progress = {'status': 'running', 'phase': name, 'controller_pid': os.getpid(),
                    'utc': utc(), 'started_monotonic': start}
        atomic(OUT / 'supplement-progress.json', progress)
        print(json.dumps(progress), flush=True)
        with (folder / 'monitor.log').open('w') as log:
            subprocess.run([PY, str(ROOT / 'stressbench/common.py'), name,
                            str(timeout), *command], stdout=log, stderr=subprocess.STDOUT)
        r = read(folder / 'resources.json')
        assert r is not None, 'No monitor evidence: ' + name
        stages.append({'name': name, **r})
        atomic(OUT / 'supplement-stages.json', stages)
        if r.get('guard_stop') == 'operator_stop_file':
            raise InterruptedError('Operator stop')
        return r

    try:
        audit=[]
        reference_scale=read(ROOT/'results/enzyme-case/functional-scale.json')
        for mixed_name in ('mixed-full','mixed-reduced'):
            work=OUT/mixed_name/'work'
            for file in sorted(work.glob('fold-*/benchmark.json')):
                b=read(file)
                if not b.get('completed'):continue
                check=validate(file.parent,Path(b['query']),b['samples_requested'])
                audit.append({'kind':'fold','folder':str(file.parent),'validation':check})
            for file in sorted(work.glob('sequence-*/timing.json')):
                scale=read(file.parent/'functional-scale.json')
                assert scale==reference_scale, ('Sequence output counts changed',str(file.parent))
                assert read(file)['total_compute_seconds']>0
                audit.append({'kind':'sequence','folder':str(file.parent),
                              'functional_counts_match_reference':True,
                              'scale':scale,'scale_sha256':sha(file.parent/'functional-scale.json')})
        atomic(OUT/'mixed-output-validation.json',{'passed':True,'utc':utc(),'checks':audit})
        name = 'cpu-extended/1x12-E64'
        r = stage(name, 1800, [VPY, str(ROOT / 'stressbench/vina_batch.py'),
                  '--folder', str(OUT / name / 'batch'), '--workers', '1', '--cpu', '12',
                  '--exhaustiveness', '64', '--duration', '0', '--panel', 'calibration',
                  '--cycles', '1', '--job-timeout', '900'])
        plan = read(OUT / 'supplement-plan.json')
        atomic(OUT / 'supplement-plan.json', {**plan, 'status': 'executed',
               'finished_utc': utc(), 'resources': r})

        name = 'gpu-extended/n768'
        query = OUT / 'inputs/queries/myosin_prefix_768.json'
        model = OUT / name / 'model'
        r = stage(name, 1800, [str(ROOT / 'openfold3/run.sh'), 'python',
                  str(ROOT / 'stressbench/fold_worker.py'), '--folder', str(model),
                  '--query', str(query), '--msa-depth', '64', '--samples', '1', '--fraction', '1.0'])
        case = {'name': name, 'query': str(query), 'msa_depth': 64, 'samples': 1,
                'fraction': 1.0, 'resources': r, 'success': False,
                'extended_time_budget_seconds': 1800}
        benchmark = read(model / 'benchmark.json')
        if benchmark: case['benchmark'] = benchmark
        if r['returncode'] == 0 and not r['guard_stop']:
            try:
                case['validation'] = validate(model, query)
                case['success'] = True
            except Exception as exc:
                case['validation_error'] = repr(exc)
        if not case['success']:
            log = (OUT / name / 'run.log').read_text(errors='replace')
            case['failure_class'] = ('mps_allocator_oom' if 'MPS backend out of memory' in log
                or 'MPS out of memory' in log else r['guard_stop'] or 'software_or_validation_error')
            case['error_tail'] = log[-2500:]
        atomic(OUT / name / 'validation.json', case)
        plan = read(OUT / 'gpu-supplement-plan.json')
        atomic(OUT / 'gpu-supplement-plan.json', {**plan, 'status': 'executed',
               'finished_utc': utc(), 'success': case['success'], 'resources': r})

        verdict = read(OUT / 'recovery-verdict.json', {})
        if not verdict or verdict.get('resume_summary', {}).get('cache_hits') != 96:
            # Same signatures and original batch directory: only missing or
            # invalidated jobs rerun. The original capped phase stays intact.
            batch = OUT / 'recovery/batch'
            assert (batch / 'fault-injection.json').exists(), 'Original recovery injection missing'
            snapshot=OUT/'recovery-before-continuation.json'
            if not snapshot.exists():
                atomic(snapshot,{'utc':utc(),'progress':read(batch/'progress.json'),
                    'jobs':[read(f) for f in sorted((batch/'jobs').glob('*/result.json'))]})
            cmd = [VPY, str(ROOT / 'stressbench/vina_batch.py'), '--folder', str(batch),
                   '--workers', '3', '--cpu', '2', '--exhaustiveness', '64',
                   '--duration', '600', '--panel', 'calibration', '--cycles', '2',
                   '--job-timeout', '900']
            continuation = stage('recovery-continued', 2400, cmd)
            total = read(batch / 'resume-summary.json', read(batch / 'summary.json'))
            assert continuation['returncode'] == 0 and not continuation['guard_stop'], continuation
            assert total['completed'] == 96 and total['failed'] == 0, total['completed']
            atomic(OUT / 'recovery-continuation-summary.json', total)
            resumed = stage('recovery-continued-resume', 120, cmd)
            cached = read(batch / 'resume-summary.json')
            assert resumed['returncode'] == 0 and cached['cache_hits'] == 96, cached['cache_hits']
            atomic(OUT / 'recovery-verdict.json', {
                'fault_was_injected': True, 'completed': 96, 'failed': 0,
                'resume_cache_hits': 96, 'resume_seconds': cached['wall_seconds'],
                'main_stage_was_truncated': bool(read(OUT / 'recovery/resources.json', {}).get('guard_stop')),
                'continuation_resources': continuation, 'resume_resources': resumed,
                'summary': total, 'resume_summary': cached})
        atomic(OUT / 'supplement-completion.json', {'status': 'completed', 'utc': utc(),
               'wall_seconds': time.monotonic() - start, 'stages': len(stages)})
        atomic(OUT / 'supplement-progress.json', {'status': 'completed', 'phase': 'completed', 'utc': utc()})
    except BaseException as exc:
        atomic(OUT / 'supplement-completion.json', {'status': 'stopped', 'utc': utc(),
               'reason': repr(exc), 'wall_seconds': time.monotonic() - start})
        raise
    finally:
        caffeine.terminate()
        subprocess.run([PY, str(ROOT / 'stressbench/report.py')], cwd=ROOT)

if __name__ == '__main__':
    main()
