"""Reproducible 7L39 OpenFold3 -> preparation -> AutoDock Vina benchmark."""
import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, required=True)
parser.add_argument('--case-dir', type=Path, help='Prepared case.json, query.json, reference/ and msa/ directory')
parser.add_argument('--single-sequence', action='store_true', help='Use raw query without MSA for a prepared case')
parser.add_argument('--reuse-prediction', action='store_true',
                    help='Use a completed OpenFold3 prediction in this output directory only')
args = parser.parse_args()
out = args.output_dir.resolve()
out.mkdir(parents=True, exist_ok=True)
if (out / 'manifest.json').exists():
    parser.error('Output has a manifest already; choose a new directory to preserve prior results.')
(out / 'logs').mkdir(exist_ok=True)
case_dir = args.case_dir.resolve() if args.case_dir else None
case = json.loads((case_dir / 'case.json').read_text()) if case_dir else None
if args.single_sequence and not case:
    parser.error('--single-sequence requires --case-dir')
if case_dir and case_dir != out:
    for name in ['case.json', 'query.json']:
        shutil.copy2(case_dir / name, out / name)
    shutil.copytree(case_dir / 'reference', out / 'reference', dirs_exist_ok=True)
msa_dir = case_dir / 'msa' if case else ROOT / 'results/openfold3-7l39/msa'
query_path = case_dir / 'query.json' if args.single_sequence else msa_dir / 'query_msa.json'
manifest = {'started_utc': datetime.now(timezone.utc).isoformat(),
            'platform': platform.platform(), 'output_dir': str(out), 'stages': [],
            'msa': 'Cached remote ColabFold MSA; inference and docking local',
            'prediction_reused_in_this_invocation': args.reuse_prediction,
            'checkpoint': 'of3-ob-2025-06-30-174k.pt',
            'protocol': 'OpenFold3-guided local docking, rigid receptor, toluene, 20 A cube',
            'seeds': [42, 2026, 20261021], 'cpu': 6, 'exhaustiveness': 32,
            'num_modes': 9, 'energy_range': 4, 'status': 'running'}
if case:
    manifest['protocol'] = case['system'] + '; OpenFold3-guided local docking; rigid receptor'
    manifest['case'] = case
if args.single_sequence:
    manifest['msa'] = 'None; single-sequence inference'
started = time.perf_counter()

def save():
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

def run(name, command):
    command = [str(x) for x in command]
    print(f'[{name}] starting', flush=True)
    t = time.perf_counter()
    log = out / 'logs' / f'{name}.log'
    with log.open('w') as stream:
        result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
    stage = {'name': name, 'command': command, 'elapsed_seconds': time.perf_counter() - t,
             'returncode': result.returncode, 'log': str(log)}
    manifest['stages'].append(stage)
    save()
    if result.returncode:
        raise RuntimeError(f'{name} failed; see {log}\n{log.read_text()[-3500:]}')
    print(f'[{name}] passed in {stage["elapsed_seconds"]:.2f} s', flush=True)

of = [ROOT / 'openfold3/run.sh', 'python']
env = [ROOT / 'run.sh']
prep = [ROOT / '.tools/prep-env/bin/python']
try:
    hashes = {}
    inputs = [ROOT / 'openfold3/t4-lysozyme-toluene.json',
              ROOT / 'openfold3/runner.yml',
              ROOT / 'results/openfold3-7l39/msa/query_msa.json',
              ROOT / 'results/openfold3-7l39/reference/7L39.cif',
              ROOT / 'results/openfold3-7l39/reference/MBN.cif',
              ROOT / '.tools/openfold3-cache/of3-ob-2025-06-30-174k.pt']
    if case:
        inputs = [case_dir / 'case.json', case_dir / 'query.json', query_path,
                  ROOT / 'openfold3/runner.yml',
                  ROOT / '.tools/openfold3-cache/of3-ob-2025-06-30-174k.pt']
        inputs += sorted((case_dir / 'reference').glob('*.cif'))
    if not args.single_sequence:
        inputs += sorted((msa_dir / 'main').glob('*.npz'))
    inputs += sorted((ROOT / 'workflow').glob('*.py'))
    inputs += [ROOT / 'openfold3/validate_7l39.py', ROOT / 'openfold3/benchmark.py']
    for path in inputs:
        with path.open('rb') as handle:
            hashes[str(path.relative_to(ROOT))] = hashlib.file_digest(handle, 'sha256').hexdigest()
    manifest['sha256'] = hashes
    save()
    if not args.reuse_prediction:
        run('openfold', of + ['-u', ROOT / 'openfold3/benchmark.py', out / 'openfold/benchmark.json',
            'predict', '--query-json', query_path,
            '--runner-yaml', ROOT / 'openfold3/runner.yml', '--num-diffusion-samples', '1',
            '--use-msa-server', 'false', '--use-templates', 'false',
            '--output-dir', out / 'openfold/prediction'])
    if case:
        run('structure_validation', of + [ROOT / 'workflow/validate_case.py', out])
    else:
        run('structure_validation', of + [ROOT / 'openfold3/validate_7l39.py',
            '--output-dir', out, '--prediction-dir', out / 'openfold'])
    run('independent_ligand', env + ['python', ROOT / 'workflow/make_ligand.py', out])
    run('ligand_pdbqt', env + ['mk_prepare_ligand.py', '-i', out / 'ligand-input.sdf',
                            '-o', out / 'ligand.pdbqt'])
    align = json.loads((out / 'alignment.json').read_text())
    for kind in ['predicted', 'reference']:
        folder = out / kind
        folder.mkdir(exist_ok=True)
        run(f'{kind}_protonation', prep + [ROOT / 'workflow/protonate.py',
            out / f'{kind}-receptor.pdb', folder / 'receptor-h.pdb', folder / 'preparation.json'])
        run(f'{kind}_receptor_pdbqt', env + ['mk_prepare_receptor.py',
            '--read_pdb', folder / 'receptor-h.pdb', '-o', folder / 'receptor', '-p', '-j', '-v',
            '--box_size'] + align['box_size_angstrom'] + ['--box_center'] + align[f'{kind}_box_center'])
        for seed in manifest['seeds']:
            base = folder / f'seed-{seed}'
            run(f'{kind}_vina_{seed}', env + ['vina', '--receptor', folder / 'receptor.pdbqt',
                '--ligand', out / 'ligand.pdbqt', '--config', folder / 'receptor.box.txt',
                '--cpu', '6', '--exhaustiveness', '32', '--seed', str(seed),
                '--num_modes', '9', '--energy_range', '4', '--out', base.with_suffix('.pdbqt')])
            run(f'{kind}_export_{seed}', env + ['mk_export.py', base.with_suffix('.pdbqt'),
                                            '-s', base.with_suffix('.sdf')])
    run('evaluate', env + ['python', ROOT / 'workflow/evaluate.py', out])
    manifest['status'] = 'completed'
except BaseException as exc:
    manifest['status'] = 'failed'
    manifest['error'] = str(exc)
    raise
finally:
    manifest['orchestrator_seconds'] = time.perf_counter() - started
    manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
    save()
