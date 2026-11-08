"""Resumable real-case docking; job engine adapted from screening/run.py."""
import argparse
import sys
import concurrent.futures
import hashlib
import json
import math
import os
import re
import signal
import subprocess
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from realcase.audit_outputs import compare_identity
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, required=True)
parser.add_argument('--stage', choices=['all','calibrate','screen','refine'], default='all')
parser.add_argument('--retry-failed', action='store_true', help='Explicitly rerun failed jobs; default resume preserves recorded failures')
args = parser.parse_args()
out = args.output_dir.resolve()
library = json.loads((out/'library.json').read_text())
records = {r['id']:r for r in library['records'] if r['status']=='prepared'}
control = 'crystal-control'
receptor = out/'receptor/receptor.pdbqt'
box = out/'receptor/receptor.box.txt'
vina = ROOT/'.tools/bin/vina'
version = subprocess.check_output([str(vina),'--version'],text=True).strip()

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def atomic(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(path)

base = {'vina':version, 'receptor':str(receptor), 'receptor_sha256':digest(receptor),
        'protocol_sha256':digest(out/'protocol.json'), 'library_sha256':digest(out/'library.json'), 'box':str(box), 'box_sha256':digest(box), 'num_modes':9, 'energy_range':4,
        'scoring':'vina', 'timeout_seconds':900}
(out/'jobs').mkdir(exist_ok=True)
file_lock = out/'runner.lock'
try:
    lock_handle = os.open(file_lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL)
except FileExistsError:
    raise SystemExit('Runner lock exists. Verify its process has exited before removing a stale lock.')
os.write(lock_handle,str(os.getpid()).encode()); os.close(lock_handle)

def job(stage, ligand_id, seed, cpu, exhaustiveness):
    entry = records[ligand_id]
    folder = out/'jobs'/stage/f'{ligand_id}-s{seed}'
    folder.mkdir(parents=True,exist_ok=True)
    spec = {**base, 'ligand_id':ligand_id,'ligand_sha256':digest(entry['pdbqt']),
            'seed':seed,'cpu':cpu,'exhaustiveness':exhaustiveness}
    signature = hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    result_path = folder/'result.json'
    if result_path.exists():
        old = json.loads(result_path.read_text())
        if old['signature'] != signature:
            raise RuntimeError(f'Job inputs changed: {folder}. Choose a new output stage.')
        if old['status']=='failed' and not args.retry_failed:
            return {**old,'cache_hit':True}
        if old['status']=='completed' and all(Path(p).exists() and digest(p)==h for p,h in old['outputs_sha256'].items()):
            return {**old,'cache_hit':True}
    started = time.perf_counter()
    dest = folder/'docked.pdbqt'
    command = ['/usr/bin/time','-l',str(vina),'--receptor',str(receptor),
        '--ligand',entry['pdbqt'],'--config',str(box),'--cpu',str(cpu),
        '--exhaustiveness',str(exhaustiveness),'--seed',str(seed),
        '--num_modes','9','--energy_range','4','--out',str(dest)]
    result = {'spec':spec,'signature':signature,'command':command,'cache_hit':False,
              'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running'}
    atomic(result_path,result)
    try:
        with (folder/'vina.log').open('w') as log:
            process = subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                code = process.wait(timeout=base['timeout_seconds'])
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGTERM)
                try: process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL); process.wait()
                raise TimeoutError('Vina exceeded per-ligand time limit')
        result['vina_wall_seconds'] = time.perf_counter()-started
        result['returncode'] = code
        if code:
            raise RuntimeError(f'Vina returned {code}')
        text = dest.read_text()
        energies = [float(v) for v in re.findall(r'^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)',text,re.M)]
        if not energies or not all(math.isfinite(x) for x in energies) or energies!=sorted(energies):
            raise ValueError('Missing, nonfinite or unsorted Vina energies')
        pdbqt_molecule = PDBQTMolecule(text,skip_typing=True)
        sdf, failures = RDKitMolCreate.write_sd_string(pdbqt_molecule)
        if failures:
            raise ValueError(f'Meeko export failed for {failures}')
        sdf_path = folder/'docked.sdf'; sdf_path.write_text(sdf)
        molecules = list(Chem.SDMolSupplier(str(sdf_path),removeHs=True))
        if len(molecules)!=len(energies) or len(energies)!=len(re.findall(r'^MODEL\s+',text,re.M)):
            raise ValueError('PDBQT / SDF / score count mismatch')
        expected = Chem.MolFromSmiles(entry['smiles'])
        normalized_sites = []
        # Check every pose, including molecular identity and finite coordinates.
        for mol in molecules:
            if mol is None:
                raise ValueError('Failed SDF molecule parse')
            same, cleared = compare_identity(mol,expected)
            if cleared: normalized_sites.append(cleared)
            if not same:
                raise ValueError('Exported ligand identity changed')
            if not all(math.isfinite(float(v)) for xyz in mol.GetConformer().GetPositions() for v in xyz):
                raise ValueError('Nonfinite coordinate')
        usage = (folder/'vina.log').read_text()
        memory = re.search(r'(\d+)\s+maximum resident set size',usage)
        result.update(status='completed',top_score=energies[0],scores=energies,
            pose_count=len(energies),pdbqt=str(dest),sdf=str(sdf_path),NH_stereo_normalization=normalized_sites,
            peak_process_rss_bytes=int(memory.group(1)) if memory else None,
            outputs_sha256={str(dest):digest(dest),str(sdf_path):digest(sdf_path)})
    except Exception as exc:
        result.update(status='failed',error=f'{type(exc).__name__}: {exc}')
    result['finished_utc'] = datetime.now(timezone.utc).isoformat()
    result['job_wall_seconds'] = time.perf_counter()-started
    atomic(result_path,result)
    return result

def stage(name, ids, seeds, workers, cpu, exhaustiveness):
    params = {'name':name,'ligand_ids':ids,'seeds':seeds,'workers':workers,'cpu_per_job':cpu,
              'exhaustiveness':exhaustiveness,**base}
    path = out/f'{name}.json'
    if path.exists() and json.loads(path.read_text())['parameters'] != params:
        raise RuntimeError(f'Stage settings changed: {name}. Refusing to overwrite prior experiment.')
    start = time.perf_counter(); results = []
    print(f'{name}: {len(ids)*len(seeds)} jobs, {workers} workers × {cpu} threads',flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(job,name,i,s,cpu,exhaustiveness) for i in ids for s in seeds]
        for future in concurrent.futures.as_completed(futures):
            result = future.result(); result['stage_elapsed_at_completion'] = time.perf_counter()-start; results.append(result)
            if len(results)%5==0 or result['status']=='failed' or len(results)==len(futures):
                progress={'stage':name,'finished':len(results),'total':len(futures),'qc_pass':sum(r['status']=='completed' for r in results),'failed':sum(r['status']=='failed' for r in results),'elapsed_seconds':time.perf_counter()-start}
                atomic(out/'progress.json',progress);print(json.dumps(progress),flush=True)
    wall = time.perf_counter()-start
    summary = {'parameters':params,'wall_seconds':wall,
        'completed':sum(r['status']=='completed' for r in results),
        'failed':sum(r['status']=='failed' for r in results),
        'cache_hits':sum(r['cache_hit'] for r in results),
        'results':sorted(results,key=lambda r:(r['spec']['ligand_id'],r['spec']['seed']))}
    if summary['cache_hits'] == len(results) and path.exists():
        # Keep measured original timing; record the resume check separately.
        atomic(out/f'{name}-resume-check.json',summary)
    else:
        atomic(path,summary)
    return json.loads(path.read_text()) if path.exists() else summary

try:
    from rdkit.Chem import rdMolAlign
    run=stage('calibrate',[control],[42,2026,20261021],3,2,32)
    ref=Chem.SDMolSupplier(str(out/'receptor/crystal-control.sdf'),removeHs=True)[0]
    checks=[]
    for r in run['results']:
        check={'seed':r['spec']['seed'],'status':r['status']}
        if r['status']=='completed':
            mol=Chem.SDMolSupplier(r['sdf'],removeHs=True)[0]
            check.update(top_rmsd_A=float(rdMolAlign.CalcRMS(mol,ref)),top_score=r['top_score'])
        checks.append(check)
    passed=all(c['status']=='completed' and c['top_rmsd_A']<=2 for c in checks)
    atomic(out/'calibration.json',{'passed':passed,'checks':checks,'note':'Performance-only panel continues regardless of pose calibration. No biological ranking validation inferred.'})
    ids=[i for i in records if i!=control];assert len(ids)==10
    result=stage('screen',ids,[42,2026,20261021],3,2,32)
    if result['failed']:raise RuntimeError('Panel contains failed jobs; see results')
finally:
    file_lock.unlink(missing_ok=True)
