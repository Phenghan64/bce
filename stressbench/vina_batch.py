"""Repeated balanced docking batches with per-job QC, resumability and fault injection."""
import argparse,concurrent.futures,datetime,hashlib,json,math,os,re,subprocess,threading,time
from pathlib import Path
from rdkit import Chem
from meeko import PDBQTMolecule,RDKitMolCreate
ROOT=Path(__file__).resolve().parents[1];OUT=Path(os.environ.get('STRESSBENCH_OUTPUT',ROOT/'results/stress-20261022-main'))
p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--workers',type=int,required=True);p.add_argument('--cpu',type=int,required=True);p.add_argument('--exhaustiveness',type=int,required=True);p.add_argument('--duration',type=float,required=True);p.add_argument('--panel',choices=['calibration','full'],default='calibration');p.add_argument('--job-timeout',type=float,default=240);p.add_argument('--cycles',type=int,default=100000);p.add_argument('--inject-after',type=float,default=-1);args=p.parse_args()
folder=Path(args.folder);folder.mkdir(parents=True,exist_ok=True)
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic(path,obj):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2)+'\n');temp.replace(path)
library=json.loads((OUT/'inputs/library.json').read_text());records={r['id']:r for r in library['records']};ids=library['calibration_panel'] if args.panel=='calibration' else list(records)
running={};mutex=threading.Lock();injected=threading.Event();done=threading.Event();start=time.monotonic();results=[];peak_workers=0
base={'workers':args.workers,'cpu':args.cpu,'exhaustiveness':args.exhaustiveness,'seed':20261022,'receptor_sha256':sha(OUT/'inputs/receptor.pdbqt'),'box_sha256':sha(OUT/'inputs/receptor.box.txt'),'vina_sha256':sha(ROOT/'.tools/bin/vina'),'num_modes':9,'energy_range':4}
def one(cycle,ligand):
    global peak_workers
    rec=records[ligand];jobfolder=folder/'jobs'/f'{cycle:05d}-{ligand}';jobfolder.mkdir(parents=True,exist_ok=True)
    spec={**base,'ligand_id':ligand,'ligand_sha256':rec['pdbqt_sha256'],'cycle':cycle};signature=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest();path=jobfolder/'result.json'
    if path.exists():
        old=json.loads(path.read_text());assert old['signature']==signature
        if old['status']=='completed' and all(Path(k).exists() and sha(k)==v for k,v in old['hashes'].items()):return {**old,'cache_hit':True}
    previous_attempts=len(list(jobfolder.glob('attempt-*.log')));attempts=[]
    for attempt in range(previous_attempts+1,previous_attempts+3):
        dest=jobfolder/f'attempt-{attempt}.pdbqt';sdf_path=jobfolder/f'attempt-{attempt}.sdf';logpath=jobfolder/f'attempt-{attempt}.log'
        command=[str(ROOT/'.tools/bin/vina'),'--receptor',str(OUT/'inputs/receptor.pdbqt'),'--ligand',rec['pdbqt'],'--config',str(OUT/'inputs/receptor.box.txt'),'--cpu',str(args.cpu),'--exhaustiveness',str(args.exhaustiveness),'--seed','20261022','--num_modes','9','--energy_range','4','--out',str(dest)]
        launched=time.monotonic();r={'signature':signature,'spec':spec,'status':'running','cache_hit':False,'started_utc':utc(),'start_offset_seconds':launched-start,'torsdof':rec['torsdof'],'complexity_group':rec['complexity_group'],'attempt':attempt,'command':command}
        atomic(path,r)
        try:
            with logpath.open('w') as log:
                process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
                with mutex:
                    running[process.pid]={'process':process,'job':str(path),'start':launched}
                    peak_workers=max(peak_workers,sum(v['process'].poll() is None for v in running.values()))
                try:code=process.wait(timeout=args.job_timeout)
                except subprocess.TimeoutExpired:process.kill();process.wait();raise TimeoutError('per_job_wall_time_limit')
                finally:
                    with mutex:running.pop(process.pid,None)
            r['docking_wall_seconds']=time.monotonic()-launched;r['returncode']=code
            if code:raise RuntimeError(f'vina_exit_{code}')
            text=dest.read_text();energies=[float(v) for v in re.findall(r'^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)',text,re.M)]
            assert energies and all(math.isfinite(v) for v in energies) and energies==sorted(energies)
            sdf,failures=RDKitMolCreate.write_sd_string(PDBQTMolecule(text,skip_typing=True));assert not failures;sdf_path.write_text(sdf)
            mols=list(Chem.SDMolSupplier(str(sdf_path),removeHs=True));assert len(mols)==len(energies)
            expected=Chem.MolFromSmiles(rec['smiles']);identity=Chem.MolToSmiles(expected,isomericSmiles=False)
            for mol in mols:
                assert mol is not None and Chem.MolToSmiles(mol,isomericSmiles=False)==identity,'molecular_connectivity_changed'
                assert all(math.isfinite(float(v)) for xyz in mol.GetConformer().GetPositions() for v in xyz)
            r.update(status='completed',scores=energies,pose_count=len(mols),hashes={str(dest):sha(dest),str(sdf_path):sha(sdf_path)})
        except Exception as error:r.update(status='failed',error=f'{type(error).__name__}: {error}')
        r.update(job_wall_seconds=time.monotonic()-launched,finish_offset_seconds=time.monotonic()-start,finished_utc=utc());attempts.append(r.copy());r['attempts']=attempts;atomic(path,r)
        if r['status']=='completed':return r
        # Only the deliberately terminated job is automatically retried.
        fault=folder/'fault-injection.json'
        if not fault.exists() or json.loads(fault.read_text()).get('job')!=str(path):return r
    return r
def injector():
    if args.inject_after<0:return
    if done.wait(args.inject_after):return
    deadline=time.monotonic()+120
    while time.monotonic()<deadline and not done.is_set():
        with mutex:
            targets=list(running.values())
            if targets:
                target=min(targets,key=lambda r:r['start']);pid=target['process'].pid
                atomic(folder/'fault-injection.json',{'pid':pid,'job':target['job'],'utc':utc(),'signal':'SIGTERM','owned_worker_only':True});target['process'].terminate();injected.set();return
        done.wait(.1)
threading.Thread(target=injector,daemon=True).start()
try:
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for cycle in range(args.cycles):
            # An explicit finite plan must contain exactly the same jobs on
            # fresh and cached runs. A duration cutoff changes that plan when
            # cached cycles complete faster; only timed workloads use it.
            if cycle and args.cycles==100000 and time.monotonic()-start>=args.duration:break
            cycle_started=time.monotonic()
            futures=[pool.submit(one,cycle,ligand) for ligand in ids]
            for future in concurrent.futures.as_completed(futures):
                r=future.result();results.append(r)
                atomic(folder/'progress.json',{'utc':utc(),'wall_seconds':time.monotonic()-start,'cycle':cycle,'completed':sum(r['status']=='completed' for r in results),'failed':sum(r['status']!='completed' for r in results),'cache_hits':sum(r['cache_hit'] for r in results),'workers':args.workers,'cpu':args.cpu,'exhaustiveness':args.exhaustiveness})
            with (folder/'cycles.jsonl').open('a') as f:f.write(json.dumps({'cycle':cycle,'start_seconds':cycle_started-start,'finish_seconds':time.monotonic()-start,'jobs':len(ids),'wall_seconds':time.monotonic()-cycle_started})+'\n')
            # A timed-out ligand must not turn the scheduled 90-minute soak into
            # a short run. Preserve failures and keep the full-panel soak loaded.
            # Stop if an entire cycle fails, which indicates a systemic problem.
            failed_this_cycle=sum(r['status']!='completed' for r in results[-len(ids):])
            if failed_this_cycle and (args.panel!='full' or failed_this_cycle==len(ids)):break
finally:done.set()
wall=time.monotonic()-start;success=sum(r['status']=='completed' for r in results);cache=sum(r['cache_hit'] for r in results)
summary={'parameters':vars(args),'base':base,'wall_seconds':wall,'completed':success,'failed':len(results)-success,'cache_hits':cache,'fresh_jobs_per_hour':(success-cache)*3600/wall,'injected':injected.is_set(),'peak_actual_active_workers':peak_workers,'results':results}
atomic(folder/('resume-summary.json' if cache else 'summary.json'),summary)
print(json.dumps({k:v for k,v in summary.items() if k!='results'}),flush=True)
if summary['failed']:raise SystemExit(1)
