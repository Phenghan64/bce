"""Bounded multi-hour capacity campaign. No worker may alter prior experiments."""
import json,math,os,random,signal,subprocess,sys,time
from pathlib import Path
import numpy as np
import psutil
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from common import ROOT,OUT,atomic,utc,sha
PY=str(ROOT/'.tools/enzyme-case/venv/bin/python');VPY=str(ROOT/'.env/bin/python')
start=time.monotonic();stages=[];fold_results=[];active=None;stop_reason=None
lock=OUT/'controller.lock'
fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
os.environ['STRESSBENCH_OUTPUT']=str(OUT);os.environ['STRESSBENCH_OWNER_PID']=str(os.getpid())
caffeine=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
def event(kind,**data):
    rec={'utc':utc(),'elapsed_seconds':time.monotonic()-start,'event':kind,**data}
    with (OUT/'events.jsonl').open('a') as f:f.write(json.dumps(rec)+'\n')
    print(json.dumps(rec),flush=True)
def progress(phase,**data):
    atomic(OUT/'progress.json',{'status':'running','phase':phase,'controller_pid':os.getpid(),'utc':utc(),'elapsed_seconds':time.monotonic()-start,'completed_stages':len(stages),'latest_stages':stages[-4:],**data})
def run_stage(name,command,timeout):
    global active
    if (OUT/'STOP').exists():raise InterruptedError('STOP file detected')
    if time.monotonic()-start>6*3600:raise TimeoutError('Six-hour campaign limit reached')
    folder=OUT/name;folder.mkdir(parents=True,exist_ok=True)
    if (folder/'resources.json').exists():raise RuntimeError('Existing evidence: '+name)
    event('stage_start',name=name,command=command);progress(name)
    with (folder/'monitor.log').open('w') as log:
        active=subprocess.Popen([PY,str(ROOT/'stressbench/common.py'),name,str(timeout),*command],stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
        while active.poll() is None:
            progress(name,stage_folder=str(folder));time.sleep(5)
    resource_file=folder/'resources.json'
    if not resource_file.exists():raise RuntimeError('Monitor failed: '+str(folder/'monitor.log'))
    result=json.loads(resource_file.read_text());result['name']=name;stages.append(result);active=None
    event('stage_end',name=name,returncode=result['returncode'],guard_stop=result['guard_stop'],wall_seconds=result['wall_seconds'])
    atomic(OUT/'stages.json',stages)
    if result['guard_stop']=='operator_stop_file':raise InterruptedError('STOP file detected')
    if result['guard_stop']:
        until=time.monotonic()+60
        while time.monotonic()<until and psutil.virtual_memory().available<2*2**30:time.sleep(2)
    return result
def vina(name,workers,cpu,e,duration,panel='calibration',cycles=None,inject=None,timeout=None):
    command=[VPY,str(ROOT/'stressbench/vina_batch.py'),'--folder',str(OUT/name/'batch'),'--workers',str(workers),'--cpu',str(cpu),'--exhaustiveness',str(e),'--duration',str(duration),'--panel',panel]
    if cycles is not None:command+=['--cycles',str(cycles)]
    if inject is not None:command+=['--inject-after',str(inject)]
    result=run_stage(name,command,timeout or duration+600)
    path=OUT/name/'batch/summary.json'
    summary=json.loads(path.read_text()) if path.exists() else None
    return result,summary,command
def validate_fold(folder,query,samples):
    q=json.loads(Path(query).read_text());expected=sum(len(c['sequence'])*len(c['chain_ids']) for x in q['queries'].values() for c in x['chains'])
    files=sorted((folder/'prediction').rglob('*_model.cif'));checks=[]
    for path in files:
        cif=MMCIF2Dict(str(path));xyz=np.array([cif[f'_atom_site.Cartn_{c}'] for c in 'xyz'],dtype=float)
        assert np.isfinite(xyz).all()
        names=cif['_atom_site.label_atom_id'];resids=cif['_atom_site.label_seq_id'];chains=cif['_atom_site.label_asym_id'];ca={(ch,r) for a,r,ch in zip(names,resids,chains) if a=='CA'}
        assert len(ca)==expected,(path,len(ca),expected)
        checks.append({'file':str(path),'sha256':sha(path),'CA_residues':len(ca),'atoms':len(names),'finite_coordinates':True})
    assert len(files)==samples,(len(files),samples)
    actual=json.loads((folder/'prediction/model_config.json').read_text())
    return {'checks':checks,'expected_residues':expected,'model_recycles':actual['architecture']['shared']['num_recycles'],'diffusion_steps':actual['architecture']['shared']['diffusion']['no_full_rollout_steps'],'actual_msa_model_max_rows':actual['architecture']['msa']['msa_module_embedder']['max_subsampled_all_msa']}
def fold(name,query,msa=64,samples=1,timeout=600,fraction=1.0):
    folder=OUT/name/'model';query=OUT/'inputs/queries'/query
    result=run_stage(name,[str(ROOT/'openfold3/run.sh'),'python',str(ROOT/'stressbench/fold_worker.py'),'--folder',str(folder),'--query',str(query),'--msa-depth',str(msa),'--samples',str(samples),'--fraction',str(fraction)],timeout)
    r={'name':name,'query':str(query),'msa_depth':msa,'samples':samples,'fraction':fraction,'resources':result,'success':False}
    if (folder/'benchmark.json').exists():r['benchmark']=json.loads((folder/'benchmark.json').read_text())
    if result['returncode']==0 and not result['guard_stop']:
        try:r['validation']=validate_fold(folder,query,samples);r['success']=True
        except Exception as exc:r['validation_error']=repr(exc)
    text=(OUT/name/'run.log').read_text(errors='replace')
    if not r['success']:
        r['failure_class']='mps_allocator_oom' if 'MPS backend out of memory' in text or 'MPS out of memory' in text else result['guard_stop'] or 'software_or_validation_error'
        r['error_tail']=text[-2500:]
    atomic(OUT/name/'validation.json',r);fold_results.append(r);atomic(OUT/'gpu-results.json',fold_results)
    event('fold_result',name=name,success=r['success'],failure_class=r.get('failure_class'),wall_seconds=result['wall_seconds'])
    return r
def stop_handler(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
signal.signal(signal.SIGTERM,stop_handler)
try:
    atomic(OUT/'controller-start.json',{'utc':utc(),'pid':os.getpid(),'protocol_sha256':sha(OUT/'protocol.json'),'maximum_hours':6})
    atomic(OUT/'software-manifest.json',[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in (ROOT/'stressbench').glob('*.py')])
    progress('baseline')
    while not (OUT/'baseline/resources.json').exists():time.sleep(2)
    base=json.loads((OUT/'baseline/resources.json').read_text());assert base['returncode']==0;stages.append({**base,'name':'baseline'});event('baseline_completed')
    configs=[(w,c,e) for w,c in [(1,12),(2,6),(3,4),(6,2),(12,1),(16,1)] for e in (8,32,64)]
    random.Random(20261022).shuffle(configs)
    calibration=[]
    for w,c,e in configs:
        result,summary,_=vina(f'cpu-calibration/{w}x{c}-E{e}',w,c,e,120,timeout=600)
        calibration.append({'workers':w,'cpu':c,'e':e,'resources':result,'summary':summary});atomic(OUT/'cpu-calibration.json',calibration)
        if summary and summary['failed']:event('calibration_qc_failure',configuration=[w,c,e])
    eligible=[r for r in calibration if r['e']==32 and r['summary'] and r['summary']['failed']==0 and r['resources']['guard_stop'] is None]
    if not eligible:raise RuntimeError('No QC-passing E32 configuration; refuse unvalidated soak')
    best=max(eligible,key=lambda r:r['summary']['fresh_jobs_per_hour']);workers,cpu=best['workers'],best['cpu']
    atomic(OUT/'selected-cpu-configuration.json',{'workers':workers,'cpu':cpu,'selection':'Highest QC-passing E32 complete-balanced-panel throughput; includes oversubscribed layout if empirically best.','calibration_jobs_per_hour':best['summary']['fresh_jobs_per_hour']})
    event('cpu_configuration_selected',workers=workers,cpu=cpu)
    gpu_start=time.monotonic();gpu_budget=4500;largest=None;failure_size=None;passed_size=0
    for n in (256,512,768,1024,1536,1964):
        if time.monotonic()-gpu_start>gpu_budget-600:break
        result=fold(f'gpu-size/n{n}',f'myosin_prefix_{n}.json')
        if result['success']:largest=result;passed_size=n
        else:failure_size=n;break
    if failure_size and passed_size and failure_size-passed_size>=256:
        midpoint=(passed_size+failure_size)//2
        if (OUT/f'inputs/queries/myosin_prefix_{midpoint}.json').exists() and time.monotonic()-gpu_start<gpu_budget-600:
            r=fold(f'gpu-size/n{midpoint}-midpoint',f'myosin_prefix_{midpoint}.json')
            if r['success']:largest=r;passed_size=midpoint
    for depth in (64,256,1024):
        if time.monotonic()-gpu_start>gpu_budget-600:break
        fold(f'gpu-msa/depth-{depth}',f'tmgh1_msa_{depth}.json',msa=depth)
    for samples in (2,4):
        if time.monotonic()-gpu_start>gpu_budget-600:break
        fold(f'gpu-samples/count-{samples}','myosin_prefix_256.json',samples=samples)
    for chains in (2,4):
        if time.monotonic()-gpu_start>gpu_budget-600:break
        r=fold(f'gpu-chains/copies-{chains}',f'tmgh1_copies_{chains}.json')
        if not r['success']:break
    # A second try at the failed size distinguishes deterministic resource failure from one transient.
    if failure_size and time.monotonic()-gpu_start<gpu_budget-600:
        fold(f'gpu-size/n{failure_size}-repeat',f'myosin_prefix_{failure_size}.json')
    if not largest:raise RuntimeError('No validated size-sweep prediction; GPU workload not ready for mixed stage')
    repeat=0
    while time.monotonic()-gpu_start<gpu_budget-60:
        remaining=gpu_budget-(time.monotonic()-gpu_start)
        expected=largest['resources']['wall_seconds']
        if remaining<expected*1.25:break
        r=fold(f'gpu-repeat/run-{repeat:03d}',Path(largest['query']).name,timeout=min(900,remaining));repeat+=1
        if not r['success']:break
    atomic(OUT/'selected-gpu-configuration.json',{'largest_success':largest,'first_failed_size':failure_size,'largest_passed_size':passed_size,'gpu_phase_seconds':time.monotonic()-gpu_start,'remaining_budget_not_filled':'Do not launch a job whose measured time exceeds remaining phase budget.'})
    vina('cpu-soak',workers,cpu,32,5400,panel='full',timeout=6000)
    mix_start=time.monotonic()
    mixed_cmd=[PY,str(ROOT/'stressbench/mixed_worker.py'),'--folder',str(OUT/'mixed-full/work'),'--duration','2700','--workers',str(workers),'--cpu',str(cpu),'--query',largest['query'],'--msa-depth','64','--sequence-mode','full']
    mixed=run_stage('mixed-full',mixed_cmd,3300)
    if (mixed['guard_stop'] or mixed['returncode']!=0) and mixed['guard_stop']!='operator_stop_file':
        remaining=max(300,2700-int(time.monotonic()-mix_start));event('mixed_fallback',cause=mixed['guard_stop'],remaining_seconds=remaining)
        fallback=[PY,str(ROOT/'stressbench/mixed_worker.py'),'--folder',str(OUT/'mixed-reduced/work'),'--duration',str(remaining),'--workers','3','--cpu','2','--query',str(OUT/'inputs/queries/myosin_prefix_256.json'),'--msa-depth','64','--sequence-mode','stream']
        run_stage('mixed-reduced',fallback,remaining+600)
    fault,summary,command=vina('recovery',3,2,64,600,cycles=2,inject=10,timeout=1000)
    # Rerun exactly the same finite plan, preserving completed outputs and recording cache hits.
    if summary and summary['failed']==0:
        resume=[x for x in command];idx=resume.index('--inject-after');del resume[idx:idx+2]
        run_stage('recovery-resume',resume,120)
        atomic(OUT/'recovery-verdict.json',{'fault_was_injected':(OUT/'recovery/batch/fault-injection.json').exists(),'summary':summary,'resume_summary':json.loads((OUT/'recovery/batch/resume-summary.json').read_text())})
    atomic(OUT/'completion.json',{'status':'completed','utc':utc(),'wall_seconds':time.monotonic()-start,'stages':len(stages)});progress('completed')
    p=json.loads((OUT/'progress.json').read_text());p['status']='completed';atomic(OUT/'progress.json',p);event('campaign_completed')
except BaseException as exc:
    stop_reason=f'{type(exc).__name__}: {exc}';event('campaign_stopped',reason=stop_reason)
    if active and active.poll() is None:
        (OUT/'STOP').write_text(stop_reason)
        try:active.wait(timeout=20)
        except subprocess.TimeoutExpired:pass
    atomic(OUT/'completion.json',{'status':'stopped','reason':stop_reason,'utc':utc(),'wall_seconds':time.monotonic()-start,'stages':len(stages)})
    atomic(OUT/'progress.json',{'status':'stopped','reason':stop_reason,'utc':utc(),'elapsed_seconds':time.monotonic()-start,'completed_stages':len(stages)})
finally:
    lock.unlink(missing_ok=True)
    caffeine.terminate()
    if (ROOT/'stressbench/report.py').exists():subprocess.run([PY,str(ROOT/'stressbench/report.py')],cwd=ROOT,check=False)
if stop_reason:raise SystemExit(1)
