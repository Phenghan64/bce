"""Use the remaining 30-minute supplemental GPU budget after n768 succeeds."""
import json, os, subprocess, time
from pathlib import Path
import psutil
from common import ROOT, OUT, atomic, utc, monitor
from supplement import read, validate, PY

def main():
    assert read(OUT/'supplement-completion.json',{}).get('status')=='completed'
    assert not (OUT/'STOP').exists()
    previous=read(OUT/'gpu-extended/n768/validation.json')
    assert previous['success']
    budget=1800-previous['resources']['wall_seconds']
    assert budget>60
    start=time.monotonic()
    os.environ['STRESSBENCH_OWNER_PID']=str(os.getpid())
    caffeine=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
    plan={'utc':utc(),'budget_seconds':budget,'reason':'768 completed with longer time budget; do not mistake it for the memory ceiling. Probe a doubled natural-sequence prefix using the unspent 1800-second supplemental GPU budget. This is an adaptive capacity bracket, not a biological claim or a precise residue threshold.','sequence_residues':1536,'same_settings_as_size_sweep':True,'repeat_failure_if_remaining_budget_allows':True}
    atomic(OUT/'boundary-followup-plan.json',plan)
    cases=[]
    try:
        for repeat in range(2):
            remaining=budget-(time.monotonic()-start)
            if repeat and (cases[-1]['success'] or remaining<max(60,cases[-1]['resources']['wall_seconds']*1.3)):
                break
            name='gpu-extended/n1536-boundary'+('-repeat' if repeat else '')
            folder=OUT/name/'model';query=OUT/'inputs/queries/myosin_prefix_1536.json'
            atomic(OUT/'boundary-progress.json',{'status':'running','phase':name,'controller_pid':os.getpid(),'started_monotonic':start,'utc':utc()})
            print(json.dumps({'phase':name,'remaining_seconds':remaining}),flush=True)
            command=[str(ROOT/'openfold3/run.sh'),'python',str(ROOT/'stressbench/fold_worker.py'),'--folder',str(folder),'--query',str(query),'--msa-depth','64','--samples','1','--fraction','1.0']
            r=monitor(command,OUT/name,remaining)
            c={'name':name,'query':str(query),'msa_depth':64,'samples':1,'fraction':1.0,'resources':r,'success':False,'adaptive_boundary_probe':True}
            b=read(folder/'benchmark.json')
            if b:c['benchmark']=b
            if r['returncode']==0 and not r['guard_stop']:
                try:c['validation']=validate(folder,query);c['success']=True
                except Exception as exc:c['validation_error']=repr(exc)
            if not c['success']:
                log=(OUT/name/'run.log').read_text(errors='replace')
                c['failure_class']='mps_allocator_oom' if 'MPS backend out of memory' in log or 'MPS out of memory' in log else r['guard_stop'] or 'software_or_validation_error'
                c['error_tail']=log[-2500:]
            atomic(OUT/name/'validation.json',c);cases.append(c)
            if r['guard_stop']=='operator_stop_file':raise InterruptedError('Operator stop')
            # Wait only for operational recovery, not for swap to be forcibly cleared.
            until=time.monotonic()+30
            while psutil.virtual_memory().available<2*2**30 and time.monotonic()<until:
                time.sleep(2)
        atomic(OUT/'boundary-completion.json',{'status':'completed','utc':utc(),'wall_seconds':time.monotonic()-start,'cases':len(cases),'supplemental_gpu_total_seconds':previous['resources']['wall_seconds']+sum(c['resources']['wall_seconds'] for c in cases)})
        atomic(OUT/'boundary-progress.json',{'status':'completed','phase':'completed','utc':utc()})
    except BaseException as exc:
        atomic(OUT/'boundary-completion.json',{'status':'stopped','utc':utc(),'reason':repr(exc),'wall_seconds':time.monotonic()-start})
        raise
    finally:
        caffeine.terminate()
        subprocess.run([PY,str(ROOT/'stressbench/report.py')],cwd=ROOT)

if __name__=='__main__':main()
