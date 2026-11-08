"""Sequential post-prediction workflow with stage timing and persisted failures."""
import subprocess,time,json,sys,os
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root/'results/tmgh1-round1';stages=[];start=time.perf_counter();env=dict(os.environ,OPENMM_CPU_THREADS='2');(out/'logs').mkdir(exist_ok=True)
def run(name,cmd):
 t=time.perf_counter();print('START',name,flush=True)
 with (out/'logs'/f'{name}.log').open('w') as log:r=subprocess.run([str(x) for x in cmd],stdout=log,stderr=subprocess.STDOUT,cwd=root,env=env)
 stages.append({'name':name,'command':list(map(str,cmd)),'seconds':time.perf_counter()-t,'returncode':r.returncode});(out/'postprediction-stages.json').write_text(json.dumps({'seconds_so_far':time.perf_counter()-start,'stages':stages},indent=2)+'\n');r.check_returncode();print('DONE',name,round(stages[-1]['seconds'],2),flush=True)
run('structure-analysis',[root/'openfold3/run.sh','python',root/'tmgh1/analyze_structures.py'])
for g in ['experimental_2J77','experimental_2J7E','predicted_42','predicted_2026','predicted_20261021']:
 for state in ['default','E166H']:run('prepare-'+g+'-'+state,[root/'.tools/prep-env/bin/python',root/'tmgh1/prepare_receptor.py',g,state])
run('receptor-audit',[root/'run.sh','python',root/'tmgh1/audit_preparation.py'])
run('docking',[root/'openfold3/run.sh','python',root/'tmgh1/monitor.py',out/'docking-resources',root/'run.sh','python',root/'tmgh1/dock.py'])
