import json,time
from pathlib import Path
from common import OUT
def read(p):
    try:return json.loads(p.read_text())
    except (FileNotFoundError,json.JSONDecodeError):return {}
p=read(OUT/'progress.json')
supp=read(OUT/'supplement-progress.json')
supp_done=read(OUT/'supplement-completion.json')
if supp_done:supp={**supp,**supp_done}
boundary=read(OUT/'boundary-progress.json')
if boundary:
    supp={**boundary,**read(OUT/'boundary-completion.json')}
if p.get('status')!='running' and supp:
    p={**p,**supp}
phase=p.get('phase','unknown');folder=OUT/phase
r=read(folder/'latest.json');b=read(folder/'batch/progress.json');mix=read(folder/'work/progress.json')
summary={'status':p.get('status'),'phase':phase,'campaign_minutes':round(p.get('elapsed_seconds',0)/60,1),'completed_stages':p.get('completed_stages'),'stage_minutes':round(r.get('elapsed_seconds',0)/60,1),'CPU_percent':round(r.get('system_cpu_percent',0),1),'GPU_percent':r.get('gpu_system',{}).get('Device Utilization %'),'RSS_GiB':round(r.get('process_tree_rss_bytes',0)/2**30,2),'available_GiB':round(r.get('system_available_bytes',0)/2**30,2),'swap_GiB':round(r.get('swap_used_bytes',0)/2**30,2),'thermal_state':r.get('thermal_state')}
if b:summary['docking']={k:b[k] for k in ['completed','failed','cycle','workers','cpu','exhaustiveness']}
if mix:summary['mixed']=mix
if p.get('reason'):summary['reason']=p['reason']
if not r:
    for key in ('CPU_percent','GPU_percent','RSS_GiB','available_GiB','swap_GiB','thermal_state'):
        summary[key]=None
if supp.get('status')=='running' and supp.get('started_monotonic'):
    summary['supplement_minutes']=round((time.monotonic()-supp['started_monotonic'])/60,1)
print(json.dumps(summary,ensure_ascii=False))
