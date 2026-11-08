"""Observe Vina processes only; system swap includes all desktop applications."""
import json,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
out=Path(__file__).resolve().parents[1]/'results/ligand-library20-preflight-20261030'
started=time.monotonic();peak_rss=0;peak_cpu=0;peak_count=0;samples=0
before=subprocess.check_output(['/usr/sbin/sysctl','-n','vm.swapusage'],text=True).strip()
while not (out/'monitor.stop').exists():
    rows=[]
    for line in subprocess.check_output(['/bin/ps','-axo','pid=,rss=,%cpu=,comm='],text=True).splitlines():
        parts=line.split(None,3)
        if len(parts)==4 and Path(parts[3]).name=='vina':rows.append((int(parts[0]),int(parts[1])*1024,float(parts[2])))
    peak_rss=max(peak_rss,sum(r[1] for r in rows));peak_cpu=max(peak_cpu,sum(r[2] for r in rows));peak_count=max(peak_count,len(rows));samples+=1
    if samples%10==0:
        (out/'vina-process-observation.json').write_text(json.dumps({'elapsed_seconds':round(time.monotonic()-started,3),'samples':samples,'peak_observed_vina_processes':peak_count,'peak_observed_sum_vina_rss_bytes':peak_rss,'peak_observed_sum_vina_cpu_percent':peak_cpu,'system_swap_before':before,'scope':'Sampled active Vina processes only, 0.5s interval; not peak total system memory.'},indent=2)+'\n')
    time.sleep(.5)
