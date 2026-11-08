"""Time a stage and sample its process-tree RSS plus system memory/swap."""
import json,subprocess,sys,time
from pathlib import Path
import psutil

root=Path(__file__).resolve().parents[1];stage=sys.argv[1];out=root/'results/abl-real-case'
if (out/f'{stage}-resources.json').exists():
    raise SystemExit('Measured resource report already exists; use run.py directly for resume checks.')
start=time.perf_counter();before=psutil.swap_memory()._asdict();peak=0;min_available=psutil.virtual_memory().available
command=[str(root/'run.sh'),'python',str(root/'realcase/run.py'),'--output-dir',str(out),'--stage',stage]
with (out/f'{stage}-console.log').open('w') as log:
    p=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT)
    while p.poll() is None:
        try:
            parent=psutil.Process(p.pid);rss=sum(q.memory_info().rss for q in [parent]+parent.children(recursive=True) if q.is_running());peak=max(peak,rss)
        except psutil.Error:pass
        min_available=min(min_available,psutil.virtual_memory().available);time.sleep(.5)
r={'stage':stage,'seconds':time.perf_counter()-start,'returncode':p.returncode,'sampled_process_tree_rss_peak_bytes':peak,
   'minimum_system_available_bytes':min_available,'swap_before':before,'swap_after':psutil.swap_memory()._asdict(),
   'note':'Process RSS sum can double-count shared pages; available RAM/swap are whole-system observations.'}
(out/f'{stage}-resources.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));sys.exit(p.returncode)
