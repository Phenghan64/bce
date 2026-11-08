"""Monitor the complete new throughput run without overwriting prior measurements."""
import json,subprocess,time,sys
from pathlib import Path
from datetime import datetime,timezone
import psutil

root=Path(__file__).resolve().parents[1];out=Path(sys.argv[1]).resolve()
assert not (out/'resources.json').exists(), 'Use run.py directly to resume; preserve original resource measurement'
assert not (out/'resources.jsonl').exists(), 'A measurement already exists; inspect before resuming'
start=time.perf_counter();before=psutil.swap_memory()._asdict();peak=0;available=psutil.virtual_memory().available
swap_peak=before['used'];samples=0;cpu_totals={}
command=[str(root/'run.sh'),'python',str(root/'largebench/run.py'),'--output-dir',str(out),'--stage','all']
with (out/'run.log').open('w') as log,(out/'resources.jsonl').open('w') as stream:
    process=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT)
    while process.poll() is None:
        rss=0;tree_cpu=0;process_count=0
        try:
            parent=psutil.Process(process.pid);children=[parent]+parent.children(recursive=True)
            for child in children:
                try:
                    rss+=child.memory_info().rss;t=child.cpu_times();cpu_totals[(child.pid,child.create_time())]=t.user+t.system;process_count+=1
                except psutil.Error:pass
        except psutil.Error:pass
        vm=psutil.virtual_memory();swap=psutil.swap_memory();peak=max(peak,rss);available=min(available,vm.available);swap_peak=max(swap_peak,swap.used)
        record={'elapsed_seconds':time.perf_counter()-start,'rss_sum_bytes':rss,'system_available_bytes':vm.available,
            'swap_used_bytes':swap.used,'swap_sout_bytes':swap.sout,'system_cpu_percent':psutil.cpu_percent(),
            'sampled_cumulative_tree_cpu_seconds':sum(cpu_totals.values()),'process_count':process_count}
        stream.write(json.dumps(record)+'\n');stream.flush();samples+=1;time.sleep(1)
r={'command':command,'seconds':time.perf_counter()-start,'returncode':process.returncode,'sample_count':samples,
    'process_tree_peak_rss_bytes':peak,'system_min_available_bytes':available,'swap_before':before,
    'swap_after':psutil.swap_memory()._asdict(),'swap_peak_used_bytes':swap_peak,'sampled_tree_cpu_seconds':sum(cpu_totals.values()),
    'note':'RSS summation can double count shared pages. CPU samples can miss short-lived processes. System memory/swap include other applications.'}
(out/'resources.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));sys.exit(process.returncode)
