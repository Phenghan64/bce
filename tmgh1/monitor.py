"""Monitor a command and preserve bounded resource-stop evidence."""
import sys,os,signal,json,time,subprocess
from pathlib import Path
import psutil
out=Path(sys.argv[1]).resolve();cmd=sys.argv[2:];out.mkdir(parents=True,exist_ok=True);assert not (out/'resources.jsonl').exists()
start=time.perf_counter();before=psutil.swap_memory()._asdict();peak=0;available=psutil.virtual_memory().available;low_since=None;reason=None;swappeak=before['used'];cpu_totals={};psutil.cpu_percent()
with (out/'run.log').open('w') as log,(out/'resources.jsonl').open('w') as stream:
 p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 while p.poll() is None:
  rss=0
  try:
   for child in [psutil.Process(p.pid),*psutil.Process(p.pid).children(recursive=True)]:
    try:
     rss+=child.memory_info().rss;t=child.cpu_times();cpu_totals[(child.pid,child.create_time())]=t.user+t.system
    except psutil.Error:pass
  except psutil.Error:pass
  elapsed=time.perf_counter()-start;vm=psutil.virtual_memory();swap=psutil.swap_memory();peak=max(peak,rss);available=min(available,vm.available);swappeak=max(swappeak,swap.used)
  rec={'elapsed_seconds':elapsed,'rss_sum_bytes':rss,'system_available_bytes':vm.available,'system_cpu_percent':psutil.cpu_percent(),'swap_used_bytes':swap.used,'swap_sout_bytes':swap.sout,'sampled_tree_cpu_seconds':sum(cpu_totals.values())};stream.write(json.dumps(rec)+'\n');stream.flush()
  if vm.available<512*2**20:low_since=low_since or time.perf_counter()
  else:low_since=None
  if elapsed>3600:reason='One-hour runtime limit'
  elif swap.used-before['used']>4*2**30:reason='System swap grew by more than4 GiB'
  elif low_since and time.perf_counter()-low_since>15:reason='System available memory below512 MiB for15 seconds'
  if reason:
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
   break
  time.sleep(1)
r={'command':cmd,'elapsed_seconds':time.perf_counter()-start,'returncode':p.returncode,'abort_reason':reason,'process_tree_peak_RSS_bytes':peak,'system_min_available_bytes':available,'swap_before':before,'swap_after':psutil.swap_memory()._asdict(),'swap_peak_used_bytes':swappeak,'sampled_tree_CPU_seconds':sum(cpu_totals.values()),'note':'RSS sum can double-count shared pages; system metrics include other apps; not a temperature measurement.'};(out/'resources.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));sys.exit(p.returncode or (1 if reason else 0))
