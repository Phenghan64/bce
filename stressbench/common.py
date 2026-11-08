"""Owned-process monitoring for the local capacity experiment."""
import ctypes,datetime,hashlib,json,os,re,signal,subprocess,time
from pathlib import Path
import psutil
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('STRESSBENCH_OUTPUT',ROOT/'results/stress-20261022-main')).resolve()
OUT.mkdir(parents=True,exist_ok=True)
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def atomic(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n');tmp.replace(path)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def thermal_state():
    try:
        ctypes.CDLL('/System/Library/Frameworks/Foundation.framework/Foundation')
        objc=ctypes.CDLL('/usr/lib/libobjc.A.dylib')
        objc.objc_getClass.argtypes=[ctypes.c_char_p];objc.objc_getClass.restype=ctypes.c_void_p
        objc.sel_registerName.argtypes=[ctypes.c_char_p];objc.sel_registerName.restype=ctypes.c_void_p
        send=ctypes.CFUNCTYPE(ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p)(('objc_msgSend',objc))
        obj=send(objc.objc_getClass(b'NSProcessInfo'),objc.sel_registerName(b'processInfo'))
        send_int=ctypes.CFUNCTYPE(ctypes.c_long,ctypes.c_void_p,ctypes.c_void_p)(('objc_msgSend',objc))
        return int(send_int(obj,objc.sel_registerName(b'thermalState')))
    except Exception:return None
def gpu_stats():
    try:
        s=subprocess.check_output(['ioreg','-r','-d','1','-c','AGXAccelerator'],text=True,timeout=3)
        line=next(x for x in s.splitlines() if '"PerformanceStatistics"' in x)
        return {k:int(v) for k,v in re.findall(r'"([^"]+)"=(\d+)',line)}
    except Exception:return {}
def terminate_group(p):
    try:os.killpg(p.pid,signal.SIGTERM)
    except ProcessLookupError:return
    try:p.wait(timeout=8)
    except subprocess.TimeoutExpired:
        try:os.killpg(p.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        p.wait()
def monitor(command,folder,timeout_seconds,allow_pressure=True):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    if (folder/'resources.jsonl').exists():raise RuntimeError('Refusing to overwrite monitoring evidence: '+str(folder))
    start=time.monotonic();initial=psutil.swap_memory();cpu_totals={};peak=0;minavail=2**63;swapmax=initial.used;low_since=None;critical_since=None;stop=None
    psutil.cpu_percent(percpu=True)
    atomic(folder/'launch.json',{'command':command,'utc':utc(),'timeout_seconds':timeout_seconds,'monitor_pid':os.getpid()})
    with (folder/'run.log').open('w') as log,(folder/'resources.jsonl').open('w') as stream:
        process=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        atomic(folder/'pid.json',{'pid':process.pid,'process_group':process.pid})
        lastcpu=0.;lastelapsed=0.;gpu={};n=0
        while process.poll() is None:
            elapsed=time.monotonic()-start;rss=0
            try:
                children=[psutil.Process(process.pid),*psutil.Process(process.pid).children(recursive=True)]
                for child in children:
                    try:
                        rss+=child.memory_info().rss;c=child.cpu_times();cpu_totals[(child.pid,child.create_time())]=c.user+c.system
                    except psutil.Error:pass
            except psutil.Error:pass
            vm=psutil.virtual_memory();swap=psutil.swap_memory();therm=thermal_state()
            if n%5==0:gpu=gpu_stats()
            cpus=psutil.cpu_percent(percpu=True);totalcpu=sum(cpu_totals.values());owncpu=100*(totalcpu-lastcpu)/max(elapsed-lastelapsed,.01)
            sample={'utc':utc(),'elapsed_seconds':elapsed,'process_tree_rss_bytes':rss,'process_tree_cpu_seconds_sampled':totalcpu,'process_tree_cpu_percent_one_core_100':owncpu,'system_cpu_per_core_percent':cpus,'system_cpu_percent':sum(cpus)/len(cpus),'system_available_bytes':vm.available,'swap_used_bytes':swap.used,'swap_sin_bytes':swap.sin,'swap_sout_bytes':swap.sout,'thermal_state':therm,'gpu_system':gpu,'disk_free_bytes':psutil.disk_usage(str(OUT)).free}
            stream.write(json.dumps(sample)+'\n');stream.flush();atomic(folder/'latest.json',sample)
            peak=max(peak,rss);minavail=min(minavail,vm.available);swapmax=max(swapmax,swap.used)
            if vm.available<384*2**20:low_since=low_since or time.monotonic()
            else:low_since=None
            if therm==3:critical_since=critical_since or time.monotonic()
            else:critical_since=None
            if elapsed>timeout_seconds:stop='stage_wall_time_limit'
            elif (swap.used>8*2**30 and swap.used-initial.used>256*2**20) or swap.used-initial.used>6*2**30:stop='system_swap_operational_limit'
            elif low_since and time.monotonic()-low_since>12:stop='system_available_below_384MiB_for_12s'
            elif critical_since and time.monotonic()-critical_since>30:stop='system_thermal_critical_30s'
            elif sample['disk_free_bytes']<20*2**30:stop='disk_free_below_20GiB'
            elif (OUT/'STOP').exists():stop='operator_stop_file'
            elif os.environ.get('STRESSBENCH_OWNER_PID') and not psutil.pid_exists(int(os.environ['STRESSBENCH_OWNER_PID'])):stop='supervisor_disappeared'
            if stop:terminate_group(process);break
            lastcpu=totalcpu;lastelapsed=elapsed;n+=1;time.sleep(2)
    result={'command':command,'started_utc':json.loads((folder/'launch.json').read_text())['utc'],'finished_utc':utc(),'wall_seconds':time.monotonic()-start,'returncode':process.returncode,'guard_stop':stop,'process_tree_peak_rss_bytes':peak,'system_min_available_bytes':minavail,'swap_before_bytes':initial.used,'swap_peak_bytes':swapmax,'swap_after_bytes':psutil.swap_memory().used,'sampled_tree_cpu_seconds':sum(cpu_totals.values()),'notes':['RSS and Metal allocations overlap: do not add.','System metrics include other apps.','CPU sampling may miss short-lived workers.','Thermal state 0 nominal, 1 fair, 2 serious, 3 critical; not a temperature.','Allocator caps and watchdog stops are operational boundaries, not proof of physical maximum.']}
    atomic(folder/'resources.json',result);return result

if __name__=='__main__':
    import sys
    r=monitor(sys.argv[3:],OUT/sys.argv[1],float(sys.argv[2]));print(json.dumps(r,indent=2));sys.exit(0 if r['returncode']==0 else 1)
