"""MPS inference worker with live allocator evidence, without changing model weights."""
import argparse,json,os,resource,sys,threading,time
from pathlib import Path
import psutil,torch
from importlib.metadata import version
p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--query',required=True);p.add_argument('--samples',type=int,default=1);p.add_argument('--msa-depth',type=int,default=64);p.add_argument('--fraction',type=float,default=1.0);args=p.parse_args()
out=Path(args.folder);out.mkdir(parents=True,exist_ok=True)
if not torch.backends.mps.is_available():raise RuntimeError('MPS is unavailable; CPU fallback refused')
torch.set_num_threads(6);torch.mps.set_per_process_memory_fraction(args.fraction)
runner=out/'runner.json'
config={'pl_trainer_args':{'accelerator':'mps','devices':1,'precision':'32-true','cluster_environment':'local'},'data_module_args':{'num_workers':0,'num_workers_validation':0},'experiment_settings':{'seeds':[42],'use_templates':False,'use_msa_server':False},'dataset_config_kwargs':{'msa':{'max_rows':args.msa_depth}},'model_update':{'presets':['predict','low_mem'],'custom':{'architecture':{'msa':{'msa_module_embedder':{'min_subsampled_all_msa':args.msa_depth,'max_subsampled_all_msa':args.msa_depth}}}}}}
# Preserve the upstream nested configuration structure.
runner.write_text(json.dumps(config,indent=2)+'\n')
start=time.monotonic();done=threading.Event();peak_driver=0;peak_tensor=0;error=None;ok=False
def sample():
    global peak_driver,peak_tensor
    with (out/'mps-memory.jsonl').open('w') as stream:
        while not done.is_set():
            try:
                driver=torch.mps.driver_allocated_memory();tensor=torch.mps.current_allocated_memory();peak_driver=max(peak_driver,driver);peak_tensor=max(peak_tensor,tensor)
                s={'seconds':time.monotonic()-start,'driver_bytes':driver,'tensor_bytes':tensor,'rss_bytes':psutil.Process().memory_info().rss,'swap_used_bytes':psutil.swap_memory().used}
                stream.write(json.dumps(s)+'\n');stream.flush()
            except Exception:pass
            done.wait(1)
thread=threading.Thread(target=sample,daemon=True);thread.start()
argv=['predict','--query-json',str(Path(args.query).resolve()),'--runner-yaml',str(runner),'--num-diffusion-samples',str(args.samples),'--use-msa-server','false','--use-templates','false','--output-dir',str(out/'prediction')]
try:
    from openfold3.run_openfold import cli
    cli.main(args=argv,standalone_mode=False);torch.mps.synchronize();ok=True
except BaseException as exc:error=f'{type(exc).__name__}: {exc}';raise
finally:
    done.set();thread.join(timeout=3)
    report={'completed':ok,'error':error,'cli_seconds':time.monotonic()-start,'openfold_version':version('openfold3'),'torch_version':torch.__version__,'backend':'PyTorch MPS (not MLX)','mps_memory_fraction':args.fraction,'mps_recommended_max_bytes':torch.mps.recommended_max_memory(),'sampled_mps_driver_peak_bytes':peak_driver,'sampled_mps_tensor_peak_bytes':peak_tensor,'process_maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'samples_requested':args.samples,'msa_depth_requested':args.msa_depth,'query':args.query,'command_args':argv}
    (out/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
