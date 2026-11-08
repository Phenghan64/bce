"""Concurrent workloads; all children belong to the externally monitored test group."""
import argparse,json,os,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=Path(os.environ.get('STRESSBENCH_OUTPUT',ROOT/'results/stress-20261022-main'))
p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--duration',type=int,required=True);p.add_argument('--workers',type=int,required=True);p.add_argument('--cpu',type=int,required=True);p.add_argument('--query',required=True);p.add_argument('--msa-depth',type=int,default=64);p.add_argument('--sequence-mode',choices=['full','stream'],default='full');args=p.parse_args();out=Path(args.folder);out.mkdir(parents=True,exist_ok=True)
start=time.monotonic();children={};rounds={'fold':0,'sequence':0};results=[]
def launch(kind):
    if kind=='vina':cmd=[str(ROOT/'.env/bin/python'),str(ROOT/'stressbench/vina_batch.py'),'--folder',str(out/'vina'),'--workers',str(args.workers),'--cpu',str(args.cpu),'--exhaustiveness','32','--duration',str(args.duration),'--panel','full']
    elif kind=='fold':
        folder=out/f'fold-{rounds[kind]:03d}';cmd=[str(ROOT/'openfold3/run.sh'),'python',str(ROOT/'stressbench/fold_worker.py'),'--folder',str(folder),'--query',args.query,'--msa-depth',str(args.msa_depth)];rounds[kind]+=1
    else:
        folder=out/f'sequence-{rounds[kind]:03d}';cmd=[str(ROOT/'.tools/enzyme-case/venv/bin/python'),str(ROOT/'enzyme_case/run.py')] if args.sequence_mode=='full' else [str(ROOT/'.tools/enzyme-case/venv/bin/python'),str(ROOT/'stressbench/sequence_stream.py'),str(folder)];rounds[kind]+=1
    env=os.environ.copy()
    if kind=='sequence':env['ENZYME_CASE_OUTPUT']=str(folder)
    log=(out/f'{kind}-launch-{rounds.get(kind,0):03d}.log').open('w');process=subprocess.Popen(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=env);children[kind]=(process,log,time.monotonic(),cmd)
for kind in ('vina','fold','sequence'):launch(kind)
while children:
    for kind,(proc,log,when,cmd) in list(children.items()):
        if proc.poll() is not None:
            log.close();results.append({'kind':kind,'returncode':proc.returncode,'seconds':time.monotonic()-when,'command':cmd});del children[kind]
            if kind=='fold' and proc.returncode!=0 and 'myosin_prefix_256.json' not in args.query:
                # Do not silently finish the stage as CPU-only after a model failure.
                original_query=args.query;args.query=str(OUT/'inputs/queries/myosin_prefix_256.json');args.msa_depth=64
                (out/'fold-fallback.json').write_text(json.dumps({'from_query':original_query,'to_query':args.query,'reason':'fold_worker_failed_under_concurrent_load','seconds':time.monotonic()-start},indent=2)+'\n')
                if time.monotonic()-start<args.duration:launch(kind)
            elif kind!='vina' and proc.returncode==0 and time.monotonic()-start<args.duration:launch(kind)
            elif proc.returncode!=0:
                (out/'summary.json').write_text(json.dumps({'wall_seconds':time.monotonic()-start,'results':results,'arguments':vars(args),'failure':'concurrent_worker_failed; terminate owned stage group'},indent=2)+'\n')
                os.killpg(os.getpgrp(),signal.SIGTERM)
    (out/'progress.json').write_text(json.dumps({'seconds':time.monotonic()-start,'rounds_launched':rounds,'completed':len(results),'failed':sum(r['returncode']!=0 for r in results),'active':list(children)})+'\n')
    if children:time.sleep(2)
(out/'summary.json').write_text(json.dumps({'wall_seconds':time.monotonic()-start,'results':results,'arguments':vars(args)},indent=2)+'\n')
