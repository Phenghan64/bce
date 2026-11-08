#!/usr/bin/env python3
"""One-command reproducible engineering and scientific verification.

Default: backend/frontend tests, TypeScript and production frontend build.
--scientific adds a fresh real 7L39 OpenFold3 -> Vina reference run.
No production service is changed by this command.
"""
from __future__ import annotations
import argparse,json,os,re,shlex,subprocess,time
from pathlib import Path
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parent
FRONTEND=Path(os.getenv('DOHO_FRONTEND_DIR','/Users/yili/Documents/dsh/doho'))
NODE_ROOT=Path('/Users/yili/Library/Application Support/BridgeDevTools/2026-10-14/node/bin')
NODE=NODE_ROOT/'node' if (NODE_ROOT/'node').is_file() else Path(os.getenv('NODE_BINARY','node'))
NPM=NODE_ROOT/'npm' if (NODE_ROOT/'npm').is_file() else Path(os.getenv('NPM_BINARY','npm'))

def run_step(name,command,cwd,log,env=None):
    started=time.perf_counter();print(f'\n[{name}] {shlex.join(map(str,command))}',flush=True)
    with log.open('w') as stream:
        process=subprocess.Popen([str(x) for x in command],cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,env=env)
        assert process.stdout
        for line in process.stdout:
            stream.write(line);print(line,end='',flush=True)
        process.stdout.close()
        code=process.wait()
    return {'name':name,'command':[str(x) for x in command],'cwd':str(cwd),'log':str(log),
            'exit_code':code,'elapsed_seconds':round(time.perf_counter()-started,3)}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--scientific',action='store_true',help='run a fresh real 7L39 OpenFold3/Vina reference calculation')
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    out=args.output_dir or ROOT/'results'/'verification'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out.mkdir(parents=True,exist_ok=False);out=out.resolve()
    backend=ROOT/'.env/bin/python3.11'
    if not backend.is_file():raise SystemExit('Missing project Python: '+str(backend))
    if not (FRONTEND/'package.json').is_file():raise SystemExit('Missing frontend: '+str(FRONTEND))
    steps=[]
    frontend_env=os.environ.copy()
    if NODE_ROOT.is_dir():frontend_env['PATH']=str(NODE_ROOT)+os.pathsep+frontend_env.get('PATH','')
    commands=[
      ('backend-tests',[backend,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'],ROOT,None),
      ('frontend-cjs-tests',[NODE,'--test',*sorted(str(p) for p in (FRONTEND/'tests').glob('*.test.cjs'))],FRONTEND,frontend_env),
      ('frontend-ts-tests',[NODE,'--import','tsx','--experimental-specifier-resolution=node','--test',*sorted(str(p) for p in (FRONTEND/'tests').glob('*.test.ts')),*sorted(str(p) for p in (FRONTEND/'tests').glob('*.test.tsx'))],FRONTEND,frontend_env),
      ('frontend-typecheck',[NPM,'run','lint'],FRONTEND,frontend_env),
    ]
    build_env=frontend_env.copy();build_env.update(VITE_BASE_PATH='/marine-agent/',VITE_MARINE_MODE='true',VITE_MARINE_AGENT_API_BASE='/doho/research-api')
    commands.append(('frontend-build',[NPM,'run','build'],FRONTEND,build_env))
    if args.scientific:
        commands.append(('real-7l39-openfold-vina',[ROOT/'run.sh','python',ROOT/'workflow/run.py','--output-dir',out/'7l39'],ROOT,None))
    failed=None
    for name,command,cwd,env in commands:
        step=run_step(name,command,cwd,out/(name+'.log'),env)
        steps.append(step)
        if step['exit_code']:
            failed=name;break
    if not failed and args.scientific:
        report=json.loads((out/'7l39'/'report.json').read_text())
        if not report.get('execution_checks_passed') or not report.get('benchmark_passed'):failed='scientific-reference-contract'
    summary={'status':'failed' if failed else 'passed','failed_step':failed,'scientific_run':args.scientific,
             'steps':steps,'created_at':datetime.now(timezone.utc).isoformat()}
    summary['backend_tests']=sum(int(x) for step in steps if step['name']=='backend-tests' for x in re.findall(r'Ran (\d+) tests',(out/(step['name']+'.log')).read_text()))
    summary['frontend_tests']=sum(int(x) for step in steps if step['name'] in {'frontend-cjs-tests','frontend-ts-tests'} for x in re.findall(r'ℹ tests (\d+)',(out/(step['name']+'.log')).read_text()))
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print('\n'+json.dumps(summary,ensure_ascii=False,indent=2))
    return 1 if failed else 0
if __name__=='__main__':raise SystemExit(main())
