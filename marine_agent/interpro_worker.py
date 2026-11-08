"""Authenticated Linux worker: pinned default analyses, durable jobs and resume.

No shell commands or filesystem paths are accepted from clients. Only one
Nextflow workflow runs at a time; member analyses retain their official memory
requirements. Docker containers are tagged for cancellation and crash recovery.
"""
from __future__ import annotations
import argparse
import base64
import csv
import hashlib
import hmac
import json
import os
from pathlib import Path
import signal
import subprocess
import threading
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from .interpro_protocol import PROTOCOL, JOB_ID, ANALYSIS_ID, digest, validate_records, validate_scan

TERMINAL = {'succeeded', 'failed', 'stopped'}

def atomic(path, value):
    tmp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2)); tmp.replace(path)

def file_hash(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576), b''): h.update(block)
    return h.hexdigest()

class Worker:
    def __init__(self, root, release):
        self.root, self.release_dir = Path(root), Path(release)
        self.release = json.loads((self.release_dir/'release.json').read_text())
        self.jobs = self.root/'jobs'; self.jobs.mkdir(parents=True,exist_ok=True)
        self.analyses = self.root/'analyses'; self.analyses.mkdir(exist_ok=True)
        self.lock=threading.RLock(); self.queue_lock=threading.Lock(); self.processes={}; self.threads={}
        # A restarted coordinator never silently repeats work. Clear only its
        # labelled containers and require a fresh attempt to resume the cache.
        for path in self.jobs.glob('*/state.json'):
            state=json.loads(path.read_text())
            if any(op['state'] not in TERMINAL for op in state.get('analyses',{}).values()):
                self._containers(state['job_id'], kill=True)
                state['cancelled']=True
                for op in state['analyses'].values():
                    if op['state'] not in TERMINAL: op.update(state='stopped',reason='worker restarted; completed stages retained')
                atomic(path,state)

    def _path(self, job):
        if not JOB_ID.fullmatch(job): raise ValueError('invalid job ID')
        directory=self.jobs/job; directory.mkdir(exist_ok=True)
        return directory/'state.json'

    def _read(self, job):
        path=self._path(job)
        return json.loads(path.read_text()) if path.exists() else {'job_id':job,'cancelled':False,'analyses':{}}

    def _save(self, state): atomic(self._path(state['job_id']),state)

    def _containers(self, job, kill=False):
        result=subprocess.run(['docker','ps','-q','--filter',f'label=eyoson.interpro.job={job}'],capture_output=True,text=True,timeout=20,check=True)
        ids=result.stdout.split()
        if ids and kill:
            subprocess.run(['docker','rm','-f',*ids],capture_output=True,timeout=30,check=True)
        return ids

    def status(self, job, analysis=None):
        with self.lock:
            state=self._read(job)
            active=any(t.is_alive() for (j,_),t in self.threads.items() if j==job)
            if analysis:
                op=dict(state['analyses'].get(analysis,{'state':'unknown'}))
                if op['state']=='running':
                    op['progress']=self.trace_summary(self.analyses/analysis/'trace.tsv')
                return {'job_id':job,'analysis_id':analysis,'active':active,**op}
            if active: status='cancelling' if state['cancelled'] else 'running'
            elif state['cancelled']: status='stopped'
            elif state['analyses']: status='failed' if any(o['state']=='failed' for o in state['analyses'].values()) else 'succeeded'
            else: status='unknown'
            return {'job_id':job,'state':status,'active':active,'analyses':state['analyses']}

    def cancel(self, job):
        with self.lock:
            state=self._read(job);state['cancelled']=True;self._save(state)
            for (j,_), process in list(self.processes.items()):
                if j==job:
                    try: os.killpg(process.pid,signal.SIGTERM)
                    except ProcessLookupError: pass
        self._containers(job,kill=True)
        return self.status(job)

    def submit(self, payload):
        if not isinstance(payload, dict) or set(payload) - {'job_id','analysis_id','release_digest','records','refresh'} or type(payload.get('refresh', False)) is not bool:
            raise ValueError('invalid scan request fields')
        job, analysis=payload.get('job_id',''),payload.get('analysis_id','')
        if not JOB_ID.fullmatch(job) or not ANALYSIS_ID.fullmatch(analysis):raise ValueError('invalid execution identifiers')
        if payload.get('release_digest') != digest(self.release):raise ValueError('execution release changed; resolve inputs again')
        records=validate_records(payload.get('records'))
        request={'records':records,'release_digest':payload['release_digest']}
        request_hash=digest(request)
        with self.lock:
            state=self._read(job)
            if state['cancelled']:raise ValueError('job was cancelled; use a new attempt ID')
            previous=state['analyses'].get(analysis)
            if previous:
                if previous['request_hash']!=request_hash:raise ValueError('same analysis with different inputs')
                return self.status(job,analysis)
            directory=self.analyses/analysis;directory.mkdir(exist_ok=True)
            spec=directory/'request.json'
            if spec.exists() and json.loads(spec.read_text())!=request:raise ValueError('analysis identity collision')
            atomic(spec,request)
            state['analyses'][analysis]={'state':'queued','request_hash':request_hash,'submitted_at':time.time()};self._save(state)
            thread=threading.Thread(target=self._run,args=(job,analysis,records,bool(payload.get('refresh',False))),daemon=True)
            self.threads[job,analysis]=thread;thread.start()
            return self.status(job,analysis)

    def _update(self, job, analysis, **fields):
        with self.lock:
            state=self._read(job);state['analyses'][analysis].update(fields);self._save(state)

    def _cancelled(self,job):
        with self.lock:return self._read(job)['cancelled']

    @staticmethod
    def trace_summary(path):
        if not path.exists():return {'completed_processes':0,'cached_processes':0}
        try:
            with path.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
            return {'completed_processes':sum(r.get('status')=='COMPLETED' for r in rows),
                    'cached_processes':sum(r.get('status')=='CACHED' for r in rows),
                    'failed_processes':sum(r.get('status')=='FAILED' for r in rows),
                    'last_processes':[{'name':r.get('name'),'status':r.get('status')} for r in rows[-6:]]}
        except (OSError,ValueError):return {'completed_processes':0,'cached_processes':0}

    def _run(self,job,analysis,records,refresh):
        started=time.monotonic();acquired=False;process=None;directory=self.analyses/analysis
        try:
            while not acquired:
                if self._cancelled(job):return self._update(job,analysis,state='stopped')
                acquired=self.queue_lock.acquire(timeout=.25)
            if self._cancelled(job):return self._update(job,analysis,state='stopped')
            images=self.release.get('images', [])
            if images:
                tags=[row['tags'][0] for row in images]
                actual=json.loads(subprocess.check_output(['docker','image','inspect',*tags],text=True,timeout=30))
                if [row['Id'] for row in actual] != [row['id'] for row in images]:
                    raise ValueError('container image identity changed; install a new verified release')
            receipt=directory/'receipt.json'
            if not refresh and receipt.exists():
                saved=json.loads(receipt.read_text())
                if all((directory/name).is_file() and file_hash(directory/name)==h for name,h in saved['files'].items()):
                    self._update(job,analysis,state='succeeded',cache_hit=True,elapsed_seconds=time.monotonic()-started);return
            fasta = ''.join(f">{r['id']}\n{r['sequence']}\n" for r in records)
            input_path = directory/'input.faa'
            # Nextflow's default file cache includes mtime. An unchanged input
            # must keep its timestamp for cancellation/resume to reuse tasks.
            if not input_path.exists() or input_path.read_text() != fasta:
                input_path.write_text(fasta)
            # Keep prior attempt logs; the work/cache directories stay in place.
            history=directory/'attempts'/job;history.mkdir(parents=True,exist_ok=True)
            for name in ('trace.tsv','nextflow.log','console.log','execution-report.html','scan.json','scan.tsv','receipt.json'):
                path=directory/name
                if path.exists():path.replace(history/name)
            cpus=int(os.getenv('INTERPRO_TASK_CPUS','2'));workers=int(os.getenv('INTERPRO_MAX_WORKERS','2'))
            config=directory/'execution.config'
            config.write_text(f"executor.cpus = {workers}\nexecutor.memory = '18 GB'\nprocess.maxForks = {workers}\ndocker.runOptions = '--label eyoson.interpro.job={job}'\n")
            env=os.environ.copy();env['JAVA_HOME']=str(self.root/'runtime/java');env['NXF_HOME']=str(self.root/'runtime/nxf');env['NXF_VER']='26.10.3';env['NXF_OPTS']='-Xms128m -Xmx1536m';env['NXF_ANSI_LOG']='false'
            env.pop('INTERPRO_WORKER_TOKEN', None)
            env['PATH']=env['JAVA_HOME']+'/bin:'+env['PATH']
            command=[str(self.root/'runtime/nextflow'),'-log','nextflow.log','run',str(self.release_dir/'interproscan6-6.0.3.0/main.nf'),'-profile','docker','-c',str(config),'-w',str(directory/'work'),'--input',str(directory/'input.faa'),'--datadir',str(self.release_dir/'databases'),'--outdir',str(directory),'--outprefix','scan','--formats','json,tsv','--interpro',self.release['interpro_version'],'--no-matches-api','--cpus',str(cpus),'--max-workers',str(workers),'-with-trace','trace.tsv','-with-report','execution-report.html']
            if not refresh and (directory/'.nextflow/cache').exists():command.append('-resume')
            atomic(directory/'command.json',{'argv':command,'scope':'official_default_available_analyses','matches_api':False})
            self._update(job,analysis,state='running',started_at=time.time(),cache_hit=False)
            with (directory/'console.log').open('w') as log:
                with self.lock:
                    if self._cancelled(job):return self._update(job,analysis,state='stopped')
                    process=subprocess.Popen(command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                    self.processes[job,analysis]=process
                cancelled_at=None
                while process.poll() is None:
                    if time.monotonic()-started > int(os.getenv('INTERPRO_TIMEOUT_SECONDS','86400')):
                        self.cancel(job)
                        self._update(job,analysis,reason='execution time budget exceeded; completed stages retained')
                    if self._cancelled(job):
                        if cancelled_at is None:cancelled_at=time.monotonic()
                        try:os.killpg(process.pid,signal.SIGKILL if time.monotonic()-cancelled_at>3 else signal.SIGTERM)
                        except ProcessLookupError:pass
                    time.sleep(.5)
                code=process.wait()
            if self._cancelled(job):
                self._containers(job,kill=True);self._update(job,analysis,state='stopped',elapsed_seconds=time.monotonic()-started);return
            if code:raise RuntimeError(f'InterProScan/Nextflow exited {code}; inspect console.log and nextflow.log')
            raw=json.loads((directory/'scan.json').read_text());validate_scan(raw,records,self.release)
            files={n:file_hash(directory/n) for n in ['scan.json','scan.tsv','trace.tsv','nextflow.log','console.log','command.json'] if (directory/n).is_file()}
            if not {'scan.json','scan.tsv','trace.tsv','command.json'}<=files.keys():raise ValueError('missing required execution artifacts')
            data={'analysis_id':analysis,'release':self.release,'files':files,'progress':self.trace_summary(directory/'trace.tsv'),'elapsed_seconds':time.monotonic()-started,'exit_code':code,'de_novo_scan_executed':True}
            atomic(receipt,data)
            self._update(job,analysis,state='succeeded',elapsed_seconds=data['elapsed_seconds'],progress=data['progress'])
        except Exception as exc:
            self._update(job,analysis,state='stopped' if self._cancelled(job) else 'failed',reason=f'{type(exc).__name__}: {exc}',elapsed_seconds=time.monotonic()-started)
        finally:
            if process and process.poll() is None:
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                process.wait()
            with self.lock:self.processes.pop((job,analysis),None)
            if acquired:self.queue_lock.release()

    def artifacts(self,job,analysis):
        status=self.status(job,analysis)
        if status['state']!='succeeded':raise ValueError('scan artifacts are not complete')
        directory=self.analyses/analysis;receipt=json.loads((directory/'receipt.json').read_text())
        if receipt['release']!=self.release or receipt['analysis_id']!=analysis:raise ValueError('receipt identity mismatch')
        files={}
        for name,h in receipt['files'].items():
            path=directory/name
            if file_hash(path)!=h:raise ValueError('artifact integrity failure')
            files[name]=base64.b64encode(path.read_bytes()).decode()
        return {'job_id':job,'analysis_id':analysis,'receipt':receipt,'files':files}

class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt,*args):pass
    def respond(self,code,value):
        body=json.dumps(value,ensure_ascii=False).encode();self.send_response(code);self.send_header('Content-Type','application/json');self.send_header('Content-Length',str(len(body)));self.end_headers()
        try:self.wfile.write(body)
        except (BrokenPipeError,ConnectionResetError):pass
    def handle_request(self):
        if not hmac.compare_digest(self.headers.get('X-InterPro-Worker-Token',''),self.server.token):return self.respond(401,{'error':'unauthorized'})
        import re
        parts=self.path.strip('/').split('/')
        try:
            worker=self.server.worker
            if self.command=='GET' and self.path=='/health':return self.respond(200,{'status':'ready','protocol':PROTOCOL,'release':worker.release,'release_digest':digest(worker.release)})
            if self.command=='POST':
                length=int(self.headers.get('Content-Length','0'))
                if not 0<=length<=40*1024*1024:raise ValueError('request exceeds transport limit')
                payload=json.loads(self.rfile.read(length) or b'{}')
                if not isinstance(payload,dict):raise ValueError('request must be a JSON object')
                if self.path=='/v1/scan':return self.respond(202,worker.submit(payload))
                if len(parts)==4 and parts[:2]==['v1','jobs'] and parts[3]=='cancel':return self.respond(200,worker.cancel(parts[2]))
            if self.command=='GET' and parts[:2]==['v1','jobs']:
                if len(parts)==3:return self.respond(200,worker.status(parts[2]))
                if len(parts) in (5,6) and parts[3]=='analyses' and ANALYSIS_ID.fullmatch(parts[4]):
                    if len(parts)==5:return self.respond(200,worker.status(parts[2],parts[4]))
                    if parts[5]=='artifacts':return self.respond(200,worker.artifacts(parts[2],parts[4]))
            self.respond(404,{'error':'not found'})
        except (ValueError,KeyError,TypeError) as exc:self.respond(400,{'error':str(exc)})
        except Exception:self.respond(503,{'error':'worker operation failed; inspect worker logs'})
    do_GET=handle_request
    do_POST=handle_request

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--release',required=True);p.add_argument('--host',default='100.64.0.15');p.add_argument('--port',type=int,default=8110);p.add_argument('--cleanup-only',action='store_true');a=p.parse_args()
    if a.cleanup_only:
        Worker(a.root,a.release)
        return
    token=os.environ.get('INTERPRO_WORKER_TOKEN','')
    if len(token)<32:raise RuntimeError('configure worker authentication')
    server=ThreadingHTTPServer((a.host,a.port),Handler);server.daemon_threads=True;server.token=token;server.worker=Worker(a.root,a.release);server.serve_forever()
if __name__=='__main__':main()
