"""Full InterProScan adapter with exact-input, release and artifact verification."""
from __future__ import annotations
import base64
import fcntl
import json
import os
import time
import uuid
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, build_opener, ProxyHandler
from .execution_control import RunCancelled, check_cancelled, current_control
from .interpro_protocol import PROTOCOL, digest, validate_scan, active_library_versions
from .protein_screening import atomic_json, file_sha, sha

FILES = {'scan.json','scan.tsv','trace.tsv','nextflow.log','console.log','command.json'}

def configured():
    return bool(os.getenv('INTERPRO_WORKER_URL') and os.getenv('INTERPRO_WORKER_TOKEN'))

class RemoteScanStopPending(RunCancelled):
    remote_stop_pending = True

def request(path, payload=None, *, timeout=20, attempts=3):
    url=os.getenv('INTERPRO_WORKER_URL','').rstrip('/')
    token=os.getenv('INTERPRO_WORKER_TOKEN','')
    if not url or not token:raise RuntimeError('InterProScan 全扫描执行节点未配置')
    headers={'X-InterPro-Worker-Token':token,'Content-Type':'application/json'}
    opener=build_opener(ProxyHandler({}))
    for attempt in range(attempts):
        try:
            with opener.open(Request(url+path,data=None if payload is None else json.dumps(payload).encode(),headers=headers),timeout=timeout) as response:
                raw=response.read(256*1024*1024+1)
                if len(raw)>256*1024*1024:raise ValueError('InterProScan artifacts exceed transfer limit')
                value=json.loads(raw)
                if not isinstance(value,dict):raise ValueError('invalid worker response')
                return value
        except HTTPError as exc:
            if exc.code<500:raise RuntimeError(f'InterProScan worker rejected request (HTTP {exc.code})') from exc
            if attempt==attempts-1:raise RuntimeError(f'InterProScan worker unavailable (HTTP {exc.code})') from exc
        except OSError as exc:
            if attempt==attempts-1:raise RuntimeError('InterProScan worker is unreachable') from exc
        time.sleep(.5)
    raise AssertionError('unreachable')

def stop_worker(job):
    for _ in range(6):
        try:
            snapshot=request(f'/v1/jobs/{job}/cancel',{},timeout=3,attempts=1)
            if snapshot.get('job_id')==job and snapshot.get('state')=='stopped' and snapshot.get('active') is False:
                return True
        except Exception:
            pass
        time.sleep(.5)
    return False

def annotate_full_scan(req, records, source, *, root, emit=None):
    from .domain_evidence import _record_evidence, sequence_md5, now, read_domain_analysis
    from .access_control import grant_analysis
    from .progress import emit as agent_progress
    started=time.perf_counter();check_cancelled()
    health=request('/health');release=health.get('release')
    if health.get('status')!='ready' or health.get('protocol')!=PROTOCOL or not isinstance(release,dict) or digest(release)!=health.get('release_digest'):
        raise ValueError('InterProScan worker did not identify a ready, pinned execution release')
    spec={'workflow_version':'interpro-fullscan-v1','execution_method':'full_scan','release_info':release,
          'records':[{'id':r['id'],'sequence_sha256':sha(r['sequence'].encode())} for r in records],
          'source':source,'expected_pfam':sorted(set(req.expected_pfam))}
    analysis=digest(spec);directory=root/'results/protein-domains'/analysis;directory.mkdir(parents=True,exist_ok=True)
    locks=directory.parent/'.locks';locks.mkdir(exist_ok=True)
    control=current_control();job=control.job_id if control else uuid.uuid4().hex
    base=f'/v1/jobs/{job}/analyses/{analysis}'
    def progress(status,detail=None):
        value={'name':'InterProScan 完整扫描','stage_id':'interpro_scan:'+analysis,'status':status,**(detail or {})}
        (emit or agent_progress)(value)
    with (locks/analysis).open('a+') as lock:
        while True:
            check_cancelled()
            try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);break
            except BlockingIOError:time.sleep(.1)
        manifest=directory/'manifest.json'
        try:state=json.loads(manifest.read_text()) if manifest.exists() else {}
        except (OSError,ValueError):state={}
        if state.get('spec')!=spec:state={'analysis_id':analysis,'spec':spec,'chunks':{},'attempts':[]}
        attempt={'started_at':now(),'status':'running','job_id':job,'refresh':req.refresh};state['attempts'].append(attempt);state['status']='running';atomic_json(manifest,state)
        submitted=False
        try:
            atomic_json(directory/'input-records.json',{'records':records,'source':source})
            atomic_json(directory/'worker-info.json',health)
            check_cancelled();submitted=True
            request('/v1/scan',{'job_id':job,'analysis_id':analysis,'records':[{'id':r['id'],'sequence':r['sequence']} for r in records],
                              'release_digest':health['release_digest'],'refresh':req.refresh})
            last=None
            while True:
                check_cancelled();snapshot=request(base)
                if snapshot.get('job_id')!=job or snapshot.get('analysis_id')!=analysis:raise ValueError('worker job identity mismatch')
                current=snapshot.get('state');detail=snapshot.get('progress',{})
                marker=json.dumps([current,detail],sort_keys=True)
                if marker!=last:
                    progress('running',{'phase':current,**detail});last=marker
                    attempt['worker_state']=snapshot;atomic_json(manifest,state)
                if current=='succeeded':break
                if current=='stopped':
                    if control:control.cancel()
                    raise RunCancelled()
                if current=='failed':raise RuntimeError(snapshot.get('reason','InterProScan failed'))
                if current not in {'queued','running'}:raise ValueError('worker returned unknown execution state')
                deadline=time.monotonic()+2
                while time.monotonic()<deadline:check_cancelled();time.sleep(.1)
            bundle=request(base+'/artifacts',timeout=120);check_cancelled()
            receipt=bundle.get('receipt',{});files=bundle.get('files',{})
            if bundle.get('job_id')!=job or bundle.get('analysis_id')!=analysis or receipt.get('analysis_id')!=analysis or receipt.get('release')!=release or receipt.get('exit_code')!=0 or receipt.get('de_novo_scan_executed') is not True:
                raise ValueError('scan artifact receipt identity mismatch')
            if set(files)!=FILES or set(receipt.get('files',{}))!=FILES:raise ValueError('missing full scan execution artifacts')
            decoded={}
            for name,encoded in files.items():
                raw=base64.b64decode(encoded,validate=True)
                if sha(raw)!=receipt['files'][name]:raise ValueError('scan artifact checksum mismatch')
                decoded[name]=raw
            by_md5=validate_scan(json.loads(decoded['scan.json']),records,release)
            command=json.loads(decoded['command.json'])
            if command.get('matches_api') is not False or command.get('scope')!='official_default_available_analyses' or '--no-matches-api' not in command.get('argv',[]):
                raise ValueError('worker did not certify de novo default-scope execution')
            evidence=[]
            for record in records:
                item=by_md5[sequence_md5(record)]
                result=_record_evidence(record,{**item,'found':True},spec['expected_pfam'])
                result['lookup_status']='full_scan_matches_available' if result['match_count'] else 'full_scan_without_matches'
                result['expected_pfam']={k:'not_observed_in_full_scan' if v=='not_observed_in_precomputed_record' else v for k,v in result['expected_pfam'].items()}
                result['full_scan_required_for_new_sequence']=False;evidence.append(result)
            for name,data in decoded.items():
                temporary=directory/(name+'.incoming');temporary.write_bytes(data);temporary.replace(directory/name)
            atomic_json(directory/'scan-receipt.json',receipt)
            report={'status':'succeeded','kind':'protein_domain_evidence','analysis_id':analysis,'workflow_version':spec['workflow_version'],
                    'method':'InterProScan 6 full scan, official default available analyses','de_novo_scan_executed':True,
                    'scan_executed_this_attempt':not bool(snapshot.get('cache_hit')),
                    'source':source,'release':release['interpro_version'],'library_versions':active_library_versions(release),
                    'analyses':release['analyses'],'record_count':len(records),'unique_sequence_count':len(by_md5),'unknown_record_count':0,
                    'match_count':sum(r['match_count'] for r in evidence),'records':evidence,
                    'quality_gate':{'decision':'full_scan_sequence_evidence_available','activity_confirmed':False,'catalytic_sites_verified':False},
                    'limits':['Domain matches support family evidence, not enzyme activity or substrate specificity.',
                              'No observed match does not establish biochemical absence.',release['excluded_optional']],
                    'compute':{'elapsed_seconds':receipt['elapsed_seconds'],'cache_hit':bool(snapshot.get('cache_hit')),'progress':receipt['progress']}}
            atomic_json(directory/'report.json',report)
            names=FILES|{'input-records.json','worker-info.json','scan-receipt.json','report.json'}
            state['artifacts']={n:file_sha(directory/n) for n in sorted(names)};state['status']=attempt['status']='succeeded'
            attempt['cache_hit']=bool(snapshot.get('cache_hit'));progress('succeeded',receipt['progress'])
        except BaseException as exc:
            state['status']=attempt['status']='stopped' if isinstance(exc,RunCancelled) else 'failed';attempt['error']=f'{type(exc).__name__}: {exc}'
            if submitted and not stop_worker(job):
                attempt['remote_stop_unconfirmed']=True
                raise RemoteScanStopPending('等待 InterProScan 节点确认停止') from exc
            raise
        finally:
            attempt['finished_at']=now();attempt['elapsed_seconds']=round(time.perf_counter()-started,4);atomic_json(manifest,state)
    grant_analysis('protein-domains',analysis)
    return read_domain_analysis(analysis,root=root)
