"""Contract tests; real scientific acceptance is recorded separately in results/."""
import base64
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import threading
import time
import unittest
import uuid
from unittest.mock import patch
from marine_agent import domain_evidence, interpro_remote
from marine_agent.interpro_protocol import PROTOCOL, digest, validate_scan
from marine_agent.interpro_worker import Worker, atomic
from marine_agent.execution_control import job_scope, RunCancelled
from marine_agent.evidence_review import render_domains, explain_activity_evidence
from marine_agent.restart_recovery import remote_services

RELEASE={'interpro_version':'110.0','interproscan_version':'6.0.3.0','analyses':[{'name':'Pfam'}],
         'library_versions':{'Pfam':'38.2','Phobius':'1.01'},'excluded_optional':'Optional licensed/ML analyses not configured.'}
RECORDS=[{'id':'new','sequence':'MALWMRLLPLLALLALWGPGPG'}, {'id':'duplicate','sequence':'MALWMRLLPLLALLALWGPGPG'}, {'id':'nohit','sequence':'ACDEFGHIKLMNPQRSTVWY'}]

def scan():
    return {'interpro-version':'110.0','interproscan-version':'6.0.3.0','results':[
        {'sequence':RECORDS[0]['sequence'],'md5':hashlib.md5(RECORDS[0]['sequence'].encode()).hexdigest().upper(),
         'xref':[{'id':'new'},{'id':'duplicate'}], 'matches':[{'signature':{'accession':'PF00082','description':'Subtilase family','signatureLibraryRelease':{'library':'Pfam','version':'38.2'}},'locations':[{'start':2,'end':10,'evalue':1e-12}]}]},
        {'sequence':RECORDS[2]['sequence'],'md5':hashlib.md5(RECORDS[2]['sequence'].encode()).hexdigest().upper(),'xref':[{'id':'nohit'}],'matches':[]}]}

def payload():return {'query_fasta':''.join(f">{r['id']}\n{r['sequence']}\n" for r in RECORDS),'source':'synthetic contract fixtures','method':'full_scan','expected_pfam':['PF00082']}

class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.calls=[];self.job=None;self.analysis=None;self.broken=None
    def remote(self,path,data=None,**kwargs):
        self.calls.append((path,data))
        if path=='/health':return {'protocol':PROTOCOL,'status':'ready','release':RELEASE,'release_digest':digest(RELEASE)}
        if path=='/v1/scan':self.job=data['job_id'];self.analysis=data['analysis_id'];return {'state':'queued'}
        if path.endswith('/cancel'):return {'job_id':self.job,'state':'stopped','active':False}
        if not path.endswith('/artifacts'):return {'job_id':self.job,'analysis_id':self.analysis,'state':'succeeded','cache_hit':False}
        raw=scan()
        if self.broken=='missing':raw['results'].pop()
        if self.broken=='coordinate':raw['results'][0]['matches'][0]['locations'][0]['end']=10000
        files={'scan.json':json.dumps(raw).encode(),'scan.tsv':b'fake\tfixture\n','trace.tsv':b'status\tname\nCOMPLETED\tSCAN\n',
               'command.json':json.dumps({'argv':['nextflow','--no-matches-api'],'matches_api':False,'scope':'official_default_available_analyses'}).encode(),
               'console.log':b'fixture log','nextflow.log':b'fixture log'}
        result={'job_id':self.job,'analysis_id':self.analysis,'receipt':{'analysis_id':self.analysis,'release':RELEASE,'exit_code':0,'de_novo_scan_executed':True,'elapsed_seconds':5,'progress':{'completed_processes':1},'files':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()}},'files':{k:base64.b64encode(v).decode() for k,v in files.items()}}
        if self.broken=='hash':result['receipt']['files']['scan.json']='a'*64
        return result
    def run_scan(self,emit=None):
        with patch.object(interpro_remote,'request',side_effect=self.remote),patch.object(domain_evidence,'_fetch',side_effect=AssertionError('Matches API must not be used')):
            return domain_evidence.annotate(payload(),root=self.root,emit=emit)
    def test_full_scan_preserves_zero_matches_duplicates_and_non_activity_scope(self):
        report=self.run_scan()
        self.assertTrue(report['de_novo_scan_executed']);self.assertEqual(report['record_count'],3)
        self.assertEqual(report['library_versions'], {'Pfam':'38.2'})
        self.assertEqual(report['unique_sequence_count'],2);self.assertEqual(report['unknown_record_count'],0)
        self.assertEqual(report['records'][2]['expected_pfam']['PF00082'],'not_observed_in_full_scan')
        self.assertFalse(report['quality_gate']['activity_confirmed'])
        self.assertIn('完整扫描完成',render_domains(report));self.assertNotIn('没有运行新的',render_domains(report))
        self.assertIn('InterProScan 完整扫描',explain_activity_evidence([report])['final_text'])
    def test_missing_invalid_coordinates_and_bad_hash_cannot_succeed(self):
        for failure in ['missing','coordinate','hash']:
            self.broken=failure
            with self.assertRaises(ValueError):self.run_scan()
            states=[json.loads(p.read_text()) for p in self.root.glob('results/protein-domains/*/manifest.json')]
            self.assertTrue(all(s['status']=='failed' for s in states))
            self.assertTrue(any(p.endswith('/cancel') for p,_ in self.calls))
    def test_tampered_local_output_is_rejected(self):
        report=self.run_scan();path=self.root/'results/protein-domains'/report['analysis_id']/'scan.json';path.write_text('{}')
        with self.assertRaises(ValueError):domain_evidence.read_domain_analysis(report['analysis_id'],root=self.root)
    def test_cancellation_calls_worker_and_saves_interrupted_attempt(self):
        with job_scope(uuid.uuid4().hex) as control:
            def emit(event):control.cancel()
            with self.assertRaises(RunCancelled):self.run_scan(emit)
        self.assertTrue(any(p.endswith('/cancel') for p,_ in self.calls))
        manifest=next(self.root.glob('results/protein-domains/*/manifest.json'))
        self.assertEqual(json.loads(manifest.read_text())['status'],'stopped')
    def test_unconfirmed_remote_stop_is_not_reported_as_stopped(self):
        with patch.object(interpro_remote, 'stop_worker', return_value=False):
            self.broken='missing'
            with self.assertRaises(interpro_remote.RemoteScanStopPending):self.run_scan()
        state=json.loads(next(self.root.glob('results/protein-domains/*/manifest.json')).read_text())
        self.assertTrue(state['attempts'][-1]['remote_stop_unconfirmed'])

    def test_protocol_rejects_wrong_xref_sequence_version_and_duplicate_results(self):
        for bad in ['xref','sequence','version','duplicate','library']:
            raw=scan()
            if bad=='xref':raw['results'][0]['xref'][0]['id']='foreign'
            if bad=='sequence':raw['results'][0]['sequence']='AAAA'
            if bad=='version':raw['interpro-version']='111.0'
            if bad=='duplicate':raw['results'].append(copy.deepcopy(raw['results'][0]))
            if bad=='library':raw['results'][0]['matches'][0]['signature']['signatureLibraryRelease']['version']='38.1'
            with self.assertRaises(ValueError):validate_scan(raw,RECORDS,RELEASE)
    def test_restart_recovery_includes_fullscan_worker(self):
        with patch.dict(os.environ,{'INTERPRO_WORKER_URL':'http://worker','INTERPRO_WORKER_TOKEN':'test'}):
            self.assertIn(('InterPro','http://worker','test'),remote_services({'execution_spec':{'workflow':'protein_domains','input':{'method':'full_scan'}}}))

    def test_multiple_recorded_scans_do_not_wait_for_unrelated_workers(self):
        env={'INTERPRO_WORKER_URL':'http://interpro','INTERPRO_WORKER_TOKEN':'test',
             'OPENFOLD_WORKER_URL':'http://unreachable-openfold','VINA_WORKER_URL':'http://unreachable-vina'}
        op={'spec':{'workflow':'protein_domains','input':{'method':'full_scan'}}}
        with patch.dict(os.environ,env,clear=True):
            for summary in ({'workflow':'multiple'}, {}):
                services=remote_services({'execution_spec':summary,'operations':[op,op]})
                self.assertEqual(services,[('InterPro','http://interpro','test')])
            self.assertEqual({s[0] for s in remote_services({})},{'OpenFold','Vina','InterPro'})

class WorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.release=self.root/'release';self.release.mkdir();atomic(self.release/'release.json',RELEASE)
        self.worker=Worker(self.root,self.release)
        self.patch=patch.object(self.worker,'_containers',return_value=[]);self.patch.start();self.addCleanup(self.patch.stop)
    def data(self,job=None,analysis=None):return {'job_id':job or uuid.uuid4().hex,'analysis_id':analysis or 'a'*64,'release_digest':digest(RELEASE),'records':RECORDS}
    def test_cancel_before_submission_persists_across_worker_restart(self):
        data=self.data();self.worker.cancel(data['job_id'])
        restarted=Worker(self.root,self.release)
        with self.assertRaises(ValueError):restarted.submit(data)
        self.assertEqual(restarted.status(data['job_id'])['state'],'stopped')
    def test_queued_job_stop_and_idempotent_submission(self):
        self.worker.queue_lock.acquire();data=self.data()
        try:
            self.worker.submit(data);self.worker.submit(data)
            self.assertEqual(len(self.worker.threads),1)
            changed=copy.deepcopy(data);changed['records'][0]['sequence']='ACDE'
            with self.assertRaises(ValueError):self.worker.submit(changed)
            self.worker.cancel(data['job_id'])
            self.worker.threads[data['job_id'],data['analysis_id']].join(3)
            state=self.worker.status(data['job_id']);self.assertFalse(state['active']);self.assertEqual(state['state'],'stopped')
        finally:self.worker.queue_lock.release()
    def test_wrong_release_is_rejected_before_creating_a_job(self):
        data=self.data();data['release_digest']='0'*64
        with self.assertRaises(ValueError):self.worker.submit(data)
        self.assertFalse(list(self.worker.jobs.glob('*/state.json')))
    def test_restart_marks_interrupted_jobs_stopped_and_cleans_only_owned_containers(self):
        data=self.data();path=self.worker._path(data['job_id'])
        atomic(path,{'job_id':data['job_id'],'cancelled':False,'analyses':{data['analysis_id']:{'state':'running'}}})
        with patch.object(Worker,'_containers',return_value=[]) as cleanup:
            recovered=Worker(self.root,self.release)
        cleanup.assert_called_once_with(data['job_id'],kill=True)
        self.assertEqual(recovered.status(data['job_id'])['state'],'stopped')

if __name__=='__main__':unittest.main()
