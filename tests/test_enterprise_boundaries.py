"""Account/session, artifact and operation recovery acceptance regressions."""
import json
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from fastapi import FastAPI, HTTPException
from fastapi.testclient import TestClient
from tests.test_vina_conversation import api
from marine_agent import access_control
from marine_agent.task_lifecycle import ConversationTasks

class AccountTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name);self.store = api.RunStore(self.root/'state.sqlite3')
        self.addCleanup(patch.stopall)
        patch.object(api,'store',self.store).start();patch.object(api,'ROOT',self.root).start();patch.object(api,'RESULTS',self.root/'results').start();patch.object(api,'_submit').start()
        def verify(cookie):
            if cookie not in ('session=a','session=b'):
                raise HTTPException(401,'login required')
            return cookie[-1]
        patch.object(access_control,'verify_session',side_effect=verify).start()
        app=FastAPI()
        app.router.routes = [route for route in api.app.routes if not route.path.startswith('/vina/structures/') and not route.path.startswith('/vina/receptors')]
        app.add_middleware(access_control.AccountMiddleware)
        from marine_agent.receptor_api import register_receptor_routes
        register_receptor_routes(app,self.store,api.executor,api._find_viewer_manifest)
        self.client=TestClient(app)
    def call(self,method,path,owner=None,**kw):
        return self.client.request(method,path,headers={} if owner is None else {'cookie':'session='+owner},**kw)
    def create(self,owner):
        response=self.call('POST','/runs',owner,json={'prompt':'你好'})
        self.assertEqual(response.status_code,202,response.text)
        return response.json()
    def test_cookie_identity_and_all_run_surfaces(self):
        a=self.create('a');b=self.create('b');rid=a['run_id']
        self.assertEqual(self.call('GET','/runs').status_code,401)
        self.assertEqual(self.client.get('/runs',headers={'X-User-Id':'a'}).status_code,401)
        self.assertEqual([r['id'] for r in self.call('GET','/runs','b').json()['runs']],[b['run_id']])
        self.assertEqual(len(self.call('GET','/conversations','b').json()['conversations']),1)
        for path in [f'/runs/{rid}',f'/runs/{rid}/report',f'/runs/{rid}/events',f'/vina/runs/{rid}',f'/vina/structures/{rid}',f'/vina/receptors/{rid}',f'/conversations/{a["conversation_id"]}']:
            self.assertEqual(self.call('GET',path,'b').status_code,404,path)
        for path in [f'/runs/{rid}/cancel',f'/runs/{rid}/continue']:
            self.assertEqual(self.call('POST',path,'b',json={}).status_code,404,path)
        self.assertEqual(self.call('POST','/runs','b',json={'prompt':'你好','conversation_id':a['conversation_id']}).status_code,404)
        self.assertEqual(self.call('GET',f'/runs/{rid}','a').status_code,200)
    def test_owned_run_cannot_read_other_runs_artifacts_or_database(self):
        a=self.create('a');b=self.create('b')
        own=self.root/'results/vina-runs'/a['run_id']/'pose.pdbqt';own.parent.mkdir(parents=True);own.write_text('own pose')
        foreign=self.root/'results/vina-runs'/b['run_id']/'pose.pdbqt';foreign.parent.mkdir(parents=True);foreign.write_text('private pose')
        for path,expected in [(own,200),(foreign,404),(self.root/'state.sqlite3',404)]:
            relative=str(path.relative_to(self.root));r=self.call('GET',f'/runs/{a["run_id"]}/artifact/{relative}','a')
            self.assertEqual(r.status_code,expected)
        fake={'agent':{'messages':[{'type':'human','content':json.dumps({'viewer_manifest':str(foreign.relative_to(self.root))})}]}}
        self.store.update(a['run_id'],result_json=json.dumps(fake))
        self.assertFalse(api._owns_artifact(self.store.get(a['run_id']),foreign))
    def test_analysis_ownership_survives_restart_and_context_does_not_bleed(self):
        with access_control.account_scope('a'):
            access_control.grant_analysis('protein-domains','a'*64)
            access_control.require_analysis('protein-domains','a'*64)
        with access_control.account_scope('b'), self.assertRaises(ValueError):
            access_control.require_analysis('protein-domains','a'*64)
        reopened=api.RunStore(self.root/'state.sqlite3')
        self.assertTrue(reopened.owns_analysis('protein-domains','a'*64,'a'))
        self.assertIsNone(access_control.principal.get())
    def test_internal_callbacks_fail_closed_without_configured_secret(self):
        a=self.create('a')
        with patch.object(api,'VINA_CALLBACK_TOKEN',''):
            self.assertEqual(self.call('POST',f'/internal/vina-progress/{a["run_id"]}',json={}).status_code,401)

class SelectiveContinuationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=api.RunStore(Path(self.temp.name)/'state.sqlite3')
        self.first=self.store.create('two independent queries','general')
        self.spec={'workflow':'protein_domains','input':{'query_fasta':'>x\n'+'A'*40,'source':'test'}}
        self.store.update(self.first['run_id'],execution_json=json.dumps({'workflow':'multiple'}),status='partial')
        for n,status in [('first','failed'),('second','succeeded')]:
            self.store.save_operation(self.first['run_id'],n,'annotate_protein_domains',self.spec)
            self.store.complete_operation(self.first['run_id'],n,status,{})
    def test_selected_failed_operation_preserves_input_without_replaying_sibling(self):
        current=self.store.create('continue first','general',conversation_id=self.first['conversation_id'])
        service=ConversationTasks(self.store,current['run_id']);r=service.prepare(self.first['run_id'],operation_index=1)
        self.assertEqual(r['status'],'handoff_ready');self.assertEqual(service.handoff['spec'],self.spec)
        self.assertEqual(self.store.get(current['run_id'])['continuation_operation_id'],'first')
        self.assertEqual(self.store.get(self.first['run_id'])['operations'][1]['status'],'succeeded')
    def test_full_scan_unavailable_is_not_a_cached_lookup_success(self):
        current=self.store.create('full scan first','general',conversation_id=self.first['conversation_id'])
        service=ConversationTasks(self.store,current['run_id'])
        self.assertEqual(service.prepare(self.first['run_id'],operation_index=1,method='full_scan')['status'],'capability_unavailable')
        self.assertIsNone(service.handoff)
    def test_completed_operation_cannot_be_selected_for_retry(self):
        current=self.store.create('continue second','general',conversation_id=self.first['conversation_id'])
        self.assertEqual(ConversationTasks(self.store,current['run_id']).prepare(self.first['run_id'],operation_index=2)['status'],'already_completed')

    def test_running_turn_can_inspect_its_persisted_operations(self):
        current=self.store.create('new analysis','general',conversation_id=self.first['conversation_id'])
        self.store.save_operation(current['run_id'],'current','annotate_protein_domains',self.spec)
        self.store.complete_operation(current['run_id'],'current','partial',{'analysis_id':'a'*64})
        service=ConversationTasks(self.store,current['run_id'])
        visible=service.snapshot(current['run_id'])['tasks']
        self.assertEqual(visible[0]['operations'][0]['status'],'partial')
        identity_tool=next(t for t in service.tools() if t.name=='inspect_saved_input_identity')
        with patch('marine_agent.input_identity.inspect_input_identity',return_value={'status':'partial'}) as inspect:
            identity_tool.invoke({'source_run_id':current['run_id'],'sequence_index':1})
            self.assertEqual(inspect.call_args.args,([self.spec],1))

class InputIdentityTests(unittest.TestCase):
    def test_retry_specs_are_deduplicated_before_input_selection(self):
        from marine_agent.input_identity import inspect_input_identity
        first={'input':{'query_fasta':'>P00782_variant\nACDE'}}
        second={'input':{'query_fasta':'>unresolved_second\nAAAA'}}
        with patch('marine_agent.input_identity.get_uniprot_sequence') as fetch:
            fetch.invoke.return_value=json.dumps({'status':'succeeded','sequence':'LCDE','accession':'P00782','metadata':{'protein_name':'Subtilisin'}})
            result=inspect_input_identity([first,second,first],1)
            self.assertEqual(len(result['records']),1)
            self.assertEqual(result['records'][0]['substitutions'],['L1A'])
            all_records=inspect_input_identity([first,second,first])['records']
            self.assertEqual(len(all_records),2)

class AuthMigrationTests(unittest.TestCase):
    def test_malformed_identity_service_response_is_unavailable(self):
        with patch.dict(os.environ,{'MARINE_AUTH_URL':'http://auth.invalid/me'}),patch('urllib.request.build_opener') as opener:
            opener.return_value.open.return_value=io.BytesIO(b'[]')
            with self.assertRaises(HTTPException) as raised:
                access_control.verify_session('test-session')
            self.assertEqual(raised.exception.status_code,503)

    def test_legacy_analysis_grants_require_real_tool_provenance(self):
        with tempfile.TemporaryDirectory() as directory,patch.dict(os.environ,{'MARINE_LEGACY_OWNER_ID':''}):
            path=Path(directory)/'state.sqlite3';store=api.RunStore(path)
            row=store.create('historic input','general')
            messages=[{'type':role,'name':'annotate_protein_domains','content':json.dumps({'analysis_id':char*64})} for role,char in [('tool','a'),('human','b'),('ai','c')]]
            store.update(row['run_id'],result_json=json.dumps({'agent':{'messages':messages}}))
            with patch.dict(os.environ,{'MARINE_LEGACY_OWNER_ID':'admin'}):
                migrated=api.RunStore(path)
            self.assertTrue(migrated.owns_analysis('protein-domains','a'*64,'admin'))
            for char in ['b','c']:
                self.assertFalse(migrated.owns_analysis('protein-domains',char*64,'admin'))

if __name__=='__main__':unittest.main()
