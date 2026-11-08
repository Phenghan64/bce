import hashlib
import json
import tempfile
import threading
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch
from langchain.agents.middleware.types import ToolCallRequest
from langchain_core.messages import ToolMessage
from langchain_core.tools import tool
from tests.test_vina_conversation import api
from marine_agent.coordination import ExecutionReceipts, domain_completion_contract, domain_contract_progress, operation_outcome
from marine_agent.task_lifecycle import ConversationTasks

@tool
def annotate_protein_domains(query_fasta: str, source: str='test', refresh: bool=False, expected_pfam: list[str]|None=None) -> str:
    """Test fixture for exact-request execution coordination."""
    return '{}'

def request(identifier='one', **args):
    return ToolCallRequest(tool_call={'name':'annotate_protein_domains','id':identifier,'args':{'query_fasta':'>q\nACDE',**args},'type':'tool_call'},tool=annotate_protein_domains,state={},runtime=None)

def result(req,status='succeeded'):
    return ToolMessage(content=json.dumps({'status':status,'analysis_id':'a'*64,'release':'110','records':[{'protein_id':'q'}]}),tool_call_id=req.tool_call['id'],name=req.tool_call['name'])

class ReceiptTests(unittest.TestCase):
    def test_parallel_duplicate_has_one_real_execution_and_distinct_receipts(self):
        middleware=ExecutionReceipts('run');calls=[]
        def handler(req):
            calls.append(req.tool_call['id']);time.sleep(.12);return result(req,'partial')
        with ThreadPoolExecutor(max_workers=2) as pool:
            responses=list(pool.map(lambda name:middleware.wrap_tool_call(request(name),handler),['one','two']))
        self.assertEqual(len(calls),1)
        data=[json.loads(r.content)['execution_reference'] for r in responses]
        self.assertEqual(sum(r['reused'] for r in data),1)
        self.assertEqual({r.tool_call_id for r in responses},{'one','two'})
        self.assertEqual({r['run_id'] for r in data},{'run'})
    def test_changed_parameter_refresh_and_new_turn_execute(self):
        middleware=ExecutionReceipts('run');calls=[]
        def handler(req):calls.append(req.tool_call['args']);return result(req)
        middleware.wrap_tool_call(request('one'),handler)
        middleware.wrap_tool_call(request('same',source='test',refresh=False,expected_pfam=None),handler)
        self.assertEqual(len(calls),1)
        middleware.wrap_tool_call(request('changed',expected_pfam=['PF00082']),handler)
        middleware.wrap_tool_call(request('refresh',refresh=True),handler)
        middleware.wrap_tool_call(request('refresh2',refresh=True),handler)
        ExecutionReceipts('new-run').wrap_tool_call(request('new'),handler)
        self.assertEqual(len(calls),4)
    def test_failed_result_can_retry_and_explicit_recompute_sets_refresh(self):
        middleware=ExecutionReceipts('run');calls=[]
        def handler(req):calls.append(req);return result(req,'failed' if len(calls)==1 else 'succeeded')
        middleware.wrap_tool_call(request('failed'),handler);middleware.wrap_tool_call(request('retry'),handler)
        self.assertEqual(len(calls),2)
        explicit=ExecutionReceipts('new',recompute=True)
        explicit.wrap_tool_call(request('fresh'),handler);explicit.wrap_tool_call(request('fresh2'),handler)
        self.assertTrue(calls[-1].tool_call['args']['refresh']);self.assertEqual(len(calls),3)
        repetitions=ExecutionReceipts('independent',allow_repeats=True)
        repetitions.wrap_tool_call(request('trial1'),handler);repetitions.wrap_tool_call(request('trial2'),handler)
        self.assertEqual(len(calls),5)

class CompletionTests(unittest.TestCase):
    def setUp(self):
        self.prompt='只查结构域\n>a\nACDE\n>b\nACDF\n'
        self.contract=domain_completion_contract(self.prompt,'provided_fasta_domain_lookup')
    def state(self, records):
        return {'messages':[{'type':'tool','name':'annotate_protein_domains','content':json.dumps({'status':'partial','analysis_id':'a'*64,'release':'110','records':records})}]}
    def row(self,name,seq):
        from marine_agent.domain_evidence import _record_evidence
        return _record_evidence({'id':name,'sequence':seq},None,[])
    def test_partial_unknown_evidence_finishes_lookup_but_not_scientific_validation(self):
        progress=domain_contract_progress(self.contract,self.state([self.row('a','ACDE'),self.row('b','ACDF')]))
        self.assertTrue(progress['complete']);self.assertEqual(progress['covered_inputs'],2)
    def test_missing_or_wrong_sequence_cannot_finish(self):
        progress=domain_contract_progress(self.contract,self.state([self.row('a','ACDE'),self.row('b','AAAA')]))
        self.assertFalse(progress['complete']);self.assertEqual(progress['missing_inputs'],['b'])
    def test_open_task_unparseable_fasta_and_recompute_do_not_get_early_stop(self):
        self.assertIsNone(domain_completion_contract(self.prompt,'open'))
        self.assertIsNone(domain_completion_contract(self.prompt+'然后再预测','provided_fasta_domain_lookup'))
        self.assertIsNone(domain_completion_contract(self.prompt,'provided_fasta_domain_lookup',True))
    def test_explicit_pfam_check_must_be_present_in_evidence(self):
        contract=domain_completion_contract('核验PF00082\n>a\nACDE','provided_fasta_domain_lookup')
        row=self.row('a','ACDE')
        self.assertFalse(domain_contract_progress(contract,self.state([row]))['complete'])
        row['expected_pfam']={'PF00082':'unknown'}
        self.assertTrue(domain_contract_progress(contract,self.state([row]))['complete'])
    def test_successful_sibling_does_not_mask_failure_but_exact_retry_can_succeed(self):
        a={'tool':'annotate_protein_domains','spec':{'input':'a'},'status':'failed'}
        b={'tool':'annotate_protein_domains','spec':{'input':'b'},'status':'succeeded'}
        self.assertEqual(operation_outcome({'operations':[a,b]}),'partial')
        self.assertEqual(operation_outcome({'operations':[a,b,{**a,'status':'succeeded'}]}),'succeeded')

class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=api.RunStore(Path(self.temp.name)/'state.sqlite3')
        self.old=self.store.create('two inputs','general');self.store.update(self.old['run_id'],status='partial')
        for n,seq,status in [('first','ACDE','partial'),('second','ACDF','succeeded')]:
            spec={'workflow':'protein_domains','input':{'query_fasta':f'>{n}\n{seq}','source':'test'}}
            self.store.save_operation(self.old['run_id'],n,'annotate_protein_domains',spec)
            self.store.complete_operation(self.old['run_id'],n,status,{'analysis_id':('a' if n=='first' else 'b')*64,'unknown_record_count':int(n=='first')})
        self.current=self.store.create('followup','general',conversation_id=self.old['conversation_id'])
        self.service=ConversationTasks(self.store,self.current['run_id'])
    def test_analysis_id_resolves_to_its_own_operation_and_full_scan_stays_unavailable(self):
        ref=self.service.resolve_reference('a'*64)
        self.assertEqual(ref['operation']['operation_id'],'first')
        self.assertEqual(self.service.prepare(operation_id='a'*64,method='full_scan')['status'],'capability_unavailable')
        self.assertIsNone(self.service.handoff)
    def test_configured_full_scan_handoff_keeps_only_selected_saved_input(self):
        with patch('marine_agent.interpro_remote.configured', return_value=True):
            result=self.service.prepare(operation_id='a'*64,method='full_scan')
        self.assertEqual(result['status'],'handoff_ready')
        self.assertEqual(self.service.handoff['spec']['input']['method'],'full_scan')
        self.assertEqual(self.service.handoff['spec']['input']['query_fasta'],'>first\nACDE')
        self.assertEqual(self.service.handoff['operation_id'],'first')
        original=self.store.get(self.old['run_id'])
        self.assertEqual(original['operations'][1]['status'],'succeeded')
        self.assertNotIn('method',original['operations'][0]['spec']['input'])
    def test_other_conversation_and_ambiguous_analysis_do_not_resolve(self):
        other=self.store.create('foreign','general')
        self.assertEqual(self.service.resolve_reference(other['run_id'])['status'],'not_found')
        self.store.complete_operation(self.old['run_id'],'second','succeeded',{'analysis_id':'a'*64})
        self.assertEqual(self.service.resolve_reference('a'*64)['status'],'ambiguous_reference')
    def test_identity_from_analysis_reads_only_that_input(self):
        tool=next(t for t in self.service.tools() if t.name=='inspect_saved_input_identity')
        with patch('marine_agent.input_identity.inspect_input_identity',return_value={'status':'partial'}) as inspect:
            tool.invoke({'source_run_id':'b'*64})
            self.assertEqual(inspect.call_args.args[0][0]['input']['query_fasta'],'>second\nACDF')
            self.assertEqual(len(inspect.call_args.args[0]),1)
    def test_fresh_turn_has_no_continuation_capability(self):
        fresh=self.store.create('new standalone query','general');service=ConversationTasks(self.store,fresh['run_id'])
        self.assertFalse(service.has_resumable_work())
        self.assertNotIn('continue_task',[t.name for t in service.tools(allow_continue=service.has_resumable_work())])


class ConversationReviewTests(unittest.TestCase):
    def test_incidental_identity_read_does_not_replace_the_requested_explanation(self):
        from marine_agent.evidence_review import review_answer
        state={'messages':[{'type':'tool','name':'inspect_saved_input_identity','content':json.dumps({'records':[{'label':'input'}],'status':'succeeded','final_text':'input identity only'})}]}
        answer='不能证明胶原酶活性；结构域证据不足，还需要胶原底物实验。'
        self.assertEqual(review_answer(state,answer,scope='conversation')['final_text'],answer)
        self.assertEqual(review_answer(state,answer)['final_text'],'input identity only')

class ActivityEvidenceTests(unittest.TestCase):
    def test_first_unknown_input_cannot_borrow_sibling_family_or_generated_annotations(self):
        from marine_agent.evidence_review import explain_activity_evidence
        reports=[{'analysis_id':'a'*64,'release':'110','status':'partial','records':[
            {'protein_id':'first','full_scan_required_for_new_sequence':True,'pfam':[]},
            {'protein_id':'second','full_scan_required_for_new_sequence':False,'pfam':[{'accession':'PF00318','description':'Ribosomal protein S2','locations':[{'start':9,'end':224}]}]}]}]
        evidence=explain_activity_evidence(reports,1)
        self.assertFalse(evidence['activity_established']);self.assertFalse(evidence['new_computation'])
        self.assertEqual(evidence['selected_protein_ids'],['first'])
        self.assertNotIn('PF00318',evidence['final_text'])
        self.assertNotIn('PF00144',evidence['final_text'])
        self.assertIn('未知',evidence['final_text'])
    def test_unavailable_records_are_not_reported_as_completed_assessment(self):
        from marine_agent.evidence_review import explain_activity_evidence
        self.assertEqual(explain_activity_evidence([{'status':'unavailable'}])['status'],'partial')



class MixedEvidenceTests(unittest.TestCase):
    def test_incidental_identity_cannot_erase_requested_domain_report(self):
        from marine_agent.evidence_review import review_answer
        from marine_agent.domain_evidence import _record_evidence
        domain={'kind':'protein_domain_evidence','status':'partial','analysis_id':'a'*64,'release':'110','record_count':1,'unknown_record_count':1,'returned_records':1,'records':[_record_evidence({'id':'mutant','sequence':'ACDE'},None,[])]}
        state={'messages':[{'type':'tool','name':'annotate_protein_domains','content':json.dumps(domain)}, {'type':'tool','name':'inspect_saved_input_identity','content':json.dumps({'records':[{'label':'mutant'}],'status':'succeeded','final_text':'identity-only'})}]}
        reviewed=review_answer(state,'draft')
        self.assertEqual(reviewed['decision'],'rendered_from_domain_evidence')
        self.assertIn('结构域证据核验',reviewed['final_text'])
        self.assertEqual(reviewed['run_status'],'partial')



class FreshEvidenceTests(unittest.TestCase):
    def test_refresh_report_uses_new_computation_not_an_incidental_old_read(self):
        from marine_agent.evidence_review import review_answer
        from marine_agent.domain_evidence import _record_evidence
        report={'kind':'protein_domain_evidence','status':'partial','release':'110','record_count':1,'unknown_record_count':1,'returned_records':1,'records':[_record_evidence({'id':'mutant','sequence':'ACDE'},None,[])]}
        state={'messages':[{'type':'tool','name':name,'content':json.dumps({**report,'analysis_id':identifier})} for name,identifier in [('annotate_protein_domains','new-result'),('inspect_domain_evidence','old-result')]]}
        answer=review_answer(state,'draft')['final_text']
        self.assertIn('new-result',answer);self.assertNotIn('old-result',answer)



class StoreResourceTests(unittest.TestCase):
    def test_commit_rollback_and_deterministic_connection_close(self):
        import sqlite3
        with tempfile.TemporaryDirectory() as directory:
            store=api.RunStore(Path(directory)/'state.sqlite3')
            with store._connect() as db:db.execute("INSERT INTO migrations VALUES ('committed')")
            with self.assertRaises(sqlite3.ProgrammingError):db.execute('SELECT 1')
            with self.assertRaises(RuntimeError):
                with store._connect() as db:
                    db.execute("INSERT INTO migrations VALUES ('rolled-back')")
                    raise RuntimeError('rollback test')
            with store._connect() as db:
                self.assertIsNotNone(db.execute("SELECT 1 FROM migrations WHERE name='committed'").fetchone())
                self.assertIsNone(db.execute("SELECT 1 FROM migrations WHERE name='rolled-back'").fetchone())
            held=[]
            for _ in range(400):
                with store._connect() as db:db.execute('SELECT 1').fetchone();held.append(db)
            self.assertEqual(len(held),400)
            for db in held:
                with self.assertRaises(sqlite3.ProgrammingError):db.execute('SELECT 1')



class FailedRefreshTests(unittest.TestCase):
    def test_failed_refresh_cannot_be_replaced_by_an_old_successful_read(self):
        from marine_agent.evidence_review import review_answer
        old={'kind':'protein_domain_evidence','status':'succeeded','analysis_id':'old-result','release':'110','record_count':0,'unknown_record_count':0,'returned_records':0,'records':[]}
        state={'messages':[{'type':'tool','name':'annotate_protein_domains','content':json.dumps({'status':'failed','reason':'upstream unavailable'})},{'type':'tool','name':'inspect_domain_evidence','content':json.dumps(old)}]}
        answer=review_answer(state,'success')
        self.assertEqual(answer['run_status'],'failed');self.assertNotIn('old-result',answer['final_text'])

if __name__=='__main__':unittest.main()
