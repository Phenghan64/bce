"""Exact sequence evidence, missing data semantics, restart and report grounding."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import uuid
from unittest.mock import patch

from marine_agent import domain_evidence as domains
from marine_agent.evidence_review import review_answer, render_domains
from marine_agent.execution_control import job_scope, RunCancelled
from marine_agent.protein_screening import parse_fasta, run_screen
from tests.test_vina_conversation import api
from marine_agent.task_lifecycle import saved_spec

FIXTURES = Path(__file__).parent / 'fixtures'
INFO = json.loads((FIXTURES / 'domain-evidence/info.json').read_text())
MATCHES = json.loads((FIXTURES / 'domain-evidence/matches.json').read_text())


def payload():
    return {'query_fasta': ''.join((FIXTURES / f'protein-screening/{a}.faa').read_text() for a in ('P00782','P0A7V0','P69905')),
            'source': 'UniProt fixtures, 2026-11-07', 'expected_pfam': ['PF00082', 'PF00040']}


def fetch(url, body=None):
    value = INFO if body is None else {'results': [r for r in MATCHES['results'] if r['md5'] in body['md5']]}
    return copy.deepcopy(value), {'url': url, 'retrieved_at': '2026-11-07T12:00:00Z'}


def tool_state(report, name='annotate_protein_domains'):
    return {'messages':[{'type':'tool','name':name,'content':json.dumps(report)}]}


class DomainTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.mock = patch.object(domains, '_fetch', side_effect=fetch)
        self.mock.start()
        self.addCleanup(self.mock.stop)

    def run_lookup(self, data=None, **kwargs):
        return domains.annotate(data or payload(), root=self.root, **kwargs)

    def test_real_api_fixture_has_correct_names_coordinates_versions_and_no_negative_activity(self):
        report = self.run_lookup()
        self.assertEqual(report['status'], 'succeeded')
        self.assertEqual(report['record_count'], 3)
        self.assertFalse(report['de_novo_scan_executed'])
        self.assertFalse(report['quality_gate']['activity_confirmed'])
        first = report['records'][0]
        hit = next(p for p in first['pfam'] if p['accession']=='PF00082')
        self.assertEqual(hit['description'], 'Subtilase family')
        self.assertEqual(hit['library_version'], '38.2')
        self.assertEqual(hit['locations'][0]['start'], 130)
        self.assertTrue(all(p['analysis_source']=='Pfam' for p in first['pfam']))
        self.assertTrue(any(p['analysis_source']=='InterPro-N' for p in first['predicted_pfam']))
        self.assertEqual(first['expected_pfam']['PF00082'], 'observed')
        self.assertEqual(first['expected_pfam']['PF00040'], 'not_observed_in_precomputed_record')
        self.assertNotIn('PF00082', [p['accession'] for r in report['records'][1:] for p in r['pfam']])

    def test_identical_sequences_preserve_all_record_ids_but_lookup_unique_md5_once(self):
        data=payload();seq=parse_fasta(data['query_fasta'])[0]['sequence']
        data['query_fasta']=f'>one\n{seq}\n>two\n{seq}\n'
        report=self.run_lookup(data)
        self.assertEqual(report['record_count'],2)
        self.assertEqual(report['unique_sequence_count'],1)
        self.assertEqual([r['protein_id'] for r in report['records']],['one','two'])
        self.assertEqual(report['records'][0]['pfam'],report['records'][1]['pfam'])

    def test_omitted_and_explicit_unknown_md5_are_unknown_not_no_domain(self):
        for mode in ('omitted','unknown'):
            data=payload();data['query_fasta']='>new\nMALWMRLLPLLALLALWGPGPG\n';data['source']=mode
            digest=domains.sequence_md5(parse_fasta(data['query_fasta'])[0])
            def missing(url,body=None):
                return (INFO if body is None else {'results':[] if mode=='omitted' else [{'md5':digest,'found':False,'matches':[]}]}),{}
            with patch.object(domains,'_fetch',side_effect=missing):report=self.run_lookup(data)
            self.assertEqual(report['status'],'partial')
            self.assertEqual(report['records'][0]['expected_pfam']['PF00082'],'unknown')
            self.assertTrue(report['records'][0]['full_scan_required_for_new_sequence'])
            self.assertIn('不能判断结构域缺失',render_domains(report))

    def test_same_release_cache_integrity_refresh_and_release_invalidation(self):
        report=self.run_lookup()
        with patch.object(domains,'_fetch',side_effect=fetch) as calls:cached=self.run_lookup()
        self.assertEqual(calls.call_count,1)  # release check only
        self.assertEqual(cached['execution']['reused_chunks'],['0'])
        chunk=self.root/'results/protein-domains'/report['analysis_id']/'chunk-0.json'
        chunk.write_text('{}')
        with self.assertRaises(ValueError):domains.read_domain_analysis(report['analysis_id'],root=self.root)
        repaired=self.run_lookup();self.assertEqual(repaired['execution']['reused_chunks'],[])
        data=payload();data['refresh']=True
        refreshed=self.run_lookup(data);self.assertEqual(refreshed['execution']['reused_chunks'],[])
        def newer(url,body=None):
            value,meta=fetch(url,body)
            if body is None:value['release']='111.0'
            return value,meta
        with patch.object(domains,'_fetch',side_effect=newer):new=self.run_lookup()
        self.assertNotEqual(new['analysis_id'],report['analysis_id'])

    def test_unexpected_md5_or_invalid_coordinates_fail_without_successful_cache(self):
        for bad in ('md5','coordinates'):
            def broken(url,body=None):
                value,meta=fetch(url,body)
                if body is not None:
                    if bad=='md5':value['results'][0]['md5']='A'*32
                    else:value['results'][0]['matches'][0]['locations'][0]['end']=999999
                return value,meta
            data=payload();data['source']=bad
            with patch.object(domains,'_fetch',side_effect=broken),self.assertRaises(ValueError):self.run_lookup(data)
        for manifest in self.root.glob('results/protein-domains/*/manifest.json'):
            state=json.loads(manifest.read_text());self.assertEqual(state['status'],'failed');self.assertEqual(state['chunks'],{})

    def test_stop_between_chunks_reuses_only_completed_chunk(self):
        alphabet='ACDEFGHIKLMNPQRSTVWY'
        data=payload();data['query_fasta']=''.join(f'>p{i}\nMALWMRLLPLLALLALWGPGPG{alphabet[i//20]}{alphabet[i%20]}\n' for i in range(101))
        with job_scope(uuid.uuid4().hex) as control:
            def stop(event):
                if event['status']=='running' and event['chunk']==2:control.cancel()
            with self.assertRaises(RunCancelled):self.run_lookup(data,emit=stop)
        resumed=self.run_lookup(data)
        self.assertEqual(resumed['execution']['reused_chunks'],['0'])
        self.assertEqual(resumed['record_count'],101)
        self.assertEqual(resumed['unknown_record_count'],101)

    def test_saved_screening_candidate_is_read_directly_and_tampering_blocks_lookup(self):
        data=payload()
        screen=run_screen({'query_fasta':data['query_fasta'], 'reference_fasta':(FIXTURES/'protein-screening/P00780.faa').read_text(),
                          'query_source':'UniProt fixtures','reference_source':'UniProt P00780','objective':'find subtilisin','threads':2},root=self.root)
        selected=self.run_lookup({'analysis_id':screen['analysis_id'],'protein_id':'sp|P00782|SUBT_BACAM'})
        self.assertEqual(selected['record_count'],1)
        self.assertEqual(selected['source']['analysis_id'],screen['analysis_id'])
        previous={'workflow':'protein_screen','input':{'query_fasta':data['query_fasta'],
            'reference_fasta':(FIXTURES/'protein-screening/P00780.faa').read_text(),
            'query_source':'UniProt fixtures','reference_source':'UniProt P00780','objective':'find subtilisin','threads':2}}
        followup={'analysis_id':screen['analysis_id'],'protein_id':'sp|P00782|SUBT_BACAM'}
        with patch.object(api,'ROOT',self.root),patch.object(api,'resolve_inputs',side_effect=lambda p:domains.resolve_inputs(p,root=self.root)):
            self.assertTrue(api._matches_screen_source(previous,followup))
            changed=copy.deepcopy(previous);changed['input']['objective']='another objective'
            self.assertFalse(api._matches_screen_source(changed,followup))
        (self.root/'results/protein-screenings'/screen['analysis_id']/'candidates.faa').write_text('modified')
        with self.assertRaises(ValueError):self.run_lookup({'analysis_id':screen['analysis_id']})

    def test_wrong_model_identifiers_names_and_activity_claims_are_not_published_as_evidence(self):
        report=self.run_lookup()
        draft='PF00040 是 Peptidase_S8，活性提高100倍。PF00082 是胶原结构域。'
        reviewed=review_answer(tool_state(report),draft)
        self.assertIn('PF00082 · Subtilase family',reviewed['final_text'])
        self.assertNotIn('100倍',reviewed['final_text'])
        self.assertNotIn('是胶原结构域',reviewed['final_text'])
        self.assertEqual(reviewed['original_text'],draft)
        # PF00040 was an explicit query, so its actual not-observed status is
        # reported; the model's invented accession/name association is absent.
        self.assertIn('PF00040：该预计算记录未观察到匹配',reviewed['final_text'])
        self.assertIn('没有运行新的 InterProScan',reviewed['final_text'])

    def test_unverified_domain_claims_after_blast_are_withheld_and_failure_is_not_success(self):
        screening={'status':'succeeded','analysis_id':'a'*64,'query_count':3,'candidate_count':0,'candidates':[]}
        reviewed=review_answer(tool_state(screening,'run_protein_homology_screen'),'确认含有PF00040胶原酶域')
        self.assertNotIn('PF00040',reviewed['final_text'])
        self.assertIn('结构域、催化位点及实验活性尚未确认',reviewed['final_text'])
        failed=review_answer(tool_state({'status':'failed','kind':'protein_domain_evidence'}),'结构域核验成功')
        self.assertEqual(failed['run_status'],'failed')
        self.assertIn('未完成',failed['final_text'])


class DomainAPITests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=api.RunStore(Path(self.temp.name)/'state.sqlite3')
        self.addCleanup(patch.stopall);patch.object(api,'store',self.store).start()
        self.submit=patch.object(api,'_submit').start()

    def test_domain_task_preserves_inputs_continues_without_openfold_and_is_not_shortcut_prediction(self):
        created=api.create_domain_run(api.DomainRunRequest(**payload()))
        spec=saved_spec(self.store.get(created.run_id));self.assertEqual(spec['workflow'],'protein_domains')
        self.store.update(created.run_id,status='stopped')
        child=api.continue_run(created.run_id);handoff=self.submit.call_args.args[2]
        with patch.object(api,'_execute_domain_spec') as execute,patch.object(api,'_execute_openfold_spec') as fold:api._execute_continuation(child.run_id,handoff)
        execute.assert_called_once();fold.assert_not_called()
        self.assertEqual(execute.call_args.args[1]['input']['query_fasta'],payload()['query_fasta'])
        self.assertIsNone(api._inline_structure_sequence('核验这个蛋白结构域\n'+(FIXTURES/'protein-screening/P00782.faa').read_text()))

    def test_ordinary_conversation_has_no_computational_evidence_claim(self):
        from marine_agent.turn_planning import TurnPlan
        record=self.store.create('你好','general')
        state={'messages':[{'type':'ai','content':'你好，很高兴与你交流。'}]}
        class FakeAgent:
            def stream(self,*args,**kwargs):yield state
        with patch.object(api,'plan_turn',return_value=TurnPlan(action='reply',reason='问候')),patch.object(api,'build_agent',return_value=FakeAgent()):api._controlled_execute(record['run_id'],api._execute_run,'你好','general')
        result=self.store.get(record['run_id'])['result']
        self.assertEqual(result['final_text'],'你好，很高兴与你交流。')
        self.assertEqual(result['scientific_review']['decision'],'conversation_interpretation')
        self.assertEqual(api.get_run_report(record['run_id'])['execution']['kind'],'dialogue')

    def test_agent_public_report_uses_real_evidence_keeps_draft_and_logs_parallel_calls(self):
        from marine_agent.turn_planning import TurnPlan
        with patch.object(domains,'_fetch',side_effect=fetch):report=domains.annotate(payload(),root=Path(self.temp.name))
        record=self.store.create('核验结构域','general')
        state={'messages':[{'type':'ai','id':'calls','tool_calls':[
            {'name':'get_uniprot_sequence','args':{'accession':'P00782'}},
            {'name':'annotate_protein_domains','args':payload()}]},
            {'type':'tool','name':'annotate_protein_domains','id':'result','content':json.dumps(report)},
            {'type':'ai','id':'draft','content':'PF00040是蛋白酶域，已证明活性提高100倍'}]}
        class FakeAgent:
            def stream(self,*args,**kwargs):yield state
        with patch.object(api,'plan_turn',return_value=TurnPlan(action='execute',reason='结构域核验')),patch.object(api,'build_agent',return_value=FakeAgent()):api._controlled_execute(record['run_id'],api._execute_run,'核验结构域','general')
        result=self.store.get(record['run_id'])['result']
        self.assertIn('PF00082 · Subtilase family',result['final_text'])
        self.assertNotIn('100倍',result['final_text'])
        self.assertIn('100倍',result['scientific_review']['original_text'])
        calls=[e['payload']['tool'] for e in self.store.events(record['run_id']) if e['event_type']=='tool_call']
        self.assertEqual(calls,['get_uniprot_sequence','annotate_protein_domains'])

    def test_failed_domain_review_keeps_public_failure_report_and_emits_failure(self):
        from marine_agent.turn_planning import TurnPlan
        record=self.store.create('核验结构域','general')
        state={'messages':[{'type':'tool','name':'annotate_protein_domains',
                            'content':json.dumps({'status':'failed','error':'lookup unavailable'})},
                           {'type':'ai','content':'核验成功，具有胶原酶活性。'}]}
        class FakeAgent:
            def stream(self,*args,**kwargs):yield state
        with patch.object(api,'plan_turn',return_value=TurnPlan(action='execute',reason='结构域核验')),patch.object(api,'build_agent',return_value=FakeAgent()):
            api._controlled_execute(record['run_id'],api._execute_run,'核验结构域','general')
        stored=self.store.get(record['run_id'])
        self.assertEqual(stored['status'],'failed')
        self.assertEqual(stored['current_stage'],'failed')
        report=api.get_run_report(record['run_id'])
        self.assertIn('未完成',report['final_text'])
        self.assertEqual(report['scientific_review']['decision'],'domain_lookup_failed')
        self.assertIn('核验成功',report['scientific_review']['original_text'])
        self.assertEqual(self.store.events(record['run_id'])[-1]['event_type'],'failed')

    def test_partial_domain_status_uses_sequence_evidence_label(self):
        from marine_agent.task_lifecycle import ConversationTasks
        source=api.create_domain_run(api.DomainRunRequest(**payload()))
        self.store.update(source.run_id,status='partial')
        record=self.store.create('目前执行怎么样','general',conversation_id=source.conversation_id)
        api._reply_task_status(record['run_id'],ConversationTasks(self.store,record['run_id']),0)
        text=self.store.get(record['run_id'])['result']['final_text']
        self.assertIn('部分序列证据未知',text)
        self.assertNotIn('配体',text)


if __name__=='__main__':unittest.main(verbosity=2)
