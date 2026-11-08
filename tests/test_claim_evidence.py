import base64,csv,io,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from fastapi import HTTPException
from tests.test_vina_conversation import api
from marine_agent.access_control import account_scope
from marine_agent.research_projects import ProjectBrief
from marine_agent.candidate_ledger import CandidateBatch,CandidateInput,ImportInput,candidate_tools
from marine_agent.claim_evidence import Claim,ClaimPlan,collect_bundle,compile_answer,render_claim,calculation_sources,public_evidence
from marine_agent.experimental_data import COLUMNS

class ClaimCompilerTests(unittest.TestCase):
    def source(self,identifier='E1',**extra):
        s={'id':identifier,'record_id':identifier,'kind':'experiment','subject':'BglB','metric':'kcat','value':880/60,'unit':'s^-1','raw_value':'880','raw_unit':'min^-1','qualifier':'=','sd':10/60,'n':None,'data_kind':'summary',
           'conditions':{'substrate':'pNPG','temperature_c':None,'ph':7.5,'protocol':'A420','control':'WT','batch':None,'construct':'C-His','source':'doi:public'},'source':'doi:public'}
        s.update(extra);return s
    def compile(self,claims,sources=None,draft='未经核对的草稿'):
        return compile_answer('核对记录',draft,{'version':'test','run_id':'run','sources':sources or [self.source()]},invoke=lambda _: {'claims':claims})
    def test_provider_json_contract_and_malformed_response(self):
        from types import SimpleNamespace
        from marine_agent.claim_evidence import _ask_plan
        class Model:
            content='```json\n{"claims":[{"op":"report","evidence_ids":["E1"]}]}\n```'
            def bind(self, **kwargs):return self
            def invoke(self, messages):return SimpleNamespace(content=self.content)
        model=Model()
        with patch('marine_agent.agent.configured_model',return_value=model):
            plan=_ask_plan({'question':'值','evidence':[self.source()]})
            self.assertEqual(plan.claims[0].op,'report')
            model.content='{"claims":[{"op":"report","evidence_ids":["E1"],"field":"metric_scope"}]}'
            with self.assertRaises(ValueError):_ask_plan({})
            model.content='{"claims":['
            r=compile_answer('值','伪造工业效果',{'sources':[self.source()]})
            self.assertEqual(r['decision'],'records_only_fallback');self.assertNotIn('伪造工业效果',r['final_text'])
    def test_provider_repairs_serialization_once_without_bypassing_validation(self):
        from types import SimpleNamespace
        from marine_agent.claim_evidence import _ask_plan
        class Model:
            calls=0
            def bind(self, **kwargs):return self
            def invoke(self, messages):
                self.calls+=1
                return SimpleNamespace(content='{"claims":[' if self.calls==1 else '{"claims":[{"op":"report","evidence_ids":["E1"]}]}')
        model=Model()
        with patch('marine_agent.agent.configured_model',return_value=model):
            plan=_ask_plan({});self.assertEqual(plan.provider_attempts,2);self.assertEqual(model.calls,2)
        r=self.compile([{'op':'report','evidence_ids':['E404']}]);self.assertEqual(r['decision'],'records_only_fallback')
    def test_wrong_number_and_unit_corrected_without_republishing_draft(self):
        r=self.compile([{'op':'report','evidence_ids':['E1'],'stated_value':880,'stated_unit':'s^-1'}],draft='催化常数为880 s^-1且工业效果优秀')
        self.assertIn('14.6667 s^-1',r['final_text']);self.assertNotIn('工业效果优秀',r['final_text']);self.assertIn('draft_numeric_value_replaced',r['claims'][0]['checks']);self.assertEqual(r['original_text'],'催化常数为880 s^-1且工业效果优秀')
    def test_unknown_batch_not_manufactured_from_source_or_draft(self):
        r=self.compile([{'op':'conditions','field':'batch','evidence_ids':['E1'],'stated_condition':'同批实验'}])
        self.assertIn('批次：未报告',r['final_text']);self.assertIn('draft_condition_replaced',r['claims'][0]['checks'])
    def test_same_batch_requires_both_known(self):
        a=self.source();b=self.source('E2',subject='N220Y')
        r=self.compile([{'op':'compare','relation':'same_batch','evidence_ids':['E1','E2']}],[a,b]);self.assertIn('不能确认同批',r['final_text'])
        for s in [a,b]:s['conditions']['batch']='B1'
        r=self.compile([{'op':'compare','relation':'same_batch','evidence_ids':['E1','E2']}],[a,b]);self.assertIn('批次一致',r['final_text']);self.assertIn('未独立核验',r['final_text'])
    def test_multi_record_batch_comparison_checks_all_records(self):
        sources=[self.source('E1'),self.source('E2'),self.source('E3')]
        for s in sources:s['conditions']['batch']='B1'
        sources[-1]['conditions']['batch']='B2'
        r=self.compile([{'op':'compare','relation':'same_batch','evidence_ids':['E1','E2','E3']}],sources)
        self.assertIn('批次不同',r['final_text']);self.assertFalse(r['rejected_claims'])
        r=self.compile([{'op':'compare','relation':'ratio','evidence_ids':['E1','E2','E3']}],sources)
        self.assertEqual(r['rejected_claims'][0]['reason'],'comparison_requires_two_records')
    def test_ratio_blocked_for_incomplete_conditions(self):
        r=self.compile([{'op':'compare','relation':'ratio','evidence_ids':['E1','E2']}],[self.source(),self.source('E2',value=2)])
        self.assertIn('暂不作倍数判断',r['final_text']);self.assertIn('批次未报告',r['final_text']);self.assertEqual(r['claims'][0]['status'],'limited')
    def test_ratio_exact_conditions_units_and_n_checked(self):
        a=self.source(value=4,n=3);b=self.source('E2',subject='variant',value=2,n=3)
        for s in [a,b]:s['conditions'].update(batch='B1',temperature_c=25)
        r=self.compile([{'op':'compare','relation':'ratio','evidence_ids':['E1','E2']}],[a,b]);self.assertIn('为 2；',r['final_text']);self.assertIn('不是显著性检验',r['final_text'])
        b['conditions']['substrate']='collagen';r=self.compile([{'op':'compare','relation':'ratio','evidence_ids':['E1','E2']}],[a,b]);self.assertIn('底物不一致',r['final_text'])
    def test_missing_not_zero_and_censored_not_exact(self):
        for source,phrase in [(self.source(value=None,qualifier='missing'),'不能按零处理'),(self.source(value=1/6,raw_value='10',qualifier='<'),'这是界限')]:
            r=self.compile([{'op':'report','evidence_ids':['E1']}],[source]);self.assertIn(phrase,r['final_text']);self.assertEqual(r['claims'][0]['status'],'limited')
    def test_grouped_report_keeps_each_value_and_missing_state(self):
        sources=[self.source('E1',subject='G12N',value=None,qualifier='missing'),self.source('E2',subject='Q19A',value=1/6,raw_value='10',qualifier='<')]
        r=self.compile([{'op':'report','evidence_ids':['E1','E2']}],sources)
        self.assertFalse(r['rejected_claims']);self.assertIn('不能按零处理',r['final_text']);self.assertIn('<0.166667',r['final_text']);self.assertIn('不是精确测量值',r['final_text'])
    def test_activity_interpretation_explains_missing_measurement(self):
        r=self.compile([{'op':'interpret','relation':'activity','evidence_ids':['E1']}],[self.source(subject='G12N',value=None,qualifier='missing')])
        self.assertIn('G12N 未报告的测量值不能按零活性',r['final_text'])
    def test_km_not_promoted_to_affinity(self):
        r=self.compile([{'op':'interpret','relation':'binding_affinity','evidence_ids':['E1']}],[self.source(metric='km')]);self.assertIn('不能直接认定',r['final_text']);self.assertIn('Kd',r['final_text']);self.assertTrue(r['claims'][0]['rule_sources'])
    def test_industrial_claim_is_bounded_to_substrate(self):
        r=self.compile([{'op':'interpret','relation':'industrial_suitability','evidence_ids':['E1']}]);self.assertIn('pNPG',r['final_text']);self.assertIn('不能单独确认',r['final_text'])
    def test_forged_reference_and_operator_not_published(self):
        r=self.compile([{'op':'report','evidence_ids':['E999']}]);self.assertEqual(r['decision'],'records_only_fallback');self.assertNotIn('E999',r['final_text']);self.assertEqual(r['rejected_claims'][0]['reason'],'unknown_evidence_id')
        r=self.compile([{'op':'assert_best_enzyme','evidence_ids':['E1']}]);self.assertEqual(r['decision'],'records_only_fallback')
    def test_provider_timeout_is_record_fallback_not_unchecked_draft(self):
        def fail(_):raise TimeoutError('private-provider-data')
        r=compile_answer('解释','这台机器比其他所有都快',{'version':'test','sources':[self.source()]},invoke=fail)
        self.assertEqual(r['fallback_reason'],'TimeoutError');self.assertNotIn('其他所有',r['final_text']);self.assertNotIn('private-provider-data',json.dumps(r))
    def test_empty_bundle_never_calls_provider_or_invents_facts(self):
        with patch('marine_agent.claim_evidence._ask_plan') as model:
            r=compile_answer('酶活是多少','42 U/mg',{'sources':[]});model.assert_not_called();self.assertIn('没有取得',r['final_text'])
    def test_repeated_selection_deduplicated_citations_are_actual_ids(self):
        r=self.compile([{'op':'report','evidence_ids':['E1']}]*2);self.assertEqual(len(r['claims']),1);self.assertIn('[1](#evidence-E1)',r['final_text'])
    def test_multiple_condition_records_and_mismatched_metric_fail_safely(self):
        r=self.compile([{'op':'conditions','field':'batch','evidence_ids':['E1','E2']}],[self.source(),self.source('E2')]);self.assertFalse(r['rejected_claims']);self.assertIn('批次：未报告',r['final_text'])
        r=self.compile([{'op':'compare','relation':'ratio','evidence_ids':['E1','E2']}],[self.source(),self.source('E2',metric='km',unit='mM')]);self.assertIn('指标或单位不一致',r['final_text'])
    def test_computational_score_does_not_support_activity_or_experimental_ratio(self):
        a=self.source(kind='calculation',metric='vina_affinity',value=-7,unit='kcal/mol');b=self.source('E2',kind='calculation',metric='vina_affinity',value=-8,unit='kcal/mol')
        r=self.compile([{'op':'interpret','evidence_ids':['E1']}],[a]);self.assertIn('不是实测',r['final_text'])
        r=self.compile([{'op':'compare','relation':'ratio','evidence_ids':['E1','E2']}],[a,b]);self.assertIn('不能用于实验效能比较',r['final_text'])

class ClaimStoreTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.store=api.RunStore(Path(self.tmp.name)/'state.sqlite3')
        self.scope=account_scope('alice');self.scope.__enter__();self.addCleanup(self.scope.__exit__,None,None,None)
        self.pid=self.store.create_project(ProjectBrief(name='证据验收'))['id'];self.cid=self.store.new_project_conversation(self.pid)['conversation_id']
        self.store.register_candidates(self.pid,CandidateBatch(candidates=[CandidateInput(key='A',name='蛋白A',sequence='ACDEFG',source='public')]))
        row=dict.fromkeys(COLUMNS,'');row.update(candidate_key='A',metric='kcat',value='60',unit='min^-1',data_kind='summary',source='doi:test')
        out=io.StringIO();w=csv.DictWriter(out,fieldnames=COLUMNS);w.writeheader();w.writerow(row)
        p=self.store.preview_experiments(self.pid,ImportInput(filename='data.csv',content_base64=base64.b64encode(out.getvalue().encode()).decode()));self.store.commit_experiments(self.pid,p['token'])
        self.run=self.store.create('读取实验','general',conversation_id=self.cid)
    def state(self):
        snapshot=self.store.run_project_context(self.run['run_id']);tool=candidate_tools(self.store,snapshot)[0];value=json.loads(tool.invoke({'candidate_key':'A'}))
        return {'messages':[{'type':'tool','name':'inspect_candidate_evidence','content':json.dumps(value)}]}
    def test_rehydrate_from_persistent_store_not_modified_tool_values(self):
        state=self.state();value=json.loads(state['messages'][0]['content']);value['measurements'][0]['value']=123456;state['messages'][0]['content']=json.dumps(value)
        bundle=collect_bundle(self.store,self.run['run_id'],state);self.assertEqual(bundle['sources'][0]['value'],1);self.assertEqual(bundle['sources'][0]['ledger_revision'],2);self.assertEqual(len(bundle['sources'][0]['record_id']),64)
    def test_exact_source_file_and_snapshot_survive_restart(self):
        bundle=collect_bundle(self.store,self.run['run_id'],self.state());review=compile_answer('值','草稿',bundle,invoke=lambda _:{'claims':[{'op':'report','evidence_ids':['E1']}]})
        self.store.update(self.run['run_id'],status='succeeded',result_json=json.dumps({'claim_evidence':review,'final_text':review['final_text']}))
        restarted=api.RunStore(self.store.path)
        with patch.object(api,'store',restarted):
            evidence=api.get_claim_evidence(self.run['run_id']);self.assertEqual(evidence['sources'][0]['file_sha256'],bundle['sources'][0]['file_sha256']);self.assertNotIn('original_text',evidence)
            self.assertTrue(api.get_run_report(self.run['run_id'])['has_claim_evidence']);self.assertTrue(restarted.get_conversation(self.cid)['runs'][0]['has_claim_evidence'])
    def test_corrupt_original_file_or_measurement_is_not_promoted(self):
        state=self.state()
        with self.store._connect() as db:db.execute("UPDATE experiment_imports SET raw=?",(b'corrupt',))
        bundle=collect_bundle(self.store,self.run['run_id'],state)
        self.assertEqual(bundle['sources'],[]);self.assertEqual(bundle['integrity_issues'][0]['reason'],'original_file_hash_mismatch')
    def test_changed_value_with_old_row_hash_is_rejected(self):
        state=self.state()
        with self.store._connect() as db:
            r=db.execute('SELECT payload FROM experimental_measurements').fetchone();value=json.loads(r[0]);value['value']=100
            db.execute('UPDATE experimental_measurements SET payload=?',(json.dumps(value),))
        bundle=collect_bundle(self.store,self.run['run_id'],state);self.assertEqual(bundle['sources'],[]);self.assertEqual(bundle['integrity_issues'][0]['reason'],'measurement_record_hash_mismatch')
    def test_owner_cannot_read_foreign_evidence(self):
        with account_scope('bob'),patch.object(api,'store',self.store):
            with self.assertRaises(HTTPException) as e:api.get_claim_evidence(self.run['run_id'])
            self.assertEqual(e.exception.status_code,404)
    def test_ordinary_chat_skips_claim_review(self):
        self.assertIsNone(collect_bundle(self.store,self.run['run_id'],{'messages':[{'type':'ai','content':'你好'}]}))
    def test_followup_can_use_prior_checked_sources_without_claiming_current_version(self):
        bundle=collect_bundle(self.store,self.run['run_id'],self.state());review=compile_answer('值','草稿',bundle,invoke=lambda _:{'claims':[{'op':'report','evidence_ids':['E1']}]})
        self.store.update(self.run['run_id'],status='succeeded',result_json=json.dumps({'claim_evidence':review}));next_run=self.store.create('同批吗','general',conversation_id=self.cid)
        next_bundle=collect_bundle(self.store,next_run['run_id'],{'messages':[]},expect_evidence=True)
        self.assertEqual(next_bundle['sources'][0]['inherited_from_run_id'],self.run['run_id']);self.assertEqual(next_bundle['sources'][0]['ledger_revision'],2)
    def test_vina_score_read_from_persisted_result_not_report_prose(self):
        run=self.store.create('dock','vina',conversation_id=self.cid);self.store.update(run['run_id'],status='succeeded',result_json=json.dumps({'final_text':'分数-999酶活高','vina':{'results':[{'id':'ligand','status':'succeeded','scores':[{'affinity_kcal_mol':-7.2}]}]}}))
        result=calculation_sources(self.store,run['run_id'],Path(self.tmp.name),lambda _:None);self.assertEqual(result[0]['value'],-7.2);self.assertEqual(result[0]['kind'],'calculation')
    def test_pipeline_persists_checked_text_and_api_citations(self):
        from marine_agent.turn_planning import TurnPlan
        state=self.state();state['messages'].append({'type':'ai','content':'同批实验，催化常数60 s^-1'})
        class Agent:
            def stream(self,*args,**kwargs):yield state
        with patch.object(api,'store',self.store),patch.object(api,'build_agent',return_value=Agent()),patch.object(api,'plan_turn',return_value=TurnPlan(action='inspect',reason='读取记录',inspection='results',context_scope='project')),patch('marine_agent.claim_evidence._ask_plan',return_value={'claims':[{'op':'report','evidence_ids':['E1'],'stated_value':60},{'op':'conditions','field':'batch','evidence_ids':['E1'],'stated_condition':'同批'}]}):
            api._controlled_execute(self.run['run_id'],api._execute_run,'读取记录','general');report=api.get_run_report(self.run['run_id']);self.assertTrue(report['has_claim_evidence']);self.assertIn('1 s^-1',report['final_text']);self.assertIn('批次：未报告',report['final_text']);self.assertEqual(self.store.get_conversation(self.cid)['messages'][-1]['content'],report['final_text'])

if __name__=='__main__':unittest.main()
