import json,tempfile,unittest,os
from pathlib import Path
from unittest.mock import patch
from marine_agent.scientific_knowledge import references,review_explanation,display_prose
from marine_agent.turn_planning import TurnPlan
from marine_agent.execution_control import RunCancelled
from tests.test_vina_conversation import api

class ScientificKnowledgeTests(unittest.TestCase):
 def test_review_limits_are_applied_when_transport_is_constructed(self):
  from marine_agent.scientific_knowledge import _review_model
  with patch.dict(os.environ,{'DEEP_AGENT_BASE_URL':'https://example.invalid/v1','DEEP_AGENT_API_KEY':'test-only'}),patch('langchain_openai.ChatOpenAI') as client:
   _review_model()
  self.assertEqual(client.call_args.kwargs['timeout'],25)
  self.assertEqual(client.call_args.kwargs['max_retries'],0)
  self.assertEqual(client.call_args.kwargs['max_completion_tokens'],2400)
 def test_reference_identifiers_stay_out_of_user_prose(self):
  self.assertEqual(display_prose('定义（plddt.definition, pldddt.granularity）。'), '定义。')
  self.assertEqual(display_prose('PAE（单位 Å），pae[i,j]。'), 'PAE（单位 Å），pae[i,j]。')
 def test_retrieval_scopes_versions_and_unknowns(self):
  b=references(['plddt','plddt']);self.assertEqual(len(b['cards']),1)
  self.assertEqual(len(b['sha256']),64);self.assertIn('eyoson_export',b['sources'])
  self.assertFalse(references(['made_up'])['cards']);self.assertNotIn('vina_faq',b['sources'])
 def test_server_sources_cannot_be_fabricated_by_reviewer(self):
  b=references(['plddt'])
  for invalid in [None,{'verdict':'supported','answer':'已证明活性','fact_ids':['invented']}, {'verdict':'supported','answer':'见 https://fake.example','fact_ids':['plddt.definition']}]:
   r=review_explanation('解释pLDDT','几何自检',b,invoke=lambda p:invalid)
   self.assertEqual(r['decision'],'reference_facts_fallback');self.assertNotIn('fake.example',r['final_text'])
   self.assertNotIn('已证明活性',r['final_text'])
 def test_revision_keeps_auditable_draft_and_exact_sources(self):
  r=review_explanation('pLDDT是什么','几何自检',references(['plddt']),invoke=lambda p:{'verdict':'revised','answer':'pLDDT 是预测的局部准确性。','fact_ids':['plddt.definition'],'issues':['预测与实际lDDT混淆']})
  self.assertEqual(r['original_text'],'几何自检');self.assertEqual([s['id'] for s in r['sources']],['ebi_plddt'])
  self.assertIn('https://www.ebi.ac.uk/',r['final_text'])
 def test_registered_link_in_review_does_not_trigger_fallback(self):
  b=references(['plddt']);url=b['sources']['ebi_plddt']['url']
  r=review_explanation('pLDDT','draft',b,invoke=lambda p:{'verdict':'supported','answer':f'pLDDT 是预测。见[官方资料]({url})。','fact_ids':['plddt.definition']})
  self.assertEqual(r['decision'],'reference_review_supported');self.assertEqual(r['final_text'].count(url),1)
  self.assertIsNone(r['fallback_reason'])
 def test_timeout_falls_back_but_cancel_propagates(self):
  def fail(p):raise TimeoutError('secret provider payload')
  r=review_explanation('PAE','错误',references(['pae']),invoke=fail)
  self.assertEqual(r['decision'],'reference_facts_fallback');self.assertNotIn('secret',r['final_text'])
  def cancel(p):raise RunCancelled('stopped')
  with self.assertRaises(RunCancelled):review_explanation('PAE','',references(['pae']),invoke=cancel)
 def test_unrelated_chat_does_not_invoke_review(self):
  with patch('marine_agent.scientific_knowledge._review_model') as model:
   self.assertIsNone(review_explanation('你好','你好',references([])));model.assert_not_called()
 def test_actual_conversation_pipeline_keeps_review_as_dialogue(self):
  with tempfile.TemporaryDirectory() as tmp,patch.object(api,'store',api.RunStore(Path(tmp)/'state.sqlite3')):
   r=api.store.create('解释pLDDT','general')
   class Agent:
    def stream(self,*args,**kwargs):yield {'messages':[{'type':'ai','content':'pLDDT是几何自检'}]}
   response={'final_text':'pLDDT是局部准确性的预测。','decision':'reference_review_revised','fact_ids':['plddt.definition'],'knowledge_version':'test'}
   with patch.object(api,'plan_turn',return_value=TurnPlan(action='reply',reason='定义',reference_topics=['plddt'])),patch.object(api,'build_agent',return_value=Agent()),patch('marine_agent.scientific_knowledge.review_explanation',return_value=response) as review:
    api._controlled_execute(r['run_id'],api._execute_run,'解释pLDDT','general')
   report=api.get_run_report(r['run_id'])
   self.assertEqual(report['final_text'],response['final_text']);self.assertEqual(report['execution']['kind'],'dialogue')
   self.assertEqual(report['explanation_evidence']['decision'],'reference_review_revised');review.assert_called_once()
   self.assertFalse(api.store.get(r['run_id']).get('operations'))
   self.assertFalse(any(e['event_type']=='stage' and e['payload'].get('name')=='核对解释与来源' for e in api.store.events(r['run_id'])))
if __name__=='__main__':unittest.main()
