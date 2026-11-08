import json,hashlib,tempfile,unittest
from pathlib import Path
from marine_agent.result_contract import verified_structure, trusted_structure_manifest, execution_receipt
from tests.test_vina_conversation import api

class ResultContractTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
  d=self.root/'results/model';d.mkdir(parents=True);self.data=d/'structure.cif';self.data.write_text('real coordinates')
  self.manifest=d/'manifest.json';self.manifest.write_text(json.dumps({'schema_version':1,'id':'test','provenance':{'model':'OpenFold3'},'structure':{'url':'structure.cif','sha256':hashlib.sha256(self.data.read_bytes()).hexdigest()}}))
  self.record={'id':'real-run','mode':'manual_openfold','status':'succeeded','result':{'agent':{'manual_openfold':{'viewer_manifest':'results/model/manifest.json'}}}}
 def test_ready_requires_real_tool_files_and_digest(self):
  self.assertEqual(verified_structure(self.record,self.root,api._find_viewer_manifest)['status'],'ready')
  self.data.write_text('tampered')
  bad=verified_structure(self.record,self.root,api._find_viewer_manifest)
  self.assertEqual(bad['status'],'unavailable');self.assertFalse(execution_receipt(self.record,bad)['artifact_verified'])
 def test_model_and_user_prose_cannot_create_result(self):
  for role in ('ai','human'):
   agent={'messages':[{'type':role,'name':'run_openfold_sequence','content':json.dumps({'viewer_manifest':'results/model/manifest.json'})}]}
   self.assertIsNone(trusted_structure_manifest(agent,api._find_viewer_manifest))
 def test_failed_confidence_file_cannot_be_ready(self):
  m=json.loads(self.manifest.read_text());m['confidence']={'url':'missing.json'};self.manifest.write_text(json.dumps(m))
  self.assertEqual(verified_structure(self.record,self.root,api._find_viewer_manifest)['status'],'unavailable')
 def test_dialogue_success_is_not_computation_or_validation(self):
  r={'id':'chat','mode':'general','status':'succeeded','result':{'requested_action':'reply','final_text':'全部验证通过'}}
  receipt=execution_receipt(r);self.assertEqual(receipt['computation_status'],'not_requested');self.assertFalse(receipt['artifact_verified'])
  r['result']['requested_action']='execute';self.assertEqual(execution_receipt(r)['computation_status'],'not_started')
 def test_method_rejection_is_not_scientific_execution(self):
  r={'id':'x','mode':'general','status':'succeeded','result':{'scientific_review':{'decision':'method_rejected'}}}
  self.assertEqual(execution_receipt(r)['computation_status'],'not_performed')
 def test_intent_analysis_is_not_scientific_execution(self):
  r={'id':'x','mode':'general','status':'partial','result':{'requested_action':'execute','agent':{'messages':[{'type':'tool','name':'analyze_request_intent','content':'{"status":"succeeded"}'}]}}}
  self.assertEqual(execution_receipt(r)['computation_status'],'not_started')

class SubmissionTests(unittest.TestCase):
 def setUp(self):
  from unittest.mock import patch
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
  self.store=api.RunStore(Path(self.tmp.name)/'state.sqlite3')
  self.patch=patch.object(api,'store',self.store);self.patch.start();self.addCleanup(self.patch.stop)
  self.submit=patch.object(api,'_submit');self.mock=self.submit.start();self.addCleanup(self.submit.stop)
 def test_same_key_survives_reopen_and_never_launches_twice(self):
  request=api.RunRequest(prompt='你好',request_id='same-request-01')
  first=api.create_run(request)
  api.store=api.RunStore(self.store.path)
  second=api.create_run(request)
  self.assertEqual(first.run_id,second.run_id);self.assertEqual(self.mock.call_count,1)
 def test_conflicting_inputs_and_parallel_conversation_are_rejected(self):
  first=api.create_run(api.RunRequest(prompt='你好',request_id='same-request-02'))
  with self.assertRaises(api.HTTPException) as err:api.create_run(api.RunRequest(prompt='不同',request_id='same-request-02'))
  self.assertEqual(err.exception.status_code,409)
  with self.assertRaises(api.HTTPException) as err:api.create_run(api.RunRequest(prompt='重复',conversation_id=first.conversation_id))
  self.assertEqual(err.exception.status_code,409);self.assertEqual(self.mock.call_count,1)

class MethodContractTests(unittest.TestCase):
 def test_language_independent_endpoint_rejects_activity_not_binding(self):
  from marine_agent.collagenase_tools import assess_collagenase_method
  for objective in ['collagenase (collagen-degrading protease) activity','胶原酶活性','activité enzymatique']:
   r=json.loads(assess_collagenase_method.invoke({'method':'Vina','objective':objective,'endpoint':'enzyme_activity'}))
   self.assertEqual(r['decision'],'reject')
  r=json.loads(assess_collagenase_method.invoke({'method':'Vina','objective':'胶原酶与小分子抑制剂结合姿态','endpoint':'binding_pose'}))
  self.assertNotEqual(r['decision'],'reject')
