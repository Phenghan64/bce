"""Task handoffs are scoped, durable, deterministic and preserve real Vina poses."""
import base64
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tests.test_vina_conversation import api, ROOT
from marine_agent.execution_control import job_scope
from marine_agent.task_lifecycle import ConversationTasks, saved_spec
from marine_agent import vina_tools


class ContinuationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT/'results')
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.store = api.RunStore(self.path/'state.sqlite3')
        self.addCleanup(patch.stopall)
        patch.object(api, 'store', self.store).start()
        self.submit = patch.object(api, '_submit').start()

    def stopped(self, accession='J7H3P9'):
        source = self.store.create('预测结构', 'manual_openfold', manual={'accession':accession,'ligands':[{'smiles':'CCO','count':1}]})
        self.store.update(source['run_id'], status='stopped', current_stage='stopped')
        return source

    def turn(self, source, text='继续'):
        record = self.store.create(text, 'general', conversation_id=source['conversation_id'])
        self.store.update(record['run_id'], status='running')
        return record, ConversationTasks(self.store, record['run_id'])

    def test_scoped_state_survives_intervening_chat_and_never_exposes_another_conversation(self):
        source=self.stopped();other=self.stopped('O74213')
        chat=self.store.create('你有执行吗','general',conversation_id=source['conversation_id'])
        self.store.update(chat['run_id'],status='succeeded')
        current,service=self.turn(source)
        self.assertEqual(service.snapshot(other['run_id'])['status'],'not_found')
        self.assertEqual(service.prepare(other['run_id'])['status'],'not_found')
        self.assertEqual(service.prepare()['source_run_id'],source['run_id'])
        self.assertEqual(service.handoff['spec']['input']['ligands'],[{'smiles':'CCO','count':1}])
        self.assertEqual(self.store.get(current['run_id'])['continuation_of'],source['run_id'])
        self.assertEqual(self.store.get(source['run_id'])['status'],'stopped')

    def test_endpoint_links_new_attempt_same_conversation_and_duplicate_active_request_is_rejected(self):
        source=self.stopped();created=api.continue_run(source['run_id'])
        self.assertNotEqual(created.run_id,source['run_id'])
        self.assertEqual(created.conversation_id,source['conversation_id'])
        self.assertEqual(self.store.get(created.run_id)['continuation_of'],source['run_id'])
        with self.assertRaises(api.HTTPException) as error:api.continue_run(source['run_id'])
        self.assertEqual(error.exception.status_code,409)
        self.assertEqual(self.submit.call_count,1)

    def test_public_actions_and_poll_payload_follow_latest_attempt(self):
        source=self.stopped()
        self.assertTrue(api.get_run(source['run_id'])['can_continue'])
        child=api.continue_run(source['run_id'])
        self.store.update(child.run_id,status='succeeded')
        rows=self.store.get_conversation(source['conversation_id'])['runs']
        self.assertFalse(any(r['can_continue'] for r in rows))
        self.assertFalse(api.get_run(source['run_id'])['can_continue'])
        self.assertNotIn('execution_spec',api.get_run(child.run_id))
        with self.assertRaises(api.HTTPException):api.continue_run(source['run_id'])
        self.assertEqual(len(self.store.get_conversation(source['conversation_id'])['runs']),2)

    def test_queued_continuation_updates_sidebar_timestamp_and_never_calls_model(self):
        source=self.stopped()
        with self.store._connect() as db:
            db.execute("UPDATE conversations SET updated_at='2000-01-01' WHERE id=?",(source['conversation_id'],))
        with patch.object(api,'plan_turn') as model:
            child=api.continue_run(source['run_id'])
        model.assert_not_called()
        self.assertEqual(self.store.get_conversation(source['conversation_id'])['conversation']['updated_at'],child.created_at)

    def test_model_rate_limit_message_does_not_expose_gateway_payload(self):
        error=RuntimeError('private gateway payload');error.status_code=429
        self.assertIn('限流',api._execution_error(error))
        self.assertNotIn('private',api._execution_error(error))

    def test_long_conversation_still_exposes_the_interrupted_scientific_task(self):
        source=self.stopped()
        for index in range(16):
            chat=self.store.create(f'讨论 {index}','general',conversation_id=source['conversation_id'])
            self.store.update(chat['run_id'],status='succeeded')
        _,service=self.turn(source)
        self.assertIn(source['run_id'],[row['run_id'] for row in service.snapshot()['tasks']])

    def test_completed_latest_attempt_does_not_restart_older_stopped_ancestor(self):
        source=self.stopped();created=api.continue_run(source['run_id'])
        self.store.update(created.run_id,status='succeeded')
        _,service=self.turn(source)
        self.assertEqual(service.prepare()['status'],'not_resumable')

    def test_cancelling_and_finished_tasks_cannot_be_continued(self):
        for status in ['running','queued','cancelling','succeeded']:
            fresh=self.store.create('任务','manual_openfold',manual={'accession':'J7H3P9'})
            self.store.update(fresh['run_id'],status=status)
            with self.assertRaises(api.HTTPException) as error:api.continue_run(fresh['run_id'])
            self.assertEqual(error.exception.status_code,409)

    def test_actual_tool_call_inputs_are_saved_before_execution_and_ambiguous_multijob_is_not_guessed(self):
        current=self.store.create('预测这个','general')
        api._capture_execution_calls(current['run_id'],{'tool_calls':[{'name':'run_openfold_sequence','args':{'sequence':'A'*30,'ligands_json':'[]'}}]})
        self.assertEqual(saved_spec(self.store.get(current['run_id']))['input']['sequence'],'A'*30)
        api._capture_execution_calls(current['run_id'],{'tool_calls':[{'name':'run_openfold_sequence','args':{'sequence':'C'*30}}]})
        self.assertIsNone(saved_spec(self.store.get(current['run_id'])))

    def test_natural_continuation_hands_off_before_model_reply_and_preserves_exact_inputs(self):
        from marine_agent.turn_planning import TurnPlan
        source=self.stopped();current=self.store.create('继续','general',conversation_id=source['conversation_id'])
        with patch.object(api,'plan_turn',return_value=TurnPlan(action='continue',reason='用户请求继续')),patch.object(api,'build_agent') as builder,patch.object(api,'_execute_openfold_spec') as execute:
            api._controlled_execute(current['run_id'],api._execute_run,'继续','general')
        execute.assert_called_once();builder.assert_not_called()
        self.assertEqual(execute.call_args.args[1]['input'],{'accession':'J7H3P9','ligands':[{'smiles':'CCO','count':1}]})
        self.assertEqual(self.store.get(current['run_id'])['continuation_of'],source['run_id'])

    def test_status_question_uses_authoritative_facts_without_starting_computation(self):
        from marine_agent.turn_planning import TurnPlan
        source=self.stopped();current=self.store.create('你有继续操作吗','general',conversation_id=source['conversation_id'])
        with patch.object(api,'plan_turn',return_value=TurnPlan(action='inspect',reason='状态询问')),patch.object(api,'build_agent') as builder,patch.object(api,'_execute_openfold_spec') as execute:
            api._controlled_execute(current['run_id'],api._execute_run,'你有继续操作吗','general')
        execute.assert_not_called();builder.assert_not_called()
        record=self.store.get(current['run_id'])
        self.assertEqual(record['status'],'succeeded')
        self.assertIn('尚未再次启动',record['result']['final_text'])
        self.assertIn('已停止',record['result']['final_text'])

    def test_missing_or_modified_pose_is_recomputed_instead_of_reused(self):
        ligand=b'ligand';path=self.path/'results/vina-runs/old/one.out.pdbqt';path.parent.mkdir(parents=True)
        raw=b'REMARK VINA RESULT: -5.0 0.0 0.0\n'
        path.write_bytes(raw)
        row={'id':'one','status':'succeeded','input_sha256':hashlib.sha256(ligand).hexdigest(),
             'output':str(path.relative_to(self.path)),'output_sha256':hashlib.sha256(raw).hexdigest(),
             'scores':[{'mode':1,'affinity_kcal_mol':-5.0,'rmsd_lb':0.0,'rmsd_ub':0.0}]}
        source={'id':'old','result':{'vina':{'results':[row]}}};config={'ligands':[{'id':'one','bytes':ligand}]}
        with patch.object(api,'ROOT',self.path):
            self.assertEqual(len(api._vina_carried_results(source,config)),1)
            path.write_bytes(raw+b'changed coordinates')
            self.assertEqual(api._vina_carried_results(source,config),[])
            path.unlink();self.assertEqual(api._vina_carried_results(source,config),[])

    def test_real_vina_continues_only_unfinished_ligands_with_same_parameters_and_keeps_3d_files(self):
        reference=ROOT/'results/vina-runs/fd3cad6e-d4c5-4777-bac5-38068fb4ea05'
        def encoded(file):return base64.b64encode(file.read_bytes()).decode()
        request={'receptor':{'name':'receptor.pdbqt','content_base64':encoded(reference/'receptor.pdbqt')},
                 'ligands':[{'id':name,'name':name+'.pdbqt','content_base64':encoded(reference/'ligands'/file)}
                     for name,file in [('caffeine','library-L001-primary-1.pdbqt'),('aspirin','library-L002-primary-1.pdbqt')]],
                 'box':{'center':[15.614,53.38,15.455],'size':[30,30,30]},'params':{'exhaustiveness':8,'cpu':1}}
        source=self.store.create('对接两种配体','vina',manual=request)
        real=vina_tools.run_batch
        with patch.object(api,'ROOT',self.path),patch.object(vina_tools,'ROOT',self.path),patch.object(api,'VINA_WORKER_URL',''),patch.object(vina_tools,'VINA_BINARY',str(ROOT/'.tools/bin/vina')):
            with job_scope(source['run_id']) as control:
                def interrupt(directory,config,emit):
                    config['batch_workers']=1
                    def event(data):
                        emit(data)
                        if data.get('succeeded')==1:control.cancel()
                    return real(directory,config,emit=event)
                with patch.object(api,'run_batch',side_effect=interrupt):api._execute_vina_run(source['run_id'],copy.deepcopy(request))
            before=self.store.get(source['run_id']);self.assertEqual(before['status'],'stopped')
            completed=[r for r in before['result']['vina']['results'] if r['status']=='succeeded']
            self.assertEqual(len(completed),1)
            created=api.continue_run(source['run_id']);handoff=self.submit.call_args.args[2]
            # Stop delivered before the worker starts its pending ligand must
            # still retain the completed pose; a later continuation follows it.
            with patch.object(api,'VINA_WORKER_URL','test-worker'),patch.object(api,'run_vina_worker',return_value={'status':'stopped','results':[]}):
                api._controlled_execute(created.run_id,api._execute_continuation,handoff)
            stopped_again=self.store.get(created.run_id)
            self.assertEqual(stopped_again['status'],'stopped')
            self.assertEqual([r['status'] for r in stopped_again['result']['vina']['results']],['succeeded','skipped'])
            created=api.continue_run(created.run_id);handoff=self.submit.call_args.args[2]
            with patch.object(api,'run_batch',wraps=real) as execute:
                api._controlled_execute(created.run_id,api._execute_continuation,handoff)
            after=self.store.get(created.run_id)
            self.assertEqual(after['status'],'succeeded',after.get('error'))
            self.assertEqual([r['id'] for r in execute.call_args.args[1]['ligands']],['aspirin'])
            self.assertEqual(execute.call_args.args[1]['params']['exhaustiveness'],8)
            vina=after['result']['vina'];self.assertEqual(len(vina['ranking']),2)
            kept=next(r for r in vina['results'] if r['id']=='caffeine')
            self.assertEqual(kept['scores'],completed[0]['scores'])
            self.assertTrue((self.path/kept['output']).is_file())
            self.assertTrue((self.path/'results/vina-runs'/created.run_id/'receptor.pdbqt').is_file())
            self.assertEqual(self.store.get(source['run_id'])['result'],before['result'])


if __name__=='__main__':unittest.main(verbosity=2)
