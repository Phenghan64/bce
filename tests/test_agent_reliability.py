"""Regression tests for failures observed in the 2026-11-07 live acceptance."""
import io
import json
import os
from pathlib import Path
import tempfile
import threading
import time
import unittest
import uuid
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from tests.test_vina_conversation import api
from marine_agent import turn_planning as planning
from marine_agent.execution_control import job_scope, job_status, control_for
from marine_agent.restart_recovery import reconcile
from marine_agent.task_lifecycle import ConversationTasks


class PlannerTests(unittest.TestCase):
    def test_method_verdict_terminates_invalid_endpoint_but_not_enzyme_ligand_question(self):
        from marine_agent.collagenase_tools import assess_collagenase_method
        from marine_agent.evidence_review import method_rejection
        verdict = assess_collagenase_method.invoke({'method':'Vina', 'objective':'预测胶原降解活性'})
        state = {'messages':[{'type':'tool','name':'assess_collagenase_method','content':verdict}]}
        self.assertEqual(method_rejection(state)['decision'], 'method_rejected')
        legitimate = json.loads(assess_collagenase_method.invoke({'method':'Vina','objective':'研究胶原酶与小分子抑制剂的结合姿态'}))
        self.assertNotEqual(legitimate['decision'],'reject')

    def test_empty_function_result_uses_validated_json_fallback(self):
        model = MagicMock()
        model.with_structured_output.return_value.invoke.return_value = None
        model.invoke.return_value = SimpleNamespace(content='{"action":"reply","reason":"讨论适用性"}')
        with patch.object(planning, 'configured_model', return_value=model):
            plan = planning.plan_turn('Vina 能预测酶活吗？', [], {})
        self.assertEqual(plan.action, 'reply')
        self.assertEqual(model.invoke.call_count, 1)

    def test_invalid_outputs_fail_closed_after_two_attempts(self):
        model = MagicMock()
        model.with_structured_output.return_value.invoke.return_value = None
        model.invoke.return_value = SimpleNamespace(content='{"action":"delete","reason":"invalid"}')
        with patch.object(planning, 'configured_model', return_value=model), self.assertRaises(planning.TurnPlanningError):
            planning.plan_turn('执行', [], {})
        self.assertEqual(model.invoke.call_count, 1)

    def test_valid_provider_response_requires_no_second_call(self):
        model = MagicMock()
        model.with_structured_output.return_value.invoke.return_value = planning.TurnPlan(action='inspect', reason='查询状态')
        with patch.object(planning, 'configured_model', return_value=model):
            self.assertEqual(planning.plan_turn('状态', [], {}).action, 'inspect')
        model.invoke.assert_not_called()


class ReliabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'state.sqlite3'
        self.store = api.RunStore(self.path)
        self.addCleanup(patch.stopall)
        patch.object(api, 'store', self.store).start()

    def domain_run(self):
        record = self.store.create('查结构域', 'general')
        api._save_execution_spec(record['run_id'], {'workflow': 'protein_domains', 'input': {'query_fasta': '>q\n' + 'A'*382, 'source':'regression'}})
        self.store.update(record['run_id'], status='running')
        return record

    def test_restart_marks_local_orphan_stopped_and_continuation_preserves_input(self):
        original = self.domain_run()
        reopened = api.RunStore(self.path)
        # Opening the store is not itself a restart of the coordinator.
        self.assertEqual(reopened.get(original['run_id'])['status'], 'running')
        records = reopened.claim_orphans()
        self.assertEqual(len(records), 1)
        self.assertTrue(reconcile(reopened, records[0], attempts=1))
        recovered = reopened.get(original['run_id'])
        self.assertEqual(recovered['status'], 'stopped')
        self.assertTrue(reopened.can_continue(recovered))
        with patch.object(api, 'store', reopened), patch.object(api, '_submit') as submit:
            created = api.continue_run(original['run_id'])
        self.assertEqual(submit.call_args.args[2]['spec']['input']['query_fasta'], '>q\n'+'A'*382)
        self.assertEqual(created.conversation_id, original['conversation_id'])

    def test_unreachable_worker_never_unlocks_duplicate_computation(self):
        record = self.store.create('结构预测', 'manual_openfold', manual={'accession':'P00782'})
        orphan = self.store.claim_orphans()[0]
        client = MagicMock();client.open.side_effect = TimeoutError()
        with patch.dict(os.environ, {'OPENFOLD_WORKER_URL':'http://unreachable'}):
            self.assertFalse(reconcile(self.store, orphan, attempts=1, client=client))
        row = self.store.get(record['run_id'])
        self.assertEqual(row['current_stage'], 'recovery_pending')
        self.assertFalse(self.store.can_continue(row))

    def test_cancel_ack_alone_is_not_worker_quiescence(self):
        record = self.store.create('结构预测', 'manual_openfold', manual={'accession':'P00782'})
        orphan = self.store.claim_orphans()[0]
        class Client:
            active = True
            def open(self, request, timeout):
                payload = {'status':'cancelling'} if request.data else {'job_id':record['run_id'], 'state':'cancelling' if self.active else 'stopped', 'active':self.active}
                return io.BytesIO(json.dumps(payload).encode())
        client = Client()
        with patch.dict(os.environ, {'OPENFOLD_WORKER_URL':'http://worker'}):
            self.assertFalse(reconcile(self.store, orphan, attempts=1, client=client))
            client.active = False
            self.assertTrue(reconcile(self.store, orphan, attempts=1, client=client))
        self.assertEqual(self.store.get(record['run_id'])['status'], 'stopped')

    def test_two_tool_operations_remain_visible_after_reopening_store(self):
        record = self.store.create('查询两条序列', 'general')
        for i, seq in enumerate(['A'*382, 'C'*241]):
            api._capture_execution_calls(record['run_id'], {'tool_calls':[{'id':f'call-{i}', 'name':'annotate_protein_domains', 'args':{'query_fasta':f'>q{i}\n{seq}', 'source':'regression'}}]})
            api._capture_operation_result(record['run_id'], {'type':'tool','tool_call_id':f'call-{i}','content':json.dumps({'status':'partial' if i==0 else 'succeeded','record_count':1})})
        self.store.update(record['run_id'], status='partial')
        current = self.store.create('做完了吗', 'general', conversation_id=record['conversation_id'])
        reopened = api.RunStore(self.path)
        facts = ConversationTasks(reopened, current['run_id']).snapshot()['tasks'][0]
        self.assertEqual([o['status'] for o in facts['operations']], ['partial','succeeded'])
        self.assertEqual(facts['operations'][0]['execution']['input_sequences'][0]['residue_count'], 382)
        api._reply_task_status(current['run_id'], ConversationTasks(reopened,current['run_id']), time.perf_counter())
        text = reopened.get(current['run_id'])['result']['final_text']
        self.assertNotIn('没有找到', text)
        self.assertIn('多个计算输入', text)
        self.assertFalse(reopened.can_continue(reopened.get(record['run_id'])))
        length_turn = reopened.create('第一条长度', 'general', conversation_id=record['conversation_id'])
        api._reply_task_status(length_turn['run_id'], ConversationTasks(reopened,length_turn['run_id']), time.perf_counter(), inspection='sequence_length', sequence_index=1)
        self.assertEqual(reopened.get(length_turn['run_id'])['result']['final_text'], '382')

    def test_latest_user_message_cannot_leak_into_earlier_queued_turn(self):
        first = self.store.create('先讨论', 'general')
        self.store.create('后面的指令', 'general', conversation_id=first['conversation_id'])
        context = self.store.conversation_context(first['conversation_id'], through_run_id=first['run_id'])
        self.assertEqual([m['content'] for m in context], ['先讨论'])

    def test_chat_finishes_while_science_pool_is_occupied(self):
        release = threading.Event();occupied = threading.Event()
        agent = MagicMock()
        agent.stream.return_value = iter([{'messages':[{'type':'ai','content':'你好'}]}])
        with ThreadPoolExecutor(max_workers=1) as science, ThreadPoolExecutor(max_workers=2) as dialogue:
            with patch.object(api,'executor',science), patch.object(api,'control_executor',dialogue), patch.object(api,'plan_turn',return_value=planning.TurnPlan(action='reply',reason='问候')), patch.object(api,'build_agent',return_value=agent):
                def block(_):
                    occupied.set();release.wait(5)
                job = self.domain_run();api._submit(job['run_id'],block)
                self.assertTrue(occupied.wait(1))
                try:
                    chat = api.create_run(api.RunRequest(prompt='你好'))
                    deadline = time.monotonic()+2
                    while time.monotonic()<deadline and self.store.get(chat.run_id)['status'] not in api.TERMINAL:
                        time.sleep(.01)
                    self.assertEqual(self.store.get(chat.run_id)['status'],'succeeded')
                    self.assertFalse(release.is_set())
                finally:
                    release.set()

    def test_execute_handoff_preserves_new_future_for_cancellation(self):
        release = threading.Event();queued = threading.Event()
        with ThreadPoolExecutor(max_workers=1) as science, ThreadPoolExecutor(max_workers=2) as dialogue:
            with patch.object(api,'executor',science), patch.object(api,'control_executor',dialogue), patch.object(api,'plan_turn',return_value=planning.TurnPlan(action='execute',reason='计算')):
                science.submit(release.wait,5)
                try:
                    chat = api.create_run(api.RunRequest(prompt='分析这个蛋白质'))
                    deadline=time.monotonic()+2
                    while time.monotonic()<deadline and self.store.get(chat.run_id)['current_stage']!='scientific_queue':
                        time.sleep(.01)
                    self.assertEqual(self.store.get(chat.run_id)['current_stage'],'scientific_queue')
                    time.sleep(.05)
                    with api.future_lock:
                        self.assertIn(chat.run_id, api.run_futures)
                    self.assertEqual(api.cancel_run(chat.run_id)['status'],'stopped')
                finally:
                    release.set()


class WorkerStatusTests(unittest.TestCase):
    def test_cancelled_active_job_is_not_reported_stopped_until_scope_exits(self):
        job=uuid.uuid4().hex
        self.assertEqual(job_status(job)['state'],'unknown')
        with job_scope(job) as control:
            control.cancel()
            self.assertEqual(job_status(job)['state'],'cancelling')
            self.assertTrue(job_status(job)['active'])
        self.assertEqual(job_status(job)['state'],'stopped')
        self.assertFalse(job_status(job)['active'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
