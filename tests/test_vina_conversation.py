"""Durable Vina task association and worker progress, isolated from live state."""
import base64
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

with tempfile.TemporaryDirectory() as initial_db:
    with patch.dict(os.environ, {'MARINE_AGENT_STATE_DB': initial_db + '/state.sqlite3'}):
        from marine_agent import api
from marine_agent import vina_tools

ROOT = Path(__file__).resolve().parents[1]


class VinaConversationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = api.RunStore(Path(self.temp.name) / 'state.sqlite3')
        self.addCleanup(patch.stopall)
        patch.object(api, 'store', self.store).start()
        patch.object(api, 'VINA_CALLBACK_TOKEN', 'test-only-callback-secret').start()
        self.executor = patch.object(api.executor, 'submit').start()
        self.payload = {
            'receptor': {'name': '1IEP.pdbqt', 'content_base64': base64.b64encode((ROOT / 'results/1iep/receptor.pdbqt').read_bytes()).decode()},
            'ligands': [{'id': 'imatinib', 'name': 'imatinib.pdbqt', 'content_base64': base64.b64encode((ROOT / 'results/1iep/ligand.pdbqt').read_bytes()).decode()}],
            'box': {'center': [15.614, 53.38, 15.455], 'size': [30, 30, 30]}, 'params': {'exhaustiveness': 1},
        }

    def test_new_vina_is_immediately_visible_in_history_without_restart(self):
        created = api.create_vina_run(api.VinaRunRequest(**self.payload))
        conversation = self.store.get_conversation(created.conversation_id)
        self.assertEqual(conversation['runs'][0]['id'], created.run_id)
        self.assertEqual(conversation['messages'][0]['run_id'], created.run_id)
        self.assertEqual(self.store.list_conversations()[0]['latest_run_status'], 'queued')

    def test_existing_conversation_is_reused_and_invalid_conversation_is_rejected(self):
        first = self.store.create('已有结构预测', 'manual_openfold')
        self.store.update(first['run_id'], status='succeeded')
        created = api.create_vina_run(api.VinaRunRequest(**self.payload, conversation_id=first['conversation_id']))
        self.assertEqual(created.conversation_id, first['conversation_id'])
        self.assertEqual(len(self.store.get_conversation(created.conversation_id)['runs']), 2)
        with self.assertRaises(api.HTTPException) as error:
            api.create_vina_run(api.VinaRunRequest(**self.payload, conversation_id='missing-conversation'))
        self.assertEqual(error.exception.status_code, 404)
        self.assertEqual(self.executor.call_count, 1)

    def test_partial_and_failed_runs_persist_their_reports_in_the_same_conversation(self):
        for worker_status in ['partial', 'failed']:
            created = api.create_vina_run(api.VinaRunRequest(**self.payload))
            result = {'status': worker_status, 'results': [{'id': 'one', 'status': 'succeeded'}, {'id': 'two', 'status': 'failed'}], 'ranking': [], 'elapsed_seconds': 1}
            with patch.object(api, 'VINA_WORKER_URL', 'unit-test-worker'), patch.object(api, 'run_vina_worker', return_value=result):
                api._execute_vina_run(created.run_id, self.payload)
            conversation = self.store.get_conversation(created.conversation_id)
            self.assertEqual(conversation['runs'][0]['status'], worker_status)
            self.assertEqual([m['role'] for m in conversation['messages']], ['user', 'assistant'])
            self.assertEqual(conversation['messages'][-1]['run_id'], created.run_id)

    def test_callback_counts_are_preserved_and_late_callbacks_do_not_reopen_finished_jobs(self):
        created = api.create_vina_run(api.VinaRunRequest(**self.payload))
        event = {'stage': '对接计算', 'name': '配体完成', 'ligand_id': 'one', 'status': 'failed', 'completed': 1, 'total': 2, 'failed': 1, 'succeeded': 0}
        api.vina_progress(created.run_id, event, api.VINA_CALLBACK_TOKEN)
        self.assertEqual(self.store.events(created.run_id)[-1]['payload']['completed'], 1)
        self.assertEqual(self.store.events(created.run_id)[-1]['payload']['status'], 'failed')
        self.store.update(created.run_id, status='partial', current_stage='completed')
        before = len(self.store.events(created.run_id))
        self.assertEqual(api.vina_progress(created.run_id, event, api.VINA_CALLBACK_TOKEN)['status'], 'ignored_terminal')
        self.assertEqual(len(self.store.events(created.run_id)), before)
        self.assertEqual(self.store.get(created.run_id)['current_stage'], 'completed')

    def test_real_vina_emits_completed_counts_including_a_failed_ligand(self):
        binary = ROOT / '.tools/bin/vina'
        if not binary.is_file(): binary = Path('/Users/yili/.marine-agent/runtime/.tools/bin/vina')
        self.assertTrue(binary.is_file())
        config = vina_tools.validate_config(self.payload)
        config['ligands'].append({'id': 'bad-input', 'name': 'bad.pdbqt', 'bytes': b'ROOT\nATOM  INVALID_FORMAT\nENDROOT\nTORSDOF 0\n'})
        config['batch_workers'] = 2
        events = []
        with tempfile.TemporaryDirectory(dir=ROOT / 'results') as directory, patch.object(vina_tools, 'VINA_BINARY', str(binary)):
            result = vina_tools.run_batch(Path(directory), config, emit=events.append)
        self.assertEqual(result['status'], 'partial')
        completed = [event for event in events if 'completed' in event]
        self.assertEqual([event['completed'] for event in completed], [1, 2])
        self.assertEqual(completed[-1]['failed'], 1)
        self.assertEqual(completed[-1]['succeeded'], 1)
        self.assertEqual({event['status'] for event in completed}, {'succeeded', 'failed'})


if __name__ == '__main__': unittest.main(verbosity=2)
