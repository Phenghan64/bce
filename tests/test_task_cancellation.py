"""Cancellation uses real owned processes; application state tests use isolated SQLite."""
import base64
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
import uuid
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

from marine_agent.execution_control import JobControl, RunCancelled, control_for, job_scope
from tests.test_vina_conversation import api, ROOT
from marine_agent.vina_tools import validate_config, run_batch


class ProcessTests(unittest.TestCase):
    def test_stop_reaps_owned_process_group_without_touching_another_job(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'child.pid'
            script = ('import subprocess,sys,time; p=subprocess.Popen([sys.executable,"-c",'
                      '"import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)"]);'
                      'open(sys.argv[1],"w").write(str(p.pid)); time.sleep(60)')
            control = JobControl(uuid.uuid4().hex)
            stopped = threading.Event()
            def run():
                try: control.run([sys.executable,'-c',script,str(path)], capture_output=True, text=True)
                except RunCancelled: stopped.set()
            other = subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'])
            try:
                thread = threading.Thread(target=run);thread.start()
                deadline=time.monotonic()+5
                while not path.exists() and time.monotonic()<deadline: time.sleep(.02)
                self.assertTrue(path.exists())
                child=int(path.read_text())
                control.cancel();thread.join(6)
                self.assertFalse(thread.is_alive());self.assertTrue(stopped.is_set())
                self.assertIsNone(other.poll())
                state=subprocess.run(['/bin/ps','-o','stat=','-p',str(child)],capture_output=True,text=True).stdout.strip()
                self.assertTrue(not state or state.startswith('Z'),state)
                self.assertEqual(len(control.processes),0)
            finally:
                other.terminate();other.wait()

    def test_cancel_before_submission_never_launches_process(self):
        control=control_for(uuid.uuid4().hex);control.cancel()
        with self.assertRaises(RunCancelled), patch('subprocess.Popen') as launch:
            control.run([sys.executable,'-c','raise Exception()'])
        launch.assert_not_called()

    def test_timeout_also_cleans_up_group(self):
        control=JobControl(uuid.uuid4().hex)
        with self.assertRaises(subprocess.TimeoutExpired):
            control.run([sys.executable,'-c','import time;time.sleep(60)'],timeout=.1,capture_output=True,text=True)
        self.assertEqual(len(control.processes),0)

    def test_real_vina_batch_stop_preserves_finished_score_and_marks_unstarted_tasks(self):
        source=ROOT/'results/vina-runs/fd3cad6e-d4c5-4777-bac5-38068fb4ea05'
        request={'receptor':{'name':'receptor.pdbqt','content_base64':base64.b64encode((source/'receptor.pdbqt').read_bytes()).decode()},
                 'ligands':[{'id':name,'name':name+'.pdbqt','content_base64':base64.b64encode((source/'ligands'/file).read_bytes()).decode()}
                            for name,file in [('fast','library-L001-primary-1.pdbqt'),('slow','imatinib-control.pdbqt'),('pending','imatinib-control.pdbqt')]],
                 'box':{'center':[15.614,53.38,15.455],'size':[30,30,30]},'params':{'exhaustiveness':32}}
        control=JobControl(uuid.uuid4().hex);events=[]
        def emit(event):
            events.append(event)
            if event.get('succeeded')==1: control.cancel()
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory, patch('marine_agent.vina_tools.VINA_BINARY',str(ROOT/'.tools/bin/vina')):
            config=validate_config(request);config['batch_workers']=1
            result=run_batch(Path(directory),config,emit=emit,control=control)
            self.assertEqual(result['status'],'stopped')
            self.assertEqual(result['results'][0]['status'],'succeeded')
            self.assertTrue(result['results'][0]['scores'])
            self.assertTrue(all(row['status'] in {'stopped','skipped'} for row in result['results'][1:]))
            self.assertEqual(len(result['ranking']),1)
            self.assertTrue(all(event.get('failed',0)==0 for event in events))
            self.assertEqual(len(control.processes),0)


class StateTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.store=api.RunStore(Path(self.temp.name)/'state.sqlite3')
        self.addCleanup(patch.stopall)
        patch.object(api,'store',self.store).start()
        patch.object(api.cancel_executor,'submit').start()

    def test_completed_task_cannot_be_cancelled_or_rewritten(self):
        record=self.store.create('complete','vina');run=record['run_id']
        self.store.update(run,status='succeeded',result_json=json.dumps({'final_text':'done'}))
        self.assertEqual(api.cancel_run(run)['status'],'succeeded')
        self.assertEqual(self.store.get(run)['result']['final_text'],'done')
        self.assertFalse(control_for(run).cancelled)

    def test_stop_is_idempotent_and_late_completion_cannot_overwrite_it(self):
        record=self.store.create('running','manual_openfold');run=record['run_id']
        self.store.update(run,status='running')
        self.assertEqual(api.cancel_run(run)['status'],'cancelling')
        self.assertEqual(api.cancel_run(run)['status'],'cancelling')
        self.store.update(run,status='succeeded',result_json='{}')
        self.store.append_message(record['conversation_id'],'assistant','late report',run)
        self.assertEqual(self.store.get(run)['status'],'cancelling')
        api._finish_stopped(run);api._finish_stopped(run)
        self.store.update(run,status='failed')
        messages=self.store.get_conversation(record['conversation_id'])['messages']
        self.assertEqual([m['role'] for m in messages],['user','assistant'])
        self.assertEqual(self.store.get(run)['status'],'stopped')

    def test_queued_job_is_removed_and_never_executed(self):
        blocker=threading.Event();called=threading.Event()
        with ThreadPoolExecutor(max_workers=1) as pool, patch.object(api,'executor',pool):
            pool.submit(blocker.wait)
            record=self.store.create('queued','vina');run=record['run_id']
            try:
                api._submit(run,lambda _:called.set())
                self.assertEqual(api.cancel_run(run)['status'],'stopped')
            finally: blocker.set()
        self.assertFalse(called.is_set())
        self.assertEqual(self.store.get(run)['status'],'stopped')

    def test_running_wrapper_waits_for_real_exit_then_records_stop(self):
        record=self.store.create('running','manual_openfold');run=record['run_id']
        launched=threading.Event()
        def execute(_):
            self.store.update(run,status='running')
            launched.set()
            control_for(run).run([sys.executable,'-c','import time;time.sleep(60)'])
        thread=threading.Thread(target=api._controlled_execute,args=(run,execute));thread.start()
        self.assertTrue(launched.wait(2));api.cancel_run(run);thread.join(5)
        self.assertFalse(thread.is_alive());self.assertEqual(self.store.get(run)['status'],'stopped')
        self.assertEqual(len(control_for(run).processes),0)


if __name__=='__main__':unittest.main()
