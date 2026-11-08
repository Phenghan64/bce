"""New scientific inputs, real BLAST and failure recovery, isolated from live state."""
import copy
import json
import os
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from marine_agent import protein_screening as screen
from marine_agent.execution_control import RunCancelled, job_scope
from tests.test_vina_conversation import api
from marine_agent.task_lifecycle import ConversationTasks, saved_spec

FIXTURES = Path(__file__).parent / "fixtures/protein-screening"


def request():
    return {"query_fasta": "".join((FIXTURES / f"{a}.faa").read_text() for a in ("P00782", "P0A7V0", "P69905")),
            "reference_fasta": (FIXTURES / "P00780.faa").read_text(),
            "query_source": "UniProt P00782, P0A7V0, P69905; sources.json hashes",
            "reference_source": "UniProt P00780; sources.json hash",
            "objective": "Find subtilisin homologs in new supplied protein records", "threads": 2}


class ProteinScreenTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_screen(self, payload=None, **kwargs):
        return screen.run_screen(payload or request(), root=self.root, **kwargs)

    def directory(self, result):
        return self.root / "results/protein-screenings" / result["analysis_id"]

    def test_real_unseen_inputs_find_homolog_and_exclude_unrelated_proteins(self):
        report = self.run_screen()
        self.assertEqual(report["query_count"], 3)
        self.assertEqual(report["candidate_count"], 1)
        candidate = report["candidates"][0]
        self.assertEqual(candidate["protein_id"], "sp|P00782|SUBT_BACAM")
        self.assertFalse(candidate["exact_reference_sequence"])
        self.assertGreater(candidate["identity_percent"], 50)
        self.assertFalse(report["quality_gate"]["activity_confirmed"])
        self.assertEqual(report["quality_gate"]["domain_verification"], "not_run")
        selected = screen.read_analysis(report["analysis_id"], root=self.root, protein_id=candidate["protein_id"])
        self.assertEqual(len(selected["selected_candidate"]["sequence"]), candidate["length"])

    def test_same_inputs_reuse_and_corrupt_hit_output_rebuilds_only_affected_stages(self):
        first = self.run_screen()
        second = self.run_screen()
        self.assertEqual(first["analysis_id"], second["analysis_id"])
        self.assertEqual(second["execution"]["reused_stages"], ["validate_inputs", "makeblastdb", "blastp", "rank_and_report"])
        (self.directory(first) / "hits.tsv").write_text("corrupt\n")
        self.assertEqual(screen.read_analysis(first["analysis_id"], root=self.root)["status"], "invalid_artifacts")
        repaired = self.run_screen()
        self.assertEqual(repaired["execution"]["reused_stages"], ["validate_inputs", "makeblastdb"])
        self.assertEqual(repaired["candidates"], first["candidates"])

    def test_stop_before_blast_preserves_database_and_resumes(self):
        import uuid
        with job_scope(uuid.uuid4().hex) as control:
            def stop(event):
                if event == {"name": "blastp", "status": "running"}:
                    control.cancel()
            with self.assertRaises(RunCancelled):
                self.run_screen(emit=stop)
        manifests = list(self.root.glob("results/protein-screenings/*/manifest.json"))
        self.assertEqual(len(manifests), 1)
        interrupted = json.loads(manifests[0].read_text())
        self.assertEqual(interrupted["status"], "stopped")
        resumed = self.run_screen()
        self.assertEqual(resumed["execution"]["reused_stages"], ["validate_inputs", "makeblastdb"])
        self.assertEqual(resumed["candidate_count"], 1)

    def test_parameter_and_input_changes_create_independent_analyses(self):
        first = self.run_screen()
        modified = request()
        modified["min_identity"] = 100
        strict = self.run_screen(modified)
        self.assertEqual(strict["candidate_count"], 0)
        self.assertNotEqual(strict["analysis_id"], first["analysis_id"])
        modified = request()
        modified["query_fasta"] = (FIXTURES / "P69905.faa").read_text()
        negative = self.run_screen(modified)
        self.assertEqual(negative["candidate_count"], 0)
        self.assertNotEqual(negative["analysis_id"], first["analysis_id"])

    def test_malformed_or_ambiguous_inputs_are_not_silently_repaired(self):
        for fasta in (">same\nACDEF\n>same\nACDEG\n", ">empty\n", "ACDEFG", ">bad\nACD*EF\n"):
            invalid = request()
            invalid["query_fasta"] = fasta
            with self.assertRaises(ValueError):
                self.run_screen(invalid)
        for field, value in (("min_identity", float("nan")), ("threads", 0), ("query_source", " ")):
            invalid = request()
            invalid[field] = value
            with self.assertRaises(ValueError):
                self.run_screen(invalid)
        self.assertFalse((self.root / "results").exists())

    def test_tool_failure_does_not_become_a_successful_artifact_and_is_recoverable(self):
        real = screen._run
        def failed(command, **kwargs):
            if "-query" in command:
                raise OSError("injected execution-node interruption")
            return real(command, **kwargs)
        with patch.object(screen, "_run", side_effect=failed):
            with self.assertRaises(OSError):
                self.run_screen()
        manifest = json.loads(next(self.root.glob("results/protein-screenings/*/manifest.json")).read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual(manifest["stages"]["blastp"]["status"], "failed")
        recovered = self.run_screen()
        self.assertEqual(recovered["candidate_count"], 1)
        self.assertEqual(recovered["execution"]["reused_stages"], ["validate_inputs", "makeblastdb"])


class ProteinAPITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = api.RunStore(Path(self.temp.name) / "state.sqlite3")
        self.addCleanup(patch.stopall)
        patch.object(api, "store", self.store).start()
        self.submit = patch.object(api, "_submit").start()

    def test_durable_submission_continuation_and_exact_input_preservation(self):
        created = api.create_protein_run(api.ProteinScreenRunRequest(**request()))
        record = self.store.get(created.run_id)
        self.assertEqual(saved_spec(record)["input"]["query_fasta"], request()["query_fasta"])
        self.store.update(created.run_id, status="stopped")
        continued = api.continue_run(created.run_id)
        handoff = self.submit.call_args.args[2]
        with patch.object(api, "_execute_protein_spec") as execute, patch.object(api, "_execute_openfold_spec") as fold:
            api._execute_continuation(continued.run_id, handoff)
        execute.assert_called_once()
        fold.assert_not_called()
        self.assertEqual(execute.call_args.args[1]["input"], record["execution_spec"]["input"])

    def test_failed_backend_has_failed_ledger_and_can_continue(self):
        created = api.create_protein_run(api.ProteinScreenRunRequest(**request()))
        spec = saved_spec(self.store.get(created.run_id))
        with patch.object(api, "execute_screen", side_effect=RuntimeError("worker unavailable")):
            api._execute_protein_spec(created.run_id, spec)
        failed = self.store.get(created.run_id)
        self.assertEqual(failed["status"], "failed")
        self.assertIsNotNone(saved_spec(failed))
        api.continue_run(created.run_id)

    def test_agent_capture_and_progress_share_the_same_task_state(self):
        current = self.store.create("筛选新蛋白", "general")
        api._capture_execution_calls(current["run_id"], {"tool_calls": [{"name": "run_protein_homology_screen", "args": request()}]})
        spec = saved_spec(self.store.get(current["run_id"]))
        self.assertEqual(spec["workflow"], "protein_screen")
        api._protein_progress(current["run_id"], {"name": "blastp", "status": "running"})
        self.assertEqual(self.store.get(current["run_id"])["current_stage"], "blastp")
        self.assertIsNone(api._inline_structure_sequence("预测候选酶\n" + request()["query_fasta"]))

    def test_new_research_tools_cannot_select_frozen_benchmark_and_intent_uses_new_inputs(self):
        from marine_agent.agent import TOOLS
        from marine_agent.intent import classify_intent
        names = {tool.name for tool in TOOLS}
        self.assertIn("run_protein_homology_screen", names)
        self.assertIn("run_openfold_sequence", names)
        self.assertNotIn("inspect_workflow", names)
        self.assertNotIn("run_collagenase_batch_screen", names)
        self.assertEqual(api._resolve_mode("寻找潜在胶原蛋白酶", "auto"), "general")
        self.assertEqual(api._resolve_mode("复查 collagenase-agent-20", "auto"), "collagenase")
        intent = classify_intent("请对以下蛋白做同源筛选，暂不预测结构\n" + request()["query_fasta"])
        self.assertEqual(intent["primary_intent"], "homology_screening")


class ProteinTransportTests(unittest.TestCase):
    def test_authenticated_worker_roundtrip_returns_verified_real_blast_artifacts(self):
        from marine_agent import vina_worker, vina_remote, protein_tools
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            server = vina_worker.ReusableThreadingHTTPServer(("127.0.0.1", 0), vina_worker.Handler)
            server.worker_token = "isolated-test-token"
            server.job_timeout = 60
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            try:
                with patch.object(screen, "ROOT", directory), patch.object(vina_worker, "ROOT", directory), patch.object(vina_remote, "ROOT", directory), patch.object(protein_tools, "ROOT", directory), patch.dict(os.environ, {
                    "PROTEIN_BLAST_BIN": str(FIXTURES.parents[2] / ".tools/blast-2.17.0/bin"),
                    "PROTEIN_WORKER_URL": f"http://127.0.0.1:{server.server_port}",
                    "PROTEIN_WORKER_TOKEN": "isolated-test-token", "VINA_CALLBACK_URL": "",
                }):
                    result = protein_tools.execute_screen(request())
                    self.assertEqual(result["status"], "succeeded")
                    self.assertEqual(result["candidate_count"], 1)
                    self.assertTrue((directory / result["artifacts"]["manifest.json"]).is_file())
            finally:
                server.shutdown()
                server.server_close()
                thread.join(timeout=5)


if __name__ == "__main__":
    unittest.main(verbosity=2)
