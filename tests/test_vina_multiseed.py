"""Real multi-seed Vina execution and per-seed continuation contracts."""
import base64
import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tests.test_vina_conversation import api, ROOT
from marine_agent import vina_tools


class VinaMultiSeedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / "results")
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.store = api.RunStore(self.path / "state.sqlite3")
        self.addCleanup(patch.stopall)
        patch.object(api, "store", self.store).start()
        patch.object(api, "_submit").start()
        reference = ROOT / "results/vina-runs/fd3cad6e-d4c5-4777-bac5-38068fb4ea05"
        def encoded(path):
            return base64.b64encode(path.read_bytes()).decode()
        self.request = {
            "receptor": {"name": "receptor.pdbqt", "content_base64": encoded(reference / "receptor.pdbqt")},
            "ligands": [{"id": "caffeine", "name": "caffeine.pdbqt",
                         "content_base64": encoded(reference / "ligands/library-L001-primary-1.pdbqt")}],
            "box": {"center": [15.614, 53.38, 15.455], "size": [30, 30, 30]},
            "params": {"exhaustiveness": 1, "num_modes": 3, "cpu": 1,
                       "seed": 42, "seeds": [42, 2026, 20261021]},
        }
        self.binary = ROOT / ".tools/bin/vina"
        if not self.binary.is_file():
            self.binary = Path("/Users/yili/.marine-agent/runtime/.tools/bin/vina")

    def test_seed_validation_rejects_duplicate_and_out_of_range_values(self):
        for seeds in ([], [1, 1], [1, -1], [1, "2"], [1, True]):
            request = copy.deepcopy(self.request)
            request["params"]["seeds"] = seeds
            with self.assertRaises(vina_tools.VinaInputError):
                vina_tools.validate_config(request)

    def test_real_three_seed_run_has_independent_commands_scores_and_artifacts(self):
        config = vina_tools.validate_config(self.request)
        with patch.object(vina_tools, "ROOT", self.path), patch.object(vina_tools, "VINA_BINARY", str(self.binary)):
            result = vina_tools.run_batch(self.path / "results/vina-runs/three", config)
        self.assertEqual(result["status"], "succeeded")
        self.assertEqual(result["ligand_count"], 1)
        self.assertEqual(result["task_count"], 3)
        self.assertEqual([row["seed"] for row in result["results"]], [42, 2026, 20261021])
        self.assertEqual(len({row["id"] for row in result["results"]}), 3)
        consistency = result["top_pose_consistency"][0]
        self.assertEqual(consistency["status"], "available")
        self.assertEqual(consistency["seed_count"], 3)
        self.assertEqual(len(consistency["pairs"]), 3)
        self.assertIn(consistency["central_top_pose_seed"], [42, 2026, 20261021])
        self.assertGreaterEqual(consistency["max_top_pose_rmsd_angstrom"], 0)
        self.assertIn("no ligand-only fitting", consistency["method"])
        self.assertIn(consistency["review_status"], {"consistent_review_required", "divergent_review_required"})
        self.assertIn(consistency["recommended_review_task_id"], {row["id"] for row in result["results"]})
        for row in result["results"]:
            self.assertEqual(int(row["command"][row["command"].index("--seed") + 1]), row["seed"])
            self.assertTrue((self.path / row["output"]).is_file())
            self.assertTrue(1 <= len(row["scores"]) <= 3)

    def test_interrupted_run_reuses_finished_seed_and_recomputes_only_missing_seeds(self):
        config = vina_tools.validate_config(self.request)
        source_id = "source-multiseed"
        complete_dir = self.path / "results/vina-runs" / source_id
        with patch.object(vina_tools, "ROOT", self.path), patch.object(vina_tools, "VINA_BINARY", str(self.binary)):
            complete = vina_tools.run_batch(complete_dir, config)
        source = self.store.create("三种子对接", "vina", manual=self.request)
        # Keep the real run ID aligned with the artifact directory.
        source_id = source["run_id"]
        source_dir = self.path / "results/vina-runs" / source_id
        source_dir.mkdir(parents=True)
        for name in ["receptor.pdbqt", "config.json"]:
            (source_dir / name).write_bytes((complete_dir / name).read_bytes())
        (source_dir / "ligands").mkdir(exist_ok=True)
        first = copy.deepcopy(complete["results"][0])
        first["output"] = f"results/vina-runs/{source_id}/ligands/{first['id']}.out.pdbqt"
        first["log"] = f"results/vina-runs/{source_id}/ligands/{first['id']}.log"
        for key in ["output", "log"]:
            (self.path / first[key]).write_bytes((complete_dir / "ligands" / Path(complete["results"][0][key]).name).read_bytes())
        partial = {**complete, "status": "stopped", "results": [first], "ranking": []}
        self.store.update(source_id, result_json=json.dumps({"vina": partial}), status="stopped")
        continued = api.continue_run(source_id)
        # _submit is patched in setUp; obtain the actual mock from the API patch.
        handoff = api._submit.call_args.args[2]
        real = vina_tools.run_batch
        with patch.object(api, "ROOT", self.path), patch.object(vina_tools, "ROOT", self.path), \
             patch.object(api, "VINA_WORKER_URL", ""), patch.object(vina_tools, "VINA_BINARY", str(self.binary)):
            with patch.object(api, "run_batch", wraps=real) as execute:
                api._controlled_execute(continued.run_id, api._execute_continuation, handoff)
        pending = execute.call_args.args[1]
        self.assertEqual(pending["task_seeds"], {"caffeine": [2026, 20261021]})
        saved = self.store.get(continued.run_id)
        self.assertEqual(saved["status"], "succeeded", saved.get("error"))
        self.assertEqual(saved["result"]["vina"]["task_count"], 3)
        self.assertEqual([row["seed"] for row in saved["result"]["vina"]["results"]], [42, 2026, 20261021])
        self.assertEqual(saved["result"]["vina"]["continuation"]["retained_tasks"], [first["id"]])
        consistency = saved["result"]["vina"]["top_pose_consistency"][0]
        self.assertEqual(consistency["status"], "available")
        self.assertEqual(consistency["seed_count"], 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
