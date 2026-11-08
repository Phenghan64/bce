"""Inline OpenFold requests preserve the complete user-declared scientific input."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tests.test_vina_conversation import api, ROOT


class OpenFoldRequestContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / "results")
        self.addCleanup(self.temp.cleanup)
        self.store = api.RunStore(Path(self.temp.name) / "state.sqlite3")
        self.addCleanup(patch.stopall)
        patch.object(api, "store", self.store).start()

    def test_explicit_smiles_is_normalized_for_openfold(self):
        ligands = api._inline_structure_ligands("预测复合物；甲苯 SMILES：Cc1ccccc1。")
        self.assertEqual(len(ligands), 1)
        self.assertEqual(ligands[0]["smiles"], "Cc1ccccc1")
        self.assertEqual(ligands[0]["chain_ids"], ["X"])

    def test_multiple_smiles_have_unique_non_protein_chains(self):
        ligands = api._inline_structure_ligands("SMILES: CCO\nSMILES=CCN")
        self.assertEqual([row["smiles"] for row in ligands], ["CCO", "CCN"])
        self.assertEqual([row["chain_ids"][0] for row in ligands], ["X", "Y"])

    def test_explicit_ccd_is_preserved(self):
        ligands = api._inline_structure_ligands("复合物配体 CCD：MBN")
        self.assertEqual(ligands[0]["ccd_codes"], ["MBN"])
        self.assertEqual(ligands[0]["chain_ids"], ["X"])

    def test_unlabelled_chemical_name_is_not_guessed_as_smiles(self):
        self.assertEqual(api._inline_structure_ligands("预测蛋白和甲苯复合物"), [])

    def test_inline_shortcut_preserves_ligands_and_records_vina_followup(self):
        record = self.store.create("复合物预测并对接", "general")
        prompt = "用OpenFold3预测复合物，再用Vina对接。SMILES：Cc1ccccc1\n" + "A" * 30
        with patch.object(api, "_execute_openfold_spec") as execute, patch.object(api, "plan_turn") as planner:
            api._execute_run(record["run_id"], prompt, "general")
        planner.assert_not_called()
        spec = execute.call_args.args[1]
        self.assertEqual(spec["input"]["sequence"], "A" * 30)
        self.assertEqual(spec["input"]["ligands"][0]["smiles"], "Cc1ccccc1")
        self.assertEqual(spec["requested_followup"], "vina")

    def test_summary_distinguishes_protein_only_from_complex(self):
        base = {"status": "succeeded", "sequence": {"length": 30},
                "confidence": {"avg_plddt": 90}, "quality_gate": {"decision": "ok"}}
        single = api._manual_openfold_summary({"source_kind": "inline_sequence", "result": {**base, "input": {"ligands": []}}})
        self.assertIn("蛋白结构预测完成", single)
        self.assertIn("无，本次为单蛋白结构预测", single)
        complex_result = {**base, "input": {"ligands": [{"smiles": "Cc1ccccc1", "chain_ids": ["X"]}]}}
        complex_text = api._manual_openfold_summary({"source_kind": "inline_sequence", "result": complex_result})
        self.assertIn("复合物结构预测完成", complex_text)
        self.assertIn("Cc1ccccc1（链 X）", complex_text)

    def test_requested_vina_is_explicitly_pending_not_reported_as_complete(self):
        record = self.store.create("结构并对接", "general")
        ligands = api._inline_structure_ligands("SMILES: Cc1ccccc1")
        payload = {"status": "succeeded", "result": {"status": "succeeded", "sequence": {"length": 30},
                   "confidence": {"avg_plddt": 90}, "quality_gate": {"decision": "ok"},
                   "input": {"ligands": ligands}}}
        spec = {"workflow": "openfold_sequence", "input": {"sequence": "A" * 30, "ligands": ligands},
                "requested_followup": "vina"}
        with patch.object(api, "run_openfold_sequence") as tool:
            tool.invoke.return_value = json.dumps(payload)
            api._execute_openfold_spec(record["run_id"], spec, 0.0)
        saved = self.store.get(record["run_id"])
        self.assertEqual(saved["status"], "succeeded")
        self.assertIn("后续 Vina：尚未运行", saved["result"]["final_text"])
        self.assertEqual(saved["result"]["pending_actions"][0]["status"], "awaiting_confirmation")


if __name__ == "__main__":
    unittest.main(verbosity=2)
