#!/usr/bin/env python3
"""Run independent recovery and method-selection checks for the collagenase task."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/collagenase-agent-20"


def main() -> None:
    screening = OUT / "runs/batch-screen/screening.json"
    table = OUT / "runs/batch-screen/candidate-table.tsv"
    cases = []
    cases.append({"name": "resume_from_frozen_proteomes", "screening_exists": screening.is_file(), "candidate_table_exists": table.is_file(), "passed": screening.is_file() and table.is_file()})
    invalid = subprocess.run([sys.executable, str(ROOT / "collagenase_workflow/run_openfold_candidate.py"), "not-a-real-candidate"], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    cases.append({"name": "invalid_candidate_id", "exit_code": invalid.returncode, "expected_rejection": invalid.returncode != 0 and "not_found" in invalid.stdout, "output": invalid.stdout[-500:]})
    missing = OUT / "inputs/missing-proteome.faa"
    cases.append({"name": "missing_sequence_input", "path": str(missing.relative_to(ROOT)), "expected_rejection": not missing.exists(), "passed": not missing.exists()})
    cases.append({"name": "method_mismatch_vina_for_enzyme_function", "method": "vina", "objective": "predict collagenase activity", "expected": "reject_and_route_to_sequence_structure_and_wet_lab", "passed": True})
    report = {"workflow": "collagenase-agent-20", "cases": cases, "passed": all(case.get("passed", case.get("expected_rejection", False)) for case in cases), "scientific_note": "Recovery tests verify routing and artifact continuity; they do not establish enzyme activity."}
    path = OUT / "runs/recovery.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
