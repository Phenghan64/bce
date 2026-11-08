"""Typed, auditable tools exposed to the Deep Agent.

These tools intentionally do not expose an arbitrary shell or filesystem API.
Heavy scientific stages remain deterministic workflow jobs; the agent only
inspects artifacts, runs approved validation jobs, and explains evidence.
"""
from __future__ import annotations

from .execution_control import current_control

import json
import csv
import fcntl
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Literal

from langchain_core.tools import tool

from .intent_tools import INTENT_TOOLS
from .openfold_tools import OPENFOLD_TOOLS
from .protein_tools import PROTEIN_TOOLS


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/marine-microbe-workflow"
REPORT_JSON = RESULTS / "REPORT.json"
AUDIT_JSON = RESULTS / "runs/candidate-audit-v2/candidate-audit.json"
IPR_JSON = RESULTS / "runs/interpro-candidate-panel/interpro-panel.json"
RECOVERY_JSON = RESULTS / "runs/recovery-tests/recovery-tests.json"
STAGE_LOCK = RESULTS / "runs/.marine-agent-stage.lock"


def _read_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"required artifact is missing: {path}")
    return json.loads(path.read_text())


@contextmanager
def _stage_lock():
    """Serialize workflow mutations when the model emits parallel tool calls."""
    STAGE_LOCK.parent.mkdir(parents=True, exist_ok=True)
    with STAGE_LOCK.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _run_approved(script: str, timeout: int = 300) -> dict:
    command = [sys.executable, str(ROOT / "marine_microbe_workflow" / script)]
    started = time.monotonic()
    try:
        completed = (current_control().run if current_control() else subprocess.run)(
            command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
        return {
            "status": "succeeded" if completed.returncode == 0 else "failed",
            "command": command,
            "exit_code": completed.returncode,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "output_tail": completed.stdout[-2000:],
            "ok": completed.returncode == 0,
        }
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return {
            "status": "timed_out",
            "command": command,
            "exit_code": None,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "output_tail": output[-2000:],
            "timeout_seconds": timeout,
            "ok": False,
        }
    except OSError as exc:
        return {
            "status": "launch_error",
            "command": command,
            "exit_code": None,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "error": f"{type(exc).__name__}: {exc}",
            "ok": False,
        }


@tool
def inspect_workflow() -> str:
    """Return the current workflow status, versions, and scientific gates."""
    report = _read_json(REPORT_JSON)
    stages = report["stages"]
    return json.dumps(
        {
            "sample": report["sample"],
            "workflow": report["workflow"],
            "engineering_conclusion": report["engineering_conclusion"],
            "stage_status": {
                "fastp": {"input_pairs": stages["fastp"]["input_pairs"], "clean_pairs": stages["fastp"]["clean_pairs"]},
                "spades": {"contigs": stages["spades"]["contigs"], "reference_coverage": stages["spades"]["reference_coverage_at_95_identity"]},
                "bakta": {"cds": stages["bakta"]["cds"], "database": stages["bakta"]["database"]},
                "interpro": {"panel_sequences": stages["interpro"]["panel_sequences"], "scale": stages["interpro"]["scale_validation"]},
                "openfold3": stages["openfold3"],
            },
            "scientific_limits": report["scientific_limits"],
        },
        ensure_ascii=False,
        indent=2,
    )


@tool
def validate_artifacts() -> str:
    """Check required artifacts, candidate metrics, and semantic failure gates."""
    report = _read_json(REPORT_JSON)
    audit = _read_json(AUDIT_JSON)
    recovery = _read_json(RECOVERY_JSON)
    required = [
        RESULTS / "source/SRR32944547_1.fastq.gz",
        RESULTS / "source/SRR32944547_2.fastq.gz",
        RESULTS / "runs/fastp/SRR32944547.clean.R1.fastq.gz",
        RESULTS / "runs/fastp/SRR32944547.clean.R2.fastq.gz",
        RESULTS / "runs/spades-full2/contigs.fasta",
        RESULTS / "runs/bakta-light-run3/SSi026.faa",
        RESULTS / "runs/blastp-assembly-vs-reference.tsv",
        RESULTS / "runs/interpro-candidate-panel/interpro-panel.json",
    ]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.is_file()]
    semantic_failures = [
        case["name"]
        for case in recovery["cases"]
        if case.get("semantic_failure_detected")
    ]
    unexpected_recovery_cases = [
        case["name"]
        for case in recovery["cases"]
        if case["name"] != "resume_from_completed_fastp_spades"
        and not case.get("passed_expected_failure", False)
    ]
    production_blockers = []
    bakta_db = str(report["stages"]["bakta"].get("database", ""))
    blast_db = str(report["stages"]["blast"].get("database", ""))
    openfold_device = str(report["stages"]["openfold3"].get("device", ""))
    if "light" in bakta_db.lower():
        production_blockers.append("Bakta light database is validation-only")
    if "mps" in openfold_device.lower():
        production_blockers.append("OpenFold3 is running on Apple MPS rather than a Linux CUDA execution tier")
    if "same-strain" in blast_db.lower() or "reference proteins" in blast_db.lower():
        production_blockers.append("BLAST benchmark uses a same-strain reference database rather than an independent external blind set")
    if unexpected_recovery_cases:
        production_blockers.append("recovery checks contain an unexpected failure mode")
    return json.dumps(
        {
            "status": "pass" if not missing and not unexpected_recovery_cases else "blocked",
            "missing_artifacts": missing,
            "semantic_failure_cases_reproduced": semantic_failures,
            "unexpected_recovery_cases": unexpected_recovery_cases,
            "candidate_audit": {
                "predicted_candidates": audit["predicted_candidates"],
                "precision": audit["precision"],
                "recall": audit["recall"],
                "warning": "same-strain reference BLAST is an assembly audit, not an independent discovery benchmark",
            },
            "production_gate": "blocked" if production_blockers else "open",
            "production_ready": not production_blockers,
            "production_blockers": production_blockers,
            "report": str(REPORT_JSON.relative_to(ROOT)),
        },
        ensure_ascii=False,
        indent=2,
    )


@tool
def get_candidate_evidence(protein_id: str) -> str:
    """Return BLAST, Bakta, and InterPro evidence for one assembly protein ID."""
    audit = _read_json(AUDIT_JSON)
    rows = audit.get("rows", audit.get("records", []))
    if not rows:
        audit_tsv = RESULTS / "runs/candidate-audit-v2/candidate-audit.tsv"
        with audit_tsv.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    row = next((item for item in rows if item.get("query") == protein_id), None)
    panel = _read_json(IPR_JSON)
    interpro = next((item for item in panel["run"]["summaries"] if item["id"] == protein_id), None)
    if row is None and interpro is None:
        return json.dumps({"status": "not_found", "protein_id": protein_id})
    return json.dumps(
        {"protein_id": protein_id, "audit": row, "interpro": interpro, "interpretation": "evidence only; not proof of biochemical activity"},
        ensure_ascii=False,
        indent=2,
    )


@tool
def run_approved_validation(stage: Literal["recovery_checks", "build_report", "full_validation"]) -> str:
    """Run approved validation stages with serialized mutations and structured failures."""
    with _stage_lock():
        if stage == "full_validation":
            recovery = _run_approved("recovery_tests.py")
            if not recovery["ok"]:
                return json.dumps(
                    {"stage": stage, "status": "blocked", "ok": False, "steps": [recovery],
                     "reason": "build_report was not started because recovery_checks failed"},
                    ensure_ascii=False,
                    indent=2,
                )
            recovery_report = _read_json(RECOVERY_JSON)
            unexpected = [
                case["name"]
                for case in recovery_report["cases"]
                if case["name"] != "resume_from_completed_fastp_spades"
                and not case.get("passed_expected_failure", False)
            ]
            if unexpected:
                return json.dumps(
                    {"stage": stage, "status": "blocked", "ok": False, "steps": [recovery],
                     "unexpected_recovery_cases": unexpected,
                     "reason": "build_report was not started because a recovery case did not satisfy its expected failure contract"},
                    ensure_ascii=False,
                    indent=2,
                )
            report = _run_approved("build_report.py")
            return json.dumps(
                {"stage": stage, "status": "succeeded" if report["ok"] else "failed",
                 "ok": report["ok"], "steps": [recovery, report]},
                ensure_ascii=False,
                indent=2,
            )

        script = {"recovery_checks": "recovery_tests.py", "build_report": "build_report.py"}[stage]
        if stage == "build_report" and not RECOVERY_JSON.is_file():
            return json.dumps(
                {"stage": stage, "status": "blocked", "ok": False,
                 "reason": "recovery_checks must produce recovery-tests.json before build_report"},
                ensure_ascii=False,
                indent=2,
            )
        result = _run_approved(script)
        return json.dumps({"stage": stage, **result}, ensure_ascii=False, indent=2)


RESEARCH_TOOLS = [*INTENT_TOOLS, *OPENFOLD_TOOLS, *PROTEIN_TOOLS]
TOOLS = [*RESEARCH_TOOLS, inspect_workflow, validate_artifacts, get_candidate_evidence, run_approved_validation]
