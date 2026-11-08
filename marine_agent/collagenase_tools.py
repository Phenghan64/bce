"""Bounded tools for the independent 20-proteome collagenase benchmark.

The tools expose deterministic workflow jobs and evidence inspection to the
Deep Agent. They deliberately do not expose a general shell or filesystem
interface. A returned candidate is a hypothesis supported by sequence/domain
evidence, never a claim of biochemical activity.
"""
from __future__ import annotations

from .execution_control import current_control

import csv
import fcntl
import json
import re
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path

from langchain_core.tools import tool

from .intent_tools import analyze_request_intent
from .protein_tools import PROTEIN_TOOLS
from .openfold_tools import get_uniprot_sequence, resolve_uniprot_sequence, run_openfold_accession, search_uniprot


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/collagenase-agent-20"
MANIFEST = OUT / "manifest.json"
SCREEN_DIR = OUT / "runs/batch-screen"
SCREENING = SCREEN_DIR / "screening.json"
INTERPRO = SCREEN_DIR / "interpro-candidates.json"
RECOVERY = OUT / "runs/recovery.json"
REPORT_JSON = OUT / "REPORT.json"
LOCK = OUT / "runs/.collagenase-agent.lock"


def _read_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"required artifact is missing: {path.relative_to(ROOT)}")
    return json.loads(path.read_text())


def _relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


@contextmanager
def _stage_lock():
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with LOCK.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _run_script(script: str, args: list[str] | None = None, timeout: int = 1800) -> dict:
    command = [sys.executable, str(ROOT / "collagenase_workflow" / script), *(args or [])]
    started = time.perf_counter()
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
            "ok": completed.returncode == 0,
            "exit_code": completed.returncode,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "command": command,
            "output_tail": completed.stdout[-3000:],
        }
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return {
            "status": "timed_out",
            "ok": False,
            "exit_code": None,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "command": command,
            "timeout_seconds": timeout,
            "output_tail": output[-3000:],
        }
    except OSError as exc:
        return {
            "status": "launch_error",
            "ok": False,
            "exit_code": None,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "command": command,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _run_interpro() -> dict:
    command = [
        sys.executable,
        str(ROOT / "collagenase_workflow" / "run_interpro_candidates.py"),
        str(SCREEN_DIR / "candidates.faa"),
        str(INTERPRO),
    ]
    started = time.perf_counter()
    try:
        completed = (current_control().run if current_control() else subprocess.run)(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900)
        return {"status": "succeeded" if completed.returncode == 0 else "failed", "ok": completed.returncode == 0, "exit_code": completed.returncode, "elapsed_seconds": round(time.perf_counter() - started, 3), "command": command, "output_tail": completed.stdout[-3000:]}
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return {"status": "timed_out", "ok": False, "exit_code": None, "elapsed_seconds": round(time.perf_counter() - started, 3), "command": command, "timeout_seconds": 900, "output_tail": output[-3000:]}
    except OSError as exc:
        return {
            "status": "launch_error",
            "ok": False,
            "exit_code": None,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "command": command,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _candidate_rows() -> tuple[list[dict], dict]:
    screening = _read_json(SCREENING)
    with (SCREEN_DIR / "candidate-table.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    interpro = _read_json(INTERPRO)
    by_id = {item["id"]: item for item in interpro["run"]["summaries"]}
    for row in rows:
        evidence = by_id.get(row["qseqid"], {})
        row["interpro_match_count"] = int(evidence.get("match_count", 0))
        row["interpro_pfam"] = evidence.get("pfam", [])
        row["interpro_interpro"] = evidence.get("interpro", [])
        row["domain_gate"] = "PF01752" in row["interpro_pfam"]
    return rows, screening


@tool
def inspect_collagenase_inputs() -> str:
    """Inspect the frozen 20-proteome benchmark and decide whether NGS is needed."""
    manifest = _read_json(MANIFEST)
    fasta_checks = []
    for item in manifest["proteomes"]:
        path = ROOT / item["fasta_path"]
        fasta_checks.append({
            "proteome_id": item["proteome_id"],
            "organism": item["organism_requested"],
            "path": _relative(path),
            "exists": path.is_file(),
            "expected_protein_count": item["protein_count"],
        })
    return json.dumps(
        {
            "status": "ready" if all(x["exists"] for x in fasta_checks) else "blocked",
            "benchmark": manifest["benchmark"],
            "objective": manifest["objective"],
            "proteome_count": len(manifest["proteomes"]),
            "total_expected_proteins": sum(x["protein_count"] for x in manifest["proteomes"]),
            "input_mode": "existing_protein_sequences",
            "ngs_decision": {
                "fastp": "skip",
                "SPAdes": "skip",
                "reason": manifest["ngs_route"]["reason"],
            },
            "traceability": "UniProt public proteomes; see manifest for URLs and SHA-256 hashes",
            "fasta_checks": fasta_checks,
            "scientific_boundary": "Input completeness and provenance do not establish enzyme activity.",
        },
        ensure_ascii=False,
        indent=2,
    )


@tool
def run_collagenase_batch_screen() -> str:
    """Run or resume BLASTP screening and verify the InterPro candidate panel."""
    with _stage_lock():
        result = _run_script("analyze_batch.py", timeout=1800)
        if not result["ok"]:
            return json.dumps({"stage": "batch_screen", **result}, ensure_ascii=False, indent=2)
        interpro_run = _run_interpro()
        if not interpro_run["ok"]:
            return json.dumps({"stage": "batch_screen", "status": "blocked", "ok": False, "blastp": result, "interpro": interpro_run, "reason": "InterPro Matches API verification failed"}, ensure_ascii=False, indent=2)
        screening = _read_json(SCREENING)
        interpro_ready = INTERPRO.is_file()
        return json.dumps(
            {
                "stage": "batch_screen",
                "status": "succeeded" if interpro_ready else "blocked",
                "ok": interpro_ready,
                "blastp": screening["blastp"],
                "interpro": interpro_run,
                "screening": {
                    "proteome_count": screening["proteome_count"],
                    "total_proteins": screening["total_proteins"],
                    "proteins_with_blast_hit": screening["proteins_with_blast_hit"],
                    "candidate_count": screening["candidate_count"],
                    "candidate_count_by_tier": screening["candidate_count_by_tier"],
                    "candidate_table": screening["candidate_table"],
                },
                "interpro_artifact": _relative(INTERPRO) if interpro_ready else None,
                "next_gate": "InterPro/Pfam evidence must be present before structure prediction" if not interpro_ready else "ready",
                "scientific_warning": screening["scientific_warning"],
            },
            ensure_ascii=False,
            indent=2,
        )


@tool
def get_collagenase_candidates(limit: int = 20, require_domain: bool = False) -> str:
    """Return ranked candidates with BLASTP and InterPro evidence."""
    if limit < 1 or limit > 100:
        return json.dumps({"status": "invalid_input", "reason": "limit must be between 1 and 100"}, ensure_ascii=False, indent=2)
    try:
        rows, screening = _candidate_rows()
    except (FileNotFoundError, json.JSONDecodeError, KeyError) as exc:
        return json.dumps({"status": "blocked", "reason": str(exc)}, ensure_ascii=False, indent=2)
    selected = [row for row in rows if not require_domain or row["domain_gate"]]
    selected.sort(key=lambda row: (-float(row["screen_score"]), int(row["rank"])))
    payload = []
    for row in selected[:limit]:
        payload.append({
            "rank": int(row["rank"]),
            "protein_id": row["qseqid"],
            "organism": row["organism"],
            "length": int(row["length"]),
            "blast": {
                "reference_id": row["reference_id"],
                "reference_title": row["reference_title"],
                "identity_percent": float(row["pident"]),
                "query_coverage_percent": float(row["qcov"]),
                "evalue": float(row["evalue"]),
                "tier": row["tier"],
            },
            "interpro": {
                "match_count": row["interpro_match_count"],
                "pfam": row["interpro_pfam"],
                "interpro": row["interpro_interpro"],
                "collagenase_domain_gate": row["domain_gate"],
            },
            "posthoc_label": row["known_collagenase_label_posthoc"].lower() == "true",
        })
    return json.dumps(
        {
            "status": "succeeded",
            "count": len(payload),
            "total_candidates": len(selected),
            "require_domain": require_domain,
            "candidates": payload,
            "screening_rule": screening["candidate_rule"],
            "scientific_warning": "Sequence similarity and domains generate hypotheses; they do not prove collagenase activity, secretion, specificity, or industrial performance.",
        },
        ensure_ascii=False,
        indent=2,
    )


@tool
def assess_collagenase_method(method: str, objective: str, endpoint: str = 'unspecified') -> str:
    """Check whether a method is scientifically suitable for the stated objective."""
    method_norm = method.strip().lower()
    objective_norm = objective.strip().lower()
    # An enzyme as the receptor does not make small-molecule docking invalid.
    # Reject the claimed endpoint (activity/degradation), not its protein name.
    enzyme_activity = any(token in objective_norm for token in (
        "enzyme activity", "catalytic activity", "collagen degradation", "collagenolysis",
        "酶活", "催化活性", "降解活性", "降解能力", "降解率", "活性排名", "耐盐", "耐胆盐"))
    is_vina = "vina" in method_norm or "autodock" in method_norm
    if is_vina and (endpoint in {'enzyme_activity','substrate_degradation','organism_phenotype'} or (endpoint == 'unspecified' and enzyme_activity)):
        return json.dumps(
            {
                "decision": "reject",
                "method": method,
                "objective": objective,
                "reason": "Vina is a small-ligand docking method. A docking score cannot establish catalytic activity, collagen substrate degradation, secretion, or salt/pH stability.",
                "route": ["BLASTP", "InterPro/Pfam", "MSA-backed OpenFold3 for fold/domain review", "expression and biochemical collagenolysis assay"],
                "allowed_vina_scope": "Only for a defined small-molecule ligand question after a validated receptor/ligand setup; never as an enzyme-activity proxy.",
            },
            ensure_ascii=False,
            indent=2,
        )
    accepted = {"blastp", "blastp screening", "interpro", "interproscan", "pfam", "openfold3", "openfold"}
    decision = "accept_with_limits" if method_norm in accepted else "review_required"
    return json.dumps(
        {
            "decision": decision,
            "method": method,
            "objective": objective,
            "reason": "This method can contribute evidence, but its output must be interpreted within the stated scientific limits.",
            "activity_confirmation": "wet-lab assay required",
        },
        ensure_ascii=False,
        indent=2,
    )


@tool
def run_collagenase_structure(candidate_id: str) -> str:
    """Run the gated MSA-backed OpenFold3 prediction for one eligible candidate."""
    if not candidate_id or len(candidate_id) > 240:
        return json.dumps({"status": "invalid_input", "reason": "candidate_id is empty or too long"}, ensure_ascii=False, indent=2)
    # Resume a completed structure job instead of spending another 20+ minutes
    # recomputing the same deterministic query.
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id)
    structure_path = OUT / "runs/openfold" / safe_id / "structure.json"
    if structure_path.is_file():
        try:
            parsed = json.loads(structure_path.read_text())
        except json.JSONDecodeError:
            parsed = {"status": "invalid_structure_artifact"}
        return json.dumps(
            {
                "stage": "openfold3_structure",
                "candidate_id": candidate_id,
                "status": "resumed_existing",
                "execution": {"status": "skipped", "reason": "completed artifact already exists"},
                "artifact": _relative(structure_path),
                "result": parsed,
                "scientific_warning": "OpenFold3 confidence supports fold/domain review; it does not establish collagenase activity or substrate specificity.",
            },
            ensure_ascii=False,
            indent=2,
        )
    with _stage_lock():
        result = _run_script("run_openfold_candidate.py", [candidate_id], timeout=2400)
    parsed = None
    if structure_path.is_file():
        try:
            parsed = json.loads(structure_path.read_text())
        except json.JSONDecodeError:
            parsed = {"status": "invalid_structure_artifact"}
    return json.dumps(
        {
            "stage": "openfold3_structure",
            "candidate_id": candidate_id,
            "execution": result,
            "artifact": _relative(structure_path) if structure_path.is_file() else None,
            "result": parsed,
            "scientific_warning": "OpenFold3 confidence supports fold/domain review; it does not establish collagenase activity or substrate specificity.",
        },
        ensure_ascii=False,
        indent=2,
    )


@tool
def run_collagenase_recovery() -> str:
    """Run the benchmark's resume, invalid-input, missing-input, and method-gate checks."""
    with _stage_lock():
        result = _run_script("recovery.py", timeout=600)
    report = _read_json(RECOVERY) if RECOVERY.is_file() else None
    return json.dumps(
        {
            "stage": "recovery",
            "execution": result,
            "artifact": _relative(RECOVERY) if RECOVERY.is_file() else None,
            "report": report,
        },
        ensure_ascii=False,
        indent=2,
    )


@tool
def build_collagenase_report() -> str:
    """Build the final JSON/Markdown report after recovery checks pass."""
    if not RECOVERY.is_file():
        return json.dumps({"status": "blocked", "reason": "run_collagenase_recovery must complete first"}, ensure_ascii=False, indent=2)
    with _stage_lock():
        result = _run_script("build_report.py", timeout=600)
    return json.dumps(
        {
            "stage": "report",
            "execution": result,
            "report_markdown": _relative(OUT / "REPORT.md") if (OUT / "REPORT.md").is_file() else None,
            "report_json": _relative(REPORT_JSON) if REPORT_JSON.is_file() else None,
        },
        ensure_ascii=False,
        indent=2,
    )


COLLAGENASE_TOOLS = [
    *PROTEIN_TOOLS,
    analyze_request_intent,
    inspect_collagenase_inputs,
    run_collagenase_batch_screen,
    get_collagenase_candidates,
    assess_collagenase_method,
    run_collagenase_structure,
    run_collagenase_recovery,
    build_collagenase_report,
    search_uniprot,
    get_uniprot_sequence,
    resolve_uniprot_sequence,
    run_openfold_accession,
]
