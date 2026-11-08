"""Input-driven protein homology screening with hashed, restartable stages.

This module has no dependency on the frozen marine/collagenase benchmarks.
Sequence similarity produces candidates, not an enzyme activity prediction.
"""
from __future__ import annotations

import csv
import fcntl
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from pydantic import BaseModel, ConfigDict, Field

from .execution_control import RunCancelled, check_cancelled, current_control

ROOT = Path(__file__).resolve().parents[1]
VERSION = "protein-homology-v1"
MAX_FASTA_BYTES = 32 * 1024 * 1024
ANALYSIS_ID = re.compile(r"^[0-9a-f]{64}$")
FIELDS = "qseqid sseqid pident length qlen slen qcovs qstart qend sstart send evalue bitscore"
LIMIT = "Homology candidates only; identity, coverage and score do not establish activity, substrate specificity, secretion or industrial performance. Domain, catalytic-site and experimental verification remain necessary."


class ScreenRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    query_fasta: str = Field(min_length=1)
    reference_fasta: str = Field(min_length=1)
    query_source: str = Field(min_length=1, max_length=4000)
    reference_source: str = Field(min_length=1, max_length=4000)
    objective: str = Field(min_length=1, max_length=2000)
    min_identity: float = Field(default=30, ge=0, le=100)
    min_query_coverage: float = Field(default=50, ge=0, le=100)
    min_reference_coverage: float = Field(default=50, ge=0, le=100)
    max_evalue: float = Field(default=1e-5, gt=0, le=10)
    threads: int = Field(default=6, ge=1, le=256)
    max_target_seqs: int = Field(default=25, ge=1, le=10000)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def parse_fasta(text: str) -> list[dict]:
    if len(text.encode("utf-8")) > MAX_FASTA_BYTES:
        raise ValueError("FASTA exceeds the 32 MiB transport limit; use a registered large-data execution backend")
    records, seen = [], set()
    header, sequence = None, []

    def finish():
        if header is None:
            return
        identifier = header.split()[0]
        if identifier in seen:
            raise ValueError(f"duplicate FASTA identifier: {identifier}")
        seq = "".join(sequence).upper()
        terminal_stop = seq.endswith("*")
        if terminal_stop:
            seq = seq[:-1]
        if not seq or any(char not in "ACDEFGHIKLMNPQRSTVWYBXZJUO" for char in seq):
            raise ValueError(f"invalid or empty protein sequence: {identifier}")
        if len(seq) > 1000000:
            raise ValueError(f"protein sequence exceeds 1,000,000 residues: {identifier}")
        seen.add(identifier)
        records.append({"id": identifier, "header": header, "sequence": seq,
                        "terminal_stop_removed": terminal_stop})

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header, sequence = line[1:].strip(), []
            if not header:
                raise ValueError("FASTA header is empty")
        elif header is None:
            raise ValueError("protein screening requires FASTA records with explicit identifiers")
        else:
            sequence.append("".join(line.split()))
    finish()
    if not records or len(records) > 200000:
        raise ValueError("FASTA must contain 1 to 200,000 proteins")
    return records


def validate_request(payload: dict) -> tuple[ScreenRequest, list[dict], list[dict]]:
    request = ScreenRequest.model_validate(payload)
    if not all(value.strip() for value in (request.query_source, request.reference_source, request.objective)):
        raise ValueError("objective and both source descriptions must be nonempty")
    return request, parse_fasta(request.query_fasta), parse_fasta(request.reference_fasta)


def _executable(name: str) -> str:
    directory = os.getenv("PROTEIN_BLAST_BIN", "")
    candidate = Path(directory) / name if directory else ROOT / ".tools/blast-2.17.0/bin" / name
    selected = str(candidate) if candidate.is_file() else shutil.which(name)
    if not selected:
        raise RuntimeError(f"{name} is not installed on the execution node; configure PROTEIN_BLAST_BIN")
    return selected


def _run(command: list[str], *, timeout: int = 1800):
    check_cancelled()
    return (current_control().run if current_control() else subprocess.run)(
        command, capture_output=True, text=True, timeout=timeout)


def _identity() -> dict:
    result = {}
    for name in ("blastp", "makeblastdb"):
        path = _executable(name)
        completed = _run([path, "-version"], timeout=30)
        if completed.returncode:
            raise RuntimeError(f"{name} version check failed")
        result[name] = {"path": path, "version": completed.stdout.splitlines()[0], "binary_sha256": file_sha(Path(path))}
    return result


def _valid_stage(directory: Path, stage: dict) -> bool:
    artifacts = stage.get("artifacts", {})
    if stage.get("status") != "succeeded" or not artifacts:
        return False
    return all((directory / name).is_file() and file_sha(directory / name) == digest
               for name, digest in artifacts.items())


def _rank(directory: Path, request: ScreenRequest, queries: list[dict], refs: list[dict]) -> dict:
    best, hit_ids = {}, set()
    with (directory / "hits.tsv").open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 13:
                raise ValueError("BLAST output has an invalid column count")
            qid, rid = int(row[0][1:]), int(row[1][1:])
            query, ref = queries[qid], refs[rid]
            identity, qcov, evalue, bitscore = map(float, (row[2], row[6], row[11], row[12]))
            rcov = (abs(int(row[10]) - int(row[9])) + 1) / len(ref["sequence"]) * 100
            if not all(math.isfinite(v) for v in (identity, qcov, evalue, bitscore, rcov)):
                raise ValueError("BLAST output contains a non-finite value")
            if int(row[4]) != len(query["sequence"]) or int(row[5]) != len(ref["sequence"]):
                raise ValueError("BLAST output sequence lengths disagree with submitted input")
            hit_ids.add(query["id"])
            if identity < request.min_identity or qcov < request.min_query_coverage or rcov < request.min_reference_coverage or evalue > request.max_evalue:
                continue
            candidate = {"protein_id": query["id"], "reference_id": ref["id"],
                         "length": len(query["sequence"]), "identity_percent": identity,
                         "query_coverage_percent": qcov, "reference_coverage_percent": round(rcov, 4),
                         "evalue": evalue, "bitscore": bitscore,
                         "exact_reference_sequence": query["sequence"] == ref["sequence"]}
            previous = best.get(query["id"])
            key = lambda value: (value["bitscore"], value["query_coverage_percent"], value["reference_coverage_percent"], value["identity_percent"])
            if previous is None or key(candidate) > key(previous):
                best[query["id"]] = candidate
    ranked = sorted(best.values(), key=lambda x: (-x["bitscore"], x["evalue"], x["protein_id"]))
    for index, item in enumerate(ranked, 1):
        item["rank"] = index
    columns = ["rank", "protein_id", "reference_id", "length", "identity_percent", "query_coverage_percent", "reference_coverage_percent", "evalue", "bitscore", "exact_reference_sequence"]
    with (directory / "candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(ranked)
    selected = set(best)
    (directory / "candidates.faa").write_text("".join(f">{q['id']}\n{q['sequence']}\n" for q in queries if q["id"] in selected))
    report = {"status": "succeeded", "objective": request.objective,
              "input_mode": "existing_protein_sequences", "ngs_stages_skipped": ["fastp", "SPAdes", "Bakta"],
              "query_count": len(queries), "reference_count": len(refs), "proteins_with_hit": len(hit_ids),
              "candidate_count": len(ranked), "candidates": ranked,
              "parameters": request.model_dump(exclude={"query_fasta", "reference_fasta", "objective", "query_source", "reference_source"}),
              "sources": {"query": request.query_source, "reference": request.reference_source, "verification": "user supplied; recorded, not independently certified"},
              "ranking_method": "best passing reference hit per query, ranked by BLAST bitscore; not predicted activity",
              "reference_overlap_count": sum(item["exact_reference_sequence"] for item in ranked),
              "quality_gate": {"decision": "candidate_evidence_only", "activity_confirmed": False, "domain_verification": "not_run"},
              "scientific_limit": LIMIT,
              "search_limit": "One HSP per reference; max_target_seqs limits retained hits. Thresholds are task-specific, not universal enzyme criteria."}
    atomic_json(directory / "report.json", report)
    return report


def run_screen(payload: dict, *, root: Path | None = None, emit: Callable | None = None) -> dict:
    """Run real BLAST, or reuse stages whose entire output set still hashes correctly."""
    started = time.perf_counter()
    request, queries, refs = validate_request(payload)
    versions = _identity()
    inputs = {"query": sha(request.query_fasta.encode()), "reference": sha(request.reference_fasta.encode())}
    spec = {"workflow_version": VERSION, "inputs": inputs, "request": request.model_dump(exclude={"query_fasta", "reference_fasta"}), "tools": versions}
    analysis_id = sha(json.dumps(spec, sort_keys=True, ensure_ascii=False).encode())
    base = (root or ROOT) / "results/protein-screenings"
    base.mkdir(parents=True, exist_ok=True)
    directory = base / analysis_id
    directory.mkdir(exist_ok=True)
    locks = base / ".locks"
    locks.mkdir(exist_ok=True)
    with (locks / analysis_id).open("a+") as lock:
        while True:
            check_cancelled()
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                time.sleep(.1)
        state_path = directory / "manifest.json"
        try:
            state = json.loads(state_path.read_text()) if state_path.is_file() else {}
        except (ValueError, OSError):
            state = {}
        if state.get("spec") != spec:
            state = {"analysis_id": analysis_id, "spec": spec, "stages": {}, "attempts": []}
        attempt = {"started_at": datetime.now(timezone.utc).isoformat(), "status": "running", "reused_stages": []}
        state["attempts"].append(attempt)
        state["status"] = "running"
        atomic_json(state_path, state)
        # Original input and normalized input are both retained. BLAST receives
        # internal IDs so punctuation in user headers cannot corrupt its IDs.
        (directory / "query.original.faa").write_text(request.query_fasta)
        (directory / "reference.original.faa").write_text(request.reference_fasta)
        for name, records, prefix in (("query", queries, "q"), ("reference", refs, "r")):
            (directory / f"{name}.faa").write_text("".join(f">{prefix}{i}\n{r['sequence']}\n" for i, r in enumerate(records)))
        atomic_json(directory / "input-records.json", {"query": queries, "reference": refs})
        dirty = False

        def stage(name, action):
            nonlocal dirty
            check_cancelled()
            previous = state["stages"].get(name, {})
            if not dirty and _valid_stage(directory, previous):
                attempt["reused_stages"].append(name)
                if emit:
                    emit({"name": name, "status": "succeeded", "cache_hit": True})
                return
            dirty = True
            state["current_stage"] = name
            state["stages"][name] = {"status": "running"}
            atomic_json(state_path, state)
            if emit:
                emit({"name": name, "status": "running"})
            tick = time.perf_counter()
            files, metadata = action()
            state["stages"][name] = {"status": "succeeded", "elapsed_seconds": round(time.perf_counter() - tick, 4),
                                      "artifacts": {p.name: file_sha(p) for p in files}, **metadata}
            atomic_json(state_path, state)
            if emit:
                emit({"name": name, "status": "succeeded", "cache_hit": False})

        def command_stage(name, command, outputs):
            completed = _run(command)
            log = directory / f"{name}.log"
            log.write_text(completed.stdout + completed.stderr)
            if completed.returncode:
                raise RuntimeError(f"{name} exited {completed.returncode}; see {log.name}")
            files = outputs()
            if not files or any(not path.is_file() for path in files):
                raise RuntimeError(f"{name} did not create required output files")
            return [log, *files], {"command": command, "exit_code": completed.returncode}

        try:
            def validate():
                qc = {"query_count": len(queries), "reference_count": len(refs),
                      "inputs": inputs, "terminal_stops_removed": sum(r["terminal_stop_removed"] for r in queries + refs),
                      "ambiguous_residue_records": sum(any(c in r["sequence"] for c in "BXZJUO") for r in queries + refs),
                      "nucleotide_like_record_count": sum(set(r["sequence"]) <= set("ACGTUN") for r in queries + refs)}
                atomic_json(directory / "input-qc.json", qc)
                return [directory / name for name in ("input-qc.json", "query.original.faa", "reference.original.faa", "query.faa", "reference.faa", "input-records.json")], {}
            stage("validate_inputs", validate)
            stage("makeblastdb", lambda: command_stage("makeblastdb", [versions["makeblastdb"]["path"], "-in", str(directory / "reference.faa"), "-dbtype", "prot", "-out", str(directory / "reference_db")], lambda: list(directory.glob("reference_db.*"))))
            stage("blastp", lambda: command_stage("blastp", [versions["blastp"]["path"], "-query", str(directory / "query.faa"), "-db", str(directory / "reference_db"), "-out", str(directory / "hits.tsv"), "-outfmt", "6 " + FIELDS, "-evalue", str(request.max_evalue), "-max_target_seqs", str(request.max_target_seqs), "-max_hsps", "1", "-seg", "yes", "-num_threads", str(request.threads)], lambda: [directory / "hits.tsv"]))

            def report():
                _rank(directory, request, queries, refs)
                return [directory / name for name in ("report.json", "candidates.tsv", "candidates.faa")], {}
            stage("rank_and_report", report)
            state["status"] = attempt["status"] = "succeeded"
            state["current_stage"] = "completed"
        except BaseException as exc:
            state["status"] = attempt["status"] = "stopped" if isinstance(exc, RunCancelled) else "failed"
            attempt["error"] = f"{type(exc).__name__}: {exc}"
            if state.get("current_stage") in state["stages"]:
                state["stages"][state["current_stage"]]["status"] = state["status"]
            raise
        finally:
            attempt["elapsed_seconds"] = round(time.perf_counter() - started, 4)
            attempt["finished_at"] = datetime.now(timezone.utc).isoformat()
            atomic_json(state_path, state)
    return read_analysis(analysis_id, root=root)


def read_analysis(analysis_id: str, *, root: Path | None = None, protein_id: str = "") -> dict:
    if not ANALYSIS_ID.fullmatch(analysis_id):
        raise ValueError("invalid protein analysis ID")
    directory = (root or ROOT) / "results/protein-screenings" / analysis_id
    state = json.loads((directory / "manifest.json").read_text())
    if state.get("analysis_id") != analysis_id or sha(json.dumps(state.get("spec", {}), sort_keys=True, ensure_ascii=False).encode()) != analysis_id:
        return {"status": "invalid_artifacts", "analysis_id": analysis_id, "reason": "Analysis manifest identity does not match its recorded specification"}
    if state.get("status") != "succeeded":
        return {"status": state.get("status"), "analysis_id": analysis_id, "current_stage": state.get("current_stage"), "reason": "Analysis is incomplete; resubmit the exact original inputs to resume validated stages"}
    if not all(_valid_stage(directory, state["stages"].get(name, {})) for name in ("validate_inputs", "makeblastdb", "blastp", "rank_and_report")):
        return {"status": "invalid_artifacts", "analysis_id": analysis_id, "reason": "An artifact is missing or modified; resubmit the exact original inputs to rebuild affected stages"}
    report = json.loads((directory / "report.json").read_text())
    candidates = report.pop("candidates")
    report.update({"analysis_id": analysis_id, "candidates": candidates[:50], "returned_candidates": min(50, len(candidates)),
                   "artifacts": {name: f"results/protein-screenings/{analysis_id}/{name}" for name in ("report.json", "manifest.json", "candidates.tsv", "candidates.faa", "input-qc.json")},
                   "execution": state["attempts"][-1], "workflow_version": VERSION})
    if protein_id:
        candidate = next((c for c in candidates if c["protein_id"] == protein_id), None)
        if candidate is None:
            raise ValueError("protein_id is not a passing candidate in this analysis")
        records = json.loads((directory / "input-records.json").read_text())
        record = next(r for r in records["query"] if r["id"] == protein_id)
        report["selected_candidate"] = {**candidate, "sequence": record["sequence"], "sequence_sha256": sha(record["sequence"].encode())}
    return report
