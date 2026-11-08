"""Versioned domain evidence via exact lookup or an explicit full InterProScan.

Lookup remains the fast default. Full scans use the authenticated Linux worker.
Missing records and unobserved signatures are not biochemical negatives.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import time
from http.client import HTTPException
from datetime import datetime, timezone
from urllib.request import Request, build_opener, ProxyHandler

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from .access_control import grant_analysis, require_analysis
from .execution_control import RunCancelled, check_cancelled, current_control
from .protein_screening import ROOT, ANALYSIS_ID, atomic_json, file_sha, parse_fasta, read_analysis, sha

API = os.getenv("INTERPRO_MATCHES_URL", "https://www.ebi.ac.uk/interpro/matches/api").rstrip("/")
VERSION = "exact-sequence-domain-evidence-v2"
PFAM = re.compile(r"^PF[0-9]{5}$")
MD5 = re.compile(r"^[0-9A-F]{32}$")


class DomainRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    query_fasta: str = ""
    source: str = Field(default="", max_length=4000)
    analysis_id: str = ""
    protein_id: str = Field(default="", max_length=1000)
    expected_pfam: list[str] = Field(default_factory=list, max_length=100)
    refresh: bool = False
    method: Literal["lookup", "full_scan"] = "lookup"


def now():
    return datetime.now(timezone.utc).isoformat()


def resolve_inputs(payload: dict, *, root: Path | None = None):
    request = DomainRequest.model_validate(payload)
    if bool(request.query_fasta) == bool(request.analysis_id):
        raise ValueError("provide either query_fasta/source or a protein screening analysis_id")
    if any(not PFAM.fullmatch(value) for value in request.expected_pfam):
        raise ValueError("expected_pfam must contain explicit PF00000-format accessions")
    if request.analysis_id:
        require_analysis('protein-screenings', request.analysis_id)
        report = read_analysis(request.analysis_id, root=root)
        if report.get("status") != "succeeded":
            raise ValueError("source protein analysis is incomplete or its artifacts failed integrity checks")
        directory = (root or ROOT) / "results/protein-screenings" / request.analysis_id
        text = (directory / "candidates.faa").read_text()
        records = parse_fasta(text) if text.strip() else []
        if request.protein_id:
            records = [record for record in records if record["id"] == request.protein_id]
        source = {"kind": "protein_screening", "analysis_id": request.analysis_id,
                  "protein_id": request.protein_id or None, "sources": report["sources"]}
    else:
        if request.protein_id or not request.source.strip():
            raise ValueError("standalone FASTA requires source; protein_id selects only a saved screening candidate")
        records = parse_fasta(request.query_fasta)
        source = {"kind": "provided_fasta", "description": request.source}
    if not records:
        raise ValueError("no proteins selected for domain lookup")
    return request, records, source


def _fetch(url: str, payload: dict | None = None) -> tuple[dict, dict]:
    client = build_opener(ProxyHandler({}))
    started = time.perf_counter()
    for attempt in range(1, 4):
        check_cancelled()
        try:
            request = Request(url, data=None if payload is None else json.dumps(payload).encode(),
                              headers={"Content-Type": "application/json", "Accept": "application/json"})
            with client.open(request, timeout=15) as response:
                raw = response.read(32 * 1024 * 1024 + 1)
                if len(raw) > 32 * 1024 * 1024:
                    raise ValueError("InterPro response exceeded the 32 MiB chunk limit")
                result = json.loads(raw)
            check_cancelled()
            if not isinstance(result, dict):
                raise ValueError("InterPro response must be a JSON object")
            return result, {"url": url, "retrieved_at": now(), "attempt": attempt,
                            "elapsed_seconds": round(time.perf_counter() - started, 4), "response_sha256": sha(raw)}
        except (OSError, ValueError, HTTPException) as exc:
            if attempt == 3:
                raise RuntimeError(f"InterPro lookup failed after {attempt} attempts: {type(exc).__name__}: {exc}") from exc
            deadline = time.monotonic() + attempt
            while time.monotonic() < deadline:
                check_cancelled()
                time.sleep(.1)
    raise AssertionError("unreachable")


def sequence_md5(record):
    return hashlib.md5(record["sequence"].encode()).hexdigest().upper()


def _check_chunk(value: dict, requested: list[str]) -> dict:
    results = value.get("results")
    if not isinstance(results, list):
        raise ValueError("InterPro response has no results list")
    found = {}
    for item in results:
        if not isinstance(item, dict):
            raise ValueError("invalid InterPro result record")
        digest = str(item.get("md5", "")).upper()
        if not MD5.fullmatch(digest) or digest not in requested or digest in found:
            raise ValueError("InterPro returned an unexpected or duplicate sequence MD5")
        if not isinstance(item.get("found"), bool) or not isinstance(item.get("matches", []), list):
            raise ValueError("InterPro result has invalid found/matches fields")
        if not item["found"] and item.get("matches"):
            raise ValueError("InterPro returned contradictory found/matches fields")
        found[digest] = item
    return found


def _record_evidence(record: dict, item: dict | None, expected: list[str]) -> dict:
    if item is None:
        status = "lookup_not_returned"
    elif not item["found"]:
        status = "not_in_precomputed_database"
    else:
        status = "matches_available" if item.get("matches") else "precomputed_record_without_matches"
    signatures, entries = [], {}
    for match in (item or {}).get("matches", []):
        signature = match.get("signature") or {}
        accession = str(signature.get("accession") or "")
        release = signature.get("signatureLibraryRelease") or {}
        entry = signature.get("entry") or {}
        if not accession or not release.get("library") or not release.get("version"):
            raise ValueError("InterPro signature lacks an accession or library version")
        locations = []
        for location in match.get("locations", []):
            start, end = location.get("start"), location.get("end")
            if type(start) is not int or type(end) is not int or not 1 <= start <= end <= len(record["sequence"]):
                raise ValueError("InterPro coordinates do not match the submitted sequence length")
            fragments = location.get("location-fragments", [])
            for fragment in fragments:
                left, right = fragment.get("start"), fragment.get("end")
                if type(left) is not int or type(right) is not int or not start <= left <= right <= end:
                    raise ValueError("InterPro fragment coordinates are invalid")
            locations.append({"start": start, "end": end, "fragments": fragments,
                              "score": location.get("score"), "evalue": location.get("evalue")})
        signatures.append({"accession": accession, "name": signature.get("name"),
                           "description": signature.get("description"), "type": signature.get("type"),
                           "library": release["library"], "library_version": release["version"],
                           "analysis_source": match.get("source", release["library"]),
                           "locations": locations, "entry": entry or None})
        if entry.get("accession"):
            entries[entry["accession"]] = {"accession": entry["accession"], "name": entry.get("name"),
                                           "description": entry.get("description"), "type": entry.get("type")}
    pfam = [sig for sig in signatures if sig["library"].lower() == "pfam" and sig["analysis_source"].lower() == "pfam"]
    predicted = [sig for sig in signatures if sig["library"].lower() == "pfam" and sig["analysis_source"].lower() != "pfam"]
    observed = {sig["accession"].split(".")[0] for sig in pfam}
    predicted_ids = {sig["accession"].split(".")[0] for sig in predicted}
    known = item is not None and item.get("found") is True
    return {"protein_id": record["id"], "sequence_length": len(record["sequence"]),
            "sequence_sha256": sha(record["sequence"].encode()), "md5": sequence_md5(record),
            "lookup_status": status, "match_count": len(signatures), "signatures": signatures,
            "pfam": pfam, "predicted_pfam": predicted, "interpro": sorted(entries.values(), key=lambda x: x["accession"]),
            "expected_pfam": {value: "observed" if value in observed else "predicted_only" if value in predicted_ids else "not_observed_in_precomputed_record" if known else "unknown" for value in expected},
            "full_scan_required_for_new_sequence": not known}


def _cached(directory: Path, entry: dict) -> bool:
    path = directory / str(entry.get("file", "missing"))
    return bool(entry.get("sha256") and path.is_file() and file_sha(path) == entry["sha256"])


def annotate(payload: dict, *, root: Path | None = None, emit=None) -> dict:
    started = time.perf_counter()
    request, records, source = resolve_inputs(payload, root=root)
    if request.method == "full_scan":
        from .interpro_remote import annotate_full_scan
        return annotate_full_scan(request, records, source, root=root or ROOT, emit=emit)
    info, info_meta = _fetch(API + "/info")
    if not info.get("release") or not isinstance(info.get("analyses"), list):
        raise ValueError("InterPro did not identify its database release and analyses")
    spec = {"workflow_version": VERSION, "api": API, "release_info": info,
            "records": [{"id": record["id"], "sequence_sha256": sha(record["sequence"].encode())} for record in records],
            "source": source, "expected_pfam": sorted(set(request.expected_pfam))}
    analysis_id = sha(json.dumps(spec, sort_keys=True, ensure_ascii=False).encode())
    directory = (root or ROOT) / "results/protein-domains" / analysis_id
    directory.mkdir(parents=True, exist_ok=True)
    lock_root = directory.parent / ".locks"
    lock_root.mkdir(exist_ok=True)
    with (lock_root / analysis_id).open("a+") as lock:
        while True:
            check_cancelled()
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                time.sleep(.1)
        manifest_path = directory / "manifest.json"
        try:
            state = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        except (ValueError, OSError):
            state = {}
        if state.get("spec") != spec or not isinstance(state.get("chunks"), dict) or not isinstance(state.get("attempts"), list):
            state = {"analysis_id": analysis_id, "spec": spec, "chunks": {}, "attempts": []}
        attempt = {"started_at": now(), "status": "running", "reused_chunks": [], "refresh": request.refresh}
        state["attempts"].append(attempt)
        state["status"] = "running"
        atomic_json(manifest_path, state)
        try:
            atomic_json(directory / "input-records.json", {"records": records, "source": source})
            atomic_json(directory / "api-info.json", {"info": info, "request": info_meta})
            md5s = sorted({sequence_md5(record) for record in records})
            representative = {sequence_md5(record): record for record in records}
            by_md5 = {}
            for offset in range(0, len(md5s), 100):
                check_cancelled()
                chunk = md5s[offset:offset + 100]
                key = str(offset // 100)
                name = f"chunk-{key}.json"
                state["current_chunk"] = key
                previous = state["chunks"].get(key, {})
                if not request.refresh and previous.get("requested_md5") == chunk and _cached(directory, previous):
                    artifact = json.loads((directory / name).read_text())
                    attempt["reused_chunks"].append(key)
                else:
                    if emit:
                        emit({"name": "domain_lookup", "status": "running", "chunk": int(key) + 1, "total": (len(md5s) + 99) // 100})
                    check_cancelled()
                    raw, meta = _fetch(API + "/matches", {"md5": chunk})
                    for digest, item in _check_chunk(raw, chunk).items():
                        _record_evidence(representative[digest], item, spec["expected_pfam"])
                    artifact = {"requested_md5": chunk, "response": raw, "request": meta}
                    atomic_json(directory / name, artifact)
                    state["chunks"][key] = {"file": name, "sha256": file_sha(directory / name), "requested_md5": chunk}
                    atomic_json(manifest_path, state)
                by_md5.update(_check_chunk(artifact["response"], chunk))
                if emit:
                    emit({"name": "domain_lookup", "status": "succeeded", "chunk": int(key) + 1,
                          "total": (len(md5s) + 99) // 100, "cache_hit": key in attempt["reused_chunks"]})
            evidence = [_record_evidence(record, by_md5.get(sequence_md5(record)), spec["expected_pfam"]) for record in records]
            unknown = sum(record["full_scan_required_for_new_sequence"] for record in evidence)
            status = "partial" if unknown else "succeeded"
            report = {"status": status, "kind": "protein_domain_evidence", "analysis_id": analysis_id,
                      "workflow_version": VERSION, "method": "InterPro Matches API exact-sequence precomputed lookup",
                      "de_novo_scan_executed": False, "source": source, "api": API,
                      "release": info["release"], "release_date": info.get("release_date"),
                      "library_versions": info["analyses"], "record_count": len(records), "unique_sequence_count": len(md5s),
                      "unknown_record_count": unknown, "match_count": sum(record["match_count"] for record in evidence),
                      "records": evidence, "quality_gate": {"decision": "incomplete_evidence" if unknown else "precomputed_sequence_evidence_available",
                      "activity_confirmed": False, "catalytic_sites_verified": False},
                      "limits": ["Precomputed records only; unknown MD5s require a full InterProScan execution backend.",
                                 "A domain match supports family evidence, not activity or substrate specificity.",
                                 "InterPro-N predictions are separate evidence from member-database analyses, with their own versions and scores.",
                                 "An unobserved domain is not an experimental negative or proof of absence."]}
            atomic_json(directory / "report.json", report)
            state["artifacts"] = {name: file_sha(directory / name) for name in ("input-records.json", "api-info.json", "report.json")}
            state["status"] = attempt["status"] = status
        except BaseException as exc:
            state["status"] = attempt["status"] = "stopped" if isinstance(exc, RunCancelled) else "failed"
            attempt["error"] = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            attempt["finished_at"] = now()
            attempt["elapsed_seconds"] = round(time.perf_counter() - started, 4)
            atomic_json(manifest_path, state)
    grant_analysis('protein-domains', analysis_id)
    return read_domain_analysis(analysis_id, root=root)


def read_domain_analysis(analysis_id: str, *, root: Path | None = None, protein_id: str = "") -> dict:
    if not ANALYSIS_ID.fullmatch(analysis_id):
        raise ValueError("invalid domain analysis ID")
    require_analysis('protein-domains', analysis_id)
    directory = (root or ROOT) / "results/protein-domains" / analysis_id
    state = json.loads((directory / "manifest.json").read_text())
    if state.get("analysis_id") != analysis_id or sha(json.dumps(state.get("spec", {}), sort_keys=True, ensure_ascii=False).encode()) != analysis_id:
        raise ValueError("domain analysis identity is invalid")
    if state.get("status") not in {"succeeded", "partial"}:
        return {"status": state.get("status"), "analysis_id": analysis_id, "reason": "lookup incomplete; resubmit exact inputs to continue"}
    expected_artifacts = {'input-records.json', 'api-info.json', 'report.json'}
    if state['spec'].get('execution_method') == 'full_scan':
        from .interpro_remote import FILES
        expected_artifacts = FILES | {'input-records.json', 'worker-info.json', 'scan-receipt.json', 'report.json'}
    if set(state.get('artifacts', {})) != expected_artifacts:
        raise ValueError('domain manifest lacks required artifact hashes; repeat original lookup to repair')
    checks = [{"file": name, "sha256": digest} for name, digest in state.get("artifacts", {}).items()] + list(state["chunks"].values())
    if not checks or not all(_cached(directory, item) for item in checks):
        raise ValueError("domain artifact is missing or modified; repeat original lookup to repair")
    report = json.loads((directory / "report.json").read_text())
    if state['spec'].get('execution_method') == 'full_scan':
        from .interpro_protocol import active_library_versions
        # The release metadata also lists optional libraries that were not run.
        # Preserve that raw artifact, but expose only the selected scan scope.
        report['library_versions'] = active_library_versions(state['spec']['release_info'])
    records = report["records"]
    if protein_id:
        records = [record for record in records if record["protein_id"] == protein_id]
        if not records:
            raise ValueError("protein_id does not belong to this domain analysis")
    report["records"] = [{k: v for k, v in record.items() if k != "signatures"} for record in records[:10]]
    report["returned_records"] = len(report["records"])
    report["execution"] = state["attempts"][-1]
    report["artifacts"] = {name: f"results/protein-domains/{analysis_id}/{name}" for name in ["manifest.json", *sorted(expected_artifacts)]}
    return report
