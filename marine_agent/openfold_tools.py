"""Bounded UniProt and OpenFold3 tools for sequence and structure jobs."""
from __future__ import annotations

import hashlib
import json
import base64
import os
import re
import subprocess
import sys
import tarfile
import tempfile
import time
import uuid
from io import BytesIO
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from langchain_core.tools import tool

from openfold3.predict_accession import fetch_uniprot, normalize_ligands

from .search_assist import assist_query, translate_candidate_names
from .execution_control import current_control, check_cancelled

ROOT = Path(__file__).resolve().parents[1]
# UniProtKB accessions are six or ten alphanumeric characters. The ten-character
# form must contain a digit so ordinary protein names cannot bypass search.
ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|(?=[A-Z0-9]{10}$)(?=.*[0-9])[A-Z0-9]{10})$",
    re.IGNORECASE,
)
PYTHON = ROOT / ".env/bin/python" if (ROOT / ".env/bin/python").is_file() else Path(sys.executable)
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_SEARCH_CACHE = ROOT / "results/openfold3-uniprot-search-cache"
OPENFOLD_WORKER_URL = re.sub(r"/+$", "", os.environ.get("OPENFOLD_WORKER_URL", ""))
OPENFOLD_WORKER_TOKEN = os.environ.get("OPENFOLD_WORKER_TOKEN", "")
OPENFOLD_WORKER_TIMEOUT = max(60, int(os.environ.get("OPENFOLD_WORKER_TIMEOUT_S", "3600")))


def _json_result(payload: dict) -> str:
    return json.dumps(payload, ensure_ascii=False, indent=2)


def _atomic_json_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, ensure_ascii=False)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def _safe_extract_bundle(encoded: str, expected_root: str) -> None:
    """Extract a worker artifact bundle below ROOT with path and link checks."""
    try:
        archive_bytes = base64.b64decode(encoded, validate=True)
    except (ValueError, TypeError) as exc:
        raise RuntimeError("OpenFold worker returned an invalid artifact bundle") from exc
    if len(archive_bytes) > 512 * 1024 * 1024:
        raise RuntimeError("OpenFold worker artifact bundle exceeds the 512 MiB limit")
    destination = (ROOT / expected_root).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=BytesIO(archive_bytes), mode="r:gz") as archive:
        members = archive.getmembers()
        if not members:
            raise RuntimeError("OpenFold worker returned an empty artifact bundle")
        for member in members:
            member_path = Path(member.name)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise RuntimeError("OpenFold worker returned an unsafe artifact path")
            if member.issym() or member.islnk() or not (member.isfile() or member.isdir()):
                raise RuntimeError("OpenFold worker returned an unsupported artifact entry")
            resolved = (ROOT / member.name).resolve()
            try:
                resolved.relative_to(ROOT.resolve())
            except ValueError as exc:
                raise RuntimeError("OpenFold worker artifact escaped the result directory") from exc
        for member in members:
            archive.extract(member, ROOT)
    if not destination.exists():
        raise RuntimeError("OpenFold worker bundle did not contain the declared output directory")


def _run_openfold_worker(accession: str | None, ligands: list[dict], *, sequence: str | None = None,
                         sequence_id: str = "sequence") -> dict:
    """Submit one deterministic accession or inline-sequence job to the Mac worker."""
    if not OPENFOLD_WORKER_URL:
        return {
            "status": "unavailable",
            "stage": "execution",
            "accession": accession,
            "reason": "OpenFold3 execution node is not configured; set OPENFOLD_WORKER_URL or run the API on a host with .tools/openfold3-env.",
        }
    job_id = uuid.uuid4().hex
    output_root = f"results/openfold3-remote/{job_id}"
    payload = {"job_id": job_id, "ligands": ligands, "output_root": output_root}
    check_cancelled()
    if control := current_control():
        payload['cancellation_id'] = control.job_id
    if sequence is not None:
        payload.update({"sequence": sequence, "sequence_id": sequence_id})
    else:
        payload["accession"] = accession
    body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if OPENFOLD_WORKER_TOKEN:
        headers["X-OpenFold-Worker-Token"] = OPENFOLD_WORKER_TOKEN
    request = Request(f"{OPENFOLD_WORKER_URL}/v1/openfold", data=body, headers=headers, method="POST")
    try:
        from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeout
        from urllib.request import build_opener, ProxyHandler
        from .progress import emit
        def submit():
            with build_opener(ProxyHandler({})).open(request, timeout=OPENFOLD_WORKER_TIMEOUT) as response:
                return json.load(response)
        previous_stage = None
        with ThreadPoolExecutor(max_workers=1, thread_name_prefix='openfold-transport') as transport:
            future = transport.submit(submit)
            while True:
                check_cancelled()
                try:
                    response_payload = future.result(timeout=1)
                    break
                except FutureTimeout:
                    if future.done():
                        raise
                    try:
                        progress_request = Request(f'{OPENFOLD_WORKER_URL}/v1/jobs/{job_id}/progress', headers=headers)
                        with build_opener(ProxyHandler({})).open(progress_request, timeout=2) as response:
                            progress = json.load(response)
                        stage = progress.get('stage')
                        if progress.get('job_id') == job_id and stage in {'sequence','msa','prediction_queue','prediction','viewer_export'} and stage != previous_stage:
                            emit(stage)
                            previous_stage = stage
                    except (OSError, ValueError):
                        pass
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")[-2000:]
        return {"status": "failed", "stage": "worker", "accession": accession, "reason": f"worker HTTP {exc.code}: {detail}"}
    except (URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
        return {"status": "unavailable", "stage": "worker", "accession": accession, "reason": str(exc)}
    if not isinstance(response_payload, dict):
        return {"status": "failed", "stage": "worker", "accession": accession, "reason": "worker returned a non-object response"}
    result = response_payload.get("result")
    if not isinstance(result, dict):
        result = {"status": response_payload.get("status", "failed"), "reason": response_payload.get("reason", "worker result is missing")}
    bundle = response_payload.get("bundle")
    if isinstance(bundle, dict) and isinstance(bundle.get("data"), str):
        declared_root = str(bundle.get("root") or output_root)
        if declared_root != output_root:
            return {"status": "failed", "stage": "worker", "accession": accession, "reason": "worker output root does not match the submitted job"}
        try:
            _safe_extract_bundle(bundle["data"], declared_root)
        except RuntimeError as exc:
            return {"status": "failed", "stage": "worker_artifacts", "accession": accession, "reason": str(exc)}
    result.setdefault("worker", {"url": OPENFOLD_WORKER_URL, "job_id": job_id})
    return result


def _parse_ligands_json(ligands_json: str) -> list[dict]:
    if not ligands_json.strip():
        return []
    parsed = json.loads(ligands_json)
    values = parsed.get("ligands", []) if isinstance(parsed, dict) else parsed
    if any(isinstance(item, dict) and (item.get("sdf_path") or item.get("sdf_file_path")) for item in values):
        raise ValueError("agent ligand input cannot reference SDF files; use inline SMILES or CCD codes")
    return normalize_ligands(parsed)


def _normalize_inline_sequence(value: str) -> str:
    """Normalize a pasted raw sequence or one-record FASTA without accepting prose."""
    text = str(value or "").strip()
    if text.startswith(">"):
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if len(lines) < 2:
            raise ValueError("FASTA sequence is missing residues")
        text = "".join(lines[1:])
    sequence = re.sub(r"\s+", "", text).upper()
    if len(sequence) < 20 or len(sequence) > 4000:
        raise ValueError("sequence must contain 20-4000 amino acids")
    invalid = sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"))
    if invalid:
        raise ValueError(f"sequence contains unsupported residues: {','.join(invalid)}")
    return sequence


def _search_cache_key(query: str, organism: str, reviewed: bool | None, limit: int) -> str:
    canonical = json.dumps(
        {"query": query, "organism": organism, "reviewed": reviewed, "limit": limit},
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _fetch_search_json(url: str) -> dict:
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "marine-scientific-agent/1.0 (UniProt search)",
        },
    )
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            with urlopen(request, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            last_error = exc
            transient = exc.code == 429 or 500 <= exc.code <= 599
            if not transient or attempt == 3:
                raise RuntimeError(f"UniProt search failed: HTTP {exc.code}") from exc
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            try:
                delay = min(8.0, max(0.5, float(retry_after))) if retry_after else min(8.0, 2**attempt)
            except ValueError:
                delay = min(8.0, 2**attempt)
            time.sleep(delay)
        except (URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt == 3:
                raise RuntimeError(f"UniProt search failed: {exc}") from exc
            time.sleep(min(8.0, 2**attempt))
    raise RuntimeError(f"UniProt search failed: {last_error}") from last_error


def _protein_name(entry: dict) -> str | None:
    description = entry.get("proteinDescription") or {}
    recommended = description.get("recommendedName") or {}
    full_name = recommended.get("fullName") or {}
    if full_name.get("value"):
        return full_name["value"]
    submission = description.get("submissionNames") or []
    if submission:
        return ((submission[0].get("fullName") or {}).get("value"))
    return None


def _search_rows(payload: dict) -> list[dict]:
    rows = []
    for entry in payload.get("results", []):
        accession = entry.get("primaryAccession")
        sequence = entry.get("sequence") or {}
        organism = entry.get("organism") or {}
        if not accession:
            continue
        rows.append(
            {
                "accession": accession,
                "entry_type": entry.get("entryType"),
                "reviewed": entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
                "entry_name": entry.get("uniProtkbId"),
                "protein_name": _protein_name(entry),
                "organism": organism.get("scientificName"),
                "organism_taxon_id": organism.get("taxonId"),
                "length": sequence.get("length"),
                "fragment": str((entry.get("proteinDescription") or {}).get("flag", "")).lower() in {"fragment", "fragments"} or bool(sequence.get("fragment")),
                "sequence_url": f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
            }
        )
    return rows


@tool
def search_uniprot(query: str, organism: str = "", reviewed: bool | None = None, limit: int = 10) -> str:
    """Search UniProt by protein/function/name and return traceable candidates.

    This is a resolver step, not an automatic accession choice. The caller
    should confirm organism, protein name, length, and reviewed status before
    requesting a sequence or starting OpenFold3.
    """
    query = query.strip()
    organism = organism.strip()
    if len(query) < 2 or len(query) > 160:
        return _json_result({"status": "invalid_input", "reason": "query must contain 2-160 characters"})
    if len(organism) > 120:
        return _json_result({"status": "invalid_input", "reason": "organism must contain at most 120 characters"})
    if limit < 1 or limit > 20:
        return _json_result({"status": "invalid_input", "reason": "limit must be between 1 and 20"})

    assistance = assist_query(query, organism)
    terms = assistance.get("search_terms") or [query]
    # UniProt's query language treats spaces as an implicit AND. For a
    # model-generated synonym list, use the first canonical phrase here; the
    # frontend can issue the remaining phrase as a separate auditable query.
    search_query = str(terms[0]).strip() or query
    terms = [search_query]
    if organism:
        terms.append(f'organism_name:"{organism.replace(chr(34), "")}"')
    if reviewed is True:
        terms.append("reviewed:true")
    if reviewed is False:
        terms.append("reviewed:false")
    uniprot_query = " AND ".join(terms)
    params = {
        "query": uniprot_query,
        "format": "json",
        "size": str(limit),
        "fields": "accession,id,protein_name,organism_name,length,reviewed,sequence",
    }
    url = f"{UNIPROT_SEARCH_URL}?{urlencode(params)}"
    cache_key = _search_cache_key(query, organism, reviewed, limit)
    cache_path = UNIPROT_SEARCH_CACHE / f"{cache_key}.json"
    started = time.perf_counter()
    cache_hit = False
    try:
        if cache_path.is_file():
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            cache_hit = True
        else:
            payload = _fetch_search_json(url)
            _atomic_json_write(cache_path, payload)
        candidates = _search_rows(payload)
        translated = translate_candidate_names(candidates)
        for candidate in candidates:
            candidate["protein_name_zh"] = translated["names"].get(candidate["accession"])
    except (OSError, json.JSONDecodeError, RuntimeError) as exc:
        return _json_result(
            {
                "status": "unavailable",
                "query": query,
                "search_query": search_query,
                "organism": organism or None,
                "reason": str(exc),
                "source": url,
            }
        )
    return _json_result(
        {
            "status": "succeeded",
            "query": query,
            "search_query": search_query,
            "normalized_queries": assistance.get("search_terms") or [query],
            "query_assistance": {
                "method": assistance.get("method"),
                "ai_used": bool(assistance.get("ai_used")),
                "model": assistance.get("model"),
            },
            "candidate_translation": {
                "ai_used": bool(translated.get("ai_used")),
                "model": translated.get("model"),
                "field": "protein_name_zh",
            },
            "organism": organism or None,
            "reviewed_filter": reviewed,
            "count": len(candidates),
            "candidates": candidates,
            "cache_hit": cache_hit,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "source": url,
            "cache_path": str(cache_path.relative_to(ROOT)),
            "next_step": "Confirm one accession before calling get_uniprot_sequence or run_openfold_accession.",
            "scientific_limit": "UniProt search results identify candidate records; they do not confirm function, activity, substrate specificity, or expression context.",
        }
    )


@tool
def get_uniprot_sequence(accession: str) -> str:
    """Fetch one validated UniProt sequence and provenance by explicit accession."""
    accession = accession.strip()
    if not ACCESSION_RE.fullmatch(accession):
        return _json_result({"status": "invalid_input", "reason": "an explicit UniProt accession is required"})
    started = time.perf_counter()
    try:
        metadata, sequence = fetch_uniprot(accession)
    except (OSError, RuntimeError, ValueError) as exc:
        return _json_result({"status": "unavailable", "accession": accession, "reason": str(exc)})
    return _json_result(
        {
            "status": "succeeded",
            "accession": metadata.get("primary_accession", accession),
            "metadata": metadata,
            "sequence": sequence,
            "fasta": f">{metadata.get('primary_accession', accession)}\n{sequence}\n",
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "scientific_limit": "The sequence and UniProt annotation are provenance evidence; they do not establish biochemical activity or suitability for a downstream method.",
        }
    )


@tool
def resolve_uniprot_sequence(query: str, organism: str = "", reviewed: bool | None = None, limit: int = 5) -> str:
    """Resolve a natural-language protein request to a sequence when safe.

    An explicit accession is fetched directly. A search with exactly one
    returned candidate fetches that sequence; multiple candidates are returned
    without choosing one so the user can confirm the intended record.
    """
    query = query.strip()
    if ACCESSION_RE.fullmatch(query):
        return _json_result(
            {
                "status": "resolved",
                "resolution": "explicit_accession",
                "result": json.loads(get_uniprot_sequence.invoke({"accession": query})),
            }
        )
    search_result = json.loads(
        search_uniprot.invoke(
            {"query": query, "organism": organism, "reviewed": reviewed, "limit": limit}
        )
    )
    if search_result.get("status") != "succeeded":
        return _json_result({"status": search_result.get("status", "failed"), "search": search_result})
    candidates = search_result.get("candidates", [])
    if len(candidates) != 1:
        return _json_result(
            {
                "status": "confirmation_required" if candidates else "not_found",
                "resolution": "ambiguous_search" if candidates else "no_search_hit",
                "query": query,
                "organism": organism or None,
                "candidates": candidates,
                "search": {
                    "cache_hit": search_result.get("cache_hit"),
                    "elapsed_seconds": search_result.get("elapsed_seconds"),
                    "source": search_result.get("source"),
                },
                "next_step": "Confirm one accession before retrieving the sequence." if candidates else "Broaden or correct the protein name, function, or organism.",
            }
        )
    accession = candidates[0]["accession"]
    sequence_result = json.loads(get_uniprot_sequence.invoke({"accession": accession}))
    return _json_result(
        {
            "status": "resolved" if sequence_result.get("status") == "succeeded" else sequence_result.get("status", "failed"),
            "resolution": "unique_search_hit",
            "candidate": candidates[0],
            "result": sequence_result,
            "search": {
                "cache_hit": search_result.get("cache_hit"),
                "elapsed_seconds": search_result.get("elapsed_seconds"),
                "source": search_result.get("source"),
            },
        }
    )


@tool
def run_openfold_accession(accession: str, ligands_json: str = "") -> str:
    """Run UniProt accession -> MSA -> registered OpenFold3 execution node.

    The accession must be explicit. This tool does not guess from a protein
    name; an upstream resolver must disambiguate names before calling it.
    Optional ligands_json is an inline JSON list of bounded ligand objects. Each
    object has one chain_id and exactly one of smiles, ccd_codes, or sdf_path;
    file paths are intentionally not accepted through the agent surface.
    """
    accession = accession.strip()
    if not ACCESSION_RE.fullmatch(accession):
        return json.dumps({"status": "invalid_input", "reason": "an explicit UniProt accession is required"}, ensure_ascii=False, indent=2)
    command = [str(PYTHON), str(ROOT / "openfold3/predict_accession.py"), "--accession", accession]
    ligands: list[dict] = []
    if ligands_json.strip():
        try:
            ligands = _parse_ligands_json(ligands_json)
        except (json.JSONDecodeError, TypeError, ValueError) as exc:
            return json.dumps({"status": "invalid_input", "reason": f"invalid ligands_json: {exc}"}, ensure_ascii=False, indent=2)
        command.extend(["--ligands-json", json.dumps(ligands, ensure_ascii=False, separators=(",", ":"))])
    if OPENFOLD_WORKER_URL:
        try:
            payload = _run_openfold_worker(accession, ligands if ligands_json.strip() else [])
        except Exception as exc:  # worker failures must become a structured tool result
            payload = {"status": "failed", "stage": "worker", "accession": accession, "reason": f"{type(exc).__name__}: {exc}"}
        payload.setdefault("accession", accession)
        return json.dumps(payload, ensure_ascii=False, indent=2)
    if not (ROOT / ".tools/bin/micromamba").is_file() or not (ROOT / ".tools/openfold3-env").is_dir():
        return json.dumps(
            {
                "status": "unavailable",
                "stage": "execution",
                "accession": accession,
                "reason": "OpenFold3 execution node is not configured; the API host has no local .tools/openfold3-env. Configure OPENFOLD_WORKER_URL for the Mac worker.",
            },
            ensure_ascii=False,
            indent=2,
        )
    try:
        completed = (current_control().run if current_control() else subprocess.run)(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=3600)
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return json.dumps({"status": "timed_out", "accession": accession, "output_tail": output[-3000:]}, ensure_ascii=False, indent=2)
    output = completed.stdout.strip()
    try:
        payload = json.loads(output)
    except json.JSONDecodeError:
        payload = {"status": "failed", "accession": accession, "exit_code": completed.returncode, "output_tail": output[-3000:]}
    payload.setdefault("accession", accession)
    payload["exit_code"] = completed.returncode
    return json.dumps(payload, ensure_ascii=False, indent=2)


@tool
def run_openfold_sequence(sequence: str, ligands_json: str = "") -> str:
    """Run a pasted protein sequence through MSA and the registered OpenFold3 worker.

    This is the explicit-sequence counterpart to run_openfold_accession. It is
    deliberately bounded to one protein sequence (20-4000 standard amino acids)
    and never treats a pasted sequence as a UniProt record. Preserve all user
    supplied ligands in ligands_json; an empty list explicitly means a
    protein-only prediction and must be reported as such.
    """
    try:
        normalized = _normalize_inline_sequence(sequence)
        ligands = _parse_ligands_json(ligands_json)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        return _json_result({"status": "invalid_input", "reason": str(exc)})
    sequence_hash = hashlib.sha256(normalized.encode("ascii")).hexdigest()
    sequence_id = f"sequence-{sequence_hash[:12]}"
    if OPENFOLD_WORKER_URL:
        try:
            payload = _run_openfold_worker(None, ligands, sequence=normalized, sequence_id=sequence_id)
        except Exception as exc:  # worker failures must remain structured
            payload = {"status": "failed", "stage": "worker", "reason": f"{type(exc).__name__}: {exc}"}
        payload.update({"source_kind": "inline_sequence", "sequence_length": len(normalized), "sequence_sha256": sequence_hash})
        return _json_result(payload)

    if not (ROOT / ".tools/bin/micromamba").is_file() or not (ROOT / ".tools/openfold3-env").is_dir():
        return _json_result({
            "status": "unavailable",
            "stage": "execution",
            "source_kind": "inline_sequence",
            "sequence_length": len(normalized),
            "sequence_sha256": sequence_hash,
            "reason": "OpenFold3 execution node is not configured; set OPENFOLD_WORKER_URL for the Mac worker.",
        })
    output_root = ROOT / "results" / "openfold3-sequences" / sequence_id
    output_root.mkdir(parents=True, exist_ok=True)
    fasta_path = output_root / "input.fasta"
    fasta_path.write_text(f">{sequence_id}\n{normalized}\n", encoding="utf-8")
    command = [str(PYTHON), str(ROOT / "openfold3/predict_accession.py"), "--fasta", str(fasta_path), "--output-dir", str(output_root)]
    if ligands:
        command.extend(["--ligands-json", json.dumps(ligands, ensure_ascii=False, separators=(",", ":"))])
    try:
        completed = (current_control().run if current_control() else subprocess.run)(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=3600)
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return _json_result({"status": "timed_out", "source_kind": "inline_sequence", "sequence_length": len(normalized), "output_tail": output[-3000:]})
    try:
        payload = json.loads(completed.stdout.strip())
    except json.JSONDecodeError:
        payload = {"status": "failed", "exit_code": completed.returncode, "output_tail": completed.stdout[-3000:]}
    payload.update({"source_kind": "inline_sequence", "sequence_length": len(normalized), "sequence_sha256": sequence_hash})
    return _json_result(payload)


OPENFOLD_TOOLS = [search_uniprot, get_uniprot_sequence, resolve_uniprot_sequence, run_openfold_accession, run_openfold_sequence]
