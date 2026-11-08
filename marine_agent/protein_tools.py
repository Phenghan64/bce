"""General protein inputs for the Agent; no arbitrary paths or benchmark fallback."""
from __future__ import annotations

import json
import os
import uuid
from typing import Literal
from urllib.request import Request, build_opener, ProxyHandler

from langchain_core.tools import tool

from .access_control import grant_analysis, require_analysis
from .execution_control import check_cancelled, current_control
from .protein_screening import ANALYSIS_ID, ROOT, read_analysis, run_screen, sha, validate_request
from .vina_remote import _safe_extract_bundle
from .domain_evidence import annotate, read_domain_analysis


def execute_screen(payload: dict, *, emit=None) -> dict:
    # Validate on the control plane as well as the execution node.
    request, _, _ = validate_request(payload)
    url = os.getenv("PROTEIN_WORKER_URL", os.getenv("VINA_WORKER_URL", "")).rstrip("/")
    if not url:
        result = run_screen(request.model_dump(), emit=emit)
        if result.get('status') == 'succeeded':
            grant_analysis('protein-screenings', result['analysis_id'])
        return result
    control = current_control()
    body = {"job_id": control.job_id if control else uuid.uuid4().hex, "request": request.model_dump()}
    headers = {"Content-Type": "application/json", "X-Vina-Worker-Token": os.getenv("PROTEIN_WORKER_TOKEN", os.getenv("VINA_WORKER_TOKEN", ""))}
    check_cancelled()
    client = build_opener(ProxyHandler({}))
    with client.open(Request(url + "/v1/proteins/screen", data=json.dumps(body).encode(), headers=headers, method="POST"), timeout=4000) as response:
        reply = json.load(response)
    check_cancelled()
    result = reply.get("result") or {}
    if result.get("status") != "succeeded":
        raise RuntimeError(reply.get("reason") or result.get("reason") or "protein execution node did not complete the analysis")
    analysis_id = result.get("analysis_id", "")
    if not ANALYSIS_ID.fullmatch(analysis_id):
        raise RuntimeError("execution node returned an invalid analysis ID")
    expected = f"results/protein-screenings/{analysis_id}"
    bundle = reply.get("bundle") or {}
    if bundle.get("root") != expected:
        raise RuntimeError("protein artifact root does not match the analysis ID")
    _safe_extract_bundle(bundle.get("data"), expected)
    manifest = json.loads((ROOT / expected / "manifest.json").read_text())
    recorded = manifest.get("spec", {})
    if recorded.get("inputs") != {"query": sha(request.query_fasta.encode()), "reference": sha(request.reference_fasta.encode())} or recorded.get("request") != request.model_dump(exclude={"query_fasta", "reference_fasta"}):
        raise RuntimeError("protein artifacts do not match submitted inputs and parameters")
    verified = read_analysis(analysis_id)
    if verified.get("status") != "succeeded":
        raise RuntimeError("protein artifact integrity verification failed")
    grant_analysis('protein-screenings', analysis_id)
    return verified


@tool
def run_protein_homology_screen(
    query_fasta: str, reference_fasta: str, query_source: str, reference_source: str,
    objective: str, min_identity: float = 30, min_query_coverage: float = 50,
    min_reference_coverage: float = 50, max_evalue: float = 1e-5,
    threads: int = 6, max_target_seqs: int = 25,
) -> str:
    """Screen NEW protein FASTA against an explicit target reference FASTA with real BLASTP.

    Both sources and the objective are required. No old benchmark is selected.
    Existing proteins skip NGS/assembly/annotation. Choose task-specific coverage,
    identity and e-value thresholds; defaults are not universal enzyme criteria.
    Outputs include analysis_id, ranked homology candidates, source/input hashes,
    versions, commands, QC and restartable stage artifacts. Repeating EXACT input
    reuses validated stages. Domain/activity validation is a separate next step.
    Accepts up to 32 MiB per FASTA; large inputs need the registered data backend.
    """
    try:
        payload = dict(query_fasta=query_fasta, reference_fasta=reference_fasta,
                       query_source=query_source, reference_source=reference_source,
                       objective=objective, min_identity=min_identity,
                       min_query_coverage=min_query_coverage,
                       min_reference_coverage=min_reference_coverage,
                       max_evalue=max_evalue, threads=threads, max_target_seqs=max_target_seqs)
        return json.dumps(execute_screen(payload), ensure_ascii=False)
    except Exception as exc:
        return json.dumps({"status": "failed", "reason": str(exc), "activity_confirmed": False}, ensure_ascii=False)


@tool
def inspect_protein_analysis(analysis_id: str, protein_id: str = "") -> str:
    """Read a new protein screening's verified evidence and stage history.

    Set protein_id to retrieve one passing candidate's exact sequence for domain
    verification or run_openfold_sequence. A BLAST candidate is not confirmed
    enzyme activity. Never substitute records from the frozen test datasets.
    """
    try:
        require_analysis('protein-screenings', analysis_id)
        return json.dumps(read_analysis(analysis_id, protein_id=protein_id), ensure_ascii=False)
    except (ValueError, OSError) as exc:
        return json.dumps({"status": "unavailable", "reason": str(exc)}, ensure_ascii=False)


@tool
def annotate_protein_domains(query_fasta: str = "", source: str = "", analysis_id: str = "",
                             protein_id: str = "", expected_pfam: list[str] | None = None,
                             refresh: bool = False, method: Literal["lookup", "full_scan"] = "lookup") -> str:
    """Get exact-sequence InterPro/Pfam evidence for NEW protein input.

    Provide either query_fasta with source OR a prior protein screening
    analysis_id (optionally select protein_id). Saved candidate sequences are
    read and hash-verified directly, avoiding model copying. expected_pfam is
    optional, explicit PF00000 accessions to check against observed matches.
    Returns real signature names, coordinates, versions and source artifacts.
    method="lookup" queries existing exact-sequence records quickly.
    method="full_scan" runs all official default available InterProScan analyses
    on Linux with the Matches API disabled; use for new/unknown sequences or an
    explicit request to scan. Progress and validated Nextflow resume are saved.
    No match does not establish absence or lack of activity.
    Repeating identical input reuses validated chunks for the same release;
    refresh=True re-queries. It does not verify catalytic sites or enzyme activity.
    """
    try:
        return json.dumps(annotate({"query_fasta": query_fasta, "source": source,
            "analysis_id": analysis_id, "protein_id": protein_id,
            "expected_pfam": expected_pfam or [], "refresh": refresh, "method": method}), ensure_ascii=False)
    except Exception as exc:
        return json.dumps({"status": "failed", "kind": "protein_domain_evidence", "reason": str(exc)}, ensure_ascii=False)


@tool
def inspect_domain_evidence(analysis_id: str, protein_id: str = "") -> str:
    """Read hash-verified domain evidence and exact signature names from a saved domain analysis.

    analysis_id refers to a domain lookup result, not a BLAST analysis. Optional
    protein_id selects one record. Unknown precomputed records need full scans.
    """
    try:
        return json.dumps(read_domain_analysis(analysis_id, protein_id=protein_id), ensure_ascii=False)
    except (ValueError, OSError) as exc:
        return json.dumps({"status": "unavailable", "reason": str(exc)}, ensure_ascii=False)


PROTEIN_TOOLS = [run_protein_homology_screen, inspect_protein_analysis,
                 annotate_protein_domains, inspect_domain_evidence]
