#!/usr/bin/env python3
"""Query the InterPro Matches API twice for the current candidate FASTA."""
from __future__ import annotations

import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
API = "https://www.ebi.ac.uk/interpro/matches/api"


def digest(value: object) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def fetch_json(url: str, payload: dict | None = None, attempts: int = 3) -> tuple[dict, dict]:
    last_error = None
    for attempt in range(1, attempts + 1):
        started = time.perf_counter()
        try:
            body = None if payload is None else json.dumps(payload).encode()
            request = urllib.request.Request(url, data=body, headers={"content-type": "application/json"} if body else {})
            with urllib.request.urlopen(request, timeout=120) as response:
                result = json.loads(response.read().decode())
            return result, {"elapsed_seconds": round(time.perf_counter() - started, 3), "attempt": attempt}
        except Exception as exc:  # network/API failures are reported to the caller
            last_error = exc
            if attempt < attempts:
                time.sleep(2 ** (attempt - 1))
    raise RuntimeError(f"InterPro API request failed after {attempts} attempts: {last_error}")


def md5(sequence: str) -> str:
    return hashlib.md5(sequence.upper().encode()).hexdigest().upper()


def parse_fasta(path: Path) -> list[dict]:
    records = []
    for record in SeqIO.parse(str(path), "fasta"):
        sequence = str(record.seq).upper()
        if not sequence or "*" in sequence or any(char not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" for char in sequence):
            raise ValueError(f"invalid protein sequence: {record.id}")
        records.append({"id": record.id, "sequence": sequence})
    if not records:
        raise ValueError("no FASTA records")
    return records


def summarize(record: dict, item: dict) -> dict:
    matches = item.get("matches") or []
    pfam, ncbifam, interpro, libraries = set(), set(), set(), set()
    for match in matches:
        signature = match.get("signature") or {}
        release = signature.get("signatureLibraryRelease") or {}
        library = release.get("library")
        if library:
            libraries.add(library)
        accession = match.get("model_ac") or signature.get("accession")
        if library == "Pfam" and accession:
            pfam.add(accession)
        if library == "NCBIFAM" and accession:
            ncbifam.add(accession)
        entry = signature.get("entry") or {}
        if entry.get("accession"):
            interpro.add(entry["accession"])
    return {
        "id": record["id"],
        "family": None,
        "length": len(record["sequence"]),
        "md5": md5(record["sequence"]),
        "found": bool(item.get("found") or matches),
        "match_count": len(matches),
        "libraries": sorted(libraries),
        "pfam": sorted(pfam),
        "ncbifam": sorted(ncbifam),
        "interpro": sorted(interpro),
    }


def lookup(records: list[dict]) -> dict:
    started = time.perf_counter()
    by_md5 = {md5(record["sequence"]): record for record in records}
    md5s = sorted(by_md5)
    summaries = []
    requests = []
    for offset in range(0, len(md5s), 100):
        chunk = md5s[offset : offset + 100]
        result, request_meta = fetch_json(f"{API}/matches", {"md5": chunk})
        requests.append({"offset": offset, "count": len(chunk), **request_meta})
        for item in result.get("results", []):
            record = by_md5.get(str(item.get("md5", "")).upper())
            if record:
                summaries.append(summarize(record, item))
    seen = {item["md5"] for item in summaries}
    for record in records:
        if md5(record["sequence"]) not in seen:
            summaries.append(summarize(record, {"matches": []}))
    summaries.sort(key=lambda item: item["id"])
    return {
        "sequence_count": len(records),
        "unique_md5_count": len(md5s),
        "found_count": sum(item["found"] and item["match_count"] > 0 for item in summaries),
        "total_match_count": sum(item["match_count"] for item in summaries),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "requests": requests,
        "summaries": summaries,
        "normalized_sha256": digest(summaries),
    }


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: run_interpro_candidates.py INPUT.faa OUTPUT.json", file=sys.stderr)
        return 2
    input_path, output_path = (Path(value).resolve() for value in sys.argv[1:])
    records = parse_fasta(input_path)
    info, info_meta = fetch_json(f"{API}/info")
    first = lookup(records)
    second = lookup(records)
    output = {
        "tool": "InterProScan 6 Matches API path",
        "workflow_version": "6.0.3.0",
        "api_url": API,
        "api_info": info,
        "input": {"path": str(input_path.relative_to(ROOT)), "sequence_count": len(records), "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest()},
        "run": first,
        "repeat": {
            "elapsed_seconds": second["elapsed_seconds"],
            "normalized_sha256": second["normalized_sha256"],
            "byte_level_normalized_match": first["normalized_sha256"] == second["normalized_sha256"],
            "requests": second["requests"],
        },
        "info_request": info_meta,
    }
    output_path.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"sequences": len(records), "found": first["found_count"], "matches": first["total_match_count"], "elapsed_seconds": first["elapsed_seconds"], "repeat_identical": output["repeat"]["byte_level_normalized_match"], "output": str(output_path.relative_to(ROOT))}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
