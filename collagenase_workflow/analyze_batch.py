#!/usr/bin/env python3
"""Run independent BLASTP screening across the frozen 20-proteome benchmark."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import time
from pathlib import Path

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/collagenase-agent-20"
BLAST = ROOT / ".tools/blast-2.17.0/bin"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], log: Path, timeout: int = 900) -> dict:
    started = time.perf_counter()
    try:
        completed = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
        log.write_text(completed.stdout)
        return {"status": "succeeded" if completed.returncode == 0 else "failed", "exit_code": completed.returncode, "elapsed_seconds": round(time.perf_counter() - started, 3), "command": command, "log": str(log.relative_to(ROOT)), "output_tail": completed.stdout[-1000:]}
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        log.write_text(output)
        return {"status": "timed_out", "exit_code": None, "elapsed_seconds": round(time.perf_counter() - started, 3), "command": command, "log": str(log.relative_to(ROOT)), "output_tail": output[-1000:], "timeout_seconds": timeout}


def main() -> None:
    manifest = json.loads((OUT / "manifest.json").read_text())
    run_dir = OUT / "runs/batch-screen"
    run_dir.mkdir(parents=True, exist_ok=True)
    combined = run_dir / "all_20_proteomes.faa"
    query_meta: dict[str, dict] = {}
    with combined.open("w") as output:
        for item in manifest["proteomes"]:
            fasta = ROOT / item["fasta_path"]
            for record in SeqIO.parse(str(fasta), "fasta"):
                query_id = f"{item['proteome_id']}::{record.id}"
                query_meta[query_id] = {
                    "proteome_id": item["proteome_id"],
                    "organism": item["organism_requested"],
                    "original_id": record.id,
                    "length": len(record.seq),
                    "header": record.description,
                    "sequence": str(record.seq),
                }
                output.write(f">{query_id} {record.description}\n{record.seq}\n")

    references = ROOT / manifest["collagenase_references"]["fasta_path"]
    database = run_dir / "collagenase_refs"
    make_db = run([str(BLAST / "makeblastdb"), "-in", str(references), "-dbtype", "prot", "-parse_seqids", "-out", str(database)], run_dir / "makeblastdb.log")
    if make_db["status"] != "succeeded":
        raise RuntimeError(json.dumps(make_db))
    out_tsv = run_dir / "all-vs-collagenase.tsv"
    outfmt = "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle"
    blast = run([
        str(BLAST / "blastp"), "-query", str(combined), "-db", str(database),
        "-out", str(out_tsv), "-outfmt", outfmt, "-evalue", "1e-5",
        "-max_target_seqs", "5", "-max_hsps", "1", "-seg", "yes", "-num_threads", "6",
    ], run_dir / "blastp.log", timeout=1800)
    if blast["status"] != "succeeded":
        raise RuntimeError(json.dumps(blast))

    hits: dict[str, dict] = {}
    with out_tsv.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 10:
                continue
            qseqid, sseqid, pident, length, qlen, slen, qcovs, evalue, bitscore, stitle = row[:10]
            if qseqid not in query_meta:
                continue
            qcov = float(qcovs) if qcovs not in {"", "-"} else float(length) / float(qlen) * 100
            candidate = {
                **query_meta[qseqid],
                "qseqid": qseqid,
                "reference_id": sseqid,
                "reference_title": stitle,
                "pident": float(pident),
                "alignment_length": int(length),
                "qcov": qcov,
                "evalue": float(evalue),
                "bitscore": float(bitscore),
            }
            previous = hits.get(qseqid)
            if previous is None or (candidate["bitscore"], candidate["qcov"], candidate["pident"]) > (previous["bitscore"], previous["qcov"], previous["pident"]):
                hits[qseqid] = candidate

    candidates = []
    for qseqid, hit in hits.items():
        # This rule is fixed before the post-hoc UniProt annotation audit.
        if hit["pident"] < 25 or hit["qcov"] < 35 or hit["evalue"] > 1e-5 or not 250 <= hit["length"] <= 1800:
            continue
        if hit["pident"] >= 40 and hit["qcov"] >= 60:
            tier = "high"
        elif hit["pident"] >= 30 and hit["qcov"] >= 45:
            tier = "medium"
        else:
            tier = "low"
        annotation = hit["header"].lower()
        hit["tier"] = tier
        hit["screen_score"] = round(0.6 * hit["pident"] + 0.4 * hit["qcov"], 4)
        hit["known_collagenase_label_posthoc"] = bool(re.search(r"collagenase|collagenolytic|peptidase m9", annotation))
        candidates.append(hit)
    candidates.sort(key=lambda x: (-x["screen_score"], x["organism"], x["qseqid"]))
    for rank, item in enumerate(candidates, 1):
        item["rank"] = rank

    candidate_fasta = run_dir / "candidates.faa"
    with candidate_fasta.open("w") as output:
        for item in candidates:
            output.write(f">{item['qseqid']} organism={item['organism']} rank={item['rank']} tier={item['tier']}\n{item['sequence']}\n")
    table = run_dir / "candidate-table.tsv"
    fields = ["rank", "qseqid", "proteome_id", "organism", "original_id", "length", "reference_id", "reference_title", "pident", "qcov", "evalue", "bitscore", "screen_score", "tier", "known_collagenase_label_posthoc"]
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(candidates)
    summary = {
        "workflow": "collagenase-agent-20",
        "input_mode": "existing_protein_sequences",
        "ngs_stages_skipped": ["fastp", "SPAdes"],
        "proteome_count": len(manifest["proteomes"]),
        "total_proteins": len(query_meta),
        "proteins_with_blast_hit": len(hits),
        "candidate_count": len(candidates),
        "candidate_count_by_tier": {tier: sum(x["tier"] == tier for x in candidates) for tier in ["high", "medium", "low"]},
        "known_collagenase_posthoc_count": sum(x["known_collagenase_label_posthoc"] for x in candidates),
        "candidates_per_strain": {item["organism_requested"]: sum(x["organism"] == item["organism_requested"] for x in candidates) for item in manifest["proteomes"]},
        "database": {"path": str(references.relative_to(ROOT)), "sha256": sha256(references), "records": manifest["collagenase_references"]["records_selected"]},
        "makeblastdb": make_db,
        "blastp": blast,
        "candidate_rule": {"pident_min": 25, "qcov_min": 35, "evalue_max": 1e-5, "length_min": 250, "length_max": 1800, "posthoc_labels_excluded_from_rule": True},
        "candidate_table": str(table.relative_to(ROOT)),
        "candidate_fasta": str(candidate_fasta.relative_to(ROOT)),
        "source_sha256": {str((ROOT / item["fasta_path"]).relative_to(ROOT)): item["fasta_sha256"] for item in manifest["proteomes"]},
        "scientific_warning": "BLAST similarity and post-hoc annotation labels are evidence layers; they do not prove collagenase activity or substrate specificity.",
    }
    (run_dir / "screening.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: summary[k] for k in ["proteome_count", "total_proteins", "candidate_count", "candidate_count_by_tier", "known_collagenase_posthoc_count", "blastp"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
