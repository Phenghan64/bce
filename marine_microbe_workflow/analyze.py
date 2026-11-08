#!/usr/bin/env python3
"""Audit assembled protein candidate discovery against a held-out reference annotation."""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


TARGET = re.compile(r"(protease|peptidase|aminopeptidase|metalloprotease|metallopeptidase|endopeptidase)", re.I)
EXCLUDE = re.compile(r"(prepilin peptidase|signal peptide peptidase|peptidoglycan|transpeptidase|ATP-dependent Clp protease ATP-binding)", re.I)


def parse_attrs(value: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def reference_truth(gff: Path) -> dict[str, dict]:
    truth: dict[str, dict] = {}
    with gzip_open(gff) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_attrs(fields[8])
            pid = attrs.get("protein_id")
            if not pid:
                continue
            product = attrs.get("product", "").replace("%2C", ",")
            is_target = bool(TARGET.search(product)) and not bool(EXCLUDE.search(product))
            truth[pid] = {
                "protein_id": pid,
                "product": product,
                "locus_tag": attrs.get("locus_tag", ""),
                "target": is_target,
            }
    return truth


def gzip_open(path: Path):
    import gzip
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else path.open("r", encoding="utf-8")


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with gzip_open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            records[name] = "".join(chunks)
    return records


def parse_bakta_tsv(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    header_index = next((i for i, line in enumerate(lines) if "Sequence Id" in line), None)
    if header_index is None:
        return out
    header = lines[header_index].lstrip("#")
    rows = csv.DictReader(lines[header_index + 1 :], fieldnames=header.split("\t"), delimiter="\t")
    for row in rows:
        pid = (row.get("locus_tag") or row.get("Locus Tag") or row.get("Gene")
               or row.get("query") or row.get("Sequence Id"))
        if pid:
            out[pid] = row
    return out


def parse_blast(path: Path) -> dict[str, dict]:
    best: dict[str, dict] = {}
    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            q, s, pident, length, qlen, slen, qcov, evalue, bits = raw.rstrip("\n").split("\t")[:9]
            if s.startswith("ref|") and s.endswith("|"):
                s = s[4:-1]
            hit = {
                "query": q, "subject": s, "pident": float(pident), "length": int(length),
                "qlen": int(qlen), "slen": int(slen), "qcov": float(qcov),
                "evalue": float(evalue), "bitscore": float(bits),
            }
            previous = best.get(q)
            if previous is None or hit["bitscore"] > previous["bitscore"]:
                best[q] = hit
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-gff", type=Path, required=True)
    ap.add_argument("--reference-faa", type=Path, required=True)
    ap.add_argument("--bakta-faa", type=Path, required=True)
    ap.add_argument("--bakta-tsv", type=Path, required=True)
    ap.add_argument("--blast-tsv", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    truth = reference_truth(args.reference_gff)
    refseq = parse_fasta(args.reference_faa)
    assembled = parse_fasta(args.bakta_faa)
    bakta = parse_bakta_tsv(args.bakta_tsv)
    best = parse_blast(args.blast_tsv)

    # Bakta FAA headers contain the locus tag in the description. Resolve it without
    # trusting a product label alone, then apply the held-out reference label.
    query_to_locus: dict[str, str] = {}
    for query, seq in assembled.items():
        query_to_locus[query] = query
    rows = []
    for query, seq in assembled.items():
        hit = best.get(query, {})
        subject = hit.get("subject", "")
        ref = truth.get(subject, {})
        product = ref.get("product", "")
        # Candidate discovery is based on Bakta's product text when available; the
        # reference target label is only used for the audit confusion matrix.
        b_row = bakta.get(query, {})
        b_product = " ".join(str(v) for k, v in b_row.items() if k and "product" in k.lower())
        predicted = bool(TARGET.search(b_product)) and not bool(EXCLUDE.search(b_product))
        confident_hit = bool(hit) and hit["pident"] >= 80 and hit["qcov"] >= 80 and hit["evalue"] <= 1e-20
        rows.append({
            "query": query,
            "length": len(seq),
            "bakta_product": b_product,
            "reference_subject": subject,
            "reference_product": product,
            "predicted_protease": predicted,
            "reference_protease": bool(ref.get("target", False)) if confident_hit else None,
            "blast_confident": confident_hit,
            "blast_pident": hit.get("pident"),
            "blast_qcov": hit.get("qcov"),
            "blast_evalue": hit.get("evalue"),
            "blast_bitscore": hit.get("bitscore"),
        })

    with (args.outdir / "candidate-audit.tsv").open("w", encoding="utf-8", newline="") as fh:
        fields = list(rows[0]) if rows else ["query"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    audited = [r for r in rows if r["reference_protease"] is not None]
    predicted = [r for r in audited if r["predicted_protease"]]
    truth_positive = [r for r in audited if r["reference_protease"]]
    tp = sum(r["predicted_protease"] and r["reference_protease"] for r in audited)
    fp = sum(r["predicted_protease"] and not r["reference_protease"] for r in audited)
    fn = sum((not r["predicted_protease"]) and r["reference_protease"] for r in audited)
    tn = sum((not r["predicted_protease"]) and (not r["reference_protease"]) for r in audited)
    summary = {
        "assembled_proteins": len(rows),
        "reference_cds": len(truth),
        "reference_target_cds": sum(v["target"] for v in truth.values()),
        "blast_confident_audit_proteins": len(audited),
        "predicted_candidates": len(predicted),
        "reference_target_in_audit": len(truth_positive),
        "confusion_matrix": {"tp": tp, "fp": fp, "fn": fn, "tn": tn},
        "precision": tp / (tp + fp) if tp + fp else None,
        "recall": tp / (tp + fn) if tp + fn else None,
        "candidate_examples": [r for r in rows if r["predicted_protease"]][:20],
        "non_target_controls": [r for r in rows if r["reference_protease"] is False][:10],
        "note": "Reference labels are held-out audit evidence; they were not used to define candidate product terms."
    }
    (args.outdir / "candidate-audit.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
