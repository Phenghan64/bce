#!/usr/bin/env python3
"""Fetch and freeze an independent set of public marine microbial proteomes."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import requests
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/collagenase-agent-20"
PROTEOMES = OUT / "inputs/proteomes"
REFERENCES = OUT / "references"

# These accessions were selected from public UniProt proteome records. None is
# the Alteromonas macleodii SSi026 dataset used by the earlier workflow.
DATASETS = [
    ("UP000001870", "Alteromonas mediterranea"),
    ("UP001152447", "Pseudoalteromonas haloplanktis"),
    ("UP000006201", "Pseudoalteromonas tunicata"),
    ("UP000054316", "Vibrio alginolyticus"),
    ("UP000002493", "Vibrio parahaemolyticus"),
    ("UP000067422", "Vibrio harveyi"),
    ("UP001569200", "Vibrio splendidus"),
    ("UP000183986", "Marinobacter hydrocarbonoclasticus"),
    ("UP000007077", "Marinobacter adhaerens"),
    ("UP000001023", "Ruegeria pomeroyi"),
    ("UP000001353", "Roseobacter litoralis"),
    ("UP000830781", "Sulfitobacter pontiacus"),
    ("UP000001062", "Marinomonas mediterranea"),
    ("UP000190064", "Oceanospirillum linum"),
    ("UP000001171", "Idiomarina loihiensis"),
    ("UP000189431", "Salinivibrio costicola"),
    ("UP000000593", "Photobacterium profundum"),
    ("UP000001557", "Shewanella baltica"),
    ("UP000000547", "Colwellia psychrerythraea"),
    ("UP000008898", "Zobellia galactanivorans"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def download_proteome(session: requests.Session, upid: str, organism: str) -> dict:
    PROTEOMES.mkdir(parents=True, exist_ok=True)
    name = slug(f"{organism}_{upid}")
    compressed = PROTEOMES / f"{name}.fasta.gz"
    fasta = PROTEOMES / f"{name}.faa"
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {"format": "fasta", "compressed": "true", "query": f"proteome:{upid}"}
    if not compressed.exists():
        response = session.get(url, params=params, timeout=(20, 180))
        response.raise_for_status()
        compressed.write_bytes(response.content)
    if not fasta.exists():
        with gzip.open(compressed, "rt") as source, fasta.open("w") as target:
            target.write(source.read())
    records = list(SeqIO.parse(str(fasta), "fasta"))
    if not records:
        raise RuntimeError(f"empty proteome: {upid}")
    return {
        "proteome_id": upid,
        "organism_requested": organism,
        "source_url": response.url if 'response' in locals() else session.get(url, params=params, timeout=(20, 60)).url,
        "compressed_path": str(compressed.relative_to(ROOT)),
        "fasta_path": str(fasta.relative_to(ROOT)),
        "compressed_bytes": compressed.stat().st_size,
        "fasta_bytes": fasta.stat().st_size,
        "compressed_sha256": sha256(compressed),
        "fasta_sha256": sha256(fasta),
        "protein_count": len(records),
        "min_length": min(len(r.seq) for r in records),
        "max_length": max(len(r.seq) for r in records),
    }


def fetch_collagenase_references(session: requests.Session) -> dict:
    REFERENCES.mkdir(parents=True, exist_ok=True)
    fasta_path = REFERENCES / "external_collagenase_refs.faa"
    meta_path = REFERENCES / "external_collagenase_refs.tsv"
    response = session.get(
        "https://rest.uniprot.org/uniprotkb/search",
        params={
            "query": "protein_name:collagenase AND taxonomy_id:2",
            "format": "tsv",
            "fields": "accession,id,protein_name,organism_name,length,sequence",
            "size": 100,
        },
        timeout=(20, 180),
    )
    response.raise_for_status()
    rows = list(csv.DictReader(response.text.splitlines(), delimiter="\t"))
    target_names = {name.lower() for _, name in DATASETS}
    selected = []
    seen = set()
    for row in rows:
        sequence = (row.get("Sequence") or "").strip().upper()
        organism = row.get("Organism") or ""
        if not sequence or len(sequence) < 300 or not re.fullmatch(r"[A-Z]+", sequence):
            continue
        if any(name in organism.lower() for name in target_names):
            continue
        if sequence in seen:
            continue
        seen.add(sequence)
        selected.append(row)
    if len(selected) < 5:
        raise RuntimeError(f"too few external collagenase references after leakage filter: {len(selected)}")
    with fasta_path.open("w") as fasta, meta_path.open("w", newline="") as meta:
        writer = csv.DictWriter(meta, fieldnames=["Entry", "Entry Name", "Protein names", "Organism", "Length"], delimiter="\t")
        writer.writeheader()
        for row in selected:
            fasta.write(f">{row['Entry']} organism={row['Organism']} protein={row['Protein names']}\n")
            seq = row["Sequence"].strip().upper()
            for i in range(0, len(seq), 80):
                fasta.write(seq[i:i + 80] + "\n")
            writer.writerow({key: row.get(key, "") for key in writer.fieldnames})
    return {
        "query": "protein_name:collagenase AND taxonomy_id:2",
        "source_url": response.url,
        "records_returned": len(rows),
        "records_selected": len(selected),
        "excluded_target_organism_records": sum(any(name in (r.get("Organism") or "").lower() for name in target_names) for r in rows),
        "fasta_path": str(fasta_path.relative_to(ROOT)),
        "fasta_sha256": sha256(fasta_path),
        "metadata_path": str(meta_path.relative_to(ROOT)),
    }


def create_direct_sequence_fixture(manifest: dict) -> dict:
    fixture = OUT / "inputs/existing_protein_sequences.faa"
    fixture.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for item in manifest["proteomes"][:3]:
        path = ROOT / item["fasta_path"]
        record = next(SeqIO.parse(str(path), "fasta"))
        records.append((f"{item['proteome_id']}|{record.id}|existing_sequence", str(record.seq)))
    with fixture.open("w") as handle:
        for ident, sequence in records:
            handle.write(f">{ident}\n{sequence}\n")
    return {"path": str(fixture.relative_to(ROOT)), "sequence_count": len(records), "sha256": sha256(fixture)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    proteomes = [download_proteome(session, upid, organism) for upid, organism in DATASETS]
    refs = fetch_collagenase_references(session)
    manifest = {
        "benchmark": "collagenase-agent-20",
        "objective": "find potential bacterial collagenases in 20 previously unused marine microbial proteomes",
        "data_policy": "public UniProt proteomes; external collagenase references exclude the 20 target organism names",
        "ngs_route": {"status": "not_required", "reason": "existing protein sequences were supplied; fastp and SPAdes must be skipped"},
        "proteomes": proteomes,
        "collagenase_references": refs,
    }
    manifest["direct_sequence_fixture"] = create_direct_sequence_fixture(manifest)
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"proteomes": len(proteomes), "total_proteins": sum(x["protein_count"] for x in proteomes), "references": refs["records_selected"], "manifest": str(OUT / "manifest.json")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
