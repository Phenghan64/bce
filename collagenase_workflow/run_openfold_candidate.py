#!/usr/bin/env python3
"""Run a gated, MSA-backed OpenFold3 prediction for one screened candidate."""
from __future__ import annotations

import csv
import json
import re
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/collagenase-agent-20"
RUN_DIR = OUT / "runs/batch-screen"


def execute(command: list[str], log: Path, timeout: int) -> dict:
    started = time.perf_counter()
    try:
        proc = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
        log.write_text(proc.stdout)
        return {"status": "succeeded" if proc.returncode == 0 else "failed", "exit_code": proc.returncode, "elapsed_seconds": round(time.perf_counter() - started, 3), "command": command, "log": str(log.relative_to(ROOT)), "output_tail": proc.stdout[-1200:]}
    except subprocess.TimeoutExpired as exc:
        text = exc.stdout or ""
        if isinstance(text, bytes):
            text = text.decode(errors="replace")
        log.write_text(text)
        return {"status": "timed_out", "exit_code": None, "elapsed_seconds": round(time.perf_counter() - started, 3), "command": command, "log": str(log.relative_to(ROOT)), "output_tail": text[-1200:], "timeout_seconds": timeout}


def main() -> None:
    candidate_id = None
    if len(__import__("sys").argv) > 1:
        candidate_id = __import__("sys").argv[1]
    screening = json.loads((RUN_DIR / "screening.json").read_text())
    with (RUN_DIR / "candidate-table.tsv").open() as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row["qseqid"]: row for row in candidates}
    if candidate_id is None:
        eligible = [row for row in candidates if row["tier"] == "high"]
        candidate_id = eligible[0]["qseqid"] if eligible else candidates[0]["qseqid"]
    if candidate_id not in by_id:
        print(json.dumps({"status": "not_found", "candidate_id": candidate_id}, ensure_ascii=False))
        raise SystemExit(2)
    candidate = by_id[candidate_id]
    interpro = json.loads((RUN_DIR / "interpro-candidates.json").read_text())
    ipr = next((x for x in interpro["run"]["summaries"] if x["id"] == candidate_id), None)
    if not ipr or "PF01752" not in ipr.get("pfam", []):
        print(json.dumps({"status": "blocked", "candidate_id": candidate_id, "reason": "candidate lacks the collagenase-domain gate PF01752"}, ensure_ascii=False))
        raise SystemExit(3)
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id)
    target = OUT / "runs/openfold" / safe_id
    target.mkdir(parents=True, exist_ok=True)
    query = {
        "seeds": [42],
        "queries": {
            safe_id: {
                "query_name": safe_id,
                "chains": [{"molecule_type": "PROTEIN", "chain_ids": ["A"], "description": candidate["organism"], "sequence": _sequence_for(candidate_id), "non_canonical_residues": None, "smiles": None, "ligand_name": None, "ccd_codes": None, "paired_msa_file_paths": None, "main_msa_file_paths": None, "template_alignment_file_path": None, "template_entry_chain_ids": None, "template_cif_paths": None, "template_cif_chain_ids": None, "sdf_file_path": None, "cyclic": False}],
                "use_msas": True,
                "use_paired_msas": False,
                "use_main_msas": True,
                "covalent_bonds": None,
                "pocket_constraint": None,
            }
        },
    }
    query_path = target / "query.json"
    query_path.write_text(json.dumps(query, indent=2) + "\n")
    msa_dir = target / "msa"
    msa = execute([str(ROOT / "openfold3/run.sh"), "python", str(ROOT / "openfold3/fetch_msa.py"), str(query_path), str(msa_dir)], target / "msa.log", timeout=900)
    if msa["status"] != "succeeded":
        result = {"status": "failed", "candidate_id": candidate_id, "candidate": candidate, "interpro": ipr, "msa": msa, "reason": "MSA preparation failed; prediction was not attempted"}
        (target / "structure.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
        print(json.dumps(result, ensure_ascii=False, indent=2))
        raise SystemExit(1)
    query_msa = msa_dir / "query_msa.json"
    prediction_dir = target / "prediction"
    benchmark = target / "benchmark.json"
    prediction = execute([str(ROOT / "openfold3/run.sh"), "python", "-u", str(ROOT / "openfold3/benchmark.py"), str(benchmark), "predict", "--query-json", str(query_msa), "--runner-yaml", str(ROOT / "openfold3/runner.yml"), "--num-diffusion-samples", "1", "--use-msa-server", "false", "--use-templates", "false", "--output-dir", str(prediction_dir)], target / "prediction.log", timeout=1800)
    confidences = list(prediction_dir.glob(f"{safe_id}/seed_*/*confidences_aggregated.json"))
    confidence = json.loads(confidences[0].read_text()) if confidences else None
    if prediction["status"] != "succeeded" or confidence is None:
        status = "failed"
        decision = "structure_unavailable"
    else:
        status = "succeeded"
        decision = "structure_review_allowed" if confidence.get("avg_plddt", 0) >= 70 and confidence.get("ptm", 0) >= 0.5 else "low_confidence_structure"
    result = {"status": status, "candidate_id": candidate_id, "candidate": candidate, "interpro": ipr, "msa": msa, "prediction": prediction, "confidence": confidence, "decision": decision, "scientific_limit": "OpenFold3 confidence supports fold/domain review; it does not establish collagenase activity or substrate specificity."}
    (target / "structure.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": status, "candidate_id": candidate_id, "decision": decision, "confidence": confidence, "elapsed_seconds": round(msa["elapsed_seconds"] + prediction["elapsed_seconds"], 3)}, ensure_ascii=False, indent=2))
    if status != "succeeded":
        raise SystemExit(1)


def _sequence_for(candidate_id: str) -> str:
    with (RUN_DIR / "candidates.faa").open() as handle:
        current = None
        seq = []
        for line in handle:
            if line.startswith(">"):
                if current == candidate_id:
                    return "".join(seq)
                current = line[1:].split()[0]
                seq = []
            elif current == candidate_id:
                seq.append(line.strip())
        if current == candidate_id:
            return "".join(seq)
    raise KeyError(candidate_id)


if __name__ == "__main__":
    main()
