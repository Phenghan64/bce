#!/usr/bin/env python3
"""Run cheap, non-destructive failure/recovery checks for the marine workflow."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/marine-microbe-workflow/runs/recovery-tests"
OUT.mkdir(parents=True, exist_ok=True)
RUN_TAG = str(int(time.time()))


def run_case(name: str, command: list[str], *, env: dict[str, str] | None = None, success_markers: list[str] | None = None) -> dict:
    started = time.perf_counter()
    proc = subprocess.run(
        command,
        cwd=ROOT,
        env={**os.environ, **(env or {})},
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=120,
    )
    log = OUT / f"{name}.log"
    log.write_text(proc.stdout)
    marker_hits = {marker: marker in proc.stdout for marker in (success_markers or [])}
    semantic_failure_detected = bool(marker_hits) and any(marker_hits.values())
    return {
        "name": name,
        "command": command,
        "exit_code": proc.returncode,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "log": str(log.relative_to(ROOT)),
        # Some tools report malformed input through a structured/log marker but
        # still exit 0. That is an expected failure once the marker is present.
        "passed_expected_failure": proc.returncode != 0 or semantic_failure_detected,
        "semantic_failure_detected": semantic_failure_detected,
        "markers": marker_hits,
        "tail": proc.stdout[-500:].replace("\n", " "),
    }


def main() -> None:
    cases: list[dict] = []
    bad = OUT / "bad.fastq"
    bad.write_text("@r1\nACGT\n+\n!!!\n")
    cases.append(run_case("invalid_fastq", [
        str(ROOT / ".tools/fastp-1.4.0/bin/fastp"),
        "-i", str(bad), "-o", str(OUT / "bad.out.fastq"), "--disable_adapter_trimming",
    ], success_markers=["ERROR: sequence and quality have different length", "Your FASTQ may be invalid"]))

    r1 = OUT / "mismatch.R1.fastq"
    r2 = OUT / "mismatch.R2.fastq"
    r1.write_text("@r1/1\nACGT\n+\n!!!!\n@r2/1\nACGT\n+\n!!!!\n")
    r2.write_text("@r1/2\nACGT\n+\n!!!!\n")
    cases.append(run_case("paired_read_mismatch", [
        str(ROOT / ".tools/fastp-1.4.0/bin/fastp"), "-i", str(r1), "-I", str(r2),
        "-o", str(OUT / "mismatch.R1.out.fastq"), "-O", str(OUT / "mismatch.R2.out.fastq"),
    ], success_markers=["WARNING: different read numbers", "Ignore the unmatched reads"]))

    cases.append(run_case("missing_spades_input", [
        str(ROOT / ".tools/spades-env/bin/python3"), str(ROOT / ".tools/spades-env/bin/spades.py"),
        "--isolate", "--pe1-1", str(OUT / "does-not-exist.R1.fastq.gz"),
        "--pe1-2", str(OUT / "does-not-exist.R2.fastq.gz"), "-o", str(OUT / "spades-missing"),
    ]))

    cases.append(run_case("missing_bakta_db", [
        str(ROOT / ".tools/bakta-env/bin/bakta"), "--db", str(OUT / "missing-bakta-db"),
        "--output", str(OUT / f"bakta-missing-{RUN_TAG}"), str(ROOT / "results/marine-microbe-workflow/runs/spades-full2/contigs.fasta"),
    ]))

    cases.append(run_case("missing_trnascan_path", [
        str(ROOT / ".tools/bakta-env/bin/bakta"), "--db", str(ROOT / "results/bakta-validation/db-light"),
        "--output", str(OUT / f"bakta-no-trnascan-{RUN_TAG}"), "--prefix", "SSi026", "--meta", "--threads", "2",
        str(ROOT / "results/marine-microbe-workflow/runs/spades-full2/contigs.fasta"),
    ], env={"PATH": "/usr/bin:/bin"}, success_markers=["dependency not found", "tRNAscan-SE"]))

    cases.append(run_case("missing_blast_db", [
        str(ROOT / ".tools/blast-2.17.0/bin/blastp"), "-query", str(ROOT / "results/marine-microbe-workflow/runs/bakta-light-run3/SSi026.faa"),
        "-db", str(OUT / "no-such-db"), "-outfmt", "6", "-out", str(OUT / "blast.out"),
    ]))

    invalid_fasta = OUT / "invalid.faa"
    invalid_fasta.write_text(">bad\nACD?EF\n")
    cases.append(run_case("invalid_interpro_fasta", [
        os.environ.get("CODEX_MCP_NODE_PATH", "node"), str(ROOT / "results/run_interproscan6_api_validation.mjs"),
        str(invalid_fasta), str(OUT / "invalid-interpro.json"),
    ]))

    timeout_script = OUT / "interpro-timeout.mjs"
    timeout_script.write_text(
        "const c = new AbortController(); setTimeout(() => c.abort(), 100); "
        "fetch('http://127.0.0.1:9/matches', {signal:c.signal}).then(() => process.exit(2))"
        ".catch(e => { console.error(e.name); process.exit(e.name === 'AbortError' ? 1 : 2); });\n"
    )
    cases.append(run_case("interpro_api_timeout_simulated", [
        os.environ.get("CODEX_MCP_NODE_PATH", "node"), str(timeout_script),
    ], success_markers=["AbortError", "TypeError"]))

    malformed_query = OUT / "malformed-openfold.json"
    malformed_query.write_text(json.dumps({"seeds": [42], "queries": {"bad": {"query_name": "bad", "chains": []}}}))
    cases.append(run_case("malformed_openfold_query", [
        str(ROOT / "openfold3/run.sh"), "python", "-u", str(ROOT / "openfold3/benchmark.py"),
        str(OUT / "openfold-malformed.json"), "predict", "--query-json", str(malformed_query),
        "--runner-yaml", str(ROOT / "openfold3/runner.yml"), "--num-diffusion-samples", "1",
        "--use-msa-server", "false", "--use-templates", "false", "--output-dir", str(OUT / f"openfold-malformed-out-{RUN_TAG}"),
    ], success_markers=["Failed Queries: bad", "prediction returned None"]))

    # A completed fastp and SPAdes result is intentionally reused; this is the
    # positive continuation check after the expected failures above.
    continuation = {
        "name": "resume_from_completed_fastp_spades",
        "fastp_clean_r1_exists": (ROOT / "results/marine-microbe-workflow/runs/fastp/SRR32944547.clean.R1.fastq.gz").is_file(),
        "fastp_clean_r2_exists": (ROOT / "results/marine-microbe-workflow/runs/fastp/SRR32944547.clean.R2.fastq.gz").is_file(),
        "spades_contigs_exists": (ROOT / "results/marine-microbe-workflow/runs/spades-full2/contigs.fasta").is_file(),
    }
    continuation["passed"] = all(v for k, v in continuation.items() if k != "name")
    cases.append(continuation)

    report = {"tool": "marine_microbe_workflow recovery checks", "cases": cases}
    (OUT / "recovery-tests.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
