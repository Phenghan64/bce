#!/usr/bin/env python3
"""Run a reproducible UniProt/FASTA -> MSA -> local OpenFold3 job.

This entry point deliberately accepts an explicit accession or FASTA file.
Name resolution belongs in a separate, confirmation-aware layer.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results/openfold3-jobs"
DEFAULT_UNIPROT_CACHE = ROOT / "results/openfold3-uniprot-cache"
DEFAULT_MSA_CACHE = ROOT / "results/openfold3-msa-cache"
DEFAULT_MPS_LOCK = ROOT / "results/openfold3-runtime/openfold3-mps.lock"
AA = set("ACDEFGHIKLMNPQRSTVWY")
OPENFOLD_CHECKPOINT = "openbind-2025-06-30-174k"
MSA_CACHE_FORMAT = 2


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def safe_id(value: str) -> str:
    result = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    if not result:
        raise ValueError("identifier is empty after sanitization")
    return result[:160]


def parse_fasta(path: Path) -> tuple[str, str]:
    header = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                raise ValueError("FASTA input must contain exactly one record")
            header = line[1:].strip() or path.stem
        else:
            chunks.append(line.replace(" ", "").upper())
    if header is None:
        raise ValueError("FASTA input has no header")
    return header.split()[0], "".join(chunks)


def validate_sequence(sequence: str) -> str:
    normalized = re.sub(r"\s+", "", sequence).upper()
    if not normalized:
        raise ValueError("protein sequence is empty")
    invalid = sorted(set(normalized) - AA)
    if invalid:
        raise ValueError(f"sequence contains unsupported residues: {','.join(invalid)}")
    if len(normalized) < 20:
        raise ValueError("protein sequence is shorter than 20 residues")
    return normalized


def atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def _uniprot_metadata(payload: dict, accession: str, url: str, *, cache_hit: bool, cache_path: Path | None,
                      retrieved_at: str | None = None, cache_created_at: str | None = None,
                      attempts: int | None = None) -> tuple[dict, str]:
    sequence = validate_sequence(payload.get("sequence", {}).get("value", ""))
    primary = payload.get("primaryAccession", accession)
    protein = payload.get("proteinDescription", {})
    recommended = protein.get("recommendedName", {}).get("fullName", {}).get("value")
    organism = payload.get("organism", {}).get("scientificName")
    metadata = {
        "accession_requested": accession,
        "primary_accession": primary,
        "entry_version": payload.get("entryVersion"),
        "dataset": payload.get("entryType") or payload.get("dataset"),
        "protein_name": recommended,
        "organism": organism,
        "length": len(sequence),
        "fragment": str(protein.get("flag", "")).lower() in {"fragment", "fragments"},
        "sequence_sha256": sha256_text(sequence),
        "source_url": url,
        "retrieved_at": retrieved_at or utc_now(),
        "cache_hit": cache_hit,
        "cache_path": str(cache_path.relative_to(ROOT)) if cache_path and cache_path.is_relative_to(ROOT) else str(cache_path) if cache_path else None,
    }
    if cache_created_at:
        metadata["cache_created_at"] = cache_created_at
    if attempts is not None:
        metadata["request_attempts"] = attempts
    return metadata, sequence


def fetch_uniprot(accession: str, *, cache_dir: Path = DEFAULT_UNIPROT_CACHE, refresh: bool = False) -> tuple[dict, str]:
    encoded = quote(accession, safe="")
    url = f"https://rest.uniprot.org/uniprotkb/{encoded}.json"
    cache_path = cache_dir / f"{safe_id(accession).upper()}.json"
    if cache_path.is_file() and not refresh:
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            payload = cached.get("payload", cached)
            metadata, sequence = _uniprot_metadata(
                payload,
                accession,
                url,
                cache_hit=True,
                cache_path=cache_path,
                retrieved_at=cached.get("retrieved_at"),
                cache_created_at=cached.get("cache_created_at"),
            )
            return metadata, sequence
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
            # A corrupt cache must never silently replace a fresh authoritative response.
            pass

    request = Request(url, headers={"Accept": "application/json", "User-Agent": "openfold3-local-pipeline/1.0"})
    last_error: Exception | None = None
    attempts = 0
    for retry in range(4):
        attempts = retry + 1
        try:
            with urlopen(request, timeout=60) as response:
                payload = json.loads(response.read().decode("utf-8"))
            break
        except HTTPError as exc:
            last_error = exc
            transient = exc.code == 429 or 500 <= exc.code <= 599
            if not transient or retry == 3:
                raise RuntimeError(f"UniProt lookup failed for {accession}: HTTP {exc.code}") from exc
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            try:
                delay = min(8.0, max(0.5, float(retry_after))) if retry_after else min(8.0, 2 ** retry)
            except ValueError:
                delay = min(8.0, 2 ** retry)
            time.sleep(delay)
        except (URLError, TimeoutError) as exc:
            last_error = exc
            if retry == 3:
                raise RuntimeError(f"UniProt lookup failed for {accession}: {exc}") from exc
            time.sleep(min(8.0, 2 ** retry))
    else:
        raise RuntimeError(f"UniProt lookup failed for {accession}: {last_error}") from last_error

    retrieved_at = utc_now()
    cache_created_at = utc_now()
    atomic_write_json(cache_path, {
        "cache_format": 1,
        "accession": accession,
        "source_url": url,
        "retrieved_at": retrieved_at,
        "cache_created_at": cache_created_at,
        "sequence_sha256": sha256_text(validate_sequence(payload.get("sequence", {}).get("value", ""))),
        "payload": payload,
    })
    return _uniprot_metadata(
        payload,
        accession,
        url,
        cache_hit=False,
        cache_path=cache_path,
        retrieved_at=retrieved_at,
        cache_created_at=cache_created_at,
        attempts=attempts,
    )


def run_command(command: list[str], log_path: Path, timeout: int) -> dict:
    started = time.perf_counter()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
        log_path.write_text(completed.stdout, encoding="utf-8")
        return {
            "status": "succeeded" if completed.returncode == 0 else "failed",
            "ok": completed.returncode == 0,
            "exit_code": completed.returncode,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "command": command,
            "log": str(log_path.relative_to(ROOT)),
            "output_tail": completed.stdout[-2000:],
        }
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        log_path.write_text(output, encoding="utf-8")
        return {
            "status": "timed_out",
            "ok": False,
            "exit_code": None,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "command": command,
            "log": str(log_path.relative_to(ROOT)),
            "timeout_seconds": timeout,
            "output_tail": output[-2000:],
        }


@contextmanager
def mps_prediction_lock(path: Path):
    """Serialize local MPS predictions across CLI and API processes."""
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    queued_at = time.perf_counter()
    with path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        acquired_at = time.perf_counter()
        owner = {"pid": os.getpid(), "acquired_at": utc_now()}
        handle.seek(0)
        handle.truncate()
        json.dump(owner, handle)
        handle.flush()
        try:
            yield {
                "lock_path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                "queue_wait_seconds": round(acquired_at - queued_at, 3),
            }
        finally:
            handle.seek(0)
            handle.truncate()
            handle.flush()
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
def normalize_ligands(value: object) -> list[dict]:
    """Validate the bounded OpenFold3 ligand shape used by this entry point."""
    if value is None:
        return []
    if isinstance(value, dict):
        value = value.get("ligands")
    if not isinstance(value, list):
        raise ValueError("ligands JSON must be a list or an object with a 'ligands' list")
    normalized: list[dict] = []
    used_chain_ids = {"A"}
    for index, raw in enumerate(value, start=1):
        if not isinstance(raw, dict):
            raise ValueError(f"ligand {index} must be an object")
        chain_ids = raw.get("chain_ids")
        if chain_ids is None and raw.get("chain_id") is not None:
            chain_ids = [raw["chain_id"]]
        if not isinstance(chain_ids, list) or len(chain_ids) != 1 or not isinstance(chain_ids[0], str):
            raise ValueError(f"ligand {index} must provide exactly one chain_id")
        chain_id = chain_ids[0].strip()
        if not re.fullmatch(r"[A-Za-z0-9]", chain_id):
            raise ValueError(f"ligand {index} chain_id must be one alphanumeric character")
        chain_id = chain_id.upper()
        if chain_id in used_chain_ids:
            raise ValueError(f"duplicate or reserved ligand chain_id: {chain_id}")
        used_chain_ids.add(chain_id)

        representations = []
        smiles = raw.get("smiles")
        if isinstance(smiles, str) and smiles.strip():
            representations.append("smiles")
        ccd_codes = raw.get("ccd_codes")
        if isinstance(ccd_codes, str) and ccd_codes.strip():
            ccd_codes = [ccd_codes.strip()]
        if isinstance(ccd_codes, list) and ccd_codes and all(isinstance(code, str) and code.strip() for code in ccd_codes):
            ccd_codes = [code.strip().upper() for code in ccd_codes]
            representations.append("ccd_codes")
        else:
            ccd_codes = None
        sdf_path = raw.get("sdf_file_path", raw.get("sdf_path"))
        if isinstance(sdf_path, str) and sdf_path.strip():
            sdf_path = str(Path(sdf_path).expanduser().resolve())
            if not Path(sdf_path).is_file():
                raise ValueError(f"ligand {index} SDF does not exist: {sdf_path}")
            representations.append("sdf_file_path")
        if len(representations) != 1:
            raise ValueError(f"ligand {index} must provide exactly one of smiles, ccd_codes, or sdf_file_path")

        ligand_name = str(raw.get("ligand_name") or raw.get("name") or f"ligand{index}").strip()
        if not re.fullmatch(r"[A-Za-z0-9]+", ligand_name):
            raise ValueError(f"ligand {index} ligand_name must contain only ASCII letters and digits")
        ligand = {
            "molecule_type": "ligand",
            "chain_ids": [chain_id],
            "ligand_name": ligand_name,
            "smiles": smiles.strip() if "smiles" in representations else None,
            "ccd_codes": ccd_codes,
            "sdf_file_path": sdf_path if "sdf_file_path" in representations else None,
        }
        normalized.append(ligand)
    return normalized


def load_ligands(value: Path | None) -> list[dict]:
    if value is None:
        return []
    raw_value = str(value)
    try:
        if raw_value.lstrip().startswith(("[", "{")):
            raw = json.loads(raw_value)
        else:
            raw = json.loads(Path(raw_value).expanduser().read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not read --ligands-json: {exc}") from exc
    return normalize_ligands(raw)


def ligand_fingerprint(ligands: list[dict]) -> str:
    canonical = json.dumps(ligands, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def protein_chain(sequence: str, description: str) -> dict:
    return {
        "molecule_type": "PROTEIN",
        "chain_ids": ["A"],
        "description": description,
        "sequence": sequence,
        "non_canonical_residues": None,
        "smiles": None,
        "ligand_name": None,
        "ccd_codes": None,
        "paired_msa_file_paths": None,
        "main_msa_file_paths": None,
        "template_alignment_file_path": None,
        "template_entry_chain_ids": None,
        "template_cif_paths": None,
        "template_cif_chain_ids": None,
        "sdf_file_path": None,
        "cyclic": False,
    }


def write_query(path: Path, query_name: str, sequence: str, description: str, ligands: list[dict] | None = None) -> None:
    chains = [protein_chain(sequence, description)]
    chains.extend(ligands or [])
    payload = {
        "seeds": [42],
        "queries": {
            query_name: {
                "query_name": query_name,
                "chains": chains,
                "use_msas": True,
                "use_paired_msas": False,
                "use_main_msas": True,
                "covalent_bonds": None,
                "pocket_constraint": None,
            }
        },
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_protein_query(path: Path, query_name: str, sequence: str, description: str) -> None:
    write_query(path, query_name, sequence, description, [])


def _rebase_paths(value: object, source_root: Path, destination_root: Path) -> object:
    if isinstance(value, str):
        try:
            path = Path(value)
            if path.is_absolute() and path.is_relative_to(source_root):
                return str(destination_root / path.relative_to(source_root))
        except ValueError:
            pass
        return value
    if isinstance(value, list):
        return [_rebase_paths(item, source_root, destination_root) for item in value]
    if isinstance(value, dict):
        return {key: _rebase_paths(item, source_root, destination_root) for key, item in value.items()}
    return value


def _copy_msa_cache(cache_root: Path, destination: Path) -> tuple[Path, str]:
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(cache_root, destination)
    cached_query_path = destination / "query_msa.json"
    payload = json.loads(cached_query_path.read_text(encoding="utf-8"))
    payload = _rebase_paths(payload, cache_root.resolve(), destination.resolve())
    cached_query_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return cached_query_path, hashlib.sha256(cached_query_path.read_bytes()).hexdigest()


def _msa_cache_valid(cache_root: Path) -> bool:
    query_path = cache_root / "query_msa.json"
    try:
        payload = json.loads(query_path.read_text(encoding="utf-8"))
        if len(payload.get("queries", {})) != 1:
            return False
        query = next(iter(payload["queries"].values()))
        proteins = [chain for chain in query.get("chains", []) if chain.get("molecule_type", "").upper() == "PROTEIN"]
        if len(proteins) != 1:
            return False
        paths = proteins[0].get("main_msa_file_paths") or []
        return bool(paths) and all(Path(path).is_file() for path in paths)
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
        return False


def _merge_cached_msa(cached_query_path: Path, current_query_path: Path, output_path: Path) -> None:
    cached = json.loads(cached_query_path.read_text(encoding="utf-8"))
    current = json.loads(current_query_path.read_text(encoding="utf-8"))
    if len(cached.get("queries", {})) != 1 or len(current.get("queries", {})) != 1:
        raise ValueError("MSA cache rebinding currently supports one query at a time")
    cached_query = next(iter(cached["queries"].values()))
    current_query = next(iter(current["queries"].values()))
    cached_proteins = [chain for chain in cached_query.get("chains", []) if chain.get("molecule_type", "").upper() == "PROTEIN"]
    current_proteins = [chain for chain in current_query.get("chains", []) if chain.get("molecule_type", "").upper() == "PROTEIN"]
    if len(cached_proteins) != 1 or len(current_proteins) != 1 or cached_proteins[0].get("sequence") != current_proteins[0].get("sequence"):
        raise ValueError("MSA cache sequence does not match the current protein query")
    for field in ("paired_msa_file_paths", "main_msa_file_paths", "template_alignment_file_path", "template_entry_chain_ids", "template_cif_paths", "template_cif_chain_ids"):
        current_proteins[0][field] = cached_proteins[0].get(field)
    current_query["use_msas"] = cached_query.get("use_msas", True)
    current_query["use_paired_msas"] = cached_query.get("use_paired_msas", False)
    current_query["use_main_msas"] = cached_query.get("use_main_msas", True)
    output_path.write_text(json.dumps(current, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def acquire_msa(query_path: Path, sequence: str, description: str, msa_dir: Path, *, cache_dir: Path,
                refresh: bool, query_name: str) -> tuple[dict, Path]:
    sequence_hash = sha256_text(sequence)
    params = {
        "format": MSA_CACHE_FORMAT,
        "service": "colabfold-public",
        "use_templates": False,
        "use_paired_msas": False,
        "checkpoint": OPENFOLD_CHECKPOINT,
    }
    params_hash = hashlib.sha256(json.dumps(params, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()[:16]
    cache_root = cache_dir / sequence_hash / params_hash
    cache_query = cache_root / "query_msa.json"
    cache_hit = not refresh and _msa_cache_valid(cache_root)
    request_started = time.perf_counter()
    if not cache_hit:
        cache_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix=f"msa-{sequence_hash[:12]}-", dir=cache_dir.parent))
        try:
            protein_query = temporary / "protein-query.json"
            write_protein_query(protein_query, query_name, sequence, description)
            raw = run_command(
                [str(ROOT / "openfold3/run.sh"), "python", "-u", str(ROOT / "openfold3/fetch_msa.py"), str(protein_query), str(temporary / "cache")],
                temporary / "msa.log",
                900,
            )
            generated = temporary / "cache"
            if not raw["ok"] or not (generated / "query_msa.json").is_file():
                msa_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(temporary / "msa.log", msa_dir / "msa.log")
                failure = {}
                for name in ["failure.json", "network.json"]:
                    if (generated / name).is_file():
                        shutil.copy2(generated / name, msa_dir / name)
                        if name == "failure.json":
                            failure = json.loads((generated / name).read_text())
                return {**raw, **failure, "log": str((msa_dir / "msa.log").relative_to(ROOT)), "cache_hit": False,
                        "cache_key": f"{sequence_hash}/{params_hash}"}, msa_dir / "query_msa.json"
            generated_query = generated / "query_msa.json"
            generated_payload = json.loads(generated_query.read_text(encoding="utf-8"))
            generated_payload = _rebase_paths(generated_payload, generated.resolve(), cache_root.resolve())
            generated_query.write_text(json.dumps(generated_payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            cache_root.parent.mkdir(parents=True, exist_ok=True)
            if cache_root.exists():
                shutil.rmtree(cache_root)
            shutil.move(str(generated), str(cache_root))
        finally:
            shutil.rmtree(temporary, ignore_errors=True)
    query_msa, query_hash = _copy_msa_cache(cache_root, msa_dir)
    _merge_cached_msa(query_msa, query_path, query_msa)
    query_hash = hashlib.sha256(query_msa.read_bytes()).hexdigest()
    elapsed = round(time.perf_counter() - request_started, 3)
    return {
        "status": "succeeded",
        "ok": True,
        "cache_hit": cache_hit,
        "cache_key": f"{sequence_hash}/{params_hash}",
        "cache_path": str(cache_root.relative_to(ROOT)),
        "query_file": str(query_msa.relative_to(ROOT)),
        "query_msa_sha256": query_hash,
        "elapsed_seconds": elapsed,
        "source": "ColabFold public service",
    }, query_msa


def find_prediction_artifacts(prediction_dir: Path, query_name: str) -> tuple[Path | None, Path | None, Path | None]:
    prefix = f"{query_name}/seed_*"
    model = next(iter(sorted(prediction_dir.glob(f"{prefix}/*_model.cif"))), None)
    aggregate = next(iter(sorted(prediction_dir.glob(f"{prefix}/*_confidences_aggregated.json"))), None)
    detail = next(iter(sorted(prediction_dir.glob(f"{prefix}/*_confidences.json"))), None)
    return model, aggregate, detail


def _percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, int(round((len(ordered) - 1) * fraction))))
    return round(ordered[index], 3)


def assess_confidence(confidence: dict, detail: dict, ligand_count: int) -> dict:
    """Turn raw OpenFold3 confidence into an explicit, conservative QC gate."""
    plddt = [float(value) for value in detail.get("plddt", []) if isinstance(value, (int, float))]
    pae = [float(value) for row in detail.get("pae", []) if isinstance(row, list) for value in row if isinstance(value, (int, float))]
    avg_plddt = confidence.get("avg_plddt")
    ptm = confidence.get("ptm")
    iptm = confidence.get("iptm")
    plddt_p10 = _percentile(plddt, 0.10)
    pae_p95 = _percentile(pae, 0.95)
    pae_mean = round(sum(pae) / len(pae), 3) if pae else None
    low_fraction = round(sum(value < 70 for value in plddt) / len(plddt), 4) if plddt else None
    reasons: list[str] = []
    if not plddt:
        reasons.append("per-atom pLDDT is missing")
    if not pae:
        reasons.append("PAE matrix is missing")
    if not isinstance(avg_plddt, (int, float)) or avg_plddt < 70:
        reasons.append("average pLDDT is below 70 or missing")
    if not isinstance(ptm, (int, float)) or ptm < 0.5:
        reasons.append("pTM is below 0.5 or missing")
    if plddt_p10 is not None and plddt_p10 < 50:
        reasons.append("the lowest-confidence 10% of model atoms have pLDDT below 50")
    review_allowed = not reasons
    docking_review_allowed = review_allowed
    if ligand_count and (not isinstance(iptm, (int, float)) or iptm < 0.5):
        docking_review_allowed = False
        reasons.append("protein-ligand interface ipTM is below 0.5 or missing")
    return {
        "status": "pass" if review_allowed else "blocked",
        "structure_review_allowed": review_allowed,
        "docking_structure_review_allowed": docking_review_allowed,
        "decision": "structure_review_allowed" if review_allowed else "low_confidence_structure",
        "metrics": {
            "avg_plddt": avg_plddt,
            "ptm": ptm,
            "iptm": iptm,
            "plddt_count": len(plddt),
            "plddt_p10": plddt_p10,
            "plddt_below_70_fraction": low_fraction,
            "pae_count": len(pae),
            "pae_mean": pae_mean,
            "pae_p95": pae_p95,
        },
        "thresholds": {
            "avg_plddt_min": 70,
            "ptm_min": 0.5,
            "plddt_p10_min": 50,
            "iptm_min_for_ligand_review": 0.5,
        },
        "reasons": reasons,
        "scientific_limit": "This gate permits structural review only; it does not establish binding, affinity, catalytic activity, substrate specificity, or industrial performance.",
    }


def load_benchmark_metrics(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    keys = ("cli_completed", "openfold3_version", "torch_version", "mps_available", "elapsed_seconds", "process_peak_rss_bytes", "sampled_peak_mps_driver_bytes", "system_swap_before_bytes", "system_swap_after_bytes", "memory_note")
    return {key: payload[key] for key in keys if key in payload}


def assess_resources(metrics: dict, lock_info: dict, max_swap_growth_gib: float) -> dict:
    before = metrics.get("system_swap_before_bytes")
    after = metrics.get("system_swap_after_bytes")
    delta_bytes = after - before if isinstance(before, (int, float)) and isinstance(after, (int, float)) else None
    delta_gib = round(delta_bytes / 1024**3, 3) if delta_bytes is not None else None
    warning = delta_gib is not None and delta_gib > max_swap_growth_gib
    return {
        "status": "warning" if warning else "pass",
        "queue_wait_seconds": lock_info.get("queue_wait_seconds"),
        "lock_path": lock_info.get("lock_path"),
        "system_swap_delta_bytes": delta_bytes,
        "system_swap_delta_gib": delta_gib,
        "max_swap_growth_gib": max_swap_growth_gib,
        "scientific_limit": "This is a host resource warning, not a model-quality assessment. A warning should reduce concurrency or move execution to a larger worker; it does not invalidate the generated structure.",
    }


def run(args: argparse.Namespace) -> dict:
    source_kind = "uniprot" if args.accession else "fasta"
    if args.accession and args.fasta:
        raise ValueError("use either --accession or --fasta, not both")
    if not args.accession and not args.fasta:
        raise ValueError("one of --accession or --fasta is required")
    ligands = load_ligands(args.ligands_json)
    ligand_hash = ligand_fingerprint(ligands)

    if args.accession:
        accession = args.accession.strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{1,63}", accession):
            raise ValueError("invalid UniProt accession")
        metadata, sequence = fetch_uniprot(accession, cache_dir=args.uniprot_cache_dir, refresh=args.refresh_uniprot)
        input_id = metadata["primary_accession"]
        description = metadata.get("protein_name") or input_id
    else:
        input_id, raw_sequence = parse_fasta(args.fasta.resolve())
        sequence = validate_sequence(raw_sequence)
        metadata = {
            "accession_requested": None,
            "primary_accession": None,
            "protein_name": input_id,
            "organism": None,
            "length": len(sequence),
            "sequence_sha256": sha256_text(sequence),
            "source_path": str(args.fasta.resolve()),
            "retrieved_at": utc_now(),
            "cache_hit": False,
        }
        description = input_id

    query_name = safe_id(input_id)
    if args.output_dir:
        output_root = args.output_dir.resolve()
    else:
        ligand_suffix = f"ligands-{ligand_hash[:12]}" if ligands else "protein-only"
        output_root = DEFAULT_OUT / query_name / ligand_suffix
    output_root.mkdir(parents=True, exist_ok=True)
    pipeline_path = output_root / "pipeline.json"
    if pipeline_path.is_file() and not args.force and not args.refresh_msa:
        try:
            existing = json.loads(pipeline_path.read_text(encoding="utf-8"))
            if existing.get("status") == "succeeded":
                return {"status": "resumed_existing", "output_dir": str(output_root), "pipeline": str(pipeline_path), "result": existing}
        except json.JSONDecodeError:
            pass
    if pipeline_path.is_file():
        try:
            existing = json.loads(pipeline_path.read_text(encoding="utf-8"))
            existing_hash = existing.get("input", {}).get("ligand_fingerprint")
            if existing_hash and existing_hash != ligand_hash:
                raise ValueError("output directory already contains a different ligand set; choose a new --output-dir")
        except json.JSONDecodeError:
            pass

    total_started = time.perf_counter()
    input_path = output_root / "input.json"
    input_path.write_text(json.dumps({"source_kind": source_kind, "metadata": metadata, "sequence": sequence, "ligands": ligands, "ligand_fingerprint": ligand_hash}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    query_path = output_root / "query.json"
    write_query(query_path, query_name, sequence, description, ligands)

    def progress(stage):
        atomic_write_json(output_root / 'progress.json', {'stage': stage, 'updated_at': datetime.now(timezone.utc).isoformat()})
    progress('msa')
    msa_dir = output_root / "msa"
    msa, query_msa = acquire_msa(
        query_path,
        sequence,
        description,
        msa_dir,
        cache_dir=args.msa_cache_dir.resolve(),
        refresh=args.refresh_msa,
        query_name=query_name,
    )
    if not msa["ok"]:
        result = {"status": "failed", "stage": "msa", "input": {**metadata, "ligands": ligands, "ligand_fingerprint": ligand_hash}, "msa": msa, "output_dir": str(output_root)}
        pipeline_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return result

    if not query_msa.is_file():
        result = {
            "status": "failed",
            "stage": "msa", "input": {**metadata, "ligands": ligands, "ligand_fingerprint": ligand_hash},
            "msa": msa,
            "reason": "MSA command completed without query_msa.json",
            "output_dir": str(output_root),
        }
        pipeline_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return result
    msa["query_file"] = str(query_msa.relative_to(ROOT))
    msa["query_msa_sha256"] = hashlib.sha256(query_msa.read_bytes()).hexdigest()
    prediction_dir = output_root / "prediction"
    progress('prediction_queue')
    with mps_prediction_lock(args.mps_lock_path) as lock_info:
        progress('prediction')
        prediction = run_command(
            [
                str(ROOT / "openfold3/run.sh"), "python", "-u", str(ROOT / "openfold3/benchmark.py"),
                str(output_root / "benchmark.json"), "predict", "--query-json", str(query_msa),
                "--runner-yaml", str(ROOT / "openfold3/runner.yml"), "--num-diffusion-samples", str(args.samples),
                "--use-msa-server", "false", "--use-templates", "false", "--output-dir", str(prediction_dir),
            ],
            output_root / "prediction.log",
            args.timeout,
        )
    prediction["queue_wait_seconds"] = lock_info["queue_wait_seconds"]
    prediction["mps_lock_path"] = lock_info["lock_path"]
    model_path, aggregate_path, detail_path = find_prediction_artifacts(prediction_dir, query_name)
    confidence = json.loads(aggregate_path.read_text(encoding="utf-8")) if aggregate_path else None
    confidence_detail = json.loads(detail_path.read_text(encoding="utf-8")) if detail_path else {}
    if not prediction["ok"] or confidence is None or model_path is None:
        result = {
            "status": "failed", "stage": "prediction", "input": {**metadata, "ligands": ligands, "ligand_fingerprint": ligand_hash}, "msa": msa, "prediction": prediction,
            "output_dir": str(output_root), "artifacts": {"model": str(model_path) if model_path else None},
        }
        pipeline_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return result
    quality_gate = assess_confidence(confidence, confidence_detail, len(ligands))
    benchmark_metrics = load_benchmark_metrics(output_root / "benchmark.json")
    resource_guard = assess_resources(benchmark_metrics, lock_info, args.max_swap_growth_gib)
    quality_path = output_root / "quality-gate.json"
    atomic_write_json(quality_path, quality_gate)

    progress('viewer_export')
    viewer_dir = output_root / "viewer"
    viewer = {"status": "skipped", "ok": True}
    if not args.no_viewer:
        seed_prediction_dir = next(iter(sorted(prediction_dir.glob(f"{query_name}/seed_*"))), None)
        if seed_prediction_dir is None:
            viewer = {"status": "failed", "ok": False, "reason": "prediction seed directory is missing"}
        else:
            viewer = run_command(
                [
                    sys.executable, str(ROOT / "openfold3/export_viewer_artifact.py"),
                    "--prediction-dir", str(seed_prediction_dir),
                    "--output-dir", str(viewer_dir),
                    "--quality-json", str(quality_path),
                    "--artifact-id", f"{query_name}-openfold",
                    "--title", description,
                    "--model-version", "openbind-2025-06-30-174k",
                    "--seed", "42",
                    "--no-templates",
                ],
                output_root / "viewer-export.log",
                300,
            )
    if not viewer["ok"]:
        result = {
            "status": "failed", "stage": "viewer_export", "input": {**metadata, "ligands": ligands, "ligand_fingerprint": ligand_hash}, "msa": msa,
            "prediction": prediction, "confidence": confidence, "viewer": viewer,
            "output_dir": str(output_root),
        }
        pipeline_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return result

    avg_plddt = confidence.get("avg_plddt")
    ptm = confidence.get("ptm")
    iptm = confidence.get("iptm")
    decision = quality_gate["decision"]
    result = {
        "status": "succeeded",
        "stage": "complete",
        "input": {**metadata, "ligands": ligands, "ligand_fingerprint": ligand_hash},
        "sequence": {"length": len(sequence), "sha256": sha256_text(sequence)},
        "msa": msa,
        "prediction": prediction,
        "confidence": {"avg_plddt": avg_plddt, "ptm": ptm, "iptm": iptm},
        "quality_gate": quality_gate,
        "viewer": viewer,
        "decision": decision,
        "artifacts": {
            "input": str(input_path.relative_to(ROOT)),
            "query": str(query_path.relative_to(ROOT)),
            "query_msa": str(query_msa.relative_to(ROOT)),
            "model_cif": str(model_path.relative_to(ROOT)),
            "confidence_aggregated": str(aggregate_path.relative_to(ROOT)),
            "confidence_detail": str(detail_path.relative_to(ROOT)) if detail_path else None,
            "viewer_manifest": str((viewer_dir / "manifest.json").relative_to(ROOT)) if not args.no_viewer else None,
            "quality_gate": str(quality_path.relative_to(ROOT)),
            "benchmark": str((output_root / "benchmark.json").relative_to(ROOT)),
        },
        "provenance": {
            "model": "OpenFold3",
            "checkpoint": OPENFOLD_CHECKPOINT,
            "device": "Apple MPS",
            "templates": False,
            "seed": 42,
            "diffusion_samples": args.samples,
            "msa_source": "ColabFold public service, cached locally",
            "uniprot_cache": metadata.get("cache_path") if source_kind == "uniprot" else None,
            "msa_cache": msa.get("cache_path"),
            "ligand_fingerprint": ligand_hash,
        },
        "resources": benchmark_metrics,
        "resource_guard": resource_guard,
        "elapsed_seconds": round(time.perf_counter() - total_started, 3),
        "scientific_limit": "Structure confidence supports structural review; it does not establish biochemical activity, substrate specificity, affinity, or industrial performance.",
    }
    pipeline_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {"status": "succeeded", "output_dir": str(output_root), "pipeline": str(pipeline_path), "result": result}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--accession", help="explicit UniProt accession")
    source.add_argument("--fasta", type=Path, help="one-record protein FASTA")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--ligands-json", help="JSON file path or inline JSON list of ligands")
    parser.add_argument("--uniprot-cache-dir", type=Path, default=DEFAULT_UNIPROT_CACHE)
    parser.add_argument("--msa-cache-dir", type=Path, default=DEFAULT_MSA_CACHE)
    parser.add_argument("--mps-lock-path", type=Path, default=DEFAULT_MPS_LOCK)
    parser.add_argument("--max-swap-growth-gib", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--refresh-uniprot", action="store_true", help="bypass and replace the UniProt response cache")
    parser.add_argument("--refresh-msa", action="store_true", help="bypass and replace the content-addressed MSA cache")
    parser.add_argument("--no-viewer", action="store_true", help="skip generation of the frontend structure manifest")
    args = parser.parse_args()
    if args.samples < 1 or args.samples > 8:
        parser.error("--samples must be between 1 and 8")
    args.uniprot_cache_dir = args.uniprot_cache_dir.resolve()
    args.msa_cache_dir = args.msa_cache_dir.resolve()
    args.mps_lock_path = args.mps_lock_path.resolve()
    if args.max_swap_growth_gib < 0:
        parser.error("--max-swap-growth-gib must be non-negative")
    try:
        print(json.dumps(run(args), ensure_ascii=False, indent=2), flush=True)
    except Exception as exc:
        print(json.dumps({"status": "failed", "error": f"{type(exc).__name__}: {exc}"}, ensure_ascii=False, indent=2), flush=True)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
