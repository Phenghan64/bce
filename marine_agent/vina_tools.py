"""Bounded, auditable AutoDock Vina execution for the docking workbench.

The workbench accepts prepared PDBQT files and keeps every input and parameter
inside a durable run directory.  Preparation of PDB/SDF files is intentionally
not guessed here: callers must supply the prepared PDBQT and can inspect the
exact preparation provenance before launching Vina.
"""
from __future__ import annotations

from .execution_control import RunCancelled, current_control

import base64
import hashlib
import itertools
import json
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
VINA_BINARY = os.environ.get("VINA_BINARY", "vina")
MAX_RECEPTOR_BYTES = 25 * 1024 * 1024
MAX_LIGAND_BYTES = 8 * 1024 * 1024
MAX_TOTAL_BYTES = 50 * 1024 * 1024


class VinaInputError(ValueError):
    """A user-supplied docking input failed a deterministic validation gate."""


def _safe_name(value: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._")
    return (cleaned[:120] or fallback)


def _decode_content(item: dict[str, Any], label: str, max_bytes: int) -> tuple[str, bytes]:
    name = _safe_name(str(item.get("name") or label), label)
    encoded = item.get("content_base64")
    if not isinstance(encoded, str) or not encoded:
        raise VinaInputError(f"{label} 缺少 PDBQT 文件内容")
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (ValueError, base64.binascii.Error) as exc:
        raise VinaInputError(f"{label} 不是有效的 base64 文件") from exc
    if not raw or len(raw) > max_bytes:
        raise VinaInputError(f"{label} 文件大小超出限制")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise VinaInputError(f"{label} 必须是 UTF-8 文本 PDBQT") from exc
    if not re.search(r"^(ATOM  |HETATM)", text, flags=re.MULTILINE):
        raise VinaInputError(f"{label} 未检测到 ATOM/HETATM 记录")
    if "ROOT" not in text and label == "配体":
        raise VinaInputError("配体 PDBQT 缺少 ROOT 记录，无法进行 Vina 对接")
    return name, raw


def _number(value: Any, label: str, low: float, high: float) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise VinaInputError(f"{label} 必须是数字") from exc
    if not low <= result <= high:
        raise VinaInputError(f"{label} 必须在 {low:g} 到 {high:g} 之间")
    return result


def validate_config(request: dict[str, Any]) -> dict[str, Any]:
    receptor = request.get("receptor")
    ligands = request.get("ligands")
    if not isinstance(receptor, dict):
        raise VinaInputError("缺少受体 PDBQT")
    if not isinstance(ligands, list) or not ligands or len(ligands) > 100:
        raise VinaInputError("配体数量必须在 1 到 100 个之间")
    total = 0
    receptor_name, receptor_bytes = _decode_content(receptor, "受体", MAX_RECEPTOR_BYTES)
    total += len(receptor_bytes)
    ligand_records = []
    seen: set[str] = set()
    for index, ligand in enumerate(ligands, 1):
        if not isinstance(ligand, dict):
            raise VinaInputError(f"第 {index} 个配体格式无效")
        name, raw = _decode_content(ligand, "配体", MAX_LIGAND_BYTES)
        ligand_id = _safe_name(str(ligand.get("id") or Path(name).stem), f"ligand-{index}")
        if ligand_id in seen:
            raise VinaInputError(f"配体 ID 重复：{ligand_id}")
        seen.add(ligand_id)
        total += len(raw)
        ligand_records.append({"id": ligand_id, "name": name, "bytes": raw})
    if total > MAX_TOTAL_BYTES:
        raise VinaInputError("本次上传文件总大小超过 50 MB")

    box = request.get("box") or {}
    center = [_number(value, f"搜索盒中心 {axis}", -1000, 1000) for axis, value in zip(("x", "y", "z"), box.get("center", []))]
    size = [_number(value, f"搜索盒尺寸 {axis}", 1, 200) for axis, value in zip(("x", "y", "z"), box.get("size", []))]
    if len(center) != 3 or len(size) != 3:
        raise VinaInputError("搜索盒必须提供 3 个中心坐标和 3 个尺寸")

    params = request.get("params") or {}
    def integer(name: str, default: int, low: int, high: int) -> int:
        try:
            value = int(params.get(name, default))
        except (TypeError, ValueError) as exc:
            raise VinaInputError(f"{name} 必须是整数") from exc
        if not low <= value <= high:
            raise VinaInputError(f"{name} 必须在 {low} 到 {high} 之间")
        return value

    primary_seed = integer("seed", 20261028, 0, 2_147_483_647)
    raw_seeds = params.get("seeds", [primary_seed])
    if not isinstance(raw_seeds, list) or not 1 <= len(raw_seeds) <= 5:
        raise VinaInputError("seeds 必须包含 1 到 5 个随机种子")
    seeds = []
    for value in raw_seeds:
        if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= 2_147_483_647:
            raise VinaInputError("每个 seed 都必须是规定范围内的整数")
        if value in seeds:
            raise VinaInputError("seeds 不能包含重复值")
        seeds.append(value)
    requested_task_seeds = request.get("task_seeds") or {}
    if not isinstance(requested_task_seeds, dict):
        raise VinaInputError("task_seeds 格式无效")
    task_seeds = {}
    for ligand in ligand_records:
        selected = requested_task_seeds.get(ligand["id"], seeds)
        if not isinstance(selected, list) or not selected or any(seed not in seeds for seed in selected):
            raise VinaInputError(f"配体 {ligand['id']} 的续跑种子不属于本次 seeds")
        task_seeds[ligand["id"]] = list(dict.fromkeys(selected))
    if sum(len(values) for values in task_seeds.values()) > 100:
        raise VinaInputError("配体与种子展开后的对接任务总数不能超过 100")

    return {
        "engine": "vina",
        "receptor": {"name": receptor_name, "bytes": receptor_bytes},
        "ligands": ligand_records,
        "box": {"center": center, "size": size},
        "task_seeds": task_seeds,
        "params": {
            "exhaustiveness": integer("exhaustiveness", 8, 1, 64),
            "num_modes": integer("num_modes", 9, 1, 20),
            "energy_range": _number(params.get("energy_range", 4), "energy_range", 0, 20),
            "seed": seeds[0],
            "seeds": seeds,
            "cpu": integer("cpu", 1, 1, 32),
            "timeout_seconds": integer("timeout_seconds", 900, 30, 7200),
        },
        "provenance": {
            "receptor_preparation": request.get('receptor_preparation'),
            "ligand_library": request.get("ligand_library"),
            "input_format": "PDBQT",
            "preparation_gate": "用户确认的受体和配体已准备为 PDBQT；系统不自动猜测质子化、互变异构体或手性",
            "scientific_limit": "Vina 亲和力估计值仅用于构象和候选排序，不等同于实验结合常数或酶活",
        },
    }


def _parse_scores(text: str) -> list[dict[str, Any]]:
    rows = []
    pattern = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)\s+([\d.]+)\s+([\d.]+)", re.MULTILINE)
    for mode, (affinity, rmsd_lb, rmsd_ub) in enumerate(pattern.findall(text), 1):
        rows.append({"mode": mode, "affinity_kcal_mol": float(affinity), "rmsd_lb": float(rmsd_lb), "rmsd_ub": float(rmsd_ub)})
    return rows


def _seed_tasks(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand ligands into independently reproducible seed tasks."""
    params = config.get("params", {})
    seeds = params.get("seeds") or [params.get("seed", 20261028)]
    tasks = []
    for ligand in config["ligands"]:
        selected = (config.get("task_seeds") or {}).get(ligand["id"], seeds)
        for seed in selected:
            task_id = ligand["id"] if len(seeds) == 1 else f"{ligand['id']}-seed-{seed}"
            tasks.append({"ligand": ligand, "seed": seed, "id": task_id, "index": len(tasks) + 1})
    return tasks


def _analyze_seed_consistency(run_dir: Path, config: dict[str, Any], results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compare top poses across seeds without ligand-only fitting.

    Meeko restores chemistry before RDKit CalcRMS enumerates symmetry matches.
    Analysis failures are explicit and never turn a failed docking into success.
    """
    groups = []
    ligand_ids = list(dict.fromkeys(row.get("ligand_id") or row["id"] for row in results))
    for ligand_id in ligand_ids:
        rows = [row for row in results if (row.get("ligand_id") or row["id"]) == ligand_id]
        successful = [row for row in rows if row.get("status") == "succeeded" and row.get("scores")]
        if len(successful) < 2:
            continue
        group = {"ligand_id": ligand_id, "name": successful[0]["name"], "seed_count": len(successful),
                 "seeds": [row["seed"] for row in successful]}
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolAlign
            exporter = Path(VINA_BINARY).resolve().parents[2] / ".env/bin/mk_export.py"
            if not exporter.is_file():
                raise RuntimeError(f"Meeko exporter unavailable: {exporter}")
            molecules = []
            for row in successful:
                source = ROOT / row["output"]
                sdf = run_dir / "ligands" / f"{row['id']}.sdf"
                completed = (control.run if (control := current_control()) else subprocess.run)(
                    [str(exporter), str(source), "-s", str(sdf)],
                    capture_output=True, text=True, timeout=60)
                if completed.returncode != 0 or not sdf.is_file():
                    raise RuntimeError(f"Meeko export failed for seed {row['seed']}: {(completed.stdout or completed.stderr or '')[-500:]}")
                molecules_found = list(Chem.SDMolSupplier(str(sdf), removeHs=True))
                if not molecules_found or molecules_found[0] is None:
                    raise RuntimeError(f"Meeko export has no readable top pose for seed {row['seed']}")
                row["sdf_output"] = str(sdf.relative_to(ROOT))
                molecules.append(molecules_found[0])
            pairs = []
            for i, j in itertools.combinations(range(len(successful)), 2):
                pairs.append({"seed_a": successful[i]["seed"], "seed_b": successful[j]["seed"],
                              "top_pose_rmsd_angstrom": float(rdMolAlign.CalcRMS(molecules[i], molecules[j]))})
            scores = [row["scores"][0]["affinity_kcal_mol"] for row in successful]
            medoid_index = min(range(len(successful)), key=lambda i: sum(
                pair["top_pose_rmsd_angstrom"] for pair in pairs
                if pair["seed_a"] == successful[i]["seed"] or pair["seed_b"] == successful[i]["seed"]))
            max_rmsd = max(pair["top_pose_rmsd_angstrom"] for pair in pairs)
            group.update({
                "status": "available",
                "method": "Meeko SDF chemistry restoration; symmetry-aware heavy-atom RDKit CalcRMS; fixed receptor frame; no ligand-only fitting",
                "pairs": pairs,
                "max_top_pose_rmsd_angstrom": max_rmsd,
                "top_score_range_kcal_mol": max(scores) - min(scores),
                "central_top_pose_seed": successful[medoid_index]["seed"],
                "recommended_review_task_id": successful[medoid_index]["id"],
                "recommended_review_seed": successful[medoid_index]["seed"],
                "review_status": "consistent_review_required" if max_rmsd <= 2.0 else "divergent_review_required",
                "review_action": "默认查看几何中心首位姿势并复核受体、口袋和接触；这不是按最低分自动选定正确姿势。",
                "review_threshold_angstrom": 2.0,
                "top_poses_within_2A": max_rmsd <= 2.0,
                "scientific_limit": "跨种子一致性只描述计算姿势是否接近；2 Å 是复核阈值，不证明结合正确性、亲和力或酶活。分数接近也不代表姿势一致。",
            })
        except Exception as exc:
            group.update({"status": "unavailable", "review_status": "analysis_unavailable",
                          "review_action": "保留各原始种子结果，先解决化学导出或 RMSD 分析问题。",
                          "reason": f"{type(exc).__name__}: {exc}"[:1000]})
        groups.append(group)
    return groups


def run_batch(run_dir: Path, config: dict[str, Any], emit: Callable[[dict[str, Any]], None] | None = None, control=None) -> dict[str, Any]:
    control = control or current_control()
    run_dir = run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = run_dir / "receptor.pdbqt"
    receptor_path.write_bytes(config["receptor"]["bytes"])
    (run_dir / "config.json").write_text(json.dumps({k: v for k, v in config.items() if k not in {"receptor", "ligands"}}, ensure_ascii=False, indent=2) + "\n")
    ligand_dir = run_dir / "ligands"
    ligand_dir.mkdir(exist_ok=True)
    for ligand in config["ligands"]:
        (ligand_dir / f"{ligand['id']}.pdbqt").write_bytes(ligand["bytes"])

    p = config["params"]
    box = config["box"]
    tasks = _seed_tasks(config)
    results: list[dict[str, Any]] = []

    def checkpoint(payload, name='checkpoint.json'):
        target = run_dir / name
        temporary = target.with_suffix(target.suffix + '.tmp')
        with temporary.open('w') as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(target)

    receptor_sha = hashlib.sha256(config['receptor']['bytes']).hexdigest()
    engine_path = Path(VINA_BINARY)
    engine_sha = hashlib.sha256(engine_path.read_bytes()).hexdigest() if engine_path.is_file() else None
    batch_started = time.perf_counter()

    def run_one(task: dict[str, Any]) -> dict[str, Any]:
        ligand, seed, task_id, index = task["ligand"], task["seed"], task["id"], task["index"]
        ligand_id = ligand["id"]
        output = run_dir / "ligands" / f"{task_id}.out.pdbqt"
        log_path = run_dir / "ligands" / f"{task_id}.log"
        command = [VINA_BINARY, "--receptor", str(receptor_path), "--ligand", str(ligand_dir / f"{ligand_id}.pdbqt"),
                   "--center_x", str(box["center"][0]), "--center_y", str(box["center"][1]), "--center_z", str(box["center"][2]),
                   "--size_x", str(box["size"][0]), "--size_y", str(box["size"][1]), "--size_z", str(box["size"][2]),
                   "--exhaustiveness", str(p["exhaustiveness"]), "--num_modes", str(p["num_modes"]),
                   "--energy_range", str(p["energy_range"]), "--seed", str(seed), "--cpu", str(p["cpu"]), "--out", str(output)]
        base = {"id": task_id, "ligand_id": ligand_id, "seed": seed, "name": ligand["name"]}
        if control and control.cancelled:
            return {**base, "status": "skipped", "scores": [], "error": "任务已停止，未开始计算"}
        if emit:
            emit({"stage": "对接计算", "status": "running", "ligand_id": task_id,
                  "source_ligand_id": ligand_id, "seed": seed, "index": index, "total": len(tasks)})
        started = time.perf_counter()
        try:
            completed = (control.run if control else subprocess.run)(command, capture_output=True, text=True, timeout=p["timeout_seconds"])
            combined = (completed.stdout or "") + (completed.stderr or "")
            log_path.write_text(combined)
            output_text = output.read_text(errors="replace") if output.is_file() else ""
            scores = _parse_scores(output_text)
            status = "succeeded" if completed.returncode == 0 and scores and output.is_file() else "failed"
            row = {**base, "status": status, "returncode": completed.returncode,
                   "elapsed_seconds": round(time.perf_counter() - started, 3), "scores": scores,
                   "output": str(output.relative_to(ROOT)), "log": str(log_path.relative_to(ROOT)),
                   "command": command, "input_sha256": hashlib.sha256(ligand["bytes"]).hexdigest()}
            if status != "succeeded":
                row["error"] = combined[-1200:] or "Vina 没有返回有效 pose"
            else:
                row["output_sha256"] = hashlib.sha256(output.read_bytes()).hexdigest()
        except RunCancelled:
            log_path.write_text("用户停止了当前对接任务\n")
            row = {**base, "status": "stopped", "scores": [], "elapsed_seconds": round(time.perf_counter() - started, 3),
                   "error": "用户已停止计算", "log": str(log_path.relative_to(ROOT))}
        except FileNotFoundError:
            raise RuntimeError(f"Vina 引擎不可用：{VINA_BINARY}")
        except subprocess.TimeoutExpired as exc:
            log_path.write_text(str(exc.stdout or "")[-4000:])
            row = {**base, "status": "timed_out", "returncode": None, "elapsed_seconds": round(time.perf_counter() - started, 3),
                   "scores": [], "output": str(output.relative_to(ROOT)), "log": str(log_path.relative_to(ROOT)),
                   "command": command, "error": f"超过 {p['timeout_seconds']} 秒"}
        return row

    requested_workers = max(1, int(config.get("batch_workers", os.environ.get("VINA_BATCH_WORKERS", "1"))))
    cpu_budget = max(1, int(os.environ.get("VINA_CPU_BUDGET", str(os.cpu_count() or 1))))
    worker_count = min(len(tasks), requested_workers, max(1, cpu_budget // max(1, int(p["cpu"]))))
    with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="vina") as pool:
        futures = [pool.submit(run_one, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            results.append(row)
            checkpoint({'status': 'partial', 'engine': 'AutoDock Vina',
                        'configuration': {'box': config['box'], 'params': config['params']},
                        'receptor_sha256': receptor_sha, 'engine_sha256': engine_sha,
                        'provenance': config['provenance'], 'results': list(results),
                        'elapsed_seconds': round(time.perf_counter() - batch_started, 3)})
            if emit:
                succeeded = sum(item['status'] == 'succeeded' for item in results)
                emit({'stage': '对接计算', 'name': '配体种子完成', 'ligand_id': row['id'],
                      'source_ligand_id': row['ligand_id'], 'seed': row['seed'], 'status': row['status'],
                      'completed': sum(item['status'] in {'succeeded', 'failed', 'timed_out'} for item in results),
                      'total': len(tasks), 'succeeded': succeeded,
                      'failed': sum(item['status'] in {'failed', 'timed_out'} for item in results)})
    order = {task['id']: task['index'] for task in tasks}
    results.sort(key=lambda item: order[item['id']])

    payload = {
        "status": "stopped" if control and control.cancelled else "succeeded" if all(item["status"] == "succeeded" for item in results) else "partial",
        "engine": "AutoDock Vina",
        "vina_binary": VINA_BINARY,
        "receptor_sha256": receptor_sha, "engine_sha256": engine_sha,
        "configuration": {"box": config["box"], "params": config["params"]},
        "provenance": config["provenance"],
        "results": results,
        "ligand_count": len(config['ligands']),
        "task_count": len(tasks),
        "top_pose_consistency": _analyze_seed_consistency(run_dir, config, results),
        "elapsed_seconds": round(time.perf_counter() - batch_started, 3),
        "execution": {"batch_workers": worker_count, "cpu_budget": cpu_budget, "per_docking_cpu": p["cpu"]},
        "ranking": sorted(({"id": item["id"], "ligand_id": item["ligand_id"], "seed": item["seed"],
                             "name": item["name"], "status": item["status"],
                             "top_affinity_kcal_mol": item["scores"][0]["affinity_kcal_mol"] if item["scores"] else None}
                            for item in results if item["status"] == "succeeded"),
                           key=lambda item: item["top_affinity_kcal_mol"]),
    }
    checkpoint(payload, "result.json")
    return payload
