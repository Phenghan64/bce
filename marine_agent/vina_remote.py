"""Transport for running validated Vina jobs on the Mac worker."""
from __future__ import annotations

import base64
import json
import os
import tarfile
from io import BytesIO
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
VINA_WORKER_URL = os.environ.get("VINA_WORKER_URL", "").rstrip("/")
VINA_WORKER_TOKEN = os.environ.get("VINA_WORKER_TOKEN", "")
VINA_WORKER_TIMEOUT = max(30, int(os.environ.get("VINA_WORKER_TIMEOUT_S", "7200")))


def _safe_extract_bundle(encoded: str, expected_root: str) -> None:
    try:
        archive_bytes = base64.b64decode(encoded, validate=True)
    except (ValueError, TypeError) as exc:
        raise RuntimeError("Vina worker returned an invalid artifact bundle") from exc
    if len(archive_bytes) > 512 * 1024 * 1024:
        raise RuntimeError("Vina worker artifact bundle exceeds the 512 MiB limit")
    destination = (ROOT / expected_root).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=BytesIO(archive_bytes), mode="r:gz") as archive:
        members = archive.getmembers()
        if not members:
            raise RuntimeError("Vina worker returned an empty artifact bundle")
        for member in members:
            member_path = Path(member.name)
            if member_path.is_absolute() or ".." in member_path.parts or not (member.name == expected_root or member.name.startswith(expected_root + '/')):
                raise RuntimeError("Vina worker returned an unsafe artifact path")
            if member.issym() or member.islnk() or not (member.isfile() or member.isdir()):
                raise RuntimeError("Vina worker returned an unsupported artifact entry")
            try:
                (ROOT / member.name).resolve().relative_to(ROOT.resolve())
            except ValueError as exc:
                raise RuntimeError("Vina worker artifact escaped the result directory") from exc
        for member in members:
            archive.extract(member, ROOT)
    if not destination.exists():
        raise RuntimeError("Vina worker bundle did not contain the declared output directory")


def run_vina_worker(run_id: str, config: dict) -> dict:
    if not VINA_WORKER_URL:
        return {"status": "unavailable", "stage": "worker", "reason": "VINA_WORKER_URL is not configured"}
    output_root = f"results/vina-runs/{run_id}"
    request_payload = {
        "receptor": {"name": config["receptor"]["name"], "content_base64": base64.b64encode(config["receptor"]["bytes"]).decode("ascii")},
        "ligands": [
            {"id": item["id"], "name": item["name"], "content_base64": base64.b64encode(item["bytes"]).decode("ascii")}
            for item in config["ligands"]
        ],
        "box": config["box"],
        "params": config["params"],
        "task_seeds": config.get("task_seeds", {}),
        "receptor_preparation": config.get('provenance',{}).get('receptor_preparation'),
        "ligand_library": config.get("provenance",{}).get("ligand_library"),
    }
    body = json.dumps({"job_id": run_id, "output_root": output_root, "request": request_payload}, ensure_ascii=False).encode("utf-8")
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if VINA_WORKER_TOKEN:
        headers["X-Vina-Worker-Token"] = VINA_WORKER_TOKEN
    request = Request(f"{VINA_WORKER_URL}/v1/vina", data=body, headers=headers, method="POST")
    try:
        with urlopen(request, timeout=VINA_WORKER_TIMEOUT) as response:
            response_payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")[-2000:]
        return {"status": "failed", "stage": "worker", "reason": f"worker HTTP {exc.code}: {detail}"}
    except (URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
        return {"status": "unavailable", "stage": "worker", "reason": str(exc)}
    if not isinstance(response_payload, dict):
        return {"status": "failed", "stage": "worker", "reason": "worker returned a non-object response"}
    result = response_payload.get("result")
    if not isinstance(result, dict):
        result = {"status": response_payload.get("status", "failed"), "reason": response_payload.get("reason", "worker result is missing")}
    bundle = response_payload.get("bundle")
    if isinstance(bundle, dict) and isinstance(bundle.get("data"), str):
        declared_root = str(bundle.get("root") or output_root)
        if declared_root != output_root:
            return {"status": "failed", "stage": "worker", "reason": "worker output root does not match submitted job"}
        try:
            _safe_extract_bundle(bundle["data"], declared_root)
        except RuntimeError as exc:
            return {"status": "failed", "stage": "worker_artifacts", "reason": str(exc)}
    result.setdefault("execution_node", {"worker_url": VINA_WORKER_URL, "run_id": run_id})
    return result


def call_receptor_worker(action: str, payload: dict) -> dict:
    if action not in {'inspect', 'prepare'} or not VINA_WORKER_URL:
        raise RuntimeError('本机受体准备服务未配置')
    headers = {'Content-Type': 'application/json', 'X-Vina-Worker-Token': VINA_WORKER_TOKEN}
    request = Request(f'{VINA_WORKER_URL}/v1/receptor/{action}', data=json.dumps(payload).encode(), headers=headers, method='POST')
    try:
        with urlopen(request, timeout=420 if action == 'prepare' else 30) as response:
            result = json.loads(response.read())
    except HTTPError as exc:
        try:
            detail = json.loads(exc.read()).get('reason', '本机准备服务请求失败')
        except (ValueError, TypeError, AttributeError):
            detail = '本机准备服务请求失败'
        raise ValueError(detail) from exc
    if result.get('status') != 'succeeded':
        raise ValueError(result.get('reason', '本机受体准备失败'))
    if action == 'prepare':
        expected_root = f"results/vina-receptors/{payload['job_id']}"
        bundle = result.get('bundle') or {}
        if bundle.get('root') != expected_root:
            raise RuntimeError('受体准备产物路径不匹配')
        _safe_extract_bundle(bundle['data'], expected_root)
    return result.get('inspection') if action == 'inspect' else result['result']


def recover_vina_artifacts(run_id: str) -> dict | None:
    """Read durable per-ligand results only after the old worker is quiescent."""
    import re
    if not VINA_WORKER_URL:
        return None
    if not re.fullmatch(r'[0-9a-f-]{32,36}', run_id):
        raise ValueError('invalid run ID')
    request = Request(f'{VINA_WORKER_URL}/v1/jobs/{run_id}/artifacts',
                      headers={'X-Vina-Worker-Token': VINA_WORKER_TOKEN})
    from urllib.request import build_opener, ProxyHandler
    try:
        with build_opener(ProxyHandler({})).open(request, timeout=60) as response:
            payload = json.load(response)
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise RuntimeError('尚未确认旧任务产物，暂不重复计算') from exc
    root = f'results/vina-runs/{run_id}'
    bundle = payload.get('bundle') or {}
    if bundle.get('root') != root or not isinstance(payload.get('result'), dict):
        raise RuntimeError('recovered Vina artifact identity mismatch')
    _safe_extract_bundle(bundle.get('data'), root)
    return payload['result']
