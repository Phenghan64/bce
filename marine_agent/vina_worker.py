#!/usr/bin/env python3
"""Tailscale-only AutoDock Vina execution worker.

The control-plane API stays on the Linux node.  This small authenticated
worker runs the validated Vina batch on the Mac that also owns OpenFold3.
Only prepared PDBQT content and bounded configuration cross the worker
boundary; no shell or arbitrary filesystem operation is exposed.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import re
import tarfile
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from pathlib import Path
from threading import Lock

from marine_agent.vina_tools import ROOT, VINA_BINARY, run_batch, validate_config
from marine_agent.execution_control import RunCancelled, control_for, job_scope, job_lock, job_status

JOB_RE = re.compile(r"^(?:[0-9a-f]{32}|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$", re.IGNORECASE)
JOB_LOCK = Lock()


def _json_bytes(payload: dict) -> bytes:
    return json.dumps(payload, ensure_ascii=False).encode("utf-8")


def _bundle(output_root: Path, root_name: str) -> str:
    buffer = BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        archive.add(output_root, arcname=root_name, recursive=True)
    return base64.b64encode(buffer.getvalue()).decode("ascii")


def _decode_request(payload: dict) -> dict:
    """Convert transport base64 fields back to the normal Vina request."""
    if not isinstance(payload, dict):
        raise ValueError("request must be an object")
    job_id = str(payload.get("job_id", ""))
    output_root = str(payload.get("output_root", ""))
    if not JOB_RE.fullmatch(job_id):
        raise ValueError("job_id is invalid")
    if output_root != f"results/vina-runs/{job_id}":
        raise ValueError("output_root is invalid")
    request = payload.get("request")
    if not isinstance(request, dict):
        raise ValueError("request payload is missing")
    return request


def _run_job(payload: dict, timeout: int) -> dict:
    job_id = str(payload.get("job_id", ""))
    root_name = str(payload.get("output_root", ""))
    try:
        request = _decode_request(payload)
        config = validate_config(request)
        output_root = (ROOT / root_name).resolve()
        output_root.relative_to((ROOT / "results" / "vina-runs").resolve())

        def emit(event: dict) -> None:
            # The control plane receives the final report.  Progress remains
            # available in the worker log and result provenance; the API job
            # itself stays in a running state while this bounded call runs.
            return None

        # A worker-level cap lets the Mac use several independent ligands at
        # once without changing the public Vina protocol.  cpu remains the
        # per-docking thread count and the cap is applied by run_batch.
        # Keep scheduling internal: users only choose the official Vina
        # per-docking CPU value. The worker fills the Mac's CPU budget with
        # independent ligand jobs and automatically reduces concurrency when
        # a docking uses more threads.
        config["batch_workers"] = max(1, int(os.environ.get("VINA_BATCH_WORKERS", "8")))
        callback_url = os.environ.get("VINA_CALLBACK_URL", "").rstrip("/")
        callback_token = os.environ.get("VINA_CALLBACK_TOKEN", "")

        def report_progress(event: dict) -> None:
            emit(event)
            if not callback_url:
                return
            body = json.dumps(event, ensure_ascii=False).encode("utf-8")
            request = urllib.request.Request(
                f"{callback_url}/internal/vina-progress/{job_id}",
                data=body,
                headers={
                    "Content-Type": "application/json",
                    "X-Vina-Worker-Callback-Token": callback_token,
                },
                method="POST",
            )
            try:
                # Bypass desktop HTTP proxy settings for the private
                # Tailscale callback address.
                opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
                opener.open(request, timeout=2).read()
            except (urllib.error.URLError, TimeoutError, OSError):
                # Progress is best-effort; a callback outage must not cancel
                # a valid Vina calculation or its final artifact bundle.
                pass

        result = run_batch(output_root, config, emit=report_progress)
        return {
            "status": result.get("status", "failed"),
            "result": result,
            "bundle": {"root": root_name, "data": _bundle(output_root, root_name)},
        }
    except Exception as exc:
        return {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}


class ReusableThreadingHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = True
    daemon_threads = True


class Handler(BaseHTTPRequestHandler):
    server_version = "MarineVinaWorker/1.0"

    def _send(self, code: int, payload: dict) -> None:
        body = _json_bytes(payload)
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def do_GET(self) -> None:  # noqa: N802
        snapshot = re.fullmatch(r'/v1/jobs/([0-9a-f-]{32,36})/artifacts', self.path)
        if snapshot:
            expected = self.server.worker_token
            if not expected or self.headers.get('X-Vina-Worker-Token', '') != expected:
                self._send(401, {'status': 'unauthorized'})
                return
            job_id = snapshot[1]
            if job_status(job_id)['active']:
                self._send(409, {'status': 'still_active'})
                return
            root_name = f'results/vina-runs/{job_id}'
            directory = ROOT / root_name
            path = directory / 'result.json'
            if not path.is_file():
                path = directory / 'checkpoint.json'
            if not path.is_file():
                self._send(404, {'status': 'no_checkpoint'})
                return
            try:
                result = json.loads(path.read_text())
                self._send(200, {'status': 'succeeded', 'result': result,
                                 'bundle': {'root': root_name, 'data': _bundle(directory, root_name)}})
            except (OSError, ValueError):
                self._send(503, {'status': 'snapshot_unavailable'})
            return
        job = re.fullmatch(r'/v1/jobs/([0-9a-f-]{32,36})', self.path)
        if job:
            expected = self.server.worker_token
            if expected and self.headers.get('X-Vina-Worker-Token', '') != expected:
                self._send(401, {'status': 'unauthorized'})
                return
            self._send(200, job_status(job[1]))
            return
        if self.path == "/health":
            self._send(200, {"status": "ok", "service": "autodock-vina-worker", "binary": VINA_BINARY,
                             "capabilities": ["vina", "receptor_preparation", "protein_homology_screen"]})
            return
        self._send(404, {"status": "not_found"})

    def do_POST(self) -> None:  # noqa: N802
        cancelling = re.fullmatch(r"/v1/jobs/([0-9a-f-]{32,36})/cancel", self.path)
        if not cancelling and self.path not in {"/v1/vina", "/v1/proteins/screen", "/v1/receptor/inspect", "/v1/receptor/prepare"}:
            self._send(404, {"status": "not_found"})
            return
        expected = self.server.worker_token  # type: ignore[attr-defined]
        if expected and self.headers.get("X-Vina-Worker-Token", "") != expected:
            self._send(401, {"status": "unauthorized"})
            return
        if cancelling:
            try:
                control_for(cancelling[1]).cancel()
                self._send(200, {"status":"cancelling"})
            except ValueError:
                self._send(400, {"status":"invalid_input"})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 90 * 1024 * 1024:
                raise ValueError("request body is too large or empty")
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("request must be a JSON object")
        except (ValueError, json.JSONDecodeError, UnicodeDecodeError) as exc:
            self._send(400, {"status": "invalid_input", "reason": str(exc)})
            return
        if self.path == '/v1/proteins/screen':
            from marine_agent.protein_screening import run_screen
            job_id = str(payload.get('job_id', ''))
            def progress(event):
                callback = os.getenv('VINA_CALLBACK_URL', '').rstrip('/')
                if not callback:
                    return
                request = urllib.request.Request(callback + '/internal/protein-progress/' + job_id,
                    data=_json_bytes(event), method='POST', headers={'Content-Type': 'application/json',
                    'X-Vina-Worker-Callback-Token': os.getenv('VINA_CALLBACK_TOKEN', '')})
                try:
                    urllib.request.build_opener(urllib.request.ProxyHandler({})).open(request, timeout=2).close()
                except (urllib.error.URLError, TimeoutError, OSError):
                    pass
            try:
                with job_scope(job_id) as control, job_lock(JOB_LOCK, control):
                    scientific = run_screen(payload.get('request', {}), emit=progress)
                    root_name = f"results/protein-screenings/{scientific['analysis_id']}"
                    result = {'status': 'succeeded', 'result': scientific,
                              'bundle': {'root': root_name, 'data': _bundle(ROOT / root_name, root_name)}}
            except RunCancelled:
                result = {'status': 'stopped', 'result': {'status': 'stopped'}}
            except Exception as exc:
                result = {'status': 'failed', 'reason': f'{type(exc).__name__}: {exc}'}
        elif self.path.startswith('/v1/receptor/'):
            from marine_agent.receptor_tools import inspect_source, prepare_receptor
            try:
                if self.path.endswith('/inspect'):
                    result = {"status": "succeeded", "inspection": inspect_source(payload)}
                else:
                    job_id = str(payload.get('job_id', ''))
                    if not JOB_RE.fullmatch(job_id):
                        raise ValueError('job_id is invalid')
                    root_name = f'results/vina-receptors/{job_id}'
                    # Do not wait behind a long Vina batch after the HTTP
                    # client has timed out: that would create an orphan prep.
                    if not JOB_LOCK.acquire(blocking=False):
                        self._send(409, {"status": "busy", "reason": "本机正在执行分子对接或受体准备，请完成后重试"})
                        return
                    try:
                        result = {"status": "succeeded", "result": prepare_receptor(payload, ROOT / root_name)}
                        result['bundle'] = {"root": root_name, "data": _bundle(ROOT / root_name, root_name)}
                    finally:
                        JOB_LOCK.release()
            except Exception as exc:
                result = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
        else:
            try:
                with job_scope(str(payload.get('job_id',''))) as control, job_lock(JOB_LOCK, control):
                    result = _run_job(payload, self.server.job_timeout)  # type: ignore[attr-defined]
            except RunCancelled:
                result = {"status":"stopped", "result":{"status":"stopped", "results":[]}}
            except ValueError as exc:
                self._send(400, {"status":"invalid_input","reason":str(exc)})
                return
        self._send(200, result)

    def log_message(self, format: str, *args: object) -> None:
        print(f"[vina-worker] {format % args}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.environ.get("VINA_WORKER_HOST", "100.64.0.22"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("VINA_WORKER_PORT", "8100")))
    parser.add_argument("--token", default=os.environ.get("VINA_WORKER_TOKEN", ""))
    parser.add_argument("--token-file", type=Path)
    parser.add_argument("--timeout", type=int, default=int(os.environ.get("VINA_WORKER_TIMEOUT_S", "7200")))
    args = parser.parse_args()
    if args.token_file:
        args.token = args.token_file.read_text(encoding="utf-8").strip()
    server = ReusableThreadingHTTPServer((args.host, args.port), Handler)
    server.worker_token = args.token  # type: ignore[attr-defined]
    server.job_timeout = max(30, args.timeout)  # type: ignore[attr-defined]
    print(f"AutoDock Vina worker listening on {args.host}:{args.port}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
