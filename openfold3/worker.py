#!/usr/bin/env python3
"""Small Tailscale-only OpenFold3 execution worker for the Marine Agent API.

The control-plane API runs on Linux, while this worker runs on the Mac that
owns the Apple MPS OpenFold3 environment. It deliberately exposes one
authenticated, bounded endpoint and returns a reproducible artifact bundle;
it does not expose a shell or arbitrary filesystem access.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import re
import subprocess
import tarfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from pathlib import Path
from threading import Lock
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from marine_agent.execution_control import RunCancelled, control_for, current_control, job_scope, job_lock, job_status


ROOT = Path(os.environ.get("OPENFOLD_PROJECT_ROOT", str(Path(__file__).resolve().parents[1]))).resolve()
JOB_RE = re.compile(r"^[0-9a-f]{32}$")
ACCESSION_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{1,63}$")
AA = set("ACDEFGHIKLMNPQRSTVWY")
JOB_LOCK = Lock()


class ReusableThreadingHTTPServer(ThreadingHTTPServer):
    """Allow launchd to restart the worker without waiting on a stale socket."""

    allow_reuse_address = True
    daemon_threads = True


def _json_bytes(payload: dict) -> bytes:
    return json.dumps(payload, ensure_ascii=False).encode("utf-8")


def _bundle(output_root: Path, root_name: str) -> str:
    buffer = BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        archive.add(output_root, arcname=root_name, recursive=True)
    return base64.b64encode(buffer.getvalue()).decode("ascii")


def _run_job(payload: dict, timeout: int) -> dict:
    job_id = str(payload.get("job_id", ""))
    accession = str(payload.get("accession", "")).strip()
    sequence = re.sub(r"\s+", "", str(payload.get("sequence", ""))).upper()
    sequence_id = str(payload.get("sequence_id", "sequence")).strip()
    root_name = str(payload.get("output_root", ""))
    ligands = payload.get("ligands", [])
    if not JOB_RE.fullmatch(job_id):
        return {"status": "invalid_input", "reason": "job_id is invalid"}
    has_accession = bool(accession)
    has_sequence = bool(sequence)
    if has_accession == has_sequence:
        return {"status": "invalid_input", "reason": "provide exactly one of accession or sequence"}
    if has_accession and not ACCESSION_RE.fullmatch(accession):
        return {"status": "invalid_input", "reason": "accession is invalid"}
    if has_sequence:
        if len(sequence) < 20 or len(sequence) > 4000:
            return {"status": "invalid_input", "reason": "sequence must contain 20-4000 amino acids"}
        invalid = sorted(set(sequence) - AA)
        if invalid:
            return {"status": "invalid_input", "reason": f"sequence contains unsupported residues: {','.join(invalid)}"}
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", sequence_id):
            return {"status": "invalid_input", "reason": "sequence_id is invalid"}
    if root_name != f"results/openfold3-remote/{job_id}":
        return {"status": "invalid_input", "reason": "output_root is invalid"}
    if not isinstance(ligands, list) or len(ligands) > 16:
        return {"status": "invalid_input", "reason": "ligands must be a list of at most 16 entries"}

    output_root = (ROOT / root_name).resolve()
    try:
        output_root.relative_to((ROOT / "results/openfold3-remote").resolve())
    except ValueError:
        return {"status": "invalid_input", "reason": "output_root escaped the worker result directory"}
    output_root.mkdir(parents=True, exist_ok=True)
    fasta_path = None
    if has_sequence:
        fasta_path = output_root / "input.fasta"
        fasta_path.write_text(f">{sequence_id}\n{sequence}\n", encoding="utf-8")
    command = [
        str(ROOT / "openfold3/run.sh"),
        "python",
        "-u",
        str(ROOT / "openfold3/predict_accession.py"),
        "--fasta" if has_sequence else "--accession",
        str(fasta_path) if has_sequence else accession,
        "--output-dir",
        str(output_root),
        "--ligands-json",
        json.dumps(ligands, ensure_ascii=False, separators=(",", ":")),
    ]
    try:
        control = current_control()
        completed = (control.run if control else subprocess.run)(
            command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
        stdout = completed.stdout or ""
        try:
            result = json.loads(stdout)
        except json.JSONDecodeError:
            result = {
                "status": "failed",
                "stage": "worker_process",
                "reason": "OpenFold3 process did not return JSON",
                "output_tail": stdout[-4000:],
            }
        if not isinstance(result, dict):
            result = {"status": "failed", "stage": "worker_process", "reason": "OpenFold3 returned a non-object result"}
        result.setdefault("exit_code", completed.returncode)
    except RunCancelled:
        result = {"status":"stopped", "stage":"stopped", "reason":"用户已停止结构预测"}
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        result = {"status": "timed_out", "stage": "worker_process", "reason": "OpenFold3 worker timeout", "output_tail": output[-4000:]}
    except OSError as exc:
        result = {"status": "failed", "stage": "worker_process", "reason": str(exc)}

    bundle = _bundle(output_root, root_name)
    return {"status": result.get("status", "failed"), "result": result, "bundle": {"root": root_name, "data": bundle}}


class Handler(BaseHTTPRequestHandler):
    server_version = "MarineOpenFoldWorker/1.0"

    def _send(self, code: int, payload: dict) -> None:
        body = _json_bytes(payload)
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            # A browser or control-plane request may be cancelled while a
            # long OpenFold3 job is still finishing. The worker must remain
            # alive so launchd does not restart it for a client disconnect.
            pass

    def do_GET(self) -> None:  # noqa: N802
        progress = re.fullmatch(r'/v1/jobs/([0-9a-f-]{32,36})/progress', self.path)
        if progress:
            expected = self.server.worker_token
            if not expected or self.headers.get('X-OpenFold-Worker-Token', '') != expected:
                self._send(401, {'status': 'unauthorized'})
                return
            path = ROOT / 'results/openfold3-remote' / progress[1] / 'progress.json'
            try:
                result = json.loads(path.read_text())
            except (ValueError, OSError):
                result = {'stage': 'sequence'}
            self._send(200, {'job_id': progress[1], **result})
            return
        job = re.fullmatch(r'/v1/jobs/([0-9a-f-]{32,36})', self.path)
        if job:
            expected = self.server.worker_token
            if expected and self.headers.get('X-OpenFold-Worker-Token', '') != expected:
                self._send(401, {'status': 'unauthorized'})
                return
            self._send(200, job_status(job[1]))
            return
        if self.path == "/health":
            self._send(200, {"status": "ok", "service": "openfold3-mps-worker", "device": "Apple MPS"})
            return
        self._send(404, {"status": "not_found"})

    def do_POST(self) -> None:  # noqa: N802
        cancelling = re.fullmatch(r"/v1/jobs/([0-9a-f-]{32,36})/cancel", self.path)
        if self.path != "/v1/openfold" and not cancelling:
            self._send(404, {"status": "not_found"})
            return
        expected = self.server.worker_token  # type: ignore[attr-defined]
        if expected and self.headers.get("X-OpenFold-Worker-Token", "") != expected:
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
            if length <= 0 or length > 128 * 1024:
                raise ValueError("request body is too large or empty")
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("request must be a JSON object")
        except (ValueError, json.JSONDecodeError, UnicodeDecodeError) as exc:
            self._send(400, {"status": "invalid_input", "reason": str(exc)})
            return
        try:
            with job_scope(str(payload.get('cancellation_id') or payload.get('job_id',''))) as control, job_lock(JOB_LOCK, control):
                result = _run_job(payload, self.server.job_timeout)  # type: ignore[attr-defined]
        except RunCancelled:
            result = {"status":"stopped", "result":{"status":"stopped"}}
        except ValueError as exc:
            self._send(400, {"status":"invalid_input","reason":str(exc)})
            return
        self._send(200, result)

    def log_message(self, format: str, *args: object) -> None:
        print(f"[openfold-worker] {format % args}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.environ.get("OPENFOLD_WORKER_HOST", "100.64.0.22"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("OPENFOLD_WORKER_PORT", "8099")))
    parser.add_argument("--token", default=os.environ.get("OPENFOLD_WORKER_TOKEN", ""))
    parser.add_argument("--token-file", type=Path, help="read the shared token from a mode-600 file")
    parser.add_argument("--timeout", type=int, default=int(os.environ.get("OPENFOLD_WORKER_TIMEOUT_S", "3600")))
    args = parser.parse_args()
    if args.token_file:
        args.token = args.token_file.read_text(encoding="utf-8").strip()
    server = ReusableThreadingHTTPServer((args.host, args.port), Handler)
    server.worker_token = args.token  # type: ignore[attr-defined]
    server.job_timeout = max(60, args.timeout)  # type: ignore[attr-defined]
    print(f"OpenFold3 MPS worker listening on {args.host}:{args.port}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
