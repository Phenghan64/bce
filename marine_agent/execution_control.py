"""Cooperative task cancellation and bounded termination of owned process groups."""
from __future__ import annotations

import contextvars
import os
import re
import signal
import subprocess
import threading
import time
from contextlib import contextmanager

JOB_ID = re.compile(r"^(?:[0-9a-f]{32}|[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12})$", re.I)
TERMINAL = {"succeeded", "partial", "failed", "stopped"}


class RunCancelled(BaseException):
    """Not a scientific failure; tool error handlers must not consume it."""


class JobControl:
    def __init__(self, job_id: str):
        self.job_id = job_id
        self.event = threading.Event()
        self.lock = threading.Lock()
        self.processes: set[subprocess.Popen] = set()
        self.active_scopes = 0

    @property
    def cancelled(self):
        return self.event.is_set()

    def check(self):
        if self.cancelled:
            raise RunCancelled()

    @staticmethod
    def _signal(process, sig):
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            pass

    def cancel(self):
        with self.lock:
            self.event.set()
            for process in self.processes:
                self._signal(process, signal.SIGTERM)

    def run(self, command, *, timeout=None, **kwargs):
        if kwargs.pop('capture_output', False):
            kwargs['stdout'] = subprocess.PIPE
            kwargs['stderr'] = subprocess.PIPE
        # Descendants (MSA and MLX included) inherit this dedicated group.
        with self.lock:
            self.check()
            process = subprocess.Popen(command, start_new_session=True, **kwargs)
            self.processes.add(process)
        started = time.monotonic()
        interrupted_at = None
        timed_out = False
        try:
            while True:
                expired = timeout is not None and time.monotonic() - started >= timeout
                if self.cancelled or expired:
                    timed_out = timed_out or expired
                    if interrupted_at is None:
                        self._signal(process, signal.SIGTERM)
                        interrupted_at = time.monotonic()
                    elif time.monotonic() - interrupted_at >= 2:
                        self._signal(process, signal.SIGKILL)
                try:
                    stdout, stderr = process.communicate(timeout=.2)
                    break
                except subprocess.TimeoutExpired:
                    continue
            if interrupted_at is not None:
                # Reap descendants that outlived the wrapper, including ones
                # which closed their output pipes before exiting.
                self._signal(process, signal.SIGKILL)
            self.check()
            if timed_out:
                raise subprocess.TimeoutExpired(command, timeout, output=stdout, stderr=stderr)
            return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
        finally:
            if process.poll() is None:
                self._signal(process, signal.SIGKILL)
                process.wait()
            with self.lock:
                self.processes.discard(process)


_lock = threading.Lock()
_controls: dict[str, JobControl] = {}
_current = contextvars.ContextVar('scientific_job_control', default=None)


def control_for(job_id: str) -> JobControl:
    if not JOB_ID.fullmatch(job_id):
        raise ValueError('invalid cancellation job ID')
    with _lock:
        # Keep cancellation tombstones: a delayed HTTP submission must not
        # start a task after its cancellation has already been acknowledged.
        return _controls.setdefault(job_id, JobControl(job_id))


def current_control() -> JobControl | None:
    return _current.get()


def job_status(job_id: str) -> dict:
    """Read without creating a tombstone; 'idle' does not certify success."""
    if not JOB_ID.fullmatch(job_id):
        raise ValueError('invalid job ID')
    with _lock:
        control = _controls.get(job_id)
    if control is None:
        return {'job_id': job_id, 'state': 'unknown', 'active': False}
    with control.lock:
        active = bool(control.active_scopes or control.processes)
        return {'job_id': job_id, 'active': active,
                'state': ('cancelling' if control.cancelled else 'running') if active else
                         ('stopped' if control.cancelled else 'idle')}


def check_cancelled():
    if control := current_control():
        control.check()


@contextmanager
def job_scope(job_id: str):
    control = control_for(job_id)
    with control.lock:
        control.active_scopes += 1
    token = _current.set(control)
    try:
        yield control
    finally:
        _current.reset(token)
        with control.lock:
            control.active_scopes -= 1


@contextmanager
def job_lock(lock, control):
    acquired = False
    try:
        while not acquired:
            control.check()
            acquired = lock.acquire(timeout=.2)
        control.check()
        yield
    finally:
        if acquired:
            lock.release()
