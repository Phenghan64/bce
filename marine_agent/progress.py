"""Execution-scoped progress transport, separate from model-written prose."""
from contextvars import ContextVar
reporter = ContextVar('scientific_progress_reporter', default=None)
def emit(stage):
    callback = reporter.get()
    if callback:
        callback(stage)
