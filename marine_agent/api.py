"""HTTP control plane for the bounded marine scientific Deep Agent.

The API owns task state and audit events. Scientific execution remains inside
the registered tools exposed by :mod:`marine_agent.agent`; no arbitrary shell
or filesystem tool is made available to the model.
"""
from __future__ import annotations

from contextlib import contextmanager

import base64
import copy
import hashlib
import json
import os
import re
import sqlite3
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from fastapi import FastAPI, Header, HTTPException, Query, status
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

from .agent import COLLAGENASE_SYSTEM_PROMPT, build_agent
from .collagenase_tools import COLLAGENASE_TOOLS
from .agent import TOOLS as CORE_TOOLS
from .intent import classify_intent
from .openfold_tools import get_uniprot_sequence, run_openfold_accession, run_openfold_sequence, search_uniprot, _parse_ligands_json
from .vina_tools import VinaInputError, run_batch, validate_config, _seed_tasks, _analyze_seed_consistency
from .ligand_library import LibraryUnavailable, public_catalog, resolve_library_ligands
from .vina_remote import VINA_WORKER_URL, run_vina_worker, recover_vina_artifacts
from .receptor_api import register_receptor_routes, resolve_prepared_receptor
from .execution_control import RunCancelled, TERMINAL, control_for, job_scope, check_cancelled
from .task_lifecycle import ConversationTasks, saved_spec, describe_spec
from .turn_planning import plan_turn, CONVERSATION_PROMPT, TurnPlan, TurnPlanningError
from .protein_tools import execute_screen
from .protein_screening import ScreenRequest, validate_request
from .domain_evidence import DomainRequest, annotate, resolve_inputs, read_domain_analysis
from .evidence_review import review_answer, render_domains, method_rejection
from .result_contract import trusted_structure_manifest, verified_structure, execution_receipt, scientific_payloads
from .restart_recovery import reconcile
from .access_control import AccountMiddleware, principal, account_scope
from .research_projects import ProjectStore, register_project_routes, project_context_message, project_tools
from urllib.request import Request, build_opener, ProxyHandler


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
STATE_DB = Path(os.getenv("MARINE_AGENT_STATE_DB", str(ROOT / "results/marine-agent-api/state.sqlite3")))
MAX_CONCURRENCY = max(1, int(os.getenv("MARINE_AGENT_MAX_CONCURRENCY", "1")))
CONTROL_CONCURRENCY = max(2, int(os.getenv('MARINE_AGENT_CONTROL_CONCURRENCY', '4')))
RECURSION_LIMIT = max(4, int(os.getenv("DEEP_AGENT_RECURSION_LIMIT", "24")))
MODEL_NAME = os.getenv("DEEP_AGENT_MODEL", "configured-provider")
VINA_CALLBACK_TOKEN = os.getenv("VINA_CALLBACK_TOKEN", "")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    model_dump = getattr(value, "model_dump", None)
    if callable(model_dump):
        return _json_safe(model_dump())
    return str(value)


def _content_text(content: Any) -> str:
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts: list[str] = []
        for item in content:
            if isinstance(item, str):
                parts.append(item)
            elif isinstance(item, dict) and isinstance(item.get("text"), str):
                parts.append(item["text"])
        return "".join(parts)
    return str(content) if content is not None else ""


def _final_text(result: Any) -> str:
    """Extract the final assistant text without assuming a provider message type."""
    messages = result.get("messages", []) if isinstance(result, dict) else []
    for message in reversed(messages):
        if isinstance(message, dict):
            role = message.get("role") or message.get("type")
            content = message.get("content")
        else:
            role = getattr(message, "role", None) or getattr(message, "type", None)
            content = getattr(message, "content", None)
        if role in ("assistant", "ai") and content:
            return _content_text(content)
    return ""


def _message_field(message: Any, field: str, default: Any = None) -> Any:
    """Read a LangChain message field from either an object or JSON-like dict."""
    if isinstance(message, dict):
        return message.get(field, default)
    return getattr(message, field, default)


def _stream_event_payload(message: Any) -> tuple[str | None, dict[str, Any]] | None:
    """Convert one newly observed agent message into a small audit event.

    The full message remains in the final run record. Live events are deliberately
    bounded so the browser can render progress without receiving giant tool
    payloads or protein sequences.
    """
    role = _message_field(message, "role") or _message_field(message, "type")
    tool_calls = _message_field(message, "tool_calls") or []
    if role in ("assistant", "ai") and tool_calls:
        call = tool_calls[0]
        if isinstance(call, dict):
            name = call.get("name") or "unknown_tool"
            args = call.get("args") or {}
        else:
            name = getattr(call, "name", None) or "unknown_tool"
            args = getattr(call, "args", {}) or {}
        return "tool_call", {
            "name": f"调用工具：{name}",
            "tool": name,
            "args_preview": str(_json_safe(args))[:600],
        }
    if role in ("tool", "tool_result"):
        name = _message_field(message, "name") or _message_field(message, "tool_name") or "tool"
        content = _content_text(_message_field(message, "content"))
        return "tool_result", {
            "name": f"工具完成：{name}",
            "tool": name,
            "preview": content[:800],
        }
    return None


class ManualOpenFoldInput(BaseModel):
    accession: str = Field(min_length=2, max_length=20)
    ligands: list[dict[str, Any]] = Field(default_factory=list, max_length=16)


class SubmissionConflict(ValueError):
    pass


class VinaRunRequest(BaseModel):
    """Explicit, browser-safe Vina workbench request.

    Files are base64 PDBQT text, or a server-verified receptor preparation ID.
    The API never receives arbitrary host paths. Ligand chemistry is prepared
    separately; OpenFold protein preparation records its standard pH treatment.
    """

    request_id: str | None = Field(default=None, min_length=8, max_length=128)
    receptor: dict[str, Any]
    ligands: list[dict[str, Any]] = Field(min_length=1, max_length=100)
    box: dict[str, Any]
    params: dict[str, Any] = Field(default_factory=dict)
    conversation_id: str | None = Field(default=None, min_length=8, max_length=64)


class RunRequest(BaseModel):
    request_id: str | None = Field(default=None, min_length=8, max_length=128)
    prompt: str = Field(min_length=1, max_length=12_000)
    mode: Literal["auto", "general", "collagenase", "manual_openfold"] = "auto"
    manual: ManualOpenFoldInput | None = None
    conversation_id: str | None = Field(default=None, min_length=8, max_length=64)


class ProteinScreenRunRequest(ScreenRequest):
    conversation_id: str | None = Field(default=None, min_length=8, max_length=64)


class DomainRunRequest(DomainRequest):
    conversation_id: str | None = Field(default=None, min_length=8, max_length=64)


class RunResponse(BaseModel):
    run_id: str
    status: str
    created_at: str
    conversation_id: str | None = None


class RunStore(ProjectStore):
    """Small durable task store suitable for the control plane.

    Scientific artifacts remain in the workflow result directories. SQLite
    stores only task metadata, event history, and the final agent response.
    """

    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._initialize()
        self.initialize_projects()
        self.initialize_candidates()

    @contextmanager
    def _connect(self):
        connection = sqlite3.connect(self.path, timeout=30)
        try:
            connection.row_factory = sqlite3.Row
            connection.execute("PRAGMA busy_timeout=30000")
            with connection:
                yield connection
        finally:
            connection.close()

    def grant_analysis(self, kind, analysis_id, owner_id):
        with self._lock, self._connect() as db:
            db.execute('INSERT OR IGNORE INTO analysis_owners VALUES (?,?,?)', (kind, analysis_id, owner_id))

    def owns_analysis(self, kind, analysis_id, owner_id):
        with self._connect() as db:
            return db.execute('SELECT 1 FROM analysis_owners WHERE kind=? AND analysis_id=? AND owner_id=?', (kind,analysis_id,owner_id)).fetchone() is not None

    def claim_orphans(self) -> list[dict]:
        # Called only by the single API process startup, never by RunStore's
        # constructor (tests/readers must not interrupt a live coordinator).
        with self._lock, self._connect() as db:
            ids = [row['id'] for row in db.execute("SELECT id FROM runs WHERE status IN ('queued','running','cancelling')")]
            for run_id in ids:
                db.execute("UPDATE runs SET status='cancelling', current_stage='recovery_pending' WHERE id=?", (run_id,))
                self._event(db, run_id, 'recovery_pending', {'name': '服务重启，核对遗留任务状态'})
        return [self.get(run_id) for run_id in ids]

    def save_recovered_vina(self, run_id, payload):
        with self._lock, self._connect() as db:
            row = db.execute('SELECT status, result_json FROM runs WHERE id=?', (run_id,)).fetchone()
            if not row or row['status'] not in ('cancelling', 'stopped', 'failed', 'partial'):
                return
            result = json.loads(row['result_json'] or '{}')
            result['vina'] = payload
            db.execute('UPDATE runs SET result_json=? WHERE id=?', (json.dumps(result, ensure_ascii=False), run_id))
            self._event(db, run_id, 'artifacts_recovered', {'name': '已核验本机保存的对接产物', 'completed': sum(r.get('status') == 'succeeded' for r in payload.get('results', []))})

    def finish_recovery(self, run_id: str) -> None:
        original = self.get(run_id) or {}
        text = '服务重启中断了本轮任务，已确认旧作业停止。已有输入和结果保留。'
        text += '可以沿用原输入继续。' if saved_spec(original) else '请明确要执行的具体操作后重新提交。'
        with self._lock, self._connect() as db:
            row = db.execute('SELECT * FROM runs WHERE id=?', (run_id,)).fetchone()
            if not row or row['current_stage'] != 'recovery_pending' or row['status'] != 'cancelling':
                return
            result = json.loads(row['result_json']) if row['result_json'] else {}
            result.setdefault('final_text', text)
            result['recovery'] = {'reason': 'coordinator_restarted', 'remote_quiescence_confirmed': True}
            db.execute("UPDATE runs SET status='stopped', current_stage='stopped', completed_at=?, error=NULL, result_json=? WHERE id=?",
                       (_now(), json.dumps(result, ensure_ascii=False), run_id))
            db.execute("UPDATE scientific_operations SET status='stopped' WHERE run_id=? AND status='running'", (run_id,))
            self._event(db, run_id, 'stopped', {'name': '重启恢复检查完成，旧作业已停止'})
        self.append_message(row['conversation_id'], 'assistant', text, run_id)

    def _initialize(self) -> None:
        with self._connect() as db:
            db.execute("PRAGMA journal_mode=WAL")
            db.executescript(
                """
                CREATE TABLE IF NOT EXISTS migrations (name TEXT PRIMARY KEY);
                CREATE TABLE IF NOT EXISTS analysis_owners (kind TEXT, analysis_id TEXT, owner_id TEXT, PRIMARY KEY(kind,analysis_id,owner_id));
                CREATE TABLE IF NOT EXISTS task_requests (owner_id TEXT, request_id TEXT, fingerprint TEXT, run_id TEXT, PRIMARY KEY(owner_id,request_id));
                CREATE TABLE IF NOT EXISTS runs (
                    id TEXT PRIMARY KEY,
                    conversation_id TEXT,
                    prompt TEXT NOT NULL,
                    mode TEXT NOT NULL,
                    status TEXT NOT NULL,
                    model TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    started_at TEXT,
                    completed_at TEXT,
                    current_stage TEXT NOT NULL DEFAULT 'queued',
                    result_json TEXT,
                    manual_json TEXT,
                    error TEXT,
                    runtime_seconds REAL
                );
                CREATE TABLE IF NOT EXISTS run_events (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    run_id TEXT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
                    created_at TEXT NOT NULL,
                    event_type TEXT NOT NULL,
                    payload_json TEXT NOT NULL
                );
                CREATE INDEX IF NOT EXISTS idx_run_events_run_id ON run_events(run_id, id);
                CREATE TABLE IF NOT EXISTS scientific_operations (
                    run_id TEXT NOT NULL,
                    operation_id TEXT NOT NULL,
                    tool TEXT NOT NULL,
                    spec_json TEXT NOT NULL,
                    status TEXT NOT NULL,
                    facts_json TEXT,
                    PRIMARY KEY (run_id, operation_id)
                );
                CREATE TABLE IF NOT EXISTS conversations (
                    id TEXT PRIMARY KEY,
                    title TEXT NOT NULL DEFAULT '新的对话',
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS conversation_messages (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    conversation_id TEXT NOT NULL REFERENCES conversations(id) ON DELETE CASCADE,
                    role TEXT NOT NULL CHECK (role IN ('user', 'assistant')),
                    content TEXT NOT NULL,
                    run_id TEXT,
                    created_at TEXT NOT NULL
                );
                CREATE INDEX IF NOT EXISTS idx_conversation_messages_conversation
                    ON conversation_messages(conversation_id, id);
                """
            )
            columns = {row[1] for row in db.execute("PRAGMA table_info(runs)").fetchall()}
            if "manual_json" not in columns:
                db.execute("ALTER TABLE runs ADD COLUMN manual_json TEXT")
            if "conversation_id" not in columns:
                db.execute("ALTER TABLE runs ADD COLUMN conversation_id TEXT")
            if "execution_json" not in columns:
                db.execute("ALTER TABLE runs ADD COLUMN execution_json TEXT")
            if "continuation_operation_id" not in columns:
                db.execute("ALTER TABLE runs ADD COLUMN continuation_operation_id TEXT")
            if "continuation_of" not in columns:
                db.execute("ALTER TABLE runs ADD COLUMN continuation_of TEXT")
            for table in ('runs', 'conversations'):
                if 'owner_id' not in {row[1] for row in db.execute(f'PRAGMA table_info({table})')}:
                    db.execute(f"ALTER TABLE {table} ADD COLUMN owner_id TEXT NOT NULL DEFAULT ''")
                db.execute(f'CREATE INDEX IF NOT EXISTS idx_{table}_owner ON {table}(owner_id)')
            legacy_owner = os.getenv('MARINE_LEGACY_OWNER_ID', '')
            migrate_owners = legacy_owner and not db.execute("SELECT 1 FROM migrations WHERE name='account_ownership_v1'").fetchone()
            if migrate_owners:
                for table in ('runs', 'conversations'):
                    db.execute(f"UPDATE {table} SET owner_id=? WHERE owner_id=''", (legacy_owner,))
            db.execute("CREATE INDEX IF NOT EXISTS idx_runs_conversation ON runs(conversation_id, created_at)")
            # Preserve pre-conversation task history. Older runs become
            # one-message conversations so they remain visible in the sidebar;
            # all new messages reuse an explicit conversation ID.
            legacy_runs = db.execute(
                "SELECT id, prompt, created_at, result_json, error FROM runs WHERE conversation_id IS NULL ORDER BY created_at"
            ).fetchall()
            for legacy in legacy_runs:
                conversation_id = str(uuid.uuid4())
                prompt = str(legacy["prompt"] or "新的对话")
                title = " ".join(prompt.strip().split())[:36] or "新的对话"
                db.execute(
                    "INSERT INTO conversations (id, title, created_at, updated_at) VALUES (?, ?, ?, ?)",
                    (conversation_id, title, legacy["created_at"], legacy["created_at"]),
                )
                db.execute("UPDATE runs SET conversation_id = ? WHERE id = ?", (conversation_id, legacy["id"]))
                db.execute(
                    "INSERT INTO conversation_messages (conversation_id, role, content, run_id, created_at) VALUES (?, 'user', ?, ?, ?)",
                    (conversation_id, prompt, legacy["id"], legacy["created_at"]),
                )
                assistant_text = None
                if legacy["result_json"]:
                    try:
                        assistant_text = json.loads(legacy["result_json"]).get("final_text")
                    except (TypeError, json.JSONDecodeError):
                        assistant_text = None
                assistant_text = assistant_text or (f"任务执行失败：{legacy['error']}" if legacy["error"] else None)
                if assistant_text:
                    db.execute(
                        "INSERT INTO conversation_messages (conversation_id, role, content, run_id, created_at) VALUES (?, 'assistant', ?, ?, ?)",
                        (conversation_id, str(assistant_text), legacy["id"], legacy["created_at"]),
                    )


            if migrate_owners:
                # Only provenance-bearing result fields are migrated. Model prose
                # and user-pasted analysis IDs never grant access.
                for row in db.execute('SELECT id, owner_id, result_json FROM runs WHERE owner_id=?', (legacy_owner,)):
                    agent = json.loads(row['result_json'] or '{}').get('agent') or {}
                    values = [(kind,agent.get(key,{})) for key,kind in [('protein_domains','protein-domains'),('protein_screen','protein-screenings')]]
                    for message in agent.get('messages', []):
                        if not isinstance(message, dict) or message.get('type', message.get('role')) not in ('tool', 'tool_result'):
                            continue
                        kind = {'annotate_protein_domains':'protein-domains', 'run_protein_homology_screen':'protein-screenings'}.get(message.get('name'))
                        if kind:
                            try:
                                value = json.loads(message.get('content') or '{}')
                            except (ValueError, TypeError):
                                continue
                            if isinstance(value, dict):
                                values.append((kind, value))
                    for op in db.execute('SELECT tool, facts_json FROM scientific_operations WHERE run_id=?',(row['id'],)):
                        kind = 'protein-domains' if op['tool']=='annotate_protein_domains' else 'protein-screenings' if op['tool']=='run_protein_homology_screen' else None
                        if kind:
                            values.append((kind,json.loads(op['facts_json'] or '{}')))
                    for kind,value in values:
                        analysis_id=value.get('analysis_id','')
                        if re.fullmatch('[a-f0-9]{64}',analysis_id):
                            db.execute('INSERT OR IGNORE INTO analysis_owners VALUES (?,?,?)',(kind,analysis_id,legacy_owner))
                db.execute("UPDATE conversations SET owner_id=? WHERE owner_id=''", (legacy_owner,))
                db.execute("INSERT INTO migrations VALUES ('account_ownership_v1')")

    def create(self, prompt: str, mode: str, manual: dict[str, Any] | None = None,
               conversation_id: str | None = None, *, persist_prompt: bool = True, request_id: str | None = None, fingerprint: str | None = None, require_idle: bool = False) -> dict[str, Any]:
        run_id = str(uuid.uuid4())
        created_at = _now()
        owner_id = principal.get() or ''
        with self._lock, self._connect() as db:
            if request_id:
                previous = db.execute('SELECT * FROM task_requests WHERE owner_id=? AND request_id=?', (owner_id,request_id)).fetchone()
                if previous:
                    if previous['fingerprint'] != fingerprint:
                        raise SubmissionConflict('重复请求编号对应不同输入，请重新提交。')
                    saved = db.execute('SELECT id AS run_id,status,created_at,conversation_id FROM runs WHERE id=? AND owner_id=?',(previous['run_id'],owner_id)).fetchone()
                    if saved: return {**dict(saved), '_reused':True}
            if require_idle and conversation_id:
                active = db.execute("SELECT id FROM runs WHERE conversation_id=? AND owner_id=? AND status IN ('queued','running','cancelling') LIMIT 1",(conversation_id,owner_id)).fetchone()
                if active: raise SubmissionConflict('当前会话已有任务在执行，请查看或停止该任务后再提交。')
            if persist_prompt:
                if conversation_id:
                    exists = db.execute("SELECT 1 FROM conversations WHERE id = ? AND owner_id = ?", (conversation_id, owner_id)).fetchone()
                    if exists is None:
                        raise ValueError("conversation not found")
                else:
                    conversation_id = str(uuid.uuid4())
                    title = " ".join(prompt.strip().split())[:36] or "新的对话"
                    db.execute(
                        "INSERT INTO conversations (id, title, created_at, updated_at) VALUES (?, ?, ?, ?)",
                        (conversation_id, title, created_at, created_at),
                    )
            db.execute(
                "INSERT INTO runs (id, conversation_id, prompt, mode, status, model, created_at, manual_json, owner_id) VALUES (?, ?, ?, ?, 'queued', ?, ?, ?, ?)",
                (run_id, conversation_id, prompt, mode, MODEL_NAME, created_at, json.dumps(manual, ensure_ascii=False) if manual else None, owner_id),
            )
            if conversation_id:
                db.execute("UPDATE conversations SET owner_id=? WHERE id=? AND owner_id=''", (owner_id, conversation_id))
            if request_id:
                db.execute('INSERT INTO task_requests VALUES (?,?,?,?)',(owner_id,request_id,fingerprint,run_id))
            self.capture_project_context(db, run_id, conversation_id)
            self._event(db, run_id, "queued", {"mode": mode})
            if persist_prompt:
                db.execute(
                    "INSERT INTO conversation_messages (conversation_id, role, content, run_id, created_at) VALUES (?, 'user', ?, ?, ?)",
                    (conversation_id, prompt, run_id, created_at),
                )
                db.execute('UPDATE conversations SET updated_at=? WHERE id=?', (created_at, conversation_id))
        return {"run_id": run_id, "status": "queued", "created_at": created_at, "conversation_id": conversation_id}

    def append_message(self, conversation_id: str | None, role: str, content: str, run_id: str | None = None) -> None:
        if not conversation_id or role not in {"user", "assistant"} or not content.strip():
            return
        now = _now()
        with self._lock, self._connect() as db:
            if role == 'assistant' and run_id:
                state = db.execute('SELECT status FROM runs WHERE id=?', (run_id,)).fetchone()
                if state and state['status'] == 'cancelling':
                    return
            db.execute(
                "INSERT INTO conversation_messages (conversation_id, role, content, run_id, created_at) VALUES (?, ?, ?, ?, ?)",
                (conversation_id, role, content, run_id, now),
            )
            db.execute("UPDATE conversations SET updated_at = ? WHERE id = ?", (now, conversation_id))

    def save_operation(self, run_id, operation_id, tool, spec):
        with self._lock, self._connect() as db:
            db.execute('INSERT OR IGNORE INTO scientific_operations VALUES (?,?,?,?,?,NULL)',
                       (run_id, operation_id, tool, json.dumps(spec, ensure_ascii=False), 'running'))

    def complete_operation(self, run_id, operation_id, status, facts):
        with self._lock, self._connect() as db:
            db.execute('UPDATE scientific_operations SET status=?, facts_json=? WHERE run_id=? AND operation_id=?',
                       (status, json.dumps(facts, ensure_ascii=False), run_id, operation_id))

    def conversation_context(self, conversation_id: str | None, max_messages: int = 40, max_chars: int = 120_000, through_run_id: str | None = None) -> list[dict[str, str]]:
        if not conversation_id:
            return []
        with self._connect() as db:
            rows = db.execute(
                """SELECT role, content FROM conversation_messages WHERE conversation_id = ?
                   AND (? IS NULL OR id <= (SELECT min(id) FROM conversation_messages WHERE run_id=? AND role='user'))
                   ORDER BY id DESC LIMIT ?""",
                (conversation_id, through_run_id, through_run_id, max(1, min(max_messages, 100))),
            ).fetchall()
        messages: list[dict[str, str]] = []
        used = 0
        for row in rows:
            content = str(row["content"] or "")
            if used + len(content) > max_chars:
                if messages:
                    break
                content = content[-max_chars:]
            messages.append({"role": str(row["role"]), "content": content})
            used += len(content)
        return list(reversed(messages))

    def list_conversations(self, limit: int = 50) -> list[dict[str, Any]]:
        bounded_limit = max(1, min(limit, 100))
        with self._connect() as db:
            rows = db.execute(
                """
                SELECT c.id, c.title, c.created_at, c.updated_at,
                       (SELECT substr(m.content, 1, 240) FROM conversation_messages m WHERE m.conversation_id = c.id ORDER BY m.id DESC LIMIT 1) AS last_message,
                       (SELECT r.id FROM runs r WHERE r.conversation_id = c.id ORDER BY r.created_at DESC LIMIT 1) AS latest_run_id,
                       (SELECT r.status FROM runs r WHERE r.conversation_id = c.id ORDER BY r.created_at DESC LIMIT 1) AS latest_run_status,
                       (SELECT r.current_stage FROM runs r WHERE r.conversation_id = c.id ORDER BY r.created_at DESC LIMIT 1) AS latest_run_stage
                FROM conversations c WHERE (? IS NULL OR c.owner_id=?) ORDER BY c.updated_at DESC LIMIT ?
                """,
                (principal.get(), principal.get(), bounded_limit),
            ).fetchall()
        return [dict(row) for row in rows]

    def get_conversation(self, conversation_id: str) -> dict[str, Any] | None:
        with self._connect() as db:
            conversation = db.execute("SELECT * FROM conversations WHERE id = ? AND (? IS NULL OR owner_id=?)", (conversation_id, principal.get(), principal.get())).fetchone()
            if conversation is None:
                return None
            messages = db.execute(
                "SELECT id, role, content, run_id, created_at FROM conversation_messages WHERE conversation_id = ? ORDER BY id",
                (conversation_id,),
            ).fetchall()
            runs = db.execute(
                "SELECT id, conversation_id, prompt, mode, status, model, created_at, started_at, completed_at, current_stage, error, runtime_seconds, continuation_of FROM runs WHERE conversation_id = ? ORDER BY created_at",
                (conversation_id,),
            ).fetchall()
        summaries = [dict(row) for row in runs]
        continued = {row['continuation_of'] for row in summaries if row.get('continuation_of')}
        for row in summaries:
            full = self.get(row['id'])
            row['has_claim_evidence'] = bool((full.get('result') or {}).get('claim_evidence'))
            row['has_structure_result'] = bool(trusted_structure_manifest((full.get('result') or {}).get('agent'), _find_viewer_manifest))
            row['can_continue'] = row['id'] not in continued and row['status'] in {'stopped','failed','partial'} and bool(saved_spec(self.get(row['id'])))
        return {
            "conversation": dict(conversation),
            "project": self.conversation_project(conversation_id),
            "messages": [dict(row) for row in messages],
            "runs": summaries,
        }

    def can_continue(self, record: dict) -> bool:
        if record['status'] not in {'stopped', 'failed', 'partial'} or not saved_spec(record):
            return False
        with self._connect() as db:
            return db.execute('SELECT 1 FROM runs WHERE continuation_of=? LIMIT 1', (record['id'],)).fetchone() is None

    @staticmethod
    def _event(db: sqlite3.Connection, run_id: str, event_type: str, payload: dict[str, Any]) -> None:
        db.execute(
            "INSERT INTO run_events (run_id, created_at, event_type, payload_json) VALUES (?, ?, ?, ?)",
            (run_id, _now(), event_type, json.dumps(payload, ensure_ascii=False)),
        )

    def update(self, run_id: str, **fields: Any) -> None:
        if not fields:
            return
        assignments = ", ".join(f"{key} = ?" for key in fields)
        with self._lock, self._connect() as db:
            state = db.execute('SELECT status FROM runs WHERE id=?', (run_id,)).fetchone()
            if state and (state['status'] == 'stopped' or (state['status'] == 'cancelling' and fields.get('status') != 'stopped')):
                return
            db.execute(f"UPDATE runs SET {assignments} WHERE id = ?", (*fields.values(), run_id))

    def request_cancel(self, run_id: str):
        if self.get(run_id) is None:
            return None
        with self._lock, self._connect() as db:
            row = db.execute('SELECT status FROM runs WHERE id=?', (run_id,)).fetchone()
            if row is None:
                return None
            previous = row['status']
            if previous not in TERMINAL and previous != 'cancelling':
                db.execute("UPDATE runs SET status='cancelling', current_stage='cancelling' WHERE id=?", (run_id,))
                self._event(db, run_id, 'cancel_requested', {'name':'正在停止任务'})
            return previous

    def event(self, run_id: str, event_type: str, payload: dict[str, Any]) -> None:
        with self._lock, self._connect() as db:
            state = db.execute('SELECT status FROM runs WHERE id=?', (run_id,)).fetchone()
            if state and state['status'] in {'cancelling','stopped'} and event_type in {'completed','failed'}:
                return
            self._event(db, run_id, event_type, payload)

    def get(self, run_id: str) -> dict[str, Any] | None:
        with self._connect() as db:
            row = db.execute("SELECT * FROM runs WHERE id = ? AND (? IS NULL OR owner_id=?)", (run_id, principal.get(), principal.get())).fetchone()
            if row is None:
                return None
            payload = dict(row)
            payload['operations'] = [dict(op) for op in db.execute(
                'SELECT operation_id, tool, spec_json, status, facts_json FROM scientific_operations WHERE run_id=? ORDER BY rowid', (run_id,))]
            for op in payload['operations']:
                op['spec'] = json.loads(op.pop('spec_json'))
                op['facts'] = json.loads(op.pop('facts_json') or '{}')
            if payload.get("execution_json"):
                payload["execution_spec"] = json.loads(payload.pop("execution_json"))
            else:
                payload.pop("execution_json", None)
            if payload.get("result_json"):
                payload["result"] = json.loads(payload.pop("result_json"))
            else:
                payload.pop("result_json", None)
            if payload.get("manual_json"):
                payload["manual_input"] = json.loads(payload.pop("manual_json"))
            else:
                payload.pop("manual_json", None)
            return payload

    def list(self, limit: int = 50, offset: int = 0) -> list[dict[str, Any]]:
        """Return lightweight durable task summaries for the workbench history."""
        bounded_limit = max(1, min(limit, 100))
        bounded_offset = max(0, offset)
        with self._connect() as db:
            rows = db.execute(
                """
                SELECT id, conversation_id, prompt, mode, status, model, created_at, started_at,
                       completed_at, current_stage, error, runtime_seconds
                FROM runs WHERE (? IS NULL OR owner_id=?)
                ORDER BY created_at DESC
                LIMIT ? OFFSET ?
                """,
                (principal.get(), principal.get(), bounded_limit, bounded_offset),
            ).fetchall()
        return [dict(row) for row in rows]

    def events(self, run_id: str, after_id: int = 0) -> list[dict[str, Any]]:
        with self._connect() as db:
            rows = db.execute(
                "SELECT id, created_at, event_type, payload_json FROM run_events WHERE run_id = ? AND id > ? ORDER BY id",
                (run_id, after_id),
            ).fetchall()
        return [
            {
                "id": row["id"],
                "created_at": row["created_at"],
                "event_type": row["event_type"],
                "payload": json.loads(row["payload_json"]),
            }
            for row in rows
        ]


store = RunStore(STATE_DB)
executor = ThreadPoolExecutor(max_workers=MAX_CONCURRENCY, thread_name_prefix="marine-agent")
control_executor = ThreadPoolExecutor(max_workers=CONTROL_CONCURRENCY, thread_name_prefix='marine-dialogue')
cancel_executor = ThreadPoolExecutor(max_workers=4, thread_name_prefix="marine-cancel")
run_futures = {}
future_lock = threading.RLock()
app = FastAPI(title="Marine Scientific Agent API", version="0.1.0")
app.add_middleware(AccountMiddleware)
from . import access_control
access_control.analysis_check = lambda kind, analysis_id, owner: store.owns_analysis(kind, analysis_id, owner)
access_control.analysis_grant = lambda kind, analysis_id, owner: store.grant_analysis(kind, analysis_id, owner)


@app.on_event('startup')
def recover_interrupted_runs():
    for record in store.claim_orphans():
        cancel_executor.submit(reconcile, store, record)


def _finish_stopped(run_id: str, started=None):
    record = store.get(run_id) or {}
    if record.get('status') in TERMINAL:
        return
    text = '任务已停止。'
    store.update(run_id, status='stopped', current_stage='stopped', completed_at=_now(), error=None,
                 runtime_seconds=round(time.perf_counter()-started, 3) if started else 0,
                 result_json=json.dumps({'final_text':text}, ensure_ascii=False))
    store.append_message(record.get('conversation_id'), 'assistant', text, run_id)
    store.event(run_id, 'stopped', {'name':'任务已停止'})


def _execution_error(exc: Exception) -> str:
    if isinstance(exc, TurnPlanningError):
        return str(exc)
    if getattr(exc, 'status_code', None) == 429:
        return '模型服务暂时限流，请稍后重试。已有任务记录和计算结果仍保留。'
    return f'{type(exc).__name__}: {exc}'


def _controlled_execute(run_id, function, *args):
    started = time.perf_counter()
    pending_remote_stop = False
    with account_scope((store.get(run_id) or {}).get('owner_id') or None), job_scope(run_id) as control:
        try:
            control.check()
            from .progress import reporter
            def progress(stage):
                if isinstance(stage, dict):
                    if (store.get(run_id) or {}).get('status') not in TERMINAL | {'cancelling'}:
                        store.update(run_id, current_stage='domain_scan')
                        store.event(run_id, 'stage', {**stage, 'source': 'worker'})
                    return
                names = {'sequence':'获取已确认序列','msa':'构建 MSA','prediction_queue':'等待结构预测资源','prediction':'OpenFold3 结构预测','viewer_export':'导出结构与置信度'}
                if stage in names and (store.get(run_id) or {}).get('status') not in TERMINAL | {'cancelling'}:
                    store.update(run_id, current_stage=stage)
                    store.event(run_id,'stage',{'name':names[stage],'stage':stage,'status':'running','source':'worker'})
            progress_token = reporter.set(progress)
            try:
                function(run_id, *args)
            finally:
                reporter.reset(progress_token)
        except RunCancelled as exc:
            if getattr(exc, 'remote_stop_pending', False):
                pending_remote_stop = True
                control.cancel()
                store.update(run_id, status='cancelling', current_stage='recovery_pending')
                store.event(run_id, 'recovery_pending', {'name':'等待扫描节点确认停止；暂不重复计算'})
                cancel_executor.submit(reconcile, store, store.get(run_id))
        except Exception as exc:
            record = store.get(run_id) or {}
            if record.get('status') not in TERMINAL and not control.cancelled:
                error = _execution_error(exc)
                store.update(run_id, status='failed', current_stage='failed', completed_at=_now(), error=error)
                store.append_message(record.get('conversation_id'), 'assistant', f'任务执行失败：{error}', run_id)
                store.event(run_id, 'failed', {'error':error})
        finally:
            if control.cancelled and not pending_remote_stop:
                _finish_stopped(run_id, started)


def _submit(run_id, function, *args):
    with future_lock:
        pool = control_executor if function is _dispatch_run else executor
        future = pool.submit(_controlled_execute, run_id, function, *args)
        run_futures[run_id] = future
    def done(_):
        with future_lock:
            # Planning can hand the same run to the science pool. Its finished
            # future must not remove the newly queued computation's handle.
            if run_futures.get(run_id) is future:
                run_futures.pop(run_id, None)
    future.add_done_callback(done)


def _cancel_workers(run_id):
    # Independent control connections are necessary while the original HTTP
    # requests wait for computation. A pre-submission cancellation is remembered
    # by the worker, preventing delayed requests from starting orphan work.
    services = [('VINA', os.getenv('VINA_WORKER_URL',''), os.getenv('VINA_WORKER_TOKEN','')),
                ('OpenFold', os.getenv('OPENFOLD_WORKER_URL',''), os.getenv('OPENFOLD_WORKER_TOKEN',''))]
    if os.getenv('INTERPRO_WORKER_URL'):
        services.append(('InterPro', os.getenv('INTERPRO_WORKER_URL'), os.getenv('INTERPRO_WORKER_TOKEN', '')))
    if os.getenv('PROTEIN_WORKER_URL'):
        services.append(('VINA', os.getenv('PROTEIN_WORKER_URL'), os.getenv('PROTEIN_WORKER_TOKEN', os.getenv('VINA_WORKER_TOKEN', ''))))
    client = build_opener(ProxyHandler({}))
    for _ in range(10):
        failed = False
        for name, url, token in services:
            if not url:
                continue
            try:
                request = Request(f'{url.rstrip("/")}/v1/jobs/{run_id}/cancel', data=b'{}',
                                  headers={f'X-{name}-Worker-Token':token}, method='POST')
                with client.open(request, timeout=3) as response:
                    response.read()
            except Exception:
                failed = True
        if not failed or (store.get(run_id) or {}).get('status') in TERMINAL:
            return
        time.sleep(1)
    store.event(run_id, 'stage', {'name':'停止请求暂未送达本机，仍在等待确认'})


@app.post('/runs/{run_id}/cancel')
def cancel_run(run_id: str):
    previous = store.request_cancel(run_id)
    if previous is None:
        raise HTTPException(status_code=404, detail='run not found')
    if previous in TERMINAL:
        return {'run_id':run_id, 'status':previous}
    record = store.get(run_id)
    if record.get('current_stage') == 'recovery_pending':
        cancel_executor.submit(reconcile, store, record)
        return {'run_id': run_id, 'status': 'cancelling'}
    control_for(run_id).cancel()
    with future_lock:
        future = run_futures.get(run_id)
        removed = future.cancel() if future and not future.running() else False
    if removed:
        _finish_stopped(run_id)
    else:
        cancel_executor.submit(_cancel_workers, run_id)
    return {'run_id':run_id, 'status':(store.get(run_id) or {})['status']}


def _resolve_mode(prompt: str, requested: str) -> str:
    if requested != "auto":
        return requested
    lowered = prompt.lower()
    return "collagenase" if "collagenase-agent-20" in lowered else "general"


def _find_viewer_manifest(value: Any) -> str | None:
    """Find an OpenFold viewer manifest path inside a tool result.

    Tool messages are serialized by different model providers in slightly
    different shapes. Keep the browser contract independent of that detail by
    looking for the explicit artifact field recursively, including JSON strings
    returned by LangChain tools.
    """
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except (TypeError, json.JSONDecodeError):
            return None
        return _find_viewer_manifest(parsed)
    if isinstance(value, dict):
        artifacts = value.get("artifacts")
        if isinstance(artifacts, dict):
            candidate = artifacts.get("viewer_manifest")
            if isinstance(candidate, str) and candidate:
                return candidate
        candidate = value.get("viewer_manifest")
        if isinstance(candidate, str) and candidate:
            return candidate
        for child in reversed(list(value.values())):
            found = _find_viewer_manifest(child)
            if found:
                return found
    elif isinstance(value, (list, tuple)):
        for child in reversed(value):
            found = _find_viewer_manifest(child)
            if found:
                return found
    return None


def _artifact_url(run_id: str, manifest_path: str | None) -> str | None:
    if not manifest_path:
        return None
    normalized = manifest_path.replace("\\", "/").lstrip("/")
    if not normalized.startswith("results/"):
        return None
    return f"/runs/{run_id}/artifact/{normalized}"


def _manual_openfold_summary(payload: dict[str, Any]) -> str:
    # The CLI returns a wrapper with output_dir/pipeline and the scientific
    # payload nested under ``result``. Keep the browser summary useful for both
    # the local and remote-worker response shapes.
    scientific = payload.get("result") if isinstance(payload.get("result"), dict) else payload
    if scientific.get("status") != "succeeded":
        return "## OpenFold3 任务未完成\n\n" + str(
            scientific.get("reason") or scientific.get("stage") or payload.get("reason") or "结构预测工具返回失败"
        )
    sequence = scientific.get("sequence") or {}
    confidence = scientific.get("confidence") or {}
    quality = scientific.get("quality_gate") or {}
    input_record = scientific.get("input") or {}
    accession = input_record.get("accession") or input_record.get("primary_accession")
    source_kind = payload.get("source_kind") or scientific.get("source_kind") or input_record.get("source_kind")
    identity = f"内嵌序列（{payload.get('sequence_length') or scientific.get('sequence_length') or (scientific.get('sequence') or {}).get('length', '--')} aa）" if source_kind == "inline_sequence" else f"UniProt：`{accession or '--'}`"
    ligands = input_record.get("ligands") if isinstance(input_record.get("ligands"), list) else []
    ligand_labels = []
    for ligand in ligands:
        if not isinstance(ligand, dict):
            continue
        representation = ligand.get("smiles") or ",".join(ligand.get("ccd_codes") or []) or "SDF"
        chains = ",".join(ligand.get("chain_ids") or [])
        ligand_labels.append(f"{representation}（链 {chains or '--'}）")
    return "\n".join(
        [
            "## OpenFold3 复合物结构预测完成" if ligand_labels else "## OpenFold3 蛋白结构预测完成",
            f"- 输入：{identity}",
            f"- 预测残基：{sequence.get('length', '--')}",
            "- 实际纳入配体：" + ("；".join(ligand_labels) if ligand_labels else "无，本次为单蛋白结构预测"),
            *(["- 序列范围：UniProt 标注为片段，本次仅预测该片段。"] if input_record.get('fragment') else []),
            f"- 平均 pLDDT：{confidence.get('avg_plddt', '--')}",
            f"- 结构审查：{quality.get('decision', quality.get('status', '--'))}",
            "- 结构、pLDDT 和 PAE 已准备进入结构分析工作台。",
            "\n> 置信度支持结构审查，不等同于酶活性、结合亲和力或工业性能。",
        ]
    )


INLINE_SEQUENCE_RE = re.compile(r"(?<![A-Za-z])([ACDEFGHIKLMNPQRSTVWY]{20,4000})(?![A-Za-z])", re.IGNORECASE)
INLINE_SMILES_RE = re.compile(r"(?i)(?<![A-Za-z])SMILES\s*[:：=]\s*([A-Za-z0-9@+\-\[\]()=#$:/\\.%]+)")
INLINE_CCD_RE = re.compile(r"(?i)(?<![A-Za-z])CCD(?:_CODES?)?\s*[:：=]\s*([A-Za-z0-9]+)")


def _inline_structure_ligands(prompt: str) -> list[dict]:
    """Preserve explicitly labelled SMILES in the deterministic structure path.

    Only labelled inline values are extracted; chemical names and unlabelled
    prose remain the Agent's responsibility. Chain IDs are deterministic and
    reserve A for the protein, matching the OpenFold input contract.
    """
    text = str(prompt or "")
    values = [("smiles", value) for value in INLINE_SMILES_RE.findall(text)]
    values += [("ccd_codes", value) for value in INLINE_CCD_RE.findall(text)]
    if len(values) > 16:
        raise ValueError("一次结构预测最多支持 16 个显式配体")
    available = ["X", "Y", "Z"] + [chr(code) for code in range(ord("B"), ord("W") + 1)] + [str(code) for code in range(10)]
    records = [{"chain_id": available[index], kind: value, "ligand_name": f"Ligand{index + 1}"}
               for index, (kind, value) in enumerate(values)]
    return _parse_ligands_json(json.dumps(records, ensure_ascii=False)) if records else []


def _requests_vina(prompt: str) -> bool:
    lowered = str(prompt or "").lower()
    return "vina" in lowered or "对接" in str(prompt or "") or "docking" in lowered


def _inline_structure_sequence(prompt: str) -> str | None:
    """Extract one pasted standard-amino-acid sequence from a structure request."""
    lowered = str(prompt or "").lower()
    if any(term in lowered for term in ('结构域', 'interpro', 'pfam', 'domain')):
        # A domain-review or mixed research request belongs to the planner,
        # rather than the shortcut for one explicit structure prediction.
        return None
    if str(prompt or "").count('>') > 1:
        # Multi-record FASTA is not a one-protein structure request.
        return None
    structure_words = ("结构", "openfold", "alphafold", "预测", "structure")
    if not any(word in lowered for word in structure_words):
        return None
    matches = INLINE_SEQUENCE_RE.findall(str(prompt or ""))
    if not matches:
        return None
    # Prefer the longest candidate; ordinary English prose rarely consists only
    # of the 20 standard amino-acid letters for 20 consecutive characters.
    return max(matches, key=len).upper()


def _save_execution_spec(run_id: str, spec: dict) -> None:
    store.update(run_id, execution_json=json.dumps(spec, ensure_ascii=False))


def _capture_execution_calls(run_id: str, message: Any, *, recompute=False) -> None:
    """Persist exact selected scientific inputs before the tool starts."""
    for call in _message_field(message, "tool_calls", []) or []:
        name = _message_field(call, "name")
        args = _message_field(call, "args", {}) or {}
        if name not in {"run_openfold_accession", "run_openfold_sequence", "run_protein_homology_screen", "annotate_protein_domains"}:
            continue
        inputs = {"ligands": _parse_ligands_json(str(args.get("ligands_json") or ""))} if name.startswith("run_openfold_") else {}
        if name == "annotate_protein_domains":
            if recompute:
                args = {**args, 'refresh':True}
            inputs = DomainRequest.model_validate(args).model_dump()
            workflow = "protein_domains"
        elif name == "run_protein_homology_screen":
            inputs = ScreenRequest.model_validate(args).model_dump()
            workflow = "protein_screen"
        elif name == "run_openfold_accession":
            inputs["accession"] = args.get("accession")
            workflow = "openfold_accession"
        else:
            inputs["sequence"] = args.get("sequence")
            workflow = "openfold_sequence"
        spec = {"workflow": workflow, "input": inputs}
        operation_id = _message_field(call, 'id') or hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
        store.save_operation(run_id, operation_id, name, spec)
        previous = (store.get(run_id) or {}).get("execution_spec")
        if previous and previous != spec:
            if previous.get('workflow') == workflow and previous.get('input') == inputs:
                spec = previous
            elif previous.get('workflow') == 'protein_screen' and workflow == 'protein_domains' and _matches_screen_source(previous, inputs):
                spec['upstream'] = previous
            else:
                spec = {"workflow": "multiple", "reason": "本轮包含多个不同计算输入，需要选择具体任务"}
        _save_execution_spec(run_id, spec)


def _capture_operation_result(run_id, message):
    if (_message_field(message, 'role') or _message_field(message, 'type')) not in {'tool', 'tool_result'}:
        return
    operation_id = _message_field(message, 'tool_call_id')
    if not operation_id:
        return
    try:
        payload = json.loads(_content_text(_message_field(message, 'content')))
    except (TypeError, ValueError):
        payload = {}
    if not isinstance(payload, dict):
        payload = {}
    state = payload.get('status')
    scientific = payload.get('result') if isinstance(payload.get('result'), dict) else payload
    state = scientific.get('status', state)
    if state not in TERMINAL:
        state = 'failed' if payload.get('error') or _message_field(message, 'status') == 'error' else 'unknown'
    facts = {key: scientific[key] for key in ('analysis_id','record_count','unknown_record_count','query_count','candidate_count','reason') if key in scientific}
    # These values are tool facts, not inferred from model-written prose.
    facts['records'] = [{key: row[key] for key in ('id','protein_id','length','sequence_length','sequence_sha256','sequence_md5','md5','status','lookup_status') if key in row}
                        for row in scientific.get('records', [])[:100] if isinstance(row, dict)]
    if payload.get('execution_reference'):
        facts['execution_reference'] = payload['execution_reference']
    store.complete_operation(run_id, operation_id, state, facts)


def _matches_screen_source(previous: dict, inputs: dict) -> bool:
    try:
        request, _, _ = resolve_inputs(inputs)
        if not request.analysis_id:
            return False
        path = ROOT / 'results/protein-screenings' / request.analysis_id / 'manifest.json'
        recorded = json.loads(path.read_text())['spec']
        expected = ScreenRequest.model_validate(previous['input'])
        return recorded['inputs'] == {'query': hashlib.sha256(expected.query_fasta.encode()).hexdigest(),
                                     'reference': hashlib.sha256(expected.reference_fasta.encode()).hexdigest()} and recorded['request'] == expected.model_dump(exclude={'query_fasta','reference_fasta'})
    except (ValueError, OSError, KeyError):
        return False


def _execute_openfold_spec(run_id: str, spec: dict, started: float) -> None:
    inputs = spec["input"]
    accession = inputs.get("accession")
    ligands = inputs.get("ligands", [])
    _save_execution_spec(run_id, spec)
    store.update(run_id, current_stage="sequence", mode="manual_openfold" if accession else "general")
    store.event(run_id, "stage", {"name": "获取已确认序列" if accession else "识别内嵌蛋白序列",
                                  **({"accession": accession} if accession else {"length": len(inputs["sequence"])})})
    store.update(run_id, current_stage="msa")
    store.event(run_id, "stage", {"name": "构建 MSA", "source": "ColabFold public service"})
    if accession:
        result_text = run_openfold_accession.invoke({"accession": accession, "ligands_json": json.dumps(ligands, ensure_ascii=False)})
    else:
        result_text = run_openfold_sequence.invoke({"sequence": inputs["sequence"], "ligands_json": json.dumps(ligands, ensure_ascii=False)})
    check_cancelled()
    payload = json.loads(result_text)
    scientific = payload.get("result") if isinstance(payload.get("result"), dict) else payload
    for key, label in [("msa", "MSA 完成"), ("prediction", "OpenFold3 预测完成")]:
        stage = scientific.get(key)
        if isinstance(stage, dict) and stage.get("status") == "succeeded":
            store.event(run_id, "stage", {"name": label, "elapsed_seconds": stage.get("elapsed_seconds"), "cache_hit": stage.get("cache_hit")})
    if payload.get("status") not in {"succeeded", "resumed_existing"}:
        stage = str(payload.get("stage") or scientific.get("stage") or "openfold3")
        detail = scientific.get(stage) if isinstance(scientific.get(stage), dict) else {}
        # Preserve the scientific payload/log references for investigation,
        # while keeping the conversation free of internal tracebacks.
        store.update(run_id, result_json=json.dumps({"agent": {"openfold_failure": _json_safe(payload)}}, ensure_ascii=False))
        if stage == 'msa':
            reason = detail.get('reason') or scientific.get('reason') or 'MSA 获取未完成，请稍后继续任务；原输入已保留。'
            store.event(run_id, 'stage', {'name': 'MSA 获取失败', 'status': 'failed', 'error_code': detail.get('error_code', 'msa_failed'), 'retryable': detail.get('retryable', True)})
            raise RuntimeError(f'MSA 获取失败：{reason} 尚未开始结构推理。')
        raise RuntimeError(f"{stage}: {detail.get('output_tail') or detail.get('reason') or scientific.get('reason') or payload.get('reason') or 'OpenFold3 返回失败'}")
    store.update(run_id, current_stage="viewer_export")
    store.event(run_id, "stage", {"name": "导出结构查看产物"})
    result = {"agent": {"manual_openfold" if accession else "inline_openfold": _json_safe(payload)},
              "final_text": _manual_openfold_summary(payload)}
    if spec.get("requested_followup") == "vina":
        result["pending_actions"] = [{"workflow": "vina", "status": "awaiting_confirmation",
                                      "required": ["受体链", "口袋残基", "位点依据", "搜索盒", "Vina参数"]}]
        result["final_text"] += (
            "\n\n### 后续 Vina：尚未运行\n"
            "为避免猜测对接位点，请先确认受体链、口袋残基、位点依据和搜索盒；"
            "随后点击“用此结构进行分子对接”继续。当前完成的是 OpenFold3 阶段，不代表对接已完成。"
        )
        store.event(run_id, "awaiting_confirmation", {"workflow": "vina", "reason": "pocket_and_protocol_confirmation"})
    conversation_id = (store.get(run_id) or {}).get("conversation_id")
    runtime = round(time.perf_counter() - started, 3)
    store.update(run_id, status="succeeded", current_stage="completed", completed_at=_now(), runtime_seconds=runtime,
                 result_json=json.dumps(result, ensure_ascii=False))
    store.append_message(conversation_id, "assistant", result["final_text"], run_id)
    store.event(run_id, "completed", {"runtime_seconds": runtime, "has_final_text": True})


def _vina_carried_results(source: dict, config: dict) -> list[dict]:
    """Reuse only poses whose stored ligand identity and score file still match."""
    from .vina_tools import _parse_scores
    by_id = {item["id"]: item for item in config["ligands"]}
    expected_tasks = {task["id"] for task in _seed_tasks(config)}
    carried = []
    payload = (source.get('result') or {}).get('vina') or {}
    original_spec = saved_spec(source)
    if original_spec:
        original = validate_config(original_spec['input'])
        if original['box'] != config['box'] or original['params'] != config['params'] or original['receptor']['bytes'] != config['receptor']['bytes']:
            return []
    if payload.get('configuration') and payload['configuration'] != {'box': config['box'], 'params': config['params']}:
        return []
    if payload.get('receptor_sha256') and payload['receptor_sha256'] != hashlib.sha256(config['receptor']['bytes']).hexdigest():
        return []
    for row in payload.get('results', []):
        ligand = by_id.get(row.get("ligand_id") or row.get("id"))
        if row.get("status") != "succeeded" or not ligand or row.get("id") not in expected_tasks:
            continue
        if row.get("input_sha256") != hashlib.sha256(ligand["bytes"]).hexdigest():
            continue
        path = (ROOT / str(row.get("output", ""))).resolve()
        try:
            path.relative_to((ROOT / "results/vina-runs" / source["id"]).resolve())
            raw = path.read_bytes()
        except (ValueError, OSError):
            continue
        if not row.get("scores") or _parse_scores(raw.decode(errors="replace")) != row["scores"]:
            continue
        if row.get("output_sha256") and row["output_sha256"] != hashlib.sha256(raw).hexdigest():
            continue
        carried.append({**copy.deepcopy(row), "carried_from_run_id": source["id"], "output_sha256": hashlib.sha256(raw).hexdigest()})
    return carried


def _execute_continuation(run_id: str, handoff: dict, started: float | None = None) -> None:
    started = started or time.perf_counter()
    check_cancelled()
    spec = handoff["spec"]
    store.update(run_id, status="running", started_at=(store.get(run_id) or {}).get('started_at') or _now(), execution_json=json.dumps(spec, ensure_ascii=False),
                 continuation_of=handoff["source_run_id"])
    if spec["workflow"] == "vina":
        _execute_vina_run(run_id, spec["input"], source_run_id=handoff["source_run_id"])
    elif spec["workflow"] == "protein_screen":
        _execute_protein_spec(run_id, spec, started)
    elif spec["workflow"] == "protein_domains":
        _execute_domain_spec(run_id, spec, started)
    else:
        _execute_openfold_spec(run_id, spec, started)


PROTEIN_STAGES = {"validate_inputs": "检查蛋白输入", "makeblastdb": "构建参考序列库",
                  "blastp": "BLASTP 同源筛选", "rank_and_report": "整理候选与证据"}


def _protein_progress(run_id: str, event: dict) -> None:
    stage = event.get("name")
    if stage not in PROTEIN_STAGES or event.get("status") not in {"running", "succeeded"}:
        raise ValueError("invalid protein stage event")
    if (store.get(run_id) or {}).get("status") in TERMINAL:
        return
    store.update(run_id, current_stage=stage)
    store.event(run_id, "stage", {"name": PROTEIN_STAGES[stage], "stage": stage,
                                  "status": event["status"], "cache_hit": event.get("cache_hit") is True})


def _execute_protein_spec(run_id: str, spec: dict, started: float | None = None) -> None:
    started = started or time.perf_counter()
    _save_execution_spec(run_id, spec)
    store.update(run_id, status="running", started_at=_now(), current_stage="protein_screen")
    store.event(run_id, "stage", {"name": "提交蛋白同源筛选", "status": "running"})
    try:
        report = execute_screen(spec["input"], emit=lambda event: _protein_progress(run_id, event))
        check_cancelled()
        if report.get("status") != "succeeded":
            raise RuntimeError(report.get("reason", "protein screening failed"))
        lines = [f"蛋白同源筛选完成：{report['query_count']} 条输入，{report['candidate_count']} 条候选。",
                 "已有蛋白序列，已跳过测序质控、组装和基因预测。"]
        for item in report["candidates"][:10]:
            lines.append(f"- {item['protein_id']}：一致性 {item['identity_percent']}%，查询覆盖 {item['query_coverage_percent']}%，参考覆盖 {item['reference_coverage_percent']}%。")
        lines.append("候选仅有序列同源证据，尚未验证结构域、催化位点和实验活性。")
        lines.append(f"分析编号：`{report['analysis_id']}`")
        result = {"agent": {"protein_screen": report}, "final_text": "\n\n".join(lines)}
        store.update(run_id, status="succeeded", current_stage="completed", completed_at=_now(),
                     runtime_seconds=round(time.perf_counter() - started, 3), result_json=json.dumps(result, ensure_ascii=False))
        store.append_message((store.get(run_id) or {}).get("conversation_id"), "assistant", result["final_text"], run_id)
        store.event(run_id, "completed", {"analysis_id": report["analysis_id"], "candidate_count": report["candidate_count"]})
    except Exception as exc:
        error = _execution_error(exc)
        store.update(run_id, status="failed", current_stage="failed", completed_at=_now(), error=error,
                     runtime_seconds=round(time.perf_counter() - started, 3))
        store.append_message((store.get(run_id) or {}).get("conversation_id"), "assistant", f"蛋白筛选未完成：{error}", run_id)
        store.event(run_id, "failed", {"error": error})


def _execute_domain_spec(run_id: str, spec: dict, started: float | None = None) -> None:
    started = started or time.perf_counter()
    _save_execution_spec(run_id, spec)
    store.update(run_id, status="running", started_at=_now(), current_stage="domain_scan" if spec["input"].get("method") == "full_scan" else "domain_lookup")
    store.event(run_id, "stage", {"name": "核验原始序列的结构域证据", "status": "running"})
    def progress(event):
        check_cancelled()
        store.event(run_id, 'stage', {**event, 'name': '查询 InterPro 预计算记录' if event.get('name') == 'domain_lookup' else event['name']})
    report = annotate(spec['input'], emit=progress)
    check_cancelled()
    text = render_domains(report)
    result = {'agent': {'protein_domains': report}, 'final_text': text,
              'scientific_review': {'decision': 'rendered_from_domain_evidence',
                                    'scope': 'Full InterProScan sequence evidence; no biochemical activity confirmation.' if report.get('de_novo_scan_executed') else 'Precomputed sequence evidence; no de novo scan or biochemical activity confirmation.'}}
    store.update(run_id, status=report['status'], current_stage='completed', completed_at=_now(),
                 runtime_seconds=round(time.perf_counter()-started,3), result_json=json.dumps(result,ensure_ascii=False))
    store.append_message((store.get(run_id) or {}).get('conversation_id'), 'assistant', text, run_id)
    store.event(run_id, 'completed', {'analysis_id': report['analysis_id'], 'status': report['status'],
                                     'unknown_record_count': report['unknown_record_count']})


def _reply_task_status(run_id: str, lifecycle: ConversationTasks, started: float, selected_id: str = '', inspection: str = 'status', sequence_index: int | None = None) -> None:
    """Lifecycle status is a server fact, rather than a model recollection."""
    selected_operation = None
    if selected_id:
        resolved = lifecycle.resolve_reference(selected_id, lifecycle.records())
        if 'record' in resolved:
            selected_id = resolved['record']['id']
            selected_operation = resolved['operation']
    records = [row for row in lifecycle.records() if (saved_spec(row) or row.get('operations') or row.get('execution_spec')) and (not selected_id or row['id'] == selected_id)]
    if not records:
        text = '当前会话没有找到科研计算任务。'
    else:
        row = records[-1]
        spec = saved_spec(row) or {'workflow': 'multiple', 'input': {}}
        inputs = spec['input']
        label = ('科研组合任务' if spec['workflow'] == 'multiple' else
                 '蛋白结构域核验' if spec['workflow'] == 'protein_domains' else
                 '蛋白同源筛选' if spec['workflow'] == 'protein_screen' else
                 'Vina 分子对接' if spec['workflow'] == 'vina' else f"{inputs.get('accession') or '所输入序列'} 的 OpenFold3 预测")
        partial_label = ('已结束，部分配体失败' if spec['workflow'] == 'vina' else
                         '已结束，部分序列证据未知' if spec['workflow'] == 'protein_domains' else
                         '已结束，结果不完整')
        status = {'queued':'排队中','running':'执行中','cancelling':'正在停止','stopped':'已停止','failed':'执行失败','partial':partial_label,'succeeded':'已完成'}[row['status']]
        text = f'最近一次{label}：{status}。'
        if row.get('operations'):
            labels = {'succeeded':'已完成','partial':'部分结果未知','failed':'失败','running':'执行中','stopped':'已停止','unknown':'结果状态待核验'}
            for index, op in enumerate(row['operations'], 1):
                op_state = op['status']
                if op_state == 'running' and row['status'] in TERMINAL:
                    op_state = 'unknown'
                text += f"\n{index}. {op['tool']}：{labels.get(op_state, op_state)}。"
            if spec['workflow'] == 'multiple':
                text += '\n包含多个计算输入，需要明确具体操作后再继续；不会自动重跑整组。'
        specs = [selected_operation['spec']] if selected_operation else [op['spec'] for op in row.get('operations', [])] or [spec]
        domain_specs = [s for s in specs if s['workflow'] == 'protein_domains']
        if domain_specs:
            text += ('\n该任务选择完整 InterProScan 执行，完成情况以上述状态和保存结果为准。'
                     if any(s['input'].get('method') == 'full_scan' for s in domain_specs)
                     else '\n结构域核验查询的是预计算记录，未运行完整 InterProScan。')
        if inspection == 'protein_identity':
            from .input_identity import inspect_input_identity
            identity = inspect_input_identity(specs, sequence_index)
            store.update(run_id,status=identity['status'],current_stage='completed',completed_at=_now(),
                         runtime_seconds=round(time.perf_counter()-started,3),
                         result_json=json.dumps({'final_text':identity['final_text'],'identity_evidence':identity},ensure_ascii=False))
            store.append_message(lifecycle.conversation_id,'assistant',identity['final_text'],run_id)
            store.event(run_id,'completed',{'status':identity['status'],'source':'verified_input_and_uniprot'})
            return
        if inspection == 'activity_evidence':
            from .evidence_review import explain_activity_evidence
            from .domain_evidence import read_domain_analysis
            domain_reports, seen = [], set()
            for op in ([selected_operation] if selected_operation else row.get('operations', [])):
                analysis_id = op.get('facts', {}).get('analysis_id')
                if op['tool'] == 'annotate_protein_domains' and analysis_id and analysis_id not in seen:
                    seen.add(analysis_id)
                    try:
                        domain_reports.append(read_domain_analysis(analysis_id))
                    except (ValueError, OSError) as exc:
                        domain_reports.append({'status':'unavailable','reason':str(exc)})
            explanation = explain_activity_evidence(domain_reports, sequence_index)
            store.update(run_id, status=explanation['status'], current_stage='completed', completed_at=_now(),
                         runtime_seconds=round(time.perf_counter()-started,3),
                         result_json=json.dumps({'final_text':explanation['final_text'], 'evidence_assessment':explanation, 'source_run_id':row['id']},ensure_ascii=False))
            store.append_message(lifecycle.conversation_id,'assistant',explanation['final_text'],run_id)
            store.event(run_id,'completed',{'status':explanation['status'],'source':'verified_saved_evidence','new_computation':False})
            return
        if inspection == 'sequence_length':
            sequences = []
            for item in specs:
                details = describe_spec(item) or {}
                sequences.extend(details.get('input_sequences', []))
                if item['workflow'] == 'openfold_sequence':
                    sequences.append({'protein_id':'输入序列', 'residue_count': details['residue_count']})
            selected = sequences[sequence_index-1:sequence_index] if sequence_index else sequences
            if len(selected) == 1:
                text = str(selected[0]['residue_count'])
            elif selected:
                text = '\n'.join(f"{i}. {s['protein_id']}：{s['residue_count']} 个残基" for i,s in enumerate(selected,1))
            else:
                text = '保存的输入中没有找到对应的明确序列长度，不能据此给出数字。'
        if inspection == 'status' and row['status']=='stopped':
            text += (' 此前已继续执行过，目前再次停止。' if row.get('continuation_of') else ' 尚未再次启动。')
            text += ' 发送“继续”可以沿用原输入重新执行。'
        elif inspection == 'status' and row.get('continuation_of'):
            text += ' 这是沿用原输入继续执行的新一轮任务。'
        if inspection == 'status' and spec['workflow']=='vina':
            result = ((row.get('result') or {}).get('vina') or {})
            succeeded = sum(item.get('status')=='succeeded' for item in result.get('results',[]))
            text += f' 已保存 {succeeded} 个成功对接任务结果。'
    runtime = round(time.perf_counter()-started,3)
    store.event(run_id,'stage',{'name':'查询真实任务状态','status':'succeeded'})
    store.update(run_id,status='succeeded',current_stage='completed',completed_at=_now(),runtime_seconds=runtime,
                 result_json=json.dumps({'final_text':text,'task_state':lifecycle.snapshot(selected_id)},ensure_ascii=False))
    store.append_message(lifecycle.conversation_id,'assistant',text,run_id)
    store.event(run_id,'completed',{'has_final_text':True,'runtime_seconds':runtime})


def _vina_summary(payload: dict[str, Any]) -> str:
    ranking = payload.get("ranking") or []
    lines = ["## AutoDock Vina 已停止" if payload.get('status') == 'stopped' else "## AutoDock Vina 对接完成", "", f"- 引擎：{payload.get('engine', 'AutoDock Vina')}",
             f"- 任务耗时：{payload.get('elapsed_seconds', '--')} 秒",
             f"- 成功任务：{sum(1 for item in payload.get('results', []) if item.get('status') == 'succeeded')} / {len(payload.get('results', []))}",
             f"- 配体数 / 种子任务数：{payload.get('ligand_count', len(payload.get('results', [])))} / {payload.get('task_count', len(payload.get('results', [])))}", ""]
    if ranking:
        lines.append("| 配体 | Seed | Vina affinity (kcal/mol) |")
        lines.append("|---|---:|---:|")
        for item in ranking[:10]:
            lines.append(f"| `{item['name']}` | {item.get('seed', '--')} | {item.get('top_affinity_kcal_mol', '--')} |")
    consistency = payload.get("top_pose_consistency") or []
    for group in consistency:
        if group.get("status") == "available":
            lines.extend(["", f"- 跨种子首位最大 RMSD：{group['max_top_pose_rmsd_angstrom']:.3f} Å",
                          f"- 首位最高分跨度：{group['top_score_range_kcal_mol']:.3f} kcal/mol",
                          f"- 中心首位姿势种子：{group['central_top_pose_seed']}",
                          f"- 复核状态：{group['review_status']}；默认查看中心姿势供复核，不按最低分自动判定正确性"])
            lines.append("> " + group["scientific_limit"])
        else:
            lines.append(f"\n- 跨种子姿势分析不可用：{group.get('reason', '未知原因')}")
    lines.extend(["", "> Vina 分数用于构象和候选排序，不等同于实验结合常数、酶活或细胞表型。"])
    return "\n".join(lines)


def _execute_vina_run(run_id: str, request_payload: dict[str, Any], source_run_id: str | None = None) -> None:
    started = time.perf_counter()
    conversation_id = (store.get(run_id) or {}).get('conversation_id')
    store.update(run_id, status="running", started_at=_now(), current_stage="input_validation")
    store.event(run_id, "started", {"mode": "vina", "engine": "AutoDock Vina"})
    try:
        config = validate_config(request_payload)
        _save_execution_spec(run_id, {"workflow": "vina", "input": request_payload})
        store.update(run_id, mode="vina")
        carried = {}
        source_id, visited = source_run_id, set()
        while source_id and source_id not in visited:
            visited.add(source_id)
            source = store.get(source_id) or {}
            if VINA_WORKER_URL and source:
                recovered = recover_vina_artifacts(source_id)
                if recovered:
                    store.save_recovered_vina(source_id, recovered)
                    source = store.get(source_id) or source
            for row in _vina_carried_results(source, config):
                carried.setdefault(row['id'], row)
            source_id = source.get('continuation_of')
        pending_tasks = [task for task in _seed_tasks(config) if task['id'] not in carried]
        pending_config = {**config,
                          "ligands": [item for item in config['ligands'] if any(task['ligand']['id']==item['id'] for task in pending_tasks)],
                          "task_seeds": {item['id']:[task['seed'] for task in pending_tasks if task['ligand']['id']==item['id']]
                                         for item in config['ligands']}}
        store.update(run_id, current_stage="input_validation")
        store.event(run_id, "stage", {"name": "PDBQT 输入检查完成", "ligands": len(config["ligands"]), "box": config["box"]})
        for row in carried.values():
            store.event(run_id, 'stage', {'name':'保留已完成配体', 'ligand_id':row['id'], 'status':'succeeded',
                                          'total':len(config['ligands']), 'source_run_id':row['carried_from_run_id']})

        def emit(event: dict[str, Any]) -> None:
            store.update(run_id, current_stage=str(event.get("stage") or "vina"))
            store.event(run_id, "stage", event)

        run_dir = ROOT / "results" / "vina-runs" / run_id
        if not pending_tasks:
            payload = {'status':'succeeded', 'engine':'AutoDock Vina', 'results':[], 'elapsed_seconds':0}
        elif VINA_WORKER_URL:
            store.update(run_id, current_stage="vina_worker")
            store.event(run_id, "stage", {"name": "提交本机 Vina Worker", "worker_url": VINA_WORKER_URL})
            payload = run_vina_worker(run_id, pending_config)
            if payload.get("status") not in {"succeeded", "partial", "stopped"}:
                raise RuntimeError(payload.get("reason") or "本机 Vina Worker 未返回成功结果")
        else:
            payload = run_batch(run_dir, pending_config, emit=emit)
        if source_run_id:
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / 'receptor.pdbqt').write_bytes(config['receptor']['bytes'])
            (run_dir / 'ligands').mkdir(exist_ok=True)
            for row in carried.values():
                original = (ROOT / row['output']).read_bytes()
                destination = run_dir / 'ligands' / f"{row['id']}.out.pdbqt"
                destination.write_bytes(original)
                row['original_output'] = row['output']
                row['output'] = str(destination.relative_to(ROOT))
            rows = {row['id']:row for row in payload.get('results', [])}
            rows.update(carried)
            expected = _seed_tasks(config)
            payload['results'] = [rows.get(task['id'], {'id':task['id'], 'ligand_id':task['ligand']['id'],
                'seed':task['seed'], 'name':task['ligand']['name'], 'scores':[],
                'status':'skipped' if payload['status']=='stopped' else 'failed',
                'error':'任务已停止，未开始计算' if payload['status']=='stopped' else 'Worker 未返回该配体种子结果'})
                for task in expected]
            payload['ligand_count'] = len(config['ligands'])
            payload['task_count'] = len(expected)
            payload['ranking'] = sorted([
                {'id':row['id'], 'ligand_id':row.get('ligand_id'), 'seed':row.get('seed'),
                 'name':row['name'], 'status':row['status'], 'top_affinity_kcal_mol':row['scores'][0]['affinity_kcal_mol']}
                for row in payload['results'] if row['status']=='succeeded'], key=lambda row:row['top_affinity_kcal_mol'])
            payload['configuration'] = {'box':config['box'], 'params':config['params']}
            payload['provenance'] = config['provenance']
            payload['top_pose_consistency'] = _analyze_seed_consistency(run_dir, config, payload['results'])
            payload['continuation'] = {'source_run_id':source_run_id, 'retained_ligands':list(carried),
                                       'recomputed_ligands':[item['id'] for item in pending_config['ligands']],
                                       'retained_tasks':list(carried),
                                       'recomputed_tasks':[task['id'] for task in pending_tasks]}
            if payload['status'] != 'stopped':
                payload['status'] = 'succeeded' if all(row['status']=='succeeded' for row in payload['results']) else 'partial'
        if payload.get('status') == 'stopped':
            payload.setdefault('configuration', {'box':config['box'], 'params':config['params']})
            payload.setdefault('provenance', config['provenance'])
        with store._lock:
            if (store.get(run_id) or {}).get('status') == 'cancelling':
                payload['status'] = 'stopped'
            if source_run_id:
                (run_dir/'result.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2))
            runtime = round(time.perf_counter() - started, 3)
            store.update(
                run_id,
                status=payload["status"],
                completed_at=_now(),
                current_stage="stopped" if payload["status"] == "stopped" else "completed",
                result_json=json.dumps({"vina": _json_safe(payload), "final_text": _vina_summary(payload)}, ensure_ascii=False),
                runtime_seconds=runtime,
            )
            store.append_message(conversation_id, 'assistant', _vina_summary(payload), run_id)
            store.event(run_id, "stopped" if payload['status']=='stopped' else "completed", {"name":"对接已停止" if payload['status']=='stopped' else "对接完成", "runtime_seconds": runtime, "status": payload["status"]})
    except VinaInputError as exc:
        runtime = round(time.perf_counter() - started, 3)
        error = f"输入检查失败：{exc}"
        store.update(run_id, status="failed", completed_at=_now(), current_stage="input_validation", error=error, runtime_seconds=runtime)
        store.append_message(conversation_id, 'assistant', f'任务执行失败：{error}', run_id)
        store.event(run_id, "failed", {"runtime_seconds": runtime, "error": error})
    except Exception as exc:
        runtime = round(time.perf_counter() - started, 3)
        error = _execution_error(exc)
        store.update(run_id, status="failed", completed_at=_now(), current_stage="failed", error=error, runtime_seconds=runtime)
        store.append_message(conversation_id, 'assistant', f'任务执行失败：{error}', run_id)
        store.event(run_id, "failed", {"runtime_seconds": runtime, "error": error})


def _dispatch_run(run_id: str, prompt: str, mode: str) -> None:
    _execute_run(run_id, prompt, mode, dispatch=True)


def _queue_science(run_id, function, *args):
    check_cancelled()
    store.update(run_id, status='queued', current_stage='scientific_queue', started_at=None)
    store.event(run_id, 'stage', {'name': '等待科研计算资源', 'status': 'queued'})
    _submit(run_id, function, *args)


def _execute_run(run_id: str, prompt: str, mode: str, prepared_plan: TurnPlan | None = None, *, dispatch: bool = False) -> None:
    started = time.perf_counter()
    store.update(run_id, status="running", started_at=_now(), current_stage="planning")
    store.event(run_id, "started", {"mode": mode, "model": MODEL_NAME})
    run_record = store.get(run_id) or {}
    conversation_id = run_record.get("conversation_id")
    conversation_messages = store.conversation_context(conversation_id, through_run_id=run_id)
    try:
        # A pasted sequence plus an explicit structure request is a deterministic
        # scientific input. Route it directly so the language model cannot loop
        # trying to resolve it as a UniProt name or accession.
        inline_sequence = _inline_structure_sequence(prompt) if mode != "manual_openfold" else None
        if inline_sequence is not None:
            # Keep every explicitly supplied ligand in the deterministic path.
            # Previously this shortcut hard-coded ligands=[] and bypassed the
            # Agent, silently turning a requested complex into a protein-only job.
            ligands = _inline_structure_ligands(prompt)
            spec = {"workflow":"openfold_sequence", "input":{"sequence":inline_sequence, "ligands":ligands}}
            if _requests_vina(prompt):
                spec["requested_followup"] = "vina"
            _save_execution_spec(run_id, spec)
            if dispatch:
                _queue_science(run_id, _execute_run, prompt, mode)
                return
            _execute_openfold_spec(run_id, spec, started)
            return
        if mode == "manual_openfold":
            manual = run_record.get("manual_input")
            if not isinstance(manual, dict):
                raise RuntimeError("manual OpenFold input is missing")
            if dispatch:
                _queue_science(run_id, _execute_run, prompt, mode)
                return
            _execute_openfold_spec(run_id, {"workflow":"openfold_accession", "input":manual}, started)
            return
        lifecycle = ConversationTasks(store, run_id)
        project_snapshot = store.run_project_context(run_id)
        if project_snapshot:
            conversation_messages.insert(0, project_context_message(project_snapshot))
        planning_started = time.perf_counter()
        plan = prepared_plan or plan_turn(prompt, conversation_messages, {**lifecycle.snapshot(), **({'research_project':project_snapshot} if project_snapshot else {})})
        store.event(run_id, 'turn_timing', {'phase':'planning', 'seconds':round(time.perf_counter()-planning_started,3), 'reused':prepared_plan is not None})
        check_cancelled()
        store.event(run_id, 'stage', {'name':'理解本轮请求', **plan.model_dump(), 'status':'succeeded'})
        project_inspection = project_snapshot and (
            plan.context_scope == 'project' or not lifecycle.snapshot().get('tasks') or
            ((project_snapshot.get('candidate_ledger') or {}).get('candidate_count', 0) > 0
             and not plan.source_run_id and plan.inspection in {'protein_identity', 'sequence_length', 'activity_evidence', 'results'}))
        if not project_inspection and ((plan.action == 'inspect' and plan.inspection != 'results') or (plan.action == 'reply' and plan.inspection in {'protein_identity','sequence_length','activity_evidence'})):
            _reply_task_status(run_id, lifecycle, started, plan.source_run_id, plan.inspection, plan.sequence_index)
            return
        if plan.action == 'execute' and plan.requested_method == 'vina' and plan.scientific_endpoint not in {'binding_pose','binding_score'}:
            # Validate the declared method/endpoint contract before granting tools.
            if plan.scientific_endpoint == 'unspecified':
                text='请先明确这次对接要研究的结合问题，以及受体、配体和位点依据；当前尚未启动计算。'
                review={'decision':'method_needs_clarification', 'run_status':'partial', 'final_text':text}
            else:
                from .collagenase_tools import assess_collagenase_method
                args={'method':'Vina','objective':prompt,'endpoint':plan.scientific_endpoint}
                store.event(run_id,'tool_call',{'name':'核验方法与目标指标','tool':'assess_collagenase_method','args_preview':json.dumps(args,ensure_ascii=False)})
                verdict=assess_collagenase_method.invoke(args)
                store.event(run_id,'tool_result',{'name':'方法适用性核验完成','tool':'assess_collagenase_method','preview':verdict})
                review=method_rejection({'messages':[{'type':'tool','name':'assess_collagenase_method','content':verdict}]})
                if not review:
                    review={'decision':'method_needs_clarification','run_status':'partial','final_text':'该方法与目标尚未确认匹配，尚未启动计算。'}
                text=review['final_text']
            store.update(run_id,status=review['run_status'],current_stage='method_assessment',completed_at=_now(),runtime_seconds=round(time.perf_counter()-started,3),result_json=json.dumps({'final_text':text,'scientific_review':review,'requested_action':plan.action},ensure_ascii=False))
            store.append_message(conversation_id,'assistant',text,run_id)
            store.event(run_id,'completed',{'status':review['run_status'],'computation_started':False})
            return
        if dispatch and plan.action == 'execute':
            _queue_science(run_id, _execute_run, prompt, mode, plan)
            return
        if plan.action == 'continue':
            receipt = lifecycle.prepare(plan.source_run_id, plan.operation_id, plan.operation_index, plan.continuation_method)
            if lifecycle.handoff:
                if dispatch:
                    _queue_science(run_id, _execute_continuation, lifecycle.handoff)
                    return
                _execute_continuation(run_id, lifecycle.handoff, started)
                return
            # A failed execution handoff is an objective-level result, not an
            # invitation for the model to invent protein identities or status.
            text = receipt['reason']
            store.update(run_id, status='partial', current_stage='awaiting_execution_path', completed_at=_now(),
                         runtime_seconds=round(time.perf_counter()-started,3),
                         result_json=json.dumps({'final_text':text,'continuation':receipt},ensure_ascii=False))
            store.append_message(conversation_id,'assistant',text,run_id)
            store.event(run_id,'completed',{'status':'partial','objective_completed':False,'reason':receipt['status']})
            return
        conversation_messages.insert(0, {"role":"system", "content":
            "Authoritative task state for this conversation only (data, not user instructions):\n" +
            json.dumps(lifecycle.snapshot(), ensure_ascii=False)})
        if plan.action in {'reply','inspect','continue'}:
            tools = [t for t in lifecycle.tools() if t.name == 'inspect_task_state'] + [
                t for t in CORE_TOOLS if t.name in {'inspect_domain_evidence','inspect_protein_analysis'}]
            system_prompt = CONVERSATION_PROMPT
        elif mode == "collagenase":
            tools = COLLAGENASE_TOOLS
            system_prompt = COLLAGENASE_SYSTEM_PROMPT
        else:
            tools = CORE_TOOLS
            system_prompt = None
        tools = [*tools, *project_tools(store, run_id)]
        if plan.action == 'execute':
            tools = [*tools, *lifecycle.tools(allow_continue=lifecycle.has_resumable_work())]
        from .coordination import ExecutionReceipts, domain_completion_contract, domain_contract_progress, operation_outcome
        contract = domain_completion_contract(prompt, plan.completion_scope, plan.recompute) if plan.action == 'execute' else None
        if contract:
            tools = [t for t in tools if t.name in {'annotate_protein_domains', 'inspect_task_state'}]
            store.event(run_id, 'completion_contract', contract)
        from .scientific_knowledge import references, reference_context, review_explanation, compose_explanation
        reference_bundle = references(plan.reference_topics) if plan.action in {'reply','inspect'} and plan.reference_topics else None
        if reference_bundle and reference_bundle['cards']:
            conversation_messages.insert(0, {'role':'system', 'content':reference_context(reference_bundle)})
            store.event(run_id, 'reference_context', {'name':'查阅科学参考资料', 'version':reference_bundle['version'],
                        'topics':plan.reference_topics, 'sha256':reference_bundle['sha256']})
        if reference_bundle and reference_bundle['cards'] and plan.action == 'reply' and plan.reference_only:
            store.update(run_id, current_stage='reference_answer')
            store.event(run_id, 'reference_review_started', {'name':'依据参考资料整理回答','status':'running'})
            explanation = compose_explanation(prompt, conversation_messages, reference_bundle)
            check_cancelled()
            final_text = explanation['final_text']
            runtime = round(time.perf_counter()-started,3)
            store.update(run_id,status='succeeded',current_stage='completed',completed_at=_now(),runtime_seconds=runtime,
                         result_json=json.dumps({'final_text':final_text,'requested_action':'reply','explanation_evidence':explanation},ensure_ascii=False))
            store.append_message(conversation_id,'assistant',final_text,run_id)
            store.event(run_id,'reference_review',{'name':'参考资料回答完成','decision':explanation['decision'],'fact_ids':explanation['fact_ids'],'version':explanation['knowledge_version']})
            store.event(run_id,'completed',{'runtime_seconds':runtime,'has_final_text':True,'computation_started':False})
            return
        receipts = ExecutionReceipts(run_id, emit=lambda kind,payload:store.event(run_id,kind,payload), recompute=plan.recompute, allow_repeats=plan.repeated_trials)
        agent = build_agent(tools=tools, system_prompt=system_prompt, middleware=[receipts])
        store.update(run_id, current_stage="agent_execution")
        store.event(run_id, "stage", {"name": "agent_execution"})
        # Stream state updates so the browser can show tool-level progress while
        # the model is working. The final state has the same shape as invoke().
        result = None
        method_review = None
        coverage = None
        seen_messages: set[str] = set()
        for state in agent.stream(
            {"messages": conversation_messages or [{"role": "user", "content": prompt}]},
            config={"recursion_limit": RECURSION_LIMIT},
            stream_mode="values",
        ):
            check_cancelled()
            if not isinstance(state, dict):
                continue
            result = state
            for index, message in enumerate(state.get("messages", [])):
                message_id = _message_field(message, "id") or f"{index}:{_message_field(message, 'type')}:{_content_text(_message_field(message, 'content'))[:80]}"
                if message_id in seen_messages:
                    continue
                seen_messages.add(message_id)
                _capture_execution_calls(run_id, message, recompute=plan.recompute)
                _capture_operation_result(run_id, message)
                calls = _message_field(message, 'tool_calls', []) or []
                event_messages = [{"type": "ai", "tool_calls": [call]} for call in calls] if calls else [message]
                for event_message in event_messages:
                    event = _stream_event_payload(event_message)
                    if event is not None:
                        event_type, payload = event
                        store.event(run_id, event_type, payload)
            if lifecycle.handoff:
                break
            coverage = domain_contract_progress(contract, _json_safe(state))
            if coverage and coverage['complete']:
                store.event(run_id, 'completion_checked', coverage)
                break
            method_review = method_rejection(_json_safe(state)) if plan.action == 'execute' else None
            if method_review:
                store.event(run_id, 'stage', {'name': '方法评估完成：不适用于所请求目标', 'status': 'succeeded'})
                break
        if lifecycle.handoff:
            _execute_continuation(run_id, lifecycle.handoff, started)
            return
        if result is None:
            raise RuntimeError("agent stream returned no final state")
        safe_result = _json_safe(result)
        final_text = _final_text(safe_result)
        review = method_review or review_answer(safe_result, final_text, scope='conversation' if plan.action in {'reply','inspect'} else 'scientific')
        if review:
            final_text = review['final_text']
        explanation_evidence = None
        if reference_bundle and reference_bundle['cards']:
            store.event(run_id, 'reference_review_started', {'name':'核对解释与来源', 'status':'running'})
            explanation_evidence = review_explanation(prompt, final_text, reference_bundle, context=conversation_messages[-6:])
            final_text = explanation_evidence['final_text']
            store.event(run_id, 'reference_review', {'name':'解释来源核对完成', 'decision':explanation_evidence['decision'],
                        'fact_ids':explanation_evidence['fact_ids'], 'version':explanation_evidence['knowledge_version']})
        claim_evidence = None
        if not method_review:
            from .claim_evidence import collect_bundle, compile_answer
            needs_records = bool(project_snapshot and (project_snapshot.get('candidate_ledger') or {}).get('measurement_count') and plan.action == 'inspect' and plan.inspection != 'status')
            evidence_bundle = collect_bundle(store, run_id, safe_result, root=ROOT, find_manifest=_find_viewer_manifest, expect_evidence=needs_records)
            if evidence_bundle:
                store.event(run_id, 'stage', {'name':'核对结论与证据', 'status':'running'})
                claim_evidence = compile_answer(prompt, final_text, evidence_bundle, context=[m for m in conversation_messages[-6:] if m['role'] in {'user','assistant'}])
                final_text = claim_evidence['final_text']
                store.event(run_id, 'stage', {'name':'核对结论与证据', 'status':'succeeded', 'decision':claim_evidence['decision'], 'claims':len(claim_evidence['claims'])})
        operation_status = operation_outcome(store.get(run_id) or {}) if plan.action == 'execute' else None
        if operation_status in {'failed','partial'}:
            review = {**review, 'run_status':operation_status, 'operation_outcome':operation_status}
        if coverage and not coverage['complete']:
            review = {**review, 'run_status':'partial', 'completion':coverage}
            final_text += '\n\n尚未取得以下输入的结构域查询结果：' + '、'.join(coverage['missing_inputs']) + '。'
        if plan.action == 'execute' and not method_review and not operation_status:
            evidence = [value for _, value in scientific_payloads(safe_result) if value.get('status') in {'succeeded','resumed_existing','partial'}]
            if not evidence and (not review or review.get('run_status') != 'failed'):
                review = {**review, 'run_status':'partial', 'decision':'awaiting_execution_evidence'}
        runtime = round(time.perf_counter() - started, 3)
        store.update(
            run_id,
            status=review['run_status'] if review else "succeeded",
            completed_at=_now(),
            current_stage="failed" if review and review['run_status'] == 'failed' else "completed",
            error=final_text if review and review['run_status'] == 'failed' else None,
            result_json=json.dumps({"agent": safe_result, "final_text": final_text, "requested_action": plan.action,
                                    **({'completion':coverage} if coverage else {}),
                                    **({'claim_evidence':claim_evidence} if claim_evidence else {}),
                                    **({'explanation_evidence':explanation_evidence} if explanation_evidence else {}),
                                    **({'scientific_review': review} if review else {})}, ensure_ascii=False),
            runtime_seconds=runtime,
        )
        store.append_message(conversation_id, "assistant", final_text or "任务已完成，但 Agent 没有返回文字报告。", run_id)
        store.event(run_id, "failed" if review and review['run_status'] == 'failed' else "completed",
                    {"runtime_seconds": runtime, "has_final_text": bool(final_text)})
    except Exception as exc:  # the task record must preserve a structured failure
        runtime = round(time.perf_counter() - started, 3)
        error = _execution_error(exc)
        store.update(
            run_id,
            status="failed",
            completed_at=_now(),
            current_stage="failed",
            error=error,
            runtime_seconds=runtime,
        )
        store.append_message(conversation_id, "assistant", f"任务执行失败：{error}", run_id)
        store.event(run_id, "failed", {"runtime_seconds": runtime, "error": error})


@app.get("/health")
def health() -> dict[str, Any]:
    return {
        "status": "ok",
        "service": "marine-scientific-agent-api",
        "state_db": str(STATE_DB),
        "max_concurrency": MAX_CONCURRENCY,
        "control_concurrency": CONTROL_CONCURRENCY,
        "openfold3_mps_policy": "serialized file lock",
        "intent_routing": "deterministic hints plus Deep Agent planning",
        "protein_screening": "new FASTA/reference inputs; versioned BLAST; hash-validated stage recovery",
        "domain_evidence": "exact-sequence InterPro lookup; release-pinned chunks; evidence-rendered reports",
        "task_lifecycle": "conversation-scoped state and linked continuation; saved inputs; Vina partial-result reuse",
        "uniprot_resolver": "search/resolve API with atomic local cache and explicit confirmation for ambiguous hits",
        "configured_model": MODEL_NAME,
        "vina_workbench": "PDBQT validation + AutoDock Vina batch execution",
        "vina_execution": "registered Mac worker over Tailscale" if VINA_WORKER_URL else "unconfigured",
        "vina_worker_url": VINA_WORKER_URL or None,
    }


@app.post("/internal/vina-progress/{run_id}")
def vina_progress(run_id: str, event: dict[str, Any], x_vina_worker_callback_token: str | None = Header(default=None)) -> dict[str, str]:
    """Receive bounded progress events from the private Mac Vina worker."""
    if not VINA_CALLBACK_TOKEN or x_vina_worker_callback_token != VINA_CALLBACK_TOKEN:
        raise HTTPException(status_code=401, detail="unauthorized")
    record = store.get(run_id)
    if record is None or record.get("mode") != "vina":
        raise HTTPException(status_code=404, detail="vina run not found")
    if record.get('status') in TERMINAL:
        return {'status': 'ignored_terminal'}
    def count(key: str) -> int:
        try:
            return max(0, min(100, int(event.get(key) or 0)))
        except (TypeError, ValueError):
            raise HTTPException(status_code=422, detail=f'invalid progress {key}')
    safe_event = {
        "stage": str(event.get("stage") or "对接计算")[:80],
        "name": str(event.get("name") or event.get("stage") or "对接计算")[:160],
        "ligand_id": str(event.get("ligand_id") or "")[:120],
        "index": count('index'),
        "total": count('total'),
    }
    if event.get('status') in {'running', 'succeeded', 'failed', 'timed_out', 'stopped', 'skipped'}:
        safe_event['status'] = event['status']
    for key in ['completed', 'succeeded', 'failed']:
        if key in event:
            safe_event[key] = count(key)
    store.update(run_id, current_stage=safe_event["name"])
    store.event(run_id, "stage", safe_event)
    return {"status": "ok"}


@app.post("/runs", response_model=RunResponse, status_code=status.HTTP_202_ACCEPTED)
def create_run(request: RunRequest) -> RunResponse:
    mode = _resolve_mode(request.prompt, request.mode)
    if mode == "manual_openfold" and request.manual is None:
        raise HTTPException(status_code=400, detail="manual_openfold requires a confirmed accession")
    manual = request.manual.model_dump() if request.manual is not None else None
    try:
        record = store.create(request.prompt, mode, manual=manual, conversation_id=request.conversation_id, request_id=request.request_id, fingerprint=hashlib.sha256(json.dumps(request.model_dump(exclude={'request_id'}),sort_keys=True).encode()).hexdigest(), require_idle=True)
    except SubmissionConflict as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    if record.pop("_reused", False): return RunResponse(**record)
    store.event(record["run_id"], "intent_classified", classify_intent(request.prompt))
    _submit(record["run_id"], _dispatch_run, request.prompt, mode)
    return RunResponse(**record)


@app.post("/protein/runs", response_model=RunResponse, status_code=status.HTTP_202_ACCEPTED)
def create_protein_run(request: ProteinScreenRunRequest) -> RunResponse:
    payload = request.model_dump(exclude={"conversation_id"})
    try:
        validate_request(payload)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    try:
        record = store.create("蛋白同源筛选：" + request.objective, "general", conversation_id=request.conversation_id)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    spec = {"workflow": "protein_screen", "input": payload}
    _save_execution_spec(record["run_id"], spec)
    _submit(record["run_id"], _execute_protein_spec, spec)
    return RunResponse(**record)


@app.post("/internal/protein-progress/{run_id}")
def protein_progress(run_id: str, event: dict[str, Any], x_vina_worker_callback_token: str | None = Header(default=None)):
    if not VINA_CALLBACK_TOKEN or x_vina_worker_callback_token != VINA_CALLBACK_TOKEN:
        raise HTTPException(status_code=401, detail="unauthorized")
    spec = (store.get(run_id) or {}).get("execution_spec") or {}
    if spec.get("workflow") != "protein_screen":
        raise HTTPException(status_code=404, detail="protein run not found")
    try:
        _protein_progress(run_id, event)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    return {"status": "ok"}


@app.post('/protein/domains/runs', response_model=RunResponse, status_code=status.HTTP_202_ACCEPTED)
def create_domain_run(request: DomainRunRequest) -> RunResponse:
    payload = request.model_dump(exclude={'conversation_id'})
    try:
        resolve_inputs(payload)
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    try:
        record = store.create('蛋白结构域证据核验', 'general', conversation_id=request.conversation_id)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    spec = {'workflow': 'protein_domains', 'input': payload}
    _save_execution_spec(record['run_id'], spec)
    _submit(record['run_id'], _execute_domain_spec, spec)
    return RunResponse(**record)


register_project_routes(app, lambda: store)


@app.get("/conversations")
def list_conversations(limit: int = Query(default=50, ge=1, le=100)) -> dict[str, Any]:
    """List persistent Agent conversations, each of which may contain many runs."""
    return {"conversations": store.list_conversations(limit)}


@app.get("/conversations/{conversation_id}")
def get_conversation(conversation_id: str) -> dict[str, Any]:
    conversation = store.get_conversation(conversation_id)
    if conversation is None:
        raise HTTPException(status_code=404, detail="conversation not found")
    return conversation


class ContinueSelection(BaseModel):
    operation_id: str = ''
    operation_index: int | None = Field(default=None, ge=1)
    method: Literal['saved','full_scan'] = 'saved'


@app.post('/runs/{source_run_id}/continue', response_model=RunResponse, status_code=status.HTTP_202_ACCEPTED)
def continue_run(source_run_id: str, selection: ContinueSelection | None = None) -> RunResponse:
    """Create a linked attempt using the same lifecycle service as the Agent."""
    with store._lock:
        source = store.get(source_run_id)
        if source is None:
            raise HTTPException(status_code=404, detail='run not found')
        if not (saved_spec(source) or source.get('operations')) or source['status'] not in {'stopped','failed','partial'}:
            raise HTTPException(status_code=409, detail='该任务没有可继续的中断计算')
        conversation = store.get_conversation(source['conversation_id'])
        if not source.get('operations') and not store.can_continue(source):
            raise HTTPException(status_code=409, detail='该任务已有后续记录，请使用最新一次任务')
        if any(row['status'] not in TERMINAL for row in conversation['runs']):
            raise HTTPException(status_code=409, detail='当前会话已有任务在执行或等待停止')
        record = store.create('继续上次任务', 'general', conversation_id=source['conversation_id'])
        lifecycle = ConversationTasks(store, record['run_id'])
        selected = selection or ContinueSelection()
        receipt = lifecycle.prepare(source_run_id, selected.operation_id, selected.operation_index, selected.method)
        if receipt['status'] != 'handoff_ready':
            store.update(record['run_id'], status='failed', error=receipt['reason'], completed_at=_now())
            raise HTTPException(status_code=409, detail=receipt['reason'])
        _submit(record['run_id'], _execute_continuation, lifecycle.handoff)
        return RunResponse(**record)


@app.post("/vina/runs", response_model=RunResponse, status_code=status.HTTP_202_ACCEPTED)
def create_vina_run(request: VinaRunRequest) -> RunResponse:
    """Queue one standards-aligned Vina batch and return a durable run ID."""
    payload = request.model_dump()
    # Validate synchronously so malformed files and boxes are shown immediately
    # in the workbench instead of becoming an opaque background failure.
    try:
        resolve_library_ligands(payload)
        prep_id = (payload.get('receptor') or {}).get('preparation_id')
        if prep_id and store.get(str(prep_id)) is None:
            raise HTTPException(404, '受体准备任务不存在')
        preparation = resolve_prepared_receptor(payload, store)
        if preparation:
            payload['receptor_preparation'] = preparation
        validate_config(payload)
    except LibraryUnavailable as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except (ValueError, OSError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    try:
        record = store.create(f"AutoDock Vina 分子对接（{len(payload['ligands'])} 个配体）", "vina", manual=payload, conversation_id=request.conversation_id, request_id=request.request_id, fingerprint=hashlib.sha256(json.dumps(request.model_dump(exclude={'request_id'}),sort_keys=True).encode()).hexdigest(), require_idle=True)
    except SubmissionConflict as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    if record.pop("_reused", False): return RunResponse(**record)
    store.event(record["run_id"], "intent_classified", {"workflow": "vina_docking", "confirmation_required": False})
    _submit(record["run_id"], _execute_vina_run, payload)
    return RunResponse(**record)


register_receptor_routes(app, store, executor, lambda agent: trusted_structure_manifest(agent, _find_viewer_manifest))


@app.get("/vina/library")
def get_ligand_library() -> dict[str, Any]:
    try:
        return public_catalog()
    except LibraryUnavailable as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc


@app.get("/vina/demo")
def get_vina_demo() -> dict[str, Any]:
    """Return a frozen, traceable PDBQT example for a first UI smoke test."""
    receptor = RESULTS / "vina-demo" / "receptor.pdbqt"
    ligand = RESULTS / "vina-demo" / "ligand.pdbqt"
    if not receptor.is_file() or not ligand.is_file():
        raise HTTPException(status_code=404, detail="Vina 示例输入尚未部署")
    return {
        "status": "succeeded",
        "source": "local validated 1IEP Vina example",
        "receptor": {"name": "1IEP-chain-A-receptor.pdbqt", "content_base64": base64.b64encode(receptor.read_bytes()).decode("ascii")},
        "ligands": [{"id": "imatinib-control", "name": "imatinib-control.pdbqt", "content_base64": base64.b64encode(ligand.read_bytes()).decode("ascii")}],
        "box": {"center": [15.614, 53.380, 15.455], "size": [30.0, 30.0, 30.0]},
        "params": {"exhaustiveness": 8, "num_modes": 9, "energy_range": 4, "seed": 20261028, "cpu": 1},
        "provenance": "PDB 1IEP 实验受体链 A；配体采用已冻结的独立 PDBQT 控制文件。",
    }


@app.get("/lookup/uniprot/search")
def lookup_uniprot_search(
    query: str = Query(min_length=2, max_length=160),
    organism: str = Query(default="", max_length=120),
    reviewed: bool | None = None,
    limit: int = Query(default=8, ge=1, le=20),
) -> dict[str, Any]:
    """Resolver endpoint for the manual structure workflow."""
    return json.loads(
        search_uniprot.invoke(
            {"query": query, "organism": organism, "reviewed": reviewed, "limit": limit}
        )
    )


@app.get("/lookup/uniprot/{accession}")
def lookup_uniprot_sequence(accession: str) -> dict[str, Any]:
    return json.loads(get_uniprot_sequence.invoke({"accession": accession}))


@app.get("/runs")
def list_runs(
    limit: int = Query(default=50, ge=1, le=100),
    offset: int = Query(default=0, ge=0),
) -> dict[str, Any]:
    """List durable runs for the authenticated marine workbench."""
    return {"runs": store.list(limit=limit, offset=offset), "limit": limit, "offset": offset}


@app.get("/runs/{run_id}")
def get_run(run_id: str) -> dict[str, Any]:
    record = store.get(run_id)
    if record is None:
        raise HTTPException(status_code=404, detail="run not found")
    record['can_continue'] = store.can_continue(record)
    record.pop('execution_spec', None)
    from .task_lifecycle import operation_facts
    record['operations'] = [operation_facts(op) for op in record.get('operations', [])]
    if record['mode']=='vina':
        record.pop('manual_input', None)
    return record


@app.get("/vina/runs/{run_id}")
def get_vina_run(run_id: str) -> dict[str, Any]:
    record = store.get(run_id)
    if record is None or record.get("mode") != "vina":
        raise HTTPException(status_code=404, detail="vina run not found")
    # Do not send base64 receptor/ligand payloads on every browser poll. The
    # report endpoint returns the compact result after completion.
    return {key: value for key, value in record.items() if key not in {"manual_input", "result", "execution_spec"}} | {
        'can_continue':store.can_continue(record)}


@app.get("/runs/{run_id}/events")
def get_run_events(run_id: str, after_id: int = 0) -> dict[str, Any]:
    if store.get(run_id) is None:
        raise HTTPException(status_code=404, detail="run not found")
    return {"run_id": run_id, "events": store.events(run_id, after_id)}


@app.get("/runs/{run_id}/report")
def get_run_report(run_id: str) -> dict[str, Any]:
    record = store.get(run_id)
    if record is None:
        raise HTTPException(status_code=404, detail="run not found")
    project_context = store.run_project_context(run_id)
    if record["status"] == "failed" and not (record.get("result") or {}).get("final_text"):
        return {"run_id": run_id, "status": "failed", "error": record.get("error"), "project_context": project_context}
    result = record.get("result") or {}
    if record.get("mode") == "vina":
        return {
            "run_id": run_id,
            "status": record["status"],
            "runtime_seconds": record.get("runtime_seconds"),
            "final_text": result.get("final_text", ""),
            "vina": result.get("vina"),
            "project_context": project_context,
            "execution": execution_receipt(record),
        }
    structure = verified_structure(record, ROOT, _find_viewer_manifest)
    manifest_path = structure.get("path") if structure and structure["status"] == "ready" else None
    return {
        "run_id": run_id,
        "status": record["status"],
        "runtime_seconds": record.get("runtime_seconds"),
        "final_text": result.get("final_text", ""),
        "agent_result": result.get("agent"),
        "scientific_review": result.get("scientific_review"),
        "explanation_evidence": result.get("explanation_evidence"),
        "has_claim_evidence": bool(result.get("claim_evidence")),
        "project_context": project_context,
        "execution": execution_receipt(record, structure),
        "structure_artifact": {
            "type": "openfold3",
            "status": "ready",
            "manifest_url": _artifact_url(run_id, manifest_path),
            "path": manifest_path,
            **structure,
        } if manifest_path else structure,
    }


@app.get('/runs/{run_id}/evidence')
def get_claim_evidence(run_id: str):
    from .claim_evidence import public_evidence
    record = store.get(run_id)
    if record is None:
        raise HTTPException(404, '任务不存在')
    evidence = (record.get('result') or {}).get('claim_evidence')
    if not evidence:
        raise HTTPException(404, '该记录没有结论核验信息')
    return public_evidence(evidence)


@app.get("/runs/{run_id}/artifact/{artifact_path:path}")
def get_run_artifact(run_id: str, artifact_path: str):
    """Serve a completed run's viewer files without exposing host paths."""
    if store.get(run_id) is None:
        raise HTTPException(status_code=404, detail="run not found")
    normalized = artifact_path.replace("\\", "/").lstrip("/")
    if not normalized.startswith("results/"):
        raise HTTPException(status_code=404, detail="artifact not found")
    path = (ROOT / normalized).resolve()
    try:
        path.relative_to(ROOT.resolve())
    except ValueError as exc:
        raise HTTPException(status_code=404, detail="artifact not found") from exc
    if not path.is_file() or not _owns_artifact(store.get(run_id), path):
        raise HTTPException(status_code=404, detail="artifact not found")
    return FileResponse(path)


def _owns_artifact(record: dict, path: Path) -> bool:
    """Bind files to trusted execution outputs, never user prose or arbitrary paths."""
    roots = [RESULTS / 'vina-runs' / record['id'], RESULTS / 'vina-receptors' / record['id']]
    agent = (record.get('result') or {}).get('agent') or {}
    trusted = [agent[key] for key in ('manual_openfold', 'inline_openfold') if key in agent]
    for message in agent.get('messages', []):
        if message.get('type', message.get('role')) in ('tool', 'tool_result') and message.get('name') in ('run_openfold_accession', 'run_openfold_sequence'):
            trusted.append(message.get('content'))
    for output in trusted:
        manifest = _find_viewer_manifest(output)
        if manifest:
            directory = (ROOT / manifest).resolve().parent
            if directory.is_relative_to(RESULTS.resolve()):
                roots.append(directory)
    for output in [agent.get('protein_domains', {}), agent.get('protein_screen', {})]:
        analysis_id = output.get('analysis_id', '')
        if re.fullmatch('[a-f0-9]{64}', analysis_id):
            roots.extend(RESULTS / kind / analysis_id for kind in ('protein-domains', 'protein-screenings'))
    for operation in record.get('operations', []):
        analysis_id = operation.get('facts', {}).get('analysis_id', '')
        if re.fullmatch('[a-f0-9]{64}', analysis_id):
            roots.extend(RESULTS / kind / analysis_id for kind in ('protein-domains', 'protein-screenings'))
    return any(path.is_relative_to(root.resolve()) for root in roots)
