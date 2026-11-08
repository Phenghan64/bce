"""Conversation-scoped task facts and deterministic continuation handoffs.

The model selects an action; the coordinator owns execution, saved inputs and
state transitions. A continuation is a new attempt, never a resurrected process.
"""
from __future__ import annotations

import copy
import json
import threading
from typing import Any

from langchain_core.tools import tool

from .execution_control import TERMINAL, check_cancelled


def saved_spec(record: dict) -> dict | None:
    spec = record.get("execution_spec")
    if isinstance(spec, dict):
        return copy.deepcopy(spec) if spec.get("workflow") in {
            "openfold_accession", "openfold_sequence", "vina", "protein_screen", "protein_domains"
        } else None
    # Older manual jobs already have authoritative, reproducible inputs.
    manual = record.get("manual_input")
    if isinstance(manual, dict):
        if record.get("mode") == "manual_openfold" and manual.get("accession"):
            return {"workflow": "openfold_accession", "input": copy.deepcopy(manual)}
        if record.get("mode") == "vina" and manual.get("receptor") and manual.get("ligands"):
            return {"workflow": "vina", "input": copy.deepcopy(manual)}
    return None


def describe_spec(spec: dict | None) -> dict | None:
    if not spec:
        return None
    inputs = spec.get("input", {})
    workflow = spec["workflow"]
    if workflow == "protein_domains":
        fasta = inputs.get('query_fasta') or ''
        lengths = []
        for entry in fasta.split('>'):
            lines = entry.strip().splitlines()
            if len(lines) > 1:
                lengths.append({'protein_id': lines[0].split()[0], 'residue_count': len(''.join(''.join(lines[1:]).split()))})
        return {"workflow": workflow, "input_sequences": lengths[:100], "source_analysis_id": inputs.get("analysis_id"),
                "protein_id": inputs.get("protein_id"), "source": inputs.get("source"),
                "method": inputs.get("method", "lookup"),
                "continuation_method": "nextflow_validated_resume" if inputs.get("method") == "full_scan" else "reuse_valid_domain_chunks_for_same_release"}
    if workflow == "protein_screen":
        return {"workflow": workflow, "objective": inputs.get("objective"),
                "query_source": inputs.get("query_source"), "reference_source": inputs.get("reference_source"),
                "continuation_method": "reuse_hash_validated_stages_retry_remaining"}
    if workflow == "openfold_accession":
        return {"workflow": workflow, "accession": inputs.get("accession"),
                "ligand_count": len(inputs.get("ligands", [])), "continuation_method": "restart_prediction_with_cache"}
    if workflow == "openfold_sequence":
        return {"workflow": workflow, "residue_count": len(inputs.get("sequence", "")),
                "ligand_count": len(inputs.get("ligands", [])), "continuation_method": "restart_prediction_with_cache"}
    return {"workflow": workflow, "ligand_count": len(inputs.get("ligands", [])),
            "params": inputs.get("params"), "continuation_method": "retain_valid_completed_ligands_retry_remaining"}


def resumable_operations(record):
    return [op for op in record.get('operations', []) if op.get('status') in {'stopped','failed','partial','unknown'}]


def canonical_record(record):
    """Retain audit events but do not number same-turn reused calls as new inputs."""
    if not record:
        return record
    return {**record, 'operations':[op for op in record.get('operations', [])
            if not op.get('facts', {}).get('execution_reference', {}).get('reused')]}


def select_operation(record, operation_id='', operation_index=None):
    operations = record.get('operations', [])
    if operation_index is not None:
        if operation_index < 1 or operation_index > len(operations):
            return None
        selected = operations[operation_index-1]
        return selected if not operation_id or operation_id == selected['operation_id'] else None
    if operation_id:
        return next((op for op in operations if op['operation_id'] == operation_id), None)
    pending = resumable_operations(record)
    return pending[0] if len(pending) == 1 else None


def operation_facts(operation):
    return {key: operation[key] for key in ('operation_id', 'tool', 'status', 'facts') if key in operation} | {
        'execution': describe_spec(operation.get('spec'))}


class ConversationTasks:
    def __init__(self, store, current_run_id: str):
        self.store = store
        self.current_run_id = current_run_id
        self.conversation_id = (store.get(current_run_id) or {}).get("conversation_id")
        self.handoff: dict | None = None
        self._lock = threading.RLock()

    def records(self) -> list[dict]:
        conversation = self.store.get_conversation(self.conversation_id) if self.conversation_id else None
        return [canonical_record(self.store.get(row["id"])) for row in (conversation or {}).get("runs", [])
                if row["id"] != self.current_run_id]

    def has_resumable_work(self):
        return any(row['status'] in {'stopped','failed','partial'} and (saved_spec(row) or resumable_operations(row)) for row in self.records())

    def resolve_reference(self, reference, records=None):
        """Resolve only recorded, conversation-owned identifiers; never guess by shape."""
        if records is None:
            current = canonical_record(self.store.get(self.current_run_id))
            records = self.records() + ([current] if current else [])
        direct = [row for row in records if row['id'] == reference]
        if len(direct) == 1:
            return {'record':direct[0], 'operation':None}
        matches = [(row,op) for row in records for op in row.get('operations', [])
                   if reference == op['operation_id'] or reference == op.get('facts', {}).get('analysis_id')]
        if len(matches) == 1:
            return {'record':matches[0][0], 'operation':matches[0][1]}
        return {'status':'ambiguous_reference' if matches else 'not_found',
                'reason':'该编号对应多个操作，请选择明确的任务和操作。' if matches else '当前会话没有对应的任务、操作或分析编号。',
                'references':[{'run_id':r['id'], 'operation_id':o['operation_id'], 'analysis_id':o.get('facts',{}).get('analysis_id')}
                              for r,o in matches]}

    def snapshot(self, run_id: str = "") -> dict:
        records = self.records()
        current = canonical_record(self.store.get(self.current_run_id))
        if current and (saved_spec(current) or current.get('operations')):
            records.append(current)
        if run_id:
            resolved = self.resolve_reference(run_id, records)
            if 'record' not in resolved:
                return resolved
            records = [resolved['record']]
        scientific = [row for row in records if saved_spec(row) or row.get('operations') or row.get('execution_spec')]
        selected = {row['id'] for row in records[-8:] + scientific[-8:]}
        superseded = {row.get("continuation_of") for row in records if row.get("continuation_of")}
        facts = []
        for row in records:
            if row['id'] not in selected:
                continue
            spec = saved_spec(row)
            facts.append({"run_id": row["id"], "status": row["status"],
                          "completion_scope": "tool_execution_only" if row in scientific else "conversation_turn_only",
                          "scientific_validation": "not_established", "mode": row["mode"],
                          "stage": row["current_stage"], "created_at": row["created_at"],
                          "continuation_of": row.get("continuation_of"), "execution": describe_spec(spec),
                          "can_continue": bool(spec and row["id"] not in superseded and row["status"] in {"stopped", "failed", "partial"}),
                          "error": row.get("error"),
                          "operations": [operation_facts(op) for op in row.get('operations', [])],
                          "completed_ligands": sum(r.get("status") == "succeeded" for r in
                              ((row.get("result") or {}).get("vina") or {}).get("results", []))})
            if row in scientific and (run_id or row is scientific[-1]):
                facts[-1]['saved_report'] = str((row.get('result') or {}).get('final_text') or '')[:6000]
        return {"status": "ok", "source": "persistent_server_task_store", "scope":"current_conversation", "conversation_id":self.conversation_id, "tasks": facts}

    def prepare(self, source_run_id: str = "", operation_id: str = "", operation_index: int | None = None, method: str = "saved") -> dict:
        check_cancelled()
        with self._lock, self.store._lock:
            if self.handoff:
                if source_run_id and source_run_id != self.handoff["source_run_id"]:
                    return {"status": "conflict", "reason": "本轮已经选定另一个继续任务"}
                return self._receipt()
            records = self.records()
            if any(row["status"] not in TERMINAL for row in records):
                return {"status": "already_running", "reason": "当前会话已有任务在执行或等待停止，请先查看其状态"}
            superseded = {row.get('continuation_of') for row in records if row.get('continuation_of')}
            candidates = [row for row in records if (saved_spec(row) or row.get('operations')) and (row['id'] not in superseded or row.get('operations'))]
            if source_run_id:
                resolved = self.resolve_reference(source_run_id, candidates)
                if 'record' not in resolved:
                    return resolved
                candidates = [resolved['record']]
                if resolved['operation'] and not operation_id and operation_index is None:
                    operation_id = resolved['operation']['operation_id']
            if operation_id:
                resolved = self.resolve_reference(operation_id, candidates)
                if 'record' not in resolved or not resolved['operation']:
                    return resolved if 'record' not in resolved else {'status':'selection_required','reason':'需要操作编号，收到的是任务编号。'}
                candidates = [resolved['record']]
                operation_id = resolved['operation']['operation_id']
            if not candidates or candidates[-1]["status"] not in {"stopped", "failed", "partial"}:
                return {"status": "not_resumable", "reason": "没有可继续的已中断科研任务；已完成结果无需重复计算"}
            source = candidates[-1]
            spec = saved_spec(source)
            selected = None
            if operation_id or operation_index is not None or not spec:
                selected = select_operation(source, operation_id, operation_index)
                if not selected:
                    return {'status':'selection_required', 'reason':'请指定一个未完成操作', 'operations':[operation_facts(op) for op in source.get('operations', [])]}
                if selected['status'] not in {'stopped','failed','partial','unknown'}:
                    return {'status':'already_completed', 'reason':'该操作已完成，不重复计算'}
                if any(row.get('continuation_of') == source['id'] and row.get('continuation_operation_id') == selected['operation_id'] for row in records):
                    return {'status':'already_continued', 'reason':'该操作已有后续记录，请使用最新一次任务'}
                spec = copy.deepcopy(selected['spec'])
            domain_unknown = selected.get('facts', {}).get('unknown_record_count', 0) if selected else (
                ((source.get('result') or {}).get('agent') or {}).get('protein_domains', {}).get('unknown_record_count', 0)
                or max([op.get('facts', {}).get('unknown_record_count', 0) for op in source.get('operations', []) if op.get('spec') == spec] or [0]))
            if method == 'full_scan' or (spec['workflow'] == 'protein_domains' and domain_unknown):
                from .interpro_remote import configured
                if spec['workflow'] != 'protein_domains':
                    return {'status':'invalid_method', 'reason':'完整 InterProScan 仅适用于蛋白结构域分析。'}
                if not configured():
                    return {'status':'capability_unavailable', 'reason':'完整 InterProScan 执行节点尚未配置；预计算查询不能替代全扫描。', 'operation_id': selected['operation_id'] if selected else None}
                spec = copy.deepcopy(spec)
                spec['input']['method'] = 'full_scan'
            self.handoff = {"source_run_id": source["id"], "spec": spec}
            if selected:
                self.handoff['operation_id'] = selected['operation_id']
                self.store.update(self.current_run_id, continuation_operation_id=selected['operation_id'])
            self.store.update(self.current_run_id, continuation_of=source["id"],
                              execution_json=json.dumps(spec, ensure_ascii=False))
            self.store.event(self.current_run_id, "stage", {
                "name": "沿用原输入继续任务", "source_run_id": source["id"],
                "execution": describe_spec(spec), "status": "succeeded"})
            return self._receipt()

    def _receipt(self) -> dict:
        return {"status": "handoff_ready", "run_id": self.current_run_id,
                "source_run_id": self.handoff["source_run_id"],
                "execution": describe_spec(self.handoff["spec"]),
                "note": "执行由任务协调器接管，进度和结果将在本轮任务中记录；不是断点续算"}

    def tools(self, *, allow_continue=True) -> list[Any]:
        @tool
        def inspect_task_state(run_id: str = "") -> str:
            """Read authoritative task states for this conversation, including stopped work and saved inputs. Empty run_id lists recent tasks. Use for progress/status questions."""
            return json.dumps(self.snapshot(run_id), ensure_ascii=False)

        @tool
        def continue_task(source_run_id: str = "", operation_id: str = "", operation_index: int | None = None, method: str = "saved") -> str:
            """Continue one saved scientific operation. operation_id or 1-based operation_index selects a subtask; method="full_scan" explicitly requests InterProScan rather than cache lookup. Continue interrupted scientific work when the user requests execution (e.g. 继续/接着跑). Reuse exact saved inputs in a new attempt; OpenFold restarts with cache, Vina retains valid finished ligands, protein screening reuses hash-validated stages. Empty source_run_id selects the latest interrupted task in this conversation. A status question alone does not request execution. The coordinator executes the returned handoff."""
            return json.dumps(self.prepare(source_run_id, operation_id, operation_index, method), ensure_ascii=False)

        @tool
        def inspect_saved_input_identity(source_run_id: str = '', sequence_index: int | None = None) -> str:
            """Verify the identity and mutation labels of saved protein inputs against UniProt and the actual input sequence. Optional 1-based sequence_index selects one sequence. Reads evidence only, does not predict a structure or activity."""
            from .input_identity import inspect_input_identity
            current = canonical_record(self.store.get(self.current_run_id))
            records = [r for r in self.records() + ([current] if current else []) if saved_spec(r) or r.get('operations')]
            if source_run_id:
                resolved = self.resolve_reference(source_run_id, records)
                if 'record' not in resolved:
                    return json.dumps(resolved,ensure_ascii=False)
                source = resolved['record']
                specs = [resolved['operation']['spec']] if resolved['operation'] else [op['spec'] for op in source.get('operations', [])] or [saved_spec(source)]
                return json.dumps(inspect_input_identity(specs, sequence_index),ensure_ascii=False)
            if not records:
                return json.dumps({'status':'not_found','reason':'当前会话没有对应的已保存输入'},ensure_ascii=False)
            source = records[-1]
            specs = [op['spec'] for op in source.get('operations', [])] or [saved_spec(source)]
            return json.dumps(inspect_input_identity(specs, sequence_index),ensure_ascii=False)

        return [inspect_task_state, *([continue_task] if allow_continue else []), inspect_saved_input_identity]
