"""Turn-scoped execution receipts, idempotency and verifiable completion contracts."""
from __future__ import annotations
import hashlib
import json
import threading
import re
from langchain.agents.middleware import AgentMiddleware
from langchain_core.messages import ToolMessage
from .access_control import account_scope, principal
from .execution_control import check_cancelled
from .protein_screening import parse_fasta

SCIENTIFIC_TOOLS = frozenset({'annotate_protein_domains', 'run_protein_homology_screen', 'run_openfold_accession', 'run_openfold_sequence'})

class ExecutionReceipts(AgentMiddleware):
    """Memoize successful exact requests only within one user's one execution turn.

    The cache never crosses a turn, database-release boundary between turns, or
    account. Refresh requests have their own key; explicit independent trials
    bypass reuse. Failed calls remain retryable.
    All arguments, including sources and scientific parameters, enter the key.
    """
    def __init__(self, run_id, emit=None, recompute=False, allow_repeats=False):
        super().__init__()
        self.run_id, self.emit, self.recompute = run_id, emit, recompute
        self.allow_repeats = allow_repeats
        self.owner = principal.get()
        self._guard = threading.Lock()
        self._locks, self._cache = {}, {}

    def wrap_tool_call(self, request, handler):
        call = request.tool_call
        if call['name'] not in SCIENTIFIC_TOOLS:
            return handler(request)
        args = call['args']
        # Materialize schema defaults; malformed calls still reach normal tool validation.
        try:
            args = request.tool.args_schema.model_validate(args).model_dump()
        except (AttributeError, ValueError, TypeError):
            pass
        if self.recompute and call['name'] == 'annotate_protein_domains':
            args = {**args, 'refresh': True}
            request = request.override(tool_call={**call, 'args': args})
        key = hashlib.sha256(json.dumps([call['name'], args], sort_keys=True, ensure_ascii=False).encode()).hexdigest()
        # Freshness refreshes the underlying data once; it is not permission
        # for the model to repeat an already completed request within this turn.
        bypass = self.allow_repeats
        with self._guard:
            lock = self._locks.setdefault(key, threading.Lock())
        # Single-flight: parallel duplicate requests wait cancellably for the owner.
        while not lock.acquire(timeout=.1):
            check_cancelled()
        try:
            check_cancelled()
            with account_scope(self.owner):
                cached = None if bypass else self._cache.get(key)
                if cached:
                    payload = json.loads(cached.content)
                    origin = payload['execution_reference']['operation_id']
                    payload['execution_reference'] = {**payload['execution_reference'], 'operation_id':call['id'], 'reused_from_operation_id':origin, 'reused':True}
                    if self.emit:
                        self.emit('tool_reused', {'tool':call['name'], 'operation_id':call['id'], 'reused_from_operation_id':origin, 'request_sha256':key})
                    return ToolMessage(content=json.dumps(payload,ensure_ascii=False),name=call['name'],tool_call_id=call['id'])
                result = handler(request)
                if not isinstance(result, ToolMessage):
                    return result
                try:
                    payload = json.loads(result.content)
                except (ValueError, TypeError):
                    return result
                if not isinstance(payload, dict):
                    return result
                payload['execution_reference'] = {'run_id':self.run_id, 'operation_id':call['id'], 'analysis_id':payload.get('analysis_id'), 'request_sha256':key, 'reused':False}
                result = result.model_copy(update={'content':json.dumps(payload,ensure_ascii=False)})
                # Stable partial domain unknowns are reusable evidence, not transport errors.
                reusable = payload.get('status') == 'succeeded' or (
                    call['name']=='annotate_protein_domains' and payload.get('status')=='partial'
                    and payload.get('analysis_id') and payload.get('release') and payload.get('records'))
                if reusable and result.status != 'error' and not bypass:
                    self._cache[key] = result
                return result
        finally:
            lock.release()


def domain_completion_contract(prompt, scope, recompute=False):
    """Only activate for a planner-declared lookup-only request with parseable input.

    Never infer an enzyme, sequence, or multi-stage objective from keywords.
    If the supplied FASTA cannot be parsed in full, leave completion to the agent.
    """
    if scope != 'provided_fasta_domain_lookup' or recompute:
        return None
    start = prompt.find('>')
    if start < 0:
        return None
    try:
        records = parse_fasta(prompt[start:])
    except ValueError:
        return None
    if not records:
        return None
    return {'kind':'provided_fasta_domain_lookup', 'expected':[(r['id'],hashlib.sha256(r['sequence'].encode()).hexdigest()) for r in records],
            'expected_pfam': sorted(set(re.findall(r'(?<![A-Z0-9])PF\d{5}(?![A-Z0-9])',prompt.upper())))}


def domain_contract_progress(contract, state):
    if not contract:
        return None
    from .evidence_review import tool_payloads
    seen = set()
    for name, value in tool_payloads(state):
        if name != 'annotate_protein_domains' or value.get('status') not in {'succeeded','partial'} or not value.get('analysis_id') or not value.get('release'):
            continue
        for row in value.get('records',[]):
            if not set(contract['expected_pfam']) <= set(row.get('expected_pfam',{})):
                continue
            seen.add((row.get('protein_id'), str(row.get('sequence_sha256','')).lower()))
    missing = [item for item in contract['expected'] if tuple(item) not in seen]
    return {'kind':contract['kind'], 'complete':not missing, 'expected_inputs':len(contract['expected']), 'covered_inputs':len(contract['expected'])-len(missing),
            'missing_inputs':[item[0] for item in missing], 'scope':'Precomputed lookup coverage only; unknown sequences remain partial and require full scan.'}


def operation_outcome(record):
    """A successful sibling cannot mask failed/unresolved required computation.

    A successful retry of identical saved inputs supersedes that input's earlier
    failure for outcome aggregation; all original operations remain in the audit.
    """
    by_input = {}
    for op in record.get('operations',[]):
        if op.get('facts',{}).get('execution_reference',{}).get('reused'):
            continue
        key = json.dumps([op['tool'],op['spec']], sort_keys=True)
        by_input[key] = op['status']
    states = list(by_input.values())
    if not states:
        return None
    if all(state == 'succeeded' for state in states):
        return 'succeeded'
    return 'failed' if all(state == 'failed' for state in states) else 'partial'
