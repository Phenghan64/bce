"""Server-owned completion facts. Dialogue is never computational validation."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from .evidence_review import tool_payloads

STRUCTURE_TOOLS = {'run_openfold_accession', 'run_openfold_sequence'}
CONTROL_TOOLS = {'analyze_request_intent', 'assess_collagenase_method', 'inspect_task_state', 'inspect_project_task', 'inspect_candidate_evidence', 'inspect_saved_input_identity', 'write_todos', 'prepare_task_continuation'}

def scientific_payloads(agent):
    return [(name, value) for name, value in tool_payloads(agent) if name not in CONTROL_TOOLS]

def trusted_structure_manifest(agent, find):
    if not isinstance(agent, dict):
        return None
    outputs = [agent[key] for key in ('manual_openfold', 'inline_openfold') if isinstance(agent.get(key), dict)]
    outputs += [value for name, value in tool_payloads(agent) if name in STRUCTURE_TOOLS]
    for output in reversed(outputs):
        found = find(output)
        if found:
            return found
    return None

def verified_structure(record, root: Path, find):
    """A real tool path is necessary but insufficient: check files and digest."""
    candidate = trusted_structure_manifest((record.get('result') or {}).get('agent'), find)
    if not candidate:
        return None
    try:
        manifest_path = (root/candidate).resolve()
        manifest_path.relative_to((root/'results').resolve())
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('schema_version') != 1 or manifest.get('provenance',{}).get('model') != 'OpenFold3':
            raise ValueError('Unsupported structure provenance')
        item = manifest['structure']; structure = (manifest_path.parent/item['url']).resolve()
        structure.relative_to(manifest_path.parent)
        if not item.get('sha256') or hashlib.sha256(structure.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError('Structure digest mismatch')
        for key in ('confidence', 'validation'):
            if manifest.get(key):
                path=(manifest_path.parent/manifest[key]['url']).resolve();path.relative_to(manifest_path.parent)
                json.loads(path.read_text())
        return {'status':'ready', 'path':str(manifest_path.relative_to(root.resolve())), 'id':manifest.get('id'),
                'verification':'files_and_structure_sha256', 'source_run_id':record['id']}
    except (OSError, ValueError, KeyError, TypeError):
        return {'status':'unavailable','source_run_id':record['id'], 'reason':'结构产物缺失或校验失败，请检查保存结果。'}

def execution_receipt(record, structure=None):
    result=record.get('result') or {};review=result.get('scientific_review') or {}
    ops=record.get('operations') or []
    payloads=scientific_payloads(result.get('agent') or {})
    successful=[name for name,value in payloads if value.get('status') in {'succeeded','resumed_existing'}]
    action=result.get('requested_action')
    scientific=bool(ops or record.get('execution_spec') or record.get('mode') in {'vina','manual_openfold'} or successful)
    if review.get('decision')=='method_rejected':
        kind,status,label='method_assessment','not_performed','方法不适用，未执行计算'
    elif scientific:
        kind,status,label='scientific',record['status'],'计算记录'
        if structure and structure['status'] != 'ready':status,label='artifact_unavailable','产物待核验'
        elif structure:label='结构产物已核验'
        elif status=='succeeded':label='工具执行完成'
    elif action=='execute':
        kind,status,label='scientific','not_started','尚未取得计算结果'
    else:
        kind,status,label='dialogue','not_requested','对话回复'
    return {'source':'server_execution_records', 'run_id':record['id'], 'turn_status':record['status'],
            'kind':kind, 'computation_status':status, 'label':label,
            'operation_count':len(ops), 'successful_tools':successful,
            'artifact_verified':bool(structure and structure['status']=='ready'),
            'scientific_validity':'not_established_by_execution_status'}
