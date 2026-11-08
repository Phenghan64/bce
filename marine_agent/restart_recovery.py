"""Conservative reconciliation after loss of the single API coordinator.

We do not label an orphaned computation successful or automatically replay it.
Remote cancellation is acknowledged and then checked for quiescence before a
saved-input continuation can be offered. Unreachable workers stay unresolved.
"""
from __future__ import annotations

import json
import os
import time
from urllib.request import Request, build_opener, ProxyHandler


def remote_services(record):
    spec = record.get('execution_spec') or {}
    workflows = {spec.get('workflow')}
    workflows.update(op.get('spec', {}).get('workflow') for op in record.get('operations', []))
    # Tool inputs are durably registered before dispatch. A `multiple` summary
    # refers to those recorded operations, not to arbitrary additional workers.
    known = {'openfold_accession', 'openfold_sequence', 'vina', 'protein_screen', 'protein_domains'}
    uncertain = not bool(workflows & known) or bool(workflows - known - {None, 'multiple'})
    services = []
    if uncertain or workflows & {'openfold_accession', 'openfold_sequence'} or record.get('mode') == 'manual_openfold':
        services.append(('OpenFold', os.getenv('OPENFOLD_WORKER_URL', ''), os.getenv('OPENFOLD_WORKER_TOKEN', '')))
    if uncertain or 'vina' in workflows or record.get('mode') == 'vina':
        services.append(('Vina', os.getenv('VINA_WORKER_URL', ''), os.getenv('VINA_WORKER_TOKEN', '')))
    if uncertain or 'protein_screen' in workflows:
        services.append(('Vina', os.getenv('PROTEIN_WORKER_URL', os.getenv('VINA_WORKER_URL', '')), os.getenv('PROTEIN_WORKER_TOKEN', os.getenv('VINA_WORKER_TOKEN', ''))))
    specs = [spec, *(op.get('spec', {}) for op in record.get('operations', []))]
    if uncertain or any(s.get('workflow') == 'protein_domains' and s.get('input', {}).get('method') == 'full_scan' for s in specs):
        services.append(('InterPro', os.getenv('INTERPRO_WORKER_URL', ''), os.getenv('INTERPRO_WORKER_TOKEN', '')))
    return list(dict.fromkeys(service for service in services if service[1]))


def reconcile(store, record, *, attempts=20, delay=1, client=None):
    client = client or build_opener(ProxyHandler({}))
    run_id = record['id']
    services = remote_services(record)
    for attempt in range(attempts):
        settled = True
        for name, url, token in services:
            headers = {f'X-{name}-Worker-Token': token, 'Content-Type': 'application/json'}
            endpoint = f'{url.rstrip("/")}/v1/jobs/{run_id}'
            try:
                # The cancellation tombstone also rejects a delayed submission.
                with client.open(Request(endpoint + '/cancel', data=b'{}', headers=headers), timeout=3) as response:
                    response.read()
                with client.open(Request(endpoint, headers=headers), timeout=3) as response:
                    state = json.load(response)
                if state.get('job_id') != run_id or state.get('state') != 'stopped' or state.get('active') is not False:
                    settled = False
            except Exception:
                settled = False
        if settled:
            if record.get('mode') == 'vina':
                try:
                    from .vina_remote import recover_vina_artifacts
                    recovered = recover_vina_artifacts(run_id)
                    if recovered:
                        store.save_recovered_vina(run_id, recovered)
                except Exception:
                    store.event(run_id, 'recovery_pending', {'name': '旧作业已停止，等待核验保存的对接产物'})
                    if attempt + 1 < attempts:
                        time.sleep(delay)
                    continue
            store.finish_recovery(run_id)
            return True
        if attempt + 1 < attempts:
            time.sleep(delay)
    store.event(run_id, 'recovery_pending', {'name': '服务已重启，等待工作节点确认停止；暂不重复计算'})
    return False
