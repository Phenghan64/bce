"""Source-bound receptor handoff, preparation jobs and docking provenance."""
from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
import time
from datetime import datetime, timezone
from typing import Callable

from fastapi import HTTPException
from pydantic import BaseModel, Field

from .vina_remote import call_receptor_worker

ROOT = Path(__file__).resolve().parents[1]


class ReceptorPreparationRequest(BaseModel):
    source_run_id: str = Field(min_length=8, max_length=64)
    chains: list[str] = Field(min_length=1, max_length=62)
    pocket_residues: str = Field(min_length=3, max_length=4000)
    pocket_evidence: str = Field(min_length=8, max_length=2000)
    ph: float = Field(default=7.4, ge=4, le=10)
    padding: float = Field(default=5, ge=2, le=15)
    confirmed: bool = False


def source_payload(record: dict, find_manifest: Callable) -> tuple[dict, dict]:
    if not record or record.get('status') != 'succeeded':
        raise ValueError('请选择已成功完成的 OpenFold3 任务')
    manifest_path = find_manifest((record.get('result') or {}).get('agent'))
    if not manifest_path:
        raise ValueError('该任务没有 OpenFold3 结构产物')
    path = (ROOT/manifest_path).resolve()
    path.relative_to((ROOT/'results').resolve())
    manifest = json.loads(path.read_text())
    if manifest.get('provenance', {}).get('model') != 'OpenFold3' or manifest.get('structure',{}).get('format') != 'mmcif':
        raise ValueError('仅支持有 OpenFold3 来源记录的 mmCIF 产物')
    relative = manifest['structure']['url']
    if not isinstance(relative, str) or '://' in relative:
        raise ValueError('结构产物必须是本地相对路径')
    structure_path = (path.parent/relative).resolve()
    structure_path.relative_to(path.parent)
    raw = structure_path.read_bytes()
    if len(raw)>25*1024*1024:
        raise ValueError('结构超过 25 MB')
    digest = hashlib.sha256(raw).hexdigest()
    if not manifest['structure'].get('sha256') or digest != manifest['structure']['sha256']:
        raise ValueError('结构文件校验和缺失或不匹配')
    return {'source_run_id':record['id'], 'structure_base64':base64.b64encode(raw).decode()}, manifest


def resolve_prepared_receptor(payload: dict, store) -> dict | None:
    receptor = payload.get('receptor') or {}
    prep_id = receptor.get('preparation_id')
    if not prep_id:
        return None
    record = store.get(str(prep_id))
    if not record or record['mode'] != 'vina_receptor' or record['status'] != 'succeeded':
        raise ValueError('受体准备尚未成功完成')
    metadata = record.get('result') or {}
    path = (ROOT/'results/vina-receptors'/record['id']/'receptor.pdbqt').resolve()
    path.relative_to((ROOT/'results/vina-receptors').resolve())
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != metadata.get('receptor_sha256'):
        raise ValueError('准备后的受体文件校验失败')
    payload['receptor'] = {'name':f"OpenFold3-{metadata['source_run_id'][:8]}.pdbqt", 'content_base64':base64.b64encode(raw).decode()}
    # These fields come from the server preparation ledger, never browser claims.
    return {'preparation_id':record['id'], 'source_run_id':metadata['source_run_id'],
            'source_sha256':metadata['source_sha256'], 'receptor_sha256':metadata['receptor_sha256'],
            'options':metadata['options'], 'pocket_quality':metadata['pocket_quality'],
            'chain_mapping':metadata['chain_mapping'], 'box':metadata['box'],
            'protonation':metadata['protonation'], 'tools':metadata['tools'],
            'excluded_nonprotein':metadata['excluded_nonprotein'], 'scientific_limit':metadata['scientific_limit']}


def register_receptor_routes(app, store, executor, find_manifest):
    @app.get('/vina/structures/{source_run_id}')
    def inspect_structure(source_run_id: str):
        if store.get(source_run_id) is None:
            raise HTTPException(404, detail='结构任务不存在')
        try:
            payload, manifest = source_payload(store.get(source_run_id),find_manifest)
            inspection = call_receptor_worker('inspect',payload)
            return {'source_run_id':source_run_id, 'title':manifest.get('title','OpenFold3 预测结构'),
                    'metrics':manifest.get('metrics',{}), 'quality':manifest.get('quality',{}),
                    'provenance':manifest.get('provenance',{}), **inspection}
        except ValueError as exc:
            raise HTTPException(422,detail=str(exc)) from exc
        except Exception as exc:
            raise HTTPException(503,detail='结构检查服务暂不可用：'+str(exc)) from exc

    def execute(preparation_id, request):
        started = time.monotonic()
        store.update(preparation_id,status='running',current_stage='receptor_preparation',started_at=datetime.now(timezone.utc).isoformat())
        store.event(preparation_id,'stage',{'name':'检查预测结构、准备氢原子与受体 PDBQT'})
        try:
            payload, manifest = source_payload(store.get(request.source_run_id),find_manifest)
            metrics, quality = manifest.get('metrics',{}), manifest.get('quality',{})
            if quality.get('structure_review_allowed') is not True or not isinstance(metrics.get('avg_plddt'),(float,int)) or metrics['avg_plddt']<70 or not isinstance(metrics.get('ptm'),(float,int)) or metrics['ptm']<.5:
                raise ValueError('预测结构未通过结构审查质量门，请先检查质量或改用实验结构')
            payload.update(job_id=preparation_id,options=request.model_dump(exclude={'source_run_id'}))
            result = call_receptor_worker('prepare',payload)
            store.update(preparation_id,status='succeeded',current_stage='completed',completed_at=datetime.now(timezone.utc).isoformat(),runtime_seconds=round(time.monotonic()-started,3),result_json=json.dumps(result,ensure_ascii=False))
            store.event(preparation_id,'completed',{'name':'受体准备完成','source_run_id':request.source_run_id})
        except Exception as exc:
            store.update(preparation_id,status='failed',current_stage='failed',completed_at=datetime.now(timezone.utc).isoformat(),runtime_seconds=round(time.monotonic()-started,3),error=str(exc))
            store.event(preparation_id,'failed',{'name':'受体准备未通过','error':str(exc)})

    @app.post('/vina/receptors',status_code=202)
    def create_preparation(request: ReceptorPreparationRequest):
        if store.get(request.source_run_id) is None:
            raise HTTPException(404, detail='结构任务不存在')
        if not request.confirmed:
            raise HTTPException(422,detail='请先确认链、辅因子、结合位点证据和准备方案')
        try:
            source_payload(store.get(request.source_run_id),find_manifest)
        except (ValueError,FileNotFoundError) as exc:
            raise HTTPException(422,detail=str(exc)) from exc
        record = store.create('OpenFold3 → Vina 受体准备','vina_receptor',manual=request.model_dump(),persist_prompt=False)
        executor.submit(execute,record['run_id'],request)
        return record

    @app.get('/vina/receptors/{preparation_id}')
    def get_preparation(preparation_id: str):
        record = store.get(preparation_id)
        if not record or record.get('mode')!='vina_receptor':
            raise HTTPException(404,detail='受体准备任务不存在')
        result = record.get('result') or {}
        return {'preparation_id':preparation_id,'status':record['status'],'stage':record['current_stage'],
                'error':record.get('error'),'runtime_seconds':record.get('runtime_seconds'),'result':result}
