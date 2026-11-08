"""Versioned, hash-verified ligand references resolved only by the API."""
from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
LIBRARY_ROOT = ROOT/'results/ligand-libraries'
VERSION = 'pubchem20-20261030-v1'
CATALOG_HASHES = {VERSION:'e5af005cfb4aaeec8a7e11f84052b45fa39e0ce4759e352d69cb1891f6d49c6d'}


class LibraryUnavailable(RuntimeError):
    pass


def _verified_file(root: Path, relative: str, expected: str) -> bytes:
    try:
        path = (root/relative).resolve()
        path.relative_to(root.resolve())
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError('hash mismatch')
        return raw
    except (OSError, ValueError) as exc:
        raise LibraryUnavailable('配体库文件校验失败，请稍后重试；仍可上传已准备配体。') from exc


def load_catalog(version: str = VERSION) -> dict[str, Any]:
    if version not in CATALOG_HASHES:
        raise ValueError('未知的配体库版本')
    raw = _verified_file(LIBRARY_ROOT/version,'catalog.json',CATALOG_HASHES[version])
    try:
        return json.loads(raw)
    except (ValueError, TypeError) as exc:
        raise LibraryUnavailable('配体库目录不可用') from exc


def public_catalog() -> dict[str, Any]:
    catalog = load_catalog()
    # No content or local filesystem paths are exposed to the browser.
    public = {k:v for k,v in catalog.items() if k != 'compounds'}
    public['compounds'] = [
        {**{k:v for k,v in c.items() if k != 'states'},'states':[
            {k:s[k] for k in ('id','display_name','formal_charge','stereo_label','smiles','chemistry_review_required','pdbqt_sha256')}
            for s in c['states']]} for c in catalog['compounds']]
    return public


def resolve_library_ligands(payload: dict[str, Any]) -> None:
    """Replace references with verified bytes; never trust browser provenance."""
    payload.pop('ligand_library', None)
    inputs = payload.get('ligands')
    if not isinstance(inputs,list):
        raise ValueError('配体格式无效')
    expanded, metadata, seen = [], {}, set()
    catalogs = {}
    for item in inputs:
        if not isinstance(item,dict):
            raise ValueError('配体格式无效')
        if 'library_compound_id' not in item and 'library_version' not in item:
            if 'ligand_library' in item or str(item.get('id','')).startswith('library-'):
                raise ValueError('上传配体不能声明配体库身份')
            expanded.append(item)
            continue
        if set(item) != {'library_compound_id','library_version'}:
            raise ValueError('配体库引用不能同时包含上传文件或自定义状态')
        version, compound_id = item['library_version'], item['library_compound_id']
        if not isinstance(version,str) or not isinstance(compound_id,str):
            raise ValueError('配体库引用格式无效')
        key = (version,compound_id)
        if key in seen:
            raise ValueError('配体库选择重复')
        seen.add(key)
        if version not in catalogs:
            catalogs[version] = load_catalog(version)
        catalog = catalogs[version]
        compound = next((c for c in catalog['compounds'] if c['id'] == compound_id),None)
        if not compound:
            raise ValueError('未知的配体库化合物')
        for state in compound['states']:
            raw = _verified_file(LIBRARY_ROOT/version,state['pdbqt'],state['pdbqt_sha256'])
            _verified_file(LIBRARY_ROOT/version,state['sdf'],state['sdf_sha256'])
            state_id = 'library-'+state['id']
            expanded.append({'id':state_id,'name':state['id']+'.pdbqt','content_base64':base64.b64encode(raw).decode('ascii')})
            metadata[state_id] = {
                **{k:compound[k] for k in ('cid','name_zh','name_en','pubchem_url','source_inchikey')},
                **{k:state[k] for k in ('display_name','smiles','source_smiles','formal_charge','stereo_label','chemistry_review_required','pdbqt_sha256','sdf_sha256')},
                'compound_id':compound_id, 'state_id':state['id'],'version':version,
                'catalog_sha256':CATALOG_HASHES[version],'versions':catalog['versions'],
                'preparation':catalog['preparation'],'state_policy':catalog['state_policy'],
                'sources':catalog['sources'],'protocol_sha256':catalog['protocol_sha256']}
    if not 1 <= len(expanded) <= 100:
        raise ValueError('展开后配体任务数量必须在 1 到 100 之间')
    payload['ligands'] = expanded
    if metadata:
        payload['ligand_library'] = metadata
