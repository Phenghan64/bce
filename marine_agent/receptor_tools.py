"""Explicit, bounded preparation of an OpenFold protein as a rigid Vina receptor."""
from __future__ import annotations

import base64
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import string
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
AA = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
MAX_CIF = 25 * 1024 * 1024


def read_source(encoded):
    import gemmi
    raw = base64.b64decode(encoded, validate=True)
    if not raw or len(raw) > MAX_CIF:
        raise ValueError('结构为空或超过 25 MB')
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(raw.decode()).sole_block())
    if len(structure) != 1:
        raise ValueError('自动准备要求单个结构模型，请手动选择模型')
    records, extras = [], []
    seen = set()
    for chain in structure[0]:
        for residue in chain:
            key = f'{chain.name}:{residue.seqid}'
            if key in seen:
                raise ValueError('存在重复残基标识，请手动修复')
            seen.add(key)
            atoms = [a for a in residue if a.element.name not in ('H', 'D')]
            if not atoms:
                continue
            if any(not math.isfinite(v) for a in atoms for v in (a.pos.x, a.pos.y, a.pos.z)):
                raise ValueError('结构含非有限坐标')
            if any(a.altloc not in ('\x00', ' ') for a in atoms):
                raise ValueError('存在多重侧链构象，需要先选择 altloc')
            if residue.name not in AA:
                extras.append({'id': key, 'name': residue.name, 'atoms': len(atoms),
                               'protein_like': any(a.name == 'CA' for a in atoms) and any(a.name == 'N' for a in atoms)})
                continue
            ca = next((a for a in atoms if a.name == 'CA'), None)
            if ca is None:
                raise ValueError(f'{key} 缺少 CA 原子')
            confidence = sum(a.b_iso for a in atoms) / len(atoms)
            records.append({'id': key, 'chain': chain.name, 'number': residue.seqid.num,
                            'name': residue.name, 'plddt': round(confidence, 2)})
    if not records or len(records) > 10000:
        raise ValueError('未检测到标准蛋白残基，或结构超过当前受体准备容量')
    chains = []
    for name in dict.fromkeys(r['chain'] for r in records):
        rows = [r for r in records if r['chain'] == name]
        chains.append({'id': name, 'residues': len(rows), 'avg_plddt': round(sum(r['plddt'] for r in rows)/len(rows), 2)})
    return raw, structure, {'chains': chains, 'residues': records, 'nonprotein': extras,
                            'source_sha256': hashlib.sha256(raw).hexdigest(),
                            'confidence_source': 'OpenFold atom_site.B_iso_or_equiv averaged per residue'}


def inspect_source(payload):
    return read_source(payload['structure_base64'])[2]


def pocket_keys(spec, records):
    chosen = set()
    available = {r['id'] for r in records}
    for token in re.split(r'[\s,，;；]+', str(spec).strip()):
        if not token:
            continue
        if token in available:
            chosen.add(token)
            continue
        match = re.fullmatch(r'([^:\s,]+):(-?\d+)-(-?\d+)', token)
        if not match:
            raise ValueError(f'口袋残基不存在或格式不正确：{token}，例如 A:50,A:60-65')
        chain, start, end = match.group(1), int(match.group(2)), int(match.group(3))
        if end < start or end-start > 1000:
            raise ValueError('口袋残基范围无效')
        requested = {f'{chain}:{i}' for i in range(start, end+1)}
        if not requested.issubset(available):
            raise ValueError('口袋范围含不存在的残基')
        chosen.update(requested)
    if not chosen or len(chosen) > 1000:
        raise ValueError('请指定 1 到 1000 个有证据支持的口袋残基')
    return chosen


def prepare_receptor(payload, output):
    import gemmi
    started = time.monotonic()
    raw, structure, inspection = read_source(payload['structure_base64'])
    options = payload.get('options') or {}
    if options.get('confirmed') is not True:
        raise ValueError('请先确认组装状态、辅因子及准备方案')
    chains = options.get('chains')
    if not isinstance(chains, list) or not chains or len(set(chains)) != len(chains) or len(chains) > 62:
        raise ValueError('请选择有效且不重复的蛋白链')
    available_chains = {r['chain'] for r in inspection['residues']}
    if not set(chains).issubset(available_chains):
        raise ValueError('所选蛋白链不存在')
    if any(x['protein_like'] for x in inspection['nonprotein'] if x['id'].split(':')[0] in chains):
        raise ValueError('所选链存在非标准氨基酸，请手动准备，系统不会静默删除')
    ph, padding = float(options.get('ph', 7.4)), float(options.get('padding', 5))
    if not 4 <= ph <= 10 or not 2 <= padding <= 15:
        raise ValueError('pH 应在 4–10，口袋边距应在 2–15 Å')
    rows = [r for r in inspection['residues'] if r['chain'] in chains]
    keys = pocket_keys(options.get('pocket_residues', ''), rows)
    confidences = [r['plddt'] for r in rows if r['id'] in keys]
    avg, minimum = sum(confidences)/len(confidences), min(confidences)
    if avg < 70 or minimum < 50:
        raise ValueError(f'口袋启发式质量门未通过：平均 pLDDT {avg:.1f}，最低 {minimum:.1f}；请审查结构或改用实验结构')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    (out/'source.cif').write_bytes(raw)
    selected = gemmi.Structure()
    selected.add_model(gemmi.Model('1'))
    mappings, pocket_xyz = [], []
    for index, original_id in enumerate(chains):
        original = structure[0].find_chain(original_id)
        chain = gemmi.Chain(string.ascii_uppercase[index] if index < 26 else (string.ascii_lowercase+string.digits)[index-26])
        previous = None
        for residue in original:
            if residue.name not in AA:
                continue
            key = f'{original_id}:{residue.seqid}'
            if previous is not None:
                c = next((a for a in previous if a.name == 'C'), None)
                n = next((a for a in residue if a.name == 'N'), None)
                if c is None or n is None or c.pos.dist(n.pos) > 2.0:
                    raise ValueError(f'{key} 附近存在主链断裂，需要手动准备')
            previous = residue
            clone = residue.clone()
            clone.seqid = gemmi.SeqId(len(chain)+1, ' ')
            clone.het_flag = 'A'
            for idx in range(len(clone)-1,-1,-1):
                if clone[idx].element.name in ('H','D'):
                    del clone[idx]
            if key in keys:
                pocket_xyz.extend((a.pos.x,a.pos.y,a.pos.z) for a in clone)
            mappings.append({'source': key, 'prepared': f'{chain.name}:{clone.seqid}'})
            chain.add_residue(clone)
        selected[0].add_chain(chain)
    lo = [min(p[i] for p in pocket_xyz)-padding for i in range(3)]
    hi = [max(p[i] for p in pocket_xyz)+padding for i in range(3)]
    box = {'center': [round((a+b)/2,3) for a,b in zip(lo,hi)], 'size': [round(b-a,3) for a,b in zip(lo,hi)]}
    if any(v>200 for v in box['size']):
        raise ValueError('口袋搜索盒超过 200 Å，请检查残基选择')
    selected.write_pdb(str(out/'selected.pdb'))
    prep_python = Path(os.environ.get('RECEPTOR_PREP_PYTHON', str(ROOT/'.tools/prep-env/bin/python')))
    if not prep_python.is_file():
        raise RuntimeError('本机受体质子化环境不可用')
    commands = [
        [str(prep_python), str(ROOT/'marine_agent/receptor_protonate.py'), str(out/'selected.pdb'), str(out/'protonated.pdb'), str(out/'protonation.json'), str(ph)],
        [sys.executable, str(Path(sys.executable).parent/'mk_prepare_receptor.py'), '--read_pdb', str(out/'protonated.pdb'), '-o', str(out/'receptor'), '-p', '-j'],
    ]
    for number, command in enumerate(commands):
        with (out/f'preparation-{number}.log').open('w') as log:
            run = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=180, cwd=ROOT)
        if run.returncode:
            message = (out/f'preparation-{number}.log').read_text()[-1800:]
            raise ValueError('受体准备失败，未删除不匹配残基：'+message)
    pdbqt = (out/'receptor.pdbqt').read_text()
    atom_rows = [line for line in pdbqt.splitlines() if line.startswith(('ATOM  ','HETATM'))]
    if not atom_rows or any(not math.isfinite(float(line[70:76])) for line in atom_rows):
        raise ValueError('PDBQT 原子或电荷检查失败')
    prepared_residues = {(line[21], line[22:26].strip()) for line in atom_rows}
    if len(prepared_residues) != len(rows):
        raise ValueError('PDBQT 残基数与所选蛋白不一致')
    result = {'status': 'succeeded', 'source_run_id': payload.get('source_run_id'),
              'source_sha256': inspection['source_sha256'], 'receptor_sha256': hashlib.sha256(pdbqt.encode()).hexdigest(),
              'options': options, 'box': box, 'residues': len(rows), 'chain_mapping': mappings,
              'excluded_nonprotein': inspection['nonprotein'],
              'pocket_quality': {'mean_plddt': round(avg,2), 'min_plddt': minimum, 'residues': sorted(keys),
                                 'thresholds': {'mean_min':70,'residue_min':50}, 'status':'review_allowed'},
              'protonation': json.loads((out/'protonation.json').read_text()),
              'tools': {name:importlib.metadata.version(name) for name in ('meeko','rdkit','gemmi')},
              'commands': commands, 'elapsed_seconds': round(time.monotonic()-started,3),
              'scientific_limit': '启发式质量门仅允许继续计算；所选链不是经验证的生物组装。标准质子化不确定催化残基状态；不证明结合或活性。'}
    (out/'preparation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    return result
