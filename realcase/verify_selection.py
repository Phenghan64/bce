"""Independently reconstruct the panel from source matrix and frozen sampling rules."""
import json,random,hashlib
from pathlib import Path
import numpy as np
from rdkit import Chem

out=Path(__file__).resolve().parents[1]/'results/abl-real-case';source=out/'source'
matrix=np.loadtxt(source/'davis-kd-nm.txt');ids=(source/'drug-cids.txt').read_text().split()
targets=(source/'target-names.txt').read_text().splitlines();smiles=json.loads((source/'ligands_iso.json').read_text())
assert matrix.shape==(len(ids),len(targets))==(68,442)
column=targets.index('ABL1-nonphosphorylated')
assert matrix[ids.index('5291'),column]==1.1
groups={'strong':[],'weak_or_not_detected':[]}
for cid,kd in zip(ids,matrix[:,column]):
    mol=Chem.MolFromSmiles(smiles[cid])
    if cid=='5291' or len(Chem.GetMolFrags(mol))!=1 or any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(mol)):continue
    if 100<kd<10000:continue
    label='strong' if kd<=100 else 'weak_or_not_detected'
    groups[label].append({'cid':cid,'kd_nm':float(kd),'label':label,'source_smiles':smiles[cid]})
rng=random.Random(20261021);panel=[]
for group in groups.values():panel.extend(rng.sample(sorted(group,key=lambda x:int(x['cid'])),8))
protocol=json.loads((out/'protocol.json').read_text());assert panel==protocol['panel']
for record in json.loads((out/'provenance.json').read_text())['files']:
    assert hashlib.sha256((source/record['file']).read_bytes()).hexdigest()==record['sha256']
print('Verified source hashes, 68x442 row/column mapping, imatinib Kd1.1nM and all16 frozen sampled labels.')
