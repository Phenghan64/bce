"""Build a deterministic, unlabelled 1000-compound throughput pilot plus control."""
import hashlib
import json
import random
import shutil
import sys
import time
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
if (out / 'library.json').exists():
    raise SystemExit('Library already exists; use the existing library or a new directory.')
start = time.perf_counter()
source = Path(RDConfig.RDDataDir) / 'NCI/first_5K.smi'
(out / 'source').mkdir(exist_ok=True)
shutil.copy2(source, out / 'source/NCI-first_5K.smi')
control_smiles = json.loads((ROOT / 'results/closed-loop-1iep/case.json').read_text())['smiles']
control_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(control_smiles))
seen = {control_canonical}
eligible, excluded = [], Counter()
for line in source.read_text().splitlines():
    fields = line.split()
    if len(fields) < 2:
        excluded['invalid_line'] += 1
        continue
    mol = Chem.MolFromSmiles(fields[0])
    if mol is None:
        excluded['invalid_smiles'] += 1
        continue
    smi = Chem.MolToSmiles(mol)
    if smi in seen:
        excluded['duplicate'] += 1
        continue
    seen.add(smi)
    mw = Descriptors.MolWt(mol)
    heavy = mol.GetNumHeavyAtoms()
    rotors = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if (len(Chem.GetMolFrags(mol)) != 1 or not 200 <= mw <= 550 or
            not 15 <= heavy <= 40 or rotors > 10 or abs(Chem.GetFormalCharge(mol)) > 1 or
            any(a.GetSymbol() not in {'C','N','O','F','P','S','Cl','Br','I','H'} for a in mol.GetAtoms()) or
            any(str(info.specified) == 'Unspecified'
                for info in Chem.FindPotentialStereo(mol))):
        excluded['outside_declared_filters'] += 1
        continue
    eligible.append({'id': f'nci-example-{fields[1]}', 'source_id': fields[1],
        'source': 'RDKit NCI first_5K example library', 'role': 'unlabelled_library',
        'smiles': smi, 'molecular_weight': mw, 'heavy_atoms': heavy,
        'rdkit_rotatable_bonds': rotors, 'formal_charge': Chem.GetFormalCharge(mol)})
assert len(eligible) >= 1000
selected = random.Random(20261021).sample(eligible, 1000)
control = Chem.MolFromSmiles(control_smiles)
selected.append({'id':'imatinib-control', 'source':'Existing official Vina example microstate',
    'role':'positive_pose_control', 'smiles':control_smiles, 'molecular_weight':Descriptors.MolWt(control),
    'heavy_atoms':control.GetNumHeavyAtoms(), 'rdkit_rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(control),
    'formal_charge':Chem.GetFormalCharge(control)})
(out / 'ligands').mkdir(exist_ok=True)
records = []
for i, entry in enumerate(selected, 1):
    t = time.perf_counter()
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(entry['smiles']))
        params = AllChem.ETKDGv3(); params.randomSeed = 20261021; params.numThreads = 1
        embedded = AllChem.EmbedMolecule(mol, params)
        entry['embedding_random_coordinate_retry'] = embedded != 0
        if embedded != 0:
            params.useRandomCoords = True
            embedded = AllChem.EmbedMolecule(mol, params)
        if embedded != 0:
            raise ValueError('3D embedding failed')
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            raise ValueError('No complete MMFF parameters; not silently changing preparation method')
        converged = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        if converged == 1:
            converged = AllChem.MMFFOptimizeMolecule(mol, maxIters=5000)
        if converged != 0:
            raise ValueError(f'MMFF minimization did not converge: {converged}')
        mol.SetProp('_Name', entry['id'])
        sdf = out / 'ligands' / f'{entry["id"]}.sdf'
        with Chem.SDWriter(str(sdf)) as writer:
            writer.write(mol)
        setups = MoleculePreparation().prepare(mol)
        if len(setups) != 1:
            raise ValueError(f'Expected one microstate, received {len(setups)} setups')
        pdbqt, valid, error = PDBQTWriterLegacy.write_string(setups[0])
        if not valid:
            raise ValueError(error)
        dest = out / 'ligands' / f'{entry["id"]}.pdbqt'
        dest.write_text(pdbqt)
        entry.update(status='prepared', pdbqt=str(dest), sdf=str(sdf),
                     pdbqt_sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
    except Exception as exc:
        entry.update(status='failed', error=f'{type(exc).__name__}: {exc}')
    entry['preparation_seconds'] = time.perf_counter()-t
    records.append(entry)
    if i % 100 == 0 or entry['status']=='failed' or i==1001: print(f'{i}/1001 {entry["id"]}: {entry["status"]}',flush=True)
manifest = {'rdkit_version':rdBase.rdkitVersion, 'selection_seed':20261021, 'protocol_sha256':hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest(),
    'source_path':str(source), 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'eligible_unique_count':len(eligible), 'exclusions':dict(excluded),
    'selection':'1000 unique eligible records sampled without replacement; imatinib added as separate control',
    'filters':'Single component; MW 200–550; 15–40 heavy atoms; <=10 RDKit rotors; charge -1..1; organic elements; no unassigned potential atom/bond stereochemistry (FindPotentialStereo)',
    'preparation':'Source protonation/tautomer state retained, no pH enumeration; ETKDGv3 + MMFF; Meeko default Gasteiger charges',
    'limitations':'Unlabelled chemical library, not an ABL active/inactive benchmark. No enrichment or AUC can be inferred.',
    'preparation_wall_seconds':time.perf_counter()-start, 'records':records}
(out/'library.json').write_text(json.dumps(manifest, indent=2)+'\n')
assert next(r for r in records if r['id']=='imatinib-control')['status']=='prepared'
print(json.dumps({'prepared':sum(r['status']=='prepared' for r in records),
                  'failed':sum(r['status']=='failed' for r in records),
                  'seconds':manifest['preparation_wall_seconds']}), flush=True)
