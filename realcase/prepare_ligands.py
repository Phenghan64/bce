"""Prepare the frozen real-case panel without consulting docking outcomes."""
import json,sys,time,hashlib
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
from meeko import MoleculePreparation,PDBQTWriterLegacy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.tools/realcase-deps'))
from dimorphite_dl import protonate_smiles
out=ROOT/'results/abl-real-case';start=time.perf_counter()
assert not (out/'library.json').exists(), 'Preserve existing prepared library'
protocol=json.loads((out/'protocol.json').read_text());meta=json.loads((out/'receptor/metadata.json').read_text())
(out/'ligands').mkdir(exist_ok=True)
entries=[];enumeration=[]
for compound in protocol['panel']:
    states=sorted(set(Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in protonate_smiles(
        compound['source_smiles'],ph_min=7.4,ph_max=7.4,precision=0,max_variants=128)))
    enumeration.append({'cid':compound['cid'],'states':states,'count':len(states)})
    if len(states)!=1:
        entries.append({'id':'cid-'+compound['cid'],'cid':compound['cid'],'status':'failed','error':'Expected exactly one nominal-pH state'});continue
    for i,s in enumerate(states):
        entries.append({'id':f'cid-{compound["cid"]}-state-{i:02}','cid':compound['cid'],'smiles':s,'role':'panel'})
entries.append({'id':'imatinib-control','cid':'5291','smiles':meta['control_smiles'],'role':'calibration'})
for entry in entries:
    if entry.get('status')=='failed':continue
    t=time.perf_counter()
    try:
        m=Chem.AddHs(Chem.MolFromSmiles(entry['smiles']));p=AllChem.ETKDGv3();p.randomSeed=20261021;p.numThreads=1
        if AllChem.EmbedMolecule(m,p)!=0:
            p.useRandomCoords=True
            if AllChem.EmbedMolecule(m,p)!=0:raise ValueError('Embedding failed')
        if not AllChem.MMFFHasAllMoleculeParams(m):raise ValueError('Missing MMFF parameters')
        c=AllChem.MMFFOptimizeMolecule(m,maxIters=1000)
        if c==1:c=AllChem.MMFFOptimizeMolecule(m,maxIters=5000)
        if c!=0:raise ValueError(f'MMFF convergence {c}')
        m.SetProp('_Name',entry['id']);sdf=out/'ligands'/f'{entry["id"]}.sdf'
        with Chem.SDWriter(str(sdf)) as w:w.write(m)
        setups=MoleculePreparation().prepare(m)
        if len(setups)!=1:raise ValueError('Unexpected Meeko setups')
        txt,valid,error=PDBQTWriterLegacy.write_string(setups[0])
        if not valid:raise ValueError(error)
        dest=sdf.with_suffix('.pdbqt');dest.write_text(txt)
        entry.update(status='prepared',sdf=str(sdf),pdbqt=str(dest),formal_charge=Chem.GetFormalCharge(m),heavy_atoms=m.GetNumHeavyAtoms())
    except Exception as exc:entry.update(status='failed',error=str(exc))
    entry['seconds']=time.perf_counter()-t
    print(entry['id'],entry['status'],entry.get('error',''),flush=True)
library={'records':entries,'enumeration':enumeration,'preparation_seconds':time.perf_counter()-start,
 'rdkit_version':rdBase.rdkitVersion,'protocol_sha256':hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest()}
(out/'library.json').write_text(json.dumps(library,indent=2)+'\n')
print('Prepared',sum(x['status']=='prepared' for x in entries),'of',len(entries),'in',library['preparation_seconds'],flush=True)
