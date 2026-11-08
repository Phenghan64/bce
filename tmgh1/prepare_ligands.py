"""Prepare four predeclared ligand states from CCD chemistry, without pose initialization."""
import sys,json,time,hashlib
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolDescriptors,rdMolAlign
from meeko import MoleculePreparation,PDBQTWriterLegacy
root=Path(__file__).resolve().parents[1];out=root/'results/tmgh1-round1';folder=out/'ligands';folder.mkdir(exist_ok=True);start=time.perf_counter();records=[]
for lig in ['NOJ','GI2']:
 crystal=Chem.SDMolSupplier(str(out/'reference'/f'{lig}.sdf'),removeHs=True)[0]
 for state in ['cation','neutral']:
  m=Chem.MolFromSmiles(Chem.MolToSmiles(crystal))
  if lig=='NOJ' and state=='cation':
   n=next(a for a in m.GetAtoms() if a.GetSymbol()=='N');n.SetFormalCharge(1);n.SetNumExplicitHs(2)
  if lig=='GI2' and state=='neutral':
   for a in m.GetAtoms():
    if a.GetSymbol()=='N':a.SetFormalCharge(0);a.SetNumExplicitHs(0)
  Chem.SanitizeMol(m);smiles=Chem.MolToSmiles(m);assert Chem.GetFormalCharge(m)==(1 if state=='cation' else 0)
  reference=AllChem.AssignBondOrdersFromTemplate(m,crystal);Chem.SanitizeMol(reference)
  refpath=folder/f'{lig}-{state}-reference.sdf'
  with Chem.SDWriter(str(refpath)) as w:w.write(reference)
  mol=Chem.AddHs(m);params=AllChem.ETKDGv3();params.randomSeed=20261021;params.numThreads=1;assert AllChem.EmbedMolecule(mol,params)==0;assert AllChem.MMFFOptimizeMolecule(mol,maxIters=5000)==0
  ident=f'{lig}-{state}';sdf=folder/f'{ident}.sdf'
  with Chem.SDWriter(str(sdf)) as w:w.write(mol)
  setup=MoleculePreparation().prepare(mol);assert len(setup)==1;qt,ok,err=PDBQTWriterLegacy.write_string(setup[0]);assert ok,err;dest=folder/f'{ident}.pdbqt';dest.write_text(qt)
  heavy=Chem.RemoveHs(mol);fitted=Chem.Mol(heavy);rms=float(rdMolAlign.GetBestRMS(fitted,reference));ringatoms=[a.GetIdx() for a in fitted.GetAtoms() if a.IsInRing()];match=reference.GetSubstructMatch(heavy);assert len(match)==heavy.GetNumAtoms();import numpy as np
  x=heavy.GetConformer().GetPositions()[ringatoms];y=reference.GetConformer().GetPositions()[[match[i] for i in ringatoms]];rx=np.linalg.norm(x[:,None]-x[None,:],axis=2);ry=np.linalg.norm(y[:,None]-y[None,:],axis=2);ringdist=float(np.sqrt(np.mean((rx-ry)**2)))
  records.append({'id':ident,'ligand':lig,'state':state,'primary':state=='cation','smiles':smiles,'formal_charge':Chem.GetFormalCharge(m),'heavy_atoms':m.GetNumHeavyAtoms(),'rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(m),'status':'prepared','pdbqt':str(dest),'sdf':str(sdf),'reference':str(refpath),'pdbqt_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'independent_input_best_fitted_RMSD_A':rms,'rigid_ring_pair_distance_RMSD_A':ringdist,'note':'Fitted RMSD used only to diagnose initial conformation mismatch; the unfitted independent conformer is the docking input.'})
(out/'library.json').write_text(json.dumps({'records':records,'preparation_seconds':time.perf_counter()-start},indent=2)+'\n');print(json.dumps(records,indent=2))
