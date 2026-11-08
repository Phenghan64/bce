"""Freeze ten ligands by flexibility, without consulting docking scores."""
import json,hashlib,shutil,time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem,Descriptors,rdMolDescriptors
from meeko import MoleculePreparation,PDBQTWriterLegacy
root=Path(__file__).resolve().parents[1];out=root/'results/large-receptors';start=time.perf_counter()
source=root/'results/vina-throughput-1000/library.json';old=json.loads(source.read_text());selected=[]
for rotors in range(10):
 eligible=[r for r in old['records'] if r['status']=='prepared' and r['role']=='unlabelled_library' and r['rdkit_rotatable_bonds']==rotors]
 r=sorted(eligible,key=lambda x:(-x['heavy_atoms'],x['id']))[0].copy();selected.append(r)
(out/'panel-selection.json').write_text(json.dumps({'rule':'One molecule for each RDKit rotor count 0 through 9; choose largest heavy-atom count, tie break ID; no docking scores consulted. Reuse frozen ligand preparations.','source_manifest_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'ids':[r['id'] for r in selected]},indent=2)+'\n')
for name in ['albumin','mtor']:
 d=out/name;folder=d/'ligands';folder.mkdir(exist_ok=True);records=[]
 for orig in selected:
  r=orig.copy()
  for key in ['pdbqt','sdf']:
   src=Path(r[key]);dest=folder/src.name;shutil.copy2(src,dest);r[key]=str(dest)
  assert hashlib.sha256(Path(r['pdbqt']).read_bytes()).hexdigest()==r['pdbqt_sha256'];records.append(r)
 meta=json.loads((d/'receptor/metadata.json').read_text());smi=meta['control_smiles'];m=Chem.AddHs(Chem.MolFromSmiles(smi));p=AllChem.ETKDGv3();p.randomSeed=20261021;p.numThreads=1
 assert AllChem.EmbedMolecule(m,p)==0;assert AllChem.MMFFOptimizeMolecule(m,maxIters=5000)==0
 control='crystal-control';sdf=folder/f'{control}.sdf'
 with Chem.SDWriter(str(sdf)) as w:w.write(m)
 setup=MoleculePreparation().prepare(m);assert len(setup)==1;qt,ok,err=PDBQTWriterLegacy.write_string(setup[0]);assert ok,err
 dest=folder/f'{control}.pdbqt';dest.write_text(qt)
 records.append({'id':control,'role':'positive_pose_control','smiles':smi,'status':'prepared','pdbqt':str(dest),'sdf':str(sdf),'heavy_atoms':m.GetNumHeavyAtoms(),'rdkit_rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(m)),'pdbqt_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
 protocol={'receptor':meta,'panel_count':10,'panel_seeds':[42,2026,20261021],'control_count':1,'control_seeds':[42,2026,20261021],'exhaustiveness':32,'workers':3,'cpu_per_job':2,'box_A':[30,30,30],'purpose':'Hardware/performance benchmark; common unlabelled ligands. Control failures are reported, never suppressed; performance runs proceed even if pose calibration fails. No activity/enrichment conclusions. No OpenFold prediction.','ligand_preparation':'Ten pre-existing ETKDGv3/MMFF conformers reused unchanged; own control newly embedded, not initialized from crystal pose. Source microstates retained; no protonation ensemble.'}
 (d/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(d/'library.json').write_text(json.dumps({'records':records,'panel_preparation':'Reused prepared ligand files; fresh control embedding','setup_seconds_so_far':time.perf_counter()-start},indent=2)+'\n')
print([(r['id'],r['heavy_atoms'],r['rdkit_rotatable_bonds']) for r in selected])
