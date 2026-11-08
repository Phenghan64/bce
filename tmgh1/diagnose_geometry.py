"""Post-hoc geometry diagnosis only; no parameter changes or selection of favorable runs."""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
out=Path(__file__).resolve().parents[1]/'results/tmgh1-round1'
def atoms(path):
 a=[]
 for l in path.read_text().splitlines():
  if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip() not in ['H','D']:a.append({'residue':int(l[22:26]),'resname':l[17:20].strip(),'atom':l[12:16].strip(),'element':l[76:78].strip(),'xyz':np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])})
 return a
radius={'C':1.7,'N':1.55,'O':1.52,'S':1.8};ref=atoms(out/'geometries/experimental_2J77/observed.pdb');refca={a['residue']:a['xyz'] for a in ref if a['atom']=='CA'};rows=[]
for geom in ['experimental_2J77','experimental_2J7E','predicted_42','predicted_2026','predicted_20261021']:
 a=atoms(out/'receptors'/f'{geom}_default'/'receptor-h.pdb');xyz=np.array([x['xyz'] for x in a]);ca={x['residue']:x['xyz'] for x in a if x['atom']=='CA'};ids=sorted(set(refca)&set(ca));absolute_rms=float(np.sqrt(np.mean(np.sum((np.array([ca[i] for i in ids])-np.array([refca[i] for i in ids]))**2,axis=1))))
 assert absolute_rms<1,'Unexpected coordinate-frame mismatch'
 for ligand,code in [('NOJ','2J77'),('GI2','2J7E')]:
  m=Chem.SDMolSupplier(str(out/'reference'/f'{ligand}.sdf'),removeHs=True)[0];lxyz=m.GetConformer().GetPositions();dist=np.linalg.norm(xyz[:,None]-lxyz[None,:],axis=2);meta=json.loads((out/'geometries'/f'experimental_{code}'/'metadata.json').read_text());assert np.all(np.abs(lxyz-np.array(meta['box_center']))<=12),'Crystal ligand outside search grid';overlap=np.array([radius.get(x['element'],1.7) for x in a])[:,None]+np.array([radius.get(x.GetSymbol(),1.7) for x in m.GetAtoms()])[None,:]-dist
  indices=np.dstack(np.unravel_index(np.argsort(overlap.ravel())[-8:][::-1],overlap.shape))[0];top=[]
  for i,j in indices:top.append({'residue':a[i]['residue'],'resname':a[i]['resname'],'protein_atom':a[i]['atom'],'ligand_atom_index':int(j),'ligand_element':m.GetAtomWithIdx(int(j)).GetSymbol(),'distance_A':float(dist[i,j]),'vdW_overlap_A':float(overlap[i,j])})
  rows.append({'geometry':geom,'ligand':ligand,'CA_RMSD_in_actual_docking_frame_to_2J77_A':absolute_rms,'crystal_pose_inside_box':True,'crystal_ligand_min_protein_heavy_distance_A':float(dist.min()),'heavy_pairs_under_2A':int(np.sum(dist<2)),'pairs_vdW_overlap_above_1A':int(np.sum(overlap>1)),'largest_overlaps':top})
(out/'posthoc-geometry-diagnosis.json').write_text(json.dumps({'scope':'Post-hoc diagnosis of observed docking failures; no input changes or new docking jobs. Geometric overlaps are not binding energy calculations.','rows':rows},indent=2)+'\n')
for r in rows:print(r['geometry'],r['ligand'],'absolute_CA_RMSD',round(r['CA_RMSD_in_actual_docking_frame_to_2J77_A'],3),'min distance',round(r['crystal_ligand_min_protein_heavy_distance_A'],3),'severe pairs',r['heavy_pairs_under_2A'],'largest',r['largest_overlaps'][:1])
