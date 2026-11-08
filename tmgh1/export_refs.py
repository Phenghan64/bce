"""Audit experimental sources and export matched native-chain references."""
import json
from pathlib import Path
import numpy as np
from biotite.structure.io import pdb,pdbx
from rdkit import Chem
from rdkit.Chem import AllChem
root=Path(__file__).resolve().parents[1];out=root/'results/tmgh1-round1'
def fit(x,y):
 xm=x.mean(0);ym=y.mean(0);u,s,vt=np.linalg.svd((x-xm).T@(y-ym));r=u@np.diag([1,1,np.linalg.det(u@vt)])@vt;return r,ym-xm@r
refarray=None;report=[]
for code,lig in [('2J77','NOJ'),('2J7E','GI2')]:
 f=pdbx.CIFFile.read(out/'source'/f'{code}.cif');a=pdbx.get_structure(f,model=1,use_author_fields=True,extra_fields=['occupancy','b_factor']);p=a[(a.chain_id=='B')&~a.hetero&~np.isin(a.element,['H','D'])].copy();l=a[(a.chain_id=='B')&(a.res_name==lig)&~np.isin(a.element,['H','D'])].copy();ids=sorted(set(map(int,p.res_id)));assert len(ids)==442 and len(set(l.res_id))==1
 ca=p[p.atom_name=='CA'];assert len(ca)==442
 if refarray is None:refarray=ca.coord.copy();rotation=np.eye(3);translation=np.zeros(3)
 else:rotation,translation=fit(ca.coord,refarray)
 original_p=p.coord.copy();original_l=l.coord.copy();p.coord=p.coord@rotation+translation;l.coord=l.coord@rotation+translation
 d=out/'geometries'/f'experimental_{code}';d.mkdir(parents=True,exist_ok=True)
 p.chain_id[:]='A';p.chain_id[p.res_id>=308]='B'
 pf=pdb.PDBFile();pf.set_structure(p);pf.write(d/'observed.pdb')
 cc=pdbx.CIFFile.read(out/'source'/f'{lig}.cif').block;des=cc['pdbx_chem_comp_descriptor'];smiles=next(s for ty,pr,s in zip(*[des[k].as_array(str) for k in ['type','program','descriptor']]) if ty=='SMILES_CANONICAL' and pr=='OpenEye OEToolkits');template=Chem.MolFromSmiles(smiles)
 b=Chem.RWMol();ind={str(x.atom_name):b.AddAtom(Chem.Atom(str(x.element))) for x in l};cat=cc['chem_comp_bond']
 for x,y,order in zip(*[cat[k].as_array(str) for k in ['atom_id_1','atom_id_2','value_order']]):
  if x in ind and y in ind:b.AddBond(ind[x],ind[y],{'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE}[order])
 m=AllChem.AssignBondOrdersFromTemplate(template,b.GetMol());c=Chem.Conformer(len(l))
 for i,xyz in enumerate(l.coord):c.SetAtomPosition(i,tuple(map(float,xyz)))
 m.AddConformer(c);Chem.AssignStereochemistryFrom3D(m);Chem.SanitizeMol(m);assert Chem.MolToSmiles(m)==Chem.MolToSmiles(template)
 ld=out/'reference';ld.mkdir(exist_ok=True)
 with Chem.SDWriter(str(ld/f'{lig}.sdf')) as w:w.write(m)
 distances=np.linalg.norm(original_p[:,None]-original_l[None,:],axis=2).min(1);pocket=sorted(set(map(int,p.res_id[distances<=6])));gapdist={str(i):float(distances[p.res_id==i].min()) for i in [306,308]}
 waters=a[(a.chain_id=='B')&(a.res_name=='HOH')];wd=np.linalg.norm(waters.coord[:,None]-original_l[None,:],axis=2).min(1);calcium=a[(a.chain_id=='B')&(a.res_name=='CA')];cd=float(np.linalg.norm(calcium.coord[:,None]-original_l[None,:],axis=2).min());assert cd>20
 attrs={'pdb':code,'author_chain':'B','observed_residues':442,'native_sequence_residues':446,'construct_residues':468,'retained_native_positions':ids,'missing_native_positions':[1,2,307,446],'gap_endpoint_distances_A':gapdist,'ligand':lig,'ccd_smiles':Chem.MolToSmiles(template),'ligand_heavy_atoms':len(l),'ligand_occupancy_min':float(l.occupancy.min()),'ligand_occupancy_max':float(l.occupancy.max()),'ligand_B_mean':float(l.b_factor.mean()),'pocket_residues_6A':pocket,'pocket_water_count_within_4A':int(sum(wd<4)),'omitted_calcium_ligand_distance_A':cd,'rotation_to_2J77':rotation.tolist(),'translation_to_2J77':translation.tolist(),'crystal_global_CA_RMSD_to_2J77_A':float(np.sqrt(np.mean(np.sum((ca.coord@rotation+translation-refarray)**2,axis=1)))),'box_center':l.coord.mean(0).tolist(),'box_size':[24,24,24],'crystal_pH':7.0,'biological_assembly':'monomer, as deposited; A/B in asymmetric unit are separate copies'}
 (d/'metadata.json').write_text(json.dumps(attrs,indent=2)+'\n');report.append(attrs);print({k:v for k,v in attrs.items() if k not in ['retained_native_positions','rotation_to_2J77','translation_to_2J77']})
(out/'reference/source-audit.json').write_text(json.dumps(report,indent=2)+'\n')
