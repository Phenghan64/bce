"""Export all observed residues of one chain; split genuine gaps without inventing loops."""
import json
from pathlib import Path
import numpy as np
from biotite.structure.io import pdb,pdbx
from rdkit import Chem
from rdkit.Chem import AllChem
root=Path(__file__).resolve().parents[1];out=root/'results/large-receptors'
for code,chain,lig,name in [('2BXD','A','RWF','albumin'),('4JSX','B','17G','mtor')]:
 d=out/name/'receptor';d.mkdir(parents=True,exist_ok=True)
 f=pdbx.CIFFile.read(out/'source'/f'{code}.cif');a=pdbx.get_structure(f,model=1,use_author_fields=True)
 p=a[(a.chain_id==chain)&~a.hetero&~np.isin(a.element,['H','D'])].copy()
 l=a[(a.chain_id==chain)&(a.res_name==lig)&~np.isin(a.element,['H','D'])]
 ligand_ids=np.unique(l.res_id);assert len(ligand_ids)==1
 ids=sorted(set(map(int,p.res_id)));gaps=[(i,j) for i,j in zip(ids,ids[1:]) if j-i>1]
 segments=[];start=ids[0]
 for i,j in gaps+[(ids[-1],None)]:
  segment_id=chr(65+len(segments));sel=(p.res_id>=start)&(p.res_id<=i);p.chain_id[sel]=segment_id
  segments.append({'chain':segment_id,'start':start,'end':i,'residues':len(set(p.res_id[sel]))})
  start=j
 pf=pdb.PDBFile();pf.set_structure(p);pf.write(d/'observed.pdb')
 ccd=pdbx.CIFFile.read(out/'source'/f'{lig}.cif').block
 desc=ccd['pdbx_chem_comp_descriptor'];rows=list(zip(*[desc[k].as_array(str) for k in ['type','program','descriptor']]))
 smi=next(x[2] for x in rows if x[0]=='SMILES_CANONICAL' and x[1]=='OpenEye OEToolkits')
 template=Chem.MolFromSmiles(smi);b=Chem.RWMol();indices={str(x.atom_name):b.AddAtom(Chem.Atom(str(x.element))) for x in l}
 cat=ccd['chem_comp_bond']
 for x,y,order in zip(*[cat[k].as_array(str) for k in ['atom_id_1','atom_id_2','value_order']]):
  if x in indices and y in indices:b.AddBond(indices[x],indices[y],{'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE}[order])
 m=AllChem.AssignBondOrdersFromTemplate(template,b.GetMol());c=Chem.Conformer(len(l))
 for i,xyz in enumerate(l.coord):c.SetAtomPosition(i,tuple(map(float,xyz)))
 m.AddConformer(c);Chem.AssignStereochemistryFrom3D(m);Chem.SanitizeMol(m)
 assert Chem.MolToSmiles(m)==Chem.MolToSmiles(template),(Chem.MolToSmiles(m),smi)
 with Chem.SDWriter(str(d/'crystal-control.sdf')) as w:w.write(m)
 gap_info=[]
 for i,j in gaps:
  neigh=p[np.isin(p.res_id,[i,j])]
  gap_info.append({'after':i,'before':j,'omitted_residues':j-i-1,'endpoint_min_ligand_distance_A':float(np.linalg.norm(neigh.coord[:,None]-l.coord[None,:],axis=2).min())})
 seqs=f.block['entity_poly']; seq=next(seq for chains,seq in zip(seqs['pdbx_strand_id'].as_array(str),seqs['pdbx_seq_one_letter_code_can'].as_array(str)) if chain in chains.split(','))
 meta={'pdb':code,'name':name,'author_chain':chain,'observed_residues':len(ids),'construct_sequence_residues':len(seq.replace('\n','')),'observed_heavy_atoms':len(p),'segments':segments,'gaps':gap_info,'control_ligand':lig,'control_smiles':Chem.MolToSmiles(template),'box_center':l.coord.mean(axis=0).tolist(),'box_size':[30,30,30], 'scope':'Single selected protein chain, all observed residues retained. Waters, other proteins and crystal ligands omitted; rigid receptor. Missing loops not built; gap ends terminalized. Performance experiment, not full biological assembly.'}
 (d/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta))
