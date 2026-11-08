"""Export the preselected 2HYY B receptor and its experimental imatinib pose."""
import json
from pathlib import Path
import numpy as np
from biotite.structure.io import pdb, pdbx
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
out = ROOT/'results/abl-real-case'
f = pdbx.CIFFile.read(out/'source/2HYY.cif')
a = pdbx.get_structure(f, model=1, use_author_fields=False)
p = a[(a.chain_id=='B') & ~a.hetero & ~np.isin(a.element,['H','D'])]
l = a[(a.chain_id=='F') & ~np.isin(a.element,['H','D'])]
assert set(l.res_name)=={'STI'}
seq = str(f.block['entity_poly']['pdbx_seq_one_letter_code_can'].as_array()[0]).replace('\n','')
assert seq in json.loads((out/'source/proteins.json').read_text())['ABL1']
file = pdb.PDBFile(); file.set_structure(p); file.write(out/'receptor/observed.pdb')
template_smiles = json.loads((ROOT/'results/closed-loop-1iep/case.json').read_text())['smiles']
template = Chem.MolFromSmiles(template_smiles)
ccd_path = ROOT/'results/closed-loop-1iep/reference/STI.cif'
(out/'source/STI.cif').write_bytes(ccd_path.read_bytes())
cat = pdbx.CIFFile.read(ccd_path).block['chem_comp_bond']
b = Chem.RWMol(); indices = {str(x.atom_name):b.AddAtom(Chem.Atom(str(x.element))) for x in l}
for x,y,order,aromatic in zip(*[cat[k].as_array(str) for k in ['atom_id_1','atom_id_2','value_order','pdbx_aromatic_flag']]):
    if x in indices and y in indices:
        b.AddBond(indices[x],indices[y],Chem.BondType.AROMATIC if aromatic=='Y' else {'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE}[order])
m = AllChem.AssignBondOrdersFromTemplate(template,b.GetMol());Chem.SanitizeMol(m)
c = Chem.Conformer(len(l))
for i,xyz in enumerate(l.coord):c.SetAtomPosition(i,tuple(map(float,xyz)))
m.AddConformer(c)
assert Chem.MolToSmiles(m)==Chem.MolToSmiles(template)
with Chem.SDWriter(str(out/'receptor/crystal-imatinib.sdf')) as w:w.write(m)
meta = {'pdb':'2HYY','protein_chain':'B','ligand_label_chain':'F','observed_residues':len(set(p.res_id)),
    'sequence':seq,'control_smiles':template_smiles,'box_center':l.coord.mean(axis=0).tolist(),'box_size':[30,30,30],
    'internal_missing_residue':{'label_id':49,'author_id':276,'name':'ASP','insert_before_residue_index':43},
    'gap_neighbor_min_ligand_distance_A':float(np.linalg.norm(p[np.isin(p.res_id,[48,50])].coord[:,None]-l.coord[None],axis=2).min()),
    'sequence_matches_davis_ABL1':True,'terminal_missing_residues_not_rebuilt':[1,2,3,4,5,272,273]}
(out/'receptor/metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2))
