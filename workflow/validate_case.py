"""Export and validate a protein/ligand prediction against one experimental chain."""
import json
import sys
from pathlib import Path

import numpy as np
from biotite.structure.io import pdb, pdbx
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

folder = Path(sys.argv[1]).resolve()
case = json.loads((folder / 'case.json').read_text())
query = json.loads((folder / 'query.json').read_text())['queries'][case['query_name']]
sequence = query['chains'][0]['sequence']
template = Chem.MolFromSmiles(case['smiles'])
seed = folder / 'openfold/prediction' / case['query_name'] / 'seed_42'
prefix = f'{case["query_name"]}_seed_42_sample_1'
pred_cif = pdbx.CIFFile.read(seed / f'{prefix}_model.cif')
ref_cif = pdbx.CIFFile.read(folder / 'reference' / f'{case["pdb_id"]}.cif')
ccd = pdbx.CIFFile.read(folder / 'reference' / f'{case["ccd"]}.cif')
pred = pdbx.get_structure(pred_cif, model=1)
# label_seq_id maps experimental residues onto the full input entity sequence.
ref_label = pdbx.get_structure(ref_cif, model=1, use_author_fields=False)
ref_auth = pdbx.get_structure(ref_cif, model=1)
protein = pred[(pred.chain_id == 'A') & ~pred.hetero & ~np.isin(pred.element, ['H', 'D'])]
ref_protein = ref_label[(ref_label.chain_id == 'A') & ~ref_label.hetero & ~np.isin(ref_label.element, ['H', 'D'])]
pred_ca = protein[protein.atom_name == 'CA']
ref_ca = ref_protein[ref_protein.atom_name == 'CA']
assert len(pred_ca) == len(sequence)
index = {int(r): i for i, r in enumerate(pred_ca.res_id)}
matching = pred_ca[[index[int(r)] for r in ref_ca.res_id]]
assert np.array_equal(matching.res_name, ref_ca.res_name)
p, q = matching.coord.astype(float), ref_ca.coord.astype(float)
pc, qc = p.mean(axis=0), q.mean(axis=0)
u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
protein_rmsd = float(np.sqrt(np.mean(np.sum(((p-pc) @ rotation+qc-q)**2, axis=1))))

def molecule(atoms, category, resname):
    builder = Chem.RWMol()
    indices = {str(a.atom_name): builder.AddAtom(Chem.Atom(str(a.element))) for a in atoms}
    assert len(indices) == len(atoms)
    kinds = {'SING': Chem.BondType.SINGLE, 'DOUB': Chem.BondType.DOUBLE, 'TRIP': Chem.BondType.TRIPLE}
    for component, a, b, order, aromatic in zip(*[category[k].as_array(str) for k in
            ['comp_id', 'atom_id_1', 'atom_id_2', 'value_order', 'pdbx_aromatic_flag']]):
        if component == resname and a in indices and b in indices:
            builder.AddBond(indices[a], indices[b], Chem.BondType.AROMATIC if aromatic == 'Y' else kinds[order])
    mol = builder.GetMol()
    # Assign the same explicitly chosen charge state to all three sources.
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(len(atoms)); conf.Set3D(True)
    for i, xyz in enumerate(atoms.coord):
        conf.SetAtomPosition(i, tuple(float(x) for x in xyz))
    mol.AddConformer(conf)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(template)
    return mol

pred_lig = pred[(pred.chain_id == 'X') & ~np.isin(pred.element, ['H', 'D'])]
ref_lig = ref_auth[(ref_auth.chain_id == 'A') & (ref_auth.res_name == case['ccd']) & ~np.isin(ref_auth.element, ['H', 'D'])]
assert len(pred_lig) == len(ref_lig) == template.GetNumHeavyAtoms()
pred_mol = molecule(pred_lig, pred_cif.block['chem_comp_bond'], str(pred_lig.res_name[0]))
ref_mol = molecule(ref_lig, ccd.block['chem_comp_bond'], case['ccd'])
aligned = Chem.Mol(pred_mol)
aligned_xyz = (pred_lig.coord.astype(float)-pc) @ rotation+qc
for i, xyz in enumerate(aligned_xyz):
    aligned.GetConformer().SetAtomPosition(i, xyz)
confidence = json.loads((seed / f'{prefix}_confidences_aggregated.json').read_text())
benchmark = json.loads((folder / 'openfold/benchmark.json').read_text())
manifest = json.loads((folder / 'manifest.json').read_text())
assert benchmark['cli_completed']
assert np.isfinite(pred.coord).all()
distances = np.linalg.norm(pred_lig.coord[:, None]-protein.coord[None], axis=2)
report = {'system': case['system'], 'protein_residues_predicted': len(sequence),
    'protein_ca_atoms_compared': len(ref_ca), 'ligand_heavy_atoms': len(pred_lig),
    'model_runtime_seconds': json.loads((seed/'timing.json').read_text())['runtime_s'],
    'cli_seconds': benchmark['elapsed_seconds'],
    'process_peak_rss_gib': benchmark['process_peak_rss_bytes']/1024**3,
    'sampled_peak_mps_driver_gib': benchmark['sampled_peak_mps_driver_bytes']/1024**3,
    'system_swap_change_mib': (benchmark['system_swap_after_bytes']-benchmark['system_swap_before_bytes'])/1024**2,
    'avg_plddt': confidence['avg_plddt'], 'iptm': confidence['iptm'], 'ptm': confidence['ptm'],
    'protein_ca_rmsd_angstrom': protein_rmsd,
    'ligand_heavy_atom_rmsd_angstrom': float(rdMolAlign.CalcRMS(aligned, ref_mol)),
    'ligand_centroid_error_angstrom': float(np.linalg.norm(aligned_xyz.mean(axis=0)-ref_lig.coord.mean(axis=0))),
    'minimum_protein_ligand_heavy_atom_distance_angstrom': float(distances.min()),
    'model_reports_clash': confidence['has_clash'],
    'residue_numbering': 'Full entity sequence positions (label_seq_id), not author residue numbers',
    'method': f'Kabsch on all {len(ref_ca)} observed CA; symmetry-aware ligand RMSD without independent ligand fitting',
    'settings': f'One sample, seed 42; {manifest["msa"]}; no templates/pocket constraints; official OpenBind; MPS float32 low_mem',
    'limitations': 'Known benchmark, chosen ligand microstate, rigid docking; no experimental affinity validation.'}
(folder/'validation.json').write_text(json.dumps(report, indent=2)+'\n')
for kind, atoms, mol in [('predicted', protein, pred_mol), ('reference', ref_protein, ref_mol)]:
    f = pdb.PDBFile(); f.set_structure(atoms); f.write(folder/f'{kind}-receptor.pdb')
    with Chem.SDWriter(str(folder/f'{kind}-{case["ligand_tag"]}.sdf')) as writer:
        writer.write(mol)
alignment = {'predicted_centroid': pc.tolist(), 'reference_centroid': qc.tolist(), 'rotation': rotation.tolist(),
    'box_size_angstrom': case['box_size'], 'predicted_box_center': pred_lig.coord.mean(axis=0).tolist(),
    'reference_box_center': ref_lig.coord.mean(axis=0).tolist(),
    'convention': '(predicted_coordinates - predicted_centroid) @ rotation + reference_centroid'}
(folder/'alignment.json').write_text(json.dumps(alignment, indent=2)+'\n')
print(json.dumps(report, indent=2))
