"""Validate the official T4 lysozyme/toluene example against PDB 7L39."""
import json
import argparse
from pathlib import Path

import numpy as np
from biotite.structure.io import pdb, pdbx
from rdkit import Chem
from rdkit.Chem import rdMolAlign

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser()
parser.add_argument('--output-dir', type=Path, default=root / 'results/openfold3-7l39')
parser.add_argument('--prediction-dir', type=Path)
parser.add_argument('--reference-dir', type=Path, default=root / 'results/openfold3-7l39/reference')
args = parser.parse_args()
folder = args.output_dir
folder.mkdir(parents=True, exist_ok=True)
prediction_dir = args.prediction_dir or folder / 'with-msa'
seed = prediction_dir / "prediction/pdb_7L39/seed_42"
prefix = "pdb_7L39_seed_42_sample_1"
pred_cif = pdbx.CIFFile.read(seed / f"{prefix}_model.cif")
ref_cif = pdbx.CIFFile.read(args.reference_dir / "7L39.cif")
ccd = pdbx.CIFFile.read(args.reference_dir / "MBN.cif")
pred = pdbx.get_structure(pred_cif, model=1)
ref = pdbx.get_structure(ref_cif, model=1)
assert np.isfinite(pred.coord).all()
protein = pred[(pred.chain_id == "A") & ~pred.hetero]
pred_ca = protein[protein.atom_name == "CA"]
ref_ca = ref[(ref.chain_id == "A") & ~ref.hetero & (ref.atom_name == "CA")]
assert len(pred_ca) == 164 and len(ref_ca) == 162
index = {int(r): i for i, r in enumerate(pred_ca.res_id)}
matching = pred_ca[[index[int(r)] for r in ref_ca.res_id]]
assert np.array_equal(matching.res_name, ref_ca.res_name)
p = matching.coord.astype(float)
q = ref_ca.coord.astype(float)
pc, qc = p.mean(axis=0), q.mean(axis=0)
u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
rotation = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
aligned_ca = (p - pc) @ rotation + qc
protein_rmsd = np.sqrt(np.mean(np.sum((aligned_ca - q) ** 2, axis=1)))


def ligand_molecule(atoms, bond_category, residue_name):
    builder = Chem.RWMol()
    atom_indices = {}
    for atom in atoms:
        atom_indices[str(atom.atom_name)] = builder.AddAtom(Chem.Atom(str(atom.element)))
    bond_types = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
                  "TRIP": Chem.BondType.TRIPLE}
    for component, a, b, order, aromatic in zip(
        *[bond_category[key].as_array(str) for key in
          ["comp_id", "atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"]]
    ):
        if component != residue_name or a not in atom_indices or b not in atom_indices:
            continue
        kind = Chem.BondType.AROMATIC if aromatic == "Y" else bond_types[order]
        builder.AddBond(atom_indices[a], atom_indices[b], kind)
    mol = builder.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(len(atoms))
    conf.Set3D(True)
    for i, xyz in enumerate(atoms.coord):
        conf.SetAtomPosition(i, tuple(float(x) for x in xyz))
    mol.AddConformer(conf)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(Chem.MolFromSmiles("Cc1ccccc1"))
    return mol


pred_lig = pred[(pred.chain_id == "X") & ~np.isin(pred.element, ["H", "D"])]
ref_lig = ref[(ref.res_name == "MBN") & ~np.isin(ref.element, ["H", "D"])]
assert len(pred_lig) == len(ref_lig) == 7
pred_mol = ligand_molecule(pred_lig, pred_cif.block["chem_comp_bond"], "LIG0")
ref_mol = ligand_molecule(ref_lig, ccd.block["chem_comp_bond"], "MBN")
aligned_mol = Chem.Mol(pred_mol)
aligned_lig_coords = (pred_lig.coord.astype(float) - pc) @ rotation + qc
for i, xyz in enumerate(aligned_lig_coords):
    aligned_mol.GetConformer().SetAtomPosition(i, tuple(float(x) for x in xyz))
# Preserve the receptor alignment: CalcRMS accounts for molecular symmetry but
# does not independently rotate/translate the ligand.
ligand_rmsd = rdMolAlign.CalcRMS(aligned_mol, ref_mol)
centroid_distance = np.linalg.norm(aligned_lig_coords.mean(axis=0) - ref_lig.coord.mean(axis=0))
distances = np.linalg.norm(pred_lig.coord[:, None, :] - protein.coord[None, :, :], axis=2)
confidence = json.loads((seed / f"{prefix}_confidences_aggregated.json").read_text())
benchmark = json.loads((prediction_dir / "benchmark.json").read_text())
config = json.loads((prediction_dir / "prediction/experiment_config.json").read_text())
assert benchmark["cli_completed"] and config["pl_trainer_args"]["accelerator"] == "mps"
report = {
    "system": "T4 lysozyme L99A / toluene, official OpenFold3 example (PDB 7L39)",
    "protein_residues_predicted": 164,
    "protein_ca_atoms_compared": 162,
    "ligand_heavy_atoms": 7,
    "reference_url": "https://files.rcsb.org/download/7L39.cif",
    "ligand_chemistry_url": "https://files.rcsb.org/ligands/download/MBN.cif",
    "model_runtime_seconds": json.loads((seed / "timing.json").read_text())["runtime_s"],
    "cli_seconds": benchmark["elapsed_seconds"],
    "process_peak_rss_gib": benchmark["process_peak_rss_bytes"] / 1024**3,
    "sampled_peak_mps_driver_gib": benchmark["sampled_peak_mps_driver_bytes"] / 1024**3,
    "system_swap_change_mib": (benchmark["system_swap_after_bytes"] - benchmark["system_swap_before_bytes"]) / 1024**2,
    "avg_plddt": confidence["avg_plddt"],
    "iptm": confidence["iptm"],
    "ptm": confidence["ptm"],
    "protein_ca_rmsd_angstrom": float(protein_rmsd),
    "ligand_heavy_atom_rmsd_angstrom": float(ligand_rmsd),
    "ligand_centroid_error_angstrom": float(centroid_distance),
    "minimum_protein_ligand_heavy_atom_distance_angstrom": float(distances.min()),
    "model_reports_clash": confidence["has_clash"],
    "method": "Kabsch alignment on 162 experimentally observed C-alpha atoms; ligand graph from CIF bond tables; symmetry-corrected ligand RMSD in the protein-aligned frame, without separate ligand fitting",
    "settings": "Cached ColabFold MSA; no templates or pocket constraints; MPS; official default OpenBind checkpoint; seed 42; one diffusion sample; float32; low_mem",
    "limitations": "Known small benchmark; not independent generalization evidence, binding-affinity validation, or a large-complex capacity test. MSA generation was remote; inference was local. RSS and MPS memory overlap.",
}
(folder / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
with Chem.SDWriter(str(folder / "predicted-toluene.sdf")) as writer:
    pred_mol.SetProp("_Name", "OpenFold3 predicted toluene; original prediction coordinates")
    writer.write(pred_mol)
pdb_file = pdb.PDBFile()
pdb_file.set_structure(protein)
pdb_file.write(folder / "predicted-receptor.pdb")
assert next(iter(Chem.SDMolSupplier(str(folder / "predicted-toluene.sdf")))) is not None
with Chem.SDWriter(str(folder / 'reference-toluene.sdf')) as writer:
    writer.write(ref_mol)
ref_protein = ref[(ref.chain_id == 'A') & ~ref.hetero]
assert len(ref_protein[ref_protein.atom_name == 'CA']) == 162
pdb_file = pdb.PDBFile()
pdb_file.set_structure(ref_protein)
pdb_file.write(folder / 'reference-receptor.pdb')
(folder / 'alignment.json').write_text(json.dumps({
    'predicted_centroid': pc.tolist(), 'reference_centroid': qc.tolist(),
    'rotation': rotation.tolist(),
    'convention': '(predicted_coordinates - predicted_centroid) @ rotation + reference_centroid',
    'box_size_angstrom': [20.0, 20.0, 20.0],
    'predicted_box_center': pred_lig.coord.mean(axis=0).tolist(),
    'reference_box_center': ref_lig.coord.mean(axis=0).tolist(),
}, indent=2) + '\n')
print(json.dumps(report, indent=2))
