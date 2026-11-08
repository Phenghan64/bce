"""Assess every docked pose in a fixed receptor frame, including top-ranked failures."""
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

out = Path(sys.argv[1])
case = json.loads((out / 'case.json').read_text()) if (out / 'case.json').exists() else {}
manifest = json.loads((out / 'manifest.json').read_text())
ligand_tag = case.get('ligand_tag', 'toluene')
system = case.get('system', '7L39 / T4 lysozyme L99A / toluene')
alignment = json.loads((out / 'alignment.json').read_text())
validation = json.loads((out / 'validation.json').read_text())
rotation = np.asarray(alignment['rotation'])
pc = np.asarray(alignment['predicted_centroid'])
qc = np.asarray(alignment['reference_centroid'])
seeds = [42, 2026, 20261021]
template = Chem.MolFromSmiles(case.get('smiles', 'Cc1ccccc1'))
expected_smiles = Chem.MolToSmiles(template)

def load_molecules(path):
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=True))
    assert molecules and all(m is not None for m in molecules), path
    for mol in molecules:
        assert mol.GetNumHeavyAtoms() == template.GetNumHeavyAtoms()
        assert Chem.MolToSmiles(mol) == expected_smiles
        assert np.isfinite(mol.GetConformer().GetPositions()).all()
    return molecules

def aligned(mol):
    result = Chem.Mol(mol)
    coordinates = (mol.GetConformer().GetPositions() - pc) @ rotation + qc
    for i, xyz in enumerate(coordinates):
        result.GetConformer().SetAtomPosition(i, xyz)
    return result

def rmsd(a, b):
    # CalcRMS enumerates symmetry matches WITHOUT rotating/translating a ligand.
    return float(rdMolAlign.CalcRMS(a, b))

def read_heavy_pdb(path):
    coordinates, residues = [], []
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')) and line[76:78].strip() not in ('H', 'D'):
            coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            residues.append(f'{line[21]}:{line[22:26].strip()}:{line[17:20].strip()}')
    assert coordinates
    return np.asarray(coordinates), np.asarray(residues)

reference = load_molecules(out / f'reference-{ligand_tag}.sdf')[0]
prediction = load_molecules(out / f'predicted-{ligand_tag}.sdf')[0]
ref_xyz, ref_residues = read_heavy_pdb(out / 'reference-receptor.pdb')
ref_distances = np.linalg.norm(reference.GetConformer().GetPositions()[:, None] - ref_xyz[None], axis=2)
reference_contacts = set(ref_residues[np.any(ref_distances < 4.0, axis=0)])
report = {'execution_checks_passed': True, 'system': system,
    'openfold': validation,
    'protocol': {'docking_engine': 'AutoDock Vina 1.2.7', 'scoring': 'vina',
        'receptor': 'rigid', 'nominal_pH': 7.4, 'box_angstrom': alignment['box_size_angstrom'],
        'cpu': 6, 'exhaustiveness': 32, 'seeds': seeds,
        'pocket': 'OpenFold3 ligand centroid for predicted receptor; crystal ligand centroid for control',
        'ligand_start': 'Independent ETKDGv3 + MMFF conformer from specified SMILES',
        'ligand_smiles': expected_smiles,
        'rmsd': f'Symmetry-aware heavy-atom CalcRMS, no ligand fitting; fixed {validation["protein_ca_atoms_compared"]}-CA receptor alignment',
        'benchmark_top_pose_threshold_angstrom': 2.0,
        'contact_distance_angstrom': 4.0},
    'reference_contact_residues': sorted(reference_contacts), 'receptors': {},
    'limitations': [
        'Known benchmark; not an unbiased novel-target or large-complex validation.',
        'Prediction-guided local docking is not blind pocket discovery or independent proof of the pocket.',
        'Vina scores are approximate docking scores, not measured binding free energies.',
        'Preparation uses standard protonation, dry rigid receptors, and no protonation/conformational ensemble.',
        manifest['msa'] + '; timing of local stages excludes prior MSA acquisition.',
        'This is a computational validation loop; no experiment or model retraining occurred.']}
top_molecules = {}
all_predicted = []
for kind in ['predicted', 'reference']:
    folder = out / kind
    protein_xyz, residue_ids = read_heavy_pdb(folder / 'receptor-h.pdb')
    summaries = []
    tops = []
    for seed in seeds:
        pdbqt = (folder / f'seed-{seed}.pdbqt').read_text()
        energies = [float(x) for x in re.findall(r'^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)', pdbqt, re.M)]
        molecules = load_molecules(folder / f'seed-{seed}.sdf')
        assert len(molecules) == len(energies) == len(re.findall(r'^MODEL\s+', pdbqt, re.M))
        assert all(np.isfinite(energies)) and energies == sorted(energies)
        poses = []
        for rank, (molecule, energy) in enumerate(zip(molecules, energies), 1):
            coords = molecule.GetConformer().GetPositions()
            center = np.asarray(alignment[f'{kind}_box_center'])
            assert np.all(np.abs(coords - center) <= np.asarray(alignment['box_size_angstrom'])/2 + 0.1), 'Pose outside docking box'
            distances = np.linalg.norm(coords[:, None] - protein_xyz[None], axis=2)
            contacts = set(residue_ids[np.any(distances < 4.0, axis=0)])
            transformed = aligned(molecule) if kind == 'predicted' else molecule
            pose = {'rank': rank, 'vina_score_kcal_mol': energy,
                    'rmsd_to_crystal_angstrom': rmsd(transformed, reference),
                    'minimum_protein_ligand_heavy_distance_angstrom': float(distances.min()),
                    'heavy_atom_pairs_below_2_angstrom': int(np.sum(distances < 2)),
                    'contact_residues': sorted(contacts),
                    'reference_contact_recall': len(contacts & reference_contacts) / len(reference_contacts)}
            if kind == 'predicted':
                pose['rmsd_to_openfold_angstrom'] = rmsd(molecule, prediction)
                all_predicted.append((energy, seed, rank, molecule, pose))
            poses.append(pose)
        tops.append(molecules[0])
        summaries.append({'seed': seed, 'pose_count': len(poses), 'poses': poses,
            'top_pose_passes_2A': poses[0]['rmsd_to_crystal_angstrom'] <= 2.0,
            'best_rmsd_pose_rank': min(poses, key=lambda p: p['rmsd_to_crystal_angstrom'])['rank']})
    consistency = [{'seeds': [seeds[i], seeds[j]], 'top_pose_rmsd_angstrom': rmsd(tops[i], tops[j])}
                   for i, j in itertools.combinations(range(len(seeds)), 2)]
    report['receptors'][kind] = {'runs': summaries, 'top_pose_consistency': consistency,
        'all_top_poses_pass_2A': all(r['top_pose_passes_2A'] for r in summaries),
        'preparation': json.loads((folder / 'preparation.json').read_text())}

selected = min(all_predicted, key=lambda entry: (entry[0], entry[1], entry[2]))
energy, seed, rank, mol, pose = selected
with Chem.SDWriter(str(out / f'selected-docked-{ligand_tag}.sdf')) as writer:
    mol.SetProp('_Name', f'Vina predicted receptor; selected by score; seed {seed}; rank {rank}')
    writer.write(mol)
with Chem.SDWriter(str(out / 'comparison-in-crystal-frame.sdf')) as writer:
    for name, molecule in [('Crystal', reference), ('OpenFold3', aligned(prediction)),
                           ('Vina_selected_by_score', aligned(mol))]:
        molecule.SetProp('_Name', name)
        writer.write(molecule)
report['selection'] = {'rule': 'Lowest Vina score across predicted-receptor runs; no crystal RMSD used in selection',
    'seed': seed, 'rank': rank, **pose,
    'output': str(out / f'selected-docked-{ligand_tag}.sdf')}
report['benchmark_passed'] = all(r['all_top_poses_pass_2A'] for r in report['receptors'].values())
report['feedback'] = ('Both receptor arms recover the crystal pose in all three top-ranked runs.'
    if report['benchmark_passed'] else
    'Execution succeeded but pose recovery did not pass in every top-ranked run. Inspect receptor/pocket geometry and sampling before using ranked poses.')
(out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
lines = [f'OpenFold3 → Meeko → AutoDock Vina：{system} 完整流程验收', '',
    f'执行检查：通过；所有种子首位姿势 RMSD ≤ 2 Å：{report["benchmark_passed"]}',
    f'OpenFold3：Cα RMSD {validation["protein_ca_rmsd_angstrom"]:.3f} Å；配体 RMSD {validation["ligand_heavy_atom_rmsd_angstrom"]:.3f} Å',
    '配体 RMSD 按蛋白质对齐并考虑分子对称性计算，不单独拟合配体。', '',
    '受体 / 随机种子 / 首位 Vina 分数(kcal/mol) / 首位晶体 RMSD(Å) / 最佳 RMSD 姿势排名']
for kind, receptor in report['receptors'].items():
    for run in receptor['runs']:
        top = run['poses'][0]
        lines.append(f'{kind} / {run["seed"]} / {top["vina_score_kcal_mol"]:.3f} / {top["rmsd_to_crystal_angstrom"]:.3f} / {run["best_rmsd_pose_rank"]}')
lines += ['', '选择结果：按照预测受体对接分数选取，不使用实验 RMSD 排名。',
          f'种子 {seed}，排名 {rank}，分数 {energy:.3f} kcal/mol，晶体 RMSD {pose["rmsd_to_crystal_angstrom"]:.3f} Å。',
          '', report['feedback'], '', '适用边界：'] + report['limitations']
(out / 'report.txt').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
