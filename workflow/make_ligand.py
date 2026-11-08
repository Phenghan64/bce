"""Create an independent starting conformer; never read a known ligand pose."""
import json
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

folder = Path(sys.argv[1])
case = json.loads((folder / 'case.json').read_text()) if (folder / 'case.json').exists() else {}
smiles = case.get('smiles', 'Cc1ccccc1')
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
params = AllChem.ETKDGv3()
params.randomSeed = 20261021
assert AllChem.EmbedMolecule(mol, params) == 0
status = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
assert status == 0
mol.SetProp('_Name', 'Independent SMILES conformer; ETKDGv3 seed 20261021; MMFF optimized')
with Chem.SDWriter(str(folder / 'ligand-input.sdf')) as writer:
    writer.write(mol)
(folder / 'ligand-input.json').write_text(json.dumps({
    'smiles': smiles, 'embedding_seed': 20261021,
    'embedding': 'ETKDGv3', 'minimization': 'MMFF', 'minimization_status': status,
    'heavy_atoms': mol.GetNumHeavyAtoms(), 'atoms_with_hydrogen': mol.GetNumAtoms(),
    'source': 'SMILES only; no crystal or predicted ligand coordinates used'
}, indent=2) + '\n')
