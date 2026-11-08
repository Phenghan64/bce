"""Run in the isolated OpenMM/PDBFixer environment; do not rebuild missing loops."""
import json
import sys
from pathlib import Path

import numpy as np
from openmm import Platform, unit
from openmm.app import PDBFile
from pdbfixer import PDBFixer

source, destination, report = map(Path, sys.argv[1:4])
ph = float(sys.argv[4])
fixer = PDBFixer(filename=str(source), platform=Platform.getPlatformByName('Reference'))

def heavy_coordinates():
    positions = np.asarray(fixer.positions.value_in_unit(unit.angstrom))
    atoms = [atom for atom in fixer.topology.atoms() if atom.element.symbol != 'H']
    result = {(atom.residue.chain.id, atom.residue.id, atom.name): positions[atom.index].copy() for atom in atoms}
    if len(result) != len(atoms):
        raise ValueError('Duplicate atom or residue identifiers')
    return result

before = heavy_coordinates()
fixer.findMissingResidues()
skipped = {str(key): value for key, value in fixer.missingResidues.items()}
if any(key[1] not in (0, len(list(list(fixer.topology.chains())[key[0]].residues()))) for key in fixer.missingResidues):
    raise ValueError('Internal missing residues require manual receptor repair')
fixer.missingResidues = {}
fixer.findMissingAtoms()
missing = {f'{res.chain.id}:{res.id}': [a.name for a in atoms] for res, atoms in fixer.missingAtoms.items() if atoms}
if missing:
    raise ValueError('Missing heavy atoms require manual repair: ' + json.dumps(missing))
terminal = {f'{res.chain.id}:{res.id}': atoms for res, atoms in fixer.missingTerminals.items()}
fixer.addMissingAtoms(seed=42)
fixer.addMissingHydrogens(pH=ph)
after = heavy_coordinates()
shift = max(float(np.linalg.norm(after[key] - xyz)) for key, xyz in before.items())
if shift > 0.002:
    raise ValueError(f'Existing heavy atoms moved by {shift:.4f} angstrom')
with destination.open('w') as handle:
    PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
report.write_text(json.dumps({
    'nominal_pH': ph, 'platform': 'Reference', 'missing_residues_not_rebuilt': skipped,
    'added_terminal_atoms': terminal, 'heavy_atoms_before': len(before), 'heavy_atoms_after': len(after),
    'hydrogens_after': sum(atom.element.symbol == 'H' for atom in fixer.topology.atoms()),
    'max_existing_heavy_atom_shift_angstrom': shift,
    'note': 'Standard pH-based protonation, no protonation ensemble or receptor minimization. Catalytic protonation requires specialist review.'
}, indent=2) + '\n')
