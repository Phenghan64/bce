"""Prepare a protein at nominal pH 7.4, retaining observed heavy atoms."""
import json
import sys
import random
from pathlib import Path

import numpy as np
from pdbfixer import PDBFixer
from openmm import Platform, unit
from openmm.app import PDBFile

source, destination, report = map(Path, sys.argv[1:])
random.seed(42)
np.random.seed(42)
fixer = PDBFixer(filename=str(source), platform=Platform.getPlatformByName('Reference'))

def coordinates():
    pos = np.asarray(fixer.positions.value_in_unit(unit.angstrom))
    heavy = [a for a in fixer.topology.atoms() if a.element.symbol != 'H']
    result = {(a.residue.chain.id, a.residue.id, a.name): pos[a.index].copy() for a in heavy}
    assert len(result) == len(heavy), 'Duplicate residue/atom identifiers in receptor'
    return result

before = coordinates()
fixer.findMissingResidues()
skipped = {str(k): v for k, v in fixer.missingResidues.items()}
fixer.missingResidues = {}  # Do not invent missing terminal/loop residues.
fixer.findMissingAtoms()
missing = {f'{r.chain.id}:{r.id}': [a.name for a in atoms]
           for r, atoms in fixer.missingAtoms.items()}
terminal = {f'{r.chain.id}:{r.id}': atoms for r, atoms in fixer.missingTerminals.items()}
fixer.addMissingAtoms(seed=42)
fixer.addMissingHydrogens(pH=7.4)
after = coordinates()
shift = max(float(np.linalg.norm(after[k] - xyz)) for k, xyz in before.items())
assert shift < 0.002, f'Existing heavy atoms moved unexpectedly: {shift} A'
with destination.open('w') as handle:
    PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
report.write_text(json.dumps({
    'source': str(source), 'nominal_pH': 7.4, 'platform': 'Reference',
    'missing_residues_not_rebuilt': skipped, 'added_missing_heavy_atoms': missing,
    'added_terminal_atoms': terminal,
    'heavy_atoms_before': len(before), 'heavy_atoms_after': len(after),
    'hydrogens_after': sum(a.element.symbol == 'H' for a in fixer.topology.atoms()),
    'max_existing_heavy_atom_shift_angstrom': shift,
    'note': 'Standard protonation assignment; no full receptor minimization or protonation ensemble.'
}, indent=2) + '\n')
