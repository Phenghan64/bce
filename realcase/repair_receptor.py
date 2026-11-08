"""Rebuild only the documented internal ASP gap; preserve observed heavy atoms."""
import json,time,random
from pathlib import Path
import numpy as np
from pdbfixer import PDBFixer
from openmm import Platform,unit,Context,VerletIntegrator,LocalEnergyMinimizer
from openmm.app import PDBFile,ForceField,NoCutoff

out=Path(__file__).resolve().parents[1]/'results/abl-real-case/receptor'
start=time.perf_counter();random.seed(42);np.random.seed(42)
fixer=PDBFixer(filename=str(out/'observed.pdb'),platform=Platform.getPlatformByName('Reference'))
def coords():
    xyz=np.asarray(fixer.positions.value_in_unit(unit.angstrom))
    return {(a.residue.id,a.name):xyz[a.index].copy() for a in fixer.topology.atoms() if a.element.symbol!='H'}
before=coords(); residues=list(fixer.topology.residues())
assert residues[42].id=='48' and residues[43].id=='50'
fixer.missingResidues={(0,43):['ASP']}
fixer.findMissingAtoms()
missing={r.id:[a.name for a in atoms] for r,atoms in fixer.missingAtoms.items()}
fixer.addMissingAtoms(seed=42)
# Initialize modeled carboxylates with planar geometry before local relaxation.
# Some generated side chains otherwise remain in a folded local minimum.
positions=np.asarray(fixer.positions.value_in_unit(unit.angstrom))
for r in fixer.topology.residues():
    if r.name not in ['ASP','GLU']:continue
    cn,attached,o1,o2=('CG','CB','OD1','OD2') if r.name=='ASP' else ('CD','CG','OE1','OE2')
    atoms={a.name:a.index for a in r.atoms()}
    if any((r.id,n) in before for n in [o1,o2]):continue
    origin=positions[atoms[cn]];u=positions[atoms[attached]]-origin;u/=np.linalg.norm(u)
    v=positions[atoms[o1]]-origin;v-=u*np.dot(u,v);v/=np.linalg.norm(v)
    theta=np.deg2rad(117)
    for name,sign in [(o1,1),(o2,-1)]:positions[atoms[name]]=origin+1.25*(np.cos(theta)*u+sign*np.sin(theta)*v)
fixer.positions=positions*unit.angstrom
fixer.addMissingHydrogens(pH=7.4)
# Relax only modeled atoms and hydrogens; the experimental heavy atoms are fixed.
system=ForceField('amber14-all.xml').createSystem(fixer.topology,nonbondedMethod=NoCutoff,constraints=None)
for a in fixer.topology.atoms():
    if a.element.symbol!='H' and (a.residue.id,a.name) in before:system.setParticleMass(a.index,0)
integrator=VerletIntegrator(.001)
context=Context(system,integrator,Platform.getPlatformByName('Reference'))
context.setPositions(fixer.positions)
energy_before=float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
LocalEnergyMinimizer.minimize(context,tolerance=1,maxIterations=15000)
energy_after=float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
fixer.positions=context.getState(getPositions=True).getPositions()
after=coords();shift=max(float(np.linalg.norm(after[k]-v)) for k,v in before.items())
assert shift<.002
residues=list(fixer.topology.residues()); xyz=np.asarray(fixer.positions.value_in_unit(unit.angstrom))
assert len(residues)==266
dist=[]
for a,b in zip(residues,residues[1:]):
    ca=next(x for x in a.atoms() if x.name=='C');nb=next(x for x in b.atoms() if x.name=='N')
    d=float(np.linalg.norm(xyz[ca.index]-xyz[nb.index]));dist.append(d)
assert min(dist)>1.1 and max(dist)<1.8, (min(dist),max(dist))
with (out/'receptor-h.pdb').open('w') as w:PDBFile.writeFile(fixer.topology,fixer.positions,w,keepIds=True)
report={'seconds':time.perf_counter()-start,'energy_before_kJ_mol':energy_before,'energy_after_kJ_mol':energy_after,'rebuilt_residue':'B:49 ASP (author B:276)',
 'added_missing_atoms':missing,'max_observed_heavy_atom_shift_A':shift,'peptide_CN_min_A':min(dist),'peptide_CN_max_A':max(dist),
 'residue_count':len(residues),'pH':7.4,'note':'Missing terminal residues omitted; internal gap modeled, not experimentally observed. Amber14 minimization of added atoms/H only, observed heavy atoms fixed.'}
(out/'repair.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
