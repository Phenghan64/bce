"""Relax modeled atoms/H only; freeze every experimentally observed heavy atom."""
import sys,time,json,shutil
from pathlib import Path
import numpy as np
from openmm import Platform,Context,VerletIntegrator,LocalEnergyMinimizer,unit
from openmm.app import PDBFile,ForceField,NoCutoff
p=Path(sys.argv[1]);t=time.perf_counter();observed=PDBFile(str(p/'observed.pdb'));keys={(a.residue.chain.id,a.residue.id,a.name) for a in observed.topology.atoms()}
shutil.copy2(p/'receptor-h.pdb',p/'before-relaxation.pdb');f=PDBFile(str(p/'before-relaxation.pdb'));before=np.asarray(f.positions.value_in_unit(unit.angstrom));system=ForceField('amber14-all.xml').createSystem(f.topology,nonbondedMethod=NoCutoff,constraints=None);frozen=[]
for a in f.topology.atoms():
 if (a.residue.chain.id,a.residue.id,a.name) in keys:system.setParticleMass(a.index,0);frozen.append(a.index)
assert len(frozen)==len(keys)
integ=VerletIntegrator(.001);ctx=Context(system,integ,Platform.getPlatformByName('CPU'),{'Threads':'2'});ctx.setPositions(f.positions)
e0=float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole));LocalEnergyMinimizer.minimize(ctx,tolerance=5,maxIterations=10000);state=ctx.getState(getPositions=True,getEnergy=True);after=np.asarray(state.getPositions().value_in_unit(unit.angstrom));shift=float(np.linalg.norm(after[frozen]-before[frozen],axis=1).max());assert shift<.002
with (p/'receptor-h.pdb').open('w') as w:PDBFile.writeFile(f.topology,state.getPositions(),w,keepIds=True)
report={'seconds':time.perf_counter()-t,'energy_before_kJmol':e0,'energy_after_kJmol':float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)),'observed_heavy_atoms_frozen':len(frozen),'max_observed_shift_A':shift,'method':'Amber14, CPU 2 threads, only added atoms and H move, tolerance 5, max 10000 iterations'};(p/'relaxation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
