"""Add missing atoms/H to already separated observed peptide segments."""
import sys,json,time,random
from pathlib import Path
import numpy as np
from pdbfixer import PDBFixer
from openmm import Platform,unit
from openmm.app import PDBFile
out=Path(sys.argv[1]);start=time.perf_counter();random.seed(42);np.random.seed(42)
f=PDBFixer(filename=str(out/'observed.pdb'),platform=Platform.getPlatformByName('CPU'))
def coords():
 xyz=np.asarray(f.positions.value_in_unit(unit.angstrom))
 return {(a.residue.chain.id,a.residue.id,a.name):xyz[a.index].copy() for a in f.topology.atoms() if a.element.symbol!='H'}
before=coords();f.findMissingResidues();skipped={str(k):v for k,v in f.missingResidues.items()};f.missingResidues={}
f.findMissingAtoms();missing={f'{r.chain.id}:{r.id}':[a.name for a in aa] for r,aa in f.missingAtoms.items()};f.addMissingAtoms(seed=42)
# Fix initial geometry of completely modeled carboxylate oxygens.
pos=np.asarray(f.positions.value_in_unit(unit.angstrom))
for r in f.topology.residues():
 if r.name not in ['ASP','GLU']:continue
 cn,attached,o1,o2=('CG','CB','OD1','OD2') if r.name=='ASP' else ('CD','CG','OE1','OE2');atoms={a.name:a.index for a in r.atoms()}
 if any((r.chain.id,r.id,n) in before for n in [o1,o2]):continue
 origin=pos[atoms[cn]];u=pos[atoms[attached]]-origin;u/=np.linalg.norm(u);v=pos[atoms[o1]]-origin;v-=u*np.dot(u,v);v/=np.linalg.norm(v)
 for name,sign in [(o1,1),(o2,-1)]:pos[atoms[name]]=origin+1.25*(np.cos(np.deg2rad(117))*u+sign*np.sin(np.deg2rad(117))*v)
f.positions=pos*unit.angstrom;f.addMissingHydrogens(pH=7.4)
after=coords();shift=max(float(np.linalg.norm(after[k]-v)) for k,v in before.items());assert shift<.002
with (out/'receptor-h.pdb').open('w') as w:PDBFile.writeFile(f.topology,f.positions,w,keepIds=True)
meta={'seconds':time.perf_counter()-start,'residues':sum(1 for _ in f.topology.residues()),'observed_heavy_atoms':len(before),'prepared_heavy_atoms':len(after),'max_observed_shift_A':shift,'added_missing_heavy_atoms':missing,'missing_residues_not_built':skipped,'pH':7.4,'platform':'CPU, 2 threads','note':'No full protein minimization. Experimental heavy coordinates fixed. Long missing loops omitted and segments independently terminalized.'}
(out/'preparation.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta),flush=True)
