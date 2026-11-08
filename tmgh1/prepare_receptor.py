"""Prepare one receptor geometry/state with original heavy atoms frozen."""
import sys,json,time,random,subprocess,hashlib
from pathlib import Path
import numpy as np
from pdbfixer import PDBFixer
from openmm import Platform,Context,VerletIntegrator,LocalEnergyMinimizer,unit
from openmm.app import PDBFile,ForceField,NoCutoff,Modeller
root=Path(__file__).resolve().parents[1];out=root/'results/tmgh1-round1';geometry,state=sys.argv[1:];source=out/'geometries'/geometry/'observed.pdb';d=out/'receptors'/f'{geometry}_{state}';d.mkdir(parents=True,exist_ok=True);assert not (d/'preparation.json').exists();start=time.perf_counter();random.seed(42);np.random.seed(42)
f=PDBFixer(filename=str(source),platform=Platform.getPlatformByName('CPU'))
def get_coords(top,positions):
 xyz=np.asarray(positions.value_in_unit(unit.angstrom));return {(a.residue.chain.id,a.residue.id,a.name):xyz[a.index].copy() for a in top.atoms() if a.element.symbol!='H'}
before=get_coords(f.topology,f.positions);f.findMissingResidues();skipped={str(k):v for k,v in f.missingResidues.items()};f.missingResidues={};f.findMissingAtoms();missing={f'{r.chain.id}:{r.id}':[a.name for a in aa] for r,aa in f.missingAtoms.items()};f.addMissingAtoms(seed=42)
pos=np.asarray(f.positions.value_in_unit(unit.angstrom))
for r in f.topology.residues():
 if r.name not in ['ASP','GLU']:continue
 cn,attached,o1,o2=('CG','CB','OD1','OD2') if r.name=='ASP' else ('CD','CG','OE1','OE2');atoms={a.name:a.index for a in r.atoms()}
 if any((r.chain.id,r.id,n) in before for n in [o1,o2]):continue
 origin=pos[atoms[cn]];u=pos[atoms[attached]]-origin;u/=np.linalg.norm(u);v=pos[atoms[o1]]-origin;v-=u*np.dot(u,v);v/=np.linalg.norm(v)
 for name,sign in [(o1,1),(o2,-1)]:pos[atoms[name]]=origin+1.25*(np.cos(np.deg2rad(117))*u+sign*np.sin(np.deg2rad(117))*v)
ff=ForceField('amber14-all.xml');model=Modeller(f.topology,pos*unit.angstrom);variants=[]
for r in model.topology.residues():
 variants.append(('GLH' if state=='E166H' else 'GLU') if r.id=='166' else ('GLU' if r.id=='351' else None))
assigned=model.addHydrogens(ff,pH=7.0,variants=variants,platform=Platform.getPlatformByName('CPU'))
system=ff.createSystem(model.topology,nonbondedMethod=NoCutoff,constraints=None);frozen=[]
for a in model.topology.atoms():
 if a.element.symbol!='H' and (a.residue.chain.id,a.residue.id,a.name) in before:system.setParticleMass(a.index,0);frozen.append(a.index)
assert len(frozen)==len(before);integ=VerletIntegrator(.001);ctx=Context(system,integ,Platform.getPlatformByName('CPU'),{'Threads':'2'});ctx.setPositions(model.positions);e0=float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole));LocalEnergyMinimizer.minimize(ctx,tolerance=5,maxIterations=10000);final=ctx.getState(getPositions=True,getEnergy=True);after=get_coords(model.topology,final.getPositions());shift=max(float(np.linalg.norm(after[k]-v)) for k,v in before.items());assert shift<.002
with (d/'receptor-h.pdb').open('w') as w:PDBFile.writeFile(model.topology,final.getPositions(),w,keepIds=True)
protein_time=time.perf_counter()-start
# Pin catalytic templates to the explicitly prepared H states.
chain166=next(r.chain.id for r in model.topology.residues() if r.id=='166');chain351=next(r.chain.id for r in model.topology.residues() if r.id=='351')
cmd=[str(root/'.env/bin/mk_prepare_receptor.py'),'--read_pdb',str(d/'receptor-h.pdb'),'-o',str(d/'receptor'),'-p','-j','--set_template',f'{chain166}:166='+('GLH' if state=='E166H' else 'GLU')+f',{chain351}:351=GLU']
t=time.perf_counter();result=subprocess.run(cmd,capture_output=True,text=True);(d/'meeko.log').write_text(result.stdout+result.stderr)
report={'geometry':geometry,'state':state,'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'protein_preparation_seconds':protein_time,'meeko_seconds':time.perf_counter()-t,'total_seconds':time.perf_counter()-start,'meeko_returncode':result.returncode,'observed_heavy_atoms':len(before),'prepared_heavy_atoms':len(after),'residues':sum(1 for _ in model.topology.residues()),'max_original_heavy_shift_A':shift,'added_missing_heavy_atoms':missing,'missing_residues_not_built':skipped,'energy_before_kJmol':e0,'energy_after_kJmol':float(final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)),'assigned_residue_variants':{f'{r.chain.id}:{r.id}':v for r,v in zip(model.topology.residues(),assigned) if v},'command':cmd,'pH':7.0};(d/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');result.check_returncode()
qt=(d/'receptor.pdbqt').read_text();lines=[x for x in qt.splitlines() if x.startswith(('ATOM  ','HETATM'))];counts=len(set((x[21],x[22:27]) for x in lines));assert counts==report['residues'];print(json.dumps({k:report[k] for k in ['geometry','state','residues','total_seconds','meeko_returncode']}),flush=True)
