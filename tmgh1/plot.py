"""Export scientific figures for structure fidelity, docking outcomes and local resource use."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from rdkit import Chem
out=Path(__file__).resolve().parents[1]/'results/tmgh1-round1';s=json.loads((out/'summary.json').read_text());st=json.loads((out/'structure-analysis.json').read_text());d=json.loads((out/'docking.json').read_text());colors=['#126b83','#b36c38','#7b6ba3']
fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
for m,color in zip(st['models'],colors):
 errors=m['comparisons']['2J77']['per_residue_CA_errors_A'];ids=sorted(map(int,errors));axes[0,0].plot(ids,[errors[str(i)] for i in ids],color=color,lw=1,label=f"OF3 seed {m['seed']}")
for i in [166,351]:axes[0,0].axvline(i,color='#777',ls=':',lw=1)
axes[0,0].set(xlabel='Native residue number',ylabel='CA deviation (A)',title='Sequence-only predictions vs 2J77 B');axes[0,0].legend(fontsize=8,frameon=False)
geoms=['experimental_2J77','experimental_2J7E','predicted_42','predicted_2026','predicted_20261021'];labels=['2J77','2J7E','OF3 42','OF3 2026','OF3 20261021']
for lid,color,offset in [('NOJ-cation','#126b83',-.13),('GI2-cation','#b36c38',.13)]:
 for i,g in enumerate(geoms):
  group=next(r for r in s['groups'] if r['geometry']==g and r['receptor_state']=='default' and r['ligand_id']==lid);ys=group['top1_RMSDs_A'];axes[0,1].scatter(np.array([i+offset]*3)+np.array([-.025,0,.025]),ys,color=color,s=28,alpha=.8,label=lid if i==0 else None)
axes[0,1].axhline(2,color='#a34444',ls='--',lw=1,label='2 A pose threshold');axes[0,1].set(xticks=range(5),xticklabels=labels,ylabel='Top-ranked ligand RMSD (A)',title='Primary states: all three docking seeds');axes[0,1].tick_params(axis='x',rotation=15);axes[0,1].legend(fontsize=8,frameon=False)
b=json.loads((out/'inference/model-benchmark.json').read_text());samples=b['samples'];minutes=np.array([r['seconds']/60 for r in samples]);axes[1,0].plot(minutes,[r['process_rss_bytes']/2**30 for r in samples],label='Process RSS',color='#126b83');axes[1,0].plot(minutes,[r['mps_driver_bytes']/2**30 for r in samples],label='MPS driver allocation',color='#b36c38');axes[1,0].set(xlabel='Minutes from model benchmark start',ylabel='GiB (overlapping measures)',title='Prediction memory: do not add the curves');axes[1,0].legend(fontsize=8,frameon=False)
r=[json.loads(x) for x in (out/'inference/resources.jsonl').read_text().splitlines()];t=[x['elapsed_seconds']/60 for x in r];base=r[0]['swap_used_bytes'];axes[1,1].plot(t,[x['system_available_bytes']/2**30 for x in r],label='System available memory',color='#126b83');axes[1,1].plot(t,[(x['swap_used_bytes']-base)/2**30 for x in r],label='System swap growth',color='#b36c38');axes[1,1].set(xlabel='Minutes from outer monitor start',ylabel='GiB',title='Whole-system memory during 3 predictions');axes[1,1].legend(fontsize=8,frameon=False)
for ax in axes.flat:
 ax.grid(alpha=.16);ax.set_axisbelow(True)
 for spine in ['top','right']:ax.spines[spine].set_visible(False)
fig.suptitle('TmGH1 round 1: can a predicted enzyme support known-site docking?\n446 native residues | no templates or ligand input | 3 OpenFold3 seeds',fontsize=15,y=.995)
fig.get_layout_engine().set(rect=(0,.06,1,.85));fig.text(.5,.018,'Primary comparison uses fixed seed 42; other seeds assess repeatability. Historical structures may have appeared in training.\nNo affinity-ranking validation. Whole-system resource metrics include other applications.',ha='center',fontsize=9,color='#5a6570');fig.savefig(out/'round1-comparison.png',dpi=180);fig.savefig(out/'round1-comparison.svg')
# Inspectable 3D pose figure in the shared protein-aligned frame.
def atom_coords(path):
 return {(int(l[22:26]),l[12:16].strip()):np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in path.read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip() not in ['H','D']}
def sticks(ax,m,color,label=None):
 xyz=m.GetConformer().GetPositions()
 for bond in m.GetBonds():
  a,b=xyz[bond.GetBeginAtomIdx()],xyz[bond.GetEndAtomIdx()];ax.plot(*np.array([a,b]).T,color=color,lw=2.5)
 ax.scatter(*xyz.T,c=color,s=17,label=label,depthshade=False)
fig=plt.figure(figsize=(12,7),layout='constrained')
for j,(lig,code) in enumerate([('NOJ','2J77'),('GI2','2J7E')],1):
 ax=fig.add_subplot(1,2,j,projection='3d');r=next(r for r in d['results'] if r['spec']['geometry']=='predicted_42' and r['spec']['receptor_state']=='default' and r['spec']['ligand_id']==lig+'-cation' and r['spec']['seed']==42);ref=Chem.SDMolSupplier(r['spec']['reference'],removeHs=True)[0];mol=Chem.SDMolSupplier(r['sdf'],removeHs=True)[0];center=ref.GetConformer().GetPositions().mean(0);pc=atom_coords(out/'geometries/predicted_42/observed.pdb');rc=atom_coords(out/'geometries'/f'experimental_{code}'/'observed.pdb');pocket=next(x for x in st['models'] if x['seed']==42)['comparisons'][code]['pocket_residues_6A']
 for coord,color in [(rc,'#9aabad'),(pc,'#75afbd')]:
  for residue in pocket:
   keys=[(residue,n) for n in ['N','CA','C'] if (residue,n) in coord];xyz=np.array([coord[k] for k in keys]);ax.plot(*xyz.T,color=color,lw=1.1,alpha=.55)
  for residue in [166,351]:
   xyz=np.array([coord[residue,n] for n in ['CA','CB','CG','CD','OE1']]);ax.plot(*xyz.T,color=color,lw=1.8)
 sticks(ax,ref,'#263a42','Crystal ligand');sticks(ax,mol,'#d2813b','Docked on predicted receptor');xyz=np.concatenate([ref.GetConformer().GetPositions(),mol.GetConformer().GetPositions()]);radius=max(7,float(np.abs(xyz-center).max()+2));ax.set(xlim=(center[0]-radius,center[0]+radius),ylim=(center[1]-radius,center[1]+radius),zlim=(center[2]-radius,center[2]+radius),xlabel='X (A)',ylabel='Y (A)',zlabel='Z (A)',title=f'{lig}: fixed OF3 seed 42 / Vina seed 42\nTop-pose RMSD = {r["top_RMSD_A"]:.2f} A');ax.view_init(elev=22,azim=-55);ax.set_box_aspect((1,1,1));ax.tick_params(labelsize=8);ax.set_title(ax.get_title(),fontsize=11)
fig.suptitle('Protein-aligned binding-pose comparison (no independent ligand fitting)',fontsize=14,y=.995);fig.legend(handles=[Line2D([0],[0],color='#263a42',lw=2,label='Crystal ligand'),Line2D([0],[0],color='#d2813b',lw=2,label='Docked ligand'),Line2D([0],[0],color='#9aabad',lw=2,label='Experimental pocket backbone'),Line2D([0],[0],color='#75afbd',lw=2,label='Predicted pocket backbone')],loc='upper center',bbox_to_anchor=(.5,.95),ncol=4,frameon=False,fontsize=9);fig.get_layout_engine().set(rect=(0,.02,1,.84));fig.savefig(out/'binding-poses.png',dpi=180);fig.savefig(out/'binding-poses.svg')
