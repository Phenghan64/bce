import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parents[1]/'results/large-receptors';s=json.loads((out/'summary.json').read_text())
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');colors=['#126b83','#bd6735'];names=['albumin','mtor']
for name,color in zip(names,colors):
 t=s['targets'][name];label=f"{name.title()} ({t['actual_residues']} residues)";samples=[json.loads(x) for x in (out/name/'resources.jsonl').read_text().splitlines()];mins=np.array([x['elapsed_seconds']/60 for x in samples])
 axes[0,1].plot(mins,[x['rss_sum_bytes']/2**30 for x in samples],label=label,color=color)
 axes[1,0].plot(mins,[x['system_cpu_percent'] for x in samples],color=color,alpha=.8,label=label)
 lib={r['id']:r for r in json.loads((out/name/'library.json').read_text())['records']}
 rows=sorted(t['per_ligand_median_seconds'].items(),key=lambda kv:lib[kv[0]]['rdkit_rotatable_bonds'])
 x=np.array([lib[k]['rdkit_rotatable_bonds'] for k,v in rows]);y=np.array([v for k,v in rows]);axes[1,1].plot(x,y,'o-',label=label,color=color)
 labels=[f"{s['targets'][n]['actual_residues']} residues" for n in names]
vals=[s['targets'][n]['panel_seconds']/60 for n in names];bars=axes[0,0].bar(labels,vals,color=colors,width=.5)
for b,v in zip(bars,vals):axes[0,0].text(b.get_x()+b.get_width()/2,v+.03,f'{v:.2f} min',ha='center')
axes[0,0].set(title='Panel elapsed time: 10 ligands x 3 seeds',ylabel='Minutes (30 jobs; controls excluded)',ylim=(0,max(vals)*1.18))
axes[0,1].set(title='Task process-tree memory',xlabel='Minutes from each batch start',ylabel='RSS sum (GiB)');axes[0,1].legend(frameon=False,fontsize=9)
axes[1,0].set(title='Whole-system CPU (including other apps)',xlabel='Minutes from each batch start',ylabel='CPU utilization (%)',ylim=(0,100))
axes[1,1].set(title='Paired per-ligand job times',xlabel='RDKit rotatable bonds (one ligand per count)',ylabel='Median of 3 seeds (seconds)',xticks=range(10));axes[1,1].legend(frameon=False,fontsize=9)
for ax in axes.flat:
 ax.grid(alpha=.15);ax.set_axisbelow(True)
 for spine in ['top','right']:ax.spines[spine].set_visible(False)
fig.suptitle('Local Vina benchmark: two large experimental receptors\nSame 10 ligands | E32 | 3 concurrent jobs x 2 threads | 30 A pocket box',fontsize=15,y=.99)
fig.get_layout_engine().set(rect=(0,.07,1,.83))
fig.text(.5,.02,'All observed residues of each selected chain retained; mTOR is a truncated construct.\nRSS may double-count shared memory. Run order: albumin, then mTOR. This is not a biological hit-rate benchmark.',ha='center',fontsize=9,color='#5a6570')
fig.savefig(out/'comparison.png',dpi=170);fig.savefig(out/'comparison.svg')
