"""Static scientific figure: experimental labels versus frozen docking ranking."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

out=Path(__file__).resolve().parents[1]/'results/abl-real-case';r=json.loads((out/'analysis.json').read_text())
assert r['metrics'] is not None
rows=r['rows'];fig,ax=plt.subplots(figsize=(11.5,8.7),layout='constrained')
for i,row in enumerate(rows):
    color='#167a9a' if row['label']=='strong' else '#bc6241'
    x=row['score'];scores=row['scores']
    ax.plot([min(scores),max(scores)],[i,i],color=color,lw=2,alpha=.6)
    ax.scatter(scores,[i]*3,color=color,s=18,alpha=.6)
    ax.scatter([x],[i],color=color,s=65,marker='D',edgecolor='white',zorder=3)
names=[row['name'] if len(row['name'])<30 else 'CID '+row['cid'] for row in rows]
ax.set_yticks(range(16),[f'{i+1:2d}  {name}' for i,name in enumerate(names)])
ax.set_ylim(15.8,-1.1);ax.axhline(3.5,color='#9aa6ad',ls='--',lw=1)
ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
ax.set_xlabel('Vina score (kcal/mol); lower is ranked better\nDiamond: median of 3 seeds. Dots/range: individual runs.')
for spine in ['top','right','left']:ax.spines[spine].set_visible(False)
ax.tick_params(axis='y',length=0)
right=ax.twinx();right.set_ylim(ax.get_ylim());right.set_yticks(range(16),[f'{x["kd_nm"]:g} nM' if x['label']=='strong' else 'capped / not detected' for x in rows]);right.tick_params(length=0,labelsize=9)
for spine in right.spines.values():spine.set_visible(False)
right.set_ylabel('Davis 2011 experimental binding label',labelpad=12)
m=r['metrics'];lo,hi=m['roc_auc_stratified_bootstrap95']
fig.suptitle(f'Human ABL1: experimental labels vs. Vina ranking\nAUC {m["roc_auc"]:.3f} (bootstrap 95% {lo:.2f}–{hi:.2f}); top 4: {m["top4_strong_count"]}/4 strong binders',fontsize=14)
ax.legend(handles=[Line2D([0],[0],marker='D',color='none',markerfacecolor='#167a9a',label='Strong: measured Kd ≤100 nM'),Line2D([0],[0],marker='D',color='none',markerfacecolor='#bc6241',label='Weak / not detected: dataset cap 10000 nM')],loc='lower center',bbox_to_anchor=(.5,1.01),frameon=False,fontsize=9)
fig.get_layout_engine().set(rect=(0,.065,1,.94))
fig.text(.5,.025,'16 preselected compounds; 2HYY B rigid receptor; one heuristic protonation state; 3 seeds at E32.\nRetrospective extreme-group pilot. Scores are not measured affinities; no OpenFold3 prediction in this run.',ha='center',fontsize=9,color='#505b64')
fig.savefig(out/'ranking.png',dpi=180)
fig.savefig(out/'ranking.svg')
