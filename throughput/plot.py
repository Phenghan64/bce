"""Create an exportable scientific plot of measured throughput and resource use."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parents[1]/'results/vina-throughput-1000'
s=json.loads((out/'summary.json').read_text());screen=json.loads((out/'screen.json').read_text())
samples=[json.loads(line) for line in (out/'resources.jsonl').read_text().splitlines()]
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
times=sorted(r['stage_elapsed_at_completion']/60 for r in screen['results'])
axes[0,0].plot([0]+times,range(len(times)+1),color='#126b83',lw=2)
axes[0,0].set(xlabel='Minutes since screening started',ylabel='Completed screening jobs',title='Actual cumulative screen completions')
t=np.array([r['elapsed_seconds']/60 for r in samples])
axes[0,1].plot(t,[r['rss_sum_bytes']/2**30 for r in samples],label='Process-tree RSS sum',color='#126b83')
axes[0,1].plot(t,[r['system_available_bytes']/2**30 for r in samples],label='System available RAM',color='#849875')
axes[0,1].set(xlabel='Minutes since monitored run started',ylabel='GiB',title='Sampled memory use');axes[0,1].legend(frameon=False,fontsize=9)
b=s['blocks_of100_completions']
axes[1,0].bar(range(1,len(b)+1),[x['submitted_per_hour'] for x in b],color='#126b83',alpha=.85)
axes[1,0].axhline(s['screen_submitted_per_hour'],ls='--',color='#bf6c35',label='Whole-screen average')
axes[1,0].set(xlabel='Successive groups of 100 completions',ylabel='Jobs/hour',title='Batch rates (compound complexity varies)');axes[1,0].set_xticks(range(1,len(b)+1));axes[1,0].legend(frameon=False,fontsize=9)
lib={r['id']:r for r in json.loads((out/'library.json').read_text())['records']}
x=[lib[r['spec']['ligand_id']]['rdkit_rotatable_bonds'] for r in screen['results']]
y=[r['job_wall_seconds'] for r in screen['results']]
axes[1,1].scatter(x,y,alpha=.25,s=14,color='#126b83');axes[1,1].set(xlabel='RDKit rotatable bonds',ylabel='Individual job seconds',title='Per-job time at 3 concurrent jobs × 2 threads')
for ax in axes.flat:
    ax.grid(alpha=.15);ax.set_axisbelow(True)
    for spine in ['top','right']:ax.spines[spine].set_visible(False)
fig.suptitle(f"Local Vina throughput: {len(times):,} screening jobs\n{s['screen_submitted_per_hour']:,.0f} jobs/hour at E8; {s['vina_jobs']:,} total docking runs including calibration/refinement",fontsize=15,y=.99)
fig.get_layout_engine().set(rect=(0,.065,1,.84))
fig.text(.5,.018,'Human ABL1 2HYY B; 999 prepared NCI library molecules + imatinib. Source ligand states retained.\nUnlabelled throughput benchmark: no experimental hit-rate claim. RSS sum may double count shared pages.',ha='center',fontsize=9,color='#5a6570')
fig.savefig(out/'throughput.png',dpi=180);fig.savefig(out/'throughput.svg')
