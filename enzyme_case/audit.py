"""Independent result checks plus explicitly post-hoc stronger length baseline.

No changes to fitted models, selected candidates or primary frozen endpoints.
"""
import gzip,json,math,time,hashlib,os
from collections import Counter
from pathlib import Path
import numpy as np
from Bio import SeqIO
from run import CONDITIONS, AA, high_conf, exact_p

ROOT=Path(__file__).resolve().parents[1];OUT=Path(os.environ.get('ENZYME_CASE_OUTPUT',ROOT/'results/enzyme-case')).resolve()
t=time.perf_counter()
rows=json.loads((OUT/'source/cod-peptides.json').read_text())
split=json.loads((OUT/'parent-split.json').read_text())
train=set(split['train_parents']);test=set(split['test_parents'])
full={r.id:str(r.seq).rstrip('*') for r in SeqIO.parse(gzip.open(OUT/'source/cod-ensembl90.fa.gz','rt'),'fasta') if r.id in train|test}
test_text='|'.join(full[p] for p in sorted(test))
train_lengths={c:Counter() for c in CONDITIONS};truth={c:set() for c in CONDITIONS}
for r in rows:
    s=r['Sequence']
    if not 7<=len(s)<=30 or not set(s)<=set(AA):continue
    parents=[p for p in (r['Proteins'] or '').split(';') if p in full and s in full[p]]
    for c in CONDITIONS:
        if high_conf(r,c):
            if any(p in test for p in parents):truth[c].add(s)
            elif s not in test_text and any(p in train for p in parents):train_lengths[c][len(s)]+=1
N=json.loads((OUT/'release-validation.json').read_text())['unique_virtual_peptides']
lengths=np.empty(N,dtype=np.uint8);labels=np.zeros((N,6),dtype=bool)
models=json.loads((OUT/'release-model.json').read_text())['fits']
idx={a:i for i,a in enumerate(AA)}
checks=0;total=0
with gzip.open(OUT/'virtual-candidates.jsonl.gz','rt') as f:
    for i,line in enumerate(f):
        r=json.loads(line);s=r['sequence'];lengths[i]=len(s)
        labels[i]=[s in truth[c] for c in CONDITIONS]
        # Check deterministic spread of candidates, plus the first entries.
        if i<10 or i%45000==0:
            for p,start,end in r['sources']:assert p in test and full[p][start-1:end]==s
            for c,m in enumerate(models):
                q=np.array(m['pair_cut_probability_surrogate']);g=m['global_cut_frequency'];values=[]
                for p,start,end in r['sources']:
                    seq=full[p];start-=1;v=0
                    for b in (start,end):
                        if 0<b<len(seq):
                            prob=q[idx[seq[b-1]],idx[seq[b]]] if seq[b-1] in idx and seq[b] in idx else g
                            v+=math.log(prob)
                    for b in range(start+1,end):v+=math.log1p(-q[idx[seq[b-1]],idx[seq[b]]])
                    values.append(v)
                assert abs(max(values)-r['scores'][c])<1e-8
            checks+=1
        total+=1
assert total==N
counts=Counter(map(int,lengths));results=[]
for j,c in enumerate(CONDITIONS):
    # Expected training peptide mass by length, divided over available candidates of each length.
    weight=np.array([(train_lengths[c][int(l)]+1)/counts[int(l)] for l in lengths])
    result={'condition':c,'training_length_counts':dict(train_lengths[c]),'heldout_positives':int(labels[:,j].sum())}
    for fraction in (.01,.1):
        k=math.ceil(N*fraction);threshold=np.partition(weight,N-k)[N-k];above=weight>threshold;ties=weight==threshold
        result[f'recall_top{int(100*fraction)}pct']=float((labels[above,j].sum()+labels[ties,j].sum()*(k-above.sum())/ties.sum())/labels[:,j].sum())
    results.append(result)
observed=json.loads((OUT/'observed-peptide-scores.json').read_text())
assays=json.loads((OUT/'experimental-assays.json').read_text())
nooverlap=[]
for score,endpoint,direction in [('FRS','DPPH_IC50',-1),('Chelator','Chelation_IC50',-1),('Emulsi_max','EAI',1)]:
    pred=[];cover=[]
    for c in CONDITIONS:
        rr=[r for r in observed if (r['conditions'][c]['intensity'] or 0)>0 and r['predictions'].get(score) is not None and not r['predictions']['training_exact_overlap']]
        weights=[r['conditions'][c]['intensity'] for r in rr]
        pred.append(float(np.average([r['predictions'][score] for r in rr],weights=weights)))
    nooverlap.append({'score':score,'endpoint':endpoint,**exact_p(pred,[a[endpoint]*direction for a in assays])})
result={'status':'Post-hoc secondary audit after primary results; no refitting or changing selection.',
        'candidate_rows_checked':total,'independently_recomputed_candidates_all_six_recipes':checks,
        'length_baseline_method':'Laplace-smoothed training observed length frequency divided by heldout candidate count of that length; tie-aware expected recall. No heldout detection labels in the baseline.',
        'length_baseline':results,'bulk_excluding_exact_training_overlap':nooverlap,'seconds':time.perf_counter()-t}
(OUT/'secondary-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
