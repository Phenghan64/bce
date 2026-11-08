"""Reproducible cod-frame screening benchmark. Run from repository root."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '2')
import gzip, hashlib, itertools, json, math, sqlite3, time
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import spearmanr, rankdata
from models import PublishedModels, ANOX, ROOT

OUT = Path(os.environ.get('ENZYME_CASE_OUTPUT', ROOT/'results/enzyme-case')).resolve()
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT/'results/enzyme-case/source'
AA = 'ACDEFGHIKLMNPQRSTVWY'
A = {a:i for i,a in enumerate(AA)}
CONDITIONS = ['Alc-MCF','Neut-MCF','Neut&Alc-MCF','Alc-HCM','Neut-HCM','Neut&Alc-HCM']
PHASES = {}
START = time.perf_counter()
def write(name, data):
    (OUT/name).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
def stage(name, started):
    PHASES[name] = time.perf_counter()-started
    print(name, round(PHASES[name],3), 's', flush=True)
def high_conf(r,c):
    return (r['PEP'] is not None and r['PEP']<=.01 and
            r['conditions'][c]['identification']=='By MS/MS' and (r['conditions'][c]['intensity'] or 0)>0)
def exact_p(x,y):
    x=rankdata(x);y=rankdata(y)
    obs=float(np.corrcoef(x,y)[0,1])
    if not np.isfinite(obs):return {'rho':None,'p_exact':None}
    extreme=sum(abs(np.corrcoef(x,p)[0,1])>=abs(obs)-1e-10 for p in itertools.permutations(y))
    return {'rho':obs,'p_exact':extreme/720}

def prepare():
    t=time.perf_counter()
    rows=json.loads((SOURCE/'cod-peptides.json').read_text())
    assert len(rows)==len({r['Sequence'] for r in rows})==11968
    full={r.id:str(r.seq).rstrip('*') for r in SeqIO.parse(gzip.open(SOURCE/'cod-ensembl90.fa.gz','rt'),'fasta')}
    mapped={}
    for r in rows:
        mappings=[]
        for p in (r['Proteins'] or '').split(';'):
            if p not in full:continue
            s=full[p];start=s.find(r['Sequence'])
            while start>=0:
                mappings.append((p,start,start+len(r['Sequence'])))
                start=s.find(r['Sequence'],start+1)
        mapped[r['Sequence']]=mappings
    active=sorted({p for maps in mapped.values() for p,_,_ in maps})
    full={p:full[p] for p in active}
    expected=[4414,6268,3702,4020,5871,4361]
    expected_lengths=[13.0,17.4,9.36,13.1,16.7,10.2]
    reproduction=[]
    for c,n,l in zip(CONDITIONS,expected,expected_lengths):
        rr=[r for r in rows if (r['conditions'][c]['intensity'] or 0)>0]
        w=np.array([r['conditions'][c]['intensity'] for r in rr]);length=float(np.average([r['Length'] for r in rr],weights=w))
        assert len(rr)==n and abs(length-l)<.051
        reproduction.append({'condition':c,'count':len(rr),'paper_count':n,'weighted_length':length,'paper_weighted_length':l,'weighted_mass':float(np.average([r['Mass'] for r in rr],weights=w))})
    audit={'peptides':len(rows),'mapped_peptides':sum(bool(x) for x in mapped.values()),'mapped_parents':len(full),'unmapped_peptides':sum(not x for x in mapped.values()),'reproduction':reproduction}
    write('data-audit.json',audit)
    stage('load_and_verify',t)
    return rows,full,mapped,audit

def split_parents(full,mapped):
    parent={p:p for p in full}
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    seen={}
    for p,s in full.items():
        for i in range(len(s)-9):
            k=s[i:i+10]
            if k in seen:union(p,seen[k])
            else:seen[k]=p
    for s,locations in mapped.items():
        if len(s)>=7 and locations:
            for p,_,_ in locations[1:]:union(p,locations[0][0])
    groups=defaultdict(list)
    for p in full:groups[find(p)].append(p)
    test=set()
    for group,members in groups.items():
        if int(hashlib.sha256(group.encode()).hexdigest(),16)%5==0:test.update(members)
    train=set(full)-test
    assert train and test
    detail={'train_parents':sorted(train),'test_parents':sorted(test),'components':dict(groups),'method':'shared 10-mer / shared annotated observed peptide union components; SHA256 modulo5 holdout'}
    write('parent-split.json',detail)
    return train,test

def release_screen(rows,full,mapped):
    t=time.perf_counter()
    train,test=split_parents(full,mapped)
    # Exact overlap exclusion includes peptides occurring in test proteins without an original annotation.
    test_text='|'.join(full[p] for p in sorted(test))
    train_rows=[r for r in rows if 7<=len(r['Sequence'])<=30 and set(r['Sequence'])<=set(AA) and r['Sequence'] not in test_text]
    probs=[];globals_=[];fit=[]
    for c in CONDITIONS:
        cuts=np.zeros((20,20));intact=np.zeros((20,20));used=0;peptides=set()
        for r in train_rows:
            if not high_conf(r,c):continue
            for p,s,e in mapped[r['Sequence']]:
                if p not in train:continue
                used+=1;peptides.add(r['Sequence']);seq=full[p]
                for b in (s,e):
                    if 0<b<len(seq) and seq[b-1] in A and seq[b] in A:cuts[A[seq[b-1]],A[seq[b]]]+=1
                for b in range(s+1,e):
                    if seq[b-1] in A and seq[b] in A:intact[A[seq[b-1]],A[seq[b]]]+=1
        g=float(cuts.sum()/(cuts.sum()+intact.sum()))
        assert 0<g<1 and len(peptides)>20
        q=(cuts+20*g)/(cuts+intact+20)
        probs.append(q);globals_.append(g)
        fit.append({'condition':c,'unique_training_peptides':len(peptides),'training_occurrences':used,'global_cut_frequency':g,'pair_cut_probability_surrogate':q.tolist()})
    write('release-model.json',{'alphabet':AA,'fits':fit,'not_a_yield_model':True})
    # Full enumeration, no enzyme union rule and no label-based candidate pruning.
    candidates={};occurrences=0
    for p in sorted(test):
        seq=full[p];n=len(seq)
        qs=np.ones((6,n+1))
        for b in range(1,n):
            qs[:,b]=[q[A[seq[b-1]],A[seq[b]]] if seq[b-1] in A and seq[b] in A else g for q,g in zip(probs,globals_)]
        cut=np.log(qs)
        intact=np.zeros_like(qs);intact[:,1:n]=np.log1p(-qs[:,1:n])
        prefix=np.cumsum(intact,axis=1)
        for length in range(7,31):
            for s in range(n-length+1):
                e=s+length;pep=seq[s:e]
                if not set(pep)<=set(AA):continue
                scores=cut[:,s]+cut[:,e]+prefix[:,e-1]-prefix[:,s]
                bases=np.array([(int(s>0)+int(e<n))*math.log(g)+(length-1)*math.log1p(-g) for g in globals_])
                occurrence=[p,s+1,e];occurrences+=1
                if pep not in candidates:
                    candidates[pep]={'sequence':pep,'sources':[occurrence],'scores':scores,'baseline':bases}
                else:
                    r=candidates[pep];r['scores']=np.maximum(r['scores'],scores);r['baseline']=np.maximum(r['baseline'],bases);r['sources'].append(occurrence)
    assert len(candidates)>=10000
    seqs=list(candidates);scores=np.array([candidates[s]['scores'] for s in seqs]);bases=np.array([candidates[s]['baseline'] for s in seqs])
    observed={c:{r['Sequence'] for r in rows if high_conf(r,c) and 7<=len(r['Sequence'])<=30 and any(p in test for p,_,_ in mapped[r['Sequence']])} for c in CONDITIONS}
    metrics=[];selected=set();ranks={};top_members={}
    for j,c in enumerate(CONDITIONS):
        truth=observed[c]&candidates.keys()
        assert truth
        labels=np.array([s in truth for s in seqs]);denom=int(labels.sum())
        order=np.argsort(-scores[:,j],kind='stable');baseorder=np.argsort(-bases[:,j],kind='stable')
        selected.update(seqs[i] for i in order[:10000]);top_members[c]={seqs[i]:int(k+1) for k,i in enumerate(order[:10000])}
        row={'condition':c,'heldout_observed_peptides':denom,'unique_candidates':len(seqs),'top10000_observed':int(labels[order[:10000]].sum())}
        for fraction in (.01,.1):
            k=math.ceil(len(seqs)*fraction)
            # Fractional allocation at a tied cutoff avoids arbitrary sequence-order gains.
            def recall(values):
                threshold=np.partition(values,len(values)-k)[len(values)-k]
                above=values>threshold;ties=values==threshold
                return float((labels[above].sum()+labels[ties].sum()*(k-above.sum())/ties.sum())/denom)
            actual=recall(scores[:,j]);baseline=recall(bases[:,j])
            row[f'recall_top{int(fraction*100)}pct']=actual
            row[f'baseline_recall_top{int(fraction*100)}pct']=baseline
            row[f'enrichment_top{int(fraction*100)}pct']=actual/(k/len(seqs))
        metrics.append(row)
    with gzip.open(OUT/'virtual-candidates.jsonl.gz','wt') as f:
        for s,r in candidates.items():
            f.write(json.dumps({'sequence':s,'sources':r['sources'],'scores':r['scores'].tolist(),'baseline':r['baseline'].tolist()},separators=(',',':'))+'\n')
    virtual=[{'Sequence':s,'sources':candidates[s]['sources'],'release_scores':candidates[s]['scores'].tolist(),'release_top10000_rank':{c:top_members[c][s] for c in CONDITIONS if s in top_members[c]},'observed_conditions':[c for c in CONDITIONS if s in observed[c]]} for s in sorted(selected)]
    result={'train_parents':len(train),'test_parents':len(test),'enumerated_occurrences':occurrences,'unique_virtual_peptides':len(seqs),'recipe_scores':len(seqs)*6,'virtual_peptides_functionally_selected':len(virtual),'metrics':metrics,'conditions_order':CONDITIONS}
    write('release-validation.json',result)
    stage('virtual_digestion_and_holdout_validation',t)
    return virtual,result

def functional_screen(rows,virtual):
    t=time.perf_counter();models=PublishedModels();stage('model_load',t)
    t=time.perf_counter()
    seqs=sorted({r['Sequence'] for r in rows+virtual})
    valid=[s for s in seqs if 3<=len(s)<=30 and set(s)<=set(AA)]
    prediction=models.anox(valid)
    training=pd.read_csv(ANOX/'Data/03_p70_AO_db.csv')
    training_sequences=set(training['Sequence'])
    training_detail={'file':'03_p70_AO_db.csv','rows':len(training),'source_counts':training['Source'].fillna('missing').value_counts().to_dict(),'observed_exact_overlap':sorted(set(r['Sequence'] for r in rows)&training_sequences),'virtual_exact_overlap':sorted(set(r['Sequence'] for r in virtual)&training_sequences)}
    write('training-overlap.json',training_detail)
    scores={}
    for s,pred in zip(valid,prediction):
        em=models.emulsify(s) if len(s)>=7 else [None,None,None]
        scores[s]={'FRS':float(pred[0]),'Chelator':float(pred[1]),'Emulsi_alpha':em[0],'Emulsi_beta':em[1],'Emulsi_gamma':em[2],'Emulsi_max':max(em) if em[0] is not None else None,'length':len(s),'GRAVY':sum(models.hydro[c] for c in s)/len(s),'training_exact_overlap':s in training_sequences}
    for r in rows+virtual:r['predictions']=scores.get(r['Sequence'],{})
    write('observed-peptide-scores.json',rows)
    with gzip.open(OUT/'virtual-shortlist-scores.jsonl.gz','wt') as f:
        for r in virtual:f.write(json.dumps(r,separators=(',',':'))+'\n')
    db=sqlite3.connect(OUT/'screen.sqlite')
    db.executescript('DROP TABLE IF EXISTS peptide; DROP TABLE IF EXISTS observed; DROP TABLE IF EXISTS virtual; CREATE TABLE peptide(sequence TEXT PRIMARY KEY,length INTEGER,FRS REAL,Chelator REAL,Emulsi_alpha REAL,Emulsi_beta REAL,Emulsi_gamma REAL,Emulsi_max REAL,GRAVY REAL,training_exact_overlap INTEGER); CREATE TABLE observed(sequence TEXT,condition TEXT,intensity REAL,relative_intensity REAL,identification TEXT,PEP REAL,proteins TEXT,source_excel_row INTEGER); CREATE TABLE virtual(sequence TEXT PRIMARY KEY,source_positions_json TEXT,release_scores_json TEXT,top10000_ranks_json TEXT);')
    for s,r in scores.items():db.execute('INSERT INTO peptide VALUES(?,?,?,?,?,?,?,?,?,?)',(s,r['length'],r['FRS'],r['Chelator'],r['Emulsi_alpha'],r['Emulsi_beta'],r['Emulsi_gamma'],r['Emulsi_max'],r['GRAVY'],int(r['training_exact_overlap'])))
    for r in rows:
        for c in CONDITIONS:
            v=r['conditions'][c]
            if (v['intensity'] or 0)>0:db.execute('INSERT INTO observed VALUES(?,?,?,?,?,?,?,?)',(r['Sequence'],c,v['intensity'],v['relative'],v['identification'],r['PEP'],r['Proteins'],r['source_row']))
    for r in virtual:db.execute('INSERT INTO virtual VALUES(?,?,?,?)',(r['Sequence'],json.dumps(r['sources']),json.dumps(r['release_scores']),json.dumps(r['release_top10000_rank'])))
    db.commit();assert db.execute('PRAGMA integrity_check').fetchone()[0]=='ok';db.close()
    result={'observed_total':len(rows),'observed_AnOx':sum(bool(r['predictions']) for r in rows),'observed_Emulsi':sum(r['predictions'].get('Emulsi_max') is not None for r in rows),'virtual_AnOx':len(virtual),'distinct_AnOx_evaluated':len(valid),'distinct_Emulsi_evaluated':sum(s['Emulsi_max'] is not None for s in scores.values()),'training_overlap_observed':len(training_detail['observed_exact_overlap']),'training_overlap_virtual':len(training_detail['virtual_exact_overlap'])}
    write('functional-scale.json',result);stage('functional_inference_and_storage',t)
    return rows,virtual,result

def bulk_validation(rows):
    t=time.perf_counter()
    # Published Table 1, in CONDITIONS order. IC50 units mg/mL, EAI m2/g, ESI min.
    values={'DPPH_IC50':[4.97,4.36,4.31,3.58,3.36,2.93],'Chelation_IC50':[.77,.53,.77,1.01,.69,1.65],'EAI':[88.18,80.75,87.54,68.35,66.41,76.48],'ESI':[29.38,33.39,36.10,65.15,46.82,45.70],'DH':[36.1,27.0,39.2,30.2,18.7,30.2]}
    sd={'DPPH_IC50':[.18,.08,.06,.21,.18,.05],'Chelation_IC50':[.03,.013,.03,.31,.02,.07],'EAI':[7.94,4.32,4.46,.72,.97,3.99],'ESI':[3.52,.56,5.12,2.03,4.75,1.94],'DH':[.03,.12,.04,.05,.09,.04]}
    groups={'DPPH_IC50':['a','b','b','c','c','d'],'Chelation_IC50':['bc','c','bc','b','c','a'],'EAI':['a','ab','a','c','c','b']}
    assays=[{'condition':c,**{k:v[i] for k,v in values.items()},'SD':{k:v[i] for k,v in sd.items()},'significance_groups':{k:v[i] for k,v in groups.items()}} for i,c in enumerate(CONDITIONS)]
    write('experimental-assays.json',assays)
    aggregate=[]
    for c in CONDITIONS:
        detected=[r for r in rows if (r['conditions'][c]['intensity'] or 0)>0]
        total=sum(r['conditions'][c]['intensity'] for r in detected)
        for mode in ('weighted','unweighted','high_confidence','length_ge10'):
            for score in ('FRS','Chelator','Emulsi_max'):
                rr=[r for r in detected if r['predictions'].get(score) is not None and (mode!='high_confidence' or high_conf(r,c)) and (mode!='length_ge10' or len(r['Sequence'])>=10)]
                w=np.array([1 if mode=='unweighted' else r['conditions'][c]['intensity'] for r in rr])
                aggregate.append({'condition':c,'mode':mode,'score':score,'n':len(rr),'intensity_coverage':sum(r['conditions'][c]['intensity'] for r in rr)/total,'value':float(np.average([r['predictions'][score] for r in rr],weights=w)),'length_baseline':float(np.average([len(r['Sequence']) for r in rr],weights=w)),'GRAVY_baseline':float(np.average([r['predictions']['GRAVY'] for r in rr],weights=w))})
    comparisons=[]
    for score,endpoint,direction in [('FRS','DPPH_IC50',-1),('Chelator','Chelation_IC50',-1),('Emulsi_max','EAI',1),('Emulsi_max','ESI',1)]:
        for mode in ('weighted','unweighted','high_confidence','length_ge10'):
            data=[next(r for r in aggregate if r['score']==score and r['mode']==mode and r['condition']==c) for c in CONDITIONS]
            target=np.array(values[endpoint])*direction
            pred=np.array([r['value'] for r in data])
            within=[]
            for material,indices in [('MCF',[0,1,2]),('HCM',[3,4,5])]:
                local=pred[indices];truth=target[indices];best=indices[int(np.argmax(local))];actual=indices[int(np.argmax(truth))]
                group=groups.get(endpoint)
                overlaps=bool(set(group[best])&set(group[actual])) if group else None
                within.append({'material':material,'predicted_best':CONDITIONS[best],'experimental_numeric_best':CONDITIONS[actual],'same_significance_group':overlaps,'within_rho':float(spearmanr(local,truth).statistic)})
            comparisons.append({'score':score,'endpoint':endpoint,'mode':mode,**exact_p(pred,target),'length_baseline':exact_p([r['length_baseline'] for r in data],target),'GRAVY_baseline':exact_p([r['GRAVY_baseline'] for r in data],target),'within_material':within})
    write('bulk-aggregates.json',aggregate);write('bulk-validation.json',comparisons)
    # A transparent shortlist: credible observed peptides ranked by the published score, without assay fitting.
    shortlist={}
    for endpoint in ('FRS','Chelator','Emulsi_max'):
        eligible=[r for r in rows if len(r['Sequence'])>=10 and r['predictions'].get(endpoint) is not None and any(high_conf(r,c) for c in CONDITIONS)]
        eligible.sort(key=lambda r:(-r['predictions'][endpoint],r['Sequence']))
        shortlist[endpoint]=[{'sequence':r['Sequence'],'score':r['predictions'][endpoint],'PEP':r['PEP'],'proteins':r['Proteins'],'source_row':r['source_row'],'training_exact_overlap':r['predictions']['training_exact_overlap'],'conditions':{c:r['conditions'][c] for c in CONDITIONS if (r['conditions'][c]['intensity'] or 0)>0},'individual_activity_experimentally_verified':False} for r in eligible[:20]]
    write('experiment-candidates.json',shortlist)
    stage('bulk_experimental_validation',t)
    return comparisons

def main():
    protocol=ROOT/'enzyme_case/protocol.json'
    write('run-start.json',{'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'protocol_sha256':hashlib.sha256(protocol.read_bytes()).hexdigest()})
    rows,full,mapped,audit=prepare()
    virtual,release=release_screen(rows,full,mapped)
    rows,virtual,scale=functional_screen(rows,virtual)
    validation=bulk_validation(rows)
    PHASES['total_compute_seconds']=time.perf_counter()-START
    write('timing.json',PHASES)
    print(json.dumps({'scale':scale,'release':release,'primary_bulk':[x for x in validation if x['mode']=='weighted'],'timing':PHASES},ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
