"""Traceable, bounded virtual digestion; no network or model downloads during compute."""
import collections, gzip, hashlib, json, math, re, sqlite3, time
from pathlib import Path
import numpy as np
from scipy.stats import rankdata, spearmanr
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from prepare import OUT, SRC, AA, save

PROTOCOL = json.loads((OUT/'protocol.json').read_text())
RULES = {k: re.compile(v) for k,v in PROTOCOL['rules'].items()}
SCHEMES = PROTOCOL['schemes']
SEED = PROTOCOL['seed']

def cuts(seq, enzymes):
    sites = {0, len(seq)}
    for enzyme in enzymes:
        s = {m.end() for m in RULES[enzyme].finditer(seq)}
        if enzyme == 'trypsin':
            s -= {m.end() for m in RULES['trypsin_exception'].finditer(seq)}
        sites.update(s)
    return sorted(sites)

def fragments(seq, enzymes, max_mc=1, bounded=True):
    sites = cuts(seq, enzymes)
    for mc in range(max_mc+1):
        for i in range(len(sites)-mc-1):
            a,b=sites[i],sites[i+mc+1]
            if not bounded or 2 <= b-a <= 20:
                yield seq[a:b],a,b,mc

def score(seq):
    return 2*int(seq[-1] in 'FWYP') + int(seq[-2] in 'AVILMFYW') + int(len(seq)<=5)

def properties(seq):
    p=ProteinAnalysis(seq)
    return dict(length=len(seq),mw_Da=p.molecular_weight(),gravy=p.gravy(),charge_pH7=p.charge_at_pH(7),
                heuristic_score=score(seq))

def references():
    rows=json.loads((SRC/'aceip-rows.json').read_text()); accepted={};rejected=[]
    for rownum,x in enumerate(rows,2):
        seq=str(x.get('Seq') or '');v=str(x.get('ic50') or '')
        m=re.fullmatch(r'\s*(\d+(?:\.\d+)?)\s*(μM|µM|uM|mM)\s*',v)
        reason=('noncanonical_or_length' if set(seq)-AA or not 2<=len(seq)<=20 else
                'no_DOI' if not x.get('doi') else 'nonexact_or_unrecognized_unit' if not m else
                'nonpositive' if float(m[1])<=0 else None)
        if reason: rejected.append(dict(row=rownum,sequence=seq,raw_ic50=v,reason=reason));continue
        if seq in accepted: raise ValueError('duplicate experimental sequence: '+seq)
        accepted[seq]=dict(sequence=seq,ic50_uM=float(m[1])*(1000 if m[2]=='mM' else 1),
                           raw_ic50=v,doi=x['doi'],reference=x['reference'],source_sheet='ACEiP_Properties',source_row=rownum)
    save('reference-exclusions.json',rejected);save('reference-values.json',list(accepted.values()))
    return accepted,collections.Counter(x['reason'] for x in rejected)

def top_weights(scores,k):
    s=np.asarray(scores,float); threshold=np.sort(s)[-k]
    w=(s>threshold).astype(float);boundary=s==threshold
    w[boundary]=(k-w.sum())/boundary.sum()
    return w

def fast_rho(a,b):
    a=rankdata(a);b=rankdata(b);a-=a.mean();b-=b.mean()
    den=np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/den) if den else None

def metrics(rows, name):
    y=np.array([-math.log10(x['ic50_uM']) for x in rows]);n=len(rows)
    k=max(1,math.ceil(.1*n));strong=np.array([x['ic50_uM']<=10 for x in rows])
    result=dict(name=name,n=n,strong_threshold_uM=10,strong_count=int(strong.sum()),top_k=k,
                caveat='Measured, already ACE-active literature set; not precision or specificity for untested candidates.',methods={})
    for method in ['heuristic','shorter','gravy']:
        s=np.array([score(x['sequence']) if method=='heuristic' else -len(x['sequence']) if method=='shorter'
                    else ProteinAnalysis(x['sequence']).gravy() for x in rows],float)
        rng=np.random.default_rng(SEED);rho=fast_rho(s,y);boot=[]
        for _ in range(1000):
            ix=rng.integers(0,n,n);r=fast_rho(s[ix],y[ix])
            if r is not None:boot.append(r)
        sr=rankdata(s)-rankdata(s).mean();yr=rankdata(y)-rankdata(y).mean()
        den=np.linalg.norm(sr)*np.linalg.norm(yr)
        null=[float(sr@rng.permutation(yr)/den) for _ in range(2000)] if den else []
        w=top_weights(s,k);hits=float(w@strong);base=float(strong.mean())
        result['methods'][method]=dict(spearman=rho,bootstrap95=np.quantile(boot,[.025,.975]).tolist() if boot else None,
            permutation_p_two_sided=(1+sum(abs(r)>=abs(rho)-1e-12 for r in null))/2001 if null else None,
            expected_topk_strong=hits,random_expected_topk_strong=k*base,
            topk_precision_within_measured=hits/k,topk_recall=hits/strong.sum() if strong.sum() else None,
            enrichment=hits/k/base if base else None,boundary_ties_fractional=True)
    return result

def process_metrics(seq, enzymes, maxmc):
    f=list(fragments(seq,enzymes,maxmc));high=sum(score(x[0])>=3 for x in f)
    return dict(products=len(f),high_score=high,density_per1000=1000*high/len(seq))

def benchmark(refs):
    proteins=json.loads((OUT/'benchmark-proteins.json').read_text());out=[]
    for alt in [False,True]:
        panel=[p for p in proteins if p['accession'] != ('P13104' if alt else 'P84335')]
        for mc in [0,1]:
            for name,enzymes in SCHEMES.items():
                parts=[process_metrics(p['sequence'],enzymes,mc) for p in panel]
                out.append(dict(scheme=name,max_missed=mc,tropomyosin='P84335' if alt else 'P13104',
                    protein_count=len(panel),products=sum(x['products'] for x in parts),
                    high_score=sum(x['high_score'] for x in parts),
                    primary_density=1000*sum(x['high_score'] for x in parts)/sum(len(p['sequence']) for p in panel),
                    equal_protein_density=float(np.mean([x['density_per1000'] for x in parts]))))
    measured=['pepsin_pH2','thermolysin','proteinase_K']
    primary=sorted([r for r in out if r['max_missed']==0 and r['tropomyosin']=='P13104' and r['scheme'] in measured],key=lambda r:-r['primary_density'])
    return dict(rows=out,experimental_winner='pepsin_pH2',compared_measured_enzymes=measured,
                predicted_primary_order=[r['scheme'] for r in primary],winner_recovered=primary[0]['scheme']=='pepsin_pH2',
                random_winner_probability=1/3,
                limits='Only3of6 actual enzymes modeled; eight cross-species homologs, unknown abundances; heuristic is theoretical release, not measured inhibition.')

def main():
    started=time.perf_counter();timings={};proteins=json.loads((OUT/'proteins.json').read_text())
    refs,rejections=references();save('reference-audit.json',dict(accepted=len(refs),rejected_by_reason=dict(rejections)))
    dbpath=OUT/'screen.sqlite'
    if dbpath.exists():raise RuntimeError('Refuse overwrite of completed/partial database; use versioned run for rerun.')
    db=sqlite3.connect(dbpath)
    db.executescript('''
    CREATE TABLE protein(id TEXT PRIMARY KEY, sequence TEXT NOT NULL, entries_json TEXT NOT NULL);
    CREATE TABLE occurrence(peptide TEXT NOT NULL, protein TEXT NOT NULL, scheme TEXT NOT NULL,
      start_1based INTEGER NOT NULL, end_1based_inclusive INTEGER NOT NULL, missed INTEGER NOT NULL);
    CREATE TABLE peptide(sequence TEXT PRIMARY KEY, length INTEGER, mw_Da REAL, gravy REAL,
      charge_pH7 REAL, heuristic_score INTEGER, reference_ic50_uM REAL, reference_doi TEXT);
    ''')
    db.executemany('INSERT INTO protein VALUES(?,?,?)',[(p['id'],p['sequence'],json.dumps(p['entries'],ensure_ascii=False)) for p in proteins]);db.commit()
    t=time.perf_counter();unique=set();zero=set();buffer=[];per_protein=[]
    for i,p in enumerate(proteins):
        for name,enzymes in SCHEMES.items():
            for pep,a,b,mc in fragments(p['sequence'],enzymes):
                unique.add(pep)
                if mc==0:zero.add(pep)
                buffer.append((pep,p['id'],name,a+1,b,mc))
                if len(buffer)>=30000:
                    db.executemany('INSERT INTO occurrence VALUES(?,?,?,?,?,?)',buffer);buffer.clear()
            for mc in [0,1]:
                per_protein.append(dict(protein=p['id'],scheme=name,max_missed=mc,**process_metrics(p['sequence'],enzymes,mc)))
        if (i+1)%100==0:print('digested',i+1,'proteins; unique peptides',len(unique),flush=True)
    if buffer:db.executemany('INSERT INTO occurrence VALUES(?,?,?,?,?,?)',buffer)
    db.commit();timings['digestion_and_occurrence_write_s']=time.perf_counter()-t
    t=time.perf_counter();buffer=[]
    for pep in sorted(unique):
        x=properties(pep);r=refs.get(pep,{})
        buffer.append((pep,x['length'],x['mw_Da'],x['gravy'],x['charge_pH7'],x['heuristic_score'],r.get('ic50_uM'),r.get('doi')))
        if len(buffer)>=10000:
            db.executemany('INSERT INTO peptide VALUES(?,?,?,?,?,?,?,?)',buffer);buffer.clear()
    if buffer:db.executemany('INSERT INTO peptide VALUES(?,?,?,?,?,?,?,?)',buffer)
    db.executescript('CREATE INDEX occurrence_peptide ON occurrence(peptide); CREATE INDEX occurrence_scheme ON occurrence(scheme,missed);')
    db.commit();timings['properties_and_index_s']=time.perf_counter()-t
    schemes=[]
    for name in SCHEMES:
        for mc in [0,1]:
            n,nunique=db.execute('SELECT count(*),count(DISTINCT peptide) FROM occurrence WHERE scheme=? AND missed<=?',(name,mc)).fetchone()
            known=db.execute('SELECT count(DISTINCT o.peptide) FROM occurrence o JOIN peptide p ON o.peptide=p.sequence WHERE o.scheme=? AND o.missed<=? AND p.reference_ic50_uM IS NOT NULL',(name,mc)).fetchone()[0]
            vals=[v for v in per_protein if v['scheme']==name and v['max_missed']==mc]
            schemes.append(dict(scheme=name,max_missed=mc,occurrences=n,unique_peptides=nunique,reference_matches=known,
                                 equal_protein_mean_density=float(np.mean([v['density_per1000'] for v in vals]))))
    # Species+name longest-entry sensitivity. Preserve species strata rather than invent mixture abundances.
    representative={}
    for p in proteins:
        for e in p['entries']:
            k=(e['species'],e['name'])
            if k not in representative or len(p['sequence'])>len(representative[k]['sequence']):representative[k]=p
    lookup={(r['protein'],r['scheme'],r['max_missed']):r for r in per_protein};species_rows=[]
    for sp in sorted({k[0] for k in representative}):
        selected={p['id']:p for k,p in representative.items() if k[0]==sp}
        for name in SCHEMES:
            for mc in [0,1]:
                vals=[lookup[(pid,name,mc)]['density_per1000'] for pid in selected]
                species_rows.append(dict(species=sp,scheme=name,max_missed=mc,representative_proteins=len(selected),
                                         mean_density=float(np.mean(vals))))
    save('schemes.json',schemes);save('species-sensitivity.json',species_rows)
    with gzip.open(OUT/'protein-scheme-results.jsonl.gz','wt') as f:
        for r in per_protein:f.write(json.dumps(r)+'\n')
    t=time.perf_counter();validation=[metrics(list(refs.values()),'curated_literature_all')]
    overlap=[x for seq,x in refs.items() if seq in unique]
    validation.append(metrics(overlap,'marine_candidate_reference_overlap'))
    sardine_values={'FIGR':(.60,.01),'FILR':(1.01,.02),'FQRL':(1.54,.14),'FRAL':(.66,.09),'KFL':(.61,.02),'KLF':(.89,.01)}
    sardine=[]
    for seq,(value,sd) in sardine_values.items():
        mw=properties(seq)['mw_Da']
        sardine.append(dict(sequence=seq,ic50_mg_mL=value,sd_mg_mL=sd,mw_Da=mw,ic50_uM=value*1e6/mw,
            heuristic_score=score(seq),marine_MC0_present=seq in zero,marine_MC01_present=seq in unique,
            doi='10.3390/foods12112216',source_section='3.1.3',
            note='Observed in Sardina hydrolysate; library coverage in other fish is not reproduction of the sardine cleavage process.'))
    validation.append(metrics(sardine,'sardine_six_primary_measurements'));save('sardine-validation.json',sardine)
    save('validation.json',validation);save('process-benchmark.json',benchmark(refs))
    timings['experimental_backtest_s']=time.perf_counter()-t
    # Two explicitly separate queues: literature-supported exact matches and untested heuristics.
    db.row_factory=sqlite3.Row
    for name,where,order in [('literature-shortlist','reference_ic50_uM IS NOT NULL','reference_ic50_uM,sequence'),
                              ('exploratory-shortlist','reference_ic50_uM IS NULL','heuristic_score DESC,length,sequence')]:
        data=[]
        for row in db.execute(f'SELECT * FROM peptide WHERE {where} ORDER BY {order} LIMIT 30'):
            rec=dict(row);q=rec['sequence']
            rec['provenance']=[dict(r) for r in db.execute('SELECT * FROM occurrence WHERE peptide=? ORDER BY missed,scheme,protein,start_1based LIMIT 12',(q,))]
            rec['total_occurrences']=db.execute('SELECT count(*) FROM occurrence WHERE peptide=?',(q,)).fetchone()[0]
            rec['selection_note']='Tie order is length/lexical bookkeeping, not validated activity order.'
            data.append(rec)
        save(name+'.json',data)
    occurrences=db.execute('SELECT count(*) FROM occurrence').fetchone()[0]
    integrity=db.execute('PRAGMA integrity_check').fetchone()[0]
    summary=dict(protein_sequences=len(proteins),schemes=len(SCHEMES),protein_scheme_jobs=len(proteins)*len(SCHEMES),
        unique_peptides_MC0=len(zero),unique_peptides_MC01=len(unique),MC1_only=len(unique-zero),
        occurrences_MC01=occurrences,known_exact_matches=len(overlap),
        known_le10uM_matches=sum(x['ic50_uM']<=10 for x in overlap),
        heuristic_score_counts=dict(collections.Counter(score(s) for s in unique)),
        stage_seconds=timings,compute_seconds=time.perf_counter()-started,sqlite_integrity=integrity,
        protocol_sha256=hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest(),
        sequence_database_bytes=dbpath.stat().st_size)
    save('summary.json',summary);db.close();print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
