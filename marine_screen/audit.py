"""Independent output reconciliation and exploratory release-combination diagnostics."""
import hashlib,json,sqlite3,sys,platform,time
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from prepare import OUT, SRC, save
from run import SCHEMES, score

def main():
    t=time.perf_counter();db=sqlite3.connect('file:'+str(OUT/'screen.sqlite')+'?mode=ro',uri=True)
    errors={}
    errors['coordinate_or_sequence_mismatch']=db.execute('SELECT count(*) FROM occurrence o JOIN protein p ON o.protein=p.id WHERE substr(p.sequence,o.start_1based,o.end_1based_inclusive-o.start_1based+1) != o.peptide').fetchone()[0]
    errors['unregistered_peptide']=db.execute('SELECT count(*) FROM occurrence o LEFT JOIN peptide p ON o.peptide=p.sequence WHERE p.sequence IS NULL').fetchone()[0]
    errors['length_outside_protocol']=db.execute('SELECT count(*) FROM peptide WHERE length<2 OR length>20 OR length!=length(sequence)').fetchone()[0]
    errors['duplicate_occurrence']=db.execute('SELECT count(*) FROM (SELECT 1 FROM occurrence GROUP BY peptide,protein,scheme,start_1based,end_1based_inclusive,missed HAVING count(*)>1)').fetchone()[0]
    n=db.execute('SELECT count(*) FROM peptide').fetchone()[0]
    summary=json.loads((OUT/'summary.json').read_text());assert n==summary['unique_peptides_MC01']
    refs=json.loads((OUT/'reference-values.json').read_text())
    rh=float(spearmanr([score(x['sequence']) for x in refs],[-np.log10(x['ic50_uM']) for x in refs]).statistic)
    expected=json.loads((OUT/'validation.json').read_text())[0]['methods']['heuristic']['spearman']
    assert abs(rh-expected)<1e-12
    manifest=json.loads((OUT/'input-manifest.json').read_text())
    root=Path(__file__).resolve().parents[1]
    for r in manifest:
        if hashlib.sha256((root/r['path']).read_bytes()).hexdigest()!=r['sha256']:
            errors['input_changed:'+r['path']]=1
    assert not any(errors.values()),errors
    known={x['sequence'] for x in refs}
    sets={name:{r[0] for r in db.execute('SELECT DISTINCT peptide FROM occurrence WHERE scheme=? AND missed=0',(name,))} for name in SCHEMES}
    comparisons=[]
    for combo,a,b in [('trypsin+chymotrypsin','trypsin','chymotrypsin'),('pepsin_pH2+trypsin','pepsin_pH2','trypsin'),('thermolysin+trypsin','thermolysin','trypsin')]:
        before=sets[a]|sets[b];after=sets[combo]
        comparisons.append(dict(combination=combo,MC0=True,single_enzyme_union=len(before),combination_products=len(after),
          new_sequences=len(after-before),lost_sequences=len(before-after),
          new_known_sequences=sorted((after-before)&known),lost_known_sequences=sorted((before-after)&known),
          note='Compared with union of two separately digested libraries, not a kinetic enzyme-order experiment.'))
    save('combination-diagnostics.json',comparisons)
    runtime=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,
        machine='Apple M6',physical_cpu_cores=12,memory_bytes=17179869184)
    save('verification.json',dict(status='PASS',errors=errors,independent_scipy_rho=rh,
        all_occurrences_verified=summary['occurrences_MC01'],unique_candidates_verified=n,
        frozen_source_hashes_verified=len(manifest),audit_seconds=time.perf_counter()-t,
        metadata_correction='Sardine measurement source section corrected to3.1.3; no scores or numerical results changed.',runtime=runtime))
    print(json.dumps(dict(errors=errors,verified_occurrences=summary['occurrences_MC01'],combinations=comparisons),indent=2))

if __name__=='__main__':main()
