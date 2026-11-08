"""Fallback streaming sequence workload; explicitly not the full retained-candidate pipeline."""
import gzip,json,sys,time
from pathlib import Path
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1];out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);start=time.monotonic();counts=0;checksum=0
for rec in SeqIO.parse(gzip.open(ROOT/'results/enzyme-case/source/cod-ensembl90.fa.gz','rt'),'fasta'):
    s=str(rec.seq)
    for length in range(7,31):
        for pos in range(len(s)-length+1):
            peptide=s[pos:pos+length];counts+=1;checksum+=peptide.count('P')+peptide.count('L')
(out/'summary.json').write_text(json.dumps({'enumerated_occurrences':counts,'checksum':checksum,'seconds':time.monotonic()-start,'mode':'streaming fallback, no retained full model results'})+'\n')
