"""Read-only Table S3 extraction, with optional verification against saved JSON.

Usage: python extract_source.py --check (requires openpyxl).
Without --check, write the extracted JSON to stdout; the workbook is never edited.
"""
import json,sys
from pathlib import Path
from openpyxl import load_workbook
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/enzyme-case/source'
workbook=load_workbook(SOURCE/'cod-supp/marinedrugs-1006687-SI-proofed/Table S3.xlsx',read_only=True,data_only=True)
sheet=workbook['All data'];iterator=sheet.iter_rows(values_only=True);header=next(iterator)
columns={name:i for i,name in enumerate(header)}
fields=['Sequence','N-term cleavage window','C-term cleavage window','Length','Mass','Proteins','Leading razor protein','Start position','End position','Unique (Proteins)','PEP','Score','Reverse','Potential contaminant','id']
conditions=['Alc-MCF','Neut-MCF','Neut&Alc-MCF','Alc-HCM','Neut-HCM','Neut&Alc-HCM']
records=[]
for excel_row,row in enumerate(iterator,2):
    if not row[columns['Sequence']]:continue
    r={field:row[columns[field]] for field in fields};r['source_row']=excel_row
    r['conditions']={c:{'intensity':row[columns[f'Intensity {c}_Crude_TopDown_3-65']], 'relative':row[columns[f'{c}_Rel_Int']], 'identification':row[columns[f'Identification type {c}_Crude_TopDown_3-65']]} for c in conditions}
    records.append(r)
workbook.close()
if '--check' in sys.argv:
    original=json.loads((SOURCE/'cod-peptides.json').read_text())
    assert records==original,'Workbook extraction differs from the saved pipeline input'
    print(f'PASS: {len(records)} records exactly match the original workbook values and row positions.')
else:json.dump(records,sys.stdout,ensure_ascii=False)
