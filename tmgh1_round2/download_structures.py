"""Download the compound-to-CCD mapping from SI Table B and official RCSB records."""
import requests,json,concurrent.futures,hashlib,time
from pathlib import Path
out=Path(__file__).resolve().parents[1]/'results/tmgh1-round2';source=out/'source';source.mkdir(exist_ok=True)
pdbs=['2J77','1OIF','2J75','1W3J','2J7H','2CBU','2CBV','1UZ1','2J78','2J79','2J7B','2CES','2CET','2J7C','2J7D','2J7E','2J7F','2J7G']
def get(url,path):
 if path.exists():return json.loads(path.read_text()) if path.suffix=='.json' else path.read_text()
 for i in range(3):
  try:
   r=requests.get(url,timeout=35);r.raise_for_status();path.write_bytes(r.content);return r.json() if path.suffix=='.json' else r.text
  except Exception:
   if i==2:raise
   time.sleep(1)
def job(it):
 i,pdb=it;entry=get(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}',source/f'{pdb}-entry.json');entities=[get(f'https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb}/{eid}',source/f'{pdb}-entity-{eid}.json') for eid in entry['rcsb_entry_container_identifiers']['non_polymer_entity_ids']]; candidates=[e for e in entities if e['pdbx_entity_nonpoly']['comp_id'] not in ['ACT','CA','HOH','CL','SO4','GOL','EDO']];assert len(candidates)==1,(pdb,candidates);j=candidates[0];code=j['pdbx_entity_nonpoly']['comp_id'];ccd=get(f'https://data.rcsb.org/rest/v1/core/chemcomp/{code}',source/f'{code}.json');get(f'https://files.rcsb.org/ligands/download/{code}.cif',source/f'{code}.cif');return {'compound':i,'pdb':pdb,'ccd':code,'name':j['pdbx_entity_nonpoly']['name'],'descriptors':ccd.get('pdbx_chem_comp_descriptor'),'chem_comp':ccd['chem_comp'],'mapping_source':'ACS SI Table B, pages S2-S5'}
with concurrent.futures.ThreadPoolExecutor(4) as ex:records=list(ex.map(job,enumerate(pdbs,1)))
(out/'compound-map.json').write_text(json.dumps(records,indent=2)+'\n')
for r in records:print(r['compound'],r['pdb'],r['ccd'],r['name'],[(d['type'],d['descriptor']) for d in r['descriptors'] if d['type']=='SMILES_CANONICAL'])
