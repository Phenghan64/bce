"""Fetch the exact user-requested PubChem panel, retaining raw provenance."""
from pathlib import Path
import hashlib,json,time,urllib.request
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/ligand-library20-preflight-20261030/source'
CIDS=[2519,2244,3672,1983,156391,338,5280343,5280863,5280443,5280445,5280961,5281708,969516,445154,65064,9064,370,445858,689043,1794427]
NAMES=['咖啡因','阿司匹林','布洛芬','对乙酰氨基酚','萘普生','水杨酸','槲皮素','山奈酚','芹菜素','木犀草素','染料木素','大豆苷元','姜黄素','反式白藜芦醇','表没食子儿茶素没食子酸酯（EGCG）','（+）-儿茶素','没食子酸','阿魏酸','咖啡酸','绿原酸']
BASE='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+','.join(map(str,CIDS))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=[];opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
    for name,suffix in [('properties.json','/property/Title,MolecularFormula,MolecularWeight,IUPACName,IsomericSMILES,InChIKey,Charge/JSON'),('structures-2d.sdf','/SDF?record_type=2d')]:
        dest=OUT/name;url=BASE+suffix
        if dest.exists():
            print('Retained source:',name,flush=True);continue
        t=time.perf_counter()
        for attempt in range(3):
            try:
                with opener.open(urllib.request.Request(url,headers={'User-Agent':'EYOSON-Ligand-Preflight/1.0'}),timeout=60) as r:raw=r.read()
                if not raw:raise ValueError('Empty PubChem response')
                dest.write_bytes(raw);break
            except Exception:
                if attempt==2:raise
                time.sleep(2*(attempt+1))
        manifest.append({'file':name,'url':url,'retrieved_at':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'elapsed_seconds':round(time.perf_counter()-t,3)})
        print('Fetched:',name,len(raw),'bytes',flush=True)
        time.sleep(1)
    records=json.loads((OUT/'properties.json').read_text())['PropertyTable']['Properties']
    assert {x['CID'] for x in records}==set(CIDS)
    (OUT.parent/'requested-panel.json').write_text(json.dumps([{'id':f'L{i:03}','cid':cid,'name_zh':NAMES[i-1]} for i,cid in enumerate(CIDS,1)],ensure_ascii=False,indent=2)+'\n')
    if manifest:(OUT/'retrieval.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
if __name__=='__main__':main()
