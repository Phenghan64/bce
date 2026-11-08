from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/marine-screen';SRC=OUT/'source'
def fetch(name,url,params=None):
 p=SRC/name
 if p.exists():return p
 r=requests.get(url,params=params,timeout=60);r.raise_for_status();p.write_bytes(r.content)
 (SRC/(name+'.meta.json')).write_text(json.dumps({'url':r.url,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest(),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
 return p
if __name__=='__main__':
 tasks=[('mackerel.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9222390/fullTextXML'),('sardine.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10252845/fullTextXML'),('aceip-repo.json','https://api.github.com/repos/Qian0324/ACEip_DB/contents'),('zenodo.json','https://zenodo.org/api/records/21708604'),('expasy.html','https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html')]
 with concurrent.futures.ThreadPoolExecutor(5) as ex:
  fs={ex.submit(fetch,*t):t for t in tasks}
  for f in concurrent.futures.as_completed(fs):
   try:p=f.result();print(p.name,p.stat().st_size,flush=True)
   except Exception as e:print(fs[f],repr(e),flush=True)
