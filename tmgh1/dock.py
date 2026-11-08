"""Run the frozen 120-job pose-transfer experiment with resumable input/output validation."""
import sys,os,json,time,hashlib,subprocess,signal,re,math,concurrent.futures
from pathlib import Path
from datetime import datetime,timezone
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule,RDKitMolCreate
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root));from realcase.audit_outputs import compare_identity
out=root/'results/tmgh1-round1';records=json.loads((out/'library.json').read_text())['records'];protocol=out/'protocol.json';vina=root/'.tools/bin/vina';geometries=['experimental_2J77','experimental_2J7E','predicted_42','predicted_2026','predicted_20261021'];states=['default','E166H'];seeds=[42,2026,20261021]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,x):
 tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
def protein(path):
 lines=[l for l in path.read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip() not in ['H','D']];return np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lines]),np.array([int(l[22:26]) for l in lines])
def contacts(xyz,res,ligand):return set(map(int,res[np.linalg.norm(xyz[:,None]-ligand[None,:],axis=2).min(1)<4]))
reference_contacts={}
for ligand,code in [('NOJ','2J77'),('GI2','2J7E')]:
 xyz,res=protein(out/'geometries'/f'experimental_{code}'/'observed.pdb');mol=Chem.SDMolSupplier(str(out/'reference'/f'{ligand}.sdf'))[0];reference_contacts[ligand]=contacts(xyz,res,mol.GetConformer().GetPositions())
version=subprocess.check_output([str(vina),'--version'],text=True).strip();base={'protocol_sha256':sha(protocol),'library_sha256':sha(out/'library.json'),'vina':version,'exhaustiveness':32,'cpu':2,'num_modes':9,'energy_range':4,'timeout_seconds':900}
boxes=out/'boxes';boxes.mkdir(exist_ok=True)
for ligand,code in [('NOJ','2J77'),('GI2','2J7E')]:
 meta=json.loads((out/'geometries'/f'experimental_{code}'/'metadata.json').read_text());(boxes/f'{ligand}.txt').write_text('\n'.join([f'center_{axis} = {value}' for axis,value in zip('xyz',meta['box_center'])]+[f'size_{axis} = 24' for axis in 'xyz'])+'\n')
lock=out/'docking.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
def job(geom,state,entry,seed):
 rid=f'{geom}_{state}';receptor=out/'receptors'/rid/'receptor.pdbqt';box=boxes/f'{entry["ligand"]}.txt';folder=out/'jobs'/rid/f'{entry["id"]}-s{seed}';folder.mkdir(parents=True,exist_ok=True)
 spec={**base,'geometry':geom,'receptor_state':state,'ligand_id':entry['id'],'seed':seed,'receptor':str(receptor),'receptor_sha256':sha(receptor),'ligand':entry['pdbqt'],'ligand_sha256':sha(entry['pdbqt']),'box':str(box),'box_sha256':sha(box),'reference':entry['reference'],'reference_sha256':sha(entry['reference'])};signature=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest();rp=folder/'result.json'
 if rp.exists():
  old=json.loads(rp.read_text());assert old['signature']==signature,'Inputs changed; use new output version'
  if old['status']=='completed' and all(Path(p).exists() and sha(p)==h for p,h in old['outputs_sha256'].items()):return {**old,'cache_hit':True}
  if old['status']=='failed':return {**old,'cache_hit':True}
 started=time.perf_counter();dest=folder/'docked.pdbqt';cmd=['/usr/bin/time','-l',str(vina),'--receptor',str(receptor),'--ligand',entry['pdbqt'],'--config',str(box),'--cpu','2','--exhaustiveness','32','--seed',str(seed),'--num_modes','9','--energy_range','4','--out',str(dest)];r={'status':'running','signature':signature,'spec':spec,'command':cmd,'started_utc':datetime.now(timezone.utc).isoformat(),'cache_hit':False};atomic(rp,r)
 try:
  with (folder/'vina.log').open('w') as log:
   p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   try:code=p.wait(timeout=900)
   except subprocess.TimeoutExpired:
    os.killpg(p.pid,signal.SIGTERM)
    try:p.wait(timeout=5)
    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
    raise TimeoutError('Vina job timeout')
  r.update(returncode=code,vina_wall_seconds=time.perf_counter()-started);assert code==0
  text=dest.read_text();scores=[float(v) for v in re.findall(r'^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)',text,re.M)];assert scores and all(math.isfinite(v) for v in scores) and scores==sorted(scores)
  sdf,fail=RDKitMolCreate.write_sd_string(PDBQTMolecule(text,skip_typing=True));assert not fail;sp=folder/'docked.sdf';sp.write_text(sdf);mols=list(Chem.SDMolSupplier(str(sp),removeHs=True));assert len(mols)==len(scores)==len(re.findall(r'^MODEL\s+',text,re.M));expected=Chem.MolFromSmiles(entry['smiles']);reference=Chem.SDMolSupplier(entry['reference'],removeHs=True)[0];rms=[]
  for m in mols:
   assert m is not None and compare_identity(m,expected)[0];assert np.isfinite(m.GetConformer().GetPositions()).all();rms.append(float(rdMolAlign.CalcRMS(m,reference)))
  xyz,res=protein(out/'receptors'/rid/'receptor-h.pdb');observed=contacts(xyz,res,mols[0].GetConformer().GetPositions());expected_c=reference_contacts[entry['ligand']];intersection=len(observed&expected_c);contact={'native_residues':sorted(expected_c),'top_pose_residues':sorted(observed),'recall':intersection/len(expected_c),'precision':intersection/len(observed) if observed else 0,'jaccard':intersection/len(observed|expected_c)}
  r.update(status='completed',scores=scores,top_score=scores[0],pose_RMSD_A=rms,top_RMSD_A=rms[0],best_RMSD_A=min(rms),best_RMSD_rank=rms.index(min(rms))+1,top_pose_pass_2A=rms[0]<=2,pose_count=len(mols),top_pose_contacts=contact,pdbqt=str(dest),sdf=str(sp),outputs_sha256={str(dest):sha(dest),str(sp):sha(sp)})
 except Exception as exc:r.update(status='failed',error=f'{type(exc).__name__}: {exc}')
 r['job_seconds']=time.perf_counter()-started;r['finished_utc']=datetime.now(timezone.utc).isoformat();atomic(rp,r);return r
try:
 start=time.perf_counter();results=[];params=[(g,s,e,seed) for g in geometries for s in states for e in records for seed in seeds];assert len(params)==120
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
  futures=[pool.submit(job,*param) for param in params]
  for future in concurrent.futures.as_completed(futures):
   r=future.result();r['batch_elapsed_at_completion']=time.perf_counter()-start;results.append(r)
   if len(results)%10==0 or r['status']=='failed':
    progress={'finished':len(results),'total':len(params),'completed':sum(x['status']=='completed' for x in results),'failed':sum(x['status']=='failed' for x in results),'elapsed_seconds':time.perf_counter()-start};atomic(out/'docking-progress.json',progress);print(json.dumps(progress),flush=True)
 summary={'wall_seconds':time.perf_counter()-start,'jobs':len(results),'failed':sum(r['status']=='failed' for r in results),'cache_hits':sum(r['cache_hit'] for r in results),'results':sorted(results,key=lambda r:(r['spec']['geometry'],r['spec']['receptor_state'],r['spec']['ligand_id'],r['spec']['seed']))};path=out/('docking-resume.json' if summary['cache_hits']==120 and (out/'docking.json').exists() else 'docking.json');atomic(path,summary)
 if summary['failed']:raise RuntimeError('Docking QC failures preserved; inspect outputs')
finally:lock.unlink(missing_ok=True)
