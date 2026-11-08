"""Analyze predeclared structure metrics without choosing a favorable prediction."""
import json,itertools
from pathlib import Path
import numpy as np
from biotite.structure.io import pdb,pdbx
from rdkit import Chem
root=Path(__file__).resolve().parents[1];out=root/'results/tmgh1-round1'
def fit(x,y):
 xm=x.mean(0);ym=y.mean(0);u,s,vt=np.linalg.svd((x-xm).T@(y-ym));r=u@np.diag([1,1,np.linalg.det(u@vt)])@vt;return r,ym-xm@r
def rms(x,y):return float(np.sqrt(np.mean(np.sum((np.asarray(x)-np.asarray(y))**2,axis=1))))
def coords(a):return {(int(x.res_id),str(x.atom_name)):x.coord.astype(float) for x in a}
def dihedral(v):
 b0=v[0]-v[1];b1=v[2]-v[1];b2=v[3]-v[2];b1/=np.linalg.norm(b1);x=b0-b1*np.dot(b0,b1);y=b2-b1*np.dot(b2,b1);return float(np.degrees(np.arctan2(np.dot(np.cross(b1,x),y),np.dot(x,y))))
chi={'SER':'OG','THR':'OG1','CYS':'SG','VAL':'CG1','ILE':'CG1'}
refs={code:pdb.PDBFile.read(out/'geometries'/f'experimental_{code}'/'observed.pdb').get_structure(model=1) for code in ['2J77','2J7E']};rbase=coords(refs['2J77']);baseids=sorted(i for i,n in rbase if n=='CA');reports=[]
for seed in [42,2026,20261021]:
 source=out/'inference/prediction/tmgh1_native'/f'seed_{seed}'/f'tmgh1_native_seed_{seed}_sample_1_model.cif';assert source.exists(),source
 a=pdbx.get_structure(pdbx.CIFFile.read(source),model=1,use_author_fields=True,extra_fields=['b_factor']);a=a[(a.chain_id=='A')&~a.hetero&~np.isin(a.element,['H','D'])].copy();assert len(set(a.res_id))==446
 pc=coords(a);r,t=fit(np.array([pc[i,'CA'] for i in baseids]),np.array([rbase[i,'CA'] for i in baseids]));a.coord=a.coord@r+t;d=out/'geometries'/f'predicted_{seed}';d.mkdir(parents=True,exist_ok=True);pf=pdb.PDBFile();pf.set_structure(a);pf.write(d/'observed.pdb');pc=coords(a)
 agg=json.loads(source.with_name(source.name.replace('_model.cif','_confidences_aggregated.json')).read_text());timing=json.loads((source.parent/'timing.json').read_text());conf={int(i):float(a.b_factor[a.res_id==i].mean()) for i in set(a.res_id)};res_names={int(x.res_id):str(x.res_name) for x in a}
 row={'seed':seed,'native_residues':446,'source':str(source),'model_runtime_seconds':timing['runtime_s'],'confidence':agg,'rotation_to_2J77':r.tolist(),'translation_to_2J77':t.tolist(),'comparisons':{}}
 for code,ref in refs.items():
  rc=coords(ref);ids=sorted(i for i,n in rc if n=='CA');pr=np.array([pc[i,'CA'] for i in ids]);rr=np.array([rc[i,'CA'] for i in ids]);rot,tr=fit(pr,rr);aligned={k:v@rot+tr for k,v in pc.items()};meta=json.loads((out/'geometries'/f'experimental_{code}'/'metadata.json').read_text());pocket=meta['pocket_residues_6A'];ca=[(i,'CA') for i in pocket if (i,'CA') in rc];back=[k for k in rc if k[0] in pocket and k[1] in ['N','CA','C','O']];heavy=[k for k in rc if k[0] in pocket and k in aligned]
  errors={str(i):float(np.linalg.norm(aligned[i,'CA']-rc[i,'CA'])) for i in ids};rows=[];sum_sq=0;count=0
  for i in pocket:
   keys=[k for k in heavy if k[0]==i];names=[k[1] for k in keys];xx=np.array([aligned[k] for k in keys]);base=np.array([rc[k] for k in keys]);name=res_names[i];swaps=[]
   if name=='ASP':swaps=[('OD1','OD2')]
   if name=='GLU':swaps=[('OE1','OE2')]
   if name=='ARG':swaps=[('NH1','NH2')]
   if name in ['PHE','TYR']:swaps=[('CD1','CD2'),('CE1','CE2')]
   alt=base.copy()
   for n1,n2 in swaps:
    if n1 in names and n2 in names:alt[names.index(n1)],alt[names.index(n2)]=base[names.index(n2)],base[names.index(n1)]
   yy=alt if np.sum((xx-alt)**2)<np.sum((xx-base)**2) else base;sum_sq+=float(np.sum((xx-yy)**2));count+=len(keys)
   sc=[j for j,n in enumerate(names) if n not in ['N','CA','C','O','OXT']];chi_names=['N','CA','CB',chi.get(name,'CG')];chierr=None
   if all((i,n) in aligned and (i,n) in rc for n in chi_names):
    xchi=dihedral(np.array([aligned[i,n] for n in chi_names]));ychi=dihedral(np.array([rc[i,n] for n in chi_names]));chierr=abs((xchi-ychi+180)%360-180)
   rows.append({'residue':i,'name':name,'pLDDT':conf[i],'CA_deviation_A':errors[str(i)],'sidechain_symmetry_corrected_RMSD_A':rms(xx[sc],yy[sc]) if sc else None,'chi1_error_degrees':chierr})
  xp=np.array([aligned[k] for k in back]);yp=np.array([rc[k] for k in back]);r2,t2=fit(xp,yp);catalytic={str(i):{'carboxylate_centroid_shift_A':float(np.linalg.norm(np.mean([aligned[i,n] for n in ['OE1','OE2']],axis=0)-np.mean([rc[i,n] for n in ['OE1','OE2']],axis=0))),'pLDDT':conf[i]} for i in [166,351]}
  ligand=Chem.SDMolSupplier(str(out/'reference'/f'{meta["ligand"]}.sdf'))[0].GetConformer().GetPositions();loop=np.array([aligned[k] for k in aligned if k[0] in [1,2,307,446]]);loopdist=float(np.linalg.norm(loop[:,None]-ligand[None,:],axis=2).min())
  comp={'common_CA_count':len(ids),'global_CA_RMSD_A':rms([aligned[i,'CA'] for i in ids],rr),'pocket_residues_6A':pocket,'pocket_CA_RMSD_A':rms([aligned[k] for k in ca],[rc[k] for k in ca]),'pocket_backbone_RMSD_A_global_alignment':rms(xp,yp),'pocket_backbone_RMSD_A_local_alignment':rms(xp@r2+t2,yp),'pocket_heavy_RMSD_A_atom_name_mapping':rms([aligned[k] for k in heavy],[rc[k] for k in heavy]),'pocket_heavy_RMSD_A_symmetry_corrected':float(np.sqrt(sum_sq/count)),'pocket_mean_pLDDT':float(np.mean([conf[i] for i in pocket])),'catalytic':catalytic,'unobserved_reference_residues_predicted_min_ligand_distance_A':loopdist,'per_residue_CA_errors_A':errors,'pocket_residue_metrics':rows,'pocket_sidechain_missing_reference_atoms':'Metrics use only experimentally observed matched atoms; no modeled reference sidechains included'}
  row['comparisons'][code]=comp
 (d/'metadata.json').write_text(json.dumps(row,indent=2)+'\n');reports.append(row)
(out/'structure-analysis.json').write_text(json.dumps({'primary_seed':42,'models':reports,'note':'All observed native CA used without trimming; separate protein alignment for each reference; no ligand input/templates.'},indent=2)+'\n')
for r in reports:print(r['seed'],r['confidence']['avg_plddt'],{code:{k:v for k,v in c.items() if k in ['global_CA_RMSD_A','pocket_backbone_RMSD_A_global_alignment','pocket_heavy_RMSD_A_symmetry_corrected','catalytic']} for code,c in r['comparisons'].items()})
