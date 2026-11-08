"""Traceable preparation and production-API preflight for the user panel.

This is a validation harness, not a production pKa predictor or library UI.
"""
from __future__ import annotations
import argparse,base64,hashlib,importlib.metadata,json,math,sys,time,urllib.error,urllib.request,urllib.parse
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,Lipinski,rdMolAlign
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers,StereoEnumerationOptions
from meeko import MoleculePreparation,PDBQTWriterLegacy,PDBQTMolecule,RDKitMolCreate

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/ligand-library20-preflight-20261030'
API='http://100.64.0.1:8098'
OPENER=urllib.request.build_opener(urllib.request.ProxyHandler({}))
SEED=20261030

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n');temp.replace(path)
def canonical(mol):return Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)
def request(path,body=None):
    req=urllib.request.Request(API+path,data=json.dumps(body).encode() if body is not None else None,headers={'Content-Type':'application/json'})
    with OPENER.open(req,timeout=45) as r:return json.load(r)
def exact_identity(mol,expected):
    if mol is None or canonical(mol)!=canonical(expected):raise ValueError('Molecular graph, charge or specified stereochemistry changed')
    xyz=mol.GetConformer().GetPositions()
    if not np.isfinite(xyz).all():raise ValueError('Nonfinite coordinates')
def decode_pdbqt(text):
    sdf,failures=RDKitMolCreate.write_sd_string(PDBQTMolecule(text,skip_typing=True))
    if failures:raise ValueError('Meeko export failures: '+str(failures))
    mols=list(Chem.ForwardSDMolSupplier(__import__('io').BytesIO(sdf.encode()),removeHs=True))
    if not mols or any(m is None for m in mols):raise ValueError('Unreadable SDF poses')
    return mols,sdf

def variants(mol):
    return sorted({canonical(m):m for m in EnumerateStereoisomers(mol,options=StereoEnumerationOptions(onlyUnassigned=True,unique=True,maxIsomers=16))}.items())
def carboxylate(mol):
    copy=Chem.RWMol(mol);query=Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]');changed=[]
    for match in mol.GetSubstructMatches(query):
        atom=copy.GetAtomWithIdx(match[2]);atom.SetFormalCharge(-1);atom.SetNumExplicitHs(0);atom.SetNoImplicit(True);changed.append(match[2])
    result=copy.GetMol();Chem.SanitizeMol(result);return result,changed

def prepare_one(entry,folder):
    folder.mkdir(parents=True,exist_ok=True);start=time.perf_counter();expected=Chem.MolFromSmiles(entry['smiles'])
    if expected is None or len(Chem.GetMolFrags(expected))!=1:raise ValueError('Invalid or disconnected ligand')
    if any(str(i.specified) in ('Unspecified','Unknown') for i in Chem.FindPotentialStereo(expected)):
        raise ValueError('Unresolved stereochemistry requires explicit enumeration or review')
    mol=Chem.AddHs(expected);params=AllChem.ETKDGv3();params.randomSeed=SEED;params.numThreads=1;params.enforceChirality=True;params.pruneRmsThresh=.35
    ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=16,params=params))
    if not ids:raise ValueError('ETKDG failed')
    if not AllChem.MMFFHasAllMoleculeParams(mol):raise ValueError('Incomplete MMFF parameters')
    energies=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=2000,mmffVariant='MMFF94s')
    accepted=[(energy,int(cid)) for cid,(status,energy) in zip(ids,energies) if status==0 and math.isfinite(energy)]
    if not accepted:raise ValueError('No MMFF conformer converged')
    energy,cid=min(accepted)
    result=Chem.Mol(mol);result.RemoveAllConformers();result.AddConformer(Chem.Conformer(mol.GetConformer(cid)),assignId=True);result.SetProp('_Name',entry['id'])
    sdf=folder/'input.sdf'
    with Chem.SDWriter(str(sdf)) as w:w.write(result)
    exact_identity(Chem.SDMolSupplier(str(sdf),removeHs=True)[0],expected)
    setups=MoleculePreparation().prepare(result)
    if len(setups)!=1:raise ValueError('Expected one Meeko setup per explicit state')
    text,valid,error=PDBQTWriterLegacy.write_string(setups[0])
    if not valid:raise ValueError(error)
    decoded,_=decode_pdbqt(text);exact_identity(decoded[0],expected)
    atoms=[line for line in text.splitlines() if line.startswith(('ATOM  ','HETATM'))]
    charges=[float(line[70:76]) for line in atoms]
    if not charges or not all(math.isfinite(x) for x in charges) or abs(sum(charges)-Chem.GetFormalCharge(expected))>.05:
        raise ValueError('PDBQT charge sum mismatch')
    dest=folder/'ligand.pdbqt';dest.write_text(text)
    return {**entry,'status':'prepared','pdbqt':str(dest.relative_to(ROOT)),'sdf':str(sdf.relative_to(ROOT)),
            'pdbqt_sha256':digest(dest),'sdf_sha256':digest(sdf),'heavy_atoms':expected.GetNumHeavyAtoms(),
            'formal_charge':Chem.GetFormalCharge(expected),'rdkit_rotatable_bonds':Lipinski.NumRotatableBonds(expected),
            'vina_torsdof':int(text.split('TORSDOF ')[-1].strip()),'embedded_conformers':len(ids),'converged_conformers':len(accepted),
            'selected_mmff_energy':energy,'preparation_seconds':round(time.perf_counter()-start,4)}

def prepare():
    if (OUT/'library.json').exists():raise RuntimeError('Prepared library exists; use audit/resume, do not overwrite')
    props=json.loads((OUT/'source/properties.json').read_text())['PropertyTable']['Properties'];bycid={r['CID']:r for r in props}
    panel=json.loads((OUT/'requested-panel.json').read_text());source_mols={int(m.GetProp('PUBCHEM_COMPOUND_CID')):m for m in Chem.SDMolSupplier(str(OUT/'source/structures-2d.sdf'),removeHs=True)}
    protocol={'frozen_at':datetime.now(timezone.utc).isoformat(),'seed':SEED,'preparation':'ETKDGv3 16 conformers, MMFF94s maxIters2000, lowest converged energy, Meeko Gasteiger PDBQT; strict identity round trips',
        'state_policy':'Primary screening reference: only carboxylic acids converted to carboxylates; phenols retained as in source. Preserve source stereochemistry, explicitly enumerate racemic ibuprofen R/S. Acid neutral source states retained as alternatives. Two intact-aromatic curcumin enol geometries retained as alternatives. This is an explicit engineering reference ensemble, not a dominant-state/pKa prediction.',
        'review_required':'Polyphenol protonation and curcumin tautomer populations require chemistry review for target assay pH; no best docking score is used to choose chemical states.',
        'receptor':'Existing official Vina 1IEP experimental receptor, imatinib site; engineering panel, no activity labels for the 20 compounds',
        'dock':{'box':{'center':[15.614,53.380,15.455],'size':[30,30,30]},'e8':'all explicit states, seed20261030, CPU1, worker default8','e32':'preselected L001,L003,L007,L013,L015,L020 primary states; same seed, CPU1','calibration':'independently prepared imatinib, E32 seeds42,2026,20261030; fixed-frame symmetry-aware heavy-atom RMSD top pose<=2A'},
        'production_api_max_ligands':100,'sources':json.loads((OUT/'source/retrieval.json').read_text())}
    dump(OUT/'protocol.json',protocol)
    sys.path.insert(0,str(ROOT/'.tools/realcase-deps'));from dimorphite_dl import protonate_smiles
    diagnostics=[];entries=[];t=time.perf_counter()
    for p in panel:
        row=bycid[p['cid']];raw=Chem.MolFromSmiles(row['SMILES']);sm=source_mols[p['cid']]
        assert canonical(raw)==canonical(sm) and Chem.MolToInchiKey(sm)==row['InChIKey'],'PubChem source mismatch'
        nominal=protonate_smiles(row['SMILES'],ph_min=7.4,ph_max=7.4,precision=0,max_variants=512)
        ranged=protonate_smiles(row['SMILES'],ph_min=6.9,ph_max=7.9,precision=1,max_variants=512)
        diagnostics.append({**p,'nominal_dimorphite_charge':[Chem.GetFormalCharge(Chem.MolFromSmiles(s)) for s in nominal],'range_state_count':len(ranged),'nominal_smiles':nominal,'range_smiles':ranged,'note':'Diagnostics only; not silently accepted as dominant chemical states'})
        primary,changed=carboxylate(raw)
        for j,(s,mol) in enumerate(variants(primary),1):
            entries.append({**p,'id':f'{p["id"]}-primary-{j}','smiles':s,'role':'primary_reference','carboxylate_atoms':changed,'source_inchikey':row['InChIKey'],'source_smiles':row['SMILES'],'chemistry_review_required':p['cid'] in {1983,5280343,5280863,5280443,5280445,5280961,5281708,969516,445154,65064,9064,370,445858,689043,1794427}})
        if changed:
            for j,(s,mol) in enumerate(variants(raw),1):entries.append({**p,'id':f'{p["id"]}-source-neutral-{j}','smiles':s,'role':'source_state_alternative','chemistry_review_required':True})
        if p['cid']==969516:
            params=rdMolStandardize.CleanupParameters();params.tautomerRemoveBondStereo=False;params.tautomerRemoveSp3Stereo=False
            all_tautomers=rdMolStandardize.TautomerEnumerator(params).Enumerate(raw)
            selected=[m for m in all_tautomers if sum(a.GetIsAromatic() for a in m.GetAtoms())==12 and canonical(m)!=canonical(raw)]
            enols=sorted({s:m for m in selected for s,m in variants(m)}.items())
            diagnostics[-1].update(raw_tautomers=len(all_tautomers),intact_aromatic_enol_variants=len(enols))
            for j,(s,mol) in enumerate(enols,1):entries.append({**p,'id':f'{p["id"]}-enol-{j}','smiles':s,'role':'tautomer_alternative','chemistry_review_required':True})
    case=json.loads((ROOT/'results/closed-loop-1iep/case.json').read_text())
    entries.append({'id':'imatinib-control','cid':5291,'name_zh':'伊马替尼（独立姿态对照）','smiles':case['smiles'],'role':'calibration','chemistry_review_required':False})
    records=[]
    for entry in entries:
        try:r=prepare_one(entry,OUT/'prepared'/entry['id'])
        except Exception as exc:r={**entry,'status':'failed','error':f'{type(exc).__name__}: {exc}'}
        records.append(r);print(entry['id'],r['status'],r.get('error',''),flush=True)
    dump(OUT/'state-diagnostics.json',{'tool':'Dimorphite-DL '+importlib.metadata.version('dimorphite-dl'),'rules_sha256':digest(ROOT/'.tools/realcase-deps/dimorphite_dl/smarts/site_substructures.smarts'),'records':diagnostics,'decision':'Do not adopt the precision=0 nominal-state policy as a general pH preparation default for this polyphenol-rich panel.'})
    dump(OUT/'library.json',{'protocol_sha256':digest(OUT/'protocol.json'),'source_cids':len(panel),'preparation_seconds':round(time.perf_counter()-t,3),'versions':{n:importlib.metadata.version(n) for n in ('rdkit','meeko','dimorphite-dl')},'records':records})
    assert all(r['status']=='prepared' for r in records),'Some states failed preparation; inspect library.json'

def audit():
    lib=json.loads((OUT/'library.json').read_text());checked=0
    for r in lib['records']:
        if r['status']!='prepared':continue
        for key in ('pdbqt','sdf'):
            if digest(ROOT/r[key])!=r[key+'_sha256']:raise ValueError('Cached artifact checksum mismatch: '+r['id'])
        mols,_=decode_pdbqt((ROOT/r['pdbqt']).read_text());exact_identity(mols[0],Chem.MolFromSmiles(r['smiles']));checked+=1
    return checked

def payload(records,exhaust=8,seed=SEED):
    protocol=json.loads((OUT/'protocol.json').read_text());raw=(ROOT/'results/1iep/receptor.pdbqt').read_bytes()
    return {'receptor':{'name':'1IEP-experimental-receptor.pdbqt','content_base64':base64.b64encode(raw).decode()},
        'ligands':[{'id':r['id'],'name':r['id']+'.pdbqt','content_base64':base64.b64encode((ROOT/r['pdbqt']).read_bytes()).decode()} for r in records],
        'box':protocol['dock']['box'],'params':{'exhaustiveness':exhaust,'num_modes':9,'energy_range':4,'seed':seed,'cpu':1}}

def submit(name,body):
    saved=OUT/'api'/f'{name}-submission.json'
    if saved.exists():
        old=json.loads(saved.read_text())
        signature=hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest()
        if old['signature']!=signature:raise ValueError('Resume request differs: '+name)
        print('Resuming existing API task:',name,old['run_id'],flush=True);return old['run_id']
    response=request('/vina/runs',body)
    dump(saved,{**response,'signature':hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest()})
    print('Submitted:',name,response['run_id'],flush=True);return response['run_id']

def wait_report(name,run_id):
    last=None;deadline=time.monotonic()+1800
    while time.monotonic()<deadline:
        r=request('/vina/runs/'+run_id)
        if r['current_stage']!=last:print(name,r['status'],r['current_stage'],flush=True);last=r['current_stage']
        if r['status'] in ('succeeded','partial','failed'):
            if r['status']=='failed':raise RuntimeError(r.get('error'))
            report=request('/runs/'+run_id+'/report');dump(OUT/'api'/f'{name}-report.json',report);return report
        time.sleep(2)
    raise TimeoutError('Run still active; preserve submission and resume later')

def audit_report(name,report,records):
    refs={r['id']:r for r in records};vina=report['vina'];rows=[]
    receptor_xyz=np.asarray([[float(line[30:38]),float(line[38:46]),float(line[46:54])] for line in (ROOT/'results/1iep/receptor.pdbqt').read_text().splitlines() if line.startswith(('ATOM  ','HETATM')) and line[77:].strip() not in ('H','HD','HS')])
    center=np.asarray(vina['configuration']['box']['center']);half=np.asarray(vina['configuration']['box']['size'])/2
    for result in vina['results']:
        row={'id':result['id'],'status':result['status'],'elapsed_seconds':result['elapsed_seconds'],'scores':result.get('scores',[])}
        if result['status']=='succeeded':
            path=result['output'].lstrip('/');url='/runs/'+report['run_id']+'/artifact/'+urllib.parse.quote(path,safe='/')
            text=OPENER.open(API+url,timeout=30).read().decode();dest=OUT/'poses'/name/(result['id']+'.pdbqt');dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text)
            mols,sdf=decode_pdbqt(text);dest.with_suffix('.sdf').write_text(sdf)
            assert len(mols)==len(result['scores'])
            expected=Chem.MolFromSmiles(refs[result['id']]['smiles'])
            minimum_distances=[];outside=0
            for mol in mols:
                exact_identity(mol,expected)
                xyz=mol.GetConformer().GetPositions()
                outside+=int(np.any(np.abs(xyz-center)>half+.1))
                minimum_distances.append(float(np.linalg.norm(xyz[:,None]-receptor_xyz[None],axis=2).min()))
            values=[s['affinity_kcal_mol'] for s in result['scores']]
            assert np.isfinite(values).all() and values==sorted(values)
            row['poses_outside_declared_box']=outside
            row['min_heavy_distance_A_by_pose']=minimum_distances
            row['geometric_review_flags']=[]
            if outside:row['geometric_review_flags'].append('One or more poses extend outside the declared search box')
            if min(minimum_distances)<1.5:row['geometric_review_flags'].append('One or more poses have protein-ligand heavy atoms closer than 1.5A')
            row['all_poses_chemistry_checked']=len(mols);row['top_score']=result['scores'][0]['affinity_kcal_mol']
            if result['id']=='imatinib-control':
                ref=Chem.SDMolSupplier(str(ROOT/'results/closed-loop-1iep/reference-imatinib.sdf'),removeHs=True)[0]
                row['crystal_top_pose_rmsd_A']=float(rdMolAlign.CalcRMS(mols[0],ref))
        rows.append(row)
    dump(OUT/'api'/f'{name}-audit.json',{'run_id':report['run_id'],'status':report['status'],'api_runtime_seconds':report['runtime_seconds'],'execution':vina['execution'],'results':rows})
    return rows

def dock(stage):
    audit();lib=json.loads((OUT/'library.json').read_text());records=[r for r in lib['records'] if r['status']=='prepared'];primary=[r for r in records if r['role']=='primary_reference'];control=[r for r in records if r['role']=='calibration']
    if stage=='calibrate':
        for seed in (42,2026,SEED):
            name=f'calibration-e32-s{seed}';body=payload(control,32,seed);rid=submit(name,body);report=wait_report(name,rid);rows=audit_report(name,report,control)
            assert rows[0]['crystal_top_pose_rmsd_A']<=2,'Known-pose calibration failed; stop panel assessment'
        return
    if stage=='e8':selected=records;exhaust=8
    elif stage=='e32':selected=[r for r in primary if r['id'].split('-')[0] in {'L001','L003','L007','L013','L015','L020'}]+control;exhaust=32
    else:raise ValueError(stage)
    name='panel-'+stage;rid=submit(name,payload(selected,exhaust));report=wait_report(name,rid);audit_report(name,report,selected)

if __name__=='__main__':
    import urllib.parse
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','audit','calibrate','e8','e32']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='audit':print('Cached artifacts verified:',audit())
    else:dock(a.stage)
