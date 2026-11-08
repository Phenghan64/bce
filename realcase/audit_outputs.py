"""Audit stored outputs, allowing only newly inferred protonated-N stereo labels.

Keep all original job statuses. Do not rerun docking or alter coordinates/scores.
"""
import json,re,math,hashlib
from pathlib import Path
from rdkit import Chem

def compare_identity(mol, expected):
    if Chem.MolToSmiles(mol)==Chem.MolToSmiles(expected):return True,[]
    mapping=mol.GetSubstructMatch(expected,useChirality=False)
    if len(mapping)!=expected.GetNumAtoms() or mol.GetNumAtoms()!=expected.GetNumAtoms():return False,[]
    copy=Chem.Mol(mol);cleared=[]
    for a in expected.GetAtoms():
        if a.GetAtomicNum()==7 and a.GetFormalCharge()==1 and a.GetTotalNumHs()>0 and a.GetChiralTag()==Chem.ChiralType.CHI_UNSPECIFIED:
            target=copy.GetAtomWithIdx(mapping[a.GetIdx()])
            if target.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:
                cleared.append(target.GetIdx());target.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignStereochemistry(copy,cleanIt=True,force=True)
    return Chem.MolToSmiles(copy)==Chem.MolToSmiles(expected),cleared

def main():
    out=Path(__file__).resolve().parents[1]/'results/abl-real-case'
    original=json.loads((out/'panel.json').read_text());library=json.loads((out/'library.json').read_text())
    entries={r['id']:r for r in library['records']};audited=[]
    for result in original['results']:
        r=dict(result);r['original_status']=result['status'];r['original_error']=result.get('error')
        folder=out/'jobs/panel'/f'{r["spec"]["ligand_id"]}-s{r["spec"]["seed"]}'
        try:
            assert r.get('returncode')==0,'Vina did not exit successfully'
            text=(folder/'docked.pdbqt').read_text();scores=[float(s) for s in re.findall(r'^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)',text,re.M)]
            assert scores and scores==sorted(scores) and all(math.isfinite(s) for s in scores)
            mols=list(Chem.SDMolSupplier(str(folder/'docked.sdf'),removeHs=True))
            assert len(mols)==len(scores)==len(re.findall(r'^MODEL\s+',text,re.M))
            expected=Chem.MolFromSmiles(entries[r['spec']['ligand_id']]['smiles']);notes=[]
            for i,m in enumerate(mols):
                assert m is not None
                same,cleared=compare_identity(m,expected)
                assert same,f'Pose {i}: graph/charge/specified stereo differs'
                assert all(math.isfinite(float(x)) for xyz in m.GetConformer().GetPositions() for x in xyz)
                if cleared:notes.append({'pose':i+1,'protonated_N_atom_indices':cleared})
            hashes={str(folder/name):hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ['docked.sdf','docked.pdbqt']}
            r.update(status='completed',scores=scores,top_score=scores[0],sdf=str(folder/'docked.sdf'),pdbqt=str(folder/'docked.pdbqt'),pose_count=len(scores),N_stereo_normalization=notes,outputs_sha256=hashes)
            r.pop('error',None)
        except Exception as exc:r.update(status='failed',error=str(exc))
        audited.append(r)
    report={**original,'results':audited,'completed':sum(r['status']=='completed' for r in audited),'failed':sum(r['status']=='failed' for r in audited),
      'audit_rule':'Only ignore newly inferred stereo at protonated NH+ whose source input had unspecified N stereo. All heavy-atom graphs, charges and other specified stereo remain checked. Original files/statuses preserved.',
      'audit_origin':'Post-docking QC correction after first masitinib export discrepancy; no scores, labels, selection or coordinates changed.',
      'original_panel_sha256':hashlib.sha256((out/'panel.json').read_bytes()).hexdigest(),
      'reclassified_jobs':sum(r['status']!=r['original_status'] for r in audited)}
    (out/'panel-validated.json').write_text(json.dumps(report,indent=2)+'\n')
    print({k:report[k] for k in ['completed','failed','reclassified_jobs']})

if __name__=='__main__':main()
