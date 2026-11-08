"""Real API rejection and isolated failed-ligand recovery checks."""
import argparse,base64,copy,json,time,urllib.error
from pathlib import Path
from ligand_library import preflight as p

def main(stage):
    library=json.loads((p.OUT/'library.json').read_text());good=next(r for r in library['records'] if r['id']=='L001-primary-1')
    if stage=='reject':
        cases=[]
        for name,body in [('duplicate',p.payload([good,good])),('empty',p.payload([])),('invalid-box',p.payload([good]))]:
            if name=='invalid-box':body['box']['size']=[-1,30,30]
            try:p.request('/vina/runs',body);raise AssertionError('Invalid request accepted: '+name)
            except urllib.error.HTTPError as exc:
                detail=json.load(exc);assert exc.code==422,(name,exc.code);cases.append({'name':name,'http_status':exc.code,'detail':detail})
        p.dump(p.OUT/'api/rejection-checks.json',cases);print('Real API rejection checks passed:',len(cases));return
    broken=p.payload([good]);broken['ligands'].append({'id':'intentional-invalid','name':'intentional-invalid.pdbqt','content_base64':base64.b64encode(b'ROOT\nATOM  INVALID_FORMAT\nENDROOT\nTORSDOF 0\n').decode()})
    rid=p.submit('mixed-invalid',broken);report=p.wait_report('mixed-invalid',rid)
    results=report['vina']['results'];assert report['status']=='partial';assert next(r for r in results if r['id']==good['id'])['status']=='succeeded';assert next(r for r in results if r['id']=='intentional-invalid')['status']=='failed'
    corrected=copy.deepcopy(good);corrected['id']='intentional-invalid'
    rid2=p.submit('corrected-input',p.payload([corrected]));report2=p.wait_report('corrected-input',rid2);assert report2['status']=='succeeded'
    p.audit_report('corrected-input',report2,[corrected])
    p.dump(p.OUT/'api/failed-input-recovery.json',{'mixed_run_id':rid,'mixed_status':report['status'],'failed_id':'intentional-invalid','unaffected_valid_id':good['id'],'corrected_run_id':rid2,'corrected_status':report2['status'],'note':'Corrected input creates a new auditable task; the original failed record remains unchanged.'})
    print('Real partial failure and corrected-input retry passed')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['reject','mixed']);main(parser.parse_args().stage)
