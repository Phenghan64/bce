"""Explicit live rollout verification. Creates three auditable Vina jobs."""
import base64, hashlib, json, time, urllib.request, urllib.error
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/ligand-library-rollout-20261030'
opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
BASE='http://100.64.0.1:8098'
def call(path,body=None):
    request=urllib.request.Request(BASE+path,data=None if body is None else json.dumps(body).encode(),headers={'Content-Type':'application/json'})
    with opener.open(request,timeout=30) as response:return json.load(response)
def save(name,value): (OUT/(name+'.json')).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
cat=call('/vina/library'); assert cat['compound_count']==20 and cat['state_count']==21
save('catalog',cat)
demo=call('/vina/demo')
body={k:demo[k] for k in ('receptor','box','params')}
body['box']={'center':[15.614,53.38,15.455],'size':[30,30,30]}
body['params'].update({'seed':20261030,'exhaustiveness':8,'cpu':1})
body['ligands']=[{'library_compound_id':c['id'],'library_version':cat['version']} for c in cat['compounds']]
cases=[]
for name,ligands in [('unknown',[{'library_compound_id':'L999','library_version':cat['version']}]),('version',[{'library_compound_id':'L001','library_version':'unknown'}]),('duplicate',[body['ligands'][0]]*2),('override',[{**body['ligands'][0],'content_base64':'forged'}])]:
    try: call('/vina/runs',{**body,'ligands':ligands});raise AssertionError(name)
    except urllib.error.HTTPError as error:
        assert error.code==422,(name,error.code)
        cases.append({'case':name,'status':error.code,'detail':json.load(error)})
save('rejection-checks',cases)
def execute(name,inputs):
    submitted=time.perf_counter(); created=call('/vina/runs',{**body,'ligands':inputs});rid=created['run_id'];save(name+'-submission',created)
    print(name,rid,flush=True)
    deadline=time.monotonic()+900
    while time.monotonic()<deadline:
        status=call('/vina/runs/'+rid)
        if status['status'] not in ('queued','running'):break
        time.sleep(1.4)
    report=call('/runs/'+rid+'/report');save(name+'-report',report)
    vina=report['vina']; metadata=vina['provenance']['ligand_library']
    for item in vina['results']:
        if item['id'] in metadata:
            assert item['input_sha256']==metadata[item['id']]['pdbqt_sha256']
            assert metadata[item['id']]['name_zh']
    audit={'run_id':rid,'status':report['status'],'api_seconds':round(time.perf_counter()-submitted,3),'worker_seconds':vina['elapsed_seconds'],'execution':vina['execution'],'result_count':len(vina['results']),'library_state_count':len(metadata),'statuses':{x['id']:x['status'] for x in vina['results']}}
    save(name+'-audit',audit);print(json.dumps(audit,ensure_ascii=False),flush=True)
    return report
report=execute('full-library-e8',body['ligands']);assert report['status']=='succeeded' and len(report['vina']['results'])==21
invalid={'id':'intentional-invalid','name':'intentional-invalid.pdbqt','content_base64':base64.b64encode(b'ROOT\nATOM  INVALID_FORMAT\nENDROOT\nTORSDOF 0\n').decode()}
report=execute('subset-upload-partial',[body['ligands'][0],body['ligands'][2],invalid]);assert report['status']=='partial' and len(report['vina']['results'])==4
report=execute('subset-corrected',[body['ligands'][0],body['ligands'][2],{**demo['ligands'][0],'id':'corrected-upload'}]);assert report['status']=='succeeded' and len(report['vina']['results'])==4
print('LIVE LIBRARY INTEGRATION PASSED',flush=True)
