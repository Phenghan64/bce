import copy,hashlib,json,os,re,shutil,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=Path(os.environ.get('STRESSBENCH_OUTPUT',ROOT/'results/stress-20261022-main'));OUT.mkdir(exist_ok=True)
def save(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
lib=json.loads((ROOT/'results/vina-throughput-1000/library.json').read_text())
groups={0:[],1:[],2:[]}
for record in lib['records']:
    if record['status']!='prepared' or record.get('role')!='unlabelled_library':continue
    record=copy.deepcopy(record);record['torsdof']=int(re.search(r'^TORSDOF\s+(\d+)',Path(record['pdbqt']).read_text(),re.M).group(1))
    group=0 if record['torsdof']<=4 else 1 if record['torsdof']<=9 else 2
    record['complexity_group']=group;groups[group].append(record)
for group in groups:groups[group].sort(key=lambda r:hashlib.sha256(r['id'].encode()).hexdigest())
selected=[groups[g][i] for i in range(32) for g in (0,1,2)]
panel=[r['id'] for r in selected[:48]]
for r in selected:
    assert sha(r['pdbqt'])==r['pdbqt_sha256']
    dest=OUT/'inputs/ligands'/Path(r['pdbqt']).name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(r['pdbqt'],dest);r['pdbqt']=str(dest)
for name in ('receptor.pdbqt','receptor.box.txt'):
    dest=OUT/'inputs'/name;shutil.copy2(ROOT/'results/vina-throughput-1000/receptor'/name,dest)
save(OUT/'inputs/library.json',{'records':selected,'calibration_panel':panel,'torsdof_bins':{'0':'0-4','1':'5-9','2':'10-14'},'source':'Existing RDKit NCI preparation; 96 ligands, balanced 32 per torsion bin; calibration uses 48 balanced ligands to actually fill all 16 workers; performance only, not activity-labelled','source_library_sha256':sha(ROOT/'results/vina-throughput-1000/library.json')})
protein=next(r for r in json.loads((ROOT/'results/marine-screen/proteins.json').read_text()) if r['id']=='A0A8C5FPM1')
sequence=protein['sequence'];assert len(sequence)==1964
save(OUT/'inputs/size-sweep-source.json',{'record':protein,'purpose':'Nested N-terminal prefixes for controlled size pressure; truncated constructs are not biological predictions'})
tm=json.loads((ROOT/'results/tmgh1-round1/query-msa-three-seeds.json').read_text())
tmseq=next(iter(tm['queries'].values()))['chains'][0]['sequence'];assert len(tmseq)==446
for n in (128,256,384,512,640,768,896,1024,1280,1536,1964):
    name=f'myosin_prefix_{n}';save(OUT/f'inputs/queries/{name}.json',{'seeds':[42],'queries':{name:{'chains':[{'molecule_type':'protein','chain_ids':['A'],'sequence':sequence[:n]}],'use_msas':False}}})
for chains in (1,2,4):
    name=f'tmgh1_copies_{chains}';save(OUT/f'inputs/queries/{name}.json',{'seeds':[42],'queries':{name:{'chains':[{'molecule_type':'protein','chain_ids':[chr(65+i) for i in range(chains)],'sequence':tmseq}],'use_msas':False}}})
source=Path(next(iter(tm['queries'].values()))['chains'][0]['main_msa_file_paths'][0]);msa=np.load(source,allow_pickle=True)['colabfold_main'].item()
for depth in (64,256,1024):
    cropped={k:v[:depth] for k,v in msa.items()};path=OUT/f'inputs/msa/tmgh1-{depth}.npz';path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,colabfold_main=cropped)
    assert ''.join(cropped['msa'][0])==tmseq
    name=f'tmgh1_msa_{depth}';save(OUT/f'inputs/queries/{name}.json',{'seeds':[42],'queries':{name:{'chains':[{'molecule_type':'protein','chain_ids':['A'],'sequence':tmseq,'main_msa_file_paths':[str(path)]}],'use_msas':True}}})
protocol={'purpose':'Measure hardware/application capacity limits, sustained degradation and recovery. Not biological validation.',
 'software':'Official OpenFold3 PyTorch/MPS, OpenBind weights; MLX is NOT installed. Vina 1.2.7.',
 'stages':{'baseline_seconds':600,'cpu_calibration':{'layouts':[[1,12],[2,6],[3,4],[6,2],[12,1],[16,1]],'exhaustiveness':[8,32,64],'minimum_seconds_per_configuration':120,'finish_whole_balanced_panel_cycles':True,'maximum_seconds_per_configuration':600},'gpu_envelope_budget_seconds':4500,'cpu_soak_seconds':5400,'mixed_seconds':2700,'recovery':'Terminate one owned Vina worker, retry same job; rerun fixed plan to verify hashes and cache resume.'},
 'cpu_fairness':'Same 48 balanced ligands and seed across calibration configurations; complete panel cycles. Soak and mixed use the same 96-ligand panel. Report peak actual active workers, censoring/failed work and per-complexity throughput. Earlier 12-ligand pilot is excluded from comparisons.',
 'gpu_method':'Same natural myosin sequence prefixes control token length. Full native TmGH1 tests MSA depth and diffusion samples. 2/4 chain copies are workload constructs, not a claim of native oligomerization. No network MSA, no templates; original model cycles/steps unless explicitly recorded.',
 'gpu_budget':'MPS allocator budget 1.0 x recommended max, compared to previous 0.7. Retain low_mem preset. Budget OOM is an operational cap, not proof of exhausting every physical byte.',
 'guardrails':{'system_available_min_MiB':384,'low_memory_grace_seconds':12,'swap_absolute_max_GiB':8,'swap_growth_per_stage_max_GiB':6,'thermal_critical_grace_seconds':30,'disk_free_min_GiB':20,'test_process_groups_only':True},
 'adaptation':'After first size failure, test midpoint where possible and repeat boundary once; do not continue arbitrarily larger failures. GPU remaining budget is spent repeating largest stable input. Mixed falls back if simultaneous pressure guard triggers. All fallback and aborted jobs remain evidence.',
 'thermal':'NSProcessInfo thermalState and ioreg GPU statistics sampled. Temperature/power unavailable without elevation; no thermal-throttling claim solely from slower runtime.',
 'reporting':'Separate task wall time, throughput, completion/QC rates, tail latency, CPU per-core use, GPU use, RSS, Metal memory, available RAM, swap and thermal pressure. Report load achieved, not merely requested.'}
save(OUT/'protocol.json',protocol)
files=list((OUT/'inputs').rglob('*'))
save(OUT/'input-manifest.json',[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files if p.is_file()])
print('Prepared',len(selected),'ligands,',len(panel),'calibration ligands; GPU input ladder and cached MSA depths.')
