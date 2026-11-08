"""Report measured throughput, resource use, failure accounting and pose controls."""
import hashlib,itertools,json,re,statistics
from collections import Counter
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

root=Path(__file__).resolve().parents[1];out=root/'results/vina-throughput-1000'
read=lambda f:json.loads((out/f).read_text())
library=read('library.json');entries={r['id']:r for r in library['records']}
stages={name:read(name+'.json') for name in ['calibrate','screen','refine']}
screen=stages['screen'];refine=stages['refine'];resources=read('resources.json');resume=read('resume-verification.json')
all_jobs=[r for s in stages.values() for r in s['results']]
for r in all_jobs:
    for key in ['receptor','box']:
        assert hashlib.sha256(Path(r['spec'][key]).read_bytes()).hexdigest()==r['spec'][key+'_sha256']
    assert hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest()==r['spec']['protocol_sha256']
    assert hashlib.sha256((out/'library.json').read_bytes()).hexdigest()==r['spec']['library_sha256']
    assert hashlib.sha256(Path(entries[r['spec']['ligand_id']]['pdbqt']).read_bytes()).hexdigest()==r['spec']['ligand_sha256']
    if r['status']=='completed':
        for path,sha in r['outputs_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
ranking=[]
for rank,r in enumerate(sorted((x for x in screen['results'] if x['status']=='completed'),key=lambda x:(x['top_score'],x['spec']['ligand_id'])),1):
    e=entries[r['spec']['ligand_id']]
    ranking.append({'rank':rank,'id':e['id'],'score':r['top_score'],'role':e['role'],'heavy_atoms':e['heavy_atoms'],
                    'rotatable_bonds':e['rdkit_rotatable_bonds'],'molecular_weight':e['molecular_weight'],'smiles':e['smiles'],'sdf':r['sdf']})
(out/'screen-ranking.json').write_text(json.dumps(ranking,indent=2)+'\n')
with Chem.SDWriter(str(out/'screen-top-poses.sdf')) as writer:
    for row in ranking:
        mol=Chem.SDMolSupplier(row['sdf'],removeHs=True)[0];mol.SetProp('_Name',row['id']);mol.SetProp('Vina_score',str(row['score']));writer.write(mol)
ref_rows=[]
for cid in read('refinement-selection.json')['ids']:
    runs=[r for r in refine['results'] if r['spec']['ligand_id']==cid and r['status']=='completed']
    mols=[Chem.SDMolSupplier(r['sdf'],removeHs=True)[0] for r in runs]
    rmsds=[float(rdMolAlign.CalcRMS(a,b)) for a,b in itertools.combinations(mols,2)]
    scores=[r['top_score'] for r in runs]
    ref_rows.append({'id':cid,'successful_seeds':len(runs),'median_score':statistics.median(scores) if scores else None,
        'score_range':max(scores)-min(scores) if scores else None,'max_pose_rmsd_A':max(rmsds) if rmsds else None,
        'consistent_within2A':len(runs)==3 and max(rmsds)<=2})
ref_rows.sort(key=lambda x:(x['median_score'] is None,x['median_score'] or 0,x['id']))
(out/'refine-ranking.json').write_text(json.dumps(ref_rows,indent=2)+'\n')
reference=Chem.SDMolSupplier(str(out/'receptor/crystal-imatinib.sdf'),removeHs=True)[0]
controls={name:[{'seed':r['spec']['seed'],'score':r['top_score'],'rmsd_A':float(rdMolAlign.CalcRMS(Chem.SDMolSupplier(r['sdf'],removeHs=True)[0],reference))}
    for r in s['results'] if r['spec']['ligand_id']=='imatinib-control' and r['status']=='completed'] for name,s in stages.items()}
ordered=sorted(screen['results'],key=lambda r:r['stage_elapsed_at_completion']);blocks=[];previous=0
for i in range(0,len(ordered),100):
    batch=ordered[i:i+100];end=batch[-1]['stage_elapsed_at_completion'];duration=end-previous
    blocks.append({'completion_start':i+1,'completion_end':i+len(batch),'seconds':duration,'submitted_per_hour':len(batch)/duration*3600,
        'qc_pass':sum(r['status']=='completed' for r in batch),'mean_heavy_atoms':statistics.mean(entries[r['spec']['ligand_id']]['heavy_atoms'] for r in batch),
        'mean_rotatable_bonds':statistics.mean(entries[r['spec']['ligand_id']]['rdkit_rotatable_bonds'] for r in batch)})
    previous=end
inputs=[e for e in entries.values() if e['role']=='unlabelled_library'];failures=[e for e in inputs if e['status']=='failed']
job_times=[r['job_wall_seconds'] for r in screen['results']]
failure_counts=Counter(r.get('error','unknown') for r in all_jobs if r['status']!='completed')
summary={'scope':'1000 selected unlabelled NCI library compounds plus imatinib; human ABL1 2HYY B; CPU Vina1.2.7',
    'selected_library_count':len(inputs),'prepared_library_count':sum(e['status']=='prepared' for e in inputs),'preparation_failures':failures,
    'preparation_seconds':library['preparation_wall_seconds'],'stages':{name:{k:s[k] for k in ['wall_seconds','completed','failed','cache_hits']} for name,s in stages.items()},
    'vina_jobs':len(all_jobs),'vina_nonzero_or_missing_exit':sum(r.get('returncode')!=0 for r in all_jobs),
    'output_qc_rejections':sum(r.get('returncode')==0 and r['status']!='completed' for r in all_jobs),'failure_reasons':dict(failure_counts),
    'validated_poses':sum(r.get('pose_count',0) for r in all_jobs),'screen_submitted_per_hour':len(screen['results'])/screen['wall_seconds']*3600,
    'screen_qc_pass_per_hour':screen['completed']/screen['wall_seconds']*3600,'screen_job_seconds':{'median':float(np.median(job_times)),'p95':float(np.percentile(job_times,95)),'max':max(job_times)},
    'mean_library_heavy_atoms':statistics.mean(e['heavy_atoms'] for e in inputs),'mean_library_rotatable_bonds':statistics.mean(e['rdkit_rotatable_bonds'] for e in inputs),
    'resources':resources,'total_preparation_and_monitored_run_seconds':library['preparation_wall_seconds']+resources['seconds'],
    'blocks_of100_completions':blocks,'imatinib_controls':controls,'refine_pose_consistent_count':sum(r['consistent_within2A'] for r in ref_rows if r['id']!='imatinib-control'),
    'refined_library_count':len(ref_rows)-1,'resume':resume,'validated_input_output_hashes':True,'openfold3_predictions':0,
    'output_file_bytes':sum(f.stat().st_size for f in out.rglob('*') if f.is_file()),
    'limitations':read('protocol.json')['limitations']+['Source protonation/tautomer retained; no experimental affinity validation.','A100-completion block is not a fixed equal-complexity benchmark; variation does not by itself establish thermal throttling.','Timing excludes existing receptor preparation, software development and analysis/resume verification.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['千级 AutoDock Vina 高通量流程实测','',
    f"抽取库分子{len(inputs)}个，准备成功{summary['prepared_library_count']}个，失败{len(failures)}个；另加伊马替尼独立对照。",
    '受体：人源ABL1晶体结构2HYY B链；复用已验证的受体；搜索盒30×30×30Å。',
    '设置：3个并发任务×每任务2线程；初筛E8；前20名＋对照E32、三个种子复筛。',
    f"Vina实际调用{len(all_jobs)}次，非零/缺失退出码{summary['vina_nonzero_or_missing_exit']}次，输出质检拦截{summary['output_qc_rejections']}次。",
    f"已通过质检姿势{summary['validated_poses']}个。所有成功任务的输入/输出哈希复核通过。",'',
    f"准备耗时：{library['preparation_wall_seconds']:.2f}秒。",
    *[f"{name}耗时：{s['wall_seconds']:.2f}秒；质检通过{s['completed']}、失败{s['failed']}。" for name,s in stages.items()],
    f"准备＋完整监测运行：{summary['total_preparation_and_monitored_run_seconds']:.2f}秒，约{summary['total_preparation_and_monitored_run_seconds']/60:.2f}分钟。",
    f"初筛实测速率：提交{summary['screen_submitted_per_hour']:.1f}个/小时；质检通过{summary['screen_qc_pass_per_hour']:.1f}个/小时。",
    f"初筛单任务耗时：中位{np.median(job_times):.2f}秒，P95={np.percentile(job_times,95):.2f}秒，最长{max(job_times):.2f}秒；并发墙钟时间另计。",'',
    f"进程树RSS采样峰值：{resources['process_tree_peak_rss_bytes']/2**30:.2f}GiB；全系统可用内存最低：{resources['system_min_available_bytes']/2**30:.2f}GiB。",
    f"系统swap使用量起止变化：{(resources['swap_after']['used']-resources['swap_before']['used'])/2**20:.2f}MiB；峰值相对开始增加{(resources['swap_peak_used_bytes']-resources['swap_before']['used'])/2**20:.2f}MiB。",
    'RSS相加可能重复计入共享页；系统内存/swap还包含其他应用，不能据此判定全部负载来源。','',
    '时间口径：不含已存在的受体准备、开发脚本、结果分析与续跑验证；包含本轮配体准备和完整校准/初筛/复筛监测运行。','',
    f"20个复筛库分子中，{summary['refine_pose_consistent_count']}个的三种子最高分姿势两两RMSD均≤2Å；这是重复一致性，不是活性命中率。",
    f"伊马替尼E8初筛最高分姿势RMSD：{[round(x['rmsd_A'],3) for x in controls['screen']]}Å；E32三种子复筛RMSD：{[round(x['rmsd_A'],3) for x in controls['refine']]}Å。",
    '即使任务正常退出并通过化学结构质检，也不保证姿势接近实验；E8初筛需要后续复筛和科学验证。',
    '伊马替尼对照：'+json.dumps(controls,ensure_ascii=False),
    '准备失败明细：'+json.dumps(failures,ensure_ascii=False),
    '输出失败原因：'+json.dumps(dict(failure_counts),ensure_ascii=False),
    '续跑验证：'+json.dumps(resume,ensure_ascii=False),'',
    f"库分子平均{summary['mean_library_heavy_atoms']:.2f}个重原子、{summary['mean_library_rotatable_bonds']:.2f}个可旋转键。",
    '本轮是千级库的实际持续运行，不是万级/十万级实测；不能把本次速率直接用于更大、更柔性的药物库。',
    '库中没有本靶点的实验活性标签，不计算AUC、富集或真实命中率；Vina分数不等同于Kd。',
    '没有重新运行OpenFold3。完整日志、构象与任务记录在jobs/，资源时序在resources.jsonl。',
    '保留所有失败，不替换分子、不自动重试。续跑验证针对已完成结果，不是断电恢复测试。']
(out/'report.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:summary[k] for k in ['vina_jobs','vina_nonzero_or_missing_exit','output_qc_rejections','screen_submitted_per_hour','refine_pose_consistent_count','total_preparation_and_monitored_run_seconds']},indent=2))
