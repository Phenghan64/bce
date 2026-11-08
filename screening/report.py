"""Summarize the unlabelled ABL throughput pilot without inventing activity labels."""
import hashlib
import itertools
import json
import statistics
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign

ROOT = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1]).resolve()
library = json.loads((out/'library.json').read_text())
entries = {r['id']:r for r in library['records']}
stages = {n:json.loads((out/f'{n}.json').read_text()) for n in ['tune-1x6','tune-3x2','screen','refine']}
screen = stages['screen']; refined = stages['refine']
control_id = 'imatinib-control'
successful = sorted((r for r in screen['results'] if r['status']=='completed'),
                    key=lambda r:(r['top_score'],r['spec']['ligand_id']))
ranked = []
for rank, result in enumerate(successful,1):
    entry = entries[result['spec']['ligand_id']]
    ranked.append({'rank_including_control':rank,'id':entry['id'],'role':entry['role'],
        'score_kcal_mol':result['top_score'],'heavy_atoms':entry['heavy_atoms'],
        'rotatable_bonds':entry['rdkit_rotatable_bonds'],'molecular_weight':entry['molecular_weight'],
        'formal_charge':entry['formal_charge'],'input_smiles':entry['smiles'],'sdf':result['sdf']})
(out/'screen-ranking.json').write_text(json.dumps(ranked,indent=2)+'\n')
with Chem.SDWriter(str(out/'screen-top-poses.sdf')) as writer:
    for row in ranked:
        mol = Chem.SDMolSupplier(row['sdf'])[0]
        mol.SetProp('_Name',row['id'])
        mol.SetProp('Vina_score_kcal_mol',str(row['score_kcal_mol']))
        mol.SetProp('Rank_including_control',str(row['rank_including_control']))
        mol.SetProp('Role',row['role'])
        writer.write(mol)
refine_ids = json.loads((out/'refinement-selection.json').read_text())['ids']
refine_rows = []
for ligand_id in refine_ids:
    runs = sorted((r for r in refined['results'] if r['spec']['ligand_id']==ligand_id and r['status']=='completed'),
                  key=lambda r:r['spec']['seed'])
    molecules = [Chem.SDMolSupplier(r['sdf'])[0] for r in runs]
    pairwise = [{'seeds':[runs[i]['spec']['seed'],runs[j]['spec']['seed']],
                 'rmsd_angstrom':float(rdMolAlign.CalcRMS(molecules[i],molecules[j]))}
                for i,j in itertools.combinations(range(len(runs)),2)]
    scores = [r['top_score'] for r in runs]
    row = {'id':ligand_id,'role':entries[ligand_id]['role'],
        'successful_seeds':len(runs), 'median_score_kcal_mol':statistics.median(scores) if scores else None,
        'score_range_kcal_mol':max(scores)-min(scores) if scores else None,
        'top_pose_pairwise_rmsd':pairwise,
        'max_top_pose_rmsd_angstrom':max((p['rmsd_angstrom'] for p in pairwise),default=None),
        'pose_consistent_within_2A':len(runs)==3 and all(p['rmsd_angstrom']<=2 for p in pairwise),
        'runs':[{'seed':r['spec']['seed'],'score_kcal_mol':r['top_score'],'sdf':r['sdf']} for r in runs]}
    refine_rows.append(row)
refine_rows.sort(key=lambda r:(r['median_score_kcal_mol'] is None,r['median_score_kcal_mol'] or 0,r['id']))
(out/'refine-ranking.json').write_text(json.dumps(refine_rows,indent=2)+'\n')
reference = Chem.SDMolSupplier(str(ROOT/'results/closed-loop-1iep/reference-imatinib.sdf'))[0]
control = {}
for stage_name in ['screen','refine']:
    control[stage_name] = []
    for run in stages[stage_name]['results']:
        if run['spec']['ligand_id']==control_id and run['status']=='completed':
            mol = Chem.SDMolSupplier(run['sdf'])[0]
            control[stage_name].append({'seed':run['spec']['seed'],'score_kcal_mol':run['top_score'],
                                       'crystal_rmsd_angstrom':float(rdMolAlign.CalcRMS(mol,reference))})
all_jobs = [r for stage in stages.values() for r in stage['results']]
for job in all_jobs:
    if job['status']=='completed':
        for path, expected in job['outputs_sha256'].items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected,path
rss = [r['peak_process_rss_bytes']/1024**2 for r in all_jobs if r['status']=='completed' and r['peak_process_rss_bytes']]
resume = json.loads((out/'screen-resume-check.json').read_text()) if (out/'screen-resume-check.json').exists() else None
summary = {'scope':'Unlabelled 100-compound NCI example library plus one known-pose imatinib control; ABL crystal receptor 1IEP',
    'library_prepared':sum(r['status']=='prepared' for r in entries.values()),
    'library_preparation_failed':sum(r['status']=='failed' for r in entries.values()),
    'preparation_wall_seconds':library['preparation_wall_seconds'],
    'stages':{name:{k:s[k] for k in ['wall_seconds','completed','failed','cache_hits']} for name,s in stages.items()},
    'scheduling':json.loads((out/'scheduling.json').read_text()),
    'screen_submitted_per_hour':len(screen['results'])/screen['wall_seconds']*3600,
    'screen_validated_per_hour':screen['completed']/screen['wall_seconds']*3600,
    'screen_effective_seconds_per_submitted_compound':screen['wall_seconds']/len(screen['results']),
    'summed_compute_stage_seconds':library['preparation_wall_seconds']+sum(s['wall_seconds'] for s in stages.values()),
    'timing_note':'Sum of measured computation stages; excludes development, source discovery, gaps and Python imports before timers. Resume check recorded separately.',
    'total_docking_jobs':len(all_jobs),
    'vina_process_failed_jobs':sum(r.get('returncode')!=0 for r in all_jobs),
    'validation_rejected_jobs':sum(r['status']!='completed' and r.get('returncode')==0 for r in all_jobs),
    'validated_pose_count':sum(r.get('pose_count',0) for r in all_jobs),
    'vina_peak_process_rss_mib':{'median':statistics.median(rss),'max':max(rss)},
    'memory_note':'Per-process maximum RSS from macOS time -l; not a measured simultaneous total or MPS usage.',
    'imatinib_control':control,
    'control_rank_in_validated':next(r['rank_including_control'] for r in ranked if r['id']==control_id),
    'ranked_compound_count':len(ranked),
    'refined_unlabelled_ligands_consistent_within_2A':sum(r['pose_consistent_within_2A'] for r in refine_rows if r['id']!=control_id),
    'refined_unlabelled_count':sum(r['id']!=control_id for r in refine_rows),
    'resume_check':{k:resume[k] for k in ['wall_seconds','completed','failed','cache_hits']} if resume else None,
    'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'screening').glob('*.py'))},
    'executed_code_snapshots_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((out/'source').glob('*-used.py'))},
    'failure_diagnosis':json.loads((out/'failure-diagnosis.json').read_text()),
    'fix_note':'Preparation now uses FindPotentialStereo to detect unresolved ring atom stereochemistry; existing pilot results retain the two rejected records. Resume now retains failures unless --retry-failed is given.',
    'limitations':[
        'No ABL activity labels for NCI compounds; enrichment/AUC/hit rate are not evaluated.',
        'Vina ranks are hypotheses, not binding affinities or evidence of activity.',
        'Source protonation and tautomer states retained; no pH/tautomer/stereoisomer ensemble.',
        'The drug-like filtered sample and fixed 1IEP receptor do not represent every library or target.',
        'Scheduling comparison used one small 9-ligand panel, not repeated or thermal-controlled benchmarking.',
        'Successful resume with completed outputs does not validate every interruption or failure scenario.',
        'Pose agreement across seeds is a consistency check, not proof of pose accuracy.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines = ['ABL / Vina 批量试筛结果', '',
    f'配体准备：{summary["library_prepared"]} 成功，{summary["library_preparation_failed"]} 失败。',
    f'对接任务：{len(all_jobs)} 个，程序失败 {summary["vina_process_failed_jobs"]} 个，结果校验拦截 {summary["validation_rejected_jobs"]} 个；通过校验的姿势 {summary["validated_pose_count"]} 个。',
    f'初筛 101 个分子：{screen["wall_seconds"]:.2f} 秒，处理 {summary["screen_submitted_per_hour"]:.1f} 个/小时，通过校验 {summary["screen_validated_per_hour"]:.1f} 个/小时。',
    f'前 10 名加对照，3 个种子复筛：{refined["wall_seconds"]:.2f} 秒。',
    f'所有计算阶段累计：{summary["summed_compute_stage_seconds"]:.2f} 秒。',
    f'复筛前 10 个未知活性候选中，{summary["refined_unlabelled_ligands_consistent_within_2A"]} 个的三次首位姿势两两 RMSD 均 ≤2 Å。',
    f'伊马替尼初筛排名（含对照）：{summary["control_rank_in_validated"]}/{len(ranked)}。',
    '', '伊马替尼对照：', json.dumps(control,ensure_ascii=False,indent=2),
    '', '初筛前 10 个库分子（分数不是活性结论）：']
for r in [r for r in ranked if r['id']!=control_id][:10]:
    lines.append(f'{r["id"]}: {r["score_kcal_mol"]:.3f} kcal/mol; {r["heavy_atoms"]} heavy atoms; {r["rotatable_bonds"]} rotatable bonds')
lines += ['', '适用边界：'] + summary['limitations']
(out/'report.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
