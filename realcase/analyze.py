"""Evaluate the frozen ABL1 panel against experimental labels; preserve failures."""
import json,itertools
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

out=Path(__file__).resolve().parents[1]/'results/abl-real-case'
read=lambda name:json.loads((out/name).read_text())
protocol=read('protocol.json');library=read('library.json');run=read('panel-validated.json');gate=read('calibration.json')
assert gate['passed'] and run['completed']+run['failed']==48
names={str(x['CID']):x['Title'] for x in read('source/pubchem-titles.json')['PropertyTable']['Properties']}
rows=[]
for c in protocol['panel']:
    entry=next(x for x in library['records'] if x.get('cid')==c['cid'])
    jobs=[r for r in run['results'] if r['spec']['ligand_id']==entry['id']]
    row={**c,'name':names[c['cid']],'formal_charge':entry['formal_charge'],'successful_runs':sum(x['status']=='completed' for x in jobs)}
    if row['successful_runs']==3:
        scores=[x['top_score'] for x in jobs]
        poses=[Chem.SDMolSupplier(x['sdf'],removeHs=True)[0] for x in jobs]
        row.update(score=float(np.median(scores)),scores=scores,score_range=float(np.ptp(scores)),
            max_top_pose_pairwise_rmsd_A=max(float(rdMolAlign.CalcRMS(a,b)) for a,b in itertools.combinations(poses,2)))
    else:row['errors']=[x.get('error') for x in jobs if x['status']!='completed']
    rows.append(row)
metrics=None
if all(x['successful_runs']==3 for x in rows):
    rows.sort(key=lambda x:(x['score'],int(x['cid'])))
    for i,x in enumerate(rows,1):x['rank']=i
    positive=np.array([-r['score'] for r in rows if r['label']=='strong'])
    negative=np.array([-r['score'] for r in rows if r['label']!='strong'])
    def auc(p,n):
        d=p[:,None]-n[None,:]
        return float(np.mean((d>0)+.5*(d==0)))
    observed=auc(positive,negative);rng=np.random.default_rng(20261021)
    bootstrap=[auc(rng.choice(positive,8),rng.choice(negative,8)) for _ in range(5000)]
    combined=np.r_[positive,negative];null=[]
    for indices in itertools.combinations(range(16),8):
        mask=np.zeros(16,dtype=bool);mask[list(indices)]=True
        null.append(auc(combined[mask],combined[~mask]))
    top4=sum(r['label']=='strong' for r in rows[:4])
    metrics={'roc_auc':observed,'roc_auc_stratified_bootstrap95':np.percentile(bootstrap,[2.5,97.5]).tolist(),
      'one_sided_exact_label_permutation_p':float(np.mean(np.array(null)>=observed-1e-12)),
      'top4_strong_count':top4,'top4_precision':top4/4,'top4_enrichment_over_half_baseline':top4/2,
      'note':'Small selected extreme-group retrospective panel. Bootstrap interval excludes dataset/selection uncertainty; no population-level efficacy claim.'}
times={'receptor_repair':read('receptor/repair.json')['seconds'],'receptor_meeko':read('receptor/meeko-timing.json')['seconds'],
       'ligand_preparation':library['preparation_seconds'],'calibration_stage':read('calibrate-resources.json')['seconds'],
       'panel_stage':read('panel-resources.json')['seconds']}
report={'case':protocol['case'],'receptor':protocol['receptor'],'calibration':gate,'completed_vina_jobs':run['completed']+sum(len(a['checks']) for a in gate['attempts']),
    'failed_panel_jobs':run['failed'],'output_audit':{k:run[k] for k in ['audit_rule','audit_origin','reclassified_jobs']},'rows':rows,'metrics':metrics,'stage_seconds':times,'sum_final_stage_seconds':sum(times.values()),
    'timing_scope':'Final retained preparation and measured docking stages only; excludes downloads, code development, failed receptor-preparation attempts and discarded ligand preparation.',
    'openfold3_in_this_run':False,'limitations':protocol['limitations'],'resources':{s:read(s+'-resources.json') for s in ['calibrate','panel']}}
(out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
metrics_text=(f"AUC={metrics['roc_auc']:.3f}；分层bootstrap95%区间={metrics['roc_auc_stratified_bootstrap95']}；前4名有{metrics['top4_strong_count']}个强结合化合物，富集倍数={metrics['top4_enrichment_over_half_baseline']:.2f}。" if metrics else '存在未通过复核的任务，不发布完整面板AUC。')
lines=['真实案例：人源 ABL1 / Davis 2011 实验标签回顾性验证','',
 '结论：重对接校准通过，工程流程可运行。当前小样本排序存在一定信号，但不足以证明可靠的生产筛选效果。',
 metrics_text,
 '本轮AUC的bootstrap区间覆盖0.5；单侧精确标签置换p值约0.052，仅作补充描述。',
 '关键误排：Quizartinib在本组实验属弱/未检出，却排第2；Bosutinib实验Kd为0.12nM，却排第14。',
 '这些结果只评价当前受体构象、配体状态和参数组合，不能概括所有Vina使用场景。','',
 '受体：2HYY B链，人源无突变、伊马替尼结合的单一晶体构象。',
 '8个实验Kd≤100nM的化合物，8个数据集中10000nM封顶/未检出结合的化合物。',
 '10000nM不解释为精确Kd；伊马替尼另作重对接校准，不进入分类指标。',
 '统一pH7.4启发式单质子化状态，原始互变异构状态；3种随机种子，固定受体，30Å搜索盒。',
 '受体缺失ASP276被补建，补建原子及氢经Amber14局部优化，实验重原子固定。',
 '本轮未重新运行OpenFold3；没有将模型置信度当作亲和力或活性证据。','',
 '重对接校准：'+json.dumps(gate,ensure_ascii=False),
 '评价：'+json.dumps(metrics,ensure_ascii=False),'',
 '排名 | 化合物 | CID | 实验Kd(nM)或标签 | Vina三次中位分数(kcal/mol) | 三次极差 | 姿势最大两两RMSD(Å)']
for r in rows:
    label=str(r['kd_nm']) if r['label']=='strong' else '10000封顶/未检出'
    lines.append(f"{r.get('rank','失败')} | {r['name']} | {r['cid']} | {label} | {r.get('score')} | {r.get('score_range')} | {r.get('max_top_pose_pairwise_rmsd_A')}")
lines+=['','阶段耗时(秒)：'+json.dumps(times,ensure_ascii=False),
    f"最终阶段合计：{sum(times.values()):.2f}秒（{sum(times.values())/60:.2f}分钟）。",
    '时间口径：最终保留的准备步骤和对接计算；不含下载、开发调试和早期失败准备，不能解释为本次对话总耗时。','',
    f"正式对接阶段进程树RSS采样峰值约{report['resources']['panel']['sampled_process_tree_rss_peak_bytes']/2**30:.2f}GiB；全系统可用内存最低约{report['resources']['panel']['minimum_system_available_bytes']/2**30:.2f}GiB。RSS相加可能重复计入共享页面。",
    f"正式对接期间系统swap占用变化：{report['resources']['panel']['swap_after']['used']-report['resources']['panel']['swap_before']['used']}字节；系统此前已有swap占用。",
    '重复稳定性：CI-1040、GDC-0879和Selumetinib最高分姿势在不同种子间偏差较大，详见上表。重复稳定本身也不证明结合真实。','',
    '局限：样本只有16个且按强弱极端分组抽取；单一失活构象可能不利于I型抑制剂；单质子化状态与互变异构状态未经逐分子验证；没有实验复测或前瞻性验证。',
    'Vina分数用于此流程内排序，不能换算为可靠Kd，也不等同于真实药效。',
    '输出复核：'+json.dumps(report['output_audit'],ensure_ascii=False),'',
    '数据原始研究：https://doi.org/10.1038/nbt.1990',
    '公开数据：https://staff.cs.utu.fi/~aatapa/data/DrugTarget/',
    '结构：https://www.rcsb.org/structure/2HYY',
    'protocol.json保存对接前冻结的规则及修改；provenance.json保存来源与哈希；jobs/保存全部日志和构象。']
(out/'report.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'metrics':metrics,'times':times,'ranks':[(r.get('rank'),r['name'],r.get('score'),r['label']) for r in rows]},indent=2))
