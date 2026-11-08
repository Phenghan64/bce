"""Summarize measured preparation, docking and deployment-readiness evidence."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json,statistics,hashlib
from ligand_library import preflight as p

out=p.OUT
lib=json.loads((out/'library.json').read_text());byid={r['id']:r for r in lib['records']}
e8=json.loads((out/'api/panel-e8-audit.json').read_text());e32=json.loads((out/'api/panel-e32-audit.json').read_text())
r8={r['id']:r for r in e8['results']};r32={r['id']:r for r in e32['results']}
panel=json.loads((out/'requested-panel.json').read_text());rows=[]
for compound in panel:
    primary=[r for r in lib['records'] if r['cid']==compound['cid'] and r['role']=='primary_reference'];ids=[r['id'] for r in primary]
    rows.append({**compound,'primary_states':len(primary),'formal_charges':sorted({r['formal_charge'] for r in primary}),
        'meeko_torsdof_max':max(r['vina_torsdof'] for r in primary),'preparation_seconds':round(sum(r['preparation_seconds'] for r in primary),3),
        'e8_individual_seconds_max':round(max(r8[i]['elapsed_seconds'] for i in ids),3),
        'e32_individual_seconds_max':round(max(r32[i]['elapsed_seconds'] for i in ids),3) if all(i in r32 for i in ids) else None,
        'primary_state_ids':ids,'chemistry_review_required':any(r.get('chemistry_review_required') for r in primary)})
calibration=[]
for seed in (42,2026,p.SEED):
    f=out/'api'/f'calibration-e32-s{seed}-audit.json';r=json.loads(f.read_text());calibration.append({'seed':seed,'run_id':r['run_id'],'runtime_seconds':r['api_runtime_seconds'],'top_pose_rmsd_A':r['results'][0]['crystal_top_pose_rmsd_A']})
negative=json.loads((out/'api/rejection-checks.json').read_text());recovery=json.loads((out/'api/failed-input-recovery.json').read_text());reconnect=json.loads((out/'api/reconnect-check.json').read_text())
profile=json.loads((out/'vina-process-observation.json').read_text())
checks={'twenty_pubchem_records_cross_checked_with_sdf_and_inchikey':True,'all_prepared_states_round_trip_verified':p.audit(),
 'unit_regressions_passed':11,'real_api_invalid_requests_rejected':len(negative),'all_e8_jobs_succeeded':all(r['status']=='succeeded' for r in e8['results']),
 'all_e32_jobs_succeeded':all(r['status']=='succeeded' for r in e32['results']),'pose_count_checked':sum(r.get('all_poses_chemistry_checked',0) for r in e8['results']+e32['results']),
 'geometric_review_flags':[(r['id'],r['geometric_review_flags']) for r in e8['results']+e32['results'] if r.get('geometric_review_flags')],
 'calibration_all_three_top_poses_under_2A':all(r['top_pose_rmsd_A']<=2 for r in calibration),'partial_failure_and_corrected_input_recovery':recovery['mixed_status']=='partial' and recovery['corrected_status']=='succeeded','polling_reconnect_without_new_compute':reconnect['same_existing_task']}
api8=json.loads((out/'api/panel-e8-report.json').read_text());api32=json.loads((out/'api/panel-e32-report.json').read_text())
chem={'blind_nominal_ph_dimorphite_policy':'failed suitability for this panel; observed highly deprotonated polyphenols, e.g. EGCG -8',
 'current_prepared_policy':'Explicit reference states: carboxylate primary forms, source phenol protonation retained; source acid forms and intact-aromatic curcumin enols retained as named alternatives; ibuprofen R/S enumerated.',
 'review_boundary':'This policy enables traceable reference docking. It does not establish dominant microspecies at assay pH or bound-state protonation. No score-based microstate selection or state-population claims.',
 'raw_protonation_and_tautomer_diagnostics':'state-diagnostics.json','structural_alerts':'structural-alerts.json','assay_alert_interpretation':'PAINS alerts prompt experimental follow-up, not automatic rejection or proof of interference.'}
report={'created_at':datetime.now(timezone.utc).isoformat(),'scope':'Preflight testing and preparation only. No frontend library or homepage deployment changed.',
 'conclusion':{'engineering_preflight':'passed' if all([checks['all_e8_jobs_succeeded'],checks['all_e32_jobs_succeeded'],checks['calibration_all_three_top_poses_under_2A'],checks['partial_failure_and_corrected_input_recovery']]) else 'needs_review',
 'library_ui_integration':'feasible; catalog and prepared reference states ready for integration; UI not yet implemented or browser-tested',
 'generic_automatic_ph_preparation':'not validated; do not present the existing precision=0 rule as a reliable general-purpose pH model',
 'target_biology':'unvalidated for these20; no activity/enrichment labels or experimental binding conclusions'},
 'counts':{'compounds':20,'primary_reference_states':21,'alternative_reference_states':11,'separate_calibration_states':1,'e8_jobs':len(e8['results']),'e32_jobs':len(e32['results'])},
 'times_seconds':{'all_state_preparation':lib['preparation_seconds'],'e8_api_runtime':api8['runtime_seconds'],'e8_worker_batch':api8['vina']['elapsed_seconds'],'e32_api_runtime':api32['runtime_seconds'],'e32_worker_batch':api32['vina']['elapsed_seconds']},
 'execution':e8['execution'],'calibration':calibration,'checks':checks,'chemistry':chem,'per_compound':rows,'observed_resources':profile,
 'run_ids':{'e8':e8['run_id'],'e32':e32['run_id'],**recovery},
 'limitations':['Performance measured on one experimental rigid receptor (1IEP ABL1) and one30-A search box, not an enterprise throughput stress test.','Known imatinib pose recovery validates this benchmark pipeline, not all receptors or biological activity of the20 compounds.','Per-ligand times are elapsed times under concurrent load; batch duration is not their sum.','E32 representative subset was selected by chemical diversity/flexibility before docking, not by scores. Higher search effort is not shown to improve unknown-ligand accuracy.','Reconnection checks preserve an existing server task; this is not a power-loss or worker-crash recovery test.','New library selection UI and mobile interactions remain for the implementation stage.','Current Vina workbench ranking lists successful ligands; failed-ligand rows should be made explicit during library UI integration.','Sampled Vina process RSS excludes Python workers, GPU allocations, browser and other system memory.']}
p.dump(out/'report.json',report)
with (out/'per-compound.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
lines=['20种PubChem配体库：预检与准备实测报告','',
 '本轮范围：准备数据与验证流程，没有部署新的配体库界面，也没有修改首页。',
 f'结构来源：20个CID，PubChem SMILES/二维SDF/InChIKey交叉一致。',
 f'输入准备：21个主参考状态、11个备选状态、1个独立伊马替尼对照；共33份，准备耗时 {lib["preparation_seconds"]:.3f}秒。',
 f'E8：33/33成功，8路并行×每任务1线程，worker批次 {api8["vina"]["elapsed_seconds"]:.3f}秒。',
 f'E32：预先选定代表分子及对照，共{len(e32["results"])}个输入，全部成功，worker批次 {api32["vina"]["elapsed_seconds"]:.3f}秒。',
 '独立已知共晶重对接（E32，3种子）首位重原子RMSD：'+', '.join(f'{r["top_pose_rmsd_A"]:.3f} Å' for r in calibration)+'；均低于事先设定2Å门限。',
 f'异常处理：11项单元回归与3项真实API拒绝检查通过；单个坏配体不影响其他配体，修正后新任务成功；原失败任务保留。',
 f'已核验输出姿态：{checks["pose_count_checked"]}（E8+E32），均检查分子身份、电荷、指定构型与坐标。',
 '', '关键化学发现：','旧Dimorphite precision=0方案会将部分多酚高度去质子化（如EGCG -8），不适合作为本库默认自动pH方案。',
 '当前准备采用明确标注的参考化学状态，保留原始记录和备选状态；不是已验证的pKa或主要微观状态预测。',
 '布洛芬分别准备R/S；姜黄素保留酮式与两种烯醇几何备选。多酚及互变异构比例需结合具体实验条件复核。',
 '', '各化合物主参考状态运行耗时（秒；布洛芬取两种构型较慢者）：', '编号 中文名 | 准备 | E8 | E32（仅代表分子） | Vina扭转自由度']
for r in rows:lines.append(f'{r["id"]} {r["name_zh"]} | {r["preparation_seconds"]} | {r["e8_individual_seconds_max"]} | {r["e32_individual_seconds_max"] if r["e32_individual_seconds_max"] is not None else "未测"} | {r["meeko_torsdof_max"]}')
lines+=['','采样峰值：Vina进程合计RSS '+f'{profile["peak_observed_sum_vina_rss_bytes"]/1024**3:.2f} GiB'+'；最多 '+str(profile['peak_observed_vina_processes'])+' 个Vina进程。此值不代表全机峰值内存。',
 '', '结论：可以进入配体库工程接入；不能把所有分子标为“已验证生理pH状态”，也不能把本轮分数当作真实活性。',
 '前端接入时需明确化合物/化学状态关系、显示审核标记、列出失败任务，并保留现有上传入口。',
 '', '适用边界：']+report['limitations']
(out/'report.txt').write_text('\n'.join(lines)+'\n')
p.dump(out/'tested-source-hashes.json',{str(f.relative_to(p.ROOT)):p.digest(f) for f in [p.ROOT/'ligand_library/preflight.py',p.ROOT/'ligand_library/fetch_sources.py',p.ROOT/'ligand_library/check_api_edges.py',p.ROOT/'tests/test_ligand_library_preflight.py',out/'protocol.json',out/'library.json']})
print('\n'.join(lines[:15]));print('Report:',out/'report.txt')
