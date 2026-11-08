"""Produce a source audit and an explicitly UNFROZEN round-2 protocol draft.
No activity values or docking jobs are generated without verified experimental data.
"""
from pathlib import Path
import json,hashlib,html,datetime
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors,Draw
root=Path(__file__).resolve().parents[1];out=root/'results/tmgh1-round2';src=out/'source';r1=root/'results/tmgh1-round1'
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
raw=json.loads((out/'compound-map.json').read_text());assert len(raw)==18
records=[];mols=[];legends=[]
expected=['NOJ','IFM','NOY','OXZ','AZF','CTS','CGB','IFL','GOX','GTL','NTZ','GIM','PGI','IDE','GI1','GI2','GI3','GI4']
for x,code in zip(raw,expected):
 assert x['ccd']==code
 ds={d['program']:Chem.MolFromSmiles(d['descriptor']) for d in x['descriptors'] if d['type']=='SMILES_CANONICAL'}
 assert all(m is not None for m in ds.values())
 m=ds['CACTVS'];smi=Chem.MolToSmiles(m);formula=rdMolDescriptors.CalcMolFormula(m);issues=[]
 same=len({Chem.MolToSmiles(z) for z in ds.values()})==1
 if not same:
  assert code=='GTL';issues.append('OpenEye omits C=N stereo. CACTVS includes Z, matching CCD bond C1-N1 Z. CCD name also contains (2E), so naming/descriptor discrepancy requires original-paper verification. Retain CACTVS as provisional, not chemically approved.')
 if code=='GI1':
  assert not m.HasSubstructMatch(Chem.MolFromSmarts('[CX3]=[OX1]'))
  issues.append('CCD title says methoxycarbonyl, but CCD graph has methoxy(hydroxy)methyl and no C=O. The graph has an additional tetrahedral stereocenter. Do not repair from the name or admit to activity comparison until original chemical drawing is checked.')
 rec={'compound':x['compound'],'pdb':x['pdb'],'ccd':code,'name':x['name'],'mapping_source':'10.1021/ja066961g.s001, Table B, S2-S5','ccd_smiles_provisional':smi,'ccd_formula':x['chem_comp']['formula'],'rdkit_formula':formula,'formal_charge_in_CCD':Chem.GetFormalCharge(m),'heavy_atoms':m.GetNumHeavyAtoms(),'rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(m),'descriptors_agree':same,'stereocenters':Chem.FindMolChiralCenters(m,includeUnassigned=True),'chemical_review_issues':issues,'chemical_review_status':'hold_original_drawing' if issues else 'CCD_graph_consistent_assay_state_pending','round1_calibration':code in ['NOJ','GI2'],'experiment':{'endpoint':None,'value':None,'unit':None,'relation':None,'uncertainty':None,'pH':None,'temperature_C':None,'buffer':None,'substrate':None,'method':None,'source_page_table':None},'screening_ready':False}
 records.append(rec);mols.append(m);legends.append(f"{x['compound']:02d} {code} / {x['pdb']}"+(' [REVIEW]' if issues else ''))
 assert not any(v is not None for v in rec['experiment'].values())
dump(out/'library-audit.json',{'status':'source_inventory_only_NOT_docking_library','records':records})
# Derived scientific 2D structure inventory. Preserve deposited states; these are not pH assignments.
Draw.MolsToGridImage(mols,molsPerRow=3,subImgSize=(370,245),legends=legends).save(out/'compound-inventory.png')
inputs=[]
for g in ['experimental_2J7E','predicted_42','experimental_2J77']:
 for state in ['default','E166H']:
  p=r1/'receptors'/f'{g}_{state}'/'receptor.pdbqt';assert p.is_file();inputs.append({'geometry':g,'state':state,'path':str(p),'sha256':sha(p)})
center=json.loads((r1/'geometries/experimental_2J7E/metadata.json').read_text())['box_center']
protocol={'status':'DRAFT_NOT_FROZEN_NOT_EXECUTED','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'target':'TmGH1 Q08638 native 446 aa','design':'Retrospective protocol-frozen evaluation; never call prospective or blinded. Freeze only after data and chemistry are verified and before docking rankings.','blocking_requirements':['Original exact Ki or Kd numerical table with units, censoring and experimental method, pH, temperature, buffer and substrate as applicable','Original compound drawings to resolve GI1 and GTL metadata discrepancies','Choose and document protonation/tautomer rules based on verified assay conditions; never optimize chemical state against experimental ranks'], 'primary_endpoint_priority':'Prefer same-condition ITC Kd if available with adequate range and precision; otherwise same-condition kinetic Ki. Never pool endpoints, conditions, per-compound optimum pH, or unverified digitized values. Endpoint choice must precede computed scores.','cohort':'All verified compounds from the same study and assay stratum; chemical exclusions before scores. Exclude round-1 NOJ and GI2 from primary novel-compound analysis; report all eligible compounds descriptively. A small single-study cohort is not independent prospective validation.','receptors':{'primary':['experimental_2J7E','predicted_42'],'sensitivity':['experimental_2J77'],'round1_selection_bias':'2J7E selected after its own GI2 pose succeeded; record this bias.','protein_states':'Default GLU166/GLU351 baseline; E166H sensitivity only after assay review. State choice not yet finalized.'},'vina':{'version':'1.2.7','exhaustiveness':[8,32],'primary_exhaustiveness':32,'seeds':[42,2026,20261021],'workers':3,'cpu_per_job':2,'num_modes':9,'energy_range':4,'one_box_for_all_compounds':{'center':center,'size':[24,24,24],'reference':'Round1 aligned GI2 center in 2J77 frame; known-site screen, not de novo pocket discovery'},'job_count_formula':'N_ligand_states * N_receptor_geometries * N_receptor_protonation_states * 2 exhaustiveness values * 3 seeds; 18*3*1*2*3=324 only if all18 and one ligand/receptor state are approved. This is not a completed or launched job count.'},'analysis':{'score':'Median top-1 Vina score over all three seeds; lower ranks first. No choosing best-performing seed.','potency':'-log10(K_in_molar); larger is stronger; Spearman(-median_Vina, pK).','failed_job_policy':'Do not silently remove failures or compute median with fewer than3 seeds; flag incomplete compounds, report coverage and apply common complete cohort for comparisons.','spearman':'Report n, rho, two-sided permutation p (100000 fixed-seed shuffles) and paired-compound bootstrap95% interval (10000 draws); intervals can be wide. No multiple-condition cherry-picking.','top_k':'K=3 primary; K=5 descriptive when cohort size permits. Fraction of experimentally strongest K recovered among computational top K; expected random recall=K/n. At ties on either boundary, use equal fractional inclusion weights. Report random distribution, not merely expectation. This is not a clinical active/inactive classifier.','baselines':['Random permutation distribution','Heavy atom count descending (predeclared simple size baseline; tie averaging)','Molecular weight descending (descriptive correlated size baseline)'],'stability':'Seed-to-seed Spearman, top-K overlap, score dispersion; compare receptors and E8/E32 on identical compounds. Repeatability does not establish validity.','cost':'Monitor per-job and batch wall time, peak process tree RSS, system available RAM and swap before/after. Paired E32/E8 time ratio and ranking changes. Separate reused OpenFold inference costs from new round2 costs.','decision':'Ranking evidence remains separate from round1 failed pose-transfer gate. Positive exploratory association alone does not authorize automatic candidate selection; weak/uncertain results explicitly fail the evidence requirement.'},'reused_receptor_inputs':inputs,'source_hashes':{str(src/'ja066961g_si_001.pdf'):sha(src/'ja066961g_si_001.pdf'),str(out/'library-audit.json'):sha(out/'library-audit.json')}}
dump(out/'protocol-draft.json',protocol)
status={'stage':'data_readiness_blocked','round2_docking_started':False,'round2_jobs_completed':0,'compound_pdb_ccd_mappings':18,'exact_same_condition_experimental_records':0,'chemical_metadata_holds':['15/GI1','10/GTL'],'all_ligand_assay_states_pending':True,'requires_user_material':'Readable main article PDF or original numerical tables + methods + compound figure for DOI10.1021/ja066961g','scientific_ranking_conclusion':None,'automatic_candidate_prioritization_supported':False,'round1_unchanged':True}
dump(out/'decision.json',status)
text='''TmGH1 第二轮：数据准入审计（尚未开始对接）

本轮目标：检验 OpenFold3 + Vina 能否依据真实实验常数优先排序候选物。
状态：已完成公开数据搜集与化合物映射，卡在实验数值/条件核验。第二轮对接任务0个；没有新的 Spearman、Top-K、运行时间或内存峰值。不能把本文件称为“第二轮筛选验证完成报告”。

已完成
1. 下载并校验 ACS 官方公开支持信息（13页，2,840,961字节，原始MD5与Figshare元数据一致）。阅读全部提取文本，目视核对S1-S5表格。
2. 根据 Table B 建立18个论文编号 → PDB → CCD 映射，下载各条目的RCSB结构元数据和CCD化学字典。逐条排除醋酸根/钙离子等结晶成分；未把它们作为候选抑制剂。
3. 规范化并交叉核对两个CCD SMILES生成器，记录立体化学、分子式、形式电荷和异常；保存原始数据及SHA256。
4. 为后续验证准备可复用受体文件清单和哈希；起草E8/E32、三随机种子、统一口袋框、连续实验常数评价及失败处理规则。协议尚未冻结，未启动Vina或重新预测OpenFold3。

为什么实验数据不合格
公开附录Table A的 Kinetics [I]、ITC [I]、ITC [E] 是实验投料/测试浓度，不是Ki或Kd，不能直接拿来做强弱标签。Table B是晶体学数据。Table C是ITC曲线和pH依赖的1/Ki图，未提供可直接转录的完整精确常数表。特别是化合物2图的1/Ki单位为nM^-1，其他多为µM^-1，不能忽略1000倍单位差异。
本轮未从曲线目测、拟合或猜测实验数值。pH曲线不能替代统一条件的精确常数表；各化合物最优pH的数据也不能混为同一实验条件。晶体pH7不能被当成抑制/ITC实验pH。温度、缓冲体系、底物及误差需从正文确认。

获取记录
官方Figshare支持信息获取成功；PubMed、EuropePMC、OpenAlex、Semantic Scholar的已核查记录未提供可读开放正文。ACS正文/PDF请求返回403；York和UWA页面提供书目信息，未取得正文。ChEMBL按准确DOI查询0条，不能误把无效过滤器返回的其他文献当本研究数据。旧Binding MOAD域名当前返回其他内容，未用其作为数据来源。尝试的其他公开链接没有取得经核验的实验表。以上仅说明本次未取得，不能推断论文数据不存在或完全无法获得。

化学异常
15/GI1：名称为 methoxycarbonyl-substituted glucoimidazole，但当前CCD键图没有C=O，而是methoxy(hydroxy)methyl，且多一个四面体手性中心。两个SMILES生成器一致也不能证明对应实验化合物正确。保留原条目、待核对论文结构图；不能凭名称自动“修复”。
10/GTL：OpenEye SMILES缺少C=N几何立体信息；CACTVS描述及CCD键表标Z，但CCD名称包含(2E)。暂存CACTVS有立体信息版本，仍需原图核验，尚未批准用于对接。
全部18种：CCD形式电荷是字典表示，不能自动视为指定实验pH下的主要微观状态。受体Glu166/Glu351、配体质子化和互变异构体规则均须在取得条件后固定，不得按实验排序表现倒选状态。

第一轮证据继续有效
第一轮结构预测全局/口袋吻合较好，但主预测结构的主条件姿态复现0/6，实验自身重对接3/6。第二轮必须分别报告姿态正确性、排序相关性和资源表现，不能相互代替。继续使用Q08638的446残基天然序列预测；468是晶体表达构建体长度。

下一步材料（不需要重新批准运行）
需要 DOI10.1021/ja066961g 的可读正文PDF，或化合物结构图 + 实验Ki/Kd数值表及表注 + 测量方法页。若只拿到部分同条件精确值，就固定符合标准的子集，明确排除原因；不随机换蛋白或掺入别的终点。
取得后：核验单位/条件/化学身份 → 固定实验队列和状态 → 冻结协议及哈希 → E8/E32三种子对接 → 独立审计输出身份/完整性 → 相关性、Top-K、随机/尺寸基线、稳定性和资源报告。候选库规模会由数据准入决定。

电脑与资源
本次仅网络下载、PDF读取和RDKit结构整理，没有进行新的OpenFold3预测或Vina批量计算。没有采集新的计算峰值，不引用第一轮内存作为本轮实测。新增PyMuPDF1.26.7仅安装在项目.tools/tmgh1-round2-deps，主对接与预测环境未修改；PDF目视核对最终使用已有的pypdfium2。
'''
(out/'report.txt').write_text(text)
rows=''.join('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in [r['compound'],r['pdb'],r['ccd'],r['formal_charge_in_CCD'],r['heavy_atoms'],'待核对原图' if r['chemical_review_issues'] else 'CCD一致；实验状态待定'])+'</tr>' for r in records)
body=''.join('<p>'+html.escape(s).replace('\n','<br>')+'</p>' for s in text.split('\n\n'))
(out/'report.html').write_text('''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>TmGH1 第二轮数据审计</title><style>body{max-width:1060px;margin:42px auto;padding:0 26px;font-family:-apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;line-height:1.85;color:#203347;background:#fafcff}h1{font-size:30px}p{white-space:normal}.status{padding:18px 24px;background:#fff1d9;border-left:5px solid #c47a09}table{width:100%;border-collapse:collapse;font-size:14px}th,td{padding:8px;border-bottom:1px solid #dde3eb;text-align:left}th{background:#e9eff5}img{width:100%;background:white}a{color:#1558a6}footer{margin-top:30px;color:#607084;font-size:13px}</style><h1>TmGH1 第二轮 · 数据准入审计</h1><div class="status"><b>实验验证尚未开始</b> · 18个化合物映射 · 0个已核实同条件精确实验记录 · 0个第二轮对接任务</div>'''+body+'<h2>化合物清单</h2><table><thead><tr><th>论文编号</th><th>PDB</th><th>CCD</th><th>CCD形式电荷</th><th>重原子数</th><th>状态</th></tr></thead><tbody>'+rows+'</tbody></table><p>下图保留CCD表示，不代表实验pH下主要状态；标记REVIEW的条目待核验。</p><img src="compound-inventory.png" alt="18种化合物的CCD结构清单"><footer>原文 DOI：<a href="https://doi.org/10.1021/ja066961g">10.1021/ja066961g</a> · <a href="https://acs.figshare.com/articles/journal_contribution/3021493">官方公开支持信息</a> · 本地：<a href="protocol-draft.json">协议草稿</a> · <a href="library-audit.json">完整数据审计</a> · <a href="decision.json">机器可读状态</a></footer></html>')
meta=json.loads((src/'si-metadata.json').read_text());pdf=src/'ja066961g_si_001.pdf';assert hashlib.md5(pdf.read_bytes()).hexdigest()==meta['files'][0]['computed_md5']
manifest={str(p.relative_to(root)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='verification.json'}
manifest.update({str(p.relative_to(root)):sha(p) for p in (root/'tmgh1_round2').glob('*.py')})
dump(out/'verification.json',{'files_sha256':manifest,'file_count':len(manifest),'SI_MD5_matches_official':True,'compound_mapping_count':18,'experimental_records_complete':0,'docking_jobs':0,'round1_receptor_hashes':inputs})
print(json.dumps(status,ensure_ascii=False,indent=2));print('Saved report and',len(manifest),'file hashes')
