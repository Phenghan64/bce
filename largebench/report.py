"""Audit frozen inputs and all poses, then report paired performance measurements."""
import json,hashlib,re,sys,statistics,subprocess,platform
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from realcase.audit_outputs import compare_identity
out=root/'results/large-receptors'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
summary={'targets':{},'protocol':{'ligands_per_target':10,'seeds_per_ligand':3,'additional_control_jobs_per_target':3,'exhaustiveness':32,'workers':3,'cpu_per_job':2,'box_A':[30,30,30]},'limitations':['Single-chain experimental receptors, not full assemblies. mTOR is a truncated experimental construct, not full-length mTOR. Missing loops not modeled; all observed selected-chain residues retained.','This measures rigid-receptor local-pocket Vina docking, not OpenFold3 inference, flexible-receptor docking, or whole-surface blind docking.','Shared unlabelled panel, no target activity validation. Source ligand protonation/tautomer states retained. Pose controls do not establish hit rates.','RSS summation may double-count shared memory; system CPU, available RAM and swap include other applications. No direct temperature or throttling measurement.','Only 10 unique ligands and fixed run order (albumin then mTOR); differences include pocket chemistry and cannot be attributed solely to residue count.']}
all_hashes={};panel_hashes=[]
for name in ['albumin','mtor']:
 d=out/name;lib={r['id']:r for r in read(d/'library.json')['records']};panel={k:v['pdbqt_sha256'] for k,v in lib.items() if k!='crystal-control'};panel_hashes.append(panel)
 jobs=[];posecount=0;cpu=0;job_rows=[]
 for stage in ['calibrate','screen']:
  st=read(d/f'{stage}.json');assert st['failed']==0 and st['cache_hits']==0
  assert len(st['results'])==(3 if stage=='calibrate' else 30)
  for r in st['results']:
   spec=r['spec'];assert spec['receptor_sha256']==sha(spec['receptor']);assert spec['box_sha256']==sha(spec['box']);assert spec['protocol_sha256']==sha(d/'protocol.json');assert spec['library_sha256']==sha(d/'library.json');assert spec['ligand_sha256']==sha(lib[spec['ligand_id']]['pdbqt']);assert r['status']=='completed' and r['returncode']==0
   for p,h in r['outputs_sha256'].items():assert sha(p)==h;all_hashes[p]=h
   mols=list(Chem.SDMolSupplier(r['sdf'],removeHs=True));assert len(mols)==r['pose_count'];expected=Chem.MolFromSmiles(lib[spec['ligand_id']]['smiles'])
   for mol in mols:
    assert mol is not None and compare_identity(mol,expected)[0];assert np.isfinite(mol.GetConformer().GetPositions()).all()
   posecount+=len(mols);text=(Path(r['sdf']).parent/'vina.log').read_text();usage=re.search(r'([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys',text);assert usage
   job_cpu=float(usage[2])+float(usage[3]);cpu+=job_cpu
   job_rows.append({'stage':stage,'id':spec['ligand_id'],'seed':spec['seed'],'seconds':r['vina_wall_seconds'],'cpu_seconds':job_cpu,'top_score':r['top_score']});jobs.append(r)
 cal=read(d/'calibration.json');ref=Chem.SDMolSupplier(str(d/'receptor/crystal-control.sdf'),removeHs=True)[0]
 for r,c in zip(read(d/'calibrate.json')['results'],cal['checks']):
  actual=rdMolAlign.CalcRMS(Chem.SDMolSupplier(r['sdf'],removeHs=True)[0],ref);assert abs(actual-c['top_rmsd_A'])<1e-7
 screen=read(d/'screen.json');samples=[json.loads(line) for line in (d/'resources.jsonl').read_text().splitlines()];resources=read(d/'resources.json');assert resources['returncode']==0
 prep=read(d/'receptor/preparation.json');relax=read(d/'receptor/relaxation.json');conversion=read(d/'receptor/conversion.json');meta=read(d/'receptor/metadata.json');counts=read(d/'receptor/pdbqt-counts.json')
 assert prep['residues']==counts['residues']==meta['observed_residues'];assert conversion['returncode']==0
 grouped={i:[r['vina_wall_seconds'] for r in screen['results'] if r['spec']['ligand_id']==i] for i in panel}
 row={'name':name,'pdb':meta['pdb'],'author_chain':meta['author_chain'],'actual_residues':counts['residues'],'construct_sequence_residues':meta['construct_sequence_residues'],'pdbqt_atoms_including_polar_H':counts['atoms'],'prepared_heavy_atoms':prep['prepared_heavy_atoms'],'receptor_preparation_seconds':prep['seconds']+relax['seconds']+conversion['seconds'],'calibration_seconds':read(d/'calibrate.json')['wall_seconds'],'panel_seconds':screen['wall_seconds'],'monitored_total_seconds':resources['seconds'],'jobs_completed':len(jobs),'poses_checked':posecount,'panel_jobs_per_hour':30/screen['wall_seconds']*3600,'per_ligand_median_seconds':{k:statistics.median(v) for k,v in grouped.items()},'individual_job_median_seconds':statistics.median([r['vina_wall_seconds'] for r in screen['results']]),'individual_job_min_seconds':min(r['vina_wall_seconds'] for r in screen['results']),'individual_job_max_seconds':max(r['vina_wall_seconds'] for r in screen['results']),'process_tree_peak_RSS_GiB':resources['process_tree_peak_rss_bytes']/2**30,'system_min_available_GiB':resources['system_min_available_bytes']/2**30,'swap_before_GiB':resources['swap_before']['used']/2**30,'swap_after_GiB':resources['swap_after']['used']/2**30,'swap_peak_GiB':resources['swap_peak_used_bytes']/2**30,'system_swapout_delta_bytes':resources['swap_after']['sout']-resources['swap_before']['sout'],'vina_CPU_seconds':cpu,'vina_average_core_equivalents':cpu/resources['seconds'],'system_cpu_mean_percent':statistics.mean(s['system_cpu_percent'] for s in samples[1:]),'system_cpu_peak_percent':max(s['system_cpu_percent'] for s in samples),'control':cal,'gaps':meta['gaps'],'receptor_coordinate_audit':read(d/'receptor/coordinate-audit.json'),'pocket_repairs_within_8A':read(d/'receptor/pocket-repairs.json'),'jobs':job_rows}
 summary['targets'][name]=row
 for path in [d/'protocol.json',d/'library.json',d/'receptor/receptor.pdbqt',d/'receptor/receptor.box.txt',d/'resources.json',d/'resources.jsonl',d/'receptor/coordinate-audit.json',d/'receptor/pocket-repairs.json']:
  all_hashes[str(path)]=sha(path)
assert panel_hashes[0]==panel_hashes[1]
summary['verified_same_panel_pdbqt_hashes']=panel_hashes[0];summary['total_jobs']=sum(t['jobs_completed'] for t in summary['targets'].values());assert summary['total_jobs']==66
summary['total_poses_checked']=sum(t['poses_checked'] for t in summary['targets'].values());summary['paired_mtor_over_albumin_median_job_ratios']={k:summary['targets']['mtor']['per_ligand_median_seconds'][k]/v for k,v in summary['targets']['albumin']['per_ligand_median_seconds'].items()}
summary['hardware']={'cpu':subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),'physical_cores':int(subprocess.check_output(['sysctl','-n','hw.physicalcpu'],text=True)),'memory_bytes':int(subprocess.check_output(['sysctl','-n','hw.memsize'],text=True)),'platform':platform.platform(),'vina':subprocess.check_output([str(root/'.tools/bin/vina'),'--version'],text=True).strip()}
summary['thermal_status_end']=subprocess.run(['pmset','-g','therm'],capture_output=True,text=True).stdout
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for p in (root/'largebench').glob('*.py'):all_hashes[str(p)]=sha(p)
for p in (out/'source').glob('*'):all_hashes[str(p)]=sha(p)
(out/'verification.json').write_text(json.dumps({'all_passed':True,'jobs':66,'poses':summary['total_poses_checked'],'same_panel_inputs':True,'hashes':all_hashes},indent=2)+'\n')
lines=['两种大受体 × 同一组10个配体：本机 Vina 实测','', '机器：'+str(summary['hardware']),'协议：每配体3个种子；E32；3并发×2线程；30×30×30 Å 局部搜索框。每受体30次正式任务＋3次独立共晶重对接。','配体：按0至9可旋转键各选一个、优先较大分子，未参考此前对接分数；复用同一批准备文件。正式库配体24–38重原子。','']
for name,t in summary['targets'].items():
 lines += [f"{name} / {t['pdb']} author chain {t['author_chain']}",f"实际残基 {t['actual_residues']}；构建序列 {t['construct_sequence_residues']}；PDBQT原子（含极性H）{t['pdbqt_atoms_including_polar_H']}",f"受体准备 {t['receptor_preparation_seconds']:.2f}s；共晶检查 {t['calibration_seconds']:.2f}s；10配体×3种子 {t['panel_seconds']:.2f}s；监测全过程 {t['monitored_total_seconds']:.2f}s",f"单任务中位数 {t['individual_job_median_seconds']:.2f}s；范围 {t['individual_job_min_seconds']:.2f}–{t['individual_job_max_seconds']:.2f}s；批次速率 {t['panel_jobs_per_hour']:.1f} 次/小时（含重复种子，不是不同化合物数）",f"任务RSS峰值 {t['process_tree_peak_RSS_GiB']:.3f}GiB；系统可用内存最低 {t['system_min_available_GiB']:.3f}GiB；Vina平均使用 {t['vina_average_core_equivalents']:.2f} 核",f"全系统CPU平均/峰值 {t['system_cpu_mean_percent']:.1f}% / {t['system_cpu_peak_percent']:.1f}%",f"swap开始/结束/峰值 {t['swap_before_GiB']:.3f}/{t['swap_after_GiB']:.3f}/{t['swap_peak_GiB']:.3f}GiB；系统swapout增加 {t['system_swapout_delta_bytes']} bytes",'共晶对照RMSD Å：'+str([round(c['top_rmsd_A'],3) for c in t['control']['checks']]),'']
lines += [f"验证：66/66任务成功，{summary['total_poses_checked']}个输出构象通过身份、坐标和文件哈希检查。",'注意：准备耗时仅含补原子/加氢、局部松弛及Meeko转换，不含下载、脚本开发、失败尝试和人工排查。监测全过程仅含对接与输出校验。','','限制：',*summary['limitations']]
(out/'report.txt').write_text('\n'.join(lines)+'\n');print('\n'.join(lines[:30]))
