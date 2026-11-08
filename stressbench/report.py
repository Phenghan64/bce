"""Evidence-first capacity report; supports interim and final generation."""
import datetime,hashlib,html,json,math,re,statistics,time
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,OUT,atomic,utc,sha

def read(path,default=None):
    try:return json.loads(Path(path).read_text())
    except (FileNotFoundError,json.JSONDecodeError):return default
def h(x):return html.escape(str(x))
def table(head,rows):return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+h(x)+'</th>' for x in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def pct(a,q):return float(np.percentile(a,q)) if a else None
def fmt(x,n=2):return '—' if x is None else f'{x:.{n}f}'
def samples(path):
    if not path.exists():return []
    result=[]
    for line in path.read_text().splitlines():
        try:result.append(json.loads(line))
        except json.JSONDecodeError:pass
    return result
def summarize_resources(folder):
    r=read(folder/'resources.json',{});s=samples(folder/'resources.jsonl')
    if not s:return r
    cpu=[x['system_cpu_percent'] for x in s];gpu=[x['gpu_system'].get('Device Utilization %',0) for x in s]
    return {**r,'sample_count':len(s),'mean_system_cpu_percent':float(np.mean(cpu)),'p95_system_cpu_percent':pct(cpu,95),'fraction_samples_CPU_ge90pct':float(np.mean(np.array(cpu)>=90)),'mean_system_gpu_percent':float(np.mean(gpu)),'p95_system_gpu_percent':pct(gpu,95),'thermal_max':max((x['thermal_state'] or 0) for x in s),'observed_wall_seconds':s[-1]['elapsed_seconds']}

progress=read(OUT/'progress.json',{});completion=read(OUT/'completion.json');machine=read(OUT/'machine.json',{})
resources=[];timeline=[]
paths=sorted(OUT.glob('*/resources.jsonl'))+sorted(OUT.glob('*/*/resources.jsonl'))
for p in paths:
    resource=summarize_resources(p.parent);resource['name']=str(p.parent.relative_to(OUT));resources.append(resource)
    for s in samples(p):s['phase']=resource['name'];timeline.append(s)
timeline.sort(key=lambda s:s['utc']);origin=datetime.datetime.fromisoformat(timeline[0]['utc']).timestamp() if timeline else time.time()
calibration=[]
for folder in sorted((OUT/'cpu-calibration').glob('*')):
    batch=read(folder/'batch/summary.json');r=summarize_resources(folder)
    partial=False
    if not batch:
        if not (folder/'resources.json').exists():continue
        partial=True
        w,c,e=map(int,re.match(r'(\d+)x(\d+)-E(\d+)',folder.name).groups())
        partial_jobs=[read(p,{}) for p in (folder/'batch/jobs').glob('*/result.json')]
        batch={'parameters':{'workers':w,'cpu':c,'exhaustiveness':e},'results':partial_jobs,'completed':sum(x.get('status')=='completed' for x in partial_jobs),'failed':sum(x.get('status')!='completed' for x in partial_jobs),'wall_seconds':r['wall_seconds'],'fresh_jobs_per_hour':None}
    jobs=batch['results'];valid=[x for x in jobs if x['status']=='completed' and not x.get('cache_hit')];latency=[x['job_wall_seconds'] for x in valid]
    # A process interrupted by a stage deadline is censored, not a failed QC
    # result. Keep never-started panel members visible as well.
    explicit_failures=sum(x.get('status')=='failed' for x in jobs)
    censored=sum(x.get('status')=='running' for x in jobs)
    unstarted=max(0,48-len(jobs)) if partial else 0
    calibration.append({'name':folder.name,'partial_or_censored':partial or batch['failed']>0,'actual_active_workers_peak':batch.get('peak_actual_active_workers'),'workers':batch['parameters']['workers'],'cpu':batch['parameters']['cpu'],'exhaustiveness':batch['parameters']['exhaustiveness'],'completed':batch['completed'],'failed':batch['failed'],'wall_seconds':batch['wall_seconds'],'observed_completion_rate_per_hour':batch['fresh_jobs_per_hour'],'jobs_per_hour':batch['fresh_jobs_per_hour'] if not batch['failed'] else None,'latency_p50':pct(latency,50),'latency_p95':pct(latency,95),'system_cpu_mean':r.get('mean_system_cpu_percent'),'CPU_ge90_fraction':r.get('fraction_samples_CPU_ge90pct'),'rss_GiB':r.get('process_tree_peak_rss_bytes',0)/2**30,'thermal_max':r.get('thermal_max'),'guard_stop':r.get('guard_stop')})
    calibration[-1].update(failed=explicit_failures,censored=censored,unstarted=unstarted,partial_or_censored=partial or explicit_failures>0 or censored>0)
    reasons={j.get('error','unknown') for j in jobs if j.get('status')=='failed'}
    calibration[-1]['error_counts']={reason:sum(j.get('error','unknown')==reason and j.get('status')=='failed' for j in jobs) for reason in reasons}
    if partial or explicit_failures or censored:calibration[-1]['jobs_per_hour']=None
calibration.sort(key=lambda x:(x['exhaustiveness'],x['workers']))
extended=[]
for folder in sorted((OUT/'cpu-extended').glob('*')):
    if not folder.is_dir():continue
    batch=read(folder/'batch/summary.json');r=summarize_resources(folder)
    if not batch:continue
    extended.append({'name':folder.name,'parameters':batch['parameters'],'completed':batch['completed'],'failed':batch['failed'],'wall_seconds':batch['wall_seconds'],'jobs_per_hour':batch['fresh_jobs_per_hour'] if not batch['failed'] and not r.get('guard_stop') else None,'resources':r})
gpu=read(OUT/'gpu-results.json',[])
for path in sorted((OUT/'gpu-extended').glob('*/validation.json')):
    case=read(path)
    if case:gpu.append(case)
for case in gpu:
    live=samples(OUT/case['name']/'model/mps-memory.jsonl')
    case['sampled_driver_peak_GiB']=max((s['driver_bytes'] for s in live),default=0)/2**30 if live else None
    case['single_sequence_input']='myosin_prefix_' in Path(case['query']).name or 'copies_' in Path(case['query']).name
soak=read(OUT/'cpu-soak/batch/summary.json');soak_resource=summarize_resources(OUT/'cpu-soak') if (OUT/'cpu-soak').exists() else {}
soak_analysis={};drift=[]
if soak:
    jobs=[r for r in soak['results'] if r['status']=='completed'];duration=soak['parameters']['duration'];third=duration/3
    for ligand in sorted({r['spec']['ligand_id'] for r in jobs}):
        early=[r['job_wall_seconds'] for r in jobs if r['spec']['ligand_id']==ligand and r['start_offset_seconds']<third]
        late=[r['job_wall_seconds'] for r in jobs if r['spec']['ligand_id']==ligand and 2*third<=r['start_offset_seconds']<duration]
        if early and late:drift.append({'ligand':ligand,'early_n':len(early),'late_n':len(late),'early_median_seconds':statistics.median(early),'late_median_seconds':statistics.median(late),'late_early_ratio':statistics.median(late)/statistics.median(early)})
    periods=[]
    for i in range(3):
        rr=[r for r in jobs if i*third<=r['finish_offset_seconds']<(i+1)*third]
        periods.append({'period':i+1,'minutes':[i*third/60,(i+1)*third/60],'completed':len(rr),'jobs_per_hour':len(rr)*3600/third,'complexity_counts':{str(g):sum(r['complexity_group']==g for r in rr) for g in (0,1,2)}})
    repeats=defaultdict(list)
    for r in jobs:repeats[r['spec']['ligand_id']].append(r['scores'][0])
    soak_analysis={'completed':len(jobs),'failed':soak['failed'],'wall_seconds':soak['wall_seconds'],'target_duration_seconds':duration,'jobs_per_hour':soak['fresh_jobs_per_hour'],'latency_p50':pct([r['job_wall_seconds'] for r in jobs],50),'latency_p95':pct([r['job_wall_seconds'] for r in jobs],95),'matched_ligands':len(drift),'median_matched_late_early_latency_ratio':statistics.median(r['late_early_ratio'] for r in drift) if drift else None,'periods':periods,'max_same_ligand_score_range_kcalmol':max((max(v)-min(v) for v in repeats.values()),default=None),'resources':soak_resource}
    atomic(OUT/'soak-ligand-drift.json',drift)
mixed=[]
for name in ('mixed-full','mixed-reduced'):
    if not (OUT/name).exists():continue
    r=summarize_resources(OUT/name);summary=read(OUT/name/'work/summary.json',{});current=read(OUT/name/'work/progress.json',{})
    vb=read(OUT/name/'work/vina/summary.json',{})
    if not vb and (OUT/name/'resources.json').exists():
        jj=[read(f,{}) for f in (OUT/name/'work/vina/jobs').glob('*/result.json')]
        last_cycle=max((j['spec']['cycle'] for j in jj),default=-1)
        vb={'results':jj,'completed':sum(j.get('status')=='completed' for j in jj),
            'failed':sum(j.get('status')=='failed' for j in jj),
            'censored':sum(j.get('status')=='running' for j in jj),
            'unstarted':max(0,(last_cycle+1)*96-len(jj)),
            'wall_seconds':r.get('wall_seconds'),'fresh_jobs_per_hour':None,
            'partial_or_censored':True}
    command=r.get('command',[])
    window=float(command[command.index('--duration')+1]) if '--duration' in command else None
    if vb and window:
        elapsed=min(window,r.get('wall_seconds',window))
        delivered=sum(j.get('status')=='completed' and j['finish_offset_seconds']<=elapsed for j in vb['results'])
        vb.update(planned_window_seconds=window,window_completed=delivered,
                  window_completion_rate_per_hour=delivered*3600/elapsed)
    fold_runs=[]
    for f in sorted((OUT/name/'work').glob('fold-*/benchmark.json')):
        b=read(f,{})
        fold_runs.append({'folder':f.parent.name,'completed':b.get('completed'),
                         'cli_seconds':b.get('cli_seconds'),'query':b.get('query'),
                         'mps_driver_peak_GiB':b.get('sampled_mps_driver_peak_bytes',0)/2**30})
    sequence_runs=[{'folder':f.parent.name,**read(f,{})} for f in sorted((OUT/name/'work').glob('sequence-*/timing.json'))]
    matched=[]
    if soak and vb:
        complete_cycles={cycle for cycle in {j['spec']['cycle'] for j in vb['results']}
                         if sum(j['status']=='completed' and j['spec']['cycle']==cycle for j in vb['results'])==96}
        for ligand in sorted({j['spec']['ligand_id'] for j in vb.get('results',[])}):
            alone=[j['job_wall_seconds'] for j in soak['results'] if j['status']=='completed' and j['spec']['ligand_id']==ligand]
            together=[j['job_wall_seconds'] for j in vb['results'] if j['status']=='completed' and j['spec']['ligand_id']==ligand and j['spec']['cycle'] in complete_cycles]
            if alone and together:matched.append(statistics.median(together)/statistics.median(alone))
    mixed.append({'name':name,'resources':r,'summary':summary,'progress':current,
                  'vina':{k:vb.get(k) for k in ('completed','failed','censored','unstarted','partial_or_censored','wall_seconds','fresh_jobs_per_hour','planned_window_seconds','window_completed','window_completion_rate_per_hour')},
                  'matched_vina_ligands':len(matched),
                  'median_vina_mixed_alone_latency_ratio':statistics.median(matched) if matched else None,
                  'fold_runs':fold_runs,'sequence_runs':sequence_runs,
                  'fold_fallback':read(OUT/name/'work/fold-fallback.json')})
recovery=read(OUT/'recovery-verdict.json',{})
if recovery:
    recovery={k:v for k,v in recovery.items() if k!='summary' and k!='resume_summary'}|{'completed':read(OUT/'recovery/batch/summary.json',{}).get('completed'),'failed':read(OUT/'recovery/batch/summary.json',{}).get('failed'),'resume_cache_hits':read(OUT/'recovery/batch/resume-summary.json',{}).get('cache_hits'),'resume_seconds':read(OUT/'recovery/batch/resume-summary.json',{}).get('wall_seconds')}
supplement_completion=read(OUT/'supplement-completion.json')
if recovery and 'completed' in read(OUT/'recovery-verdict.json',{}):
    verdict=read(OUT/'recovery-verdict.json');recovery={k:v for k,v in verdict.items() if k not in ('summary','resume_summary')}
raw_recovery=read(OUT/'recovery-verdict.json',{})
if recovery:
    resume_wall=raw_recovery.get('resume_resources',{}).get('wall_seconds')
    if resume_wall is None:resume_wall=read(OUT/'recovery-resume/resources.json',{}).get('wall_seconds')
    reused=raw_recovery.get('summary',{}).get('cache_hits',0)
    recovery_text=f"工作进程终止后重试成功；恢复计划完成{recovery.get('completed')}个任务、失败{recovery.get('failed')}个。续算复用了{reused}个已完成任务。最后重新提交时，{recovery.get('resume_cache_hits')}/96个任务通过输出文件哈希核对并命中缓存；恢复核验阶段墙钟{fmt(resume_wall)}秒，批处理内部计时{fmt(recovery.get('resume_seconds'),3)}秒。原主阶段的时间截断记录保留。"
else:recovery_text='尚未完成。'
now=datetime.datetime.now(datetime.timezone.utc)
boundary_completion=read(OUT/'boundary-completion.json')
end_record=boundary_completion or supplement_completion or completion
summary={'generated_utc':now.isoformat(),'campaign_status':completion or progress,'supplement_status':supplement_completion,'boundary_followup_status':boundary_completion,'machine':machine,'elapsed_from_baseline_start_seconds':(datetime.datetime.fromisoformat(end_record['utc']).timestamp() if end_record else now.timestamp())-origin,'CPU_calibration':calibration,'CPU_extended':extended,'CPU_extended_matched_comparison':read(OUT/'cpu-extended/matched-e64-comparison.json'),'GPU_cases':[{k:v for k,v in r.items() if k!='error_tail'} for r in gpu],'soak':soak_analysis,'mixed':mixed,'mixed_output_validation':read(OUT/'mixed-output-validation.json'),'recovery':recovery,'stage_resources':resources,'post_run_health':read(OUT/'post-run-health.json'),'scope':'Application capacity on this host and fixed workload, not universal physical maximum or biological validation.'}
atomic(OUT/'summary.json',summary)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def save(name):plt.savefig(OUT/f'{name}.svg',bbox_inches='tight');plt.savefig(OUT/f'{name}.png',bbox_inches='tight',dpi=160);plt.close()
if timeline:
    x=np.array([(datetime.datetime.fromisoformat(s['utc']).timestamp()-origin)/3600 for s in timeline]);fig,axs=plt.subplots(3,1,figsize=(11,8),sharex=True)
    x[np.where(np.diff(x,prepend=x[0])>30/3600)[0]]=np.nan
    axs[0].plot(x,[s['system_cpu_percent'] for s in timeline],label='System CPU',lw=.7);axs[0].plot(x,[s['gpu_system'].get('Device Utilization %',0) for s in timeline],label='System GPU',lw=.7);axs[0].set_ylabel('Utilization (%)');axs[0].set_ylim(-3,103);axs[0].legend(loc='upper left')
    axs[1].plot(x,[s['process_tree_rss_bytes']/2**30 for s in timeline],label='Test tree RSS',lw=.8);axs[1].plot(x,[s['system_available_bytes']/2**30 for s in timeline],label='System available RAM',lw=.8);axs[1].plot(x,[s['swap_used_bytes']/2**30 for s in timeline],label='System swap used',lw=.8);axs[1].set_ylabel('GiB');axs[1].legend(loc='upper left')
    axs[2].step(x,[s['thermal_state'] if s['thermal_state'] is not None else np.nan for s in timeline],where='post',color='#b36724');axs[2].set_yticks([0,1,2,3],['Nominal','Fair','Serious','Critical']);axs[2].set_xlabel('Hours since baseline start');axs[2].set_ylabel('Thermal pressure');fig.suptitle('Measured load over time — system metrics include other applications');fig.tight_layout();save('resource-timeline')
if calibration:
    fig,ax=plt.subplots(figsize=(10,4));layouts=[(1,12),(2,6),(3,4),(6,2),(12,1),(16,1)]
    for e in (8,32,64):
        ax.plot(range(6),[next((r['jobs_per_hour'] for r in calibration if (r['workers'],r['cpu'],r['exhaustiveness'])==(w,c,e)),np.nan) for w,c in layouts],marker='o',label=f'E{e}')
    ax.set_xticks(range(6),[f'{w} jobs × {c} threads' for w,c in layouts],rotation=15,ha='right');ax.set_ylabel('QC-passing jobs / hour');ax.set_title('Fixed balanced ligand panel; complete cycles; identical seeds');ax.legend();fig.tight_layout();save('cpu-scaling')
if drift:
    fig,ax=plt.subplots(figsize=(10,4));ratios=sorted(r['late_early_ratio'] for r in drift);ax.plot(range(1,len(ratios)+1),ratios,'o',ms=3);ax.axhline(1,color='#888',ls='--');ax.set_xlabel('Matched ligand (sorted by latency ratio)');ax.set_ylabel('Late / early median job latency');ax.set_title('Same-ligand comparison during sustained CPU load; >1 means slower');fig.tight_layout();save('soak-drift')
cpu_table=table(['布局','实际并发峰值','强度','完成/失败/截断/未开始','任务/小时','延迟中位/P95 秒','CPU均值','CPU≥90%比例','RSS峰值GiB'],[[f"{r['workers']}×{r['cpu']}",r.get('actual_active_workers_peak','—'),r['exhaustiveness'],f"{r['completed']}/{r['failed']}/{r['censored']}/{r['unstarted']}",fmt(r['jobs_per_hour'],0),f"{fmt(r['latency_p50'])}/{fmt(r['latency_p95'])}",fmt(r['system_cpu_mean'],1)+'%',fmt((r['CPU_ge90_fraction'] or 0)*100,1)+'%',fmt(r['rss_GiB'])] for r in calibration])
extended_table=table(['补充配置','完成/失败','实际秒数','任务/小时','阶段停止原因'],[[h(r['name']),f"{r['completed']}/{r['failed']}",fmt(r['wall_seconds']),fmt(r['jobs_per_hour'],0),h(r['resources'].get('guard_stop') or '正常完成')] for r in extended])
gpu_table=table(['任务','结果','耗时 秒','MPS采样峰值 GiB','RSS采样峰值 GiB','swap峰值 GiB','备注'],[[h(r['name']),'通过' if r['success'] else '失败/中止',fmt(r['resources']['wall_seconds']),fmt(r.get('sampled_driver_peak_GiB')),fmt(r['resources']['process_tree_peak_rss_bytes']/2**30),fmt(r['resources']['swap_peak_bytes']/2**30),h(r.get('failure_class','输出坐标及残基数通过检查'))] for r in gpu])
totalhours=summary['elapsed_from_baseline_start_seconds']/3600;maxcpu=max((r.get('p95_system_cpu_percent',0) or 0 for r in resources),default=0);maxrss=max((r.get('process_tree_peak_rss_bytes',0) for r in resources),default=0)/2**30
findings=[]
eligible=[r for r in calibration if r['exhaustiveness']==32 and r['jobs_per_hour'] is not None]
if eligible:
    best=max(eligible,key=lambda r:r['jobs_per_hour'])
    findings.append(f"独立Vina：本面板E32实测最高完整吞吐为{best['jobs_per_hour']:.0f}次/小时，布局{best['workers']}×{best['cpu']}。多个低并发布局差距较小，不能宣称通用最优；12/16进程的高并发增加了尾部超时。")
if soak_analysis:
    findings.append(f"持续承载：{soak_analysis['wall_seconds']/60:.1f}分钟完成{soak_analysis['completed']}次对接、失败{soak_analysis['failed']}次，约{soak_analysis['jobs_per_hour']:.0f}次/小时。同配体后/前耗时比中位数{soak_analysis['median_matched_late_early_latency_ratio']:.3f}。这是96个分子的重复负载。")
sizes=[r for r in gpu if r['success'] and 'myosin_prefix_' in Path(r['query']).name and r['samples']==1]
if sizes:
    largest=max(sizes,key=lambda r:r['validation']['expected_residues'])
    findings.append(f"结构预测：本次最大完成长度构建体为{largest['validation']['expected_residues']}残基，阶段耗时{largest['resources']['wall_seconds']/60:.1f}分钟（{largest['name']}）。这不是所有序列的残基上限；需结合失败原因、MSA及链数判断。")
for m in mixed:
    if m['name']=='mixed-full' and m['vina'].get('window_completion_rate_per_hour') is not None:
        findings.append(f"混合争用：计划45分钟窗口交付{m['vina']['window_completed']}次对接，折合{m['vina']['window_completion_rate_per_hour']:.0f}次/小时；完整周期中同配体耗时约为独立运行的{fmt(m['median_vina_mixed_alone_latency_ratio'],2)}倍。整批排队导致收尾超过55分钟时限，原配置不能记为全批次通过，也不能据此宣称硬件无法完成。")
        findings.append(f"混合内存：全机swap峰值{m['resources'].get('swap_peak_bytes',0)/2**30:.2f}GiB，最低可用内存{m['resources'].get('system_min_available_bytes',0)/2**30:.2f}GiB。能运行伴随着明显资源争用，不能称为没有压力。")
ooms=[r['name'] for r in gpu if r.get('failure_class')=='mps_allocator_oom']
if ooms:findings.append('已观测MPS分配器内存失败：'+', '.join(ooms)+'。此限制属于当前分配预算与软件设置，未禁用系统保护。')
elif boundary_completion:
    findings.append('边界仍有未决项：本轮未观测到MPS OOM，1536残基是在约16分钟预算到期时被截断。尚未定位物理内存的绝对上限，不能把最大已完成的768残基当成硬件上限。')
findings_html='<ul>'+''.join('<li>'+h(x)+'</li>' for x in findings)+'</ul>'
timing_rows=[]
for label,names in [('环境基线（复用）',('baseline',)),('CPU校准',('cpu-calibration/',)),('GPU主阶段',('gpu-size/','gpu-msa/','gpu-samples/','gpu-chains/','gpu-repeat/')),('CPU持续负载',('cpu-soak',)),('混合原配置',('mixed-full',)),('混合降载',('mixed-reduced',)),('恢复及续算',('recovery',)),('CPU延时复测',('cpu-extended/',)),('GPU延时复测',('gpu-extended/',))]:
    sec=sum(r.get('wall_seconds',0) for r in resources if any(r['name'].startswith(n) for n in names))
    if sec:timing_rows.append([label,fmt(sec/60), '各阶段监测墙钟之和'])
timing_table=table(['项目','分钟','计时口径'],timing_rows)
isdone=completion is not None;state='已结束' if isdone else '运行中：阶段性报告';selected=read(OUT/'selected-cpu-configuration.json',{});gpuselected=read(OUT/'selected-gpu-configuration.json',{})
soak_text=f"持续负载完成 {soak_analysis.get('completed',0)} 个任务，实际进程时间 {soak_analysis.get('wall_seconds',0)/60:.1f} 分钟，失败 {soak_analysis.get('failed',0)}。同配体后段/前段延迟比中位数为 {fmt(soak_analysis.get('median_matched_late_early_latency_ratio'),3)}（大于1表示变慢）。" if soak_analysis else '持续90分钟负载尚未完成，暂不下稳定性结论。'
mix_table=table(['阶段','墙钟 分钟','退出码','保护中止','峰值RSS GiB','最低可用 GiB','swap峰值 GiB'],[[h(m['name']),fmt(m['resources'].get('wall_seconds',m['resources'].get('observed_wall_seconds',0))/60),m['resources'].get('returncode','运行中'),h(m['resources'].get('guard_stop') or '无'),fmt(m['resources'].get('process_tree_peak_rss_bytes',0)/2**30),fmt(m['resources'].get('system_min_available_bytes',0)/2**30),fmt(m['resources'].get('swap_peak_bytes',0)/2**30)] for m in mixed])
mix_work_table=table(['阶段','Vina完成/失败/截断/未开始','窗口内完成/小时','同配体混合/独立耗时比','完成预测数','预测CLI中位秒','完整序列流程次数','序列计算中位秒'],[[h(m['name']),f"{m['vina'].get('completed')}/{m['vina'].get('failed')}/{m['vina'].get('censored') or 0}/{m['vina'].get('unstarted') or 0}",fmt(m['vina'].get('window_completion_rate_per_hour'),0),fmt(m.get('median_vina_mixed_alone_latency_ratio'),3),sum(bool(f['completed']) for f in m['fold_runs']),fmt(statistics.median([f['cli_seconds'] for f in m['fold_runs'] if f['completed']])) if any(f['completed'] for f in m['fold_runs']) else '—',len(m['sequence_runs']),fmt(statistics.median([r['total_compute_seconds'] for r in m['sequence_runs']])) if m['sequence_runs'] else '—'] for m in mixed])
body=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>本机性能极限测试</title><style>body{{background:#f2f5f5;color:#183039;font:16px/1.75 -apple-system,"PingFang SC",sans-serif;margin:0}}main{{max-width:1120px;padding:40px 26px;margin:auto}}h1{{font-size:32px}}h2{{font-size:23px}}section{{background:white;border:1px solid #d6e0e0;border-radius:10px;padding:24px;margin:24px 0}}.tag{{color:#0b7777}}p,li{{overflow-wrap:anywhere}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}}.cards div{{padding:20px;background:#e1eeeb;border-radius:8px}}.cards strong{{font-size:25px;display:block}}.cards small{{font-size:13px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:10px;border-bottom:1px solid #d6e0e0;text-align:left}}th{{white-space:nowrap;background:#eaf2ef}}img{{width:100%;height:auto}}a{{color:#087c78}}.note{{background:#fff3df;border-left:4px solid #c38a39;padding:16px}}code{{font-size:12px;overflow-wrap:anywhere}}@media(max-width:650px){{main{{padding:16px}}.cards{{grid-template-columns:1fr 1fr}}section{{padding:16px}}}}</style></head><body><main><p class="tag">{h(state)} · {h(progress.get('phase',''))}</p><h1>本机性能极限与持续承载测试</h1><p>Apple M6 · 12 核 · 16 GiB · 官方 OpenFold3 / PyTorch MPS + Vina 1.2.7</p><div class="cards"><div><strong>{totalhours:.2f} 小时</strong><small>从基线开始的实际时间</small></div><div><strong>{len(calibration)}/18</strong><small>CPU配置执行数</small></div><div><strong>{len(gpu)}</strong><small>结构预测容量试验</small></div><div><strong>{maxcpu:.1f}%</strong><small>阶段最高CPU P95利用率</small></div></div>
<section><h2>实测结论与使用建议</h2>{findings_html}<p>建议按资源分别排队：独立对接可从3×4、E32起步；大结构预测保持单个GPU任务，混合运行预留内存并按实际交付率安排时限。给对接队列设置有界派发与尾部收尾预算，避免截止前一次排入整个96任务面板。降载方案只完成了较短并行窗口，且使用流式简化序列处理，不等于完整功能流程的长时验证。</p><p>这些结果验证本机应用负载与调度边界，不改变功能肽实验排序验证的科学结论，也不证明企业业务SLA或任意蛋白/配体的吞吐。</p></section>
<section><h2>实际耗时</h2>{timing_table}<p>主控制器墙钟：{fmt((completion or {}).get('wall_seconds',progress.get('elapsed_seconds',0))/60)}分钟；补充流程墙钟：{fmt((supplement_completion or {}).get('wall_seconds',0)/60)}分钟；追加边界探测墙钟：{fmt((boundary_completion or {}).get('wall_seconds',0)/60)}分钟。页面总时间从10分钟基线开始计算，包含被排除的先导测试、重启间隔与阶段切换；不把这些时间加入正式配置吞吐。主流程6小时检查在阶段边界执行，已启动阶段可能越过该时点；这里展示实际值。</p></section>
<section><h2>这次要回答什么</h2><p>测量实际施加的负载、吞吐平台、内存边界、长时间性能变化及混合任务失败点。对接库为96个真实NCI小分子，校准取其中48个；持续负载重复计算这些输入，不代表同等数量的唯一候选物或万级功能肽对接。数据规模和启动的线程数不能代替真实资源利用率；同一硬件也没有脱离软件与任务类型的单一“极限”。</p><p class="note">{'所有已执行阶段的结果如下。完成计划不表示所有配置成功，必须同时查看失败与保护中止。' if isdone else '测试仍在进行，下面只报告已获得的证据。短时满载不能提前作为持续稳定性的结论。'}</p><p>使用现有官方 OpenFold3 的 MPS 后端，未安装 MLX。初始 MPS 分配预算提高到推荐上限的100%；OOM可能首先触及分配器预算，报告不将它等同于物理内存绝对上限。基线包含正常后台应用及程序准备活动，不是隔离空载环境。</p></section>
<section><h2>实际资源曲线</h2><img src="resource-timeline.svg" alt="Resource timeline"><p>系统/RSS约每2秒、GPU利用率约每10秒、MPS分配约每1秒采样，瞬时峰值可能漏采。CPU为全机12核平均；GPU来自ioreg，其空闲时的读数可能保留末次统计值，不能单凭此值判断是否仍有后台任务；热压力来自NSProcessInfo。未获取需要管理员权限的温度和功耗。RSS与MPS驱动内存有重叠，不能相加。系统swap和可用内存包含其他应用；超过30秒无采样的区间不连线。</p></section>
<section><h2>CPU：满载是否换来更多吞吐</h2><p>固定48个分层配体、同一个受体/搜索框和种子；各配置至少120秒并完成完整面板周期。测试1×12、2×6、3×4、6×2、12×1及16×1超额并发，分别E8/E32/E64。每个完成任务检查分子连接关系、有限坐标、姿态数和分数。失败原因逐配置写入汇总数据的error_counts；阶段截止仍在运行的任务记为截断，与输出检查失败分开。吞吐空值表示未完整通过，不用删掉慢任务后的完成率参与配置排名。</p>{cpu_table}{'<img src="cpu-scaling.svg" alt="CPU scaling">' if calibration else ''}<p>选定布局：{selected.get('workers','—')}个任务×{selected.get('cpu','—')}线程（独立E32校准{fmt(selected.get('calibration_jobs_per_hour'),0)}次/小时）。它是在该校准面板上的实测选择；几个通过的E32布局差距较小，未证明统计显著优势，不是所有分子和所有受体的通用最佳配置。</p></section>
<section><h2>E64延长时限复测</h2><p>主测试中600秒阶段上限截断的结果保留。另用全新输出、同一48配体面板、1×12布局及E64，将单任务预算放宽至900秒、阶段预算至1800秒，检验时间预算截断是否能够通过延时完成。以下数据不混入原18组校准。</p>{extended_table if extended else '<p>尚未完成。</p>'}<p>同配体与原E64测试的<a href="cpu-extended/matched-e64-comparison.json">逐项匹配比较</a>单独保留，原测试未完成的那一个配体不参与速度比计算。</p></section>
<section><h2>GPU与统一内存：在哪里首先失败</h2><p>长度阶梯使用真实鱼肌球蛋白的嵌套序列前缀，前缀及多链复制任务均仅输入单序列。MSA深度使用完整446残基TmGH1的缓存比对，输入分别核实为64、256、1024个不同序列行，数据与模型采样上限核实生效；本轮未单独记录模型前向张量行数。扩散样本数在256残基前缀上变化，多链复制使用TmGH1。前缀和复制体是性能构建体，不作为真实蛋白功能或天然复合物结论。模型保留3次recycle、200步diffusion及low_mem设置。</p>{gpu_table}<p>主阶段首次未完成长度：{h(gpuselected.get('first_failed_size','尚未确定'))}；主阶段最大完成长度：{h(gpuselected.get('largest_passed_size','尚未确定'))}。表中gpu-extended任务是主流程结束后的独立补测，不混入原GPU阶段预算。768残基若完成，则用30分钟GPU补测预算中剩余的时间增加1536残基压力点；属于自适应边界探测，未完成的任务不提供硬件残基上限。长度之外，MSA、样本数、链数和当时背景内存都会改变边界，不能只给一个残基上限。</p></section>
<section><h2>90分钟持续计算</h2><p>{h(soak_text)}</p>{'<img src="soak-drift.svg" alt="Same-ligand latency drift">' if drift else ''}<p>同配体前后段比较用于减小不同分子复杂度的影响。重复种子用于性能重复性，不是独立实验样本；不凭速度下降单独宣称热降频。</p><pre>{h(json.dumps(soak_analysis.get('periods',[]),ensure_ascii=False,indent=2))}</pre></section>
<section><h2>混合业务负载与恢复</h2><p>同时运行Vina、OpenFold3和此前完整虚拟酶解/功能筛选。若触及保护边界，保留失败配置，再测试缩小预测和流式序列处理的降载方案；降载通过不能冒充原配置通过。</p>{mix_table}{mix_work_table}<p>窗口内完成率按计划并行窗口内实际交付数计算，包含排队与未完成任务的影响；不把超时配置当作完整面板通过。混合/独立耗时比仅采用完整完成的96配体周期，与独立持续负载匹配；大于1代表变慢。序列流程耗时为流水线自身记录的计算时间，预测CLI耗时包含模型加载；不与另一种计时口径直接比较。</p><p>混合输出核验：<a href="mixed-output-validation.json">逐个预测的坐标/残基数检查，以及完整序列流程的输出计数核对</a>（补充流程开始时生成）。</p><p>恢复检查：{h(recovery_text)} <a href="recovery-verdict.json">查看恢复明细</a></p><p>批处理内部计时不包含Python与依赖加载；完整阶段墙钟包含这些成本。故障注入仅向本测试拥有的一个Vina工作进程发送SIGTERM；同一任务重试后，再重新提交固定计划，核对完成文件哈希并检查缓存命中。</p></section>
<section><h2>停止条件与可复现证据</h2><p>可用内存低于384 MiB持续12秒、阶段swap增长超过6 GiB、全机swap高于8 GiB且继续增长、严重到critical的热压力持续30秒、磁盘不足20 GiB或阶段超时，会停止该测试进程组。单配体预算240秒、CPU校准配置预算600秒，单步GPU测试通常预算600秒；超时表示本轮时间预算，不代表硬件绝对无法完成。不会终止其他应用，也不通过禁用系统保护追求死机。</p><p>运行后进程清理及系统快照：<a href="post-run-health.json">post-run-health.json</a>。</p><p>运行状态：<a href="progress.json">progress.json</a>；阶段明细：<a href="stages.json">stages.json</a>；<a href="summary.json">汇总数据</a>；<a href="protocol.json">协议</a>；<a href="input-manifest.json">输入哈希</a>；<a href="software-manifest.json">主流程启动时程序哈希</a>；<a href="software-final-manifest.json">最终程序快照哈希</a>；<a href="campaign-amendment.json">先导排除说明</a>；<a href="pre-recovery-adjustment.json">恢复计划修正</a>；<a href="../../stressbench/README.txt">操作与复现说明</a>。</p><p>所有容量建议仅限本轮软件、输入、参数和后台环境。未测温度及功耗，未做操作系统交互延迟测试；热压力正常不能单独证明没有频率调整。</p></section><footer>报告生成：{h(now.isoformat())}。主流程状态：{h((completion or {}).get('status','running'))}；补充状态：{h((supplement_completion or {}).get('status','pending'))}。{h((completion or {}).get('reason',''))}</footer></main></body></html>'''
summary['findings']=findings
summary['timing_breakdown']=timing_rows
atomic(OUT/'summary.json',summary)
(OUT/'report.html').write_text(body)
print(json.dumps({'report':str(OUT/'report.html'),'status':summary['campaign_status'],'elapsed_hours':totalhours,'cpu_configurations':len(calibration),'gpu_cases':len(gpu)},ensure_ascii=False,default=str)[:1600])
