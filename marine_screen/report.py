"""Chinese offline report and scientific figure built from saved results."""
import datetime,html,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare import OUT

def load(name):return json.loads((OUT/name).read_text())
def esc(x):return html.escape(str(x))
def table(headers,rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in r)+'</tr>' for r in rows)+'</tbody></table></div>'

def main():
    s=load('summary.json');a=load('input-audit.json');v=load('validation.json');b=load('process-benchmark.json')
    res=load('resources/resources.json');schemes=load('schemes.json');sardine=load('sardine-validation.json')
    combos=load('combination-diagnostics.json');qa=load('verification.json');short=load('literature-shortlist.json')
    log=(OUT/'resources/run.log').read_text();maxrss=int(re.search(r'(\d+)\s+maximum resident set size',log)[1])
    names=[x['scheme'] for x in schemes if x['max_missed']==0]
    names_zh={'trypsin':'胰蛋白酶','chymotrypsin':'胰凝乳蛋白酶','pepsin_pH1.3':'胃蛋白酶 pH 1.3规则',
        'pepsin_pH2':'胃蛋白酶 pH 2.0规则','thermolysin':'嗜热菌蛋白酶','proteinase_K':'蛋白酶 K',
        'trypsin+chymotrypsin':'胰蛋白酶＋胰凝乳蛋白酶','pepsin_pH2+trypsin':'胃蛋白酶＋胰蛋白酶',
        'thermolysin+trypsin':'嗜热菌蛋白酶＋胰蛋白酶'}
    short_names=['Trypsin','Chymotrypsin','Pepsin pH1.3','Pepsin pH2','Thermolysin','Proteinase K','Tryp + Chymo','Pepsin + Tryp','Thermo + Tryp']
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160})
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    ax=axes[0,0];y=np.arange(len(names));n0=[x['unique_peptides'] for x in schemes if x['max_missed']==0];n1=[x['unique_peptides'] for x in schemes if x['max_missed']==1]
    ax.barh(y+.17,n0,.33,color='#087f8c',label='Complete cleavage')
    ax.barh(y-.17,n1,.33,color='#9cd6d5',label='0 or 1 missed site')
    ax.set_yticks(y,short_names);ax.invert_yaxis();ax.set_xlabel('Unique sequences per scheme (overlap across schemes)');ax.set_title('A. Screening scale');ax.legend(fontsize=8)
    ax=axes[0,1]
    for i,(item,label) in enumerate(zip(v,['Literature (n=494)','Candidate overlap (n=86)','Sardine primary (n=6)'])):
        m=item['methods']['heuristic'];rho=m['spearman'];lo,hi=m['bootstrap95']
        ax.errorbar(rho,2-i,xerr=[[rho-lo],[hi-rho]],fmt='o',color='#c65f35',capsize=4)
    ax.set_yticks([2,1,0],['Literature (n=494)','Candidate overlap (n=86)','Sardine primary (n=6)'])
    ax.axvline(0,ls='--',color='#aaa');ax.set_xlim(-1.1,1.1);ax.set_xlabel('Spearman rho vs -log10(IC50); 95% peptide bootstrap');ax.set_title('B. Heuristic ranking has no reliable advantage')
    ax=axes[1,0];rows=[r for r in b['rows'] if r['max_missed']==0 and r['tropomyosin']=='P13104' and r['scheme'] in b['compared_measured_enzymes']]
    ax.bar(['Pepsin','Thermolysin','Proteinase K'],[r['primary_density'] for r in rows],color=['#087f8c','#a5b7c1','#a5b7c1'])
    ax.set_ylabel('Score >=3 fragments / 1,000 source residues');ax.set_title('C. Mackerel process benchmark: wrong winner')
    ax.text(.03,.94,'Experimental winner: pepsin\nHeuristic winner: proteinase K',transform=ax.transAxes,va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    ax=axes[1,1];data=[json.loads(line) for line in (OUT/'resources/resources.jsonl').read_text().splitlines()]
    ax.plot([x['elapsed_seconds'] for x in data],[x['rss_sum_bytes']/2**20 for x in data],'o-',color='#087f8c')
    ax.set_xlabel('Seconds since monitor launch');ax.set_ylabel('Process-tree RSS sum (MiB), 1-second samples');ax.set_title('D. Resource use of the compute run')
    fig.savefig(OUT/'overview.png');fig.savefig(OUT/'overview.svg');plt.close(fig)
    resources_table=table(['项目','实测值','解释'],[
        ['电脑','Apple M6 · 12物理核 · 16 GiB · macOS27','本次只读核验'],
        ['核心脚本总耗时',f"{s['compute_seconds']:.2f}秒",'含酶解、性质、索引、回测和结果写出'],
        ['命令墙钟耗时','4.87秒','/usr/bin/time -l；不含资料下载、开发与人工核查'],
        ['进程树峰值 RSS',f"{res['process_tree_peak_RSS_bytes']/2**20:.1f} MiB",'每1秒采样；共享页可能重复计数、短峰可能漏采'],
        ['time记录最高 RSS',f'{maxrss/2**20:.1f} MiB','time独立口径，不能与上行相加'],
        ['系统最低可用内存',f"{res['system_min_available_bytes']/2**30:.2f} GiB",'包含其他应用影响'],
        ['系统 swap',f"{res['swap_before']['used']/2**20:.0f} → {res['swap_after']['used']/2**20:.0f} MiB",'本次监测期间无增长'],
        ['完整候选与位置数据库',f"{s['sequence_database_bytes']/2**20:.1f} MiB",'SQLite，1144604条出现记录可逐条追溯'],
        ['结构预测／对接任务','0','此轮检验序列筛选，不将前轮OpenFold3/Vina表现外推到本轮']])
    comparison=[]
    for name in names:
        x=next(x for x in schemes if x['scheme']==name and x['max_missed']==0)
        y=next(x for x in schemes if x['scheme']==name and x['max_missed']==1)
        comparison.append([names_zh[name],f"{x['unique_peptides']:,}",f"{y['unique_peptides']:,}",x['reference_matches'],y['reference_matches']])
    validation_rows=[]
    for item,title in zip(v,['494条文献数据库记录','其中86条与本轮候选重合','沙丁鱼论文6条实测肽']):
        for key,label in [('heuristic','末端＋长度规则'),('shorter','优先短肽'),('gravy','疏水性 GRAVY')]:
            m=item['methods'][key];ci=m['bootstrap95']
            validation_rows.append([title,label,f"{m['spearman']:.3f}",f'[{ci[0]:.3f}, {ci[1]:.3f}]',
                f"{m['enrichment']:.2f}" if m['enrichment'] is not None else '不适用：无≤10µM记录'])
    benchmark_rows=[]
    for row in rows:
        benchmark_rows.append([names_zh[row['scheme']],row['products'],row['high_score'],f"{row['primary_density']:.2f}",
                              '实验最强，抑制率86.14%' if row['scheme']=='pepsin_pH2' else '其余5种实验酶解物无显著差异'])
    seqrows=[]
    for x in sardine:seqrows.append([x['sequence'],f"{x['ic50_mg_mL']:.2f} ± {x['sd_mg_mL']:.2f}",f"{x['ic50_uM']:.1f}",x['heuristic_score'],'是' if x['marine_MC01_present'] else '否'])
    selected=[]
    proteins={x['id']:x for x in load('proteins.json')}
    for x in short[:10]:
        p=x['provenance'][0];entry=proteins[p['protein']]['entries'][0]
        selected.append([x['sequence'],x['reference_ic50_uM'],p['protein'],entry['species'],names_zh[p['scheme']],
                         f"{p['start_1based']}–{p['end_1based_inclusive']}",p['missed']])
    text=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>海洋蛋白万级虚拟酶解与ACE功能肽筛选实测报告</title><style>
body{{margin:0;background:#f3f6f8;color:#20333e;font:16px/1.8 -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif}}main{{max-width:1120px;margin:auto;padding:38px 26px 80px}}h1{{font-size:32px;line-height:1.4}}h2{{font-size:23px;margin-top:38px}}h3{{font-size:18px}}.lead{{background:#fff;border-left:5px solid #087f8c;padding:20px 26px;border-radius:8px}}.warn{{background:#fff1e9;border-left:4px solid #c65f35;padding:14px 20px}}.muted{{color:#586c77;font-size:14px}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:22px 0}}.card{{background:white;padding:17px;border-radius:8px}}.big{{display:block;font-size:28px;font-weight:650;color:#087f8c}}table{{border-collapse:collapse;width:100%;font-size:14px;background:white}}td,th{{padding:10px 12px;border-bottom:1px solid #dde6eb;text-align:left;vertical-align:top}}th{{background:#e6eef2}}.scroll{{overflow-x:auto}}a{{color:#006c7a}}code{{background:#e7edf0;padding:2px 4px}}img{{width:100%;height:auto;background:white}}li{{margin:7px 0}}@media(max-width:760px){{.cards{{grid-template-columns:repeat(2,1fr)}}main{{padding:20px 14px}}h1{{font-size:25px}}}}@media print{{body{{background:white}}main{{padding:0}}table{{font-size:10px}}h2{{break-after:avoid}}.scroll{{overflow:visible}}}}
</style><main><p class="muted">实测日期：2026-10-22（北京时间） · 协议版本1.0 · 回顾性验证</p>
<h1>海洋蛋白万级虚拟酶解与 ACE 功能肽筛选</h1>
<div class="lead"><b>计算规模已达到目标，当前规则未通过研发决策验证。</b><br>本机完成699条真实数据库蛋白序列的9种方案筛选。完全切割已产生32,454条独立候选肽；加入一次漏切后为88,637条。计算资源充足，但末端残基与长度规则未可靠预测活性，也未找回鲭鱼案例中实验最好的胃蛋白酶。现阶段可作为候选生成与数据追溯工具，不能自动推荐最佳酶制剂或生产工艺。</div>
<div class="cards"><div class="card"><span class="big">32,454</span>完全切割独立肽</div><div class="card"><span class="big">88,637</span>0–1次漏切独立肽</div><div class="card"><span class="big">6,291</span>蛋白×方案任务</div><div class="card"><span class="big">4.87秒</span>计算命令墙钟时间</div></div>
<img src="overview.png" alt="候选规模、排序相关性、选酶失败和资源占用四项实测结果">
<h2>1. 本次跑了什么</h2>
<p>目标是食品与海洋副产物研发中的 ACE 抑制肽候选发现：真实蛋白序列 → 酶切与漏切情景 → 全库去重 → 分子量、电荷、疏水性与规则评分 → 文献标注 → 实验回测 → 候选和工艺审查。全流程在本地 CPU 上完成，没有将万条肽送去结构预测或对接。</p>
<p>UniProt完整流式查询返回3,644个条目。按蛋白名称限定肌肉结构、相关代谢蛋白与Ⅰ型胶原，去除2,871条不相关记录、22条片段；保留751个条目，其中仅9条为 reviewed。完全相同蛋白序列合并为699条，共347,906个氨基酸。未审阅条目多为基因组推断，不能称作699种经实验确认的原料蛋白。</p>
{table(['物种','保留条目数（序列去重前）'],sorted(a['species_entry_counts'].items()))}
<p>库包含大西洋鲑、鳕、鲭鱼和日本鲭。没有足够完整的欧洲沙丁鱼或澳洲鲭目标家族序列进入主库。多物种、异构体和旁系同源序列均保留来源；它们不代表同一份有已知含量的鱼肉。另做“物种＋蛋白名称选最长代表序列”敏感性分析，降低数据库异构体数量造成的偏置，但这仍不是序列聚类或丰度加权。</p>
<h2>2. 酶与组合方案</h2>
<p>5种不同蛋白酶、胃蛋白酶2种pH规则，共6个单酶条件和3个组合。规则来自冻结的Pyteomics／ExPASy源码快照；启用胰蛋白酶例外位点。候选长度限定2–20aa，主分析MC0为完全切割；MC0–1加入恰好一个内部漏切位点。两种情景不能相加当作实际产率。</p>
{table(['方案','MC0独立肽','MC0–1独立肽','MC0文献匹配','MC0–1文献匹配'],comparison)}
<p>各行有重叠，不能相加当独立候选数。全库记录了1,144,604条“肽—蛋白—位置—方案—漏切数”出现关系。组合采用原序列切点的并集；不能模拟加酶先后、酶失活或pH兼容性。胃蛋白酶＋胰蛋白酶应理解为切点情景，不能直接作为同槽配方。</p>
{table(['组合','两种单酶产物并集','组合产物','新出现序列','被进一步切碎而消失'],[[names_zh[x['combination']],x['single_enzyme_union'],x['combination_products'],x['new_sequences'],x['lost_sequences']] for x in combos])}
<p>增加酶会同时产生和破坏候选肽，不能按“酶越多、效果越好”选择组合。蛋白酶K在此是切割参照，不意味着它是合适的食品级商业酶制剂。</p>
<h2>3. 排序规则有没有用</h2>
<p>冻结的启发式评分为：C端属于F/W/Y/P加2分，倒数第二位属于A/V/I/L/M/F/Y/W加1分，肽长≤5加1分，总分0–4。没有训练或调参，也没有把已知IC50输入该分数。对照为优先短肽、优先疏水肽和均匀随机选择。评分只是可解释基线，不是活性概率或预测IC50。</p>
<p>ACEip_DB原始1,807条记录中，仅494条满足标准氨基酸、2–20aa、单个精确正数µM或mM值且有DOI。剔除缺DOI658条、单位或范围不明确590条、零或非正值14条，以及长度／字符不符51条。多篇文献共用一条数据库数值的归属和不同实验条件未全部逐篇复核，因此属于异质文献回测，不是统一实验验证集。</p>
{table(['集合','排序方法','Spearman ρ','95%肽级bootstrap区间','Top10%富集倍数'],validation_rows)}
<div class="warn">494条记录上，主要规则ρ={v[0]['methods']['heuristic']['spearman']:.3f}，置换检验p={v[0]['methods']['heuristic']['permutation_p_two_sided']:.3f}；前50条按并列分数折算包含13.24条≤10µM肽，随机期望12.25条，富集仅1.08倍。没有足够证据认定该规则能可靠提高实验命中率。</div>
<p>本轮候选中86条与清洗后的参考记录精确匹配，其中17条被数据库标为IC50≤10µM。这个数不是新实验验证。86条上的规则ρ=0.171、95%区间跨零。GRAVY在该小子集中相关性较高，但Top10%富集仅0.56倍；排序相关和前列命中不是同一指标。</p>
<p>所有集合都由已报道活性肽构成。≤10µM仅是较强活性阈值，其他肽不能当作无活性阴性。这里没有报告二分类AUC，也不能由86条已知肽推断88,637条未测候选的真实精确率。区间采用肽级重采样，未按研究或相似序列聚类，可能低估不确定性；三类集合并非互相独立。</p>
<h2>4. 已发表实验对照</h2>
<h3>鲭鱼蒸煮汁：选酶未通过</h3>
<p>Foods2022，DOI10.3390/foods11121785。按补充表S1取得八种蛋白的同源序列：B5DGQ7、Q5XJ10、O93484、B5DGM7、Q90339、P13104、P49055、P82159。它们来自鲑、斑马鱼、虹鳟、鲤、金鱼与鲻鱼等。P84335替代P13104作为原论文另一个原肌球蛋白匹配的敏感性分析。</p>
<p>实验以1g蛋白、10mg酶粉、4h酶解；胃蛋白酶pH2/37°C，嗜热菌蛋白酶pH8/70°C，蛋白酶K pH8/50°C。相同酶粉质量不等于相同酶活单位。FAPGG法在0.1mg/mL样品条件下，胃蛋白酶组最强，正文86.14%（摘要86.15%），其余五组没有显著差异。</p>
{table(['比较方案','2–20aa产物出现数','规则分≥3出现数','每1000来源残基密度','公开实验'],benchmark_rows)}
<p><b>主分析把蛋白酶K排在最前，胃蛋白酶排在最后，未找回实验赢家。</b>该密度按每种来源蛋白一分子计算，并另检验每蛋白长度归一化。更换原肌球蛋白或允许一次漏切仍未改变蛋白酶K领先的结论。该规则容易偏好释放更多短肽的酶，不能体现混合物中真实活性、丰度、抑制机制和回收损失。</p>
<p>仅覆盖实际6种酶中的3种。木瓜蛋白酶、菠萝蛋白酶和Alcalase缺少本轮可追溯的产品级切割模型，未伪造规则补齐。胃蛋白酶pH2.0规则与论文BIOPEP pH&gt;2模型不完全相同；本轮并非重现原论文BIOPEP结果。单案例三选一的随机正确概率为1/3，不能据此估计跨项目选酶准确率。</p>
<h3>沙丁鱼：六条合成肽活性排序未通过</h3>
<p>Foods2023，DOI10.3390/foods12112216，正文3.1.3节的六个精确IC50值。该研究为dispase酶解3h后加入碱性蛋白酶再3h，HHL底物／HPLC测定。与鲭鱼FAPGG测定不同，不能直接混合数值。按未修饰游离肽分子量换算：µM = mg/mL × 10⁶ ÷ MW。</p>
{table(['肽序列','实验IC50 ± SD（mg/mL）','换算IC50（µM）','本轮规则分','其他鱼类候选库中出现'],seqrows)}
<p>规则ρ=-0.185；样本仅6条且都为经筛选的阳性，结果不足以支持排序有效。六条的活性都在毫摩尔量级，不能因数值约0.6mg/mL就误当作亚微摩尔强抑制剂。原文定性强弱描述与部分数字不完全一致，本报告保留数字原值。由于主库缺少完整沙丁鱼对应蛋白，也没有产品级dispase／碱性蛋白酶规则，未把其他鱼类中能生成这些序列声称为沙丁鱼酶解复现。</p>
<h2>5. 可追溯候选，而非已确认产品</h2>
<p>以下是数据库已有数值的精确序列匹配示例，按数据库IC50排序，供查证与合成可行性评估；它与未测候选的启发式队列分开保存。每行只展示一个来源，完整来源见SQLite。</p>
{table(['序列','数据库IC50 µM','来源蛋白','物种','示例释放方案','残基位置','漏切数'],selected)}
<p>排首的MDLA=0.01µM引用指向一篇综述，而非本轮已核实的原始测定；应先回溯原始实验及单位。VFPS等也仅核验了DOI元数据，未因此升级为本轮实验证据。禁止把这些条目直接当采购或投产清单。未测队列同分时按长度／字典序整理，其先后没有额外科学意义。</p>
<h2>6. 数据异常与范围限制</h2>
<ul><li>Zenodo食品数据集21708604下载与MD5校验成功，审查了24条推导集与16条外部评价摘要；部分曲线明确注明由抑制率反算。保留其作者报告值供后续核查，没有把它当原始仪器数据，也没有将作者Rule5的相关性冒充本流程成绩。这不是对数据真实性的否定。</li>
<li>所有主库序列为数据库完整条目，未统一移除信号肽／前肽，也未模拟胶原羟脯氨酸等修饰、构象可及性、蒸煮或加工。部分候选可能并不存在于原料的成熟蛋白中。</li>
<li>没有原料蛋白质组丰度、批次差异、酶活、时间曲线、LC-MS肽释放强度和回收率，因此不能给出真实肽产率、酶成本收益、最佳温度／pH或工业放大结论。</li>
<li>源蛋白含鱼类过敏原相关家族。肽级毒性、致敏性、苦味、消化稳定性及生物利用度没有验证。ACE体外抑制不等于人体功能或可作食品功能宣称。</li>
<li>本轮成功的是可重复的数据处理与压力测试，失败的是当前活性优先级与选酶规则。未在看见结果后调规则以制造阳性结论。</li></ul>
<h2>7. Mac表现与核验</h2>{resources_table}
<p>核心阶段：酶解与出现位置写入{s['stage_seconds']['digestion_and_occurrence_write_s']:.2f}秒；性质计算与索引{s['stage_seconds']['properties_and_index_s']:.2f}秒；实验回测{s['stage_seconds']['experimental_backtest_s']:.2f}秒；其余为汇总和文件写出。资料获取、方案整理、代码编写、核查与报告另耗人工／交互时间，没有被计入5秒。</p>
<p>6项科学不变量测试通过，覆盖切点例外、N端切割方向、质量守恒、坐标、漏切、组合顺序不变性、分子量及并列排名。独立核验全部1,144,604条出现关系与来源序列一致，无重复出现记录或丢失候选；SQLite完整性通过，冻结输入哈希通过；SciPy独立复核相关系数一致。</p>
<h2>8. 对精准制肽／酶制剂业务的结论</h2>
<p><b>不需要为了这类万级序列筛选升级电脑；需要升级的是实验数据和过程模型。</b>该流程可以进入研发前端，承担原料序列整理、候选生成、来源追踪和实验设计辅助。当前评分不应作为自动淘汰候选或选择最佳酶制剂的依据。</p>
<p>下一步应锁定一种实际原料及3–5种食品级酶制剂，取得批次蛋白质组、具体酶活与pH／温度／时间条件，测量同批酶解物的肽释放谱、ACE活性、肽含量与回收率。用这些数据检验实际释放，再做候选排序。优先比较固定少量方案的实验成本与命中率，暂不以扩大库规模或增加结构模型代替验证。</p>
<p>面试可这样描述：“我搭建并实测了近9万独立海洋蛋白来源候选肽的可追溯筛选流程。16GB Mac计算耗时约5秒。通过真实实验回测发现简单序列规则不能可靠选酶，因此把系统定位为候选生成与数据审查平台，并提出以LC-MS释放数据校准的下一阶段方案。”</p>
<h2>9. 交付与重现</h2><ul>
<li><a href="screen.sqlite">完整SQLite候选／来源／切割位置数据库</a>；<a href="literature-shortlist.json">30条文献标注候选</a>；<a href="exploratory-shortlist.json">30条探索候选</a>。</li>
<li><a href="protocol.json">冻结协议</a>；<a href="summary.json">运行摘要</a>；<a href="validation.json">统计结果</a>；<a href="verification.json">独立核验</a>；<a href="resources/run.log">原始运行日志</a>。</li>
<li><a href="species-sensitivity.json">按物种／代表蛋白的敏感性分析</a>；<a href="process-benchmark.json">选酶回测</a>；<a href="combination-diagnostics.json">组合产生／破坏候选分析</a>；<a href="README.txt">数据库查询与重现说明</a>。</li></ul>
<p class="muted">代码：marine_screen/prepare.py、run.py、audit.py、report.py。输入与协议哈希均保存。重新执行必须使用新的结果目录或明确归档旧结果，脚本拒绝覆盖已有数据库和冻结协议。</p>
<h2>来源</h2><ol>
<li><a href="https://www.uniprot.org/">UniProtKB</a>：完整API查询与每条序列快照存于source目录。</li>
<li><a href="https://doi.org/10.3390/foods11121785">Panjaitan等，Foods2022，11:1785</a>：正文及补充表S1，鲭鱼蒸煮汁酶解实验。</li>
<li><a href="https://doi.org/10.3390/foods12112216">Foods2023，12:2216</a>：沙丁鱼来源合成肽IC50，3.1.3节。</li>
<li><a href="https://github.com/Qian0324/ACEip_DB">ACEip_DB</a>：二次整理文献数据，原始xlsx只读保存；其“学术用途”声明不等于已获得商业再分发授权。</li>
<li><a href="https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html">ExPASy PeptideCutter</a>与<a href="https://github.com/levitsky/pyteomics/blob/master/pyteomics/parser.py">Pyteomics规则快照</a>。</li>
<li><a href="https://doi.org/10.5281/zenodo.21708604">Zenodo21708604</a>：辅助审查，未作为本流程主要实验证据。</li></ol>
</main></html>'''
    (OUT/'report.html').write_text(text)
    readme='''海洋蛋白虚拟酶解与ACE肽筛选实测

先打开report.html。结论：规模可行，当前启发式排序和鲭鱼选酶回测未通过。

screen.sqlite包含：
protein(id,sequence,entries_json)：完整源序列及所有原始UniProt条目映射。
peptide(sequence,length,mw_Da,gravy,charge_pH7,heuristic_score,reference_ic50_uM,reference_doi)。
occurrence(peptide,protein,scheme,start_1based,end_1based_inclusive,missed)：来源位置。位置从1起、两端都包含。
missed=0为完全切割；missed=1仅为恰好一次漏切。MC0–1需查询missed<=1。
NULL参考值表示未纳入参考标注，绝不是无活性。

示例（在当前结果目录运行）：
sqlite3 screen.sqlite "SELECT * FROM peptide WHERE sequence='VFPS';"
sqlite3 screen.sqlite "SELECT * FROM occurrence WHERE peptide='VFPS' LIMIT 20;"
sqlite3 screen.sqlite "SELECT count(DISTINCT peptide) FROM occurrence WHERE missed=0;"

从项目根目录核验：
./run.sh python marine_screen/test_core.py
./run.sh python marine_screen/audit.py

完整重现：先复制marine_screen脚本并将prepare.py中的OUT改为新结果目录；复制冻结source、protocol、proteins、benchmark-proteins输入快照；在新目录执行run.py。
不得覆盖当前已完成结果。网络源可能更新，优先使用保存的快照。

监控命令（实际执行）：
./openfold3/run.sh python tmgh1/monitor.py results/marine-screen/resources /usr/bin/time -l ./run.sh python marine_screen/run.py
仅借用已有环境中的监控依赖，没有运行OpenFold3模型。

完整数据含公开来源资料；对外商业使用或重新发布前按各源条款确认许可。报告不构成食品安全、临床效果、产率或工业工艺认证。
'''
    (OUT/'README.txt').write_text(readme)
    print(OUT/'report.html')

if __name__=='__main__':main()
