鳕鱼骨架：真实工艺参照的虚拟酶解与功能肽筛选评估
================================================

已完成结果：../results/enzyme-case/report.html
这是回顾性公开实验复现和计算能力评估，不是企业生产验证。
本轮不调用 OpenFold 或 AutoDock/Vina；肽组学、虚拟切割与功能分数
的评价终点与小分子对接不同，没有为凑流程加入不适用的软件。

复现（在当前机器与保存的环境中）
  cd /Users/yili/Documents/ChatGPT/autodock
  bash enzyme_case/reproduce.sh
脚本创建带时间的新结果目录，不覆盖已交付结果。
该脚本沿用本次机器信息与已核验模型锚点文件，适用于本机原环境复现。

核心命令与环境
  .tools/enzyme-case/venv/bin/python enzyme_case/models.py
  .tools/enzyme-case/venv/bin/python enzyme_case/run.py
  .tools/enzyme-case/venv/bin/python enzyme_case/audit.py
  .tools/enzyme-case/venv/bin/python enzyme_case/report.py
上述直接命令默认写入原结果目录；保存原结果时请用 reproduce.sh。
Python 3.11、TensorFlow 2.15.1，依赖精确版本见 environment.txt。
模型架构与权重没有重新训练；检查点六个张量直接加载。

数据重提取（需 openpyxl，可使用 Codex 自带 Python）
  python enzyme_case/extract_source.py --check
重读 Table S3.xlsx，逐值、逐原始行核对 JSON。没有改动原始工作簿。
去掉 --check 将提取结果写至标准输出，不覆盖文件。

原始输入和下载地址
  研究：10.3390/md18120599 / PMC7759894
  附件：https://mdpi-res.com/d_attachment/marinedrugs/marinedrugs-18-00599/article_deploy/marinedrugs-18-00599-s001.zip
  肽组学数据：10.17632/pc5h9drk6d.2
  蛋白：https://ftp.ensembl.org/pub/release-90/fasta/gadus_morhua/pep/Gadus_morhua.gadMor1.pep.all.fa.gz
  AnOxPePred：https://github.com/TobiasHeOl/AnOxPePred
  EmulsiPred：https://github.com/MarcatiliLab/EmulsiPred
下载归档留存在 results/enzyme-case/source/，解压的模型在
.tools/enzyme-case/；归档 SHA256、源文件 SHA256、代码和权重
校验信息见结果 manifest.json。未重跑原始 .raw 质谱鉴定。

方案与主要限制
  protocol.json 是主计算前冻结的回顾性方案；不是预注册或盲测。
  secondary-audit.json 是看到主结果后追加的长度基线和训练重叠审计。
  不能用其选择新参数，也没有重拟合/改变本次筛选结果。
  六条件按 Alc-MCF, Neut-MCF, Neut&Alc-MCF,
  Alc-HCM, Neut-HCM, Neut&Alc-HCM 排列。
  MCF=绞碎骨架；HCM=加热去骨鱼肉。二者有处理与成分双重差别。
  单酶3h、串联6h，不是等工时对照；阶段投酶量未擅自推算。
  释放是该研究终态的相邻残基切割替代模型，不是酶学动力学。
  未检出不等于不释放，质谱强度不等于产率。
  测试蛋白在该研究中已被检出，故测试群体存在检出选择偏差。
  共享片段分组与精确肽排除降低泄漏，但未完全消除远缘同源。
  单肽模型输出不能直接当水解物实验活性；本次该聚合验证失败。
  相同研究的六组配方不能替代独立原料批次验证。

结果格式
  screen.sqlite：peptide 为公开模型分数；observed 为每条件实测
  强度、鉴定方式、PEP和源行号；virtual 为入选虚拟肽来源与释放排名。
  SQL 示例：
    SELECT sequence,FRS FROM peptide ORDER BY FRS DESC LIMIT 20;
  这只按模型分数排序，不等于已验证功能肽。
  observed-peptide-scores.json：11968条实测肽，含完整证据。
  virtual-candidates.jsonl.gz：4502080条去重虚拟肽，含所有源位置，
  六条件最大来源释放分数及简单基线。起止为1-based、闭区间。
  一条肽可能多个来源；分数按各条件来源最大值聚合。
  virtual-shortlist-scores.jsonl.gz：35691条虚拟入选肽及功能分数。
  experiment-candidates.json：各目标20条高置信检出、>=10aa候选。
  assay/bulk/release/model/audit JSON：全部统计与限制。
  资源：compute/ 与 audit-compute/ 下的资源逐秒采样与日志。
  图：SVG可编辑导出，PNG可直接放入面试材料。

结论：本机吞吐可行，酶解释放排序局部有效；当前功能配方排序
未通过验证。没有做新湿实验，没有证明产率、活性或工业可放大。
