# OpenFold3 本机测试

2026-10-21 在 Apple M6 / 16 GiB / macOS 27.0 上完成官方 OpenFold3 的两次泛素预测。实际使用 Apple GPU 的 MPS 后端，无 CUDA，无模板输入。

## 安装来源

- 官方代码：[aqlaboratory/openfold-3](https://github.com/aqlaboratory/openfold-3)，2026-10-21 下载的 main 源码快照，保存在项目 `.tools/openfold3/`。
- 元数据版本为 0.5.0，这是源码归档中的 setuptools-scm fallback 版本，不是已验证的 Git release 标签。归档 SHA-256 记录在本目录 `downloads.sha256`。
- 使用当时官方默认的 **OpenBind 权重** `of3-ob-2025-06-30-174k.pt`，位于 `.tools/openfold3-cache/`。这不是 DeepMind 原版 AlphaFold3 权重，也不是第三方 OpenFold3-MLX。
- 环境位于 `.tools/openfold3-env/`，Python 3.14.6 / PyTorch 2.12.1；Metal 矩阵运算及真实模型推理均已通过。
- 由于 Pixi 下载器连接失败，导出了上游 `pixi.lock` 中 `openfold3-base/osx-arm64` 的精确包列表，使用 micromamba 安装，再安装官方源码。记录在 `environment-osx-arm64.lock.txt`、`python-packages.lock.txt`。
- 没有修改已有 Vina 环境或接受 Xcode 许可。

## 实测结果

输入为官方泛素示例，76 个残基；1 个随机种子（42）、1 个扩散样本、float32、low_mem 模式、0 个数据加载工作进程。MPS 分配预算设为其推荐上限的 70%，本机约 8.29 GiB。

| 指标 | 无 MSA 首次运行 | 带缓存 MSA 第二次运行 |
|---|---:|---:|
| OpenFold3 模型阶段耗时 | 26.58 s | 11.37 s |
| CLI 加载与执行耗时 | 97.26 s | 18.65 s |
| 进程峰值 RSS | 4.03 GiB | 4.07 GiB |
| 采样峰值 MPS 驱动分配 | 2.47 GiB | 2.57 GiB |
| 全机交换空间增量 | 268.19 MiB | 0 MiB |
| 平均 pLDDT | 85.31 | 90.28 |
| pTM | 0.842 | 0.890 |
| 对 1UBQ 的 Cα RMSD | 2.40 Å | 1.50 Å |

CLI 时间不包含 Python/PyTorch 启动、软件安装、权重下载和 MSA 获取。第二次运行受益于预热缓存，不能把速度差异解释成 MSA 加速。RSS 与 MPS 显存指标在统一内存架构下有重叠，不能相加；交换空间属于全机所有应用。

参考结构来自上游测试数据 `openfold3/tests/test_data/mmcifs/1ubq.cif`。RMSD 对全部 76 个匹配的 Cα 原子做 Kabsch 刚体对齐计算，没有去除末端或离群点。输出均能解析、坐标有限且残基序列匹配；模型自身未报告碰撞。泛素是已知小蛋白，结果不能证明对未见靶点的泛化准确性或大复合物的内存容量。

## 运行

```bash
cd /Users/yili/Documents/ChatGPT/autodock
./openfold3/test.sh
```

自动复用缓存 MSA，并创建新的 `results/openfold3-ubiquitin/repeat-时间戳/` 目录。可传入一个尚不存在的输出路径：

```bash
./openfold3/test.sh results/my-ubiquitin-test
```

运行任意 OpenFold3 命令：

```bash
./openfold3/run.sh run_openfold predict --help
```

`runner.yml` 显式指定 MPS。`benchmark.py` 在 MPS 不可用时会报错，避免误把 CPU 运行当作 GPU 测试。

## 明确 accession 的通用任务入口

当用户已经提供 UniProt accession 时，可以运行完整的确定性链路：UniProt 序列获取与校验、MSA、Apple MPS 上的 OpenFold3 推理、置信度闸门和前端结构 manifest。

```bash
./.env/bin/python openfold3/predict_accession.py \
  --accession P0CG47 \
  --output-dir results/openfold3-jobs/P0CG47
```

也支持一个只含单条记录的 FASTA：

```bash
./.env/bin/python openfold3/predict_accession.py \
  --fasta input/protein.faa \
  --output-dir results/openfold3-jobs/protein
```

任务输出包括 `input.json`、`query.json`、`msa/query_msa.json`、OpenFold3 原始预测、`pipeline.json`，以及 `viewer/manifest.json`、mmCIF、置信度 JSON。默认使用一个固定 seed、一个 diffusion sample、关闭运行时模板输入；相同输出目录再次运行会复用已完成任务，`--force` 才会重新计算。

该入口故意不接受模糊蛋白名称。名称解析、物种确认和候选选择必须在上游完成后，再把明确 accession 交给这里，避免把同名蛋白、片段或错误物种静默送进结构预测。

### 缓存与多配体输入

UniProt 响应会写入 `results/openfold3-uniprot-cache/`，默认对 429、5xx 和网络超时做有限指数退避。`pipeline.json` 记录来源 URL、序列哈希、缓存命中和获取时间；`--refresh-uniprot` 才会绕过缓存。

MSA 按蛋白序列 SHA-256 与 MSA 参数写入 `results/openfold3-msa-cache/`。任务的 `query_msa.json` 会被复制到当前任务并重写绝对路径，再把当前任务的配体链绑定回蛋白 MSA。因此 `--force` 只重跑预测，`--refresh-msa` 才重新请求公共 ColabFold 服务。

本机 MPS 推理使用 `results/openfold3-runtime/openfold3-mps.lock` 做跨进程排队。CLI、Deep Agent 和 HTTP API 同时提交时，MSA 可以并行复用，但 OpenFold3 模型阶段会串行执行；`pipeline.json` 会记录排队等待时间及系统交换空间增量。交换空间增长超过 `--max-swap-growth-gib`（默认 1 GiB）会标记资源警告，提示降低并发或迁移到更大执行节点。

可通过 `--ligands-json` 传入 JSON 文件或内联 JSON 数组。每个配体必须使用唯一的单字符链 ID，并且只能提供一种化学来源：`smiles`、`ccd_codes` 或 `sdf_file_path`。默认输出目录包含配体集合指纹，避免相同 accession 的不同配体任务相互覆盖；所有预测仍应被解释为结构审查输入，不能把配体置信度或模型分数当作实验证据。

例如：

```bash
./.env/bin/python openfold3/predict_accession.py \
  --accession P0CG47 \
  --ligands-json '[{"chain_id":"X","smiles":"Cc1ccccc1"},{"chain_id":"Y","smiles":"Oc1ccccc1"}]'
```

### MSA 网络恢复

仅对官方 `https://api.colabfold.com` 请求使用有限重试。Mac 直连失败时，可以通过 `~/.config/marine-agent/msa-network.json` 中的 `fallback_proxy` 指定现有 SSH 的回环 SOCKS5h 出口（例如 `socks5h://127.0.0.1:18199`）；环境变量 `OPENFOLD_MSA_PROXY` 可覆盖此设置。该设置不修改系统代理，TLS 校验保持开启，结构推理仍在本机。当前 SSH 出口由 `chat.marine.msa-egress` LaunchAgent 维护。

每个网络请求最多尝试三次，整个客户端网络等待预算为十分钟，外层 MSA 子进程另有十五分钟超时。提交响应读取超时或网关错误可能意味着服务端已经受理，因此不会自动重复该 POST。失败时保留 `msa/msa.log`、`failure.json`、`network.json`（若已生成）；API 对话仅显示简洁原因，原始诊断留在任务记录。没有关闭 MSA 或回退到单序列预测。

UniProt `proteinDescription.flag` 的片段标记同时进入候选列表、确认页、输入记录和结果摘要；获取到条目序列不再表述为“完整蛋白序列”。

## MSA 和结果

结构推理全部在本机进行；带 MSA 测试的序列比对由公共 ColabFold 服务生成，所以首次完整流程需要联网。默认下载超时触发重试，`fetch_msa.py` 将该服务的读取超时延长为 120 秒后下载成功；未改变模型计算。已有缓存后可复用，无需再次提交序列。

- 综合验收：`../results/openfold3-ubiquitin/validation.json`
- 结果比较：`../results/openfold3-ubiquitin/comparison.csv`
- 带 MSA 的结构：`../results/openfold3-ubiquitin/with-msa/prediction/ubiquitin/seed_42/ubiquitin_seed_42_sample_1_model.cif`
- 运行日志与逐秒资源采样：各次运行目录的 `run.log` 和 `benchmark.json`
- 重新验证这两次验收结果：`./openfold3/run.sh python openfold3/validate.py`（在项目根目录运行）

泛素测试只验证单体蛋白结构预测，没有把泛素当作小分子对接靶点。

## 第二个案例：T4 溶菌酶与甲苯

随后使用上游 `query_single_protein_single_ligand.json` 示例，完成了 **164 残基的 T4 溶菌酶 L99A 与甲苯**的联合结构预测。输入保存在 `t4-lysozyme-toluene.json`，参考实验结构为 [PDB 7L39](https://www.rcsb.org/structure/7L39)。采用缓存 MSA、无模板、无口袋约束、种子 42、1 个扩散样本，其余参数与泛素测试一致。

| 指标 | 实测结果 |
|---|---:|
| 模型阶段耗时 | 30.92 s |
| CLI 加载与执行耗时 | 38.64 s |
| 进程峰值 RSS | 3.98 GiB |
| 采样峰值 MPS 驱动分配 | 2.49 GiB |
| 全机交换空间增量 | 69.81 MiB |
| 平均 pLDDT | 97.95 |
| 蛋白–配体 ipTM | 0.934 |
| 蛋白 Cα RMSD | 0.350 Å |
| 配体重原子 RMSD | 0.391 Å |
| 配体中心位置误差 | 0.319 Å |

蛋白比较使用实验结构中可见的全部 162 个 Cα 原子；末端两个残基没有实验坐标。配体比较先按蛋白做刚体对齐，再对甲苯 7 个重原子计算化学图匹配和对称性校正后的 RMSD，没有单独拟合配体位置。参考配体键级来自 RCSB 的 MBN 化学组分定义。此结果证明该小型复合物可在本机运行，不能推广为未见靶点准确率或大型复合物容量结论。

完整记录位于 `../results/openfold3-7l39/validation.json`，原始 CIF 与置信度文件位于 `../results/openfold3-7l39/with-msa/prediction/pdb_7L39/seed_42/`。另已导出同一预测坐标系下的 `predicted-receptor.pdb` 和 `predicted-toluene.sdf`，供后续接入 Meeko/Vina；尚未做对接前质子化准备或运行这一案例的 Vina 对接。

重新验证结果（在项目根目录）：

```bash
./openfold3/run.sh python openfold3/validate_7l39.py
```

### 导出到 DOHO 结构查看器

预测完成后，可将原始 mmCIF、置信度和验证指标导出为前端统一的 manifest：

```bash
./.env/bin/python openfold3/export_viewer_artifact.py \
  --prediction-dir results/openfold3-7l39/with-msa/prediction/pdb_7L39/seed_42 \
  --output-dir /Users/yili/Documents/dsh/doho/public/structure-preview/live-7l39 \
  --artifact-id 7l39-openfold-live \
  --title 'T4 溶菌酶 L99A' \
  --validation-json results/openfold3-7l39/validation.json \
  --model-version openbind-2025-06-30-174k \
  --seed 42 --no-templates
```

本地查看时使用：

```text
http://localhost:4178/doho/?structure-preview=1&structure-manifest=/doho/structure-preview/live-7l39/manifest.json
```

导出器会复制原始 mmCIF、写入 `validation.json` 和 `manifest.json`，并记录结构文件 SHA-256。前端加载前会再次检查坐标和置信度记录；任务失败或产物不完整时不会进入就绪状态。

本次推理命令（重新测试时请换一个输出目录和 benchmark 路径，避免覆盖已有验收记录）：

```bash
./openfold3/run.sh python -u openfold3/benchmark.py \
  results/openfold3-7l39/with-msa/benchmark.json predict \
  --query-json results/openfold3-7l39/msa/query_msa.json \
  --runner-yaml openfold3/runner.yml --num-diffusion-samples 1 \
  --use-msa-server false --use-templates false \
  --output-dir results/openfold3-7l39/with-msa/prediction
```
