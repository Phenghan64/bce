# OpenFold3 → AutoDock Vina 完整流程验收

较复杂的 ABL 激酶–伊马替尼（293 aa、37 个配体重原子）使用 `results/closed-loop-1iep/` 中的 `case.json`、`query.json`、`reference/` 与 `msa/`。复用这些输入重跑时使用新输出目录：

该案例本机全流程 **4 分 29 秒**，预测受体出现首位姿势排序歧义，详见[复杂案例实测报告](1iep-results.md)。

```bash
./run.sh python workflow/run.py --case-dir results/closed-loop-1iep \
  --output-dir results/closed-loop-1iep-repeat
```

`--single-sequence` 可显式取消 MSA，作为单独的消融测试。不能将无 MSA 结果与带 MSA 结果混记；受体准备报错时流程停止并保留失败日志。`case.json` 固定配体的电荷状态和对接盒大小。蛋白实验编号通过 mmCIF 的 `label_seq_id` 对应完整输入序列；不是直接比较不同编号体系的原子。

本目录提供已在本机实际运行的 **7L39 / T4 lysozyme L99A–甲苯**基准测试。它覆盖结构预测、蛋白及配体准备、对接、实验结构对照、姿势验证和按分数选择结果，是计算流程闭环，不包含实验验证或模型再训练。本次使用官方 OpenFold3 的 OpenBind 权重和 AutoDock Vina 1.2.7；不是 DeepMind AlphaFold 3、OpenFold3-MLX 或 AutoDock4。

## 一条命令重跑

在项目根目录执行，选择一个尚未运行过的新目录：

```bash
./run.sh python workflow/run.py --output-dir results/closed-loop-7l39-new
```

脚本固定用于这个已知基准，不是可直接套用所有新靶点的通用筛选程序。它依赖本项目已安装的三个独立环境、OpenBind 权重，以及 `results/openfold3-7l39/` 中的缓存 MSA 和实验参考。重跑无需再次访问 MSA 服务。默认会重新执行 OpenFold3；`--reuse-prediction` 仅用于使用指定结果目录内已完成的预测，不会去别处复制预测。已有 `manifest.json` 时拒绝覆盖，请换新目录。

## 实际执行的步骤

1. 记录输入、代码和模型权重的 SHA-256。
2. 输入 164 aa 蛋白序列和甲苯 SMILES，使用缓存 MSA，在 Apple MPS 上运行新的 OpenFold3 预测：seed 42、1 个扩散样本、无模板、float32、low_mem。
3. 导出预测蛋白/配体和实验参考。将 162 个实验可见 Cα 与预测对应残基对齐；实验末端缺失的两个残基不补建。
4. 蛋白使用 PDBFixer/OpenMM 在标称 pH 7.4 下加氢、补末端 OXT，检查已有重原子不移动。本例无须补其他重原子。实验对照删除水和结晶添加物，保持与预测受体一致的干燥受体策略。
5. **仅从甲苯 SMILES** 使用 ETKDGv3（seed 20261021）和 MMFF 生成独立初始构象。Meeko 准备配体及受体 PDBQT，不以已知姿势初始化配体。
6. 预测受体使用 OpenFold3 配体质心定义 20 × 20 × 20 Å 对接盒；实验受体使用实验配体质心定义同样大小的盒。因此这是预测引导的局部对接，不是盲对接或独立口袋发现。
7. 两种受体各运行 Vina 三次：seeds 42 / 2026 / 20261021、6 CPU 线程、exhaustiveness 32、num_modes 9、energy_range 4、刚性受体、vina 评分函数。受 energy_range 过滤影响，实际每次输出 8–9 个姿势。
8. Meeko 导出 SDF；检查全部 51 个姿势的化学图、7 个重原子、有限坐标、盒内坐标、PDBQT/SDF/评分记录数量一致。报告全部姿势，不只报告最接近实验的那个。
9. 用固定蛋白对齐和考虑分子对称性的 `CalcRMS` 计算配体 RMSD，**不对配体独立旋转或平移拟合**。评价首位姿势、最佳 RMSD 姿势排名、种子间一致性、近距离接触和口袋残基接触。
10. 按预测受体上的最低 Vina 分数选择交付姿势；实验 RMSD 只作验收，不参与候选选择。把选中姿势、OpenFold3 和实验配体导出到同一实验坐标系，形成对比反馈。

## 2026-10-21 验收结果

正式通过目录：[`results/closed-loop-7l39-validated`](../results/closed-loop-7l39-validated/)。所有 21 个记录阶段退出码为 0，六次对接首位姿势均满足本次设定的重原子 RMSD ≤ 2 Å 基准标准。

| 指标 | 实测 |
|---|---:|
| 整条计算流程墙钟时间 | 51.86 s |
| OpenFold3 模型内部运行 | 29.59 s |
| OpenFold3 CLI 导入和执行 | 36.32 s |
| OpenFold3 进程峰值 RSS | 3.93 GiB |
| 采样到的峰值 MPS driver 内存 | 2.49 GiB |
| 预测期间整机 swap 增量 | 62.31 MiB |
| OpenFold3 平均 pLDDT / ipTM | 97.95 / 0.934 |
| 蛋白 162 Cα RMSD | 0.350 Å |
| OpenFold3 配体 RMSD | 0.391 Å |

RSS 与 MPS 内存有重叠，不能相加。swap 是整机指标，包含其他应用。总耗时不包含软件安装、模型下载和之前的远程 MSA 获取。本例规模很小，不能据此外推大型复合物的内存需求。

| 受体 | 种子 | 首位 Vina 分数（kcal/mol） | 首位配体 RMSD（Å） |
|---|---:|---:|---:|
| OpenFold3 预测 | 42 | −6.176 | 0.578 |
| OpenFold3 预测 | 2026 | −6.173 | 0.522 |
| OpenFold3 预测 | 20261021 | −6.170 | 0.572 |
| 实验结构对照 | 42 | −5.968 | 1.585 |
| 实验结构对照 | 2026 | −5.975 | 1.588 |
| 实验结构对照 | 20261021 | −5.958 | 0.622 |

预测受体的三个首位姿势互相 RMSD 为 0.075–0.103 Å，结果一致；首位姿势均无小于 2 Å 的蛋白–配体重原子对。实验对照首位姿势出现两个方向，其相互 RMSD 最高 1.48 Å，说明分数排序仍有歧义。≤2 Å 是本例声明的验收阈值，不意味着姿势完全相同。

**本例 Vina 没有比 OpenFold3 原始配体更接近实验。** OpenFold3 为 0.391 Å，预测受体对接首位为 0.522–0.578 Å。Vina 在这里提供另一种采样、评分与复核，不能默认视作精修。不同受体分数也不能直接证明哪种结构具有更强真实亲和力。

这是已知小基准，不是未知靶点的泛化准确性证明，也不检验蛋白柔性、水、辅因子或质子化状态的多种可能性。

## 交付文件与复现记录

- [`report.txt`](../results/closed-loop-7l39-validated/report.txt)：简明结果。
- [`report.json`](../results/closed-loop-7l39-validated/report.json)：51 个姿势的完整指标、接触、种子一致性和选择规则。
- [`manifest.json`](../results/closed-loop-7l39-validated/manifest.json)：命令、阶段耗时、退出码和输入/代码/权重 SHA-256。
- [`selected-docked-toluene.sdf`](../results/closed-loop-7l39-validated/selected-docked-toluene.sdf)：按预测受体 Vina 分数选中的姿势，位于预测受体原始坐标系，应与该目录下 `predicted/receptor-h.pdb` 一起使用。
- [`comparison-in-crystal-frame.sdf`](../results/closed-loop-7l39-validated/comparison-in-crystal-frame.sdf)：实验、OpenFold3、选中 Vina 配体三个构象，已统一到实验坐标系，应与 `reference-receptor.pdb` 一起显示。
- `predicted/`、`reference/`：加氢结构、Meeko 参数、对接盒、六次 PDBQT 和 SDF。
- `logs/`：每一步日志；`openfold/benchmark.json`：推理内存采样。
- [`software-versions.json`](software-versions.json)：实际工具版本；[`prep-osx-arm64.lock.txt`](prep-osx-arm64.lock.txt)：新增加氢环境的精确 Conda 包。

首次整合在 `results/closed-loop-7l39/` 保留了失败记录：复用 Biotite PDB 对象写第二个结构导致结构追加。修正为独立 PDB 对象，并增加重复原子检查后，上述正式目录从 OpenFold3 起完整重跑通过；失败目录不作为成功结果。

需要重建本项目独立加氢环境时：

```bash
.tools/bin/micromamba create -y --no-rc -r "$PWD/.tools/mamba" \
  -p "$PWD/.tools/prep-env" -f workflow/prep-osx-arm64.lock.txt
```
