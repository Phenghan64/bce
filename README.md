# 本机分子对接环境

本项目使用 Apple Silicon 原生 AutoDock Vina 命令行程序，以及独立的 Python / Meeko 环境，用于蛋白–小分子对接。

另已安装独立的官方 OpenFold3 环境，并通过本机 Apple GPU 泛素预测测试，见 [OpenFold3 使用与实测记录](openfold3/README.md)。

已完成 **OpenFold3 → 结构准备 → AutoDock Vina → 实验结构对照 → 结果选择** 的完整本机验收，见 [闭环流程、实测结果与一键重跑](workflow/README.md)。7L39 小基准完整运行约 52 秒，预测受体和实验受体共六次对接全部通过本次姿势恢复标准。

更复杂的 [ABL 激酶–伊马替尼测试](workflow/1iep-results.md)已完成：293 aa、37 个配体重原子，本机全流程约 4 分 29 秒。程序执行通过，但预测受体三次对接中有两次首位姿势偏差大，不能只凭首位分数判断结果可靠。

另已完成 [101 个分子的 Vina 批量试筛](screening/README.md)：初筛约 2 分 04 秒，前 10 名加对照的三种子复筛约 2 分 43 秒。99 个分子通过初筛结果校验，2 个因未明确的环系立体化学隔离；已验证缓存续跑并保存完整排名。该库无 ABL 活性标签，不代表活性富集验证。

## 使用

在终端进入项目并检查版本：

```bash
cd /Users/yili/Documents/ChatGPT/autodock
./run.sh vina --version
./run.sh python --version
./run.sh mk_prepare_receptor.py --help
```

`run.sh` 自动选择本项目的环境，无需激活 Conda 或配置全局 PATH。Vina 安装的是命令行程序，不包含 Python 的 `import vina` 接口。

重新运行官方 1IEP 示例：

```bash
./example.sh
```

脚本依次准备配体、准备受体及对接盒、执行 Vina 对接、导出 SDF。结果写入 `results/1iep/`；重复运行会覆盖该目录中同名的示例输出。参数固定为 6 个 CPU 线程、`exhaustiveness=32`、随机种子 `20261021`。

- `docking.log`：评分及运行日志。
- `docked.pdbqt`：Vina 输出构象。
- `docked.sdf`：通过 Meeko 恢复键级的分子构象，可供分子可视化工具打开。
- `receptor.pdbqt`、`ligand.pdbqt`、`receptor.box.txt`：准备后的输入。

本次安装验收：完整示例执行成功，最佳 Vina 评分为 −13.261 kcal/mol。日志列出 9 个候选模式；默认 `energy_range=3` 筛选后，PDBQT 和 SDF 各保存 4 个构象，SDF 已通过 RDKit 读取及三维坐标检查。验收记录见 `results/1iep/validation.json`。ProDy 会提示 `pkg_resources` 弃用警告，但当前已固定兼容版本，执行成功。

## 软件与文件位置

- `.tools/bin/vina`：Vina 1.2.7，macOS ARM64 预编译程序。
- `.tools/bin/micromamba`：项目环境管理器。
- `.env/`：Python 3.11、Meeko、RDKit、SciPy、Gemmi 及其依赖。
- `environment.yml`：环境依赖；固定 `setuptools<81` 以支持 ProDy 的 `pkg_resources` 导入。
- `environment-osx-arm64.lock.txt`：本次安装的精确 Conda 包列表，仅适用于 macOS ARM64。
- `downloads.sha256`：本次下载文件的 SHA-256 记录，用于后续完整性比对，不代表上游签名认证。
- `examples/1iep/`：官方 Vina v1.2.7 仓库中的示例输入。

环境使用预编译包，无需系统 Python 或本地编译。本次未修改系统 Python、Xcode 许可及 shell 启动配置。PyMOL 尚未安装；需要交互查看结构时可另外配置。

## 用于自己的分子

通过 `./run.sh mk_prepare_ligand.py` 和 `./run.sh mk_prepare_receptor.py` 准备文件，再用 `./run.sh vina` 对接。应为自己的受体指定合适的对接盒，检查缺失原子、质子化状态、保留的水及辅因子，不能直接套用 1IEP 示例参数。示例通过仅验证工具链可用，不能证明新靶点的预测准确性。

默认 Vina 评分不需要 AutoDock4 或 AutoGrid4；若使用 `--scoring ad4`，则需要另外配置 AutoGrid4 并生成网格。

## 来源

- [Vina 1.2.7 发布页](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7)
- [官方基础对接教程](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)
- [示例输入来源](https://github.com/ccsb-scripps/AutoDock-Vina/tree/v1.2.7/example/basic_docking/data)
- [Meeko 文档](https://meeko.readthedocs.io/)
- [Micromamba 文档](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)

Python 环境中的包来自 conda-forge。示例数据和软件的使用遵循各上游项目许可。

## 一键验收与发布基线

使用项目环境运行统一的工程验收：

```bash
./.env/bin/python3.11 verify.py
```

该命令依次运行后端 `unittest`、前端原生/TypeScript 测试、TypeScript 检查和生产前端构建，并在 `results/verification/<时间>/` 保存日志与 `summary.json`。它不修改生产服务。

需要同时执行一次新的真实 7L39 OpenFold3→Vina 基准时：

```bash
./.env/bin/python3.11 verify.py --scientific
```

科学基准会新建输出目录，检查真实产物、实验结构对照及姿势恢复标准；不会覆盖已有结果。

通用发布工具位于 `deploy/marine-agent/publish_backend.py` 和 `publish_frontend.py`。两者都支持 `--dry-run`，先验证指纹、验收证据和当前版本，再决定是否需要正式发布。发布参数和私有配置不写入仓库。
