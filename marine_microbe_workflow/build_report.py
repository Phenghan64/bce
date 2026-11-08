#!/usr/bin/env python3
"""Build the reproducible marine workflow report from recorded run artifacts."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/marine-microbe-workflow"
RUNS = RESULTS / "runs"


def read_json(path: str):
    return json.loads((ROOT / path).read_text())


def gib(value: int | float) -> float:
    return round(value / (1024**3), 2)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    audit = read_json("results/marine-microbe-workflow/runs/candidate-audit-v2/candidate-audit.json")
    assembly = read_json("results/marine-microbe-workflow/runs/spades-full2/assembly-audit.json")
    ipr = read_json("results/marine-microbe-workflow/runs/interpro-candidate-panel/interpro-panel.json")
    recovery = read_json("results/marine-microbe-workflow/runs/recovery-tests/recovery-tests.json")
    of_batch = read_json("results/marine-microbe-workflow/runs/openfold3-candidates/benchmark-batch.json")
    of_msa = read_json("results/marine-microbe-workflow/runs/openfold3-candidates/benchmark-msa-EOHHJP_01311.json")
    panel_meta = read_json("results/marine-microbe-workflow/runs/interpro-candidate-panel/panel.json")
    panel_by_id = {x["query"]: x for x in panel_meta["records"]}
    of_conf = {}
    for p in (RUNS / "openfold3-candidates/prediction-batch").glob("*/seed_*/*confidences_aggregated.json"):
        of_conf[p.parent.parent.name] = json.loads(p.read_text())
    msa_conf_path = next((RUNS / "openfold3-candidates/prediction-msa-EOHHJP_01311").glob("*/seed_*/*confidences_aggregated.json"))
    msa_conf = json.loads(msa_conf_path.read_text())

    sources = {}
    for rel in [
        "source/SAMN47740168.xml", "source/SRR32944547.xml", "source/ena-filereport.tsv",
        "source/SRR32944547_1.fastq.gz", "source/SRR32944547_2.fastq.gz",
        "source/reference/reference.fna", "source/reference/reference.faa", "source/reference/reference.gff",
    ]:
        p = RESULTS / rel
        sources[rel] = {"bytes": p.stat().st_size, "sha256": sha256(p)}
    sources["source/SRR32944547_1.fastq.gz"]["ena_md5"] = "7449be587413dc3ee634521064759001"
    sources["source/SRR32944547_2.fastq.gz"]["ena_md5"] = "3b67d6f1685535d45147f988f67324f7"

    fastp = {
        "version": "1.4.0", "input_pairs": 6733603, "clean_pairs": 6514384,
        "retained_fraction": round(6514384 / 6733603, 4), "runtime_s": 17.92,
        "peak_rss_gib": 1.34, "swap_delta_bytes": 0,
        "outputs": ["fastp.json", "fastp.html", "SRR32944547.clean.R1.fastq.gz", "SRR32944547.clean.R2.fastq.gz"],
    }
    spades = {
        "version": "4.3.0", "runtime_s": 597.96, "peak_rss_gib": 5.25,
        "swap_delta_bytes": 0, "contigs": assembly["assembly_contigs"], "assembly_bp": assembly["assembly_bp"],
        "max_contig_bp": assembly["assembly_max_bp"], "n50_bp": assembly["n50"]["n"],
        "reference_coverage_at_95_identity": assembly["reference_coverage_at_95_identity"],
        "warning": "macOS setrlimit could not enforce the requested 12 Gb SPAdes limit; use Linux cgroups/HPC scheduler in production",
    }
    bakta = {
        "version": "1.12.2", "database": "light 6.0 (validation only)", "runtime_s": 96.02,
        "peak_rss_gib": 2.80, "cds": 3951, "trna": 99, "tmrna": 1, "rrna": 8, "ncrna": 8,
        "hypothetical_proteins": 421,
        "warning": "full production DB and fixed AMRFinderPlus/dependency closure are still required",
    }
    blast = {
        "version": "2.17.0+", "query_proteins": 3951, "database": "SSi026 reference proteins",
        "cutoffs": {"identity_percent": 80, "query_coverage_percent": 80, "evalue": 1e-20},
        "output": "runs/blastp-assembly-vs-reference.tsv",
    }
    interpro = {
        "workflow": "InterProScan 6.0.3.0",
        "api": "InterPro Matches API 0.6.0, release 110.0 (2026-09-03)",
        "panel_sequences": ipr["input"]["sequence_count"], "panel_found": ipr["run"]["found_count"],
        "panel_matches": ipr["run"]["total_match_count"], "runtime_s": round(ipr["run"]["elapsed_seconds"], 3),
        "repeat_runtime_s": round(ipr["repeat"]["elapsed_seconds"], 3),
        "repeat_identical": ipr["repeat"]["byte_level_normalized_match"],
        "scale_validation": {"input_sequences": 1000, "unique_sequences": 981, "matches": 64751, "first_runtime_s": 24.3, "repeat_runtime_s": 14.5, "repeat_unique_sequences": 981},
        "local_workflow_smoke": {"input": "A0A194WDG1", "workflow": "InterProScan 6 Nextflow, MobiDB-lite + COILS path", "runtime_s": 8.6, "outputs": ["JSON", "JSONL", "TSV"], "full_local_blocker": "idrpred command missing when Matches API disabled; another local path lacked groovy"},
        "panel": [{"id": x["id"], "role": panel_by_id[x["id"]]["panel_role"], "length": x["length"], "pfam": x["pfam"], "ncbifam": x["ncbifam"], "interpro": x["interpro"], "match_count": x["match_count"]} for x in ipr["run"]["summaries"]],
    }
    openfold = {
        "version": "0.5.0", "torch": "2.12.1", "device": "Apple MPS",
        "no_msa_batch": {"queries": 4, "successful": 4, "runtime_s": of_batch["elapsed_seconds"], "peak_rss_gib": gib(of_batch["process_peak_rss_bytes"]), "peak_mps_driver_gib": gib(of_batch["sampled_peak_mps_driver_bytes"]), "swap_delta_gib": gib(of_batch["system_swap_after_bytes"] - of_batch["system_swap_before_bytes"]), "confidence": {k: {"avg_plddt": round(v["avg_plddt"], 2), "ptm": round(v["ptm"], 3), "gpde": round(v["gpde"], 3)} for k, v in of_conf.items()}},
        "msa_control": {"id": "EOHHJP_01311", "msa_fetch_s": 153, "runtime_s": of_msa["elapsed_seconds"], "peak_rss_gib": gib(of_msa["process_peak_rss_bytes"]), "peak_mps_driver_gib": gib(of_msa["sampled_peak_mps_driver_bytes"]), "swap_delta_gib": gib(of_msa["system_swap_after_bytes"] - of_msa["system_swap_before_bytes"]), "avg_plddt": round(msa_conf["avg_plddt"], 2), "ptm": round(msa_conf["ptm"], 3), "gpde": round(msa_conf["gpde"], 3)},
        "interpretation": "no-MSA structures were emitted but low confidence; MSA materially improved the same 458-aa candidate and is required before catalytic-site or docking decisions",
    }
    report = {
        "title": "真实海洋微生物候选蛋白酶发现工作流验证报告",
        "date": "2026-10-23",
        "workflow": ["public ENA FASTQ", "fastp", "SPAdes", "Bakta", "BLAST+", "InterProScan/Matches API", "OpenFold3"],
        "sample": {"organism": "Alteromonas macleodii", "strain": "SSi026", "biosample": "SAMN47740168", "ena_run": "SRR32944547", "isolation": "Sargassum ilicifolium, St John's Island wharf, Singapore", "collection_date": "2022-08-02", "reference": "NCBI GCF_049379645.1 / ASM4937964v1"},
        "source_files": sources,
        "stages": {"fastp": fastp, "spades": spades, "bakta": bakta, "blast": blast, "interpro": interpro, "openfold3": openfold},
        "candidate_audit": audit,
        "recovery_tests": recovery,
        "scientific_limits": [
            "The held-out reference labels are annotation/homology evidence, not biochemical activity measurements.",
            "The BLAST database in this benchmark is the same-strain reference protein set; it validates assembly/annotation consistency and is not an independent discovery benchmark. Production discovery must use versioned external MEROPS/UniProt/NR databases and a blind holdout.",
            "The sample is a cultured marine isolate, not a mixed metagenome; metagenomic projects require metaSPAdes/MEGAHIT, binning and contamination checks.",
            "InterPro Matches API results are precomputed evidence; complete member scans require Linux containers/HPC and full licensed databases.",
            "OpenFold3 confidence supports fold/domain and structure review only; it does not prove catalytic activity, substrate specificity, secretion, thermostability or industrial performance.",
            "A protease score or structural confidence must never be converted directly into enzyme activity or host salt tolerance.",
        ],
        "engineering_conclusion": {"fastp": "retain with paired-read and log/output QC", "SPAdes": "retain; Linux cgroups/HPC for production", "Bakta": "retain; full DB and pinned dependency image", "BLAST+": "retain; versioned databases and taxonomy", "InterProScan6": "retain API-first plus Linux full scan fallback", "OpenFold3": "retain as gated structural evidence; GPU/Linux for production", "overall": "scientifically credible as a candidate-generation and evidence-integration pipeline, not an autonomous enzyme-activity predictor"},
        "deployment": {"mac_mini": "integration, small panels, cached API lookups, and low-throughput structure triage", "linux_hpc": "full Bakta/InterProScan databases, cgroup resource limits, GPU OpenFold3, Slurm/Nextflow", "cloud": "object-store inputs, container digests, API/cache layer, batch queues and immutable provenance", "gates": ["input schema/QC", "database checksum and release", "per-stage output count", "log semantic error scan", "retry/idempotency", "candidate evidence threshold", "human/experimental confirmation"]},
    }
    (RESULTS / "REPORT.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")

    audit_cm = audit["confusion_matrix"]
    panel_lines = []
    for x in interpro["panel"]:
        panel_lines.append(f"| {x['id']} | {x['role']} | {x['length']} | {', '.join(x['pfam']) or '-'} | {', '.join(x['ncbifam']) or '-'} | {x['match_count']} |")
    recovery_lines = []
    for x in recovery["cases"]:
        if x["name"] == "resume_from_completed_fastp_spades":
            recovery_lines.append(f"| {x['name']} | continuation | {'PASS' if x['passed'] else 'FAIL'} |")
        else:
            state = "semantic error detected" if x.get("semantic_failure_detected") else ("non-zero exit" if x.get("exit_code") else "exit 0")
            recovery_lines.append(f"| {x['name']} | exit {x.get('exit_code')} | {state} |")
    md = f"""# 真实海洋微生物候选蛋白酶发现工作流验证报告

日期：2026-10-23。目标链路为：**公开原始 reads → fastp → SPAdes → Bakta → BLAST+ / InterProScan 6 → OpenFold3**。本报告记录真实运行结果、批量和异常恢复，不把注释、结构置信度或对接分数等同于真实酶活。

## 数据与可追溯性

样本是海洋分离株 *Alteromonas macleodii* SSi026，来自新加坡 St John's Island wharf 的海藻 *Sargassum ilicifolium*，BioSample `SAMN47740168`，ENA run `SRR32944547`，Illumina paired-end 151 bp。参考组装为 NCBI `GCF_049379645.1`，23 contigs、4,709,881 bp。ENA MD5：R1 `7449be587413dc3ee634521064759001`，R2 `3b67d6f1685535d45147f988f67324f7`；本地 SHA-256 和文件大小见 `REPORT.json`。

这是真实分离株基因组，不是混合宏基因组，因此本轮没有把 binning 或 MAG 完整性误报为已验证能力。

## 阶段结果

| 阶段 | 实测结果 | 生产解释 |
|---|---|---|
| fastp 1.4.0 | 6,733,603 对 → 6,514,384 对；17.92 s；峰值 RSS 1.34 GiB；swap 0 | 保留；必须检查 R1/R2 配对、gzip、read 数量和 JSON，而不是只看退出码 |
| SPAdes 4.3.0 | 597.96 s；峰值 RSS 5.25 GiB；169 contigs；4,731,335 bp；N50 326,264 bp；参考覆盖 100%（≥95% identity） | 保留；macOS 的 `setrlimit` 不能提供严格内存隔离，生产用 Linux cgroup/HPC |
| Bakta 1.12.2 | light DB 6.0；96.02 s；峰值 RSS 2.80 GiB；3,951 CDS、99 tRNA、1 tmRNA、8 rRNA、8 ncRNA、421 hypothetical | 保留；企业必须 full DB、固定 AMRFinderPlus 和容器依赖 |
| BLAST+ 2.17.0+ | 3,951 条 Bakta 蛋白对参考蛋白库；固定阈值 identity ≥80%、coverage ≥80%、E-value ≤1e-20 | 保留；生产需大规模、版本化数据库和 taxonomy 证据 |
| InterProScan 6 / API | 12 条候选/对照全部命中；228 个匹配；首次 2.038 s，重复 0.961 s，规范化结果完全一致；另有 1,000 条 reviewed enzyme 批量：981 个唯一序列、64,751 匹配、24.3/14.5 s；另完成 A0A194WDG1 的 Nextflow 局部 workflow smoke（8.6 s，JSON/JSONL/TSV） | 保留 API-first；未命中或新序列送 Linux full scan，不能把 API 命中称为本地重算；完整本地成员扫描在本机被 idrpred/groovy 依赖阻断 |
| OpenFold3 0.5.0 | MPS；无 MSA 4/4 成功但低置信度；同一 458 aa 序列加 MSA 后 pLDDT 39.6→78.8、pTM 0.223→0.808 | 保留为结构证据层；生产需要 MSA 缓存/服务和 Linux GPU |

## 候选识别审计

候选发现只使用 Bakta 产品词和组装蛋白对参考蛋白的 BLAST 证据；参考 GFF 标签在审计时才揭示。3,951 条组装蛋白中 3,856 条达到 BLAST 置信阈值，参考注释中有 112 条预先定义的蛋白酶/肽酶相关目标。规则识别 113 条候选，审计为 TP=92、FP=21、FN=19、TN=3,724，precision **0.814**、recall **0.829**。这组混淆矩阵只覆盖 3,856 条 BLAST 置信蛋白；其余 95 条应进入生产系统的低证据/待复核分支。

这里的参考蛋白库来自同一菌株，目的是验证 reads→assembly→annotation 的一致性以及候选规则是否能回收已知注释；它不是未知样本上的独立发现评测。生产系统应替换为版本化的 MEROPS、UniProtKB、NCBI RefSeq/NR 或企业自建非冗余库，并用没有参与规则制定的外部蛋白集进行盲测。

例子包括 M12B、M23、M24、S1、S8 和 SohB 家族候选；非目标对照包括短链脂肪酸 CoA 连接酶、组氨酸激酶和 DNA 相关蛋白。一个带有 “Regulator of protease activity HflC” 产品词的序列在参考中是 SPFH 蛋白，说明产品词本身会造成假阳性，必须由跨库家族证据和人工规则复核。

## InterPro 证据面板

| 序列 | 角色 | 长度 | Pfam | NCBIFAM | 匹配数 |
|---|---|---:|---|---|---:|
{chr(10).join(panel_lines)}

所有 12 条序列均能从 Matches API 获得结果，重复查询归一化哈希一致。面板中已知目标通常显示与其家族一致的 Pfam/NCBIFAM/IPR 组合；HflC 假阳性仅显示 SPFH 域，组氨酸激酶和 CoA 连接酶也没有被误作肽酶家族。InterPro 是功能证据层，不是活性测定。

## OpenFold3 结构批量与 MSA 对照

无 MSA、无模板的 4 条串行预测全部生成结构，分别为 EOHHJP_00526、EOHHJP_00314、EOHHJP_01311 和非目标 EOHHJP_00006。结构文件存在并可解析，但 avg pLDDT 仅约 36.1–46.1，pTM 约 0.223–0.303；这类结果不能直接用于催化位点或对接结论。批量 CLI 运行 657.9 s，峰值 RSS 4.07 GiB，峰值 MPS driver 8.55 GiB，swap 增加约 1.60 GiB。

同一 458 aa EOHHJP_01311 通过 ColabFold MSA 服务获取 MSA（约 153 s），OpenFold3 推理 260.5 s；avg pLDDT 78.8、pTM 0.808、GPDE 0.563。MSA 显著改变结构置信度，说明企业系统应把 MSA 生成、缓存和版本作为结构任务的一等输入。即便置信度提升，也不能单凭结构证明催化活性、底物特异性、分泌、耐盐或工业性能。

## 批量与异常恢复

| 测试 | 观察 | 工程含义 |
|---|---|---|
{chr(10).join(recovery_lines)}

fastp 对格式错误 FASTQ 和 R1/R2 数量不一致会出现明确日志错误或 “Ignore the unmatched reads”，但退出码仍为 0；OpenFold3 对空链查询也可能退出 0，同时摘要报告 Failed Queries。接口必须执行语义校验：日志错误扫描、输入/输出记录数、成功查询数、输出文件可解析性和 retry/idempotency。SPAdes 缺失输入、Bakta 缺失数据库或 tRNAscan-SE、BLAST 缺失数据库、InterPro 非法 FASTA 和模拟 API 连接异常均被捕获；完成的 fastp/SPAdes 产物可被后续阶段直接续跑。

## 科学准确性边界

这条链路已经证明：公开 reads 可以经过组装和注释，候选蛋白可以用 BLAST 与 InterPro 的多源证据进行筛选，OpenFold3 可以对候选生成结构并输出置信度。它没有证明候选的真实酶活、底物特异性、催化效率、表达/分泌、耐胆盐或菌株表型。参考 GFF 的 92/112 识别结果是注释层面的审计，不是生化金标准；同菌株参考 BLAST 也不能代表未知样本的独立发现准确率。下一步需要外部盲测、纯化蛋白、底物实验和必要的结构实验验证；宏基因组项目还需要 binning、污染率和完整度评估。

## 企业部署结论

整条链路的工程完成度为**可运行的研究原型，尚未达到无人值守生产服务**。建议采用 Linux + 容器/Apptainer + Nextflow：fastp、SPAdes、Bakta 和 BLAST 在 CPU 队列中运行；InterProScan 6 使用 API-first 缓存，未命中序列进入带完整数据库的 Linux 节点；OpenFold3 使用 Linux GPU 节点，并把 MSA、模型权重、数据库和镜像 digest 固定记录。Mac mini 适合接口联调、缓存 API 和小批量结构 triage，不应作为企业吞吐或内存隔离基准。

系统必须对每一步暴露标准输入/输出、版本、数据库 release、校验和、资源、日志错误、部分失败和恢复状态；Deep Agents 只能调用通过这些门禁的任务。最终输出应是“候选及证据等级、需要实验验证的清单”，不能把任何单一分数自动写成“已确认蛋白酶”。

详细机器可读记录见 [`REPORT.json`](REPORT.json)，原始运行产物位于 `results/marine-microbe-workflow/runs/`。
"""
    (RESULTS / "REPORT.md").write_text(md)
    print(json.dumps({"report": str(RESULTS / "REPORT.md"), "json": str(RESULTS / "REPORT.json")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
