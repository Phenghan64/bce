#!/usr/bin/env python3
"""Build an auditable report for the independent collagenase Agent benchmark."""
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/collagenase-agent-20"
RUN = OUT / "runs/batch-screen"


def main() -> None:
    manifest = json.loads((OUT / "manifest.json").read_text())
    screening = json.loads((RUN / "screening.json").read_text())
    interpro = json.loads((RUN / "interpro-candidates.json").read_text())
    recovery = json.loads((OUT / "runs/recovery.json").read_text())
    candidates = list(csv.DictReader((RUN / "candidate-table.tsv").open(), delimiter="\t"))
    by_ipr = {x["id"]: x for x in interpro["run"]["summaries"]}
    for row in candidates:
        row["interpro_match_count"] = by_ipr.get(row["qseqid"], {}).get("match_count", 0)
        row["interpro_pfam"] = by_ipr.get(row["qseqid"], {}).get("pfam", [])
        row["interpro_ipr"] = by_ipr.get(row["qseqid"], {}).get("interpro", [])
    structure_files = list((OUT / "runs/openfold").glob("*/structure.json"))
    structures = [json.loads(path.read_text()) for path in structure_files]
    report = {
        "title": "20 株未参与调试的海洋微生物潜在胶原蛋白酶自主科研 Agent 测试",
        "workflow": "collagenase-agent-20",
        "data_source": "UniProt public proteomes",
        "proteomes": manifest["proteomes"],
        "input_route": {"mode": "existing_protein_sequences", "ngs_stages_skipped": ["fastp", "SPAdes"], "reason": "input already contained protein sequences"},
        "reference_database": manifest["collagenase_references"],
        "screening": screening,
        "interpro": {"sequence_count": interpro["input"]["sequence_count"], "found_count": interpro["run"]["found_count"], "matches": interpro["run"]["total_match_count"], "repeat_identical": interpro["repeat"]["byte_level_normalized_match"], "runtime_seconds": interpro["run"]["elapsed_seconds"]},
        "top_candidates": candidates[:20],
        "structures": structures,
        "recovery": recovery,
        "method_selection": {"vina_for_collagenase_activity": "rejected", "reason": "Vina is a small-ligand docking method and cannot establish enzyme activity or collagen substrate degradation", "preferred_route": ["BLASTP", "InterPro/Pfam", "MSA-backed OpenFold3 for structural review", "expression and biochemical collagenolysis assay"]},
        "scientific_limits": [
            "The 20 proteomes are independent of the earlier SSi026 workflow, but the post-hoc labels are UniProt annotations rather than wet-lab activity measurements.",
            "The external collagenase reference set excludes target organism names, but homology leakage and annotation transfer remain possible.",
            "BLASTP and InterPro identify family-level evidence; they do not prove catalytic activity, substrate specificity, secretion, stability or industrial performance.",
            "OpenFold3 confidence supports fold/domain review only. It is not an activity classifier.",
            "The sequence-input route intentionally skips NGS; it does not validate raw-read assembly or genome completeness for these 20 organisms.",
        ],
        "enterprise_assessment": {"agent_autonomy": "validated for bounded planning and tool selection", "human_in_loop": True, "production_ready": False, "blockers": ["full external database/version governance", "Linux/HPC/GPU execution for scale", "wet-lab activity confirmation", "RBAC/queue/audit service layer"]},
    }
    (OUT / "REPORT.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    lines = ["# 20 株海洋微生物潜在胶原蛋白酶 Agent 测试", "", f"输入：{screening['proteome_count']} 株真实公开蛋白组，{screening['total_proteins']} 条蛋白；由于已有蛋白序列，跳过 fastp 和 SPAdes。", f"BLASTP 筛选候选：{screening['candidate_count']} 条（high={screening['candidate_count_by_tier']['high']}，medium={screening['candidate_count_by_tier']['medium']}，low={screening['candidate_count_by_tier']['low']}）。", f"候选中后验 UniProt collagenase 标签：{screening['known_collagenase_posthoc_count']} 条；该标签没有参与筛选规则。", f"InterPro：{interpro['input']['sequence_count']} 条，{interpro['run']['total_match_count']} 个匹配，重复结果一致={interpro['repeat']['byte_level_normalized_match']}。", f"OpenFold3 结构任务：{len(structures)} 条；只对通过 PF01752 结构域门的候选执行。", "", "## 方法选择", "Vina 被拒绝用于预测胶原蛋白酶活性；它不能替代底物降解实验。推荐路径是 BLASTP + InterPro/Pfam + MSA/OpenFold3 结构复核 + 湿实验。", "", "## 结论", "Agent 能从自然语言目标选择序列证据、批量筛选、结构复核、异常恢复和方法拒绝；当前结论仍是候选生成，不是活性确认。"]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"report": str((OUT / "REPORT.md").relative_to(ROOT)), "json": str((OUT / "REPORT.json").relative_to(ROOT)), "structures": len(structures)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
