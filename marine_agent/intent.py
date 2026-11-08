"""Deterministic request-intent hints for the scientific Deep Agent.

This is a routing aid, not a replacement for model planning. It deliberately
returns uncertainty and confirmation requirements instead of silently choosing
an accession, tool, or scientific interpretation.
"""
from __future__ import annotations

import re


ACCESSION_RE = re.compile(
    r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|(?=[A-Z0-9]{10}\b)(?=[A-Z0-9]*[0-9])[A-Z0-9]{10})\b",
    re.IGNORECASE,
)
RAW_SEQUENCE_RE = re.compile(r"(?<![A-Za-z])([ACDEFGHIKLMNPQRSTVWY]{20,4000})(?![A-Za-z])", re.IGNORECASE)

INTENT_RULES = {
    "homology_screening": ("同源", "同源筛选", "homolog", "protein screening"),
    "sequence_search": ("查序列", "找序列", "搜索序列", "查找", "查", "蛋白名称", "protein name", "find sequence", "search sequence", "uniprot search"),
    "sequence_retrieval": ("获取序列", "下载序列", "取序列", "fasta", "sequence", "序列", "accession"),
    "structure_prediction": ("结构预测", "预测结构", "openfold", "alphafold", "三维结构", "structure prediction"),
    "docking": ("对接", "vina", "autodock", "小分子结合", "docking"),
    "enzyme_discovery": ("找酶", "发现酶", "候选酶", "酶筛选", "蛋白酶", "胶原酶", "enzyme discovery", "protease"),
    "ngs_workflow": ("测序", "fastq", "ngs", "组装", "spades", "fastp", "基因组分析", "metagenome"),
    "annotation": ("注释", "blast", "interpro", "pfam", "bakta", "功能注释", "annotation", "结构域", "domain"),
    "peptide_screening": ("功能肽", "精准制肽", "虚拟酶解", "肽筛选", "peptide", "in silico digestion"),
}


RECOMMENDED_TOOLS = {
    "homology_screening": ["run_protein_homology_screen", "inspect_protein_analysis"],
    "sequence_search": ["search_uniprot", "resolve_uniprot_sequence"],
    "sequence_retrieval": ["get_uniprot_sequence", "resolve_uniprot_sequence"],
    "structure_prediction": ["resolve_uniprot_sequence", "run_openfold_accession"],
    "docking": ["method suitability review", "validated receptor/ligand preparation", "AutoDock Vina"],
    "enzyme_discovery": ["BLASTP", "InterPro/Pfam", "get_uniprot_sequence", "OpenFold3 gated review"],
    "ngs_workflow": ["fastp", "SPAdes", "Bakta", "BLASTP", "InterPro/Pfam"],
    "annotation": ["annotate_protein_domains", "inspect_domain_evidence", "BLASTP", "Bakta"],
    "peptide_screening": ["protease-specific cleavage rules", "peptide scoring", "experimental validation"],
}


def classify_intent(prompt: str) -> dict:
    text = " ".join(str(prompt or "").strip().lower().split())
    if not text:
        return {"status": "invalid_input", "reason": "prompt is empty"}
    scores = {
        intent: sum(1 for token in tokens if token.lower() in text)
        for intent, tokens in INTENT_RULES.items()
    }
    matched = sorted(((score, intent) for intent, score in scores.items() if score), key=lambda item: -item[0])
    accession = ACCESSION_RE.search(text)
    inline_sequence = RAW_SEQUENCE_RE.search(str(prompt or ""))
    if str(prompt or "").count('>') > 1 and scores["homology_screening"]:
        return {"status": "succeeded", "primary_intent": "homology_screening",
                "candidate_intents": ["homology_screening"], "confidence": "high",
                "requires_confirmation": False, "accession": None,
                "recommended_tools": RECOMMENDED_TOOLS["homology_screening"],
                "safety_flags": [], "reason": "Multiple protein FASTA records accompany a homology task; use supplied query/reference inputs, not a single-protein prediction or a frozen benchmark."}
    if inline_sequence and any(token in text for token in ('结构域','interpro','pfam','domain')):
        return {"status":"succeeded", "primary_intent":"annotation", "candidate_intents":["annotation"],
                "confidence":"high", "requires_confirmation":False, "accession":None,
                "recommended_tools":["annotate_protein_domains"], "safety_flags":[],
                "reason":"Provided protein input with domain-review intent; use sequence evidence rather than the single-structure shortcut."}
    if inline_sequence and any(token in text for token in ("结构", "预测", "openfold", "alphafold", "structure")):
        return {
            "status": "succeeded",
            "primary_intent": "structure_prediction",
            "candidate_intents": ["structure_prediction"],
            "confidence": "high",
            "requires_confirmation": False,
            "accession": None,
            "sequence_length": len(inline_sequence.group(1)),
            "recommended_tools": ["run_openfold_sequence"],
            "safety_flags": [],
            "reason": "A bounded inline protein sequence was detected; route directly to MSA/OpenFold3 without UniProt name resolution.",
        }
    if not matched:
        return {
            "status": "succeeded",
            "primary_intent": "unknown",
            "candidate_intents": [],
            "confidence": "low",
            "requires_confirmation": True,
            "accession": accession.group(0).upper() if accession else None,
            "recommended_tools": [],
            "reason": "No registered intent pattern matched; the Agent should ask a focused clarification question.",
        }
    top_score = matched[0][0]
    candidates = [intent for score, intent in matched if score == top_score or score >= max(1, top_score - 1)]
    primary = matched[0][1]
    confidence = "high" if top_score >= 2 and len(candidates) == 1 else "medium" if top_score >= 1 else "low"
    enzyme_activity = any(token in text for token in ("酶活", "活性", "降解", "collagenase", "胶原酶", "enzyme activity", "degradation"))
    requires_confirmation = (
        "sequence_search" in candidates and accession is None
    ) or len(candidates) > 1
    safety_flags = []
    if enzyme_activity and ("docking" in candidates or "docking" in primary):
        safety_flags.append("docking_score_cannot_establish_enzyme_activity")
    if "structure_prediction" in candidates and accession is None:
        safety_flags.append("explicit_accession_or_confirmed_search_result_required")
    return {
        "status": "succeeded",
        "primary_intent": primary,
        "candidate_intents": candidates,
        "confidence": confidence,
        "requires_confirmation": requires_confirmation,
        "accession": accession.group(0).upper() if accession else None,
        "recommended_tools": sorted({tool for intent in candidates for tool in RECOMMENDED_TOOLS[intent]}),
        "safety_flags": safety_flags,
        "reason": "Deterministic keyword hints are advisory; the model must inspect tool evidence and preserve scientific gates.",
    }
