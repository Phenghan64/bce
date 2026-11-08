"""Render bounded scientific facts from actual tool results, retaining drafts.

This validates structured domain evidence, not arbitrary scientific reasoning.
Ordinary conversation is not a deterministic computational workflow.
"""
from __future__ import annotations

import json
import re

IDENTIFIER = re.compile(r"(?<![A-Za-z0-9])(?:PF[0-9]{5}|IPR[0-9]{6})(?![A-Za-z0-9])", re.I)


def tool_payloads(state: dict) -> list[tuple[str, dict]]:
    payloads = []
    for message in state.get("messages", []):
        if message.get("type", message.get("role")) not in {"tool", "tool_result"}:
            continue
        content = message.get("content", "")
        if isinstance(content, list):
            content = "".join(block.get("text", "") for block in content if isinstance(block, dict))
        try:
            value = json.loads(content) if isinstance(content, str) else content
        except (ValueError, TypeError):
            continue
        if isinstance(value, dict):
            payloads.append((message.get("name", ""), value))
    return payloads


def method_rejection(state: dict) -> dict | None:
    for name, value in tool_payloads(state):
        if name == 'assess_collagenase_method' and value.get('decision') == 'reject':
            return {
                'decision': 'method_rejected', 'run_status': 'succeeded', 'evidence': value,
                'final_text': '不能按这个方案计算或排名。Vina 对接分数不能用于推断胶原降解活性、酶活或菌株耐胆盐能力。\n\n'
                              '可先用序列同源性和结构域证据筛选候选，必要时用 OpenFold3 检查结构，再通过实验测定活性。'
                              '请提供候选蛋白序列或明确的数据来源；本轮到此完成方法评估，尚未启动替代筛选或结构预测。',
                'scope': 'Method suitability assessment completed; requested scientific computation not performed.'}
    return None


def label(value) -> str:
    return str(value).replace('`', 'ˋ').replace('\n', ' ')


def render_screen(report: dict) -> str:
    lines = [f"蛋白同源筛选完成：{report['query_count']} 条输入，{report['candidate_count']} 条候选。"]
    for candidate in report.get("candidates", [])[:10]:
        lines.append(f"- `{label(candidate['protein_id'])}`：一致性 {candidate['identity_percent']}%，查询覆盖 {candidate['query_coverage_percent']}%，参考覆盖 {candidate['reference_coverage_percent']}%。")
    lines.extend(["当前结果支持同源候选筛选，结构域、催化位点及实验活性尚未确认。",
                  f"分析编号：`{report['analysis_id']}`"])
    return '\n\n'.join(lines)


def render_domains(report: dict) -> str:
    scanned = bool(report.get("de_novo_scan_executed"))
    lines = [f"结构域证据核验：{report['record_count']} 条蛋白，{report['unknown_record_count']} 条预计算记录未知。",
             f"来源：InterPro Matches API，数据库 {report['release']}；按原始序列指纹匹配。"]
    if scanned:
        lines = [f"InterProScan 完整扫描完成：{report['record_count']} 条蛋白，数据库 {report['release']}。",
                 "已运行官方默认可用分析，未使用预计算查询；可选授权与 ML 分析未启用。"]
    if scanned and report.get("scan_executed_this_attempt") is False:
        lines[0] = f"已载入核验过的完整扫描结果：{report['record_count']} 条蛋白，数据库 {report['release']}。"
        lines[1] = "本次复用同一输入与版本的完整扫描产物，未重新计算。"
    for record in report.get("records", []):
        observed = []
        for signature in record.get("pfam", []):
            ranges = ', '.join(f"{loc['start']}–{loc['end']}" for loc in signature.get('locations', []))
            name = signature.get('description') or signature.get('name') or '数据库未提供名称'
            values = ', '.join(str(loc['evalue']) for loc in signature.get('locations', []) if loc.get('evalue') is not None)
            observed.append(f"{signature['accession']} · {label(name)}" + (f"（{ranges}）" if ranges else '') + (f"，E-value {values}" if values else ''))
        if record["full_scan_required_for_new_sequence"]:
            result = "预计算库未给出可用记录，需要完整扫描；不能判断结构域缺失。"
        else:
            result = '；'.join(observed) or ('本次扫描没有 Pfam 匹配，不能据此推断没有相关功能。' if scanned else '该预计算记录没有 Pfam 匹配，不能据此推断没有相关功能。')
        lines.append(f"- `{label(record['protein_id'])}`：{result}")
        for accession, status in record.get("expected_pfam", {}).items():
            decision = {'observed': '已观察到匹配', 'predicted_only': '仅有模型预测，尚无对应的 Pfam 分析匹配', 'unknown': '证据未知',
                        'not_observed_in_precomputed_record': '该预计算记录未观察到匹配', 'not_observed_in_full_scan': '本次全扫描未观察到匹配'}[status]
            lines.append(f"  核验 {accession}：{decision}。")
    if report.get('record_count', 0) > report.get('returned_records', 0):
        lines.append(f"以上展示 {report['returned_records']} 条，完整记录保存在报告中。")
    lines.extend([("" if scanned else "本次查询了已有预计算结果，没有运行新的 InterProScan。") + "结构域证据不等同于催化位点、酶活或底物特异性确认。",
                  f"结构域分析编号：`{report['analysis_id']}`"])
    return '\n\n'.join(lines)


def review_answer(state: dict, draft: str, *, scope='scientific') -> dict:
    if scope == 'conversation':
        # Interpretation is not a new computation. An incidental read-only tool
        # must not replace the user's requested explanation with its raw report.
        return {'final_text':draft, 'run_status':'succeeded', 'decision':'conversation_interpretation',
                'scope':'Explanation using saved evidence, not a new scientific validation.'}
    payloads = tool_payloads(state)
    identities = [value for name,value in payloads if name == 'inspect_saved_input_identity' and value.get('records') and value.get('final_text')]
    if identities and not any(name == 'annotate_protein_domains' for name, _ in payloads):
        report = identities[-1]
        return {'final_text':report['final_text'], 'run_status':report['status'], 'decision':'rendered_from_verified_input_identity',
                'original_text':draft,'scope':'Input sequence compared to the explicitly labeled UniProt reference; not an activity determination.'}
    screens = [value for name, value in payloads if name in {'run_protein_homology_screen', 'inspect_protein_analysis'}
               and value.get('status') == 'succeeded' and value.get('analysis_id')]
    domains = [value for name, value in payloads if name in {'annotate_protein_domains', 'inspect_domain_evidence'}
               and value.get('kind') == 'protein_domain_evidence' and value.get('status') in {'succeeded', 'partial'}]
    new_domains = [value for name,value in payloads if name == 'annotate_protein_domains'
                   and value.get('kind') == 'protein_domain_evidence' and value.get('status') in {'succeeded','partial'}]
    if any(name == 'annotate_protein_domains' for name, _ in payloads):
        domains = new_domains  # Prior evidence read for context cannot replace a newly requested result.
    if domains:
        # Use actual tool outputs; a suggested PF accession in model prose is
        # never promoted into sequence evidence. Save prose separately for audit.
        unique = {value['analysis_id']: value for value in domains}
        allowed = set()
        texts = []
        for report in unique.values():
            for record in report.get('records', []):
                allowed.update(sig['accession'].split('.')[0] for sig in record.get('pfam', []))
                allowed.update(entry['accession'] for entry in record.get('interpro', []))
                allowed.update(record.get('expected_pfam', {}))
            texts.append(render_domains(report))
        return {'final_text': '\n\n'.join(texts), 'run_status': 'partial' if any(value['status'] == 'partial' for value in unique.values()) else 'succeeded',
                'decision': 'rendered_from_domain_evidence', 'original_text': draft,
                'unsupported_domain_identifiers': sorted(set(IDENTIFIER.findall(draft.upper())) - allowed),
                'scope': 'Structured domain evidence only; general scientific claims and enzyme activity are not validated.'}
    attempted = [value for name, value in payloads if name == 'annotate_protein_domains']
    if attempted:
        text = (render_screen(screens[-1]) + '\n\n') if screens else ''
        return {'final_text': text + '结构域查询未完成，本次没有可核验的结构域结果。原输入已保留，可稍后继续。',
                'run_status': 'partial' if screens else 'failed', 'decision': 'domain_lookup_failed',
                'original_text': draft, 'scope': 'A failed domain lookup does not certify protein function.'}
    if IDENTIFIER.search(draft) and screens:
        return {'final_text': render_screen(screens[-1]), 'run_status': 'succeeded',
                'decision': 'domain_claims_withheld_without_sequence_evidence', 'original_text': draft,
                'unsupported_domain_identifiers': sorted(set(IDENTIFIER.findall(draft.upper()))),
                'scope': 'The homology calculation succeeded; unverified domain explanations were not published as results.'}
    return {'final_text': draft, 'run_status': 'succeeded', 'decision': 'no_structured_domain_report',
            'scope': 'No domain evidence certification was performed.'}


def explain_activity_evidence(reports, sequence_index=None):
    """Activity sufficiency is bounded by tool evidence, not generated family lore.

    This answers whether current computation establishes activity. It does not
    estimate the real activity or synthesize new residue/family annotations.
    """
    records = [(report, record) for report in reports for record in report.get('records', [])]
    if sequence_index:
        records = records[sequence_index-1:sequence_index]
    lines = ['不能。当前计算证据尚不足以确认该输入具有目标酶活性或底物特异性，也不能据此断言它没有活性。']
    missing = not records or any(report.get('status')=='unavailable' for report in reports)
    if not records:
        lines.append('本轮未取得可核验的结构域明细。已有同源性、结构置信度或对接分数都不能代替目标底物的活性实验。')
    for report, record in records:
        lines.append(f"`{label(record['protein_id'])}`：")
        if record.get('full_scan_required_for_new_sequence'):
            lines.append('预计算库没有返回这条序列的可用记录，因此结构域证据仍是未知；这不表示结构域缺失。')
        else:
            signatures = []
            for signature in record.get('pfam',[]):
                ranges = ', '.join(f"{loc['start']}–{loc['end']}" for loc in signature.get('locations',[]))
                signatures.append(f"{signature['accession']} · {label(signature.get('description') or signature.get('name') or '名称未提供')}" + (f"（{ranges}）" if ranges else ''))
            lines.append('已记录的 Pfam 匹配：'+'；'.join(signatures)+'。这些匹配支持家族/结构域层面的判断，不能直接证明活性。' if signatures else '该结果没有可展示的 Pfam 匹配；不能据此推断缺乏目标功能。')
        evidence_source = "InterProScan 完整扫描" if report.get("de_novo_scan_executed") else "InterPro Matches API"
        lines.append(f"证据来源：{evidence_source}，版本 {report.get('release','未提供')}；分析编号 `{report['analysis_id']}`。")
    if any(record.get('full_scan_required_for_new_sequence') for _,record in records):
        from .interpro_remote import configured
        lines.append('下一步可运行完整 InterProScan 补齐新序列的结构域证据，仍不能代替活性实验。' if configured() else '下一步需要完整 InterProScan 补齐新序列的结构域证据；当前 Agent 尚未配置该执行通道，不能用重复预计算查询代替。')
    lines.append('要确认活性，需要对实际输入对应的蛋白开展目标底物实验，设置空白、已知阳性及适用的亲本对照，再评价突变或候选之间的差异。催化残基、底物谱和突变影响尚未在本轮验证。')
    return {'status':'partial' if missing else 'succeeded', 'activity_established':False, 'new_computation':False,
            'final_text':'\n\n'.join(lines), 'analysis_ids':list(dict.fromkeys(report['analysis_id'] for report,_ in records)),
            'selected_protein_ids':[record['protein_id'] for _,record in records],
            'scope':'Evidence sufficiency assessment, not biological activity prediction.'}
