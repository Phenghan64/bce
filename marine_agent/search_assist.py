"""Bounded bilingual assistance for UniProt search.

The model is used only to normalize a natural-language query and translate
returned UniProt protein names. It never invents accessions or replaces the
traceable UniProt result set.
"""
from __future__ import annotations

import json
import os
import re
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

CJK_RE = re.compile(r"[\u3400-\u9fff]")

# High-confidence domain terms avoid an unnecessary model call and make common
# Chinese scientific queries deterministic.
LOCAL_QUERY_TERMS: dict[str, tuple[str, ...]] = {
    "胆盐水解酶": ("bile salt hydrolase", "choloylglycine hydrolase"),
    "胆盐水解": ("bile salt hydrolase",),
    "胶原蛋白酶": ("collagenase",),
    "胶原酶": ("collagenase",),
    "β-葡萄糖苷酶": ("beta-glucosidase",),
    "β葡萄糖苷酶": ("beta-glucosidase",),
    "蛋白酶": ("protease",),
    "肽酶": ("peptidase",),
}

LOCAL_NAME_ZH: tuple[tuple[str, str], ...] = (
    ("bile salt hydrolase", "胆盐水解酶"),
    ("choloylglycine hydrolase", "胆盐水解酶"),
    ("conjugated bile acid hydrolase", "共轭胆汁酸水解酶"),
    ("collagenase", "胶原蛋白酶"),
    ("beta-glucosidase", "β-葡萄糖苷酶"),
    ("glucosidase", "葡萄糖苷酶"),
    ("protease", "蛋白酶"),
    ("peptidase", "肽酶"),
    ("lysozyme", "溶菌酶"),
    ("lipase", "脂肪酶"),
    ("amylase", "淀粉酶"),
    ("chitinase", "几丁质酶"),
)


def _model_name() -> str:
    return os.getenv("DEEP_AGENT_MODEL", "").split(":", 1)[-1]


def _extract_json(text: str) -> dict[str, Any]:
    cleaned = text.strip()
    if cleaned.startswith("```"):
        cleaned = re.sub(r"^```(?:json)?\s*|\s*```$", "", cleaned, flags=re.IGNORECASE | re.DOTALL)
    try:
        value = json.loads(cleaned)
        return value if isinstance(value, dict) else {}
    except json.JSONDecodeError:
        match = re.search(r"\{.*\}", cleaned, flags=re.DOTALL)
        if not match:
            return {}
        try:
            value = json.loads(match.group(0))
            return value if isinstance(value, dict) else {}
        except json.JSONDecodeError:
            return {}


def _call_model(system: str, user: str, max_tokens: int = 1200) -> dict[str, Any]:
    base_url = os.getenv("DEEP_AGENT_BASE_URL", "").rstrip("/")
    api_key = os.getenv("DEEP_AGENT_API_KEY") or os.getenv("OPENAI_API_KEY")
    model = _model_name()
    if not base_url or not api_key or not model:
        return {}
    body = json.dumps(
        {
            "model": model,
            "temperature": 0,
            "max_tokens": max_tokens,
            "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}],
        },
        ensure_ascii=False,
    ).encode("utf-8")
    request = Request(
        f"{base_url}/chat/completions",
        data=body,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    try:
        with urlopen(request, timeout=20) as response:
            payload = json.loads(response.read().decode("utf-8"))
        choices = payload.get("choices") or []
        message = choices[0].get("message") if choices and isinstance(choices[0], dict) else None
        content = message.get("content") if isinstance(message, dict) else ""
        return _extract_json(content) if isinstance(content, str) else {}
    except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError, ValueError):
        return {}


def local_query_terms(query: str) -> list[str]:
    """Return high-confidence English search terms for known Chinese terms."""
    lowered = query.strip().lower()
    terms: list[str] = []
    for zh, candidates in LOCAL_QUERY_TERMS.items():
        if zh in query or zh.lower() in lowered:
            for candidate in candidates:
                if candidate not in terms:
                    terms.append(candidate)
    return terms


def assist_query(query: str, organism: str = "") -> dict[str, Any]:
    """Normalize Chinese natural-language input without inventing identifiers."""
    query = query.strip()
    local = local_query_terms(query)
    if local:
        return {"search_terms": local, "ai_used": False, "method": "controlled_dictionary"}
    if not CJK_RE.search(query):
        return {"search_terms": [query], "ai_used": False, "method": "original_query"}
    result = _call_model(
        """You normalize biological search text for UniProt. Return JSON only:
{"english_queries":["..."],"query_zh":"..."}.
Translate the user's Chinese protein/function description into at most two
specific English UniProt search phrases. Do not invent UniProt accessions,
organisms, experimental claims, or gene names. Keep queries concise.""",
        json.dumps({"query": query, "organism": organism}, ensure_ascii=False),
        max_tokens=500,
    )
    terms = [str(item).strip() for item in result.get("english_queries", []) if str(item).strip()][:2]
    if not terms:
        terms = [query]
    return {
        "search_terms": terms,
        "query_zh": str(result.get("query_zh") or query),
        "ai_used": bool(result),
        "method": "model_assist" if result else "original_query_fallback",
        "model": _model_name() if result else None,
    }


def _local_name_zh(name: str | None) -> str | None:
    if not name:
        return None
    lowered = name.lower()
    for needle, translation in LOCAL_NAME_ZH:
        if needle in lowered:
            return translation
    return None


def translate_candidate_names(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    """Attach Chinese display names while retaining the original UniProt text."""
    names: dict[str, str] = {}
    missing: list[dict[str, str]] = []
    for candidate in candidates[:20]:
        accession = str(candidate.get("accession") or "")
        english = str(candidate.get("protein_name") or candidate.get("entry_name") or "")
        local = _local_name_zh(english)
        if accession and local:
            names[accession] = local
        elif accession and english:
            missing.append({"accession": accession, "name": english})
    ai_used = False
    if missing:
        result = _call_model(
            """Translate UniProt protein names into concise Simplified Chinese.
Return JSON only in the form {"translations":{"ACCESSION":"中文名称"}}.
Translate literally and conservatively. Do not infer activity, substrate,
organism, or an accession that is not supplied.""",
            json.dumps({"candidates": missing}, ensure_ascii=False),
            max_tokens=1000,
        )
        translations = result.get("translations") if isinstance(result.get("translations"), dict) else {}
        for accession, value in translations.items():
            if accession in {item["accession"] for item in missing} and isinstance(value, str) and value.strip():
                names[accession] = value.strip()[:120]
        ai_used = bool(translations)
    return {"names": names, "ai_used": ai_used, "model": _model_name() if ai_used else None}
