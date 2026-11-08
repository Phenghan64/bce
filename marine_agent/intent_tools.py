"""Agent-visible deterministic intent analysis."""
from __future__ import annotations

import json

from langchain_core.tools import tool

from .intent import classify_intent


@tool
def analyze_request_intent(prompt: str) -> str:
    """Return auditable intent hints, ambiguity, and confirmation requirements."""
    return json.dumps(classify_intent(prompt), ensure_ascii=False, indent=2)


INTENT_TOOLS = [analyze_request_intent]
