"""Command line entry point for the controlled marine Deep Agent."""
from __future__ import annotations

import argparse
import json
import os
import sys

from .agent import COLLAGENASE_SYSTEM_PROMPT, TOOLS, build_agent
from .collagenase_tools import COLLAGENASE_TOOLS
from .tools import TOOLS as VALIDATION_TOOLS


def _json_safe(value):
    """Preserve message/tool-call structure without leaking provider objects."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    model_dump = getattr(value, "model_dump", None)
    if callable(model_dump):
        return _json_safe(model_dump())
    return str(value)


def offline_report() -> dict:
    """Run a deterministic smoke test without an LLM API key."""
    inspect = next(t for t in VALIDATION_TOOLS if t.name == "inspect_workflow")
    validate = next(t for t in VALIDATION_TOOLS if t.name == "validate_artifacts")
    return {
        "mode": "offline",
        "agent_runtime": "deepagents 0.7.30",
        "tool_count": len(TOOLS),
        "tool_names": [t.name for t in TOOLS],
        "inspect_workflow": json.loads(inspect.invoke({})),
        "validate_artifacts": json.loads(validate.invoke({})),
        "next_action": "Configure a provider API key and DEEP_AGENT_MODEL for natural-language planning; keep execution on the approved tools.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true", help="run artifact/QC smoke test without an LLM")
    parser.add_argument("--prompt", help="natural-language task for the live Deep Agent")
    parser.add_argument("--model", help="provider:model, e.g. openai:gpt-5.4")
    args = parser.parse_args()
    if args.offline or not args.prompt:
        print(json.dumps(offline_report(), ensure_ascii=False, indent=2))
        return 0
    try:
        prompt_lower = args.prompt.lower()
        collagenase_task = "collagenase-agent-20" in prompt_lower
        agent = build_agent(
            args.model,
            tools=COLLAGENASE_TOOLS if collagenase_task else None,
            system_prompt=COLLAGENASE_SYSTEM_PROMPT if collagenase_task else None,
        )
        recursion_limit = int(os.getenv("DEEP_AGENT_RECURSION_LIMIT", "24"))
        result = agent.invoke(
            {"messages": [{"role": "user", "content": args.prompt}]},
            config={"recursion_limit": recursion_limit},
        )
    except Exception as exc:  # clear CLI failure, no traceback required for users
        print(json.dumps({
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
            "model": args.model or os.getenv("DEEP_AGENT_MODEL", "openai:gpt-5.4"),
            "recursion_limit": int(os.getenv("DEEP_AGENT_RECURSION_LIMIT", "24")),
        }, ensure_ascii=False, indent=2))
        return 2
    print(json.dumps(_json_safe(result), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
