"""Semantic action planning separates conversation from execution authority."""
from __future__ import annotations

import json
from typing import Literal

from langchain_core.messages import HumanMessage, SystemMessage
from pydantic import BaseModel, Field

from .agent import configured_model
from .execution_control import check_cancelled
from .scientific_knowledge import ReferenceTopic


class TurnPlanningError(RuntimeError):
    """Provider output could not establish execution authority."""


class TurnPlan(BaseModel):
    action: Literal["reply", "inspect", "execute", "continue"] = Field(
        description="Requested action: discussion/explanation; inspect existing work; perform new research; or continue interrupted research.")
    reference_topics: list[ReferenceTopic] = Field(default_factory=list, max_length=4, description='Reference topics needed to answer this turn: plddt (local confidence and residue/atom granularity), pae (relative-position uncertainty/heatmap), structure_confidence (pTM/ipTM, prediction quality vs experimental activity), vina_interpretation (docking scores, pose RMSD and experimental meaning). Resolve follow-ups using context. Empty for greetings, unrelated discussion, status queries or tool-only requests.')
    context_scope: Literal['conversation','project'] = Field(default='conversation', description='Which saved context is being asked about: project for research goals, constraints, hypotheses, candidate identities, imported experiments, project-wide results or other linked conversations; conversation for the current chat and its specific tasks.')
    reference_only: bool = Field(default=False, description='True only for a conceptual question fully answerable from the selected reference topics, including conceptual follow-ups. False when needing task results, project facts, measurements, sequence lookup, new literature, status or computation.')
    requested_method: Literal['unspecified','vina','openfold','sequence_analysis'] = Field(default='unspecified', description='Method explicitly requested in this turn, canonicalized across languages; vina includes AutoDock docking.')
    scientific_endpoint: Literal['unspecified','binding_pose','binding_score','enzyme_activity','substrate_degradation','organism_phenotype','protein_structure','sequence_evidence'] = Field(default='unspecified', description='The outcome the user wants to infer. Enzyme activity or activity ranking is enzyme_activity even when docking scores are proposed as its proxy; a ligand binding pose is binding_pose.')
    source_run_id: str = Field(default="", description="An unambiguous existing task selected for continuation; otherwise empty.")
    operation_id: str = Field(default='', description='For continuation: explicit saved operation ID, or empty.')
    operation_index: int | None = Field(default=None, ge=1, description='For continuation: explicit 1-based operation number selected by the user, or null.')
    continuation_method: Literal['saved','full_scan'] = Field(default='saved', description='full_scan when the user requires actual InterProScan, not precomputed lookup.')
    reason: str = Field(description="One short sentence describing the user's requested action.")
    completion_scope: Literal['open', 'provided_fasta_domain_lookup'] = Field(default='open', description='Use provided_fasta_domain_lookup ONLY when the entire requested work is InterPro/Pfam precomputed domain lookup of supplied FASTA inputs, with no additional identity, structure, activity, literature, or full-scan task. Otherwise open.')
    recompute: bool = Field(default=False, description='True only when the user explicitly requests fresh recomputation or repeated independent executions of the same inputs; ordinary new requests use false.')
    repeated_trials: bool = Field(default=False, description='True only when the user explicitly requests multiple independent executions of identical inputs in this turn. Refreshing once is false.')
    inspection: Literal['status', 'sequence_length', 'protein_identity', 'activity_evidence', 'results'] = Field(
        default='status', description='For inspect: progress/status; exact input sequence_length; protein_identity for saved input identity/mutation; activity_evidence when asking whether existing computed evidence establishes enzyme activity or substrate specificity; results for other interpretation.')
    sequence_index: int | None = Field(default=None, ge=1, description='For sequence_length, protein_identity or activity_evidence: explicit 1-based input sequence ordinal selected by the user (e.g. 第一条 = 1); null lists all inputs when no sequence was selected.')


PLANNING_PROMPT = """Classify the action requested by the latest user message using conversation context and server task facts.
reply: greetings, conceptual explanations, discussion, or asking about capability.
inspect: asking what was executed, whether work continued, input identity, or what existing evidence establishes; select the corresponding inspection operation.
execute: requesting new scientific work with tools, including sequence lookup.
continue: requesting that interrupted scientific work actually run again. '继续' refers to the latest unresolved work when context is clear; '继续解释' continues discussion.
Classify the user's present request, rather than treating past commands or past assistant suggestions as new authorization.
Distinguish '你有继续操作吗' (a status inquiry) from '继续操作' (an execution request).
Return a structured plan. An ambiguous continuation may leave source_run_id empty for the coordinator to resolve or clarify.
"""


CONVERSATION_PROMPT = """You are EYOSON, a conversational scientific assistant.
Answer the current question naturally and concisely using the conversation and authoritative task facts scoped to this conversation.
Keep the answer focused on that question; use task details when they help answer it.
The inspect_task_state tool reads actual task status in this conversation. Computational capabilities are available through the task coordinator when execution is requested.
Explain stopped work accurately: OpenFold prediction can restart with saved inputs and existing caches; Vina can retain validated completed ligand results and retry the remainder. This is distinct from resuming an interrupted numerical process.
"""


def plan_turn(prompt: str, messages: list[dict], task_state: dict) -> TurnPlan:
    model = configured_model()
    if isinstance(model, str):
        from langchain.chat_models import init_chat_model
        model = init_chat_model(model)
    dialogue = [{"role":m["role"], "content":m["content"][-2500:]} for m in messages[-8:] if m["role"] in {"user","assistant"}]
    inputs = [
        SystemMessage(content=PLANNING_PROMPT),
        HumanMessage(content=json.dumps({"latest_message":prompt,"recent_dialogue":dialogue,"server_task_state":task_state}, ensure_ascii=False)),
    ]
    # Some OpenAI-compatible providers return None when no function call was
    # emitted. Never infer execution authority from that empty response.
    try:
        result = model.with_structured_output(TurnPlan, method="function_calling").invoke(inputs)
        return TurnPlan.model_validate(result)
    except Exception as first_error:
        check_cancelled()
        # One alternate serialization attempt, with the same action schema.
        # This also supports providers with incomplete tool-calling support.
        try:
            response = model.invoke([*inputs, SystemMessage(content=
                'Return only one JSON object conforming to this schema: ' + json.dumps(TurnPlan.model_json_schema()))])
            text = response.content
            if not isinstance(text, str):
                raise ValueError('non-text planning response')
            text = text.strip()
            if text.startswith('```') and text.endswith('```'):
                text = text.split('\n', 1)[1].rsplit('```', 1)[0].strip()
            return TurnPlan.model_validate_json(text)
        except Exception as second_error:
            raise TurnPlanningError('暂时无法可靠理解本轮操作，请重试。没有因此启动新的科研计算。') from second_error
