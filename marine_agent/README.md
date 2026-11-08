# Marine Protease Discovery Deep Agent

This package provides a constrained Deep Agents control plane for the validated
marine microbe workflow. It exposes typed inspection, evidence, QC, and bounded
validation tools. It deliberately does not expose arbitrary shell or filesystem
tools to the model.

`analyze_request_intent` provides deterministic routing hints for sequence
search, sequence retrieval, structure prediction, docking, enzyme discovery,
NGS, annotation, and peptide screening. It reports ambiguity and confirmation
requirements; it is advisory and does not replace model planning or scientific
quality gates.

## Offline smoke test

```bash
./.env/bin/python -m marine_agent.cli --offline
```

This runs without a model API key and validates the current artifacts and
semantic failure gates.

## Live agent

For an OpenAI-compatible internal relay, configure credentials outside the
repository. The adapter uses the Chat Completions endpoint when
`DEEP_AGENT_BASE_URL` is set:

```bash
export DEEP_AGENT_BASE_URL='https://your-relay.example/v1'
export DEEP_AGENT_API_KEY='...'
export DEEP_AGENT_MODEL='your-model-name'
export DEEP_AGENT_REQUEST_TIMEOUT_S='180'
export DEEP_AGENT_MAX_RETRIES='2'
export DEEP_AGENT_MAX_COMPLETION_TOKENS='12000'
export DEEP_AGENT_RECURSION_LIMIT='24'
./.env/bin/python -m marine_agent.cli --prompt \
  '审查当前海洋候选蛋白酶结果，列出证据不足的候选和下一步验证。'
```

`run_approved_validation` supports `full_validation`, which serializes recovery
checks before rebuilding the report. The tool also returns structured timeout,
launch, and production-readiness fields. `production_ready=false` is expected
for the current validation dataset because it uses the Bakta light database,
Apple MPS for OpenFold3, and a same-strain BLAST audit rather than an external
blind benchmark.

The agent is a coordinator. The deterministic scientific tools remain the
source of truth, and no output is allowed to imply biochemical activity without
experimental evidence.

## New protein project inputs

`run_protein_homology_screen` accepts explicit query/reference protein FASTA,
their source descriptions, an objective and task-specific BLAST thresholds.
It does not read the old marine or collagenase benchmark. Protein inputs skip
fastp/SPAdes/Bakta; missing project data or references need clarification or
sequence retrieval, not substitution with the frozen validation data.

The default research Agent exposes general input tools only. Legacy benchmark
inspection/recovery remains available through explicit benchmark routes;
`collagenase-agent-20` selects its opt-in tools. Mentioning collagenase alone
does not select the 20-proteome benchmark or inject its results into new work.

`POST /protein/runs` accepts the same fields plus optional `conversation_id`.
It creates an ordinary durable conversation task without a model-planning step.
Use the existing run/events/report/stop/continue endpoints. Agent tool calls
save the same exact inputs, and the Mac worker reports four live stages:
input QC, reference database build, BLASTP, candidate ranking/report.

Outputs live in `results/protein-screenings/<analysis_id>/`. The analysis ID
hashes original inputs, source descriptions, objective, parameters, workflow
version, BLAST versions and binary hashes. Original/normalized FASTA, commands,
QC, raw hits, candidate TSV/FASTA and an attempt ledger are retained. Repeating
the exact request verifies all stage output hashes and retries only invalid
stages and their dependents. Changing inputs or parameters creates a separate
analysis. `inspect_protein_analysis` retrieves evidence and an exact candidate
sequence for downstream domain review or `run_openfold_sequence`.

Execution uses the authenticated Mac worker's `/v1/proteins/screen` route when
`VINA_WORKER_URL` is configured (`PROTEIN_WORKER_URL/TOKEN` can override it).
The existing cancellation channel also cancels protein jobs. A standalone
execution host can use local BLAST via `PROTEIN_BLAST_BIN`; missing executables
are an explicit failure. This does not deploy a reduced scientific database.

The current JSON transport supports 32 MiB per FASTA and 200,000 records per
input, not arbitrary-size uploads. The tool returns the first 50 candidates;
the full table and report are retained. Duplicate IDs, malformed FASTA and
internal stop codons fail validation; a terminal `*` is removed and recorded.
Ambiguous residues and nucleotide-like records are counted in QC. Both query
and reference coverage are checked. Sources are recorded user assertions,
not independently certified provenance; exact reference matches are flagged.
Ranking by bitscore can favor long alignments and is not an activity ranking.
No de novo InterProScan or experimental activity validation is implied.

Checks: `./.env/bin/python3.11 -m unittest tests.test_protein_screening -v`.
The fixtures contain fresh UniProt subtilisins and unrelated proteins with
source URLs/hashes. These test new-input handling and homology behavior, not
independent enzyme-activity accuracy or production-scale throughput.

## Domain evidence for new inputs

`annotate_protein_domains` accepts either `query_fasta` plus a source, or a
previous protein-screening `analysis_id` and optional `protein_id`. The latter
reads the original candidate sequences after artifact integrity checks. Optional
`expected_pfam` accessions are compared with actual results, never filled in
from the language model. `inspect_domain_evidence` reads the saved evidence.
`POST /protein/domains/runs` exposes the same inputs with conversation tracking,
progress, stop and continuation. A chained BLAST -> domain request saves the
domain continuation only after verifying it refers to the actual BLAST inputs.

The default backend queries the InterPro Matches API using exact-sequence MD5s,
100 unique sequences per request, preserving duplicate sequence identifiers.
This is precomputed evidence, not a new InterProScan run. Unknown or omitted MD5
records produce `partial` and require a full-scan backend; they are not negative
function labels. Member-database signatures and InterPro-N predictions retain
their different sources/versions/scores. Coordinates are checked against the
submitted sequence length. A signature reported by a member database retains
its original score/e-value; no universal extra e-value cutoff is invented.

Artifacts live in `results/protein-domains/<analysis_id>/`. The specification
pins sequence hashes, source, API release information, parser version and
requested accessions. Raw API chunks, retrieval times and hashes are kept.
Retries reuse intact completed chunks for that release; `refresh=true` queries
again and a database release change creates a separate analysis. No missing
response is silently treated as an empty domain list.

Executed domain reports are rendered from actual tool records. Model prose is
retained as a separate draft in `scientific_review`, rather than publishing
invented accession/name associations or unmeasured activity claims as results.
After a homology-only task, unsupported PF/IPR identifiers trigger a factual
homology summary instead. This checks a bounded evidence type, not arbitrary
scientific reasoning, catalytic residues or industrial function. Ordinary chat
still uses the conversational path. Regression checks:
`./.env/bin/python3.11 -m unittest tests.test_domain_evidence -v`.

The registered `run_openfold_accession` tool accepts an explicit UniProt accession and runs the deterministic sequence -> MSA -> local OpenFold3 -> quality gate -> viewer manifest path. `run_openfold_sequence` provides the same bounded path for a pasted one-record protein sequence (20-4000 standard amino acids), so a chat request does not need to resolve an accession first. Ambiguous protein names are handled by the separate confirmation-aware `search_uniprot` resolver before this tool is called. Both structure tools accept an optional inline `ligands_json` list for unique chain IDs with SMILES or CCD codes. SDF paths are available only to the local CLI, so the model cannot introduce arbitrary filesystem reads.

The agent also exposes `search_uniprot`, `get_uniprot_sequence`, and
`resolve_uniprot_sequence`. A name,
function, or organism request is resolved through the UniProt REST search API;
the search returns traceable candidates and a cache record, while sequence
retrieval requires an explicit accession. The agent must ask for confirmation
when multiple candidates are plausible before starting OpenFold3. Search and
sequence metadata are provenance evidence, not proof of biochemical function.

For a direct tool smoke test:

```bash
./.env/bin/python - <<'PY'
import json
from marine_agent.openfold_tools import search_uniprot
print(search_uniprot.invoke({
    "query": "beta-glucosidase",
    "organism": "Thermotoga maritima",
    "reviewed": True,
    "limit": 5,
}))
PY
```

## AutoDock Vina standards workbench

The browser workbench at `/marine-agent/vina-workbench.html` follows the
explicit Vina sequence: prepared receptor PDBQT, one or more prepared ligand
PDBQT files, search-box center and size, then `exhaustiveness`, `num_modes`,
`energy_range`, `seed`, and CPU. `POST /vina/runs` accepts base64 PDBQT inputs,
validates them, records the exact command and provenance, and runs each ligand
as an auditable job. `GET /vina/runs/{run_id}` is a lightweight polling
endpoint; `/runs/{run_id}/report` returns ranking, poses, logs, and scientific
limits. Vina affinity values are computational ranking hypotheses, never
experimental binding constants or enzyme-activity labels.

### Execution topology

The public API and task ledger run on the Linux control node. OpenFold3 and
AutoDock Vina execution run on the Mac worker over the Tailscale address
`100.64.0.22`; the worker accepts only authenticated, bounded scientific job
payloads and returns an artifact bundle under the run ID. The worker uses a
small parallel batch pool for independent ligands, while `cpu` remains the
Vina per-docking thread count. The browser URL and API route remain unchanged.
If the worker is unavailable, the API reports a structured failure rather than
silently running the job on the control node.

### OpenFold3 receptor handoff

`GET /vina/structures/{source_run_id}` accepts a successful OpenFold3 task with
a local mmCIF viewer artifact and matching SHA-256. `POST /vina/receptors`
creates a durable asynchronous preparation job; its explicit inputs are selected
chains, original-numbered pocket residues (`A:50,A:60-65`), site evidence, nominal
pH, padding and confirmation. Poll `GET /vina/receptors/{preparation_id}`.
The rigid receptor and full preparation record are downloadable run artifacts.
`POST /vina/runs` accepts `receptor: {preparation_id: "..."}`; it checks the
successful ledger entry and receptor hash and carries provenance into Vina.

The Mac uses its existing Meeko/RDKit/Gemmi environment plus the isolated
`.tools/prep-env/bin/python` PDBFixer/OpenMM environment (override with
`RECEPTOR_PREP_PYTHON`). Standard pH-based hydrogens and terminal OXT may be
added; existing heavy atoms must stay within 0.002 Å. Missing heavy atoms,
backbone breaks, altlocs and nonstandard protein residues are rejected rather
than silently removed. The source file remains immutable; chain/residue
renumbering is recorded. Preparation and Vina share the worker lock; preparation
returns a busy error instead of silently waiting past its request timeout.

Current heuristic gates require a review-allowed manifest, mean model pLDDT
>=70, pTM>=0.5 and selected-pocket mean pLDDT>=70/minimum>=50. They permit
review, not proof of docking accuracy. This protein-only route excludes
nonprotein components after explicit confirmation. It does not infer biological
assembly, catalytic protonation, metal/cofactor parameters or missing loops;
those cases require specialist preparation. Ligand PDBQT chemistry is prepared
independently. Automated handoff regression checks use a real T4 lysozyme
prediction: `PYTHONPATH=. .env/bin/python3.11 tests/test_receptor_handoff.py`.

### PubChem reference ligand library

`GET /vina/library` exposes the frozen `pubchem20-20261030-v1` catalog (20
compounds, 21 primary reference states because ibuprofen is expanded to R/S).
Build it with `./.env/bin/python3.11 ligand_library/build_catalog.py` from the
validated preflight artifacts; publish `results/ligand-libraries/<version>/` to
the API results volume. The catalog SHA is pinned in `ligand_library.py`.

`POST /vina/runs` accepts a mixture of uploaded PDBQT records and references:
`{"library_compound_id":"L003","library_version":"pubchem20-20261030-v1"}`.
The API verifies the catalog, PDBQT and SDF hashes, expands both ibuprofen
states, checks the expanded 100-job limit and records PubChem sources, charges,
stereochemistry and preparation versions. The Mac transport preserves this
metadata in `provenance.ligand_library`. References cannot override prepared
content. Missing/corrupt library data returns 503; uploads remain available.

These are explicit reference states (carboxylates, source-state phenols), not
predictions of dominant microspecies at an assay pH. Alternative states are
retained in the preflight dataset and are not silently chosen by docking score.
Polyphenol protonation and curcumin tautomers require chemistry review for the
experimental conditions. Source data and preparation protocol remain with the
versioned library. Scores do not establish experimental activity.

Regression checks: `./.env/bin/python3.11 -m unittest tests.test_ligand_library
 tests.test_ligand_library_preflight tests.test_receptor_handoff` (on one line).
Rollout integration evidence is in `results/ligand-library-rollout-20261030/`.

### Coordinator reliability (2026-11-07)

The API uses a separate control pool (`MARINE_AGENT_CONTROL_CONCURRENCY`, default
4) for planning, conversation and status reads. Explicit scientific execution is
handed to the existing science pool (`MARINE_AGENT_MAX_CONCURRENCY`). Increasing
control concurrency does not increase OpenFold or Vina compute parallelism.
Each queued turn reads conversation messages only through its own user message.
Cancellation tracks the future across the planning-to-execution handoff.

Structured plans are validated; an invalid function-call result gets one JSON
serialization retry. If both fail, the coordinator fails explicitly without
granting tool execution authority. A rejected scientific method verdict ends
that execution attempt with the recorded assessment, rather than automatically
starting an alternative study. Input length inspection reads stored sequence
facts; result interpretation remains distinct from a progress query.

Scientific tool calls are persisted in `scientific_operations` with tool-call
IDs, exact inputs, statuses and compact result facts. Multiple operations remain
visible to conversation status tools. Automatic continuation of an ambiguous
multi-input run is deliberately not inferred; selecting and resuming arbitrary
individual operations is not yet a complete public workflow.

Run **one API coordinator process per SQLite database**. Startup claims orphaned
queued/running/cancelling runs. Local interrupted work becomes stopped with its
inputs preserved. For remote work, recovery sends a job-scoped cancellation and
checks the authenticated worker `GET /v1/jobs/{id}` until no active execution
remains. Only then is a saved-input continuation allowed. An unreachable worker
leaves `current_stage=recovery_pending`; repeating Stop retries reconciliation.
Worker status `idle` is not evidence of scientific success. This recovery path
does not import unreturned remote artifacts or resume a numerical process.

Deploy the worker status endpoints before the new API. Multi-process/replicated
coordinators require leases and distributed ownership and are not supported by
this startup recovery design. Regression evidence, deployment backups and known
remaining limits are in `results/agent-reliability-20261107/`.

### Versioned scientific explanation references (2026-11-08)

`knowledge/scientific_references.json` contains scoped facts for pLDDT, PAE,
structure-confidence interpretation and Vina output. Each fact has a stable ID
and source; source receipts include retrieval date and SHA-256. Local viewer
aggregation is explicitly distinguished from AlphaFold documentation.

The turn planner selects `reference_topics` for explanations, including contextual
follow-ups. `scientific_knowledge.py` supplies those references to the answerer
and performs one additional model consistency pass for result interpretation. Pure
reference-only conceptual replies use a single reference-grounded composition
pass, selected by `reference_only`, rather than generating and reviewing a draft.
This path is labeled `reference_grounded_composition` in audit data; it is not a
second review. The parser also accepts schema-valid JSON from the same response
when a compatible provider returns it in text instead of a tool call, without
issuing another provider request. Both paths validate reference IDs
and source URLs and attaches source links in the final Markdown. Invalid review
output or a provider failure returns labeled curated facts rather than the
unchecked draft. Transport limits are set at model client construction; ordinary
chat and scientific computation paths keep their existing configuration.

Run reports include `explanation_evidence`: knowledge version/hash, used fact IDs,
source receipts, original draft, decision, issues, review duration and sanitized
fallback reason. Reference events are separate from scientific execution events.
No new structure prediction or docking is authorized by a reference lookup.

This is **same-model reference consistency checking**, not independent scientific
validation. It covers a bounded knowledge set, does not prove task-specific
measurements, and can miss unsupported embellishments. New references need
source review and new factual evaluation cases, not additions to prohibition
lists in the general system prompt. Source collection and real-model evidence
for this release are under `results/scientific-explanations-20261108/`.


## Research projects (2026-11-08)

The project panel in doho links multiple conversations to an account-owned research
brief: name, objective, constraints and unverified hypotheses. User edits use
optimistic revision checks. There is no model write tool and no automatic promotion
of an assistant answer to a confirmed fact. A brief is context, not experimental
evidence or permission to execute its objective.

`research_projects.py` creates additive SQLite tables, preserving historical runs.
Each new run captures an immutable brief revision and up to 20 recent project task
references in the same transaction as submission. Older-reference omission is
explicit. Existing runs are not retroactively rewritten when a conversation is
linked. A read-only `inspect_project_task` tool can retrieve only snapshot-listed,
account-owned task reports. Reports remain model-written reports with provenance,
not independent validation. The run report includes `project_context`.

API: GET/POST `/projects`, GET/PUT `/projects/{id}`, POST
`/projects/{id}/conversations`, PUT `/conversations/{id}/project`. PUT of a brief
requires its current revision. Linking requires `expected_project_id` and an idle
conversation. Project brief changes during a run affect future submissions only.
Project-level questions bypass conversation-only status shortcuts; ordinary
conversation and project task permissions remain distinct.

Fixed acceptance: `tests/test_research_projects.py` covers ownership, version
conflicts, restart persistence, immutable snapshots, recent-context retention,
read-only task inspection, API authentication and source-validated composition.
Real-provider evaluations and scoped deployment scripts are under
`results/project-foundation-20261108/`. These cases are engineering acceptance,
not a held-out scientific-accuracy benchmark or end-to-end autonomous research
certification. There is no experiment data import, shared project membership,
automatic candidate entity registry or unbounded archival search in this release.
