#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
export ENZYME_CASE_OUTPUT="$PWD/results/enzyme-case-repeat-$(date +%Y%m%d-%H%M%S)"
mkdir "$ENZYME_CASE_OUTPUT"
ln -s "$PWD/results/enzyme-case/source" "$ENZYME_CASE_OUTPUT/source"
cp results/enzyme-case/model-compatibility.json results/enzyme-case/environment.txt results/enzyme-case/machine.json "$ENZYME_CASE_OUTPUT/"
date -u +%Y-%m-%dT%H:%M:%SZ > "$ENZYME_CASE_OUTPUT/started-utc.txt"
.tools/enzyme-case/venv/bin/python tmgh1/monitor.py "$ENZYME_CASE_OUTPUT/compute" .tools/enzyme-case/venv/bin/python enzyme_case/run.py
.tools/enzyme-case/venv/bin/python tmgh1/monitor.py "$ENZYME_CASE_OUTPUT/audit-compute" .tools/enzyme-case/venv/bin/python enzyme_case/audit.py
.tools/enzyme-case/venv/bin/python enzyme_case/report.py
echo "Results: $ENZYME_CASE_OUTPUT/report.html"
