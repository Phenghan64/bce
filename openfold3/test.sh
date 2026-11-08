#!/bin/bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
query="$project_dir/results/openfold3-ubiquitin/msa-long-timeout/query_msa.json"
if [ ! -f "$query" ]; then
    ./openfold3/run.sh python -u openfold3/fetch_msa.py openfold3/ubiquitin.json \
        results/openfold3-ubiquitin/msa-long-timeout
fi
output="${1:-$project_dir/results/openfold3-ubiquitin/repeat-$(date +%Y%m%d-%H%M%S)}"
if [ -e "$output" ]; then
    printf 'Output already exists; choose a new directory: %s\n' "$output" >&2
    exit 2
fi
mkdir -p "$output"
./openfold3/run.sh python -u openfold3/benchmark.py "$output/benchmark.json" predict \
    --query-json "$query" --runner-yaml openfold3/runner.yml \
    --num-diffusion-samples 1 --use-msa-server false --use-templates false \
    --output-dir "$output/prediction" 2>&1 | tee "$output/run.log"
