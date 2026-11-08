#!/bin/bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$project_dir/.tools/bin:$PATH"
if [ "$#" -eq 0 ]; then
    printf 'Usage: ./run.sh <command> [arguments...]\nExample: ./run.sh vina --version\n' >&2
    exit 2
fi
exec "$project_dir/.tools/bin/micromamba" run -p "$project_dir/.env" -- "$@"
