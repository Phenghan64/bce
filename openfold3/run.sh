#!/bin/bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export OPENFOLD_CACHE="$project_dir/.tools/openfold3-cache"
# Conda's copied Python/OpenSSL build can retain the original environment
# prefix in its compiled default paths.  Pin the runtime paths to the active
# project so LaunchAgent jobs outside the Documents directory do not block
# while OpenSSL tries to open an inaccessible config/certificate file.
export OPENSSL_CONF="$project_dir/.tools/openfold3-env/ssl/openssl.cnf"
export SSL_CERT_FILE="$project_dir/.tools/openfold3-env/ssl/cert.pem"
export SSL_CERT_DIR="$project_dir/.tools/openfold3-env/ssl/certs"
# Prefer the relocated source tree over the editable-install path embedded in
# the copied conda environment.  This keeps launchd jobs independent of the
# original Documents checkout.
export PYTHONPATH="$project_dir/.tools/openfold3${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=6
export VECLIB_MAXIMUM_THREADS=6
exec "$project_dir/.tools/bin/micromamba" run -p "$project_dir/.tools/openfold3-env" -- "$@"
