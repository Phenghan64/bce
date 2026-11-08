#!/bin/bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"
mkdir -p results/1iep
./run.sh mk_prepare_ligand.py -i examples/1iep/1iep_ligand.sdf -o results/1iep/ligand.pdbqt
./run.sh mk_prepare_receptor.py -i examples/1iep/1iep_receptorH.pdb -o results/1iep/receptor -p -v \
    --box_size 20 20 20 --box_center 15.190 53.903 16.917
./run.sh vina --receptor results/1iep/receptor.pdbqt --ligand results/1iep/ligand.pdbqt \
    --config results/1iep/receptor.box.txt --exhaustiveness 32 --cpu 6 --seed 20261021 \
    --out results/1iep/docked.pdbqt 2>&1 | tee results/1iep/docking.log
./run.sh mk_export.py results/1iep/docked.pdbqt -s results/1iep/docked.sdf
