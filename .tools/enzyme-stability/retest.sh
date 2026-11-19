#!/bin/bash
set -eu
cd "$(dirname "$0")/../.."
base="$PWD/.tools/enzyme-stability"
out="$PWD/results/enzyme-stability"
python="$PWD/.tools/openfold3-env/bin/python"
export OMP_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4
PYTHONPATH="$PWD/.tools/foundry-eval/extras" "$python" "$base/thermo_scan.py"
for variant in hyper baseline; do
  weights="$base/HyperMPNN-weights"
  model=v48_020_epoch300_hyper
  if [[ "$variant" == baseline ]]; then
    weights="$base/ThermoMPNN-main/vanilla_model_weights"
    model=v_48_020
  fi
  "$python" "$base/proteinmpnn/protein_mpnn_run.py" --pdb_path "$out/1BLI.pdb" --pdb_path_chains A --out_folder "$out/$variant" --path_to_model_weights "$weights" --model_name "$model" --seed 557 --sampling_temp 0.2 --num_seq_per_target 4 --batch_size 4 --fixed_positions_jsonl "$out/fixed.jsonl"
done
