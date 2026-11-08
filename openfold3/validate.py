"""Check ubiquitin predictions and compare all 76 C-alpha atoms with 1UBQ."""
import csv
import json
from pathlib import Path

import numpy as np
from biotite.structure.io import pdbx

root = Path(__file__).resolve().parent.parent
results = root / "results/openfold3-ubiquitin"
reference = root / ".tools/openfold3/openfold3/tests/test_data/mmcifs/1ubq.cif"


def ca_atoms(path):
    atoms = pdbx.get_structure(pdbx.CIFFile.read(path), model=1)
    assert np.isfinite(atoms.coord).all(), "Non-finite coordinates"
    return atoms[(atoms.atom_name == "CA") & (atoms.chain_id == "A") & ~atoms.hetero]


ref = ca_atoms(reference)
assert len(ref) == 76
rows = []
for mode in ["single-sequence", "with-msa"]:
    folder = results / mode
    seed = folder / "prediction/ubiquitin/seed_42"
    prefix = "ubiquitin_seed_42_sample_1"
    predicted = ca_atoms(seed / f"{prefix}_model.cif")
    assert len(predicted) == 76
    assert np.array_equal(predicted.res_id, ref.res_id)
    assert np.array_equal(predicted.res_name, ref.res_name)
    p = predicted.coord.astype(float)
    q = ref.coord.astype(float)
    p -= p.mean(axis=0)
    q -= q.mean(axis=0)
    u, _, vt = np.linalg.svd(p.T @ q)
    correction = np.diag([1.0, 1.0, np.linalg.det(u @ vt)])
    aligned = p @ (u @ correction @ vt)
    rmsd = np.sqrt(np.mean(np.sum((aligned - q) ** 2, axis=1)))
    confidence = json.loads((seed / f"{prefix}_confidences_aggregated.json").read_text())
    benchmark = json.loads((folder / "benchmark.json").read_text())
    assert benchmark["cli_completed"] and benchmark["mps_available"]
    config = json.loads((folder / "prediction/experiment_config.json").read_text())
    assert config["pl_trainer_args"]["accelerator"] == "mps"
    rows.append({
        "mode": mode,
        "residues": len(predicted),
        "cli_seconds": benchmark["elapsed_seconds"],
        "model_runtime_seconds": json.loads((seed / "timing.json").read_text())["runtime_s"],
        "process_peak_rss_gib": benchmark["process_peak_rss_bytes"] / 1024**3,
        "sampled_peak_mps_driver_gib": benchmark["sampled_peak_mps_driver_bytes"] / 1024**3,
        "system_swap_change_mib": (benchmark["system_swap_after_bytes"] - benchmark["system_swap_before_bytes"]) / 1024**2,
        "avg_plddt": confidence["avg_plddt"],
        "ptm": confidence["ptm"],
        "ca_rmsd_to_1ubq_angstrom": float(rmsd),
        "model_reports_clash": confidence["has_clash"],
    })
report = {
    "target": "Ubiquitin, chain A, 76 residues",
    "source": "Official aqlaboratory/openfold-3 main source archive, downloaded 2026-10-21; archive digest in openfold3/downloads.sha256",
    "checkpoint": "of3-ob-2025-06-30-174k.pt (current upstream default OpenBind checkpoint)",
    "settings": "MPS, float32, 1 seed (42), 1 diffusion sample, low_mem, no templates",
    "rmsd_method": "All 76 matching C-alpha atoms, Kabsch rigid alignment, no trimming or outlier removal",
    "caveats": [
        "Small known protein; this does not establish unseen-target accuracy or large-complex capacity",
        "MSA was obtained from the public ColabFold server; model inference was local",
        "Second run benefits from warmed caches; speed difference is not evidence that MSA accelerates inference",
        "Model timing excludes setup/MSA downloads; process RSS and MPS driver memory overlap and must not be added",
    ],
    "runs": rows,
}
(results / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
with (results / "comparison.csv").open("w") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(report, indent=2))
