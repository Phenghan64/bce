"""Export an OpenFold3 prediction as a versioned structure-viewer artifact.

The exporter keeps the raw mmCIF and writes a small manifest consumed by the
DOHO viewer. It is intentionally independent of the model runtime so it can
be used by a worker after a successful prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--artifact-id", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--organism")
    parser.add_argument("--validation-json", type=Path)
    parser.add_argument("--quality-json", type=Path)
    parser.add_argument("--model", default="OpenFold3")
    parser.add_argument("--model-version")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--templates", action=argparse.BooleanOptionalAction, default=False)
    return parser.parse_args()


def first_file(root: Path, pattern: str) -> Path | None:
    if root.is_file() and root.match(pattern):
        return root
    return next(iter(sorted(root.rglob(pattern))), None)


def load_json(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def residue_plddt_from_cif(cif_path: Path, plddt: list[float]) -> list[float]:
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("gemmi is required to export residue-level pLDDT") from exc
    block = gemmi.cif.read(str(cif_path)).sole_block()
    loop = block.find_loop("_atom_site.id").get_loop()
    tags = [str(tag) for tag in loop.tags]
    rows = [[loop[row, column] for column in range(len(tags))] for row in range(loop.length())]
    if len(rows) != len(plddt):
        raise ValueError(f"pLDDT length {len(plddt)} does not match atom count {len(rows)}")
    chain_index = tags.index("_atom_site.label_asym_id")
    seq_index = tags.index("_atom_site.label_seq_id")
    comp_index = tags.index("_atom_site.label_comp_id")
    groups: dict[tuple[str, str, str], list[float]] = {}
    order: list[tuple[str, str, str]] = []
    for atom_index, row in enumerate(rows):
        key = (row[chain_index], row[seq_index], row[comp_index])
        if key[1] in {".", "?"}:
            continue
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(float(plddt[atom_index]))
    return [sum(groups[key]) / len(groups[key]) for key in order]


def atomic_write(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        os.replace(temp_name, path)
        os.chmod(path, 0o644)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)


def main() -> None:
    args = parse_args()
    if not args.artifact_id.replace("-", "").replace("_", "").isalnum():
        raise SystemExit("artifact-id may contain only letters, numbers, '-' and '_'")
    prediction_dir = args.prediction_dir.resolve()
    source_cif = first_file(prediction_dir, "*_model.cif")
    if source_cif is None:
        raise FileNotFoundError(f"no *_model.cif found below {prediction_dir}")
    confidence_path = first_file(prediction_dir, "*_confidences_aggregated.json")
    confidence_detail_path = first_file(prediction_dir, "*_confidences.json")
    confidence = load_json(confidence_path)
    confidence_detail = load_json(confidence_detail_path)
    validation = load_json(args.validation_json.resolve() if args.validation_json else None)
    quality = load_json(args.quality_json.resolve() if args.quality_json else None)
    metrics = {
        key: confidence[key]
        for key in ("avg_plddt", "ptm", "iptm")
        if isinstance(confidence.get(key), (int, float))
    }
    metrics.update({
        key: validation[key]
        for key in (
            "protein_residues_predicted",
            "model_runtime_seconds",
            "cli_seconds",
            "protein_ca_rmsd_angstrom",
        )
        if isinstance(validation.get(key), (int, float))
    })
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    structure_name = f"{args.artifact_id}.cif"
    structure_path = output_dir / structure_name
    shutil.copy2(source_cif, structure_path)
    residue_plddt = residue_plddt_from_cif(source_cif, confidence_detail.get("plddt", []))
    confidence_payload = {
        "plddt": confidence_detail.get("plddt", []),
        "pae": confidence_detail.get("pae", []),
        "residue_plddt": residue_plddt,
    }
    atomic_write(output_dir / "confidence.json", confidence_payload)
    digest = hashlib.sha256(structure_path.read_bytes()).hexdigest()
    validation_name = "validation.json"
    atomic_write(output_dir / validation_name, {**validation, **metrics})
    manifest = {
        "schema_version": 1,
        "id": args.artifact_id,
        "title": args.title,
        "organism": args.organism,
        "structure": {"url": structure_name, "format": "mmcif", "sha256": digest},
        "validation": {"url": validation_name},
        "confidence": {
            "url": "confidence.json",
            "protein_residue_count": len(residue_plddt),
            "token_count": len(confidence_payload["pae"]),
        },
        "metrics": metrics,
        "quality": quality or None,
        "provenance": {
            "model": args.model,
            "model_version": args.model_version,
            "seed": args.seed,
            "templates": args.templates,
            "source": source_cif.name,
        },
    }
    manifest["organism"] = args.organism
    manifest["provenance"] = {key: value for key, value in manifest["provenance"].items() if value is not None}
    manifest = {key: value for key, value in manifest.items() if value is not None}
    atomic_write(output_dir / "manifest.json", manifest)
    print(json.dumps({
        "artifact_dir": str(output_dir),
        "manifest": str(output_dir / "manifest.json"),
        "structure": str(structure_path),
        "source_prediction": str(source_cif),
        "sha256": digest,
        "metrics": metrics,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
