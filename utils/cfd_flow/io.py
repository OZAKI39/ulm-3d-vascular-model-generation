"""Input, provenance, and run-layout I/O for the CFD-flow stage."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class FlowInputs:
    run_root: Path
    tagged_surface_vtp: Path
    meter_surface_stl: Path
    boundary_conditions_json: Path
    boundary_manifest_csv: Path
    boundary_conditions: dict[str, Any]
    boundary_manifest: tuple[dict[str, str], ...]


@dataclass(frozen=True, slots=True)
class RunLayout:
    root: Path
    input: Path
    geometry: Path
    geometry_solver_m: Path
    seeder: Path
    seeder_mesh: Path
    musubi: Path
    flow: Path
    qc: Path
    figures: Path
    proteus: Path


class FlowError(RuntimeError):
    """Expected production failure with a stable public status."""

    def __init__(self, status: str, message: str) -> None:
        self.status = status
        super().__init__(f"{status}: {message}")


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing required JSON: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_flow_inputs(source_run: Path) -> FlowInputs:
    """Load exactly the accepted surface, BC JSON, and boundary manifest."""

    root = Path(source_run).resolve()
    tagged = root / "geometry" / "cfd_surface_vmtk_tps_boundarynormal_crossseam_um.vtp"
    meter = root / "geometry" / "cfd_surface_vmtk_tps_boundarynormal_crossseam_m.stl"
    bc_path = root / "bc" / "boundary_conditions_vmtk_boundarynormal_crossseam.json"
    manifest_path = root / "boundaries" / "boundary_manifest.csv"
    for path in (tagged, meter, bc_path, manifest_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing frozen CFD-flow input: {path}")
    with manifest_path.open(encoding="utf-8-sig", newline="") as stream:
        manifest = tuple(csv.DictReader(stream))
    if len(manifest) != 4:
        raise ValueError(f"Expected 4 cap records, found {len(manifest)}")
    bc = read_json(bc_path)
    return FlowInputs(root, tagged, meter, bc_path, manifest_path, bc, manifest)


def create_run_layout(
    output_root: Path,
    *,
    timestamp: datetime | None = None,
    recovery: bool = False,
    solver_recovery: bool = False,
) -> RunLayout:
    """Create the one concise production or explicitly authorized recovery tree."""

    stamp = (timestamp or datetime.now()).strftime("%Y%m%d_%H%M%S")
    if recovery and solver_recovery:
        raise ValueError("A run cannot be both full recovery and solver-only recovery")
    if solver_recovery:
        prefix = "musubi_solver_recovery_anchor003274"
    else:
        prefix = "musubi_recovery_anchor003274" if recovery else "musubi_anchor003274"
    root = Path(output_root).resolve() / f"{prefix}_{stamp}"
    if root.exists():
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", f"Run directory already exists: {root}")
    names = ("input", "geometry", "seeder", "musubi", "flow", "qc", "figures", "proteus")
    directories = {name: root / name for name in names}
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=False)
    solver_geometry = directories["geometry"] / "geometry_solver_m"
    solver_geometry.mkdir()
    seeder_mesh = directories["seeder"] / "mesh"
    if not solver_recovery:
        seeder_mesh.mkdir()
    return RunLayout(
        root=root,
        input=directories["input"],
        geometry=directories["geometry"],
        geometry_solver_m=solver_geometry,
        seeder=directories["seeder"],
        seeder_mesh=seeder_mesh,
        musubi=directories["musubi"],
        flow=directories["flow"],
        qc=directories["qc"],
        figures=directories["figures"],
        proteus=directories["proteus"],
    )


def save_input_provenance(layout: RunLayout, inputs: FlowInputs, config_path: Path) -> dict[str, Any]:
    """Copy only small source-of-truth records and save immutable hashes."""

    shutil.copy2(config_path, layout.input / "cfd_flow.yaml")
    shutil.copy2(inputs.boundary_conditions_json, layout.input / "source_boundary_conditions.json")
    shutil.copy2(inputs.boundary_manifest_csv, layout.input / "boundary_manifest.csv")
    reference = {
        "source_run": str(inputs.run_root),
        "surface_geometry_modified": False,
        "files": {
            "tagged_surface_vtp": str(inputs.tagged_surface_vtp),
            "meter_surface_stl": str(inputs.meter_surface_stl),
            "boundary_conditions_json": str(inputs.boundary_conditions_json),
            "boundary_manifest_csv": str(inputs.boundary_manifest_csv),
        },
        "sha256": {
            "tagged_surface_vtp": sha256_file(inputs.tagged_surface_vtp),
            "meter_surface_stl": sha256_file(inputs.meter_surface_stl),
            "boundary_conditions_json": sha256_file(inputs.boundary_conditions_json),
            "boundary_manifest_csv": sha256_file(inputs.boundary_manifest_csv),
        },
    }
    write_json(layout.input / "source_surface_reference.json", reference)
    return reference
