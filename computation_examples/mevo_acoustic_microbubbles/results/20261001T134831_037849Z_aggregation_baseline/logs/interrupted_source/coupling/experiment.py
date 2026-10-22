"""One-way field/response/transport coupling, independent of the flow implementation."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import yaml

from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.acoustics.kwave_backend import acoustic_identity
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.control.pzt import PZTArray
from mevo_mb.flow.fields import PoiseuilleCylinderFlowField, SimVascularFlowField, ZeroFlowField
from mevo_mb.io.artifacts import git_commit, write_json


def flow_from_config(config):
    f = config["flow"]
    if f["model"] == "zero":
        return ZeroFlowField()
    if f["model"] == "poiseuille":
        return PoiseuilleCylinderFlowField(config["vessel"]["radius_m"], f["maximum_velocity_m_s"])
    if f["model"] == "simvascular":
        return SimVascularFlowField(f["file"], f["velocity_array"], f["coordinate_scale_to_m"], f["velocity_scale_to_m_s"])
    raise ValueError("Unknown flow provider")


def run_particles(config: dict, basis_path: Path, run: Path, pattern: str = "baseline", field_override=None) -> dict:
    db = AcousticBasisDatabase.load_basis(basis_path)
    if acoustic_identity(config) != db.metadata["acoustic_identity"]:
        raise ValueError("Basis physics/grid/frequency identity differs from requested config")
    population = BubblePopulation.initialize(config)
    array = PZTArray.from_config(config)
    array.set_phases(config["controls"]["phase_patterns_rad"][pattern])
    u = array.get_control_vector()
    field = db.combine(u)
    scale = 1.
    if "aggregation" in config:
        initial = np.abs(field.interpolate_complex_pressure(population.positions_m))
        scale = config["aggregation"]["target_initial_local_PNP_Pa"] / np.median(initial)
        u *= scale
        field = db.combine(u)
    if field_override is not None:
        field = field_override
    local_initial = np.abs(field.interpolate_complex_pressure(population.positions_m))
    write_json(run/"acoustic/CONTROL.json", {"pattern": pattern, "source_control_real_Pa": u.real, "source_control_imag_Pa": u.imag,
               "source_amplitudes_Pa": np.abs(u), "aggregation_source_scale": scale,
               "initial_local_PNP_min_median_max_Pa": [local_initial.min(), np.median(local_initial), local_initial.max()],
               "note": "Source amplitude is not local pressure. Aggregation targets median initial cloud PNP only; no response formula adjusted.",
               "basis_path": str(basis_path), "basis_sha256": hashlib.file_digest(basis_path.open("rb"), "sha256").hexdigest()})
    resolved = {"control_pattern": pattern, "source_amplitudes_Pa": np.abs(u).tolist(), "source_phases_rad": np.angle(u).tolist()}
    (run/"acoustic/resolved_controls.yaml").write_text(yaml.safe_dump(resolved))
    summary = ParticleTracker(config, field, flow_from_config(config)).run(population, run)
    summary.update(pattern=pattern, basis=str(basis_path), git_commit=git_commit(), flow_model=config["flow"]["model"])
    write_json(run/"metadata.json", summary)
    from mevo_mb.visualization.plots import plot_experiment
    from mevo_mb.io.vtk_export import export_run
    plot_experiment(run, field, config)
    export_run(run, config["output"]["vtk_stride"], field)
    return summary
