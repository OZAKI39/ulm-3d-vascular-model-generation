"""One-way field/response/transport coupling, independent of the flow implementation."""

from __future__ import annotations

import hashlib
import time
from pathlib import Path

import numpy as np
import yaml

from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.acoustics.kwave_backend import acoustic_identity
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.control.pzt import PZTArray
from mevo_mb.flow.fields import PoiseuilleCylinderFlowField, SimVascularFlowField, VTKFlowField, ZeroFlowField
from mevo_mb.io.artifacts import git_commit, write_json


def flow_from_config(config):
    f = config["flow"]
    backend = f.get("backend", f.get("model"))
    if backend == "zero":
        return ZeroFlowField()
    if backend == "poiseuille":
        return PoiseuilleCylinderFlowField(config["vessel"]["radius_m"], f["maximum_velocity_m_s"])
    if backend == "simvascular":
        return SimVascularFlowField(
            f["file"], f["velocity_array"], f["coordinate_scale_to_m"], f["velocity_scale_to_m_s"]
        )
    if backend == "vtk":
        return VTKFlowField(
            f["files"],
            f.get("times_s"),
            velocity_array=f.get("velocity_array", "Velocity"),
            coordinate_scale_to_m=f.get("coordinate_scale_to_m", 1.0),
            velocity_scale_to_m_s=f.get("velocity_scale_to_m_s", 1.0),
        )
    raise ValueError("Unknown flow provider")


def run_particles(
    config: dict, basis_path: Path, run: Path, pattern: str = "baseline", field_override=None
) -> dict:
    gate_metadata = None
    if config['bubble_response']['model'] != 'linear_surrogate' or 'validation_gates' in config:
        from mevo_mb.validation.gates import ValidationGateManager, ValidationGateError
        gate = config.get('validation_gates', {})
        try:
            gate_metadata = ValidationGateManager(gate.get('registry', 'results/stage3/audit/PHYSICS_GATES.json')).check(
                gate.get('command', 'run-cylinder'), scope=gate.get('scope', 'sonovue_cylinder'),
                force_demo=gate.get('force_demo', False))
        except ValidationGateError as error:
            write_json(run/'GATE_DECISION.json', error.metadata)
            raise
        write_json(run/'GATE_DECISION.json', gate_metadata)
    basis_start = time.perf_counter()
    db = AcousticBasisDatabase.load_basis(basis_path)
    basis_elapsed = time.perf_counter()-basis_start
    if acoustic_identity(config) != db.metadata["acoustic_identity"]:
        raise ValueError("Basis physics/grid/frequency identity differs from requested config")
    population = BubblePopulation.initialize(config)
    array = PZTArray.from_config(config)
    array.set_phases(config["controls"]["phase_patterns_rad"][pattern])
    u = array.get_control_vector()
    field = db.combine(u)
    scale = 1.0
    if "aggregation" in config:
        initial = np.abs(field.interpolate_complex_pressure(population.positions_m))
        scale = config["aggregation"]["target_initial_local_PNP_Pa"] / np.median(initial)
        u *= scale
        field = db.combine(u)
    if field_override is not None:
        field = field_override
    local_initial = np.abs(field.interpolate_complex_pressure(population.positions_m))
    field.save(run / "acoustic/combined_field.h5")
    write_json(
        run / "acoustic/CONTROL.json",
        {
            "pattern": pattern,
            "source_control_real_Pa": u.real,
            "source_control_imag_Pa": u.imag,
            "source_amplitudes_Pa": np.abs(u),
            "aggregation_source_scale": scale,
            "initial_local_PNP_min_median_max_Pa": [
                local_initial.min(),
                np.median(local_initial),
                local_initial.max(),
            ],
            "note": "Source amplitude is not local pressure. Aggregation targets median initial cloud PNP only; no response formula adjusted.",
            "basis_path": str(basis_path),
            "field_override": field_override is not None,
            "exact_combined_field": "acoustic/combined_field.h5",
            "basis_sha256": hashlib.file_digest(basis_path.open("rb"), "sha256").hexdigest(),
        },
    )
    resolved = {
        "control_pattern": pattern,
        "source_amplitudes_Pa": np.abs(u).tolist(),
        "source_phases_rad": np.angle(u).tolist(),
    }
    (run / "acoustic/resolved_controls.yaml").write_text(yaml.safe_dump(resolved))
    summary = ParticleTracker(config, field, flow_from_config(config)).run(population, run)
    summary.update(
        pattern=pattern, basis=str(basis_path), git_commit=git_commit(), flow_model=config["flow"].get("backend", config["flow"].get("model"))
    )
    summary.setdefault("timings", {})["acoustic_basis_loading_s"] = basis_elapsed
    if gate_metadata:
        summary.update(gate_metadata)
    write_json(run / "metadata.json", summary)
    from mevo_mb.visualization.plots import plot_experiment
    from mevo_mb.io.vtk_export import export_run

    plot_experiment(run, field, config)
    export_run(run, config["output"]["vtk_stride"], field)
    return summary
