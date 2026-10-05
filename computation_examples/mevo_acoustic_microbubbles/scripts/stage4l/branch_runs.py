"""New-model branch runs from immutable legacy inputs, with explicit S0 gating."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import time

import numpy as np

from scripts.stage4l.dense_probe import BASE, sha, write


def run_case(legacy, case, dt_factor=1.0, secondary=False, duration=6.0, max_wall_s=1800):
    from mevo_mb.acoustics.basis import ComplexAcousticField
    from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
    from mevo_mb.bubbles.population import BubblePopulation
    from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
    from mevo_mb.flow.fields import VTKFlowField
    from mevo_mb.io.artifacts import new_run
    from scripts.stage4.synthetic import AcousticOff
    import mevo_mb.bubbles.tracker as tracking

    baseline = json.loads((legacy / "results/stage4/primary_only/STATE.json").read_text())
    row = next(row for row in baseline["rows"] if row["case"] == case and row["seed"] == 20261002)
    source = legacy / row["run"]
    previous = json.loads((source / "RUN_IDENTITY.json").read_text())
    config = copy.deepcopy(previous["config"])
    inputs = {str(source / "RUN_IDENTITY.json"): sha(source / "RUN_IDENTITY.json")}

    def frozen(path):
        path = legacy / str(path).removeprefix("/workspace/mevo_acoustic_microbubbles/")
        if not path.is_file():
            raise FileNotFoundError(path)
        inputs[str(path)] = sha(path)
        return str(path)

    config["geometry"]["mesh_surface"] = frozen(config["geometry"]["mesh_surface"])
    config["bubble_response"]["path"] = frozen(config["bubble_response"]["path"])
    f = config["flow"]
    f["snapshot"] = frozen(f["snapshot"])
    flow = VTKFlowField(
        f["snapshot"],
        velocity_array=f["velocity_array"],
        coordinate_scale_to_m=f["coordinate_scale_to_m"],
        velocity_scale_to_m_s=f["velocity_scale_to_m_s"],
        spatial_method=f["spatial_method"],
    )
    if case == "B0":
        field = AcousticOff()
    else:
        config["effective_acoustic_field"]["path"] = frozen(config["effective_acoustic_field"]["path"])
        field = ComplexAcousticField.load(config["effective_acoustic_field"]["path"])
    initial = Path(frozen(source / "INITIAL_POPULATION.npz"))
    with np.load(initial) as p:
        population = BubblePopulation(p["positions_m"].copy(), p["radii_m"].copy(), p["particle_id"].copy())
    config["particles"].update(
        translation_backend="lubrication_coupled",
        secondary_enabled=secondary,
        duration_s=duration,
        dt_s=config["particles"]["dt_s"] * dt_factor,
        pair_benchmark_path=frozen("results/stage4/performance_refactor/protocol/PAIR_BACKEND_STATE.json"),
        lubrication=dict(activation=0.1, floor=0.001, residual_tolerance=1e-11),
    )
    kind = (
        ("secondary_full" if duration == 6 else "secondary_short")
        if secondary
        else ("primary_baseline" if dt_factor == 1 else "timestep_sensitivity")
    )
    out = new_run(BASE / kind, f"{case}_dtfactor{dt_factor:g}", config)
    write(
        out / "RUN_IDENTITY.json",
        dict(
            config=config,
            legacy_reference=row["run"],
            initial_population_sha256=previous["initial_population_sha256"],
            inputs_sha256=inputs,
            translation_backend="lubrication_coupled",
            secondary_enabled=secondary,
            case=case,
            dt_factor=dt_factor,
        ),
    )
    np.savez_compressed(
        out / "INITIAL_POPULATION.npz",
        positions_m=population.positions_m,
        radii_m=population.radii_m,
        particle_id=population.ids,
    )
    tracker = GlobalSynchronousIntegrator(config, field, flow)
    original = tracking.project_contacts
    sweep = OrderedContactSweep(BASE / "unit_tests/ordered_contact_cache")
    tracking.project_contacts = sweep.project
    start = time.perf_counter()
    step = tracker.step

    def finite_step(*args, **kwargs):
        if time.perf_counter() - start > max_wall_s:
            raise RuntimeError("DEVELOPMENT_RUNTIME_BUDGET_REACHED")
        return step(*args, **kwargs)

    tracker.step = finite_step
    try:
        result = tracker.run(population, out)
    finally:
        tracking.project_contacts = original
        write(out / "INPUTS_POSTCHECK.json", dict(unchanged=all(sha(p) == v for p, v in inputs.items())))
    write(
        out / "SUMMARY.json", dict(**result, run=str(out), legacy_reference=row["run"], dt_factor=dt_factor)
    )
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--case", choices=["B0", "B1", "B2", "B3"], required=True)
    parser.add_argument("--dt-factor", type=float, choices=[1.0, 0.5], default=1.0)
    parser.add_argument("--secondary", action="store_true")
    parser.add_argument("--duration", type=float, choices=[0.1, 0.5, 1.0, 6.0], default=6.0)
    parser.add_argument("--gate", type=Path, required=True)
    args = parser.parse_args()
    gate = json.loads(args.gate.read_text())
    if gate.get("dense_probe_gate") != "PASS_WITH_SCOPE":
        raise ValueError("S0 dense probe gate has not passed")
    if args.secondary and args.duration in (0.5, 1.0) and gate.get("secondary_0p1_gate") != "PASS_WITH_SCOPE":
        raise ValueError("Secondary 0.1 s gate has not passed")
    if args.secondary and args.duration == 6.0:
        required = [
            "secondary_0p5_or_1_gate",
            "full_secondary_runtime_gate",
            "primary_baseline_gate",
            "timestep_gate",
        ]
        if any(gate.get(key) != "PASS_WITH_SCOPE" for key in required):
            raise ValueError("Full secondary staged gates have not all passed")
    run_case(args.legacy_root, args.case, args.dt_factor, args.secondary, args.duration)


if __name__ == "__main__":
    main()
