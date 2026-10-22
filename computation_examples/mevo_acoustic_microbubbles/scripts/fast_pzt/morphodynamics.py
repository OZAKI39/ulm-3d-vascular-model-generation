"""One continuous real-basis trajectory, with a same-state fixed-control branch."""

import argparse
import json
import numpy as np
import h5py
from mevo_mb.fast_sim.models import Physics, ResponseSection, Basis, pair_forces_numpy
from mevo_mb.bubbles.forces import primary_bjerknes
from .common import BASE, inputs, write
from .experiment import run_case


def main(label, performance):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    chosen = json.loads((BASE / performance / "SELECTED_STEP.json").read_text())
    p = Physics(contact_iterations=chosen["contact_iterations"])
    basis = Basis(inputs()[1], p.frequency_Hz)
    response = ResponseSection(inputs()[0], p)
    center = np.array([0, 80e-6, 40e-6])
    controls = []
    for t, amp, phase, on in [
        (0, [1, 1, 1, 1], [0, 0, 0, 0], [True] * 4),
        (0.08, [1, 0.6, 1.2, 1], [0, np.pi / 2, np.pi / 3, 0], [True, True, True, False]),
    ]:
        field, _, _ = basis.combine(amp, phase, on)
        scale = 850 / abs(field.interpolate_complex_pressure(center)[0])
        controls.append(
            dict(time_s=t, amplitude_Pa=(scale * np.array(amp)).tolist(), phase_rad=phase, enabled=on)
        )
    cfg = dict(
        count=32,
        seed=20261008,
        cloud_radius_m=30e-6,
        center_m=center.tolist(),
        basis=True,
        duration_s=0.20,
        dt_s=chosen["dt_s"],
        sample_interval_s=0.0002,
        physics={"contact_iterations": chosen["contact_iterations"]},
        events=controls,
        control_response_criteria={
            "centroid_difference_over_radius": 0.01,
            "shape_tensor_relative_difference": 0.005,
            "scope": "Observable response in this exploratory numerical model, not calibrated physical accuracy",
        },
    )
    mainrun = run_case(out / "controlled", cfg)
    with h5py.File(out / "controlled/trajectory.h5") as h:
        times = h["time_s"][:]
        idx = int(np.argmin(abs(times - 0.08)))
        assert abs(times[idx] - 0.08) < 1e-14
        switch_x = h["positions_m"][idx]
        controlled_x = h["positions_m"][idx:]
    fixed_cfg = dict(
        cfg,
        positions_m=switch_x.tolist(),
        duration_s=0.12,
        events=[dict(controls[0], time_s=0)],
        absolute_time_offset_s=0.08,
        initial_condition="Actual continuous controlled-run state at 80 ms",
    )
    fixed = run_case(out / "fixed_from_switch", fixed_cfg)
    forces = []
    for event in controls:
        field, _, _ = basis.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])
        pressure = field.interpolate_complex_pressure(switch_x[0])
        harmonics = response.query(pressure)
        primary = primary_bjerknes(harmonics[:, 0], field.interpolate_complex_gradient(switch_x[0]))
        secondary = pair_forces_numpy(switch_x[0], harmonics, p)
        forces.append(dict(pressure_Pa=pressure, primary_N=primary, secondary_N=secondary))
    np.savez_compressed(
        out / "CONTROL_SWITCH_FORCES.npz",
        positions_m=switch_x,
        **{f"{key}_{side}": value for side, record in enumerate(forces) for key, value in record.items()},
    )
    with h5py.File(out / "fixed_from_switch/trajectory.h5") as h:
        fixed_x = h["positions_m"][:]
        fixed_times = h["time_s"][:]
    comparison = None
    if mainrun["status"] == "COMPLETE" and fixed["status"] == "COMPLETE":
        np.testing.assert_allclose(times[idx:], 0.08 + fixed_times, rtol=0, atol=2e-14)
        np.testing.assert_array_equal(controlled_x[0], fixed_x[0])
        delta = controlled_x[-1, 0].mean(axis=0) - fixed_x[-1, 0].mean(axis=0)
        tensors = []
        for x in [controlled_x[-1, 0], fixed_x[-1, 0]]:
            xx = x - x.mean(axis=0)
            tensors.append(xx.T @ xx / len(x))
        shape = float(np.linalg.norm(tensors[0] - tensors[1]) / np.linalg.norm(tensors[1]))
        eigen = float(
            np.linalg.norm(np.linalg.eigvalsh(tensors[0]) - np.linalg.eigvalsh(tensors[1]))
            / np.linalg.norm(np.linalg.eigvalsh(tensors[1]))
        )
        comparison = dict(
            final_centroid_delta_m=delta.tolist(),
            centroid_distance_m=float(np.linalg.norm(delta)),
            shape_tensor_relative_difference=shape,
            eigenvalue_relative_difference=eigen,
            same_switch_state=True,
            no_interpolation=True,
        )
    force_before = forces[0]["primary_N"] + forces[0]["secondary_N"]
    force_after = forces[1]["primary_N"] + forces[1]["secondary_N"]
    force_difference = float(
        np.linalg.norm(force_after - force_before) / max(np.linalg.norm(force_before), 1e-300)
    )
    outcome = "NOT_EVALUABLE"
    if comparison is not None:
        outcome = (
            "PASS"
            if mainrun["numerical_quality"] == "PASS"
            and fixed["numerical_quality"] == "PASS"
            and comparison["centroid_distance_m"] >= 0.01 * p.radius_m
            and comparison["shape_tensor_relative_difference"] >= 0.005
            and force_difference > 1e-6
            else "FAIL"
        )
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            controlled=mainrun,
            fixed_from_actual_switch=fixed,
            comparison=comparison,
            relative_force_change_at_identical_state=force_difference,
            pzt_morphodynamics=outcome,
            scope="EXPLORATORY; same-state counterfactual, no claims of physical calibration",
        ),
    )


if __name__ == "__main__":
    a = argparse.ArgumentParser()
    a.add_argument("--label", required=True)
    a.add_argument("--performance", required=True)
    q = a.parse_args()
    main(q.label, q.performance)
