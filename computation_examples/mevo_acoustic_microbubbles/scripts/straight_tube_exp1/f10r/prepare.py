"""Pre-register finite experiments before any F10R force/trajectory calculation."""
import copy
from datetime import datetime, timezone
import numpy as np
from mevo_mb.control.scheduled_pzt import Segment, Schedule, one_hot
from .common import BASE, OLD, Q, CONFIG, ROOT, original, policy, write, read, sha, context


def main():
    assert not CONFIG.exists() and not (BASE / "PREREGISTRATION.json").exists()
    p = dict(public_baseline="467b17411deeb13a77370565b7382f071b41a774",
             complete_F10Q_source="41ae8ed2b88f2cbbed834fffccaa0b5de6e91589",
             labels=["EXPLORATORY", "HOMOGENEOUS_WATER_CONTROL_PROXY", "PHYSICS_UNCALIBRATED", "ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED"],
             training_ready=False, channel_order=["Y-", "Y+", "X-", "X+"],
             total_scientific_budget_s=600, remote_budget_s=480, per_run_budget_s=45,
             max_trajectories=32, max_full_navigation_groups=3, frequency_Hz=2.5e6,
             source_range_Pa_equivalent=[0, 1500], flow="STRICT_ZERO",
             targets_m={"X": [2e-7, 0, 0], "Y": [0, 2e-7, 0]}, target_tolerance_m=7.5e-8,
             target_definition="Fixed absolute world positions, not shifted with the measured result",
             dt_s=1e-6, dt_half_s=0.5e-6, duration_s=0.002, short_s=0.0005,
             phase_perturbation_deg=[-10, -6, -3, -1, 0, 1, 3, 6, 10],
             phase_full_trajectory_deg=[-6, 6, 10],
             phase_prediction="User hypothesis: about 6% change per degree in center primary Fx; beyond roughly +/-6 degrees target tolerance may fail. No claim before measured results.",
             phase_changed_channel="X+ only; every other full-precision F10Q parameter frozen",
             original_X_control=original("X"), original_Y_control=original("Y"),
             initial_states=["original", "rotate_z_plus5", "rotate_z_minus5"],
             perturbation_scope="Two deterministic rigid rotations about unchanged COM; not independent experimental samples; no projection",
             pair_static_search=dict(x_delta_deg=[150, 160, 170, 175, 185, 190, 200, 210], amplitude_Pa=[500, 1000, 1500], max_initial_pressure_Pa=1200),
             pair_peak_phases_deg=list(range(-60, 61, 10)), pair_peak_scan_halfwidth_m=0.0003,
             pair_peak_scope="Field-only pressure maxima on axis; potential traps, not proven 3D force equilibria or long-time transport",
             Y_static_candidates=["Y_ONE_HOT", "Y_PAIR_90", "Y_FOUR_X025", "Y_FOUR_X050"],
             Y_screen_gate="Complete .5ms, positive +Y displacement, |Rg change| <=5%, one 2.1R component; select smallest Rg change among candidates with >=1nm projected 2ms motion, otherwise retain failure/weak-control result",
             Y_scope_rule="If near-field zero/efficiency loss confirmed, all Y statements HOMOGENEOUS_WATER_PROXY_SPECIFIC",
             initial_legality="Original sphere nonpenetration and cavity queries; reject invalid state without projection",
             admission_policy=policy(),
             switches="IDEAL_SCHEDULED_SEGMENTS; exact shortened step at nonaligned boundary; global coherent phase clock; no measured USB latency/quantization",
             hardware_executed=False, voltage_conversion=None,
             sequence_stop="Compatibility before all new trajectories; short screens before full candidates; dt/cache repeats for admitted candidate groups only; stop on quality or wall-budget failure")
    write(CONFIG, p)
    c = context()
    x = c["x0"]
    states = {"original": x}
    for name, degree in [("rotate_z_plus5", 5), ("rotate_z_minus5", -5)]:
        a = np.deg2rad(degree)
        R = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        states[name] = (x - x.mean(0)) @ R.T + x.mean(0)
    for n, row in states.items():
        assert all(c["cavity"].query(v, p["admission_policy"]["radius_m"])["valid"] for v in row), n
        i, j = np.triu_indices(32, 1)
        assert np.linalg.norm(row[i] - row[j], axis=1).min() >= 2 * p["admission_policy"]["radius_m"] - 1e-15, n
    np.savez_compressed(BASE / "INITIAL_STATES.npz", **states)
    plans = {}
    plans["X_ORIGINAL"] = Schedule((Segment.from_legacy(original("X")),)).record()
    off = dict(original("X"), amplitude_Pa=[0.0]*4, phase_rad=[0.0]*4, enabled=[False]*4)
    plans["OFF"] = Schedule((Segment.from_legacy(off),)).record()
    pair = []
    for i, shift in enumerate([-3, 3]):
        e = copy.deepcopy(original("X")); e["phase_rad"][3] += np.deg2rad(shift)
        pair.append(Segment.from_legacy(e, segment_id=i, start_s=i*0.001, duration_s=0.001))
    plans["X_SEGMENTED"] = Schedule(tuple(pair)).record()
    for shift in p["phase_full_trajectory_deg"]:
        e = copy.deepcopy(original("X")); e["phase_rad"][3] += np.deg2rad(shift)
        plans[f"X_PHASE_{shift:+d}"] = Schedule((Segment.from_legacy(e),)).record()
    plans["ONE_HOT_X_SEQUENCE"] = Schedule((one_hot("X-", 400, duration_s=.001),
                                             one_hot("X+", 400, duration_s=.001, start_s=.001, segment_id=1))).record()
    write(BASE / "REGISTERED_SCHEDULES.json", plans)
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in [CONFIG, BASE/"INITIAL_STATES.npz", BASE/"REGISTERED_SCHEDULES.json",
        Q/"NAVIGATION_SELECTION.json", Q/"b_navigation_v1/runs/NAV_X_DT1/INPUT.json", Q/"b_navigation_v1/runs/NAV_X_DT1/trajectory.h5",
        Q/"b_navigation_v1/runs/NAV_Y_DT1/trajectory.h5", Q/"b_navigation_v1/runs/NAV_Y_DT05/trajectory.h5", Q/"b_navigation_v1/runs/NAV_Y_FINE/trajectory.h5"]}
    inputs.update({str((OLD/n).relative_to(ROOT)): sha(OLD/n) for n in read(Q/"PREREGISTRATION.json")["input_sha256"]})
    write(BASE / "PREREGISTRATION.json", dict(registered_utc=datetime.now(timezone.utc).isoformat(),
          protocol=p, input_sha256=inputs, force_diagnostics_before_Y_trajectories=True,
          predictions_are_hypotheses_not_results=True, source_scope="Full source controls imported without rounding"))
    print("Registered fixed targets, 3 initial states and finite phase/candidate budgets")


if __name__ == "__main__":
    main()
