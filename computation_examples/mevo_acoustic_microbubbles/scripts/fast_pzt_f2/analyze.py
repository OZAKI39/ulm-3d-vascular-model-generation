"""Offline F2 evidence aggregation. Reads immutable runs; never edits trajectories."""

import json
import hashlib
import subprocess
import h5py
import numpy as np
from .common import BASE, ROOT, PROTOCOL, write

BRANCHES = {
    name: BASE / "mature_v1" / name
    for name in ["unchanged", "amplitude", "unchanged_repeat", "F1_near_field", "mobility4"]
}
BRANCHES["phase"] = BASE / "phase_v1/phase"
RADIUS = 1.78e-6


def read(path):
    with h5py.File(path / "trajectory.h5") as h:
        d = {k: h[k][:] for k in ["time_s", "positions_m", "particle_id", "radii_m", "control_index"]}
        d.update({k: v[:] for k, v in h["metrics"].items()})
        d["quality"] = h["quality_samples/raw_F1_diagonal"][:, 0]
        d["edge_times"] = h["edge_events/time_s"][:]
        d["edge_ids"] = h["edge_events/batch_i_j_linked"][:]
    return d


def js(path):
    return json.loads(path.read_text())


def adjacency(d):
    x = d["positions_m"][:, 0]
    i, j = np.triu_indices(x.shape[1], 1)
    linked = np.linalg.norm(x[:, i] - x[:, j], axis=-1) <= 2.1 * RADIUS
    union = np.sum(linked | linked[0], axis=1)
    common = np.sum(linked & linked[0], axis=1)
    return dict(
        time_s=d["time_s"],
        pairs=np.column_stack([i, j]),
        linked=linked,
        jaccard_with_initial=np.divide(common, union, out=np.ones_like(common, dtype=float), where=union > 0),
        edge_changes_from_initial=np.sum(linked != linked[0], axis=1),
    )


def edge_dwells(d):
    """Classify raw pair-link dwell, without changing the recorded edge stream."""
    x = d["positions_m"][0, 0]
    graph = np.linalg.norm(x[:, None] - x[None, :], axis=-1) <= 2.1 * RADIUS
    starts = {}
    records = []
    for t, (_, i, j, on) in zip(d["edge_times"], d["edge_ids"]):
        key = (int(i), int(j))
        assert graph[i, j] != bool(on)
        records.append([*key, starts.get(key, 0.0), float(t), int(graph[i, j]), int(key not in starts), 0])
        starts[key] = float(t)
        graph[i, j] = graph[j, i] = bool(on)
    for (i, j), start in starts.items():
        records.append([i, j, start, float(d["time_s"][-1]), int(graph[i, j]), 0, 1])
    a = np.asarray(records, dtype=float).reshape(-1, 7)
    dwell = a[:, 3] - a[:, 2]
    summary = dict(
        raw_flips=len(d["edge_times"]),
        observed_changed_pair_episodes=len(a),
        short_completed_linked_episodes=int(np.sum((a[:, 4] == 1) & (a[:, 6] == 0) & (dwell < 0.002))),
        short_completed_unlinked_episodes=int(np.sum((a[:, 4] == 0) & (a[:, 6] == 0) & (dwell < 0.002))),
        interpretation="Dwell below 2 ms is a brief pair-threshold episode, not a component split/merge or physical fusion/breakup",
    )
    return a, summary


def time_windows(t, mask):
    ids = np.flatnonzero(mask)
    groups = np.split(ids, np.flatnonzero(np.diff(ids) > 1) + 1) if len(ids) else []
    intervals = [
        dict(
            start_s=float(t[g[0]]),
            end_s=float(t[g[-1]]),
            sample_count=len(g),
            observed_span_s=float(t[g[-1]] - t[g[0]]),
        )
        for g in groups
    ]
    return dict(
        sample_count=int(mask.sum()),
        longest_sampled_span_s=max([w["observed_span_s"] for w in intervals], default=0.0),
        intervals=intervals,
    )


def brief(path):
    row = js(path / "SUMMARY.json")
    q = row["quality"]
    d = read(path)
    n = d["positions_m"].shape[2]
    formation = {}
    for fraction in [0.5, 0.9, 1.0]:
        ix = np.flatnonzero(d["largest_size"][:, 0] >= int(np.ceil(n * fraction)))
        formation[str(fraction)] = float(d["time_s"][ix[0]]) if len(ix) else None
    return dict(
        formation_metric_definition="Legacy formation_time_s is first sampled Nmax >= 16; fractional first-hit times are also provided, not maturity",
        first_component_fraction_times_s=formation,
        source=str(path.relative_to(BASE)),
        execution=js(path / "IDENTITY.json"),
        dt_s=row["nominal_dt_s"],
        wall_s=row["wall_s"],
        macro_steps=row["macro_step_count"],
        physical_substeps=row["actual_substep_count"],
        contact_iterations=row["model"]["contact_iterations"],
        quality=row["numerical_quality"],
        max_accepted_overlap_nm=q["max_accepted_penetration_m"] * 1e9,
        max_contact_pass_correction_nm=q["max_contact_pass_correction_m"] * 1e9,
        cumulative_correction_m=q["cumulative_correction_m"],
        mobility_path_m=q["cumulative_physical_displacement_m"],
        free_Stokes_path_m=q["cumulative_free_Stokes_displacement_m"],
        correction_over_mobility_path=q["correction_fraction"],
        correction_over_free_Stokes_path=q["correction_over_free_Stokes_path"],
        max_speed_m_s=q["max_speed_m_s"],
        error_bits=q["error_bits"],
        formation_time_s=row["formation_time_s"],
        final_cluster_count=row["final_cluster_count"],
        final_largest_size=row["final_largest_size"],
    )


def main():
    adoption = js(BASE / "adoption_v1/RESULT.json")
    near = dict(
        status="PASS_IMPLEMENTATION",
        physical_scope="EXPLORATORY",
        comparison=adoption["near_field"],
        comparison_dt_s=1e-5,
        model="Signed pairwise fundamental monopole pressure feedback; local table secant compliance; cubic blend over 1–4 radius sums; harmonics 2–8 and Primary unchanged",
        limitations="Force-only pair closure, not a collective radial solution, multipole translation theory, or quantitative near-field calibration",
        coarse_sensitivity_thresholds=PROTOCOL["sensitivity"],
        interpretation="Coarse formation metrics remain below registered sensitivity bounds, but fixed-ID packing differs by more than two radii; single-seed comparison does not establish physical accuracy",
        runs={n: brief(BASE / "sensitivity_v1" / n) for n in ["f1_10us", "near_10us"]},
        two_bubble_checks=js(BASE / "core_v1/GPU_CHECKS.json"),
        two_bubble_runs={
            n: brief(BASE / "core_v1" / n) for n in ["f1_smooth_near_field", "f2_near_field_candidate"]
        },
    )
    write(BASE / "NEAR_FIELD_COMPARISON.json", near)
    contact = dict(
        status="CONTACT_MODEL_SENSITIVE",
        comparison=adoption["contact"],
        comparison_dt_s=5e-6,
        candidate_disposition=adoption["contact_audit"],
        selected=adoption["selected"],
        runs={
            n: brief(BASE / "sensitivity_v1" / n)
            for n in [
                "near_20us",
                "contact_20us",
                "f1_10us",
                "near_10us",
                "contact_10us",
                "f1_5us",
                "near_5us",
                "contact_5us",
                "secondary_off",
            ]
        },
        interpretation="Both matched 5 us runs satisfy overlap budgets, but eigenvalues differ by 12.45%; candidate angular bias makes it unsuitable for control conclusions. Retain central F1 PBD; final packing is unvalidated.",
    )
    for name in ["near_5us", "contact_5us"]:
        path = BASE / "sensitivity_v1" / name
        a = adjacency(read(path))
        np.savez_compressed(path / "ADJACENCY.npz", **a)
        contact["runs"][name]["final_neighbor_jaccard_to_initial"] = float(a["jaccard_with_initial"][-1])
    write(BASE / "CONTACT_SENSITIVITY.json", contact)
    mature = js(BASE / "mature_v1/RESULT.json")
    phase = js(BASE / "phase_v1/RESULT.json")
    comparisons = dict(mature["comparisons"], phase=phase["comparison"])
    noise = comparisons["unchanged_repeat"]
    thresholds = dict(
        rg_relative=PROTOCOL["observation"]["rg_relative"],
        eigenvalues_relative=PROTOCOL["observation"]["eigenvalues_relative"],
        aligned_internal_rmsd_m=RADIUS * PROTOCOL["observation"]["aligned_rmsd_over_radius"],
        rigid_rotation_deg=PROTOCOL["observation"]["rigid_rotation_deg"],
    )
    effects = {}
    for name in ["phase", "amplitude"]:
        c = comparisons[name]
        detected = {k: bool(c[k] > value and c[k] > 3 * noise[k]) for k, value in thresholds.items()}
        envelope = {k: max(comparisons[b][k] for b in ["F1_near_field", "mobility4"]) for k in thresholds}
        effects[name] = dict(
            endpoint_contrast_to_unchanged=c,
            above_fixed_and_repeat_thresholds=detected,
            exceeds_tested_model_envelope={k: detected[k] and c[k] > envelope[k] for k in thresholds},
            centroid_shift_above_repeat_and_tested_model_envelope=bool(
                c["centroid_distance_m"]
                > 3
                * max(
                    comparisons[b]["centroid_distance_m"]
                    for b in ["unchanged_repeat", "F1_near_field", "mobility4"]
                )
            ),
            status="PASS"
            if any(detected[k] for k in ["rg_relative", "eigenvalues_relative", "aligned_internal_rmsd_m"])
            else "NO_SIGNIFICANT_SHAPE_RESPONSE",
        )
    transient = {}
    repeat_series = np.load(BRANCHES["unchanged_repeat"] / "CONTRAST.npz")
    model_series = [np.load(BRANCHES[k] / "CONTRAST.npz") for k in ["F1_near_field", "mobility4"]]
    for name in ["phase", "amplitude"]:
        c = np.load(BRANCHES[name] / "CONTRAST.npz")
        transient[name] = {}
        for k, limit in thresholds.items():
            numerical = (c[k] > limit) & (c[k] > 3 * repeat_series[k])
            robust = numerical & (c[k] > np.maximum(model_series[0][k], model_series[1][k]))
            transient[name][k] = dict(
                above_fixed_repeat_and_tested_model=time_windows(c["time_s"], robust),
                scope="Descriptive time-resolved screening with unchanged thresholds; contact robustness unavailable; no multiple-comparison inference",
            )
    observed = {}
    for name, path in BRANCHES.items():
        d = read(path)
        a = adjacency(d)
        np.savez_compressed(path / "ADJACENCY.npz", **a)
        m = np.load(path / "MORPHOLOGY.npz")
        initial = m["fixed_centroid_m"][0]
        dwells, dwell_summary = edge_dwells(d)
        np.savez_compressed(
            path / "EDGE_DWELLS.npz",
            columns=["i", "j", "start_s", "end_s", "linked", "left_censored", "right_censored"],
            episodes=dwells,
        )
        observed[name] = dict(
            run=brief(path),
            raw_pair_dwell=dwell_summary,
            absolute_centroid_displacement_m=(m["fixed_centroid_m"][-1] - initial),
            initial_rg_m=float(m["fixed_rg_m"][0]),
            final_rg_m=float(m["fixed_rg_m"][-1]),
            initial_aspect_ratio=float(m["fixed_aspect_ratio"][0]),
            final_aspect_ratio=float(m["fixed_aspect_ratio"][-1]),
            endpoint_kabsch_rotation_deg=float(m["kabsch_rotation_deg"][-1]),
            endpoint_principal_axis_angle_deg=float(m["principal_axis_angle_deg"][-1])
            if np.isfinite(m["principal_axis_angle_deg"][-1])
            else None,
            principal_axis_identifiable_frame_fraction=float(m["principal_axis_identifiable"].mean()),
            endpoint_internal_rmsd_m=float(m["aligned_internal_rmsd_m"][-1]),
            any_membership_change=bool(np.any(m["dynamic_member_mask"] != m["dynamic_member_mask"][0])),
            final_neighbor_jaccard_to_mature=float(a["jaccard_with_initial"][-1]),
            sampled_neighbor_change_count=int(np.sum(a["linked"][1:] != a["linked"][:-1])),
            topology=js(path / "TOPOLOGY_PERSISTENCE.json"),
        )
    pzt = dict(
        status="NO_SIGNIFICANT_SHAPE_RESPONSE"
        if all(v["status"] == "NO_SIGNIFICANT_SHAPE_RESPONSE" for v in effects.values())
        else "PASS",
        endpoint_status_scope="120 ms continuation endpoint; does not describe the whole time series",
        transient_screening=transient,
        transient_interpretation="Short intervals of eigenvalue contrast exceed the same-repeat and tested force/mobility envelopes. ID-aligned internal contrast never exceeds the tested model envelope. These are descriptive transient responses, not validated lasting shape control.",
        maturity=mature["maturity"],
        continuation_s=PROTOCOL["continuation_s"],
        effects=effects,
        endpoint_unchanged_repeat=noise,
        endpoint_model_comparisons={k: comparisons[k] for k in ["F1_near_field", "mobility4"]},
        observation_criteria=PROTOCOL["observation"],
        branches=observed,
        controls=dict(
            unchanged=mature["field_selection"]["baseline"],
            amplitude=mature["field_selection"]["selected"]["amplitude"],
            phase=phase["field_selection"],
        ),
        field_diagnostics=["mature_v1/FIELD_DIAGNOSTICS.npz", "phase_v1/FIELD_DIAGNOSTICS.npz"],
        interpretation="Control produces resolved centroid translation. Endpoint shape, orientation and ID-level rearrangement do not exceed the registered repeat envelope. Large best-fit rotations coexist with deformation and must not be interpreted as controlled rigid rotation. No membership change or raw component split/merge in mature branches.",
        limits="One seed and one same-control repeat: descriptive sensitivity, not a statistical confidence interval. Contact robustness is not established; rejected contact candidate omitted from mature controls. Near-isotropic principal axes are masked. Negative endpoint result does not exclude transient responses.",
    )
    write(BASE / "PZT_MORPHODYNAMICS.json", pzt)
    write(BASE / "PERFORMANCE.json", js(BASE / "performance_v1/RESULT.json"))
    # Join only the observed common prefix and its actual continuation, with exact equality at restart.
    d = read(BASE / "mature_v1/formation")
    e = read(BRANCHES["amplitude"])
    idx = mature["maturity"]["frame_index"]
    np.testing.assert_array_equal(d["positions_m"][idx], e["positions_m"][0])
    t = np.concatenate([d["time_s"][: idx + 1], d["time_s"][idx] + e["time_s"][1:]])
    np.testing.assert_allclose(np.diff(t), 0.0002, rtol=0, atol=1e-14)
    vis = BASE / "visualization"
    vis.mkdir(exist_ok=True)
    with h5py.File(vis / "continuous_actual_states.h5", "w") as h:
        h["time_s"] = t
        for key in ["positions_m", "labels", "cluster_count", "largest_size", "global_rg_m"]:
            h.create_dataset(key, data=np.concatenate([d[key][: idx + 1], e[key][1:]]), compression="gzip")
        h["particle_id"] = d["particle_id"]
        h["radii_m"] = d["radii_m"]
        h["control_index"] = (np.arange(len(t)) >= idx).astype(int)
        h["source_frame_index"] = np.concatenate([np.arange(idx + 1), np.arange(1, len(e["time_s"]))])
        h["source_run_index"] = np.concatenate(
            [np.zeros(idx + 1, dtype=int), np.ones(len(e["time_s"]) - 1, dtype=int)]
        )
        h.attrs["sources"] = json.dumps(
            ["mature_v1/formation/trajectory.h5", "mature_v1/amplitude/trajectory.h5"]
        )
        h.attrs["position_rule"] = (
            "Exact shared-prefix branch, equal restart coordinates, no interpolation, no synthetic states"
        )
        h.attrs["switch_time_s"] = float(t[idx])
    write(
        vis / "STATE_PROVENANCE.json",
        dict(
            source_prefix="mature_v1/formation/trajectory.h5",
            prefix_frames=[0, idx],
            source_continuation="mature_v1/amplitude/trajectory.h5",
            continuation_frames=[1, len(e["time_s"]) - 1],
            restart_bitwise_equal=True,
            actual_radius_m=RADIUS,
            frame_count=len(t),
            time_step_s=0.0002,
            absolute_switch_time_s=float(t[idx]),
            terminal_time_s=float(t[-1]),
            initial_control=mature["field_selection"]["baseline"],
            changed_control=mature["field_selection"]["selected"]["amplitude"]["event"],
        ),
    )
    # Verify historical simulator code remains byte-identical to the exact F1 archive.
    old = "6de76d508ea088e137069abf058b86b7f6d47a1c"
    kept = {}
    for name in ["__init__.py", "engine.py", "kernels.py", "models.py", "metrics.py", "storage.py"]:
        p = ROOT / "src/mevo_mb/fast_sim" / name
        expected = subprocess.check_output(
            [
                "git",
                "show",
                f"{old}:computation_examples/mevo_acoustic_microbubbles/src/mevo_mb/fast_sim/{name}",
            ],
            cwd=ROOT,
        )
        assert expected == p.read_bytes()
        kept[name] = hashlib.sha256(expected).hexdigest()
    write(
        BASE / "F1_PRESERVATION.json",
        dict(
            archive_commit=old,
            historical_execution_commits=[
                "9829156893339efff713483f3125848b84adbc98",
                "3721bebd597164262691af766f723d245e043116",
            ],
            unchanged_original_simulator_sha256=kept,
        ),
    )
    print(
        json.dumps(
            dict(
                status=pzt["status"], maturity_ms=mature["maturity"]["absolute_time_s"] * 1e3, effects=effects
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
