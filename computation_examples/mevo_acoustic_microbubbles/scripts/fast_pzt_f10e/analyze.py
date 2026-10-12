"""Offline analysis of retained original attempts, accepted prefixes and contact traces."""

import gzip
import json
import hashlib
import h5py
import numpy as np
from mevo_mb.fast_sim.nearfield import F2Physics, coupled_fundamental, feedback_weight, pair_reference
from mevo_mb.fast_sim.morphodynamics import rigid_fit
from .common import BASE, ROOT, HISTORY, CAMPAIGN, js, sha, write, identity


def trace(path):
    with gzip.open(path / "CONTACT_TRACE.json.gz", "rt") as stream:
        return json.load(stream)


def overlaps(x, radius):
    d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    return np.triu(np.maximum(2 * radius - d, 0), 1)


def contact_report(path):
    data = trace(path)
    if not data["traces"]:
        return dict(
            path=str(path.relative_to(BASE)), audited_steps=0, first_rejection=data["first_rejection_step"]
        )
    t = data["traces"][-1]
    radius = data["radius_m"]
    prediction = np.array(t["predicted_positions_m"])
    previous = overlaps(prediction, radius)
    rows = []
    chain = []
    cumulative = np.zeros(len(prediction))
    for r in t["rounds"]:
        after = overlaps(np.array(r["positions_m"]), radius)
        delta = after - previous
        increased = np.argwhere(delta > 1e-12)
        decreased = np.argwhere(delta < -1e-12)
        correction = np.linalg.norm(r["correction_m"], axis=1)
        cumulative += correction
        for i, j in increased:
            matching = [(int(a), int(b)) for a, b in decreased if len({int(i), int(j), int(a), int(b)}) == 3]
            if matching and len(chain) < 12:
                chain.append(
                    dict(
                        round=r["round"],
                        increased_pair=[int(i), int(j)],
                        overlap_before_nm=float(previous[i, j] * 1e9),
                        overlap_after_nm=float(after[i, j] * 1e9),
                        neighboring_decreased_pair=list(matching[0]),
                        interpretation="Shared-particle overlap redistribution in the same simultaneous Jacobi pass; not a serial causal ordering",
                    )
                )
        pair = np.unravel_index(np.argmax(after), after.shape)
        rows.append(
            dict(
                round=r["round"],
                max_residual_nm=r["max_overlap_m"] * 1e9,
                max_single_pass_nm=r["max_correction_m"] * 1e9,
                max_pair=list(pair),
                increases=len(increased),
                decreases=len(decreased),
                max_increase_nm=max(0, float(delta.max() * 1e9)),
            )
        )
        previous = after
    x = np.array(t["initial_positions_m"])
    p = F2Physics()
    rawh = np.array(t["harmonics_complex_volume_m3"])
    harmonics = rawh[..., 0] + 1j * rawh[..., 1]
    rawp = np.array(t["pressure_complex_Pa"])
    pressure = rawp[:, 0] + 1j * rawp[:, 1]
    chi = complex(*data["chi_zero"])
    pair_force = np.array(t["pair_force_N"])
    secondary = np.array(t["secondary_force_N"])
    reference = pair_reference(x, harmonics, pressure, p, chi)
    pair_rows = []
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            d = float(np.linalg.norm(x[j] - x[i]))
            ui, uj, pi, pj, det = coupled_fundamental(
                harmonics[i, 0], harmonics[j, 0], pressure[i], pressure[j], d, p, chi
            )
            pair_rows.append(
                dict(
                    pair=[i, j],
                    distance_over_sum=d / (2 * radius),
                    gap_nm=(d - 2 * radius) * 1e9,
                    feedback_weight=float(feedback_weight(d, radius)),
                    determinant_abs=abs(det),
                    coupled_pressure_abs_Pa=[abs(pi), abs(pj)],
                    force_i_due_to_j_N=pair_force[i, j],
                    force_norm_N=float(np.linalg.norm(pair_force[i, j])),
                    radial_attraction_N=float(np.dot(pair_force[i, j], (x[j] - x[i]) / d)),
                )
            )
    # Residual of the coupled normal-resistance equations, evaluated with actual
    # local velocity. This diagnostic does NOT perform extra solver iterations.
    free = np.array(t["free_velocity_m_s"])
    local = np.array(t["local_velocity_m_s"])
    matrix = np.repeat(np.eye(3)[None], len(x), axis=0)
    coupling = np.zeros((len(x), len(x), 3, 3))
    activation_pairs = 0
    floor_pairs = 0

    def phi(e):
        return 0.5 / e - 0.45 * np.log(e) - (3 / 56) * e * np.log(e)

    for i in range(len(x)):
        for j in range(len(x)):
            if i == j:
                continue
            delta = x[j] - x[i]
            d = np.linalg.norm(delta)
            e = (d - 2 * radius) / radius
            if e < p.lubrication_activation:
                n = delta / d
                c = 0.5 * (phi(max(e, p.lubrication_floor)) - phi(p.lubrication_activation))
                coupling[i, j] = c * np.outer(n, n)
                matrix[i] += coupling[i, j]
                if i < j:
                    activation_pairs += 1
                    floor_pairs += e < p.lubrication_floor
    residual = np.einsum("ikl,il->ik", matrix, local) - np.einsum("ijkl,jl->ik", coupling, local) - free
    iterations = [free, *[np.array(v) for v in t["lubrication"]]]
    pass_delta = [float(np.linalg.norm(b - a, axis=1).max()) for a, b in zip(iterations[:-1], iterations[1:])]
    iteration_checks = []
    for before, after in zip(iterations[:-1], iterations[1:]):
        rhs = free + np.einsum("ijkl,jl->ik", coupling, before)
        expected = 0.5 * (before + np.linalg.solve(matrix, rhs[..., None])[..., 0])
        error = float(np.max(np.abs(expected - after)))
        # NumPy LAPACK solve and Warp's explicit 3x3 inverse have different
        # rounding near the lubrication floor. This is a read-only algebra
        # cross-check (1 pm/s), not a physical acceptance-budget change.
        assert error < 1e-12, "Local normal-resistance diagnostic must reproduce the actual kernel"
        iteration_checks.append(error)
    preconditioned = np.linalg.solve(matrix, residual[..., None])[..., 0]
    worst = sorted(pair_rows, key=lambda r: r["force_norm_N"], reverse=True)[:8]
    report = dict(
        path=str(path.relative_to(BASE)),
        source_sha256=sha(path / "CONTACT_TRACE.json.gz"),
        mode=data["numerical_mode"],
        selected_step=dict(
            step=t["step_number_one_based"],
            start_s=t["start_s"],
            end_s=t["attempted_end_s"],
            accepted=t["accepted"],
        ),
        first_rejection=data["first_rejection_step"],
        audited_steps=len(data["traces"]),
        pre_pbd_penetration_nm=t["prediction"]["max_overlap_m"] * 1e9,
        post_pbd_residual_nm=t["final"]["max_overlap_m"] * 1e9,
        max_single_pass_nm=max(r["max_single_pass_nm"] for r in rows),
        residual_budget_nm=17.8,
        single_pass_budget_nm=178,
        rounds=rows,
        residual_violating_pairs=t["residual_violating_pairs"],
        residual_monotone=bool(np.all(np.diff([r["max_residual_nm"] for r in rows]) <= 1e-8)),
        last_to_previous_residual_ratio=rows[-1]["max_residual_nm"] / rows[-2]["max_residual_nm"],
        last_to_first_residual_ratio=rows[-1]["max_residual_nm"] / rows[0]["max_residual_nm"],
        redistribution_examples=chain,
        top4_correction_share=float(np.sort(cumulative)[-4:].sum() / cumulative.sum()),
        correction_path_per_id_nm=cumulative * 1e9,
        wall_correction_max_nm=max(np.linalg.norm(w["correction_m"], axis=1).max() for w in t["wall_passes"])
        * 1e9,
        force_audit=dict(
            primary_max_N=float(np.linalg.norm(t["primary_force_N"], axis=1).max()),
            secondary_max_N=float(np.linalg.norm(secondary, axis=1).max()),
            pair_sum_max_difference_N=float(np.abs(pair_force.sum(axis=1) - secondary).max()),
            numpy_reference_max_difference_N=float(np.abs(reference - secondary).max()),
            strongest_pairs=worst,
            min_closure_determinant=min(r["determinant_abs"] for r in pair_rows),
            max_coupled_pressure_Pa=max(max(r["coupled_pressure_abs_Pa"]) for r in pair_rows),
            physical_calibration=False,
        ),
        lubrication=dict(
            iterations=len(t["lubrication"]),
            active_pairs=activation_pairs,
            floor_limited_pairs=int(floor_pairs),
            matrix_condition_max=float(np.linalg.cond(matrix).max()),
            pass_max_velocity_change_m_s=pass_delta,
            coupled_equation_residual_max_m_s=float(np.linalg.norm(residual, axis=1).max()),
            residual_over_max_free_speed=float(
                np.linalg.norm(residual, axis=1).max() / np.linalg.norm(free, axis=1).max()
            ),
            unpreconditioned_residual_warning="Large resistance coefficients amplify this residual; use the preconditioned residual and reconstructed-iteration checks alongside it",
            preconditioned_residual_max_m_s=float(np.linalg.norm(preconditioned, axis=1).max()),
            preconditioned_residual_over_max_free_speed=float(
                np.linalg.norm(preconditioned, axis=1).max() / np.linalg.norm(free, axis=1).max()
            ),
            reconstructed_iteration_max_errors_m_s=iteration_checks,
            local_velocity_matches_last_iteration=bool(np.array_equal(local, iterations[-1])),
            scope="Algebraic diagnostic of two fixed Jacobi passes, not calibrated bubble hydrodynamics",
        ),
        model_scope="Pairwise fundamental closure, unchanged harmonic regularizer, two local lubrication iterations, diagonal-metric PBD; not multipole/collective radial/contact calibration",
    )
    return report


def morphology_compare(reference, candidate):
    with h5py.File(reference / "trajectory.h5") as a, h5py.File(candidate / "trajectory.h5") as b:
        n = min(len(a["time_s"]), len(b["time_s"]))
        ta, tb = a["time_s"][:n], b["time_s"][:n]
        assert np.array_equal(ta, tb)
        x, y = a["positions_m"][:n, 0], b["positions_m"][:n, 0]
        delta = np.linalg.norm(x - y, axis=-1)
        cdelta = np.linalg.norm(x.mean(1) - y.mean(1), axis=-1)
        aligned = [rigid_fit(xx, yy)[1] for xx, yy in zip(x, y)]
        events = []
        for h in [a, b]:
            et = h["edge_events/time_s"][:]
            ids = h["edge_events/batch_i_j_linked"][:]
            mask = et <= ta[-1] + 1e-14
            events.append(
                dict(
                    count=int(mask.sum()),
                    times_sha256=hashlib.sha256(et[mask].tobytes()).hexdigest(),
                    ids_sha256=hashlib.sha256(ids[mask].tobytes()).hexdigest(),
                )
            )
        return dict(
            reference=str(reference.relative_to(BASE)),
            candidate=str(candidate.relative_to(BASE)),
            common_end_s=float(ta[-1]),
            reference_after_failure="UNAVAILABLE",
            max_fixed_id_position_difference_um=float(delta.max() * 1e6),
            max_centroid_difference_nm=float(cdelta.max() * 1e9),
            max_aligned_rmsd_um=float(max(aligned) * 1e6),
            max_rg_difference_um=float(
                np.abs(a["metrics/global_rg_m"][:n, 0] - b["metrics/global_rg_m"][:n, 0]).max() * 1e6
            ),
            max_aspect_difference=float(
                np.abs(
                    a["morphology/fixed_aspect_ratio"][:n, 0] - b["morphology/fixed_aspect_ratio"][:n, 0]
                ).max()
            ),
            component_labels_exact=bool(
                np.array_equal(a["metrics/labels"][:n, 0], b["metrics/labels"][:n, 0])
            ),
            largest_sizes=[np.unique(h["metrics/largest_size"][:n, 0]).tolist() for h in [a, b]],
            accepted_edge_events=events,
            numerical_equivalence_admitted=False,
            scope="Sensitivity only; no preregistered physics-equivalence tolerance. Unchanged component count does not imply unchanged internal arrangement.",
        )


def yield_report():
    old = js(CAMPAIGN / "RESULTS.json")
    isolation = js(BASE / "isolation_v1/RESULTS.json")
    rows = [r for tree in isolation["trees"] for r in tree["attempts"]]
    # Production includes original failures/preparation/static and every isolated retry.
    iso_production = isolation["static_prepare_s"] + sum(r["timing"]["complete_run_wall_s"] for r in rows)
    iso_validation = isolation["process_body_s"] - iso_production
    stages = [
        ("F10D original", 2, old["total_generation_s"]),
        ("Original + isolation", 2, old["total_generation_s"] + iso_production),
    ]
    extension = js(BASE / "extension_v1/RESULTS.json")
    extra_valid = sum(r["status"] == "COMPLETE" and r["target_end_s"] == 0.01 for r in extension["rows"])
    ext_production = extension["static_prepare_s"] + sum(
        r["timing"]["complete_run_wall_s"] for r in extension["rows"]
    )
    diagnostics = sum(
        js(BASE / label / "RESULTS.json")["process_body_s"] for label in ["contact_v1", "alternatives_v1"]
    )
    total_all = (
        old["total_generation_s"] + isolation["process_body_s"] + diagnostics + extension["process_body_s"]
    )
    stages.append(("All attempts + diagnostics", 2 + extra_valid, total_all))
    envelope = js(BASE / "isolation_v1_SUPERVISOR.json")
    d = envelope["stop"] - envelope["start"]
    cold_old = next(
        r["full_process_wall_s"] for r in js(HISTORY / "campaign_v2/JOBS.json") if r["module"] == "diversity"
    )
    cold_others = sum(
        js(BASE / (label + "_PROCESS.json"))["full_process_wall_s"]
        for label in ["contact_v1", "alternatives_v1", "extension_v1"]
    )
    timings = [js(p) for p in BASE.glob("*_v1/*/PIPELINE_TIMING.json")]
    historical = [js(p) for p in CAMPAIGN.glob("*/PIPELINE_TIMING.json")]
    breakdown = {
        key: sum(t.get(key, 0) for t in timings + historical)
        for key in [
            "model_prepare_s",
            "graph_capture_s",
            "integration_quality_wall_s",
            "sync_readback_wall_s",
            "cuda_event_trajectory_span_s",
            "morphology_quality_output_prepare_s",
            "hdf5_output_s",
            "observer_output_s",
        ]
    }
    return dict(
        independent_seed_families=4,
        prepared_valid=2,
        baseline_complete_control_trajectories=2,
        recovered_by_isolation=0,
        independent_failed_retries=4,
        alternative_only_complete_trajectories=extra_valid,
        complete_same_mode_control_triplets=0,
        complete_control_triplets_per_gpu_hour=0,
        original_generation_s=old["total_generation_s"],
        isolation_production_replay_s=iso_production,
        isolation_validation_s=iso_validation,
        diagnostic_body_s=diagnostics,
        extension_production_s=ext_production,
        all_attempts_body_s=total_all,
        body_rates=[
            dict(
                stage=name,
                valid_completed_trajectories=n,
                complete_control_triplets=0,
                workflow_wall_s=wall,
                valid_completed_trajectories_per_gpu_hour=3600 * n / wall,
            )
            for name, n, wall in stages
        ],
        timing_breakdown=breakdown,
        all_process_wall_interval_s=[cold_old + max(0, d - 1) + cold_others, cold_old + d + 1 + cold_others],
        all_process_rate_interval_per_gpu_hour=[
            3600 * (2 + extra_valid) / (cold_old + d + 1 + cold_others),
            3600 * (2 + extra_valid) / (cold_old + max(0, d - 1) + cold_others),
        ],
        isolation_cold_timing=dict(
            source="Supervisor integer-second start/stop; exact cold timer absent for initial isolation job, do not invent precision",
            measured_difference_s=d,
            true_duration_interval_s=[max(0, d - 1), d + 1],
        ),
        production_only_with_alternative_s=old["total_generation_s"] + iso_production + ext_production,
        production_only_rate=3600
        * (2 + extra_valid)
        / (old["total_generation_s"] + iso_production + ext_production),
        cost_scope="Cross-stage historical workflow accounting, not a paired speed benchmark. CUDA event spans include host gaps and are not pure kernel time; GPU-hour uses exclusive device reservation including static/CPU/output work. Diagnostic and failed replay costs retained.",
        training_ready=False,
    )


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    contact = [
        contact_report(p.parent)
        for folder in ["contact_v1", "alternatives_v1"]
        for p in sorted((BASE / folder).glob("*/CONTACT_TRACE.json.gz"))
    ]
    write(out / "CONTACT_ANALYSIS.json", contact)
    reference = BASE / "contact_v1/O_1000_1102_C0"
    morphology = [
        morphology_compare(reference, p.parent)
        for folder in ["alternatives_v1", "extension_v1"]
        for p in sorted((BASE / folder).glob("*/trajectory.h5"))
    ]
    write(out / "MORPHOLOGY_SENSITIVITY.json", morphology)
    yield_data = yield_report()
    write(out / "EFFECTIVE_YIELD.json", yield_data)
    catalog = []
    paths = [*sorted(CAMPAIGN.glob("*/trajectory.h5")), *sorted(BASE.glob("*_v1/*/trajectory.h5"))]
    for p in paths:
        s = js(p.parent / "SUMMARY.json")
        cfg = js(p.parent / "CONFIG.json")
        with h5py.File(p) as h:
            accepted = float(h["time_s"][-1])
            interval = dict(h["rejected_interval"].attrs) if "rejected_interval" in h else None
        for i, sid in enumerate(s["sample_id"]):
            x = np.array(cfg["initial_positions_m"][i], dtype=np.float64)
            catalog.append(
                dict(
                    sample_id=sid,
                    initial_state_family=s["initial_state_family"][i],
                    split_group=s["split_groups"][i],
                    initial_positions_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
                    path=str(p.parent.relative_to(ROOT)),
                    hdf5_sha256=sha(p),
                    source_commit=s["source"]["execution_commit"],
                    mode=s.get("original_physics_mode"),
                    status=s["status"],
                    quality=s["per_sample_quality"][i],
                    full_scheduled_duration=accepted == cfg["duration_s"],
                    scheduled_duration_s=cfg["duration_s"],
                    full_10ms=accepted == 0.01 and s["status"] == "COMPLETE",
                    accepted_saved_end_s=accepted,
                    rejected_interval=interval,
                    attempt_cost=js(p.parent / "PIPELINE_TIMING.json"),
                    training_ready=False,
                )
            )
    write(
        out / "SAMPLE_CATALOG.json",
        dict(independent_seeds=4, records=catalog, correlated_controls_variants_retries_share_split=True),
    )
    write(out / "ANALYSIS_IDENTITY.json", identity())
    print(json.dumps(dict(yield_data=yield_data, morphology=morphology), indent=2))


if __name__ == "__main__":
    main()
