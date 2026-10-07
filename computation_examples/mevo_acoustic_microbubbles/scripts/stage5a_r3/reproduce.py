"""Full-system R2 boundary replays, original states and exact geometry evidence."""

import copy
from concurrent.futures import ProcessPoolExecutor
from time import perf_counter
import h5py
import numpy as np
from mevo_mb.bubbles.motion_r2 import SemiImplicit
from mevo_mb.bubbles.motion_sweep_r2 import segment_triangle2
from mevo_mb.bubbles.geometry_exact_r3 import segment_triangle, record
from mevo_mb.bubbles.constrained_motion import SubdivideInterval, AcceptedCache
from scripts.stage5a_r2.run import restore, checkpoint
from mevo_mb.bubbles.motion_recovery_r2 import begin_recovery, accepted_recovery
from scripts.stage5a_r3.common import BASE, R2, B2_RUNS, load, read, write, sha, identity


def reconstruct_controller(a, rejected, index, cfg):
    controller = dict(
        next_dt=cfg["initial_dt_s"],
        last_error=None,
        last_reason=None,
        consecutive_rejections=0,
        elapsed_s=0.0,
    )
    events = [(r["wall_s"], "r", r) for r in rejected if r["elapsed_s"] < a[index - 1]["elapsed_s"] - 1e-15]
    events += [(r["wall_s"], "a", r) for i, r in enumerate(a[:index]) if i % 2 == 1]
    for _, kind, row in sorted(events, key=lambda z: z[0]):
        if kind == "r":
            dt = row["dt_s"]
            controller.update(
                next_dt=dt / 2,
                last_reason=row["reason"],
                consecutive_rejections=controller["consecutive_rejections"] + 1,
            )
            if row["trial_method"] == "S1":
                exc = SubdivideInterval(row["reason"], row["diagnostics"])
                ceiling = min(
                    cfg["max_dt_s"],
                    cfg["duration_s"] - row["elapsed_s"],
                    (np.floor((row["elapsed_s"] + 1e-14) / 1e-5) + 1) * 1e-5 - row["elapsed_s"],
                )
                begin_recovery(controller, exc, dt, ceiling, cfg["tolerance_factor"])
        else:
            dt = 2 * row["dt_s"]
            factor = min(2.0, max(0.5, 0.85 / max(row["error_ratio"], 1e-8) ** 0.5))
            controller.update(
                next_dt=accepted_recovery(
                    controller, row["substep_method"], dt, dt * factor, cfg["max_dt_s"]
                ),
                last_error=row["error_ratio"],
                last_reason="ACCEPTED",
                consecutive_rejections=0,
                elapsed_s=row["elapsed_s"],
            )
            controller["total_accepted_steps"] = controller.get("total_accepted_steps", 0) + 2
    return controller


def geometry_audit(sweep, x, y, r, ids):
    g0, n0, _ = sweep.geometry.wall_gap_and_normal(x, r)
    g1, n1, _ = sweep.geometry.wall_gap_and_normal(y, r)
    length = np.linalg.norm(y - x, axis=1)
    rows = []
    for i in np.flatnonzero(g0 - length < 1e-10):
        low = np.nextafter(np.minimum(x[i], y[i]) - r[i] - 1e-11, -np.inf)
        high = np.nextafter(np.maximum(x[i], y[i]) + r[i] + 1e-11, np.inf)
        sweep.locator.FindCellsWithinBounds(np.column_stack([low, high]).ravel(), sweep.ids)
        cells = np.array([sweep.ids.GetId(k) for k in range(sweep.ids.GetNumberOfIds())], int)
        cells = cells[
            np.all(sweep.triangles[cells].max(1) >= low, axis=1)
            & np.all(sweep.triangles[cells].min(1) <= high, axis=1)
        ]
        values = []
        for cell in cells:
            tri = sweep.triangles[cell]
            exact = segment_triangle(x[i], y[i], *tri)
            row = record(exact, r[i], 1e-11)
            row.update(
                wall_triangle_id=int(cell),
                wall_vertex_ids=sweep.mesh.faces.reshape(-1, 4)[cell, 1:].tolist(),
                triangle_vertices_m=tri.tolist(),
                float_distance_gap_m=float(np.sqrt(segment_triangle2(x[i], y[i], *tri)) - r[i]),
            )
            values.append(row)
        worst = min(values, key=lambda z: z["gap_m"]) if values else None
        rows.append(
            dict(
                particle_id=int(ids[i]),
                start_m=x[i].tolist(),
                end_m=y[i].tolist(),
                radius_m=float(r[i]),
                signed_start_gap_m=float(g0[i]),
                signed_endpoint_gap_m=float(g1[i]),
                start_normal=n0[i].tolist(),
                end_normal=n1[i].tolist(),
                normal_jump=float(np.linalg.norm(n1[i] - n0[i])),
                coordinate_changed=bool(np.any(x[i] != y[i])),
                actual_displacement_m=float(length[i]),
                candidate_triangles=len(cells),
                minimum=worst,
                all_triangles=values,
            )
        )
    return rows


def case_job(kind):
    started = perf_counter()
    run = R2 / "runs" / B2_RUNS[kind]
    out = BASE / "diagnosis/reproduction_v1" / kind
    out.mkdir(parents=True, exist_ok=False)
    tracker, sweep, cfg, initial, inputs = load("B2", "numpy")
    data, meta, saved_cache = restore(run / "CHECKPOINT.npz")
    runid = read(run / "IDENTITY.json")
    assert sha(run / "IDENTITY.json") == meta["identity"]
    accepted, rejected = read(run / "ACCEPTED.jsonl"), read(run / "REJECTED.jsonl")
    with h5py.File(run / "trajectories.h5") as h:
        times = h["time_s"][:]
        masks = h["active"][:]
        assert np.array_equal(h["positions_m"][-1], data["positions_m"], equal_nan=True)
        near = next(
            z
            for z in rejected
            if z["reason"] == "SWEEP_NOT_CERTIFIED" and z["elapsed_s"] >= accepted[-1]["elapsed_s"] - 2e-6
        )
        event_index = int(np.argmin(abs(times - near["time_s"])))
        assert abs(times[event_index] - near["time_s"]) < 1e-15 and event_index % 2 == 0
        selected = [
            ("before_persistent_sweep", event_index),
            ("last_successful_macro_start", len(times) - 3),
            ("terminal", len(times) - 1),
        ]
        fixtures = []
        for label, index in selected:
            full = h["positions_m"][index]
            mask = masks[index]
            ids = h["particle_id"][:]
            rr = h["radii_m"][:]
            controller = (
                copy.deepcopy(meta["controller"])
                if label == "terminal"
                else reconstruct_controller(accepted, rejected, index, runid)
            )
            velocity = np.zeros_like(full)
            cache = saved_cache if label == "terminal" else AcceptedCache()
            scope = (
                "saved exact checkpoint"
                if label == "terminal"
                else "actual coordinates; controller reconstructed from chronological ledgers; warm cache independently reconstructed from preceding actual substep"
            )
            mismatch = None
            if label == "terminal":
                velocity = data["velocity_m_s"]
            elif index:
                prev = h["positions_m"][index - 1][mask]
                actual = accepted[index - 1]
                tracker.translation.cache = AcceptedCache()
                step = SemiImplicit(
                    tracker.translation,
                    tracker.cylinder,
                    tracker.evaluate,
                    sweep,
                    method=actual["substep_method"],
                    nonlinear_strategy="quasi_newton",
                )
                replay, u, cache, diag = step.advance(
                    prev,
                    rr[mask],
                    ids[mask],
                    float(times[index - 1]),
                    actual["dt_s"],
                    "fixture_reconstruction",
                )
                mismatch = float(np.max(np.linalg.norm(replay - full[mask], axis=1)))
                velocity[mask] = u
            np.savez_compressed(
                out / f"{label}_source.npz",
                positions_m=full,
                radii_m=rr,
                particle_id=ids,
                active=mask,
                time_s=times[index],
            )
            fixture_identity = dict(
                source=str(run),
                source_trajectory_sha256=sha(run / "trajectories.h5"),
                index=index,
                cache_scope=scope,
                reconstructed_endpoint_difference_m=mismatch,
                original_state_unchanged=True,
                run_identity_sha256=sha(run / "IDENTITY.json"),
            )
            write(out / f"{label}_IDENTITY.json", fixture_identity)
            checkpoint(
                out / f"{label}.npz",
                full,
                velocity,
                rr,
                ids,
                mask,
                float(times[index]),
                controller,
                cache,
                sha(out / f"{label}_IDENTITY.json"),
            )
            fixtures.append(dict(label=label, path=str(out / f"{label}.npz"), **fixture_identity))
    # Replay the final failed full S0 candidate using the actual saved cache.
    mask = data["active"]
    x = data["positions_m"][mask]
    r = data["radii_m"][mask]
    ids = data["particle_id"][mask]
    t = float(data["time_s"])
    trials = []
    for dt in [rejected[-1]["dt_s"], 1e-9, 5e-8]:
        tracker.translation.cache = copy.deepcopy(saved_cache)
        captured = []
        original = tracker.translation.constrained_stage
        last_qp = {}

        def observe_qp(result, constraints, pids, stage, **kw):
            u, c, d = original(result, constraints, pids, stage, **kw)
            last_qp.clear()
            last_qp.update(
                diag=d,
                A=result.matrix,
                G=constraints.matrix,
                lower=constraints.lower,
                U=u,
                cache=c,
                b=result.matrix @ result.velocity.ravel(),
                constraint_ids=constraints.ids,
            )
            return u, c, d

        tracker.translation.constrained_stage = observe_qp

        def observe_sweep(start, end, radii, tol):
            old = sweep(start, end, radii, tol)
            captured.append(
                dict(
                    start=start.copy(),
                    end=end.copy(),
                    old=old,
                    exact=geometry_audit(sweep, start, end, radii, ids),
                )
            )
            return old

        step = SemiImplicit(
            tracker.translation,
            tracker.cylinder,
            tracker.evaluate,
            observe_sweep,
            method="S0",
            nonlinear_strategy="quasi_newton",
        )
        begin = perf_counter()
        outcome = None
        try:
            step.advance(x, r, ids, t, dt, "full")
            outcome = "ACCEPTED"
        except SubdivideInterval as exc:
            outcome = exc.reason
        finally:
            tracker.translation.constrained_stage = original
        assert len(captured) == 1
        c = captured[0]
        name = f"candidate_{dt:.12g}"
        arrays = dict(
            start_m=c["start"],
            end_m=c["end"],
            radii_m=r,
            particle_id=ids,
            requested_dt_s=dt,
            time_s=t,
            A_data=last_qp["A"].data,
            A_indices=last_qp["A"].indices,
            A_indptr=last_qp["A"].indptr,
            A_shape=last_qp["A"].shape,
            G_data=last_qp["G"].data,
            G_indices=last_qp["G"].indices,
            G_indptr=last_qp["G"].indptr,
            G_shape=last_qp["G"].shape,
            lower=last_qp["lower"],
            b=last_qp["b"],
            U=last_qp["U"],
        )
        np.savez_compressed(out / f"{name}.npz", **arrays)
        trials.append(
            dict(
                dt_s=dt,
                actual_time_change_s=float((t + dt) - t),
                outcome=outcome,
                old_safety=c["old"],
                exact_geometry=c["exact"],
                qp=last_qp["diag"],
                constraint_ids=last_qp["constraint_ids"],
                wall_s=perf_counter() - begin,
                path=str(out / f"{name}.npz"),
            )
        )
    result = dict(
        kind=kind,
        source_run=str(run),
        fixtures=fixtures,
        trials=trials,
        input_sha256=inputs,
        identity=identity(),
        wall_s=perf_counter() - started,
    )
    write(out / "RESULT.json", result)
    return result


def main():
    with ProcessPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(case_job, B2_RUNS))
    write(
        BASE / "final/B2_FAILURE_REPRODUCTION.json",
        dict(
            results=results,
            scope="Full-system local replays; no particle deletion or change to physical inputs",
        ),
    )
    print(
        [
            (
                r["kind"],
                [
                    (
                        t["dt_s"],
                        t["outcome"],
                        [p["minimum"]["classification"] for p in t["exact_geometry"] if p["minimum"]],
                    )
                    for t in r["trials"]
                ],
            )
            for r in results
        ],
        flush=True,
    )


if __name__ == "__main__":
    main()
