"""Common-state equations, identical-QP comparisons and legacy source attribution."""

import copy
import importlib.util
import sys
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import perf_counter
import h5py
import numpy as np
from scipy import sparse
from scipy.optimize import nnls
from mevo_mb.bubbles.constrained_motion import Constraints, sparse_rows
from mevo_mb.bubbles.motion_r2 import MotionR2, SemiImplicit
from mevo_mb.bubbles.motion_r3 import MotionR3, SemiImplicitR3
from mevo_mb.bubbles.motion_sweep_r3 import WallSweepR3
from mevo_mb.bubbles.translation_backends import _pairs, _gap_vectors
from mevo_mb.bubbles.motion_qp_r2 import certificate
from scripts.stage5a_r3.common import BASE, ROOT, load, read, write, sha, identity


def legacy_module(kind, archive):
    rel = "src/mevo_mb/bubbles/translation_backends.py"
    if kind == "archived":
        target = ROOT / ".runtime/archived_stage4l/translation_backends.py"
        target.parent.mkdir(parents=True, exist_ok=True)
        data = zipfile.ZipFile(archive / "source_snapshot.zip").read(rel)
        target.write_bytes(data)
        assert sha(target) == read(archive / "IDENTITY.json")["source_sha256"][rel]
    else:
        target = Path("/workspace/mevo_acoustic_microbubbles_stage4l") / rel
    name = "mevo_mb.bubbles._r3_audit_" + kind
    spec = importlib.util.spec_from_file_location(name, target)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module, target


def legacy_constraints(result, x, r, velocity, tolerance, horizon, ids, backend):
    radius = backend.mono_radius(r)
    if horizon is None:
        cutoff = 2 * radius + tolerance
    else:
        centered = velocity - np.mean(velocity, axis=0)
        reach = 2 * np.max(np.linalg.norm(centered, axis=1)) * horizon
        cutoff = 2 * radius + backend.floor * radius + reach
    pairs = _pairs(x, cutoff)
    gaps, normal = _gap_vectors(x, pairs, radius)
    wg, wn, _ = result.wall
    wi = np.flatnonzero(
        wg
        <= (
            tolerance
            if horizon is None
            else np.maximum(backend.floor * radius, np.linalg.norm(velocity, axis=1) * horizon)
        )
    )
    matrix = sparse_rows(pairs, normal, wi, wn[wi], len(r))
    h = np.r_[gaps, wg[wi]]
    lower = np.zeros(len(h)) if horizon is None else -np.maximum(h, 0) / horizon
    lower[h <= tolerance] = 0
    keys = tuple(("pair", int(ids[i]), int(ids[j])) for i, j in pairs)
    keys += tuple(("wall", int(ids[i])) for i in wi)
    return Constraints(matrix, lower, keys, pairs)


def independent_velocity_certificate(result, c, velocity):
    vs = max(float(np.max(abs(result.velocity))), 1e-30)
    a = result.matrix / result.scale
    b = a @ (result.velocity.ravel() / vs)
    u = velocity.ravel() / vs
    residual = a @ u - b
    slack = c.matrix @ u - c.lower / vs
    active = np.flatnonzero(slack <= 1e-7)
    lam = np.zeros(len(c.lower))
    if len(active):
        lam[active], _ = nnls(c.matrix[active].T.toarray(), residual, maxiter=10000)
    cert = certificate(a, b, c.matrix, c.lower / vs, u, lam, velocity_scale=vs, force_scale=result.scale * vs)
    return cert, (c.matrix.T @ lam) * result.scale * vs


def job(item):
    case, source = item
    started = perf_counter()
    out = BASE / "diagnosis/common_state_v1" / f"{case}_{source}"
    out.mkdir(parents=True, exist_ok=True)
    lineage = read(BASE / "diagnosis" / f"{case}_LINEAGE_V1.json")
    archived = Path(lineage["paths"]["archived"]).parent
    trajectory = Path(lineage["paths"]["archived" if source == "legacy" else "nominal"])
    with h5py.File(trajectory) as h:
        k = int(np.searchsorted(h["time_s"][:], 3.00003, side="right") - 1)
        mask, full, t = h["active"][k], h["positions_m"][k], float(h["time_s"][k])
        ids, r = h["particle_id"][:][mask], h["radii_m"][:][mask]
    x = full[mask]
    tracker, old_sweep, cfg, state, inputs = load(case)
    forcing = tracker.evaluate(x, r, t)
    wall = tracker.cylinder.wall_gap_and_normal(x, r)
    cold = copy.deepcopy(wall)
    tracker.cylinder._distance_cache = None
    cold2 = tracker.cylinder.wall_gap_and_normal(x, r)
    np.savez_compressed(
        out / "INPUT.npz",
        positions_m=x,
        radii_m=r,
        particle_id=ids,
        time_s=t,
        **forcing,
        wall_gap_m=wall[0],
        wall_normal=wall[1],
    )
    force = forcing["primary_force_N"] + forcing["secondary_force_N"]
    flow = forcing["flow_velocity_m_s"]
    backends, source_hashes = {}, {}
    for kind in ["archived", "fresh"]:
        module, path = legacy_module(kind, archived)
        backends[kind] = module.LubricationCoupledTranslation(
            cfg["medium"]["viscosity_Pa_s"], **cfg["particles"]["lubrication"]
        )
        source_hashes[kind] = dict(path=str(path), sha256=sha(path))
    for kind, cls in [("R2", MotionR2), ("R3", MotionR3)]:
        backends[kind] = cls(
            cfg["medium"]["viscosity_Pa_s"],
            **cfg["particles"]["lubrication"],
            **cfg["particles"]["constrained_solver"],
        )
    results = {key: backend.solve(x, r, force, flow, wall) for key, backend in backends.items()}
    base = results["archived"]
    matrices = []
    for key, a in results.items():
        diff = (a.matrix - base.matrix).toarray()
        matrices.append(
            dict(
                method=key,
                A_max_abs_N_s_m=float(np.max(abs(diff))),
                free_velocity_max_abs_m_s=float(np.max(abs(a.velocity - base.velocity))),
                b_max_abs_N=float(
                    np.max(abs(a.matrix @ a.velocity.ravel() - base.matrix @ base.velocity.ravel()))
                ),
                independent_linear_residual=a.diagnostics["relative_residual"],
            )
        )
    qp_rows = []
    for horizon in [None, 1e-7, 1e-6]:
        for legacy in ["archived", "fresh"]:
            a = results[legacy]
            backend = backends[legacy]
            c = legacy_constraints(a, x, r, a.velocity, 1e-11, horizon, ids, backend)
            hname = "instant" if horizon is None else f"{horizon:g}"
            label = legacy + "_" + hname
            sparse.save_npz(out / (label + "_A.npz"), a.matrix)
            sparse.save_npz(out / (label + "_G.npz"), c.matrix)
            try:
                u, d = backend.contact(a, x, r, a.velocity, 1e-11, horizon=horizon)
                cert, reaction = independent_velocity_certificate(a, c, u)
                old = dict(status="RETURNED", diagnostic=d, independent_certificate=cert)
            except Exception as exc:
                u, reaction = None, None
                old = dict(status="FAILED", error=repr(exc))
            sols = []
            save = dict(lower=c.lower, free_velocity_m_s=a.velocity, particle_id=ids)
            if u is not None:
                save.update(legacy_velocity_m_s=u, legacy_reaction_N=reaction)
            for method in ["R2", "R3"]:
                v, cache, diag = backends[method].constrained_stage(a, c, ids, "common_" + label)
                lam = np.array([cache.dual[k] for k in c.ids])
                react = c.matrix.T @ lam
                cert_new, _ = independent_velocity_certificate(a, c, v)
                sols.append(
                    dict(
                        method=method,
                        kkt=diag["kkt"],
                        independent_certificate=cert_new,
                        max_velocity_difference_m_s=None
                        if u is None
                        else float(np.max(np.linalg.norm(v - u, axis=1))),
                        reaction_difference_l2_N=None
                        if reaction is None
                        else float(np.linalg.norm(react - reaction)),
                    )
                )
                save[method + "_velocity_m_s"] = v
                save[method + "_reaction_N"] = react
            np.savez_compressed(out / (label + "_data.npz"), **save)
            qp_rows.append(dict(legacy=legacy, horizon_s=horizon, rows=len(c.lower), old=old, shared_QP=sols))
    steps = []
    for name, cls in [("R2", SemiImplicit), ("R3", SemiImplicitR3)]:
        sweep = old_sweep if name == "R2" else WallSweepR3(tracker.cylinder, old_sweep.path)
        for method in ["S0", "S1"]:
            step = cls(
                backends[name],
                tracker.cylinder,
                tracker.evaluate,
                sweep,
                method=method,
                nonlinear_strategy="quasi_newton",
            )
            try:
                y, u, cache, diag = step.advance(x, r, ids, t, 1e-7, "common_one_step")
                np.savez_compressed(out / f"{name}_{method}_step.npz", positions_m=y, velocity_m_s=u)
                steps.append(dict(method=name + "_" + method, status="PASS", diagnostic=diag))
            except Exception as exc:
                steps.append(
                    dict(
                        method=name + "_" + method,
                        status="FAIL",
                        error=repr(exc),
                        diagnostic=getattr(exc, "diagnostics", None),
                    )
                )
    # Use the unchanged tracker legacy branch with the exact archived/fresh
    # translation module. This is a controlled adapter, not a new frozen reference.
    import mevo_mb.bubbles.tracker as tracking
    from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep

    tracker.motion_constrained = False
    tracker.lubrication_enabled = True
    ordered = OrderedContactSweep(out / "ordered_cache")
    original = tracking.project_contacts
    tracking.project_contacts = ordered.project
    try:
        for name in ["archived", "fresh"]:
            tracker.translation = backends[name]
            try:
                f = tracker.evaluate(x, r, t)
                discrepancy = {k: float(np.max(abs(f[k] - forcing[k]))) for k in forcing}
                y, wallhit, contacts, correction, iterations = tracker.step(x, r, t, 1e-7, f)
                np.savez_compressed(out / f"{name}_step.npz", positions_m=y)
                steps.append(
                    dict(
                        method=name,
                        status="PASS",
                        force_query_differences=discrepancy,
                        final_projection_m=correction,
                        contact_iterations=iterations,
                    )
                )
            except Exception as exc:
                steps.append(dict(method=name, status="FAIL", error=repr(exc)))
    finally:
        tracking.project_contacts = original
    row = dict(
        case=case,
        state_source=source,
        actual_time_s=t,
        frame=k,
        active_count=len(r),
        source=str(trajectory),
        source_sha256=sha(trajectory),
        source_modules=source_hashes,
        force_scope="Identical full-system acoustic/response/flow inputs for all matrix solves; legacy adapter independently re-queries them",
        cold_cache_max_gap_difference_m=float(np.max(abs(cold[0] - cold2[0]))),
        cold_cache_max_normal_difference=float(np.max(abs(cold[1] - cold2[1]))),
        assembly=matrices,
        identical_qp_comparisons=qp_rows,
        one_step=steps,
        wall_s=perf_counter() - started,
    )
    write(out / "RESULT.json", row)
    return row


def main():
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(job, [(c, s) for c in ["B2", "B3"] for s in ["legacy", "R2"]]))
    write(
        BASE / "diagnosis/COMMON_STATE_REPLAY_V1.json",
        dict(rows=rows, identity=identity(), charged_experiment_wall_s=sum(r["wall_s"] for r in rows)),
    )
    print([(r["case"], r["state_source"], r["wall_s"]) for r in rows], flush=True)


if __name__ == "__main__":
    main()
