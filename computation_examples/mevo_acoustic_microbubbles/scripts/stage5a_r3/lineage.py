"""All-active, fixed-ID reference lineage audit before any R3 tuning."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import perf_counter
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, R2, read, write, sha, identity


def sample_all(path, clock):
    with h5py.File(path) as h:
        times = h["time_s"][:]
        if np.any(np.asarray(clock) < times[0] - 1e-14) or np.any(np.asarray(clock) > times[-1] + 1e-14):
            raise ValueError("Requested comparison clock exceeds actual stored interval")
        ids = h["particle_id"][:]
        radii = h["radii_m"][:]
        masks = h["active"][:]
        assert np.all(masks == masks[0])
        mask = masks[0]
        indices = []
        for t in clock:
            j = int(np.searchsorted(times, t, side="right") - 1)
            nxt = min(j + 1, len(times) - 1)
            if abs(times[nxt] - t) <= 16 * np.finfo(float).eps * max(1, abs(t)):
                j = nxt
            indices.append(max(0, j))
        indices = np.array(indices)
        unique, inverse = np.unique(indices, return_inverse=True)
        positions = h["positions_m"][unique][:, mask][inverse]
        uncertainty = np.zeros(len(indices))
        # A neighboring stored-frame chord is a rigorous hold-last-state bound
        # only for the explicitly recorded straight accepted path. Legacy
        # projection trajectories do not supply that continuous-path guarantee.
        straight = bool(h.attrs.get("all_accepted_substeps_recorded", False)) and str(
            h.attrs.get("actual_path", "")
        ).startswith("piecewise linear accepted half steps")
        certified = np.ones(len(indices), bool)
        for k, j in enumerate(indices):
            if j + 1 < len(times) and clock[k] - times[j] > 1e-14:
                uncertainty[k] = np.linalg.norm(
                    h["positions_m"][j + 1][mask] - h["positions_m"][j][mask], axis=1
                ).max()
                certified[k] = straight
        return dict(
            positions=positions,
            ids=ids[mask],
            radii=radii[mask],
            time=times[indices],
            index=indices,
            uncertainty_m=uncertainty,
            uncertainty_is_bound=certified,
            initial_full=h["positions_m"][0],
            active=mask,
            all_ids=ids,
        )


def summarize(norms):
    return dict(
        rms=float(np.sqrt(np.mean(norms**2))),
        p50=float(np.quantile(norms, 0.5)),
        p95=float(np.quantile(norms, 0.95)),
        p99=float(np.quantile(norms, 0.99)),
        max=float(norms.max()),
    )


def metrics(left, right):
    raw = right - left
    lc = left - left.mean(1, keepdims=True)
    rc = right - right.mean(1, keepdims=True)
    centered = rc - lc
    rigid = []
    for a, b in zip(lc, rc):
        u, _, vt = np.linalg.svd(b.T @ a)
        correction = np.eye(3)
        correction[-1, -1] = np.linalg.det(u @ vt)
        rigid.append(b @ (u @ correction @ vt) - a)
    return {
        k: np.linalg.norm(v, axis=2)
        for k, v in [("raw", raw), ("translation_removed", centered), ("rigid_fit", np.array(rigid))]
    }


def compare_all(case, name, left, right, group):
    with h5py.File(left) as a, h5py.File(right) as b:
        assert np.array_equal(a["particle_id"][:], b["particle_id"][:])
        assert np.array_equal(a["positions_m"][0], b["positions_m"][0], equal_nan=True)
        assert np.array_equal(a["active"][0], b["active"][0])
        start = max(float(a["time_s"][0]), float(b["time_s"][0]))
        end = min(float(a["time_s"][-1]), float(b["time_s"][-1]))
    clock = start + np.arange(int(np.floor((end - start) / 1e-5 + 1e-7)) + 1) * 1e-5
    if end - clock[-1] > 1e-14:
        clock = np.r_[clock, end]
    a, b = sample_all(left, clock), sample_all(right, clock)
    assert np.array_equal(a["ids"], b["ids"])
    subsets = {
        "all_active": np.ones(len(a["ids"]), bool),
        "original_group": np.isin(a["ids"], group),
        "complement": ~np.isin(a["ids"], group),
    }
    rows = {}
    arrays = dict(
        clock_s=clock,
        particle_id=a["ids"],
        radii_m=a["radii"],
        left_positions_m=a["positions"],
        right_positions_m=b["positions"],
        left_actual_time_s=a["time"],
        right_actual_time_s=b["time"],
        left_frame_index=a["index"],
        right_frame_index=b["index"],
        left_sampling_bound_m=a["uncertainty_m"],
        right_sampling_bound_m=b["uncertainty_m"],
        left_sampling_bound_certified=a["uncertainty_is_bound"],
        right_sampling_bound_certified=b["uncertainty_is_bound"],
    )
    for label, mask in subsets.items():
        errors = metrics(a["positions"][:, mask], b["positions"][:, mask])
        pids = a["ids"][mask]
        rows[label] = dict(particle_ids=pids.tolist(), count=len(pids), errors_m={})
        for mode, v in errors.items():
            j, i = np.unravel_index(np.argmax(v), v.shape)
            rows[label]["errors_m"][mode] = dict(
                endpoint=summarize(v[-1]),
                whole_sampled_interval=summarize(v),
                worst_particle_id=int(pids[i]),
                worst_time_s=float(clock[j]),
            )
            arrays[f"{label}_{mode}_errors_m"] = v
        norm = errors["raw"]
        threshold = 0.05 * float(a["radii"][0])
        crossed = np.flatnonzero(norm.max(1) > threshold)
        rows[label]["first_over_resolution_target"] = (
            None
            if not len(crossed)
            else dict(
                time_s=float(clock[crossed[0]]),
                prior_sample_time_s=float(clock[max(0, crossed[0] - 1)]),
                target_m=threshold,
                worst_particle_id=int(pids[np.argmax(norm[crossed[0]])]),
                max_error_m=float(norm[crossed[0]].max()),
            )
        )
    output = BASE / "comparisons" / f"{case}_{name}_all_active.npz"
    np.savez_compressed(output, **arrays)
    return dict(
        case=case,
        name=name,
        left=str(left),
        right=str(right),
        left_sha256=sha(left),
        right_sha256=sha(right),
        duration_s=end - start,
        all_active=len(a["ids"]),
        groups=rows,
        samples=str(output),
        sample_sha256=sha(output),
        sample_age_max_s=[float(np.max(clock - a["time"])), float(np.max(clock - b["time"]))],
        sampling_bound_max_m=[float(a["uncertainty_m"].max()), float(b["uncertainty_m"].max())],
        sampling_audit_version=2,
        sampling_bound_certified=bool(a["uncertainty_is_bound"].all() and b["uncertainty_is_bound"].all()),
        sampling_semantics="Exact stored common times need no hold error. At unmatched times, the neighboring-frame chord is a bound only for explicitly stored straight accepted paths; otherwise it is a descriptive displacement scale, not a rigorous legacy-path bound.",
        status="NOT_EVALUABLE",
        reason="Descriptive pre-R3 lineage audit; no independent all-active trajectory error bound established",
        no_coordinate_interpolation=True,
        original_ID_mapping=True,
    )


def flatten(value, prefix=""):
    if isinstance(value, dict):
        result = {}
        for key, v in value.items():
            result.update(flatten(v, prefix + "." + key if prefix else key))
        return result
    return {prefix: value}


def case_job(case):
    begin = perf_counter()
    campaign = read(R2 / "final/CAMPAIGN_V4A/RESULT.json")
    entries = {
        r["kind"]: r
        for r in campaign["ledger"]
        if r["case"] == case and r["duration"] == 0.005 and r.get("repeat") == 0 and r.get("summary")
    }
    archived = next(r for r in read(R2 / "fixtures/HOLDOUT_SELECTION.json") if r["case"] == case)
    paths = {
        "archived": Path(archived["source"]),
        **{k: Path(v["run"]) / "trajectories.h5" for k, v in entries.items()},
    }
    assert sha(paths["archived"]) == archived["source_sha256"]
    original_group = np.load(R2 / "comparisons" / f"{case}_0.001_v4a_legacy_reference.npz")["particle_ids"]
    pairings = [
        ("archived_vs_fresh", "archived", "reference"),
        ("fresh_vs_r2", "reference", "nominal"),
        ("archived_vs_r2", "archived", "nominal"),
        ("r2_refinement", "refined", "nominal"),
    ]
    comparisons = [
        compare_all(case, name, paths[left], paths[right], original_group) for name, left, right in pairings
    ]
    identities = {k: read(p.parent / "IDENTITY.json") for k, p in paths.items()}
    archive, fresh = identities["archived"], identities["reference"]
    sa, sb = archive["source_sha256"], fresh["source_sha256"]
    source_difference = [
        dict(path=k, archived=sa.get(k), fresh=sb.get(k))
        for k in sorted(sa.keys() & sb.keys())
        if k.startswith("src/") and sa[k] != sb[k]
    ]
    ca, cb = flatten(archive["config"]), flatten(fresh["config"])
    configuration_difference = [
        dict(key=k, archived=ca.get(k), fresh=cb.get(k))
        for k in sorted(ca.keys() | cb.keys())
        if ca.get(k) != cb.get(k)
    ]
    ia, ib = archive["input_sha256"], fresh["input_sha256"]
    common_inputs = [
        dict(path=k, sha256=ia[k], matches=ia[k] == ib[k]) for k in sorted(ia.keys() & ib.keys())
    ]
    row = dict(
        case=case,
        paths={k: str(v) for k, v in paths.items()},
        identity_sha256={k: sha(p.parent / "IDENTITY.json") for k, p in paths.items()},
        common_source_differences=source_difference,
        configuration_differences=configuration_difference,
        common_inputs=common_inputs,
        archived_only_input_paths=sorted(ia.keys() - ib.keys()),
        fresh_only_input_paths=sorted(ib.keys() - ia.keys()),
        archived_driver_sha256=archive["source_sha256"].get("scripts/stage4l/dense_probe.py"),
        fresh_driver_sha256=fresh["driver_sha256"],
        archived_output_clock_s=archive["config"]["particles"]["dt_s"],
        fresh_output_clock_s=fresh["output_clock_s"],
        reference_equivalence="UNRESOLVED until driver/time-horizon semantics and observed full-ID differences are audited; code/input hashes alone do not establish identical trajectories",
        comparisons=comparisons,
        wall_s=perf_counter() - begin,
    )
    write(BASE / "diagnosis" / f"{case}_LINEAGE_V1.json", row)
    return row


def main():
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(case_job, ["B2", "B3"]))
    write(
        BASE / "final/REFERENCE_LINEAGE_AUDIT.json",
        dict(rows=rows, identity=identity(), stage="Pre-R3 development audit"),
    )
    write(
        BASE / "final/ALL_ACTIVE_TRAJECTORY_AUDIT.json",
        dict(
            rows=[r for c in rows for r in c["comparisons"]],
            status="NOT_EVALUABLE",
            reason="Independent all-active precision bound not yet established; descriptive errors and source arrays retained",
        ),
    )
    print(
        [
            (
                r["case"],
                [
                    (c["name"], c["groups"]["all_active"]["errors_m"]["raw"]["endpoint"])
                    for c in r["comparisons"]
                ],
            )
            for r in rows
        ],
        flush=True,
    )


if __name__ == "__main__":
    main()
