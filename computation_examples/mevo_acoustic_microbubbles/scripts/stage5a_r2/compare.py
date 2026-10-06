"""Same-ID comparisons on actual states, with absolute and sampling-age bounds."""

from pathlib import Path
import json
import h5py
import numpy as np
from scripts.stage5a.compare import sample
from scripts.stage4l.dense_probe import sha, write


def stitch(paths, out):
    times = []
    positions = []
    active = []
    identity = None
    for path in paths:
        with h5py.File(path) as h:
            t = h["time_s"][:]
            x = h["positions_m"][:]
            m = h["active"][:]
            current = (h["particle_id"][:], h["radii_m"][:])
            if identity is None:
                identity = current
            else:
                assert all(np.array_equal(a, b) for a, b in zip(identity, current))
                assert abs(t[0] - times[-1][-1]) <= 1e-15
                assert np.array_equal(x[0], positions[-1][-1], equal_nan=True)
                t, x, m = t[1:], x[1:], m[1:]
            times.append(t)
            positions.append(x)
            active.append(m)
    with h5py.File(out, "x") as h:
        h["time_s"] = np.concatenate(times)
        h.create_dataset("positions_m", data=np.concatenate(positions), compression="gzip")
        h.create_dataset("active", data=np.concatenate(active), compression="gzip")
        h["particle_id"] = identity[0]
        h["radii_m"] = identity[1]
        h.attrs["segments"] = json.dumps([dict(path=str(p), sha256=sha(p)) for p in paths])
        h.attrs["coordinate_interpolation"] = False


def actual_sample(path, clock, members=None):
    # Allow only rounding-sized clock error, never advance to a future physical
    # frame farther than 16 eps of the absolute timestamp.
    with h5py.File(path) as h:
        t = h["time_s"][:]
    chosen = []
    for at in clock:
        near = int(np.argmin(abs(t - at)))
        chosen.append(t[near] if abs(t[near] - at) <= 16 * np.finfo(float).eps * max(1, abs(at)) else at)
    data = sample(path, np.array(chosen), members)
    data["requested_time_s"] = clock
    data["roundoff_clock_adjustment_s"] = np.array(chosen) - clock
    data["sample_age_s"] = clock - data["actual_time_s"]
    return data


def compare(reference, candidate, out, label):
    out = Path(out)
    out.mkdir(exist_ok=True, parents=True)
    with h5py.File(reference) as h, h5py.File(candidate) as k:
        t0 = max(float(h["time_s"][0]), float(k["time_s"][0]))
        end = min(float(h["time_s"][-1]), float(k["time_s"][-1]))
    clock = t0 + np.arange(int(np.floor((end - t0) / 1e-5 + 1e-7)) + 1) * 1e-5
    if end - t0 > 1e-15 and end - clock[-1] > 1e-14:
        clock = np.r_[clock, end]
    a = actual_sample(reference, clock)
    b = actual_sample(candidate, clock, a["particle_ids"])
    meaningful = len(clock) > 1
    rg_end = float(abs(b["rg_m"][-1] / a["rg_m"][-1] - 1))
    rg_l2 = float(np.linalg.norm(b["rg_m"] - a["rg_m"]) / np.linalg.norm(a["rg_m"]))
    quart = abs(b["pair_quantiles_m"][-1, 1:4] / a["pair_quantiles_m"][-1, 1:4] - 1)
    travel = float(np.linalg.norm(np.diff(a["centroid_m"], axis=0), axis=1).sum())
    centroid = float(np.linalg.norm(b["centroid_m"][-1] - a["centroid_m"][-1]))
    budget = max(0.1 * float(a["radius_m"]), 0.02 * travel)

    # Bound the reference sampling uncertainty with actual next accepted-frame
    # displacements along its recorded linear path; do not interpolate coordinates.
    def sampling_bound(path, data):
        with h5py.File(path) as h:
            ts = h["time_s"][:]
            ids = h["particle_id"][:]
            idx = np.array([np.flatnonzero(ids == i)[0] for i in a["particle_ids"]])
            bounds = []
            for time_actual, time_requested in zip(data["actual_time_s"], clock):
                j = int(np.searchsorted(ts, time_actual))
                j = min(j, len(ts) - 1)
                if j == len(ts) - 1 or time_requested - time_actual < 1e-14:
                    bounds.append(0.0)
                    continue
                displacement = h["positions_m"][j + 1][idx] - h["positions_m"][j][idx]
                bounds.append(float(np.linalg.norm(displacement, axis=1).max()))
        return np.array(bounds)

    ab = sampling_bound(reference, a)
    bb = sampling_bound(candidate, b)
    uncertainty = float(ab[-1] + bb[-1])
    rg_series_uncertainty = float(np.linalg.norm(ab + bb) / np.linalg.norm(a["rg_m"]))
    passed = meaningful and rg_end <= 0.02 and rg_l2 <= 0.02 and max(quart) <= 0.02 and centroid <= budget
    # Interval uncertainty cannot silently decide a marginal gate.
    sampling_ambiguous = (
        meaningful
        and (uncertainty > 0 or rg_series_uncertainty > 0)
        and (
            abs(centroid - budget) <= uncertainty
            or abs(rg_end - 0.02) <= uncertainty / a["rg_m"][-1]
            or np.any(abs(quart - 0.02) <= 2 * uncertainty / a["pair_quantiles_m"][-1, 1:4])
            or abs(rg_l2 - 0.02) <= rg_series_uncertainty
        )
    )
    for key, data in [("reference", a), ("candidate", b)]:
        np.savez_compressed(out / f"{label}_{key}.npz", **data)
    report = dict(
        label=label,
        reference=str(reference),
        candidate=str(candidate),
        reference_sha256=sha(reference),
        candidate_sha256=sha(candidate),
        compared_duration_s=end - t0,
        status="INCOMPLETE_SAMPLING_UNCERTAINTY"
        if sampling_ambiguous
        else ("PASS" if passed else "FAIL")
        if meaningful
        else "NOT_EVALUABLE",
        pass_scientific=bool(passed and not sampling_ambiguous),
        Rg_endpoint_relative=rg_end if meaningful else None,
        Rg_series_relative_L2=rg_l2 if meaningful else None,
        pair_quartiles_relative=quart.tolist() if meaningful else None,
        centroid_difference_m=centroid if meaningful else None,
        centroid_budget_m=budget,
        Rg_endpoint_absolute_m=float(abs(b["rg_m"][-1] - a["rg_m"][-1])) if meaningful else None,
        position_endpoint_max_difference_m=float(
            np.linalg.norm(b["positions_m"][-1] - a["positions_m"][-1], axis=1).max()
        )
        if meaningful
        else None,
        pair_quartiles_absolute_m=abs(
            b["pair_quantiles_m"][-1, 1:4] - a["pair_quantiles_m"][-1, 1:4]
        ).tolist()
        if meaningful
        else None,
        same_ID_particles=a["particle_ids"].tolist(),
        reference_max_sample_age_s=float(a["sample_age_s"].max()),
        candidate_max_sample_age_s=float(b["sample_age_s"].max()),
        reference_sampling_displacement_bound_max_m=float(ab.max()),
        candidate_sampling_displacement_bound_max_m=float(bb.max()),
        endpoint_sampling_uncertainty_m=uncertainty,
        Rg_series_sampling_relative_L2_uncertainty=rg_series_uncertainty,
        rigid_fit_reference_end_m=float(a["rigid_fit_m"][-1]),
        rigid_fit_candidate_end_m=float(b["rigid_fit_m"][-1]),
        graph_change_mean_reference=float(a["neighbor_turnover"].mean()) if meaningful else None,
        graph_change_mean_candidate=float(b["neighbor_turnover"].mean()) if meaningful else None,
        hard_contact_chatter_is_sole_failure_criterion=False,
        coordinate_interpolation=False,
    )
    write(out / f"{label}.json", report)
    return report
