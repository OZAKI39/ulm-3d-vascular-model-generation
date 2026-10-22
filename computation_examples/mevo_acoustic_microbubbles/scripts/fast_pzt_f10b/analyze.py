"""Offline shadow comparison only; never supplies a predictor with future states."""

import h5py
import numpy as np
from .common import ROOT, BASE, F10A, js, sha, write, identity

KEYS = ["OFF", "S_1000", "O_1000"]


def read_proxy(mode, key):
    path = BASE / "rollout_v1" / (mode + "_" + key) / "centroid.h5"
    with h5py.File(path) as h:
        return h["time_s"][:], h["centroid_m"][:]


def read_reference(key, phase="extended_5_v1"):
    with h5py.File(F10A / phase / key / "trajectory.h5") as h:
        return h["time_s"][:], h["metrics/global_centroid_m"][:, 0]


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    proto = js(ROOT / "configs/fast_pzt/f10b/protocol.json")
    limits = proto["error"]
    results = {}
    curves = {}
    sources = {}
    for mode in ["G0", "G1"]:
        rows = {}
        increments = {}
        for key in KEYS:
            t, c = read_proxy(mode, key)
            tr, r = read_reference(key)
            np.testing.assert_allclose(t, tr, rtol=0, atol=2e-15)
            _, coarse = read_reference(key, "extended_10_v1")
            diff = c - r
            row = dict(
                endpoint_error_m=diff[-1],
                endpoint_Z_error_m=float(abs(diff[-1, 2])),
                maximum_Z_error_m=float(abs(diff[:, 2]).max()),
                endpoint_3D_error_m=float(np.linalg.norm(diff[-1])),
                maximum_3D_error_m=float(np.linalg.norm(diff, axis=1).max()),
                prediction_displacement_m=c[-1] - c[0],
                reference_displacement_m=r[-1] - r[0],
                coarse_maximum_Z_error_m=float(abs(c[:, 2] - coarse[:, 2]).max()),
            )
            row["accuracy"] = all(
                row[a] <= limits[b]
                for a, b in [
                    ("endpoint_Z_error_m", "endpoint_Z_m"),
                    ("maximum_Z_error_m", "max_Z_m"),
                    ("endpoint_3D_error_m", "endpoint_3D_m"),
                ]
            )
            rows[key] = row
            curves[mode + "_" + key + "_error_m"] = diff
            curves[mode + "_" + key + "_centroid_m"] = c
            curves["reference_" + key + "_centroid_m"] = r
            curves["time_s"] = t
            for path in [
                BASE / "rollout_v1" / (mode + "_" + key) / "centroid.h5",
                F10A / "extended_5_v1" / key / "trajectory.h5",
                F10A / "extended_10_v1" / key / "trajectory.h5",
            ]:
                sources[str(path.relative_to(ROOT))] = sha(path)
        for comparator in ["OFF", "S_1000"]:
            a = rows["O_1000"]
            b = rows[comparator]
            pred = float(a["prediction_displacement_m"][2] - b["prediction_displacement_m"][2])
            ref = float(a["reference_displacement_m"][2] - b["reference_displacement_m"][2])
            increments["O_minus_" + comparator] = dict(
                predicted_m=pred,
                reference_m=ref,
                error_m=pred - ref,
                relative_error=(pred - ref) / ref,
                direction_correct=pred > 0 and ref > 0,
            )
        results[mode] = dict(
            cases=rows,
            increments=increments,
            surrogate_numerical_admission="PASS"
            if all(row["accuracy"] for row in rows.values())
            and all(v["direction_correct"] for v in increments.values())
            else "FAIL",
        )
    timing = js(BASE / "rollout_v1/TIMING.json")
    full = timing["rows"][0]
    for mode in ["G0", "G1"]:
        proxy = next(r for r in timing["rows"] if r["mode"] == mode and r["control"] == "O_1000")
        speed = full["body_s"] / proxy["body_s"]
        results[mode]["performance"] = dict(
            body_speedup=speed,
            subprocess_speedup=full["subprocess_wall_s"] / proxy["subprocess_wall_s"],
            integration_speedup=full["integration_s"] / proxy["integration_s"],
            body_gate="PASS" if speed >= proto["performance"]["min_body_speedup"] else "FAIL",
            engineering_value_gate="PASS"
            if speed >= 2 and results[mode]["surrogate_numerical_admission"] == "PASS"
            else "FAIL",
            full_body_s=full["body_s"],
            surrogate_body_s=proxy["body_s"],
            replicates=1,
        )
    contacts = []
    for folder, labels in [("contact_v1", ["A0", "A1"]), ("contact_v2", ["A2", "A3"])]:
        for label in labels:
            d = js(BASE / folder / (label + ".json"))
            trace = d["traces"]
            last = trace[-1]
            contacts.append(
                dict(
                    case=label,
                    dt_s=d["dt_s"],
                    rounds=d["contact_rounds"],
                    accepted_duration_s=d["accepted_duration_s"],
                    first_rejection_s=d["first_rejection_s"],
                    classification=d["classification"],
                    max_predicted_overlap_m=max(t["prediction"]["max_overlap_m"] for t in trace),
                    max_final_overlap_m=max(t["final"]["max_overlap_m"] for t in trace),
                    max_pass_correction_m=max(t["quality"]["max_contact_pass_correction_m"] for t in trace),
                    last_error_bits=last["quality"]["error_bits"],
                    failure_reason="residual_overlap>.01R"
                    if last["quality"]["error_bits"] & 2048
                    else "single_pass_correction>.1R"
                    if d["first_rejection_s"]
                    else None,
                    source=f"{folder}/{label}.json",
                    last_free_speed_m_s=np.linalg.norm(last["free_velocity_m_s"], axis=1).max(),
                    last_secondary_rms_N=np.sqrt(
                        np.mean(np.sum(np.array(last["secondary_force_N"]) ** 2, axis=1))
                    ),
                    physics="PHYSICAL_MODEL_UNCALIBRATED",
                )
            )
    # Only one new full trajectory: require exact reproduction of archived10us states.
    positions = np.load(BASE / "rollout_v1/FULL_O_1000/reference_validation_positions.npz")["positions_m"]
    with h5py.File(F10A / "extended_10_v1/O_1000/trajectory.h5") as h:
        np.testing.assert_array_equal(positions, h["positions_m"][:, 0])
    np.savez_compressed(out / "curves.npz", **curves)
    write(
        out / "RESULTS.json",
        dict(
            models=results,
            contacts=contacts,
            thresholds=limits,
            full_reference_reproduction="BITWISE_IDENTICAL_ALL201_SAVED_STATES",
            one_step_status=js(BASE / "checks_v1/RESULTS.json")["status"],
            support_status=js(BASE / "support_v1/RESULTS.json")["status"],
            timing=timing,
            training_ready=False,
            physics_calibrated=False,
            morphology_output=False,
            engineering_scope="All requested bounded experiments completed; accuracy and speed gates reported separately; no batch campaign after standalone2x gate failed",
            failure_record="contact_v1 A2 construction refused F3 frozen initializer before integrating; A0/A1 preserved, only missingA2/A3 run incontact_v2",
        ),
    )
    write(out / "SOURCES.json", dict(identity=identity(), files=sources))
    print(
        {
            m: (r["surrogate_numerical_admission"], r["performance"]["body_speedup"])
            for m, r in results.items()
        }
    )


if __name__ == "__main__":
    main()
