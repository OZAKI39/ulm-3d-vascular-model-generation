"""One seed, sequential single-factor comparisons; smaller steps only on quality failure."""

import argparse
import json
import h5py
import numpy as np
from .common import BASE, F1, PROTOCOL, write
from .experiment import run_case


def compare(a, b):
    with h5py.File(a / "trajectory.h5") as ha, h5py.File(b / "trajectory.h5") as hb:
        np.testing.assert_array_equal(ha["positions_m"][0], hb["positions_m"][0])
        np.testing.assert_allclose(ha["time_s"][:], hb["time_s"][:], atol=1e-14, rtol=0)
        rg_a, rg_b = ha["metrics/global_rg_m"][:, 0], hb["metrics/global_rg_m"][:, 0]
        g_a, g_b = ha["metrics/global_shape_tensor_m2"][-1, 0], hb["metrics/global_shape_tensor_m2"][-1, 0]
        e_a, e_b = ha["metrics/global_eigenvalues_m2"][-1, 0], hb["metrics/global_eigenvalues_m2"][-1, 0]
        x = ha["positions_m"][-1, 0]
        y = hb["positions_m"][-1, 0]
        x = x - x.mean(0)
        y = y - y.mean(0)
        u, _, vt = np.linalg.svd(y.T @ x)
        fix = np.eye(3)
        fix[-1, -1] = np.linalg.det(u @ vt)
        rotation = u @ fix @ vt
        result = dict(
            same_initial_positions=True,
            final_cluster_count=[
                int(ha["metrics/cluster_count"][-1, 0]),
                int(hb["metrics/cluster_count"][-1, 0]),
            ],
            final_largest_size=[
                int(ha["metrics/largest_size"][-1, 0]),
                int(hb["metrics/largest_size"][-1, 0]),
            ],
            final_cluster_sizes=[
                ha["metrics/cluster_sizes"][-1, 0].tolist(),
                hb["metrics/cluster_sizes"][-1, 0].tolist(),
            ],
            rg_series_relative_RMS=float(np.linalg.norm(rg_b - rg_a) / np.linalg.norm(rg_a)),
            final_rg_relative_difference=float(abs(rg_b[-1] - rg_a[-1]) / rg_a[-1]),
            final_shape_tensor_relative_difference=float(np.linalg.norm(g_b - g_a) / np.linalg.norm(g_a)),
            final_eigenvalues_relative_difference=float(np.linalg.norm(e_b - e_a) / np.linalg.norm(e_a)),
            final_ID_aligned_rmsd_m=float(np.sqrt(np.mean(np.sum((y @ rotation - x) ** 2, axis=1)))),
            raw_edge_changes=[len(ha["edge_events/time_s"]), len(hb["edge_events/time_s"])],
        )
    sa = json.loads((a / "SUMMARY.json").read_text())
    sb = json.loads((b / "SUMMARY.json").read_text())
    result["formation_times_s"] = [sa["formation_time_s"][0], sb["formation_time_s"][0]]
    result["numerical_quality"] = [sa["numerical_quality"], sb["numerical_quality"]]
    ta, tb = result["formation_times_s"]
    result["formation_time_relative_difference"] = None if ta is None or tb is None else abs(tb - ta) / ta
    limits = PROTOCOL["sensitivity"]
    result["statistically_model_sensitive"] = bool(
        result["final_rg_relative_difference"] > limits["rg_relative"]
        or result["final_eigenvalues_relative_difference"] > limits["shape_eigenvalues_relative"]
        or (result["formation_time_relative_difference"] or 0) > limits["formation_time_relative"]
    )
    return result


def main(label):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    models = {
        "f1": dict(force_model="f1_smooth_near_field", contact_model="f1_pbd_contact"),
        "near": dict(force_model="f2_near_field_candidate", contact_model="f1_pbd_contact"),
        "contact": dict(force_model="f2_near_field_candidate", contact_model="f2_mobility_aware_contact"),
    }
    cfg = dict(
        count=32,
        seed=20261007,
        cloud_radius_m=20e-6,
        duration_s=0.05,
        sample_interval_s=0.0005,
        events=[{"time_s": 0, "pressure_Pa": 1000}],
    )
    attempts = []
    selected = None
    for dt in [PROTOCOL["dt_s"], *PROTOCOL["fallback_dt_s"]]:
        paths = {}
        rows = {}
        for name, physics in models.items():
            if name == "f1" and dt == PROTOCOL["dt_s"]:
                paths[name] = F1 / "balanced_steps_v1/dt_20us"
                rows[name] = json.loads((paths[name] / "SUMMARY.json").read_text())
            else:
                paths[name] = out / f"{name}_{dt * 1e6:g}us"
                rows[name] = run_case(paths[name], dict(cfg, dt_s=dt, physics=physics))
        attempt = dict(
            dt_s=dt,
            paths={k: str(v) for k, v in paths.items()},
            qualities={k: v["numerical_quality"] for k, v in rows.items()},
            rows=rows,
        )
        attempts.append(attempt)
        write(out / "PARTIAL.json", dict(attempts=attempts))
        if all(v["numerical_quality"] == "PASS" for v in rows.values()):
            selected = attempt
            break
    if selected is None:
        # Honest component fallback; no change to numerical thresholds or force coefficients.
        for attempt in attempts:
            if attempt["qualities"]["near"] == "PASS":
                selected = attempt
                selected["fallback_contact"] = True
                break
        if selected is None:
            selected = attempts[0]
            selected["fallback_force_and_contact"] = True
    paths = {k: __import__("pathlib").Path(v) for k, v in selected["paths"].items()}
    force_result = compare(paths["f1"], paths["near"])
    contact_result = compare(paths["near"], paths["contact"])
    selected_name = (
        "f1"
        if selected.get("fallback_force_and_contact")
        else "near"
        if selected.get("fallback_contact")
        else "contact"
    )
    chosen = dict(
        dt_s=selected["dt_s"],
        physics=models[selected_name],
        selection_basis="Original F1 numerical quality; no morphology-target tuning",
        comparison_run=label,
        selected_name=selected_name,
    )
    off = run_case(
        out / "secondary_off",
        dict(cfg, dt_s=chosen["dt_s"], physics=dict(chosen["physics"], secondary=False)),
    )
    with h5py.File(out / "secondary_off/trajectory.h5") as h:
        drift = float(np.max(np.linalg.norm(h["positions_m"][:] - h["positions_m"][0], axis=-1)))
    assert drift == 0
    write(out / "SELECTED.json", chosen)
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            attempts=attempts,
            selected=chosen,
            near_field=force_result,
            contact=contact_result,
            secondary_off=off,
            secondary_off_max_displacement_m=drift,
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    main(a.label)
