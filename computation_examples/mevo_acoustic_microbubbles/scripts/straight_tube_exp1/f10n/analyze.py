"""Fixed-ID effects separated from translation, replay, timestep and cache error."""

import numpy as np
import h5py
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from .common import BASE, configs, read, write, sha, identity


def load_run(out, name):
    p = out / "runs" / name
    with h5py.File(p / "trajectory.h5") as h:
        arrays = {
            k: h[k][:]
            for k in [
                "positions_m",
                "time_s",
                "particle_id",
                "pressure_Pa",
                "primary_force_N",
                "secondary_force_N",
            ]
        }
    return read(p / "RESULT.json"), arrays, dict(np.load(p / "metrics.npz"))


def compare(out, a, b):
    ra, aa, ma = load_run(out, a)
    rb, ab, mb = load_run(out, b)
    if not (ra["eligible"] and rb["eligible"]):
        return dict(status="NOT_COMPARABLE_FAILED_TRAJECTORY")
    np.testing.assert_array_equal(aa["particle_id"], ab["particle_id"])
    np.testing.assert_allclose(aa["time_s"], ab["time_s"], rtol=0, atol=1e-14)
    xa, xb = aa["positions_m"][:, 0], ab["positions_m"][:, 0]
    ca, cb = xa.mean(1), xb.mean(1)
    rela = xa - ca[:, None]
    relb = xb - cb[:, None]
    return dict(
        status="PASS",
        max_fixed_id_position_difference_m=float(np.linalg.norm(xa - xb, axis=-1).max()),
        max_COM_difference_m=float(np.linalg.norm(ca - cb, axis=-1).max()),
        max_centered_fixed_id_difference_m=float(np.linalg.norm(rela - relb, axis=-1).max()),
        max_Rg_difference_m=float(np.max(abs(ma["shape_rg_m"] - mb["shape_rg_m"]))),
        max_directional_xy_difference=float(
            np.max(abs(ma["shape_directional_xy"] - mb["shape_directional_xy"]))
        ),
        max_aspect_difference=float(np.max(abs(ma["shape_aspect"] - mb["shape_aspect"]))),
        max_tensor_difference_m2=float(
            np.linalg.norm(ma["shape_tensor_m2"] - mb["shape_tensor_m2"], axis=(-2, -1)).max()
        ),
        cluster_count_identical=bool(np.array_equal(ma["cluster_count"], mb["cluster_count"])),
        scope="observed sensitivity over registered short run; not a proof of asymptotic convergence",
    )


def main():
    p = configs()["protocol"]
    accepted = read(BASE / "ACCEPTED_CAMPAIGN.json")
    out = BASE / accepted["path"]
    assert sha(out / "CAMPAIGN.json") == accepted["sha256"]
    nav = read(out / "NAV_SELECTION.json")["selected"]
    shape = read(out / "SHAPE_SELECTION.json")["selected"]
    cmp = {
        "navigation_replay": compare(out, nav, "OPTIMIZED_REPLAY"),
        "navigation_dt": compare(out, "OPTIMIZED_REPLAY", "OPTIMIZED_HALF_DT"),
        "navigation_cache": compare(out, "OPTIMIZED_REPLAY", "OPTIMIZED_FINE_CACHE"),
        "baseline_dt": compare(out, "BASELINE", "BASELINE_HALF_DT"),
        "shape_replay": compare(out, shape, "SHAPE_REPLAY"),
        "shape_dt": compare(out, "SHAPE_REPLAY", "SHAPE_HALF_DT"),
        "shape_cache": compare(out, "SHAPE_REPLAY", "SHAPE_FINE_CACHE"),
    }
    write(BASE / "analysis/NUMERICAL_SENSITIVITY.json", cmp)
    off, ao, mo = load_run(out, "OFF")
    base, ab, mb = load_run(out, "BASELINE")
    opt, an, mn = load_run(out, "OPTIMIZED_REPLAY")
    sh, ass, ms = load_run(out, "SHAPE_REPLAY")
    numerical = max(
        cmp["navigation_dt"].get("max_COM_difference_m", np.inf),
        cmp["navigation_cache"].get("max_COM_difference_m", np.inf),
        cmp["navigation_replay"].get("max_COM_difference_m", np.inf),
    )
    particle_noise = max(
        cmp["navigation_dt"].get("max_fixed_id_position_difference_m", np.inf),
        cmp["navigation_cache"].get("max_fixed_id_position_difference_m", np.inf),
    )
    shape_noise = max(
        cmp["shape_dt"].get("max_directional_xy_difference", np.inf),
        cmp["shape_cache"].get("max_directional_xy_difference", np.inf),
        cmp["shape_replay"].get("max_directional_xy_difference", np.inf),
    )
    rg_noise = max(
        cmp["shape_dt"].get("max_Rg_difference_m", np.inf),
        cmp["shape_cache"].get("max_Rg_difference_m", np.inf),
    )
    table = []
    for name in [
        "OFF",
        "BASELINE",
        "OPTIMIZED_REPLAY",
        "SHAPE_REPLAY",
        "DISPERSED_OFF",
        "DISPERSED_OPTIMIZED",
    ]:
        rec, a, m = load_run(out, name)
        x = a["positions_m"][:, 0]
        center = x.mean(1)
        rel = x - center[:, None]
        gap = squareform(pdist(x[-1])) - 2 * p["radius_m"]
        np.fill_diagonal(gap, np.inf)
        count = {
            str(f): int(connected_components(gap <= f * p["radius_m"], directed=False)[0])
            for f in [0.02, 0.05, 0.1, 0.2]
        }
        table.append(
            dict(
                name=name,
                endpoint_error_um=rec["endpoint_error_m"] * 1e6,
                displacement_um=np.array(rec["displacement_m"]) * 1e6,
                displacement_norm_um=float(np.linalg.norm(rec["displacement_m"]) * 1e6),
                rg_um=m["shape_rg_m"][[0, -1]] * 1e6,
                aspect=m["shape_aspect"][[0, -1]],
                directional_xy=m["shape_directional_xy"][[0, -1]],
                axis_initial_final=m["shape_principal_axis"][[0, -1]],
                axis_identifiable=m["shape_axis_identifiable"][[0, -1]],
                angle_to_x_deg=m["shape_angle_to_x_deg"][[0, -1]],
                tensor_initial_final_m2=m["shape_tensor_m2"][[0, -1]],
                max_centered_fixed_id_displacement_um=float(
                    np.linalg.norm(rel[-1] - rel[0], axis=-1).max() * 1e6
                ),
                clusters=int(m["cluster_count"][-1, 0]),
                cluster_threshold_sensitivity=count,
                actual_particle_pressure_range_Pa=[
                    float(np.abs(a["pressure_Pa"]).min()),
                    float(np.abs(a["pressure_Pa"]).max()),
                ],
                peak_primary_force_N=float(np.linalg.norm(a["primary_force_N"], axis=-1).max()),
                peak_secondary_force_N=float(np.linalg.norm(a["secondary_force_N"], axis=-1).max()),
                mean_primary_force_initial_N=a["primary_force_N"][0].mean(0),
                quality=rec["quality"],
            )
        )
    improvement = off["endpoint_error_m"] - opt["endpoint_error_m"]
    xy_change = sh["directional_xy_final"] - sh["directional_xy_initial"]
    rg_change = sh["rg_final_m"] - sh["rg_initial_m"]
    gate = dict(
        engineering_status="PASS",
        position_target_m=p["target_position_m"],
        target_unchanged=True,
        target_reached=bool(opt["eligible"] and opt["endpoint_error_m"] <= p["target_tolerance_m"]),
        optimized_control=nav,
        shape_control=shape,
        initial_error_um=off["initial_error_m"] * 1e6,
        optimized_error_um=opt["endpoint_error_m"] * 1e6,
        error_reduction_um=improvement * 1e6,
        acoustic_extra_displacement_um=(np.array(opt["displacement_m"]) - off["displacement_m"]) * 1e6,
        position_effect_exceeds_tested_numerical_difference=bool(improvement > 3 * numerical),
        navigation_numerical_COM_envelope_um=numerical * 1e6,
        navigation_fixed_id_numerical_envelope_um=particle_noise * 1e6,
        navigation_particle_resolution_pass=bool(particle_noise <= p["position_dt_tolerance_m"]),
        navigation_particle_resolution_tolerance_um=p["position_dt_tolerance_m"] * 1e6,
        navigation_morphology_status="UNRESOLVED_SAME_ID_DT_AND_CACHE_SENSITIVITY"
        if particle_noise > p["position_dt_tolerance_m"]
        else "WITHIN_TESTED_TOLERANCE",
        shape_directional_xy_change=xy_change,
        shape_directional_xy_numerical_envelope=shape_noise,
        shape_directional_effect_exceeds_tested_numerical_difference=bool(
            abs(xy_change) > p["shape_significance_factor"] * shape_noise
        ),
        shape_target_reached=bool(
            sh["eligible"]
            and sh["directional_xy_final"] >= p["target_shape"]["directional_aspect_x_over_y"]
            and ms["shape_axis_identifiable"][-1]
            and ms["shape_angle_to_x_deg"][-1] <= p["target_shape"]["angle_tolerance_deg"]
        ),
        Rg_change_um=rg_change * 1e6,
        Rg_numerical_envelope_um=rg_noise * 1e6,
        Rg_effect_exceeds_tested_numerical_difference=bool(abs(rg_change) > 3 * rg_noise),
        connected=opt["cluster_count_max"] == 1,
        Z_scope="local symmetric prototype: no independently demonstrated signed Z translation",
        calibration="EXPLORATORY / PHYSICS_UNCALIBRATED / HOMOGENEOUS_WATER_CONTROL_PROXY",
        training_ready=False,
        experiment_matched=False,
        morphology_interpretation="fixed-ID geometry, not fusion; instantaneous graph events and dense contact remain uncalibrated",
        identity=identity(),
    )
    write(BASE / "analysis/CONTROL_COMPARISON.json", table)
    write(BASE / "analysis/FINAL_GATE.json", gate)
    print(gate)


if __name__ == "__main__":
    main()
