"""Accepted fixed-ID trajectories and explicit rejected-attempt evidence."""

import copy
import time
import numpy as np
import h5py
from .common import BASE, ROOT, js, write, sha, identity
from .validation import validate_ledger


def read_run(path):
    s = js(path / "SUMMARY.json")
    with h5py.File(path / "trajectory.h5") as h:
        t = h["time_s"][:]
        x = h["positions_m"][:, 0]
        c = x.mean(1)
        d = h["motion_ledger/displacement_m"][:]
        m = h["motion_ledger/step_meta"][:]
        np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
        dist = np.linalg.norm(x[:, :, None] - x[:, None, :], axis=-1)
        dist[:, np.arange(32), np.arange(32)] = np.inf
        pressure = abs(h["control/pressure_Pa"][:])
        fraction = h["metrics/largest_size"][:, 0] / 32
        rg = h["metrics/global_rg_m"][:, 0]
        row = dict(
            status=s["status"],
            quality=s["numerical_quality"],
            reason=s.get("reason"),
            accepted_time_s=t[-1],
            dt_s=js(path / "CONFIG.json")["dt_s"],
            displacement_m=c[-1] - c[0],
            minimum_saved_pair_distance_m=dist.min(),
            last_largest_fraction=fraction[-1],
            last_Rg_m=rg[-1],
            pressure_min_Pa=pressure.min(),
            pressure_max_Pa=pressure.max(),
            ledger_signed_displacement_m=d.sum(0),
            mean_particle_contact_path_m=m[:, 4].sum(),
            mean_particle_wall_path_m=m[:, 5].sum(),
            quality_raw=s["quality"],
            ledger_validation=s["ledger_validation"],
            rejected_intervals=s["rejected_audit_intervals"],
            wall_s=s["wall_s"],
            process_body_wall_s=s["process_body_wall_s"],
            output_s=s["output_and_postprocessing_s"],
            execution_commit=js(path / "IDENTITY.json")["execution_commit"],
            trajectory_sha256=sha(path / "trajectory.h5"),
        )
        arrays = dict(
            time_s=t,
            centroid_m=c,
            positions_m=x,
            largest_fraction=fraction,
            Rg_m=rg,
            ledger=d,
            meta=m,
            minimum_pair_distance_m=dist.min((1, 2)),
        )
    return row, arrays


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    start = time.perf_counter()
    rows = {}
    arrays = {}
    for campaign in ["pilot_v3", "extended_10_v1", "extended_5_v1"]:
        for path in sorted((BASE / campaign).glob("*/SUMMARY.json")):
            key = campaign + "/" + path.parent.name
            rows[key], arrays[key] = read_run(path.parent)
    effects = {}
    steps = {}
    for cond in ["OFF", "S_1000", "O_1000"]:
        co = arrays["extended_10_v1/" + cond]
        fi = arrays["extended_5_v1/" + cond]
        np.testing.assert_allclose(co["time_s"], fi["time_s"], atol=2e-15)
        diff = co["centroid_m"] - fi["centroid_m"]
        steps[cond] = dict(
            max_abs_Z_difference_m=abs(diff[:, 2]).max(),
            final_Z_difference_m=diff[-1, 2],
            max_centroid_norm_difference_m=np.linalg.norm(diff, axis=1).max(),
            max_particle_norm_difference_m=np.linalg.norm(
                co["positions_m"] - fi["positions_m"], axis=2
            ).max(),
            max_Rg_relative_difference=np.max(abs(co["Rg_m"] / fi["Rg_m"] - 1)),
        )
    for campaign in ["extended_10_v1", "extended_5_v1"]:
        c = {key: arrays[campaign + "/" + key]["centroid_m"] for key in ["OFF", "S_1000", "O_1000"]}
        effects[campaign] = {other: float((c["O_1000"][-1] - c[other][-1])[2]) for other in ["OFF", "S_1000"]}
    limits = {
        other: max(
            5e-8, 3 * (steps["O_1000"]["max_abs_Z_difference_m"] + steps[other]["max_abs_Z_difference_m"])
        )
        for other in ["OFF", "S_1000"]
    }
    passed = all(v["quality"] == "PASS" for k, v in rows.items() if k.startswith("extended_"))
    resolved = bool(
        passed and all(effects[cam][other] > limits[other] for cam in effects for other in limits)
    )
    key = rows["extended_5_v1/O_1000"]
    effect = effects["extended_5_v1"]["OFF"]
    contact = key["mean_particle_contact_path_m"]
    wall = key["ledger_signed_displacement_m"][5][2]
    flags = dict(
        RADIAL_NUMERICAL_PASS=all(js(BASE / "radial_v1/RESULTS.json")["pressure_admission"].values()),
        SHELL_BRANCH_CROSSED=False,
        shell_branch_scope="production sigma=.036 at four tested pressures; auxiliary sigma0/.072 cases do cross branches",
        DIRECTION_FEASIBLE=True,
        HIGH_DRIVE_EXPLORATORY=True,
        ACOUSTIC_EFFECT_RESOLVED=resolved,
        CONTACT_DOMINATED=bool(contact >= abs(effect)),
        WALL_DOMINATED=bool(abs(wall) >= 0.25 * abs(effect) or abs(wall) > 5e-8),
        PHYSICS_UNCALIBRATED=True,
        NAVIGATION_NOT_ESTABLISHED=True,
        TRAINING_READY=False,
    )
    # Use accepted actual data to challenge the ledger validator; both corruptions must be rejected.
    path = BASE / "extended_10_v1/O_1000/trajectory.h5"
    with h5py.File(path) as h:
        ledger = {
            name: h["motion_ledger/" + name][:]
            for name in ["displacement_m", "step_meta", "wall_contact_events"]
        }
        ledger["rejected_intervals"] = []
        x = h["positions_m"][:]
        t = h["time_s"][:]
    validate_ledger(ledger, x, t)
    tests = {}
    for name in ["changed_position", "changed_term"]:
        q = copy.deepcopy(ledger)
        xx = x.copy()
        if name == "changed_position":
            xx[-1, 0, :, 2] += 1e-9
        else:
            q["displacement_m"][0, 1, 2] += 1e-9
        try:
            validate_ledger(q, xx, t)
            tests[name] = "UNEXPECTED_ACCEPT"
        except AssertionError:
            tests[name] = "REJECTED_AS_REQUIRED"
    assert all(v == "REJECTED_AS_REQUIRED" for v in tests.values())
    rad = js(BASE / "radial_v1/RESULTS.json")
    prod = [r for r in rad["rows"] if r["name"].startswith("production_") and r["samples_per_cycle"] == 512]
    reference = prod[0]
    scaling = []
    for r in prod:
        v = np.array(r["harmonics_real_m3"]) + 1j * np.array(r["harmonics_imag_m3"])
        v0 = np.array(reference["harmonics_real_m3"]) + 1j * np.array(reference["harmonics_imag_m3"])
        ratio = r["pressure_Pa"] / reference["pressure_Pa"]
        scaling.append(
            dict(
                pressure_Pa=r["pressure_Pa"],
                complex_primary_response_over_quadratic=abs(v[0] / v0[0]) / ratio,
                standing_gradient_primary_over_quadratic=v[0].real / v0[0].real / ratio,
                quadrature_gradient_primary_over_quadratic=v[0].imag / v0[0].imag / ratio,
                equal_in_phase_far_pair_over_quadratic=float(
                    np.sum((np.arange(1, 9) * abs(v)) ** 2)
                    / np.sum((np.arange(1, 9) * abs(v0)) ** 2)
                    / ratio**2
                ),
            )
        )
    write(
        out / "RESULTS.json",
        dict(
            identity=identity(),
            flags=flags,
            runs=rows,
            effects=effects,
            resolution_threshold_m=limits,
            timestep=steps,
            direct_radial_scaling=scaling,
            ledger_corruption_checks=tests,
            short_case="10ms interior location; fixed source commands, initial median1kPa; pressure varies along true transport",
            contact_interpretation="Large unsigned individual contact path changes internal packing; central PBD signed group correction remains near zero. Do not subtract unsigned path from centroid motion.",
            static_mobility_limit="No full near-wall RPY calibration; no production RPY trajectories",
            analysis_wall_s=time.perf_counter() - start,
        ),
    )
    write(
        out / "SOURCE_MANIFEST.json",
        {
            str(p.relative_to(ROOT)): dict(sha256=sha(p), bytes=p.stat().st_size)
            for p in BASE.rglob("*")
            if p.is_file()
            and p.suffix in [".json", ".h5", ".npz"]
            and "analysis" not in p.parts
            and "figures" not in p.parts
        },
    )
    print("F10A flags", flags, "paired effects", effects, "sensitivity", steps, flush=True)


if __name__ == "__main__":
    main()
