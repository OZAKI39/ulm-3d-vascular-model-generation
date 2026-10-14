"""Read-only postprocessing of every accepted F10Q force-evaluation state.

No new integration; exact saved positions, original cubic cache/APECSS/F2.
Pair-attraction diagnostics are checked against original saved GPU force sums.
The archived F10Q delivery and its manifest are never rewritten.
"""

import csv
import time
import subprocess
import numpy as np
import h5py
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics, pair_reference, feedback_weight
from .common import ROOT, BASE, OLD, read, write, sha

OUT = ROOT / "results/straight_tube_exp1/f10q_parameter_summary_20261014"


def pairs(x, vol, p, physics, chi_zero):
    """Batched algebra of unchanged F2 pair_reference; positive scalar attracts."""
    i, j = np.triu_indices(32, 1)
    delta = x[:, j] - x[:, i]
    d = np.linalg.norm(delta, axis=-1)
    vi, vj = vol[:, i], vol[:, j]
    pi, pj = p[:, i], p[:, j]
    ci = np.full(pi.shape, chi_zero, dtype=complex)
    cj = ci.copy()
    np.divide(vi[:, :, 0], pi, out=ci, where=abs(pi) > 0)
    np.divide(vj[:, :, 0], pj, out=cj, where=abs(pj) > 0)
    k = physics.density_kg_m3 * (2 * np.pi * physics.frequency_Hz) ** 2 / (4 * np.pi * d)
    a, b = k * ci, k * cj
    det = 1 - a * b
    assert abs(det).min() >= 1e-6
    ui = (vi[:, :, 0] - a * vj[:, :, 0]) / det
    uj = (vj[:, :, 0] - b * vi[:, :, 0]) / det
    w = feedback_weight(d, physics.radius_m, physics.feedback_end_over_sum)
    corr = np.sum(np.arange(1, 9) ** 2 * np.real(vi * vj.conj()), axis=-1)
    corr += w * (np.real(ui * uj.conj()) - np.real(vi[:, :, 0] * vj[:, :, 0].conj()))
    q = np.minimum(d / (physics.transition_over_sum * 2 * physics.radius_m), 1)
    core = physics.soft_core_over_sum * 2 * physics.radius_m * (1 - q) ** 2
    scalar = (
        physics.density_kg_m3
        * (2 * np.pi * physics.frequency_Hz) ** 2
        / (8 * np.pi)
        * corr
        * d
        / (d * d + core * core) ** 1.5
    )
    vector = scalar[:, :, None] * delta / d[:, :, None]
    net = np.zeros_like(x)
    per_bubble = np.zeros(x.shape[:2])
    partner = np.full(x.shape[:2], -1, dtype=int)
    for pair_index, (ii, jj) in enumerate(zip(i, j)):
        net[:, ii] += vector[:, pair_index]
        net[:, jj] -= vector[:, pair_index]
        val = np.maximum(scalar[:, pair_index], 0)
        for left, right in [(ii, jj), (jj, ii)]:
            better = val > per_bubble[:, left]
            per_bubble[better, left] = val[better]
            partner[better, left] = right
    effective = np.r_[abs(pi - k * uj)[w > 0], abs(pj - k * ui)[w > 0]]
    return net, scalar, per_bubble, partner, effective


def main():
    start = time.perf_counter()
    OUT.mkdir(exist_ok=False)
    physics = F2Physics()
    basis = CachedWaterBasis(OLD / "cache/coarse.h5")
    response = ResponseSection(OLD / "inputs/response.h5", physics)
    chi = response.query(np.array([response.knots[0] + 0j]))[0, 0] / response.knots[0]
    geometry = read(OLD / "CONFIG_SNAPSHOT.json")["geometry"]
    all_rows = []
    reports = {}
    hashes = {}
    paths = {
        "OFF": BASE / "b_navigation_v1/runs/OFF_DT1",
        "X": BASE / "b_navigation_v1/runs/NAV_X_DT1",
        "Y": BASE / "b_navigation_v1/runs/NAV_Y_DT1",
        "OLD_STRONG_VALID_PREFIX": BASE / "a_strong_v1/runs/STRONG_DT1",
    }
    for tag, folder in paths.items():
        started = time.perf_counter()
        e = read(folder / "INPUT.json")["control"]
        u = basis.vector(e["amplitude_Pa"], e["phase_rad"], e["enabled"])
        result = read(folder / "RESULT.json")
        with h5py.File(folder / "trajectory.h5") as h:
            xs = h["positions_m"][:]
            times = h["time_s"][:]
            extra = h["extra_metrics"][:]
            sample = {k: h["force_samples"][k][:] for k in h["force_samples"]}
        states = xs[:-1]
        n = len(states)
        forces = []
        sec = []
        pressures = []
        max_pair = []
        partners = []
        pairmax = []
        pairid = []
        effective = [np.inf, 0.0]
        i, j = np.triu_indices(32, 1)
        for first in range(0, n, 128):
            x = states[first : first + 128]
            pb, gb = basis.query_basis(x.reshape(-1, 3))
            p = (u @ pb).reshape(-1, 32)
            g = np.einsum("i,inj->nj", u, gb).reshape(-1, 32, 3)
            v = response.query(p)
            primary = -0.5 * np.real(v[:, :, 0, None] * g.conj())
            secondary, scalar, attraction, partner, eff = pairs(x, v, p, physics, chi)
            assert np.all((eff == 0) | (eff >= response.knots[0])) and np.all(eff <= response.knots[-1])
            effective = [min(effective[0], float(eff.min())), max(effective[1], float(eff.max()))]
            if first == 0:
                ref = pair_reference(x[0], v[0], p[0], physics, chi)
                assert np.max(abs(ref - secondary[0])) < max(np.max(abs(ref)) * 1e-12, 1e-24)
            forces.append(primary)
            sec.append(secondary)
            pressures.append(p)
            max_pair.append(attraction)
            partners.append(partner)
            pairmax.append(np.maximum(scalar, 0).max(axis=1))
            pairid.extend(np.argmax(scalar, axis=1).tolist())
        F = np.concatenate(forces)
        S = np.concatenate(sec)
        P = np.concatenate(pressures)
        attraction = np.concatenate(max_pair)
        partner = np.concatenate(partners)
        pairmax = np.concatenate(pairmax)
        valid = sample["time_s"] < times[-1] - 1e-13
        ids = np.rint(sample["time_s"][valid] / read(folder / "INPUT.json")["dt_s"]).astype(int)
        checks = {}
        for key, new in [
            ("pressure_Pa", P),
            ("primary_force_N", F),
            ("secondary_force_N", S),
            ("total_force_N", F + S),
        ]:
            diff = float(abs(new[ids] - sample[key][valid]).max())
            tol = max(float(abs(sample[key][valid]).max()) * 1e-9, 1e-24)
            assert diff <= tol, (tag, key, diff, tol)
            checks[key] = dict(max_absolute_difference=diff, tolerance=tol)
        # Validate pressure extrema at every accepted step, beyond the sparse force output.
        assert (
            np.max(abs(abs(P).min(1) - extra[:n, 3])) < 1e-8
            and np.max(abs(abs(P).max(1) - extra[:n, 4])) < 1e-8
        )
        meanvec = F.mean((0, 1))
        magnitudes = np.linalg.norm(F, axis=-1)
        secondarynorm = np.linalg.norm(S, axis=-1)
        center = xs[[0, -1]].mean(1)
        centerpb, _ = basis.query_basis(center)
        centerp = abs(u @ centerpb)
        peak = int(np.argmax(pairmax))
        pair_index = pairid[peak]
        single = np.unravel_index(np.argmax(secondarynorm), secondarynorm.shape)
        record = dict(
            source_path=str(folder.relative_to(ROOT)),
            control=e,
            phase_deg_signed=np.degrees(e["phase_rad"]),
            phase_deg_0_360=np.degrees(e["phase_rad"]) % 360,
            completed_s=times[-1],
            force_evaluation_time_s=[times[0], times[-2]],
            evaluated_step_count=n,
            pressure_mean_peak_amplitude_Pa=float(abs(P).mean()),
            pressure_min_max_Pa=[float(abs(P).min()), float(abs(P).max())],
            centroid_pressure_initial_final_Pa=centerp,
            secondary_pair_effective_pressure_min_max_Pa=effective,
            mean_primary_vector_per_bubble_N=meanvec,
            mean_primary_magnitude_per_bubble_N=float(magnitudes.mean()),
            mean_primary_magnitude_range_across_IDs_N=[
                float(magnitudes.mean(0).min()),
                float(magnitudes.mean(0).max()),
            ],
            max_primary_per_bubble_N=float(magnitudes.max()),
            group_mean_primary_vector_N=meanvec * 32,
            mean_push_unit_vector=meanvec / np.linalg.norm(meanvec)
            if np.linalg.norm(meanvec) > 0
            else np.zeros(3),
            max_pair_attraction_N=float(pairmax[peak]),
            max_pair_ids=[int(i[pair_index]), int(j[pair_index])] if pairmax[peak] > 0 else None,
            max_pair_time_s=times[peak] if pairmax[peak] > 0 else None,
            max_secondary_net_per_bubble_N=float(secondarynorm[single]),
            max_secondary_net_id=single[1],
            max_secondary_net_time_s=times[single[0]],
            max_primary_plus_secondary_per_bubble_N=float(np.linalg.norm(F + S, axis=-1).max()),
            secondary_internal_force_cancellation_N=float(np.linalg.norm(S.sum(1), axis=1).max()),
            original_saved_GPU_checks=checks,
            wall_s=time.perf_counter() - started,
            trajectory_status=result["status"],
            scope="Arithmetic average over all accepted dt=1us force evaluations and 32 fixed IDs. Forces are cycle-averaged; pressures are peak harmonic amplitudes. No rejected proposal forces, no new trajectory. Y shape/ID convergence failed; its extrema are discretization-specific.",
        )
        reports[tag] = record
        for bubble in range(32):
            peak_step = int(np.argmax(attraction[:, bubble]))
            r = dict(
                run=tag,
                particle_id=bubble,
                primary_mean_x_pN=F[:, bubble, 0].mean() * 1e12,
                primary_mean_y_pN=F[:, bubble, 1].mean() * 1e12,
                primary_mean_z_pN=F[:, bubble, 2].mean() * 1e12,
                primary_mean_magnitude_pN=magnitudes[:, bubble].mean() * 1e12,
                primary_max_magnitude_pN=magnitudes[:, bubble].max() * 1e12,
                max_pair_attraction_pN=attraction[peak_step, bubble] * 1e12,
                max_pair_partner_id=int(partner[peak_step, bubble]),
                max_pair_time_ms=times[peak_step] * 1e3,
                secondary_net_max_pN=secondarynorm[:, bubble].max() * 1e12,
                pressure_mean_Pa=abs(P[:, bubble]).mean(),
                pressure_min_Pa=abs(P[:, bubble]).min(),
                pressure_max_Pa=abs(P[:, bubble]).max(),
            )
            all_rows.append(r)
        np.savez_compressed(
            OUT / (tag + "_force_evaluation_data.npz"),
            time_s=times[:-1],
            primary_force_N=F,
            secondary_force_N=S,
            pressure_Pa=P,
            max_attractive_pair_per_bubble_N=attraction,
            partner_id=partner,
        )
        for f in ["trajectory.h5", "INPUT.json", "RESULT.json"]:
            hashes[str((folder / f).relative_to(ROOT))] = sha(folder / f)
        print(
            tag,
            {
                k: record[k]
                for k in [
                    "phase_deg_0_360",
                    "pressure_mean_peak_amplitude_Pa",
                    "pressure_min_max_Pa",
                    "centroid_pressure_initial_final_Pa",
                    "mean_primary_vector_per_bubble_N",
                    "mean_primary_magnitude_per_bubble_N",
                    "max_pair_attraction_N",
                    "max_pair_ids",
                    "max_pair_time_s",
                    "max_secondary_net_per_bubble_N",
                    "wall_s",
                ]
            },
            flush=True,
        )
    with (OUT / "per_bubble_forces.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0]))
        w.writeheader()
        w.writerows(all_rows)
    with (OUT / "pzt_controls.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "run",
                "PZT",
                "center_x_mm",
                "center_y_mm",
                "center_z_mm",
                "normal_x",
                "normal_y",
                "normal_z",
                "diameter_mm",
                "enabled",
                "source_amplitude_Pa",
                "phase_deg_0_360",
                "phase_deg_signed",
            ]
        )
        for tag, v in reports.items():
            for k, element in enumerate(geometry["pzt"]):
                w.writerow(
                    [
                        tag,
                        element["id"],
                        *(np.array(element["center_local_m"]) * 1e3),
                        *element["normal_local"],
                        element["active_diameter_m"] * 1e3,
                        v["control"]["enabled"][k],
                        v["control"]["amplitude_Pa"][k],
                        v["phase_deg_0_360"][k],
                        v["phase_deg_signed"][k],
                    ]
                )
    hashes.update(
        {
            str((OLD / f).relative_to(ROOT)): sha(OLD / f)
            for f in ["CONFIG_SNAPSHOT.json", "cache/coarse.h5", "inputs/response.h5"]
        }
    )
    write(
        OUT / "SUMMARY.json",
        dict(
            source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            script_sha256=sha(__file__),
            input_sha256=hashes,
            pzt_geometry=geometry["pzt"],
            source_amplitude_definition="rho*c*normal surface velocity in Pa; not calibrated transducer-face pressure and not voltage",
            harmonic_convention="Re[P exp(+i omega t)]",
            frequency_Hz=physics.frequency_Hz,
            radius_m=physics.radius_m,
            initial_COM_m=[0, 0, 0],
            runs=reports,
            postprocess_wall_s=time.perf_counter() - start,
            training_ready=False,
            labels=["EXPLORATORY", "HOMOGENEOUS_WATER_CONTROL_PROXY", "PHYSICS_UNCALIBRATED"],
        ),
    )
    write(
        OUT / "MANIFEST.json",
        {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != "MANIFEST.json"},
    )


if __name__ == "__main__":
    main()
