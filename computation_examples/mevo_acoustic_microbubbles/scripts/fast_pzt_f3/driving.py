"""First F3 experiment: six instantaneous controls at the frozen F2 mature state."""

import argparse
import numpy as np
import h5py
from mevo_mb.fast_sim.models import Basis, ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics, pair_reference
from mevo_mb.fast_sim.f2_engine import F2Simulator
from mevo_mb.fast_sim.kinematics import affine_velocity, force_summary, rms_vectors
from .common import BASE, F2, PROTOCOL, write, js, state, identity
from .snapshot import evaluate


def main(label):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    p = F2Physics()
    x = state()[0]
    basis = Basis(BASE / "inputs/basis.h5", p.frequency_Hz)
    table = ResponseSection(BASE / "inputs/response.h5", p)
    f2 = js(F2 / "mature_v1/FIELD_SELECTION.json")
    base = f2["baseline"]
    phase = js(F2 / "phase_v1/FIELD_SELECTION.json")["event"]
    margin = PROTOCOL["screening"]["roi_margin_radii"] * p.radius_m
    axes = [
        np.linspace(x[:, a].min() - margin, x[:, a].max() + margin, PROTOCOL["screening"]["roi_axis_samples"])
        for a in range(3)
    ]
    roi = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)

    def field(event):
        return basis.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])[0]

    f0 = field(base)
    target = np.mean(abs(f0.interpolate_complex_pressure(roi)) ** 2)
    # One analytic Z-pair phase pattern, selected without trajectory information.
    a = field(
        dict(base, phase_rad=[0, 0, np.pi / 6, 0], enabled=[True, True, True, False])
    ).interpolate_complex_pressure(roi)
    b = field(dict(base, enabled=[False, False, False, True])).interpolate_complex_pressure(roi)
    cross = np.mean(a * np.conj(b))
    offset = np.mean(abs(a) ** 2 + abs(b) ** 2)
    cosine = float((target - offset) / (2 * abs(cross)))
    assert -1 <= cosine <= 1, "Analytic phase candidate unavailable; retain failure, no search"
    roots = [float(np.angle(cross) + s * np.arccos(cosine)) for s in [-1, 1]]
    phase_z = dict(base, phase_rad=[0, 0, float(np.pi / 6), min(roots, key=abs)])
    candidates = [
        ("unchanged", "baseline", base),
        ("phase_F2", "phase", phase),
        ("phase_Z", "phase", phase_z),
    ]
    candidates += [
        (v["name"], "amplitude", v["event"]) for v in f2["candidates"] if v["family"] == "amplitude"
    ]
    assert len(candidates) == 6
    cells = np.floor((x - basis.origin) / basis.spacing).astype(int)
    unique, counts = np.unique(cells, axis=0, return_counts=True)
    grid = dict(
        spacing_m=basis.spacing,
        cluster_center_span_m=np.ptp(x, axis=0),
        sphere_extent_m=np.ptp(x, axis=0) + 2 * p.radius_m,
        particle_cell_indices=cells,
        occupied_cells=unique,
        occupancy=counts,
        origin_m=basis.origin,
        grid_shape=basis.shape,
        gradient_source="Stored nodal numpy.gradient on k-Wave grid, separately trilinearly interpolated from pressure",
    )
    force_scale = 6 * np.pi * p.viscosity_Pa_s * p.radius_m * PROTOCOL["force_scale"]["reference_speed_m_s"]
    sim = F2Simulator(x, p, table, basis)
    rows = {}
    raw = {}
    initial = None
    for name, family, event in candidates:
        snap = evaluate(sim, event)
        f = sim.last_field
        pressure = f.interpolate_complex_pressure(x)
        g = f.interpolate_complex_gradient(x)
        rp = f.interpolate_complex_pressure(roi)
        rg = f.interpolate_complex_gradient(roi)
        np.testing.assert_allclose(pressure, snap["pressure_complex_Pa"], rtol=1e-11, atol=1e-9)
        hp = table.query(pressure)
        primary = -0.5 * np.real(hp[:, 0, None] * np.conj(g))
        secondary = pair_reference(x, hp, pressure, p, sim.chi_zero)
        np.testing.assert_allclose(primary, snap["primary_N"], rtol=1e-10, atol=1e-24)
        np.testing.assert_allclose(secondary, snap["secondary_N"], rtol=1e-10, atol=1e-24)
        table.query(rp)  # no silent out-of-domain ROI pressure
        frac = (x - basis.origin) / basis.spacing
        distance = np.minimum(frac - np.floor(frac), np.ceil(frac) - frac) * basis.spacing
        eps = np.minimum(5e-8, distance * 0.25)
        assert (eps > 1e-14).all(), (
            "Finite difference stencil touches a cell plane; needs explicit one-sided audit"
        )
        fd = []
        for scale in [1, 0.5]:
            gg = np.zeros_like(g)
            for a in range(3):
                delta = np.zeros_like(x)
                delta[:, a] = eps[:, a] * scale
                gg[:, a] = (
                    f.interpolate_complex_pressure(x + delta) - f.interpolate_complex_pressure(x - delta)
                ) / (2 * eps[:, a] * scale)
            fd.append(gg)
        primary_fd = -0.5 * np.real(hp[:, 0, None] * np.conj(fd[0]))
        local = dict(
            particle_pressure_RMS_Pa=float(np.sqrt(np.mean(abs(pressure) ** 2))),
            roi_pressure_RMS_Pa=float(np.sqrt(np.mean(abs(rp) ** 2))),
            particle_pressure_min_Pa=float(abs(pressure).min()),
            particle_pressure_max_Pa=float(abs(pressure).max()),
            roi_pressure_min_Pa=float(abs(rp).min()),
            roi_pressure_max_Pa=float(abs(rp).max()),
            gradient_relative_L2_vs_pressure_derivative=float(
                np.linalg.norm(g - fd[0]) / np.linalg.norm(fd[0])
            ),
            derivative_step_halving_relative_L2=float(np.linalg.norm(fd[0] - fd[1]) / np.linalg.norm(fd[0])),
            primary_relative_L2_vs_pressure_derivative=float(
                np.linalg.norm(primary - primary_fd) / max(np.linalg.norm(primary_fd), 1e-300)
            ),
            gradient_error_action="Retain frozen production gradient; inconsistency is a field-discretization limitation, not a force correction",
        )
        row = dict(
            family=family,
            event=event,
            force=force_summary(primary, secondary, force_scale),
            field=local,
            kinematics={
                key: affine_velocity(x, snap[key])
                for key in ["free_velocity_m_s", "local_velocity_m_s", "mobility_velocity_m_s"]
            },
            primary_only_Stokes_kinematics=affine_velocity(
                x, primary / (6 * np.pi * p.viscosity_Pa_s * p.radius_m)
            ),
            coupled_pressure_quality=sim.quality(),
        )
        if initial is None:
            initial = (snap, row)
        delta = {
            key: snap[key] - initial[0][key]
            for key in ["primary_N", "secondary_N", "total_N", "free_velocity_m_s", "mobility_velocity_m_s"]
        }
        for key in ["primary_N", "secondary_N", "total_N"]:
            delta["centered_" + key] = delta[key] - delta[key].mean(0)
        row["control_difference"] = dict(
            centered_primary_RMS_N=rms_vectors(delta["centered_primary_N"]),
            centered_secondary_RMS_N=rms_vectors(delta["centered_secondary_N"]),
            centered_total_RMS_N=rms_vectors(delta["centered_total_N"]),
            mean_primary_delta_N=delta["primary_N"].mean(0),
            mean_total_delta_N=delta["total_N"].mean(0),
            free_kinematics=affine_velocity(x, delta["free_velocity_m_s"]),
            mobility_kinematics=affine_velocity(x, delta["mobility_velocity_m_s"]),
        )
        row["pressure_particle_relative_change"] = abs(
            local["particle_pressure_RMS_Pa"] / initial[1]["field"]["particle_pressure_RMS_Pa"] - 1
        )
        row["pressure_ROI_relative_change"] = abs(
            local["roi_pressure_RMS_Pa"] / initial[1]["field"]["roi_pressure_RMS_Pa"] - 1
        )
        row["eligible"] = bool(
            max(event["amplitude_Pa"]) <= PROTOCOL["screening"]["source_amplitude_cap_Pa"] + 1e-9
            and max(row["pressure_particle_relative_change"], row["pressure_ROI_relative_change"])
            <= PROTOCOL["screening"]["pressure_RMS_relative_tolerance"]
            and row["control_difference"]["mobility_kinematics"]["status"] == "IDENTIFIABLE"
        )
        row["score_per_s"] = row["control_difference"]["mobility_kinematics"].get(
            "deviatoric_strain_norm_per_s", 0.0
        )
        ds = max(initial[1]["force"]["D_S_N"], force_scale)
        row["delta_primary_over_baseline_secondary"] = (
            row["control_difference"]["centered_primary_RMS_N"] / ds
        )
        row["delta_total_over_baseline_secondary"] = row["control_difference"]["centered_total_RMS_N"] / ds
        row["strong_driving_diagnostic"] = bool(
            row["delta_primary_over_baseline_secondary"] >= 0.01
            or (row["delta_total_over_baseline_secondary"] >= 0.01 and row["score_per_s"] * 0.12 >= 0.05)
        )
        rows[name] = row
        raw[name] = dict(
            **snap,
            **{"control_delta_" + k: v for k, v in delta.items()},
            stored_gradient_Pa_m=g,
            pressure_derivative_Pa_m=fd[0],
            fd_epsilon_m=eps,
            primary_from_pressure_derivative_N=primary_fd,
            roi_pressure_Pa=rp,
            roi_gradient_Pa_m=rg,
        )
        print(
            name, "D", row["force"]["D"], "score", row["score_per_s"], "eligible", row["eligible"], flush=True
        )
    selected = {
        family: max(
            [name for name in rows if rows[name]["family"] == family and rows[name]["eligible"]],
            key=lambda name: rows[name]["score_per_s"],
        )
        for family in ["phase", "amplitude"]
    }
    limited = all(
        r["force"]["D"] < 0.01 and r["force"]["primary_nonuniformity"] < 0.1
        for r in rows.values()
        if r["eligible"]
    )
    result = dict(
        complete=True,
        identity=identity(),
        inputs=dict(basis=basis.metadata, response=table.metadata),
        grid=grid,
        force_scale_N=force_scale,
        candidate_rows=rows,
        selected=selected,
        phase_Z_analytic_roots_rad=roots,
        differential_driving_status="FIELD_DIFFERENTIAL_DRIVE_LIMITED"
        if limited
        else "DIFFERENTIAL_DRIVE_MEASURED",
        screening_status="STRONG_DRIVING_DIAGNOSTIC_FOUND"
        if any(rows[n]["strong_driving_diagnostic"] for n in selected.values())
        else "NO_STRONG_SHAPE_DRIVING_CONTROL_FOUND",
        meaning="Instantaneous diagnostics at identical actual mature positions; neither force ratio nor strain predicts validated sustained deformation",
    )
    with h5py.File(out / "FORCE_FIELD_DATA.h5", "w") as h:
        h["positions_m"] = x
        h["particle_id"] = np.arange(len(x))
        h["radii_m"] = np.full(len(x), p.radius_m)
        h["roi_points_m"] = roi
        for name, arrays in raw.items():
            for key, value in arrays.items():
                h.create_dataset(name + "/" + key, data=value)
    write(out / "RESULT.json", result)
    write(
        out / "SELECTED.json",
        dict(
            selected={f: rows[n]["event"] for f, n in selected.items()},
            names=selected,
            selection_before_trajectories=True,
            score=PROTOCOL["screening"]["score"],
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    main(a.label)
