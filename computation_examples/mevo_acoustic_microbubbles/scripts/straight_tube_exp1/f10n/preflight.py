"""Preregister before rollouts; validate piston integral, gradients and caches on CPU."""

import time
import json
import numpy as np
import h5py
from scipy.spatial.distance import pdist
from mevo_mb.acoustics.water_disc_proxy import WaterDiscs, CachedWaterBasis, cubic_query
from mevo_mb.coupling.pdms_cavity import apparatus
from mevo_mb.fast_sim.models import ResponseSection, random_cloud
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.metrics import frame_metrics
from .common import BASE, configs, write, sha, identity


def main():
    started = time.perf_counter()
    c = configs()
    p = c["protocol"]
    r = p["radius_m"]
    (BASE / "analysis").mkdir(exist_ok=True)
    (BASE / "cache").mkdir(exist_ok=True)
    assert not (BASE / "PREREGISTRATION.json").exists(), "Preserve completed preregistration"
    cavity, _ = apparatus(c["geometry"])
    proxy = WaterDiscs(c["geometry"]["pzt"], p["frequency_Hz"])
    grid = np.stack(
        np.meshgrid(np.arange(4) - 1.5, np.arange(4) - 1.5, np.arange(2) - 0.5, indexing="ij"), -1
    ).reshape(-1, 3)
    aggregate = grid * (2 * r + p["initial_surface_gap_over_radius"] * r)
    dispersed = random_cloud(32, r, p["dispersed_cloud_radius_m"], p["seed"])
    dispersed -= dispersed.mean(0)
    np.savez_compressed(
        BASE / "initial_states.npz",
        aggregate_m=aggregate,
        dispersed_m=dispersed,
        particle_id=np.arange(32),
        radius_m=r,
    )
    checks = {}
    for key, x in [("aggregate", aggregate), ("dispersed", dispersed)]:
        m = frame_metrics(x, r)
        checks[key + "_finite_radius_valid"] = all(cavity.query(q, r)["valid"] for q in x)
        checks[key + "_nonoverlap"] = pdist(x).min() >= 2 * r
        checks[key + "_network_correct"] = m["cluster_count"] == (1 if key == "aggregate" else 32)
        write(BASE / "analysis" / ("INITIAL_" + key + ".json"), m)
    response = ResponseSection(
        BASE / "inputs/response.h5", F2Physics(radius_m=r, frequency_Hz=p["frequency_Hz"])
    )
    assert sha(BASE / "inputs/response.h5") == p["response_sha256"]
    write(BASE / "analysis/RESPONSE_SCOPE.json", response.metadata)
    rng = np.random.default_rng(p["seed"])
    q = rng.uniform(-45e-6, 45e-6, (48, 3))
    a, g = proxy.evaluate(q, 48)
    a2, g2 = proxy.evaluate(q, 96)
    ar, gr = proxy.rayleigh_reference(q[:12])
    err = dict(
        quadrature_pressure=float(np.max(abs(a - a2))),
        quadrature_gradient_relative=float(np.max(abs(g - g2)) / np.max(abs(g2))),
        area_reference_pressure=float(np.max(abs(a[:, :12] - ar))),
        area_reference_gradient_relative=float(np.max(abs(g[:, :12] - gr)) / np.max(abs(gr))),
    )
    checks["FNM_matches_independent_Rayleigh_area"] = (
        err["area_reference_pressure"] < 1e-9 and err["area_reference_gradient_relative"] < 1e-9
    )
    checks["48_to_96_quadrature_convergence"] = (
        err["quadrature_pressure"] < 1e-10 and err["quadrature_gradient_relative"] < 1e-10
    )
    fd = []
    for j in range(3):
        dx = np.eye(3)[j] * 1e-8
        fd.append((proxy.evaluate(q + dx)[0] - proxy.evaluate(q - dx)[0]) / 2e-8)
    err["analytic_gradient_vs_finite_difference_relative"] = float(
        np.max(abs(g - np.stack(fd, -1))) / np.max(abs(g))
    )
    checks["analytic_gradient_consistent"] = err["analytic_gradient_vs_finite_difference_relative"] < 1e-6
    axial = []
    for i, e in enumerate(c["geometry"]["pzt"]):
        n = np.array(e["normal_local"])
        z = -np.dot(e["center_local_m"], n)
        rad = e["active_diameter_m"] / 2
        exact = np.exp(-1j * proxy.k * z) - np.exp(-1j * proxy.k * np.sqrt(z * z + rad * rad))
        axial.append(abs(proxy.evaluate([[0, 0, 0]])[0][i, 0] - exact))
    checks["on_axis_closed_form"] = max(axial) < 1e-11
    for name, spacing in [("coarse", p["cache_spacing_m"]), ("fine", p["cache_refinement_spacing_m"])]:
        t = time.perf_counter()
        n = round(2 * p["cache_half_width_m"] / spacing) + 1
        axes = tuple(np.linspace(-p["cache_half_width_m"], p["cache_half_width_m"], n) for _ in range(3))
        points = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
        pp = np.empty((4, len(points)), complex)
        for i in range(0, len(points), 4096):
            pp[:, i : i + 4096] = proxy.evaluate(points[i : i + 4096], p["quadrature_order"])[0]
        metadata = dict(
            model=p["model"],
            frequency_Hz=p["frequency_Hz"],
            source_amplitude_meaning=p["source_amplitude_meaning"],
            elements=c["geometry"]["pzt"],
            harmonic_convention="Re[P exp(+i omega t)]",
            cache_spacing_m=spacing,
            quadrature_order=p["quadrature_order"],
            literature="https://pmc.ncbi.nlm.nih.gov/articles/PMC2117898/ Eq.9",
            training_ready=False,
            identity=identity(),
        )
        path = BASE / "cache" / f"{name}.h5"
        with h5py.File(path, "x") as h:
            h.attrs["metadata_json"] = json.dumps(metadata)
            for s, axis in zip("xyz", axes):
                h.create_dataset("axes/" + s + "_m", data=axis)
            h.create_dataset(
                "pressure_basis", data=pp.reshape((4, n, n, n)), compression="gzip", compression_opts=1
            )
        basis = CachedWaterBasis(path)
        pc, gc = basis.query_basis(q)
        ep = float(np.max(abs(pc - a)))
        eg = float(np.max(abs(gc - g)) / np.max(abs(g)))
        checks[name + "_cache_pressure_gradient_accuracy"] = ep < 1e-4 and eg < 0.002
        err[name] = dict(
            pressure_error_per_source_Pa=ep,
            gradient_relative=eg,
            build_wall_s=time.perf_counter() - t,
            bytes=path.stat().st_size,
            sha256=sha(path),
        )
        u = np.array([300, 250j, -100, 75 - 120j])
        field = np.einsum("i,ixyz->xyz", u, basis.pressures)
        v, d = cubic_query(field, basis.origin, basis.spacing, q)
        checks[name + "_linear_controls"] = np.allclose(v, u @ pc, rtol=1e-12, atol=1e-10)
        fdc = []
        for j in range(3):
            dx = np.eye(3)[j] * 1e-9
            fdc.append(
                (
                    cubic_query(field, basis.origin, basis.spacing, q + dx)[0]
                    - cubic_query(field, basis.origin, basis.spacing, q - dx)[0]
                )
                / 2e-9
            )
        checks[name + "_interpolant_gradient_consistent"] = (
            np.max(abs(np.stack(fdc, -1) - d)) / np.max(abs(d)) < 1e-6
        )
        plane = q.copy()
        plane[:, 2] = 0
        _, zgrad = basis.query_basis(plane)
        checks[name + "_Z_symmetry"] = np.max(abs(zgrad[:, :, 2])) < 1e-7
        try:
            basis.vector([1] * 3, [0] * 3)
            checks["reject_three_sources"] = False
        except ValueError:
            checks["reject_three_sources"] = True
    checks["target_legal"] = cavity.query(np.array(p["target_position_m"]), r)["valid"]
    checks["response_no_extrapolation"] = response.metadata["extrapolation"] is False
    write(BASE / "CONFIG_SNAPSHOT.json", c)
    write(
        BASE / "analysis/CPU_CHECKS.json",
        dict(
            status="PASS" if all(checks.values()) else "FAIL",
            checks=checks,
            errors=err,
            wall_s=time.perf_counter() - started,
            identity=identity(),
        ),
    )
    assert all(checks.values()), checks
    write(
        BASE / "PREREGISTRATION.json",
        dict(
            protocol=p,
            initial_sha256=sha(BASE / "initial_states.npz"),
            cache_sha256={n: sha(BASE / "cache" / f"{n}.h5") for n in ["coarse", "fine"]},
            target_locked_before_any_GPU_trajectory=True,
            position_objective="endpoint fixed-ID COM target error + registered contact/wall penalties; no target relocation",
            shape_objective="abs(sqrt(Qxx/Qyy)-1.25), separate from centroid displacement; preserve all IDs",
            optimizer="two bounded static force predictor starts, then at most eight trajectory evaluations and bounded deterministic local perturbations; shape at most four",
            stop="budget, pressure invalid, contact/wall violation, cache exit; failed prefixes retained; no silent clipping",
            hypotheses="Z reflection symmetry forbids independent signed Z translation at symmetric center; measured response decides",
            identity=identity(),
        ),
    )
    print("CPU preflight PASS", len(checks), "checks", time.perf_counter() - started, flush=True)


if __name__ == "__main__":
    main()
