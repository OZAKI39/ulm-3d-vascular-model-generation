"""Bounded F6 admission: interpolation, RPY and exact archived F5 prefixes."""

import argparse
import h5py
import numpy as np
import warp as wp
from mevo_mb.fast_sim.consistent_field import query_numpy, consistent_primary
from mevo_mb.fast_sim.rpy import matrix_numpy, velocity
from mevo_mb.fast_sim.f6_engine import F6Simulator
from .common import BASE, F5, sha, write, js, identity, objects


def gpu_query(field, origin, spacing, x, harmonics=None):
    n = len(x)
    shape = field.shape
    g = np.zeros((field.size, 4, 2))
    g[:, 0, 0] = field.real.ravel()
    g[:, 0, 1] = field.imag.ravel()
    pp = wp.zeros(n, dtype=wp.vec2d)
    grad = wp.zeros((n, 3), dtype=wp.vec2d)
    f = wp.zeros(n, dtype=wp.vec3d)
    v = np.zeros((n, 8, 2))
    if harmonics is not None:
        v[..., 0] = harmonics.real
        v[..., 1] = harmonics.imag
    error = wp.zeros(1, dtype=int)
    wp.launch(
        consistent_primary,
        n,
        [
            wp.array(x, dtype=wp.vec3d),
            wp.array(g, dtype=wp.vec2d),
            wp.vec3d(*origin),
            wp.vec3d(*spacing),
            wp.vec3d(*(origin + spacing * (np.array(shape) - 1))),
            *shape,
            wp.ones(n, dtype=int),
            wp.array(v, dtype=wp.vec2d),
            f,
            pp,
            grad,
            error,
        ],
    )
    p = pp.numpy()
    gg = grad.numpy()
    return p[:, 0] + 1j * p[:, 1], gg[..., 0] + 1j * gg[..., 1], f.numpy(), int(error.numpy()[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    wp.init()
    wp.set_device("cuda:0")
    rows = {}
    origin = np.array([-2.0, -1.0, -0.5])
    spacing = np.array([0.7, 0.4, 0.3])
    shape = (5, 4, 3)
    axes = [origin[d] + np.arange(shape[d]) * spacing[d] for d in range(3)]
    xx, yy, zz = np.meshgrid(*axes, indexing="ij")
    # Exact multilinear polynomial includes mixed xyz derivatives.
    field = (
        (1 + 2j)
        + (2 - 3j) * xx
        + (4 + 1j) * yy
        + (-3 + 0.2j) * zz
        + (1 + 2j) * xx * yy
        + (3 - 1j) * xx * yy * zz
    )
    q = np.random.default_rng(6).uniform(0.12, np.array(shape) - 1.12, (23, 3))
    x = origin + q * spacing
    p, g = query_numpy(field, origin, spacing, x)
    exact = np.column_stack(
        [
            (2 - 3j) + (1 + 2j) * x[:, 1] + (3 - 1j) * x[:, 1] * x[:, 2],
            (4 + 1j) + (1 + 2j) * x[:, 0] + (3 - 1j) * x[:, 0] * x[:, 2],
            (-3 + 0.2j) + (3 - 1j) * x[:, 0] * x[:, 1],
        ]
    )
    np.testing.assert_allclose(g, exact, rtol=2e-13, atol=2e-13)
    gp, gg, _, error = gpu_query(field, origin, spacing, x)
    assert error == 0
    np.testing.assert_allclose(gp, p, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(gg, g, rtol=2e-13, atol=2e-13)
    # Finite differences never cross cells.
    for eps in (1e-5, 5e-6):
        fd = np.column_stack(
            [
                (
                    query_numpy(field, origin, spacing, x + np.eye(3)[d] * eps)[0]
                    - query_numpy(field, origin, spacing, x - np.eye(3)[d] * eps)[0]
                )
                / (2 * eps)
                for d in range(3)
            ]
        )
        np.testing.assert_allclose(fd, g, rtol=2e-9, atol=2e-9)
    rows["analytic_numpy_warp_and_cell_interior_FD"] = "PASS"
    # Boundary convention explicitly tested on piecewise x^2 nodal samples.
    field = xx**2 + 0j
    xb = np.array([origin + spacing * np.array([2.0, 0.4, 0.3]), origin + spacing * (np.array(shape) - 1)])
    _, gb = query_numpy(field, origin, spacing, xb)
    expected = [axes[0][2] + axes[0][3], axes[0][-2] + axes[0][-1]]
    np.testing.assert_allclose(gb[:, 0], expected, atol=1e-13)
    _, gg, _, err = gpu_query(field, origin, spacing, xb)
    assert err == 0
    np.testing.assert_allclose(gg, gb, atol=1e-13)
    for bad in (
        origin - np.array([1e-10, 0, 0]),
        origin + spacing * (np.array(shape) - 1) + 1e-10,
        np.array([np.nan, 0, 0]),
    ):
        try:
            query_numpy(field, origin, spacing, bad)
        except ValueError:
            pass
        else:
            raise AssertionError("CPU accepted invalid query")
        assert gpu_query(field, origin, spacing, np.atleast_2d(bad))[-1] & 1
    rows["one_sided_faces_maximum_and_invalid_rejection"] = "PASS"
    p, basis, response, geometry, flow = objects()
    for rel, entry in js(F5 / "INPUT_MANIFEST.json").items():
        assert sha(F5 / rel) == entry["sha256"], rel
    rows["current_input_identity"] = "PASS"
    with h5py.File(F5 / "science_v1/main/trajectory.h5") as h:
        configs = js(F5 / "ACOUSTIC_RELEASE_PROTOCOL.json")
        for name, start in [("ON", 0), ("OFF", 600)]:
            e = dict(configs["events"][int(start > 0)], time_s=0.0)
            sim = F6Simulator(
                h["positions_m"][start, 0],
                p,
                response,
                basis,
                vessel=geometry,
                flow=flow,
                gradient_mode="legacy_nodal",
                acoustic_enabled=start == 0,
            )
            run = sim.run(0.001, 1e-5, [e], 0.0001)
            assert run["summary"]["status"] == "COMPLETE"
            np.testing.assert_array_equal(run["positions_m"][::5], h["positions_m"][start : start + 3])
            rows["F5_" + name + "_archived_prefix"] = "BITWISE_PASS"
        positions = [
            h["positions_m"][0, 0],
            h["positions_m"][600, 0],
            np.array([[1.1e-6, 80e-6, 40e-6], [2.2e-6, 74e-6, 35e-6]]),
        ]
    force_rows = []
    for region, x in zip(["F5_initial", "F5_switch", "F3_region"], positions):
        e = (
            configs["events"][0]
            if region != "F3_region"
            else dict(time_s=0.0, amplitude_Pa=[151.340514] * 4, phase_rad=[0.0] * 4)
        )
        e = dict(e, flow_enabled=False)
        e.pop("mode", None)
        field, _, _ = basis.combine(e["amplitude_Pa"], e["phase_rad"])
        pressure, grad = query_numpy(field.pressure_Pa, basis.origin, basis.spacing, x)
        vv = response.query(pressure)
        sim = F6Simulator(x, p, response, basis, vessel_enabled=False, gradient_mode="consistent_trilinear")
        sim.control(e)
        sim.fields()
        wanted = -0.5 * np.real(vv[:, 0, None] * grad.conj())
        np.testing.assert_allclose(sim.force.numpy(), wanted, rtol=3e-11, atol=1e-26)
        old = F6Simulator(x, p, response, basis, vessel_enabled=False, gradient_mode="legacy_nodal")
        old.control(e)
        old.fields()
        of = old.force.numpy()
        force_rows.append(
            dict(
                region=region,
                pressure_range_Pa=[abs(pressure).min(), abs(pressure).max()],
                legacy_mean_N=of.mean(0),
                consistent_mean_N=wanted.mean(0),
                relative_L2_force_change=np.linalg.norm(wanted - of) / np.linalg.norm(of),
                mean_direction_cosine=np.dot(of.mean(0), wanted.mean(0))
                / (np.linalg.norm(of.mean(0)) * np.linalg.norm(wanted.mean(0))),
            )
        )
    rows["real_pressure_primary_numpy_warp"] = "PASS"
    write(out / "GRADIENT_FORCE_COMPARISON.json", force_rows)
    rr = p.radius_m
    mu = p.viscosity_Pa_s
    # N=1, separated pair, regularized overlap pair, irregular small group.
    states = [
        np.zeros((1, 3)),
        np.array([[0, 0, 0], [4 * rr, 0, 0]]),
        np.array([[0, 0, 0], [1.9 * rr, 0, 0]]),
        rr * np.array([[0, 0, 0], [2.05, 0, 0], [0.4, 2.1, 0.1], [0.4, 0.5, 2.2]]),
    ]
    eigs = []
    for x in states:
        n = len(x)
        M = matrix_numpy(x, rr, mu)
        np.testing.assert_allclose(M, M.T, rtol=0, atol=0)
        eigs.append(float(np.linalg.eigvalsh(M).min()))
        assert eigs[-1] > 0
        forces = np.random.default_rng(9).normal(size=x.shape) * 1e-13
        u = wp.zeros(n, dtype=wp.vec3d)
        err = wp.zeros(1, dtype=int)
        wp.launch(
            velocity,
            n,
            [
                wp.array(x, dtype=wp.vec3d),
                wp.array(forces, dtype=wp.vec3d),
                wp.zeros(n, dtype=wp.vec3d),
                wp.ones(n, dtype=int),
                u,
                n,
                rr,
                mu,
                0,
                0.0,
                wp.zeros(n, dtype=wp.float64),
                err,
            ],
        )
        assert not err.numpy()[0]
        np.testing.assert_allclose(u.numpy().ravel(), M @ forces.ravel(), rtol=2e-13, atol=1e-20)
        if n == 1:
            np.testing.assert_allclose(u.numpy(), forces / (6 * np.pi * mu * rr), rtol=2e-15)
    # Equal and opposite force pair has zero centroid drift, including overlap.
    for d in [1.9, 2.0, 2.01, 10.0]:
        x = np.array([[0, 0, 0], [d * rr, 0, 0]])
        M = matrix_numpy(x, rr, mu)
        vel = (M @ np.array([1.0, 0, 0, -1.0, 0, 0])).reshape(2, 3)
        np.testing.assert_allclose(vel.sum(0), 0, atol=1e-8)
    rows["RPY_limit_reciprocity_symmetry_SPD_GPU"] = dict(status="PASS", minimum_eigenvalues=eigs)
    write(
        out / "CHECKS.json", dict(status="PASS", checks=rows, physical_validation=False, identity=identity())
    )
    print("F6 CHECKS PASS", flush=True)


if __name__ == "__main__":
    main()
