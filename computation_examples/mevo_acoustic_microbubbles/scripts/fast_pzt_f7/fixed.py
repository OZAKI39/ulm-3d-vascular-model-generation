"""Actual-state force/mobility audit; no integration or synthetic particle states."""

import argparse
import h5py
import numpy as np
from mevo_mb.fast_sim.f6_engine import F6Simulator
from mevo_mb.fast_sim import kernels as k
from mevo_mb.fast_sim.consistent_field import query_numpy
from mevo_mb.fast_sim.nearfield import pair_reference
from mevo_mb.fast_sim.rpy import matrix_numpy
from .common import ROOT, BASE, MODES, objects, write, js, identity, sha


def independent_rpy(x, a, eta):
    """Independent vectorized tensor, axes (i,j,component,component)."""
    dvec = x[:, None] - x[None, :]
    distance = np.linalg.norm(dvec, axis=-1)
    safe = np.where(distance == 0, 1.0, distance)
    nn = dvec[..., :, None] * dvec[..., None, :] / safe[..., None, None] ** 2
    self_m = 1 / (6 * np.pi * eta * a)
    iso = np.where(
        distance >= 2 * a,
        (1 + 2 * a * a / (3 * safe**2)) / (8 * np.pi * eta * safe),
        self_m * (1 - 9 * distance / (32 * a)),
    )
    axial = np.where(
        distance >= 2 * a,
        (1 - 2 * a * a / safe**2) / (8 * np.pi * eta * safe),
        self_m * 3 * distance / (32 * a),
    )
    tensor = iso[..., None, None] * np.eye(3) + axial[..., None, None] * nn
    return tensor.transpose(0, 2, 1, 3).reshape(3 * len(x), 3 * len(x))


def local_numpy(x, free, p):
    n = len(x)
    blocks = np.zeros((n, n, 3, 3))
    parent = np.arange(n)

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    def phi(e):
        return 0.5 / e - 0.45 * np.log(e) - 3 / 56 * e * np.log(e)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dx = x[j] - x[i]
            d = np.linalg.norm(dx)
            gap = (d - 2 * p.radius_m) / p.radius_m
            if gap < p.lubrication_activation:
                c = 0.5 * (phi(max(gap, p.lubrication_floor)) - phi(p.lubrication_activation))
                blocks[i, j] = c * np.outer(dx, dx) / d**2
                a, b = root(i), root(j)
                parent[max(a, b)] = min(a, b)
    old = free.copy()
    for _ in range(p.mobility_iterations):
        rhs = free + np.einsum("ijab,jb->ia", blocks, old)
        old = 0.5 * (old + np.linalg.solve(np.eye(3) + blocks.sum(axis=1), rhs[..., None])[..., 0])
    labels = np.array([root(i) for i in range(n)])
    for label in np.unique(labels):
        take = labels == label
        old[take] += (free[take] - old[take]).mean(axis=0)
    return old


class StopProbe(Exception):
    pass


class Probe(F6Simulator):
    def fields(self):
        super().fields()
        self.primary = self.force.numpy().copy()

    def launch(self, kernel, dim, args):
        if kernel is k.predict:
            parent = self.component_parent.numpy()
            roots = parent.copy()
            while np.any(roots != parent[roots]):
                roots = parent[roots]
            self.velocity = (
                args[1].numpy()
                + self.component_delta.numpy()[roots] / self.component_count.numpy()[roots, None]
            )
            raise StopProbe()
        return super().launch(kernel, dim, args)


def probe(x, event, obj, mode, flow_enabled=False, wall=False, gradient="consistent_trilinear"):
    p, b, r, g, f = obj
    sim = Probe(
        x,
        p,
        r,
        b,
        flow=f if flow_enabled else None,
        flow_enabled=flow_enabled,
        vessel=g,
        vessel_enabled=wall,
        mobility_mode=MODES[mode],
        gradient_mode=gradient,
    )
    ev = dict(event, flow_enabled=flow_enabled)
    ev.pop("mode", None)
    sim.control(ev)
    try:
        sim.step(0.0)
    except StopProbe:
        pass
    assert not sim.error.numpy()[0], sim.error.numpy()
    np.testing.assert_array_equal(x, sim.x.numpy())
    return dict(
        primary_N=sim.primary,
        secondary_N=sim.force.numpy() - sim.primary,
        flow_m_s=sim.flow_query.velocity.numpy() if flow_enabled else np.zeros_like(x),
        velocity_m_s=sim.velocity,
        free_m_s=sim.free.numpy(),
        pressure_Pa=sim.pressure.numpy(),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f7/protocol.json")
    obj = objects()
    p, b, r, _, _ = obj
    e0, _, ey = proto["environments"]
    cases = [
        (
            "mature",
            "results/fast_pzt/f6/navigation_v1/zero_flow_feedback/trajectory.h5",
            0,
            e0,
            proto["frozen_F6_events"][0],
            False,
        ),
        (
            "dispersed",
            "results/fast_pzt/f6/navigation_v1/Y_CFD_off/trajectory.h5",
            -1,
            ey,
            ey["events"]["unchanged"],
            True,
        ),
        (
            "near_wall",
            "results/fast_pzt/f6/navigation_v1/Y_CFD_unchanged/trajectory.h5",
            0,
            ey,
            ey["events"]["unchanged"],
            True,
        ),
    ]
    rows = {}
    arrays = {}
    math = []
    for name, path, frame, env, ev, flow in cases:
        with h5py.File(ROOT / path) as h:
            x = h["positions_m"][frame, 0]
            t = float(h["time_s"][frame])
        target = np.array(env["target_m"])
        direction = target - x.mean(axis=0)
        direction /= np.linalg.norm(direction)
        m = independent_rpy(x, p.radius_m, p.viscosity_Pa_s)
        np.testing.assert_allclose(m, matrix_numpy(x, p.radius_m, p.viscosity_Pa_s), rtol=2e-14, atol=1e-7)
        eig = np.linalg.eigvalsh(m)
        assert eig.min() > 0
        np.testing.assert_allclose(m, m.T, atol=1e-9)
        np.testing.assert_allclose(np.diag(m), 1 / (6 * np.pi * p.viscosity_Pa_s * p.radius_m), rtol=1e-14)
        math.append(
            dict(
                state=name,
                min_eigenvalue=eig.min(),
                max_eigenvalue=eig.max(),
                symmetry_max=float(abs(m - m.T).max()),
            )
        )
        values = {mode: probe(x, ev, obj, mode, flow, flow) for mode in MODES}
        field, _, _ = b.combine(ev["amplitude_Pa"], ev["phase_rad"], ev["enabled"])
        pp, gg = query_numpy(field.pressure_Pa, b.origin, b.spacing, x)
        vv = r.query(pp)
        primary = -0.5 * np.real(vv[:, 0, None] * gg.conj())
        secondary = pair_reference(x, vv, pp, p, r.query(np.array([r.knots[0] + 0j]))[0, 0] / r.knots[0])
        for mode, v in values.items():
            np.testing.assert_allclose(v["primary_N"], primary, rtol=2e-11, atol=1e-24)
            np.testing.assert_allclose(v["secondary_N"], secondary, rtol=2e-10, atol=1e-23)
        force = primary + secondary
        uf = values["M0"]["flow_m_s"]
        free = uf + force / (6 * np.pi * p.viscosity_Pa_s * p.radius_m)
        refs = {"M0": local_numpy(x, free, p), "M1": free, "M2": uf + (m @ force.ravel()).reshape(-1, 3)}
        rows[name] = dict(
            source=path,
            source_sha256=sha(ROOT / path),
            frame=frame,
            time_s=t,
            event=ev,
            target_m=target,
            models={},
        )
        arrays[name + "_positions_m"] = x
        arrays[name + "_target_m"] = target
        for mode, v in values.items():
            np.testing.assert_allclose(v["velocity_m_s"], refs[mode], rtol=5e-10, atol=1e-16)
            u = v["velocity_m_s"]
            mean = u.mean(axis=0)
            relative = u - mean
            dx = x - x.mean(axis=0)
            rows[name]["models"][mode] = dict(
                mean_velocity_m_s=mean,
                target_component_m_s=mean @ direction,
                transverse_speed_m_s=np.linalg.norm(mean - (mean @ direction) * direction),
                relative_rms_m_s=np.sqrt(np.mean(np.sum(relative**2, axis=1))),
                rg_rate_m_s=np.mean(np.sum(dx * relative, axis=1)) / np.sqrt(np.mean(np.sum(dx**2, axis=1))),
                cpu_gpu_max_m_s=np.abs(u - refs[mode]).max(),
            )
            for key, val in v.items():
                arrays[name + "_" + mode + "_" + key] = val
        for ma, mb in [("M0", "M1"), ("M1", "M2"), ("M0", "M2")]:
            diff = np.linalg.norm(values[ma]["velocity_m_s"] - values[mb]["velocity_m_s"], axis=1)
            i = int(diff.argmax())
            rows[name][ma + "_" + mb] = dict(
                max_difference_bubble_id=i, position_m=x[i], max_velocity_difference_m_s=diff[i]
            )
        old = probe(x, ev, obj, "M2", flow, flow, gradient="legacy_nodal")["velocity_m_s"].mean(axis=0)
        rows[name]["legacy_gradient_mean_velocity_m_s"] = old
    # Few meaningful regularization cases, not a random stress suite.
    for d in [0.0, 1.5 * p.radius_m, 2 * p.radius_m, 6 * p.radius_m]:
        x = np.zeros((1 if d == 0 else 2, 3))
        x[-1, 0] = d
        m = independent_rpy(x, p.radius_m, p.viscosity_Pa_s)
        np.testing.assert_allclose(m, matrix_numpy(x, p.radius_m, p.viscosity_Pa_s), rtol=1e-14)
        assert np.linalg.eigvalsh(m).min() > 0
    np.savez_compressed(out / "fixed_states.npz", **arrays)
    write(
        out / "CHECKS.json",
        dict(
            status="PASS",
            identity=identity(),
            states=rows,
            rpy_math=math,
            audit="Shared F6 fields/pairs/PBD/integration. M1/M2 intercept local_mobility and component_links; self mobility once; no old lubrication in M2. Instantaneous velocity before PBD, not post-contact displacement.",
            limitation="RPY is free-space; near-wall state checks implementation, not wall hydrodynamic accuracy.",
        ),
    )
    print("Fixed-state CPU/GPU and mode checks PASS", flush=True)


if __name__ == "__main__":
    main()
