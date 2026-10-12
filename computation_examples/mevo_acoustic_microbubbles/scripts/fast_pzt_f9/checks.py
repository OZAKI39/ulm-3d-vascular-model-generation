"""Bounded mathematical checks, real-Y noninterference and failure retention."""

import argparse
import numpy as np
import h5py
import warp as wp
from mevo_mb.fast_sim.f9_ledger import AuditSimulator
from mevo_mb.fast_sim.f8_gate import GateSimulator
from mevo_mb.fast_sim.f4_engine import GPUVessel
from mevo_mb.fast_sim.flow_gpu import AffineFlow
from .common import BASE, ROOT, js, write, identity, objects
from .validation import validate_ledger
from .runner import run_case


def box_geometry(radius):
    # Explicit synthetic planar-wall unit fixture. Original Y geometry is untouched.
    import pyvista as pv

    mesh = pv.Box(bounds=(-30 * radius, 30 * radius, -30 * radius, 30 * radius, 0, 30 * radius)).triangulate()
    vertices = np.asarray(mesh.points, dtype=float)
    faces = mesh.faces.reshape(-1, 4)[:, 1:].astype(np.int32)
    g = GPUVessel.__new__(GPUVessel)
    g.device = "cuda:0"
    g.origin = wp.vec3d(0.0)
    g.padding = 2e-8
    g.vertices = wp.array(vertices, dtype=wp.vec3d, device=g.device)
    g.faces = wp.array(faces, dtype=int, device=g.device)
    g.wall_faces = g.faces
    g.tags = wp.ones(len(faces), dtype=int, device=g.device)
    g.points32 = wp.array(vertices.astype(np.float32), dtype=wp.vec3, device=g.device)
    g.indices = wp.array(faces.ravel(), dtype=int, device=g.device)
    g.closed = wp.Mesh(points=g.points32, indices=g.indices)
    g.wall = g.closed
    g.cap_origins = wp.array([[1.0, 0, 0]] * 3, dtype=wp.vec3d, device=g.device)
    g.cap_normals = wp.array([[1.0, 0, 0]] * 3, dtype=wp.vec3d, device=g.device)
    g.metadata = dict(fixture="synthetic closed box with planar lower wall; not a new vessel model")
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    write(out / "REGISTERED_PROTOCOL.json", js(ROOT / "configs/fast_pzt/f9/protocol.json"))
    p, basis, response, geometry, flow = objects()
    rows = {}
    event = dict(
        time_s=0.0,
        amplitude_Pa=[0.0] * 4,
        phase_rad=[0.0] * 4,
        enabled=[False] * 4,
        flow_enabled=True,
        secondary_enabled=False,
    )
    for name in ["free", "plane"]:
        x = np.array([[0.0, -4 * p.radius_m, 5 * p.radius_m], [0.0, 4 * p.radius_m, 5 * p.radius_m]])
        velocity = np.array([1e-5, 0.0, 0.0])
        g = None
        if name == "plane":
            x[:, 2] = p.radius_m + 0.5e-9
            velocity[2] = -1e-4
            g = box_geometry(p.radius_m)
        provider = AffineFlow(velocity, np.zeros((3, 3)), np.zeros(3))
        results = []
        for audit in [False, True]:
            cls = AuditSimulator if audit else GateSimulator
            extra = dict(audit_enabled=True) if audit else {}
            sim = cls(
                x,
                p,
                response,
                basis=None,
                flow=provider,
                flow_enabled=True,
                acoustic_enabled=False,
                vessel=g,
                vessel_enabled=g is not None,
                **extra,
            )
            r = sim.run(1e-4, 1e-5, [event], sample_interval_s=5e-5)
            assert r["summary"]["status"] == "COMPLETE", r["summary"]
            results.append(r)
        assert np.array_equal(results[0]["positions_m"], results[1]["positions_m"])
        r = results[1]
        check = validate_ledger(r["motion_ledger"], r["positions_m"], r["time_s"])
        expected = x + velocity * 1e-4
        if name == "plane":
            expected[:, 2] = p.radius_m
        np.testing.assert_allclose(r["positions_m"][-1, 0], expected, atol=1e-18, rtol=0)
        ev = r["motion_ledger"]["wall_contact_events"]
        if name == "free":
            assert len(ev) == 0
        else:
            assert len(ev) > 0 and np.all(ev[:, 13] > 0) and np.max(np.abs(ev[:, 11:13])) < 1e-18
        with h5py.File(out / (name + ".h5"), "w") as h:
            h["positions_m"] = r["positions_m"]
            h["time_s"] = r["time_s"]
            for key, v in r["motion_ledger"].items():
                if key != "rejected_intervals":
                    h.create_dataset("motion_ledger/" + key, data=v, compression="gzip")
        rows[name] = dict(**check, audit_bitwise_equal=True, known_solution_pass=True)
    f8 = ROOT / "results/fast_pzt/f8/y_formal_v1/Y1"
    cfg = js(f8 / "CONFIG.json")
    proto = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    for audit in [False, True]:
        s = run_case(
            out / ("Y_short_on" if audit else "Y_short_off"),
            cfg["environment"],
            cfg["events"],
            proto,
            (p, basis, response, geometry, flow),
            dt=cfg["dt_s"],
            duration=0.01,
            audit=audit,
        )
        assert s["numerical_quality"] == "PASS"
    with (
        h5py.File(out / "Y_short_on/trajectory.h5") as on,
        h5py.File(out / "Y_short_off/trajectory.h5") as off,
        h5py.File(f8 / "trajectory.h5") as old,
    ):
        x = on["positions_m"][:]
        assert np.array_equal(x, off["positions_m"][:])
        assert np.array_equal(x, old["positions_m"][: len(x)])
    rows["real_Y"] = dict(
        audit_on_off_bitwise_equal=True, F8_first_21_frames_bitwise_equal=True, duration_s=0.01
    )
    write(out / "CHECKS.json", dict(status="PASS", checks=rows, identity=identity()))
    print("F9 independent and Y noninterference checks PASS", flush=True)


if __name__ == "__main__":
    main()
