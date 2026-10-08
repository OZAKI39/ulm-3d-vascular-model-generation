"""Small finite-edge and terminal-rim geometry checks, no acoustic trajectory campaign."""

import argparse
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f4_engine import GPUVessel, F4VesselSimulator
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.vessel_kernels import swept_distance2
from mevo_mb.bubbles.geometry_exact_r3 import segment_triangle
from .common import BASE, write, identity


@wp.kernel
def distance_test(data: wp.array2d(dtype=wp.vec3d), result: wp.array(dtype=wp.float64)):
    i = wp.tid()
    result[i] = swept_distance2(data[i, 0], data[i, 1], data[i, 2], data[i, 3], data[i, 4])


def main(label):
    out = BASE / label
    out.mkdir(exist_ok=False)
    g = GPUVessel(BASE / "inputs/prepared_geometry.npz", BASE / "inputs/labelled_lumen_m.vtp")
    prep = np.load(BASE / "inputs/prepared_geometry.npz")
    vertices = prep["vertices_m"]
    faces = prep["faces"]
    region = np.load(BASE / "inputs/parent_initial.npz")["positions_m"].reshape(-1, 3).mean(axis=0)
    face = int(g.query(region[None])["face"][0])
    tri = vertices[faces[prep["labels"] == 1][face]]
    a, b, c = tri
    center = tri.mean(axis=0)
    normal = np.cross(b - a, c - a)
    normal /= np.linalg.norm(normal)
    tangent = (b - a) / np.linalg.norm(b - a)
    edge = (a + b) / 2
    paths = [
        (center - normal * 2e-6, center + normal * 2e-6),
        (center + normal * 2e-6 - tangent * 1e-6, center + normal * 2e-6 + tangent * 1e-6),
        (edge + normal * 2e-6, edge + normal * 3e-6),
        (a + normal * 2e-6, a + normal * 2e-6),
        (a - tangent * 2e-6 + normal * 1e-6, a - tangent * 2e-6 - normal * 1e-6),
    ]
    data = np.array([[p, q, a, b, c] for p, q in paths])
    gpu = wp.array(data, dtype=wp.vec3d, device=g.device)
    output = wp.empty(len(paths), dtype=wp.float64, device=g.device)
    wp.launch(distance_test, dim=len(paths), inputs=[gpu, output], device=g.device)
    reference = np.array([float(segment_triangle(*row)["distance2"]) for row in data])
    measured = output.numpy()
    errors = abs(np.sqrt(reference) - np.sqrt(measured))
    assert errors.max() < 1e-12, errors
    np.savez_compressed(
        out / "finite_triangle_paths.npz",
        data_m=data,
        reference_distance2_m2=reference,
        gpu_distance2_m2=measured,
    )
    p = F2Physics()
    response = ResponseSection(BASE / "inputs/response.h5", p)
    capfaces = faces[prep["labels"] == 2]
    rim = vertices[np.unique(capfaces)]
    origin = prep["cap_origins_m"][0]
    n = prep["cap_normals"][0]
    radial = rim - origin
    radial -= np.outer(radial @ n, n)
    v = radial[np.argmax(np.linalg.norm(radial, axis=1))]
    radial = v / np.linalg.norm(v)
    start = origin + v - radial * (p.radius_m + 1e-6) - n * 100e-9
    q = g.query(start[None])
    assert q["inside"].all() and q["distance_m"][0] > p.radius_m
    # A cap-plane crossing that simultaneously cuts the real rim must be rejected,
    # not turned into a disappearing particle/terminal event.
    velocity = (n * 200e-9 + radial * 5 * p.radius_m) / 1e-5
    sim = F4VesselSimulator(start[None], p, response, vessel=g, geometry_unit_test_velocity=velocity)
    run = sim.run(1e-5, 1e-5, [dict(time_s=0, pressure_Pa=0)], 1e-5)
    assert run["summary"]["status"] == "FAILED"
    assert len(run["boundary_events"]) == 0 and run["active"][-1].all()
    np.testing.assert_array_equal(run["positions_m"][-1, 0, 0], start)
    write(
        out / "RESULT.json",
        dict(
            status="PASS",
            identity=identity(),
            scope="GEOMETRY_UNIT_TEST",
            finite_path_count=len(paths),
            max_distance_error_m=float(errors.max()),
            reference="Existing exact rational segment_triangle on original binary64 input coordinates",
            rim_crossing_rejected_without_exit=True,
            rim_rejection_summary=run["summary"],
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    main(ap.parse_args().label)
