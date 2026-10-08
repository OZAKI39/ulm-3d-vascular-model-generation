"""Original labelled Stage-4 mesh preparation and finite-radius initialization in SI.

VTK is used only before/after GPU integration. The prepared arrays retain every
original triangle; they are not an implicit-cylinder replacement.
"""

import hashlib
import numpy as np

ORIGINAL_SHA256 = "e14b06f1164be3beab261ec052dbe7edad6c2394d80077787cb318cbece15f4b"
LABELS = {1: "WALL", 2: "INLET", 3: "OUTLET_A", 4: "OUTLET_B"}


class VesselGeometry:
    def __init__(self, path):
        import pyvista as pv
        import vtk
        from mevo_mb.coupling.warp_wall import mesh_audit

        self.path = str(path)
        self.sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
        if self.sha256 != ORIGINAL_SHA256:
            raise ValueError("F4 requires the original archived Stage-4 mesh")
        self.surface = pv.read(path)
        self.vertices = np.asarray(self.surface.points, np.float64)
        self.faces = self.surface.faces.reshape(-1, 4)[:, 1:].astype(np.int32)
        self.labels = np.asarray(self.surface["ModelFaceID"], np.int32)
        self.audit = mesh_audit(self.surface)
        if not self.audit["closed_consistent_outward"] or set(self.labels) != set(LABELS):
            raise ValueError("Closed outward labelled lumen required")
        self.wall = self.surface.extract_cells(self.labels == 1).extract_surface().triangulate()
        self.closed_distance = vtk.vtkImplicitPolyDataDistance()
        self.closed_distance.SetInput(self.surface)
        self.wall_distance = vtk.vtkImplicitPolyDataDistance()
        self.wall_distance.SetInput(self.wall)

    def query(self, positions):
        """Signed WALL distance uses closed-lumen sign, never open-wall sign."""
        x = np.atleast_2d(positions).astype(float)
        distance, closest, normal, inside = [], [], [], []
        for p in x:
            q = [0.0, 0.0, 0.0]
            self.wall_distance.EvaluateFunctionAndGetClosestPoint(p, q)
            d = np.linalg.norm(p - q)
            flag = self.closed_distance.EvaluateFunction(p) <= 0
            sign = 1 if flag else -1
            distance.append(sign * d)
            closest.append(q)
            normal.append(sign * (np.asarray(q) - p) / max(d, 1e-30))
            inside.append(flag)
        return np.asarray(distance), np.asarray(closest), np.asarray(normal), np.asarray(inside)

    def prepare(self, path, planes):
        tri = self.vertices[self.faces]
        normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        names = ["INLET", "OUTLET_A", "OUTLET_B"]
        origin = self.vertices.mean(0)
        quantization = np.max(
            np.linalg.norm(
                (self.vertices - origin).astype(np.float32).astype(float) - (self.vertices - origin), axis=1
            )
        )
        np.savez_compressed(
            path,
            vertices_m=self.vertices,
            faces=self.faces,
            labels=self.labels,
            outward_normals=normal,
            local_origin_m=origin,
            cap_origins_m=np.array([planes[n]["origin_m"] for n in names]),
            cap_normals=np.array([planes[n]["normal"] for n in names]),
            original_sha256=np.array(self.sha256),
            broadphase_padding_m=np.array(max(2e-8, 16 * quantization)),
        )
        return dict(
            original_sha256=self.sha256,
            original_coordinate_dtype=str(self.surface.points.dtype),
            **self.audit,
            labels={LABELS[k]: int((self.labels == k).sum()) for k in LABELS},
            bounds_m=np.asarray(self.surface.bounds).reshape(3, 2),
            local_origin_m=origin,
            float32_vertex_quantization_max_m=float(quantization),
            geometry_scope="Original parametric smooth Y lumen; not patient-specific anatomy",
            physical_boundary="WALL only; caps are open terminals",
            units="m",
            prepared_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        )


def initialize_particles_in_vessel(
    geometry, field_bounds, target, count=32, seed=20261008, radius=1.78e-6, wall_gap=0.2e-6, pair_gap=0.2e-6
):
    """Rejection sampling near a real wall patch; never projects invalid centers."""
    _, closest, normals, _ = geometry.query(np.array(target)[None])
    anchor, normal = closest[0], normals[0]
    tangent = np.cross(normal, [1.0, 0.0, 0.0])
    if np.linalg.norm(tangent) < 0.1:
        tangent = np.cross(normal, [0.0, 1.0, 0.0])
    tangent /= np.linalg.norm(tangent)
    other = np.cross(normal, tangent)
    rng = np.random.default_rng(seed)
    accepted = []
    trials = 0
    while len(accepted) < count and trials < 20000:
        trials += 1
        a, b = rng.uniform(-18e-6, 18e-6, 2)
        if a * a + b * b > (18e-6) ** 2:
            continue
        depth = radius + wall_gap + rng.uniform(0, 7e-6)
        p = anchor + a * tangent + b * other - depth * normal
        if np.any(p < field_bounds[:, 0]) or np.any(p > field_bounds[:, 1]):
            continue
        distance, _, _, inside = geometry.query(p[None])
        if not inside[0] or distance[0] < radius + wall_gap:
            continue
        if accepted and np.min(np.linalg.norm(np.array(accepted) - p, axis=1)) < 2 * radius + pair_gap:
            continue
        accepted.append(p)
    if len(accepted) != count:
        raise RuntimeError("Requested wall patch cannot accommodate valid initialization")
    x = np.asarray(accepted)
    gaps = geometry.query(x)[0] - radius
    distances = np.linalg.norm(x[:, None] - x[None], axis=-1)
    np.fill_diagonal(distances, np.inf)
    return x, dict(
        seed=seed,
        trials=trials,
        wall_anchor_m=anchor,
        outward_normal=normal,
        tangent_1=tangent,
        tangent_2=other,
        count=count,
        radius_m=radius,
        minimum_wall_gap_m=float(gaps.min()),
        minimum_pair_gap_m=float(distances.min() - 2 * radius),
        bounds_m=np.column_stack([x.min(0), x.max(0)]),
        initialization="Finite-radius rejection sampling on actual mesh; no correction or projection",
    )
