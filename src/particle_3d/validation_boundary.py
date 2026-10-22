"""Validation-only segment/surface classification. No wall force or correction."""
from dataclasses import dataclass
import numpy as np
import pyvista as pv
import vtk
from .microbubble import vector3

EPS = np.finfo(np.float64).eps


@dataclass(frozen=True)
class BoundaryEvent:
    role: str
    segment_fraction: float
    position_m: np.ndarray
    triangle_id: int
    role_triangle_id: int
    simultaneous_roles: tuple[str, ...]


class ValidationBoundaryClassifier:
    """Every segment is tested against all five original boundary surfaces.

    VTK only supplies AABB candidates. Plane intersection and triangle weights
    are computed here. Geometry tolerance is 64*eps*coordinate_scale (m),
    with triangle Gram conditioning for its dimensionless weights. At a shared
    rim/tie, WALL takes priority conservatively. Coplanar contacts are reported.
    This class does not move, reflect, project or exert any force on a particle.
    """

    def __init__(self, boundaries):
        xyz, roles, local_ids = [], [], []
        for role, surface in sorted(boundaries.items()):
            faces = surface.faces.reshape(-1, 4)
            if not np.all(faces[:, 0] == 3):
                raise ValueError("boundary classifier requires triangles")
            triangles = np.asarray(surface.points)[faces[:, 1:]]
            xyz.extend(triangles); roles.extend([role] * len(triangles)); local_ids.extend(range(len(triangles)))
        self.xyz = np.asarray(xyz, dtype=np.float64)
        self.roles = np.asarray(roles)
        self.local_ids = np.asarray(local_ids)
        self.edge1 = self.xyz[:, 1] - self.xyz[:, 0]
        self.edge2 = self.xyz[:, 2] - self.xyz[:, 0]
        normal = np.cross(self.edge1, self.edge2)
        norm = np.linalg.norm(normal, axis=1)
        if np.any(norm == 0): raise ValueError("degenerate boundary triangle")
        self.normal = normal / norm[:, None]
        self.h = np.maximum(np.linalg.norm(self.edge1, axis=1), np.linalg.norm(self.edge2, axis=1))
        basis = np.stack((self.edge1, self.edge2), axis=2) / self.h[:, None, None]
        gram = basis.transpose(0, 2, 1) @ basis
        self.inverse = np.linalg.inv(gram) @ basis.transpose(0, 2, 1)
        condition = np.linalg.cond(gram, np.inf)
        self.coord = np.maximum(np.max(np.abs(self.xyz), axis=(1, 2)), self.h)
        self.tolerance_m = 64 * EPS * self.coord
        self.weight_tolerance = 64 * EPS * condition * (1 + self.coord / self.h)
        cells = np.column_stack((np.full(len(self.xyz), 3), np.arange(len(self.xyz) * 3).reshape(-1, 3)))
        self._surface = pv.PolyData(self.xyz.reshape(-1, 3), cells.ravel())
        self._locator = vtk.vtkStaticCellLocator()
        self._locator.SetDataSet(self._surface); self._locator.BuildLocator()

    def _weights(self, point, ids):
        tail = np.einsum("nij,nj->ni", self.inverse[ids], (point - self.xyz[ids, 0]) / self.h[ids, None])
        return np.column_stack((1 - tail.sum(1), tail))

    def first_event(self, start_m, end_m):
        start, end = vector3(start_m, "start_m"), vector3(end_m, "end_m")
        direction = end - start
        length = np.linalg.norm(direction)
        if length == 0: return None
        pad = float(self.tolerance_m.max())
        bounds = np.column_stack((np.minimum(start, end) - pad, np.maximum(start, end) + pad)).ravel()
        candidate = vtk.vtkIdList(); self._locator.FindCellsWithinBounds(bounds, candidate)
        ids = np.array([candidate.GetId(i) for i in range(candidate.GetNumberOfIds())], dtype=np.int64)
        if not len(ids): return None
        d0 = np.einsum("ij,ij->i", self.normal[ids], start - self.xyz[ids, 0])
        d1 = np.einsum("ij,ij->i", self.normal[ids], end - self.xyz[ids, 0])
        tol = self.tolerance_m[ids]
        selected = ((d0 >= -tol) & (d1 <= tol)) | ((d0 <= tol) & (d1 >= -tol))
        ids, d0, d1, tol = ids[selected], d0[selected], d1[selected], tol[selected]
        hits = []
        for cell, s0, s1, tolerance in zip(ids, d0, d1, tol):
            tau = self.weight_tolerance[cell]
            if abs(s0) <= tolerance and abs(s1) <= tolerance:
                # Closed triangle half-plane clipping of a coplanar segment,
                # only to determine first contact fraction, not to alter state.
                w0 = self._weights(start, np.array([cell]))[0]
                w1 = self._weights(end, np.array([cell]))[0]
                lower, upper = 0., 1.
                for a, slope in zip(w0, w1 - w0):
                    if abs(slope) <= tau:
                        if a < -tau: upper = -1.; break
                    elif slope > 0: lower = max(lower, -a / slope)
                    else: upper = min(upper, -a / slope)
                if lower > upper + tau: continue
                fraction = min(1., max(0., lower))
            else:
                fraction = float(s0 / (s0 - s1))
                t_tol = tolerance / max(abs(s0 - s1), np.finfo(float).tiny)
                if fraction < -t_tol or fraction > 1 + t_tol: continue
                fraction = min(1., max(0., fraction))  # roundoff-only endpoint representation
                weights = self._weights(start + fraction * direction, np.array([cell]))[0]
                if np.min(weights) < -tau or np.max(weights) > 1 + tau: continue
            hits.append((fraction, int(cell)))
        if not hits: return None
        first = min(t for t, _ in hits)
        tied = [(t, c) for t, c in hits if abs(t - first) <= 2 * pad / length]
        priority = lambda item: (0 if self.roles[item[1]] == "WALL" else 1 if self.roles[item[1]] == "INLET" else 2, self.roles[item[1]], item[1])
        fraction, cell = min(tied, key=priority)
        return BoundaryEvent(str(self.roles[cell]), fraction, vector3(start + fraction * direction, "event_position"),
                             cell, int(self.local_ids[cell]), tuple(sorted({str(self.roles[c]) for _, c in tied})))
