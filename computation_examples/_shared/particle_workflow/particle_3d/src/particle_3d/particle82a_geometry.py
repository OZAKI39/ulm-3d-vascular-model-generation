"""P8.2A diagnostics: open aperture is not a solid cap.

Distances are to finite mesh features in SI. A full-domain representation adds
open-cap distance constraints, explicitly separate from the P7 wall contract.
No force, radius, field, or historical admission implementation is modified.
"""
from dataclasses import dataclass
import numpy as np
from .wall_geometry import WallGeometry
from .nearfield_regularization import NearFieldRegularizationV1

FAILURES = ('WALL_INTERSECTION', 'INLET_PERIMETER_CLEARANCE_FAIL',
    'OPEN_CAP_ONLY_OUTSIDE_DOMAIN', 'NEARFIELD_HANDOFF_VIOLATION',
    'POINT_LOCATION_FAIL', 'NO_OWNER_TETRA', 'PAIR_CONFLICT',
    'NUMERICAL_GEOMETRY_UNRESOLVED', 'OTHER')


def segment_distance(point, segments):
    p = np.asarray(point); a = segments[:, 0]; v = segments[:, 1]-a
    t = np.clip(np.einsum('ij,ij->i', p-a, v)/np.einsum('ij,ij->i', v, v), 0., 1.)
    return float(np.linalg.norm(p-a-t[:, None]*v, axis=1).min())


def perimeter_edges(points, faces):
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    if np.any(counts > 2):
        raise ValueError('Non-manifold inlet edge')
    edge_ids = unique[counts == 1]
    if not len(edge_ids):
        raise ValueError('No open aperture perimeter')
    vertices, degree = np.unique(edge_ids, return_counts=True)
    if np.any(degree != 2):
        raise ValueError('Inlet perimeter is not a collection of closed loops')
    return np.asarray(points)[edge_ids], edge_ids


def lower_gap(radius):
    return np.maximum(2e-9, 1e-3*np.asarray(radius))


def maximum_handoff_radius(clearance):
    # Solve a + max(h0, chi*a) <= clearance exactly, without clipping raw gaps.
    c = np.asarray(clearance)
    return np.maximum(0., np.minimum(c-2e-9, c/1.001))


def classify_failure(*, radius, wall_distance, perimeter_distance=None,
                     owner=True, finite=True, full_domain=False, cap_distance=np.inf,
                     pair_conflict=False, tolerance=0.):
    if not finite:
        return 'NUMERICAL_GEOMETRY_UNRESOLVED'
    if not owner:
        return 'NO_OWNER_TETRA'
    if pair_conflict:
        return 'PAIR_CONFLICT'
    if perimeter_distance is not None and perimeter_distance-radius < -tolerance:
        return 'INLET_PERIMETER_CLEARANCE_FAIL'
    if wall_distance-radius < -tolerance:
        return 'WALL_INTERSECTION'
    if wall_distance-radius-lower_gap(radius) < -tolerance:
        return 'NEARFIELD_HANDOFF_VIOLATION'
    if full_domain and cap_distance-radius < -tolerance:
        return 'OPEN_CAP_ONLY_OUTSIDE_DOMAIN'
    return 'ACCEPTED'


class InletGeometryAudit:
    def __init__(self, env):
        self.env = env
        surf = env.boundaries['INLET']
        self.points = np.asarray(surf.points, float)
        self.faces = surf.faces.reshape(-1, 4)[:, 1:]
        self.perimeter, self.edge_ids = perimeter_edges(self.points, self.faces)
        self.cap = WallGeometry(self.points[self.faces], provenance={'role': 'OPEN_CAP_DISTANCE_ONLY_NOT_SOLID'})
        other = [np.asarray(s.points)[s.faces.reshape(-1, 4)[:, 1:]]
                 for name, s in env.boundaries.items() if name.startswith('OUTLET_')]
        self.other_caps = WallGeometry(np.concatenate(other))
        self.origin = self.points.mean(axis=0)
        _, _, vh = np.linalg.svd(self.points-self.origin, full_matrices=False)
        self.basis = vh[:2]; self.normal = vh[2]
        mean_v = env.field.velocity_nodes_m_s[np.asarray(surf.point_data['GlobalNodeID'], int)-1].mean(axis=0)
        if self.normal@mean_v < 0:
            self.normal *= -1
        self.area = env.audit['INLET']['area_m2']
        self.equivalent_diameter = 2*np.sqrt(self.area/np.pi)
        # Search neighborhood is declared from aperture geometry, not MB radius.
        self.search_horizon = 2*self.equivalent_diameter

    def coordinates(self, points):
        return (np.asarray(points)-self.origin)@self.basis.T

    def distances(self, point, full=False):
        wall = self.env.wall.nearest_center_triangle(point)[1]
        perimeter = segment_distance(point, self.perimeter)
        cap = self.cap.nearest_center_triangle(point)[1] if full else 0.
        return wall, perimeter, cap

    def aperture_radius(self, point):
        wall, edge, _ = self.distances(point)
        return min(wall, edge)

    def full_margin(self, point, radius):
        if self.env.field.locate(point)[0] < 0:
            return -np.inf
        wall = self.env.wall.nearest_center_triangle(point)[1]-radius-lower_gap(radius)
        cap = self.cap.nearest_center_triangle(point)[1]-radius
        other = self.other_caps.nearest_center_triangle(point)[1]-radius
        return float(min(wall, cap, other))

    def semantics(self):
        wall_nodes = set(self.env.wall.global_node_ids.ravel().tolist())
        cap_ids = np.asarray(self.env.boundaries['INLET'].point_data['GlobalNodeID'], int)-1
        wall_faces = {tuple(sorted(x)) for x in self.env.wall.global_node_ids}
        cap_faces = {tuple(sorted(x)) for x in cap_ids[self.faces]}
        return dict(solid_boundaries=self.env.wall.provenance['solid_boundaries'],
            inlet_cap_wall_face_overlap=len(wall_faces & cap_faces),
            perimeter_vertices_on_wall=all(int(cap_ids[i]) in wall_nodes for i in np.unique(self.edge_ids)),
            cap_planarity_max_residual_m=float(np.max(np.abs((self.points-self.origin)@self.normal))),
            inlet_area_m2=self.area, equivalent_diameter_m=self.equivalent_diameter,
            search_horizon_m=self.search_horizon,
            search_horizon_rule='TWO_INLET_EQUIVALENT_DIAMETERS_PREDECLARED_GEOMETRY_ONLY',
            current_checker_requires_center_inside=True,
            current_checker_requires_whole_sphere_inside=False,
            current_checker_cap_is_solid=False,
            aperture_definition='MINIMUM_DISTANCE_TO_TRUE_WALL_AND_OPEN_APERTURE_PERIMETER',
            full_domain_definition='CENTER_INSIDE_AND_WALL_GAP_GE_HANDOFF_AND_DISTANCE_TO_EVERY_OPEN_CAP_GE_RADIUS',
            crossing_limitation='LEGAL_BIRTH_CENTER_DOES_NOT_CERTIFY_FINITE_SPHERE_SWEEP_FROM_ANCHOR')


def first_admissible(path, margin, *, horizon, tolerance=1e-9, max_evaluations=4096):
    """Find earliest legal point on saved streamline polyline, with bounds.

The signed distance margin is 1-Lipschitz on each inside-lumen segment. An
interval is discarded only if its Lipschitz upper bound is negative. Search is
left-first; ambiguous intervals at tolerance are reported, never skipped.
"""
    path = np.asarray(path, float)
    if path.ndim != 2 or path.shape[1] != 4 or not len(path):
        raise ValueError('Expected [time,x,y,z] path')
    calls = 0
    def evaluate(p):
        nonlocal calls
        calls += 1
        if calls > max_evaluations:
            raise RuntimeError('SEARCH_EVALUATION_GUARD')
        return float(margin(p))
    def bracket(a, b, ga, gb, sa, sb, ta, tb):
        # A nonfinite margin invalidates the Lipschitz certificate.
        if not np.isfinite([ga, gb]).all():
            raise RuntimeError('NONFINITE_GEOMETRY_MARGIN')
        if ga >= 0:
            return (a, sa, ta, ga, sa, sa)
        if (ga+gb+(sb-sa))/2 < 0:
            return None
        if sb-sa <= tolerance:
            if gb >= 0:
                return (b, sb, tb, gb, sa, sb)
            raise RuntimeError('SUBTOLERANCE_ADMISSIBILITY_UNRESOLVED')
        mid = (a+b)/2; sm = (sa+sb)/2; tm = (ta+tb)/2; gm = evaluate(mid)
        left = bracket(a, mid, ga, gm, sa, sm, ta, tm)
        return left if left is not None else bracket(mid, b, gm, gb, sm, sb, tm, tb)
    s = 0.; a = path[0, 1:]; ta = path[0, 0]
    try:
        ga = evaluate(a)
        for row in path[1:]:
            b = row[1:]; tb = row[0]; ds = float(np.linalg.norm(b-a))
            if ds == 0:
                continue
            if s+ds > horizon:
                f = (horizon-s)/ds; b = a+f*(b-a); tb = ta+f*(tb-ta); ds = horizon-s
            gb = evaluate(b)
            hit = bracket(a, b, ga, gb, s, s+ds, ta, tb)
            if hit is not None:
                point, distance, time, gap, lo, hi = hit
                return dict(status='ACCEPTED', center_m=point.tolist(), s_birth_m=distance,
                    entry_transition_time_s=time, full_margin_m=gap, bracket_m=[lo, hi],
                    evaluations=calls, minimality_tolerance_m=tolerance)
            s += ds; a = b; ta = tb; ga = gb
            if s >= horizon:
                break
        return dict(status='NO_ADMISSIBLE_INWARD_LOCATION', evaluations=calls, searched_distance_m=s)
    except RuntimeError as exc:
        return dict(status='NUMERICAL_GEOMETRY_UNRESOLVED', detail=str(exc), evaluations=calls,
                    searched_distance_m=s)


def cylinder_triangles(radius=5e-6, length=20e-6, sides=128):
    theta = np.linspace(0, 2*np.pi, sides, endpoint=False)
    ring = np.column_stack([radius*np.cos(theta), radius*np.sin(theta), np.zeros(sides)])
    triangles = []; cap = []
    for i in range(sides):
        a = ring[i]; b = ring[(i+1) % sides]; c = b+[0, 0, length]; d = a+[0, 0, length]
        triangles.extend([[a, b, c], [a, c, d]])
        cap.append([[0, 0, 0], b, a])
    return np.asarray(triangles), np.asarray(cap)
