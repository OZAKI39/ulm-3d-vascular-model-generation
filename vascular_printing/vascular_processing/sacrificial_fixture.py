"""Local box/port design around a frozen mesh. No upstream modelling is imported.

CadQuery/OCC constructs solids and splines; Manifold3D performs mesh Booleans;
FCL computes triangle-surface distances. SWC supplies endpoint identity; the
provenance-matched cap can supply the actual attachment contour and normal.
All coordinates are in the already selected print frame, in millimetres.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
import itertools
import json
import math
from pathlib import Path

import cadquery as cq
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
import trimesh
import yaml

from .swc_export import read_source
from .topbrain_qc import sha256

ROOT = Path(__file__).resolve().parents[1]
SIDE_FACES = ('-X', '+X', '-Y', '+Y')
FACES = SIDE_FACES
NORMALS = {name: np.eye(3)[axis] * sign for name, axis, sign in
           [('-X', 0, -1), ('+X', 0, 1), ('-Y', 1, -1), ('+Y', 1, 1)]}


def unit(value):
    value = np.asarray(value, dtype=float)
    length = np.linalg.norm(value)
    if not np.isfinite(length) or length <= 1e-12:
        raise ValueError('Undefined geometric tangent')
    return value / length


def angle(a, b):
    return float(np.degrees(np.arccos(np.clip(unit(a) @ unit(b), -1, 1))))


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows({k: '' if isinstance(v, (float, np.floating)) and not math.isfinite(v) else v
                          for k, v in row.items()} for row in rows)


def load_config(path):
    cfg = yaml.safe_load(Path(path).read_text(encoding='utf-8'))
    for group in ('box', 'ports', 'geometry'):
        for key, value in cfg[group].items():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                if not math.isfinite(value) or value <= 0:
                    raise ValueError(f'{group}.{key} must be positive and finite')
    if cfg['box']['open_top'] is not True:
        raise ValueError('This stage supports an open-top box only')
    if cfg['routing']['max_outlet_faces'] not in (1, 2):
        raise ValueError('At most two adjacent outlet faces are supported')
    faces = enabled_faces(cfg)
    if not faces or len(set(faces)) != len(faces) or any(face not in FACES for face in faces):
        raise ValueError('Only +/-X and +/-Y side-wall ports are allowed; Z ports are prohibited')
    fallback = cfg.get('clearance_fallback', {})
    if fallback.get('enabled'):
        if not fallback.get('authorization') or not cfg.get('native_restoration', {}).get('require_original_topology_audit_before_relaxation'):
            raise ValueError('Clearance fallback requires authorization and original topology audit')
        step, floor = fallback['decrement_mm'], fallback['minimum_mm']
        if not math.isfinite(step) or step <= 0 or not math.isfinite(floor) or not 0 <= floor < cfg['ports']['route_clearance_mm']:
            raise ValueError('Invalid clearance fallback search range')
        if fallback.get('retain_port_to_port_clearance') is not True:
            raise ValueError('Port-to-port clearance must remain unchanged')
    attachment = cfg.get('attachment_alignment', {})
    if attachment:
        for key, value in attachment.items():
            if isinstance(value, (int, float)) and (not math.isfinite(value) or value <= 0):
                raise ValueError('Attachment alignment values must be positive and finite')
        if not 0 < attachment['collar_inner_radius_fraction'] < 1 or attachment['section_offset_mm'] >= attachment['straight_lead_mm']:
            raise ValueError('Invalid internal collar or attachment inspection offset')
    return cfg


def enabled_faces(cfg):
    # Old configurations retain their original four-wall search.
    return tuple(cfg['routing'].get('allowed_faces', SIDE_FACES))


def boundary_kind(face):
    face_axis(face)
    return 'SIDE_WALL'


def boundary_thickness(face, cfg):
    face_axis(face)
    return cfg['box']['wall_thickness_mm']


@dataclass
class Endpoint:
    endpoint_id: str
    role: str
    swc_id: int
    original_swc_id: int
    branch_id: str
    position: np.ndarray
    radius: float
    tangent: np.ndarray
    tangent_window_mm: float
    sample_count: int
    original_radius_mm: float
    compensated_radius_mm: float
    neighborhood_ids: list[int] = field(default_factory=list)
    cap_ring: np.ndarray | None = None
    attachment_qc: dict = field(default_factory=dict)
    attachment_collar: object = None


@dataclass
class BoxBounds:
    inner: np.ndarray
    outer: np.ndarray

    @property
    def inner_size(self):
        return self.inner[1] - self.inner[0]

    @property
    def outer_size(self):
        return self.outer[1] - self.outer[0]


@dataclass
class Route:
    endpoint: Endpoint
    face: str
    kind: str
    target: np.ndarray
    end: np.ndarray
    direction: np.ndarray
    wire: cq.Wire
    shape: cq.Shape
    mesh: trimesh.Trimesh
    points: np.ndarray
    length: float
    turn: float
    minimum_bend_radius: float
    edge_clearance: float
    vessel_clearance: float = math.inf
    score: float = math.inf
    option_id: str = ''
    nearest_other_port: float = math.inf
    wall_spacing: float = math.inf
    required_vessel_clearance: float = 3.0
    preferred_vessel_clearance: float = 3.0


def endpoint_tangent(graph, node, cfg):
    """Arc-length regression over a terminal branch, with a clipped 5 mm window.

    Uniform arc-length sampling avoids bias from uneven native node spacing.
    The number of original samples, not interpolated samples, is reported.
    """
    inlet = graph.in_degree(node) == 0
    ids, pts, distances = [node], [graph.nodes[node]['coords'][:3]], [0.0]
    window = cfg['routing']['tangent_window_mm']
    while distances[-1] < window:
        neighbors = list(graph.successors(ids[-1]) if inlet else graph.predecessors(ids[-1]))
        if len(neighbors) != 1:
            break
        nxt = neighbors[0]
        point = graph.nodes[nxt]['coords'][:3]
        length = float(np.linalg.norm(point - pts[-1]))
        if length <= 1e-12:
            raise ValueError('Duplicate SWC centerline positions')
        ids.append(nxt)
        pts.append(point)
        distances.append(distances[-1] + length)
        if graph.out_degree(nxt) != 1 or graph.in_degree(nxt) != 1:
            break
    if len(ids) < 2:
        raise ValueError(f'Insufficient samples at SWC endpoint {node}')
    actual = min(window, distances[-1])
    sample_s = np.linspace(0, actual, max(3, math.ceil(actual / cfg['routing']['tangent_resample_step_mm']) + 1))
    points = np.asarray(pts)
    resampled = np.column_stack([np.interp(sample_s, distances, points[:, k]) for k in range(3)])
    slope = np.linalg.lstsq(np.column_stack([sample_s, np.ones(len(sample_s))]), resampled, rcond=None)[0][0]
    return unit(-slope), actual, len(ids), ids


def decode_endpoints(graph, mapping, edge_rows, transform, cfg):
    if not nx.is_arborescence(graph):
        raise ValueError('Companion SWC must be a single directed tree')
    node_rows = {int(row['swc_id']): row for row in mapping}
    branches = {(int(r['parent_id']), int(r['child_id'])): r['original_branch'] for r in edge_rows}
    endpoints = []
    outlet_index = 0
    for n in sorted(graph):
        inlet = graph.in_degree(n) == 0
        if not inlet and graph.out_degree(n) != 0:
            continue
        if inlet and graph.out_degree(n) != 1:
            raise ValueError('Inlet must have one incident branch')
        edge = (n, next(graph.successors(n))) if inlet else (next(graph.predecessors(n)), n)
        t, length, count, ids = endpoint_tangent(graph, n, cfg)
        p = graph.nodes[n]['coords']
        row = node_rows[n]
        if not inlet:
            outlet_index += 1
        endpoints.append(Endpoint('I1' if inlet else f'O{outlet_index}', 'inlet' if inlet else 'outlet',
            n, int(row['original_swc_id']), branches[edge], transform[:3, :3] @ p[:3] + transform[:3, 3],
            float(p[3]), unit(transform[:3, :3] @ t), length, count,
            float(row['original_radius_mm']), float(row['manufacturing_radius_mm']), ids))
    endpoints.sort(key=lambda e: (e.role != 'inlet', e.swc_id))
    return endpoints, branches


def load_inputs(cfg):
    base = (ROOT / cfg['input']['compact_root']).resolve()
    paths = {key: base / value for key, value in cfg['input'].items() if key not in ('compact_root', 'candidate')}
    manifest = json.loads(paths['manifest'].read_text())
    transform = np.asarray(json.loads(paths['transform'].read_text())['transform_4x4'])
    np.testing.assert_allclose(transform[:3, :3].T @ transform[:3, :3], np.eye(3), atol=1e-8)
    if not np.isclose(np.linalg.det(transform[:3, :3]), 1):
        raise ValueError('Saved orientation must be a rigid rotation, without scaling/reflection')
    compact = read_source(paths['compact_swc']).graph
    graph = read_source(paths['fitted_swc']).graph
    if set(graph.edges) != set(compact.edges) or set(graph) != set(compact):
        raise ValueError('Fitted/compact SWC topology mismatch')
    mapping = read_csv(paths['node_mapping'])
    for n in graph:
        if sum(int(row['swc_id']) == n for row in mapping) != 1:
            raise ValueError('CAP_PROVENANCE_UNRESOLVED_NODE_MAPPING_' + str(n))
    endpoints, branches = decode_endpoints(graph, mapping, read_csv(paths['edge_mapping']), transform, cfg)
    candidate = manifest['candidates'][cfg['input']['candidate']]
    caps = candidate['surface_qc']['port_caps']
    if len(caps) != len(endpoints) or Path(manifest['selected_stl']).name != paths['source_stl'].name:
        raise ValueError('Frozen STL selection or endpoint count does not match the manifest')
    np.testing.assert_allclose(transform, manifest['selected_orientation']['transform_4x4'], atol=1e-10)
    cap_matches = {}
    for endpoint in endpoints:
        native = graph.nodes[endpoint.swc_id]['coords']
        matches = [j for j, c in enumerate(caps)
            if np.linalg.norm(native[:3] - c['center_mm']) < cfg['geometry']['coordinate_tolerance_mm']
            and abs(native[3] - c['median_radius_mm']) < cfg['geometry']['coordinate_tolerance_mm']]
        if len(matches) != 1 or matches[0] in cap_matches.values():
            raise ValueError('CAP_PROVENANCE_UNRESOLVED_' + endpoint.endpoint_id)
        cap_matches[endpoint.endpoint_id] = matches[0]
    mesh = trimesh.load_mesh(paths['source_stl'], process=True)
    if not mesh.is_volume:
        raise ValueError('SOURCE_CORE_NOT_CLOSED_POSITIVE_VOLUME')
    aligned_ids = set(cfg.get('attachment_alignment', {}).get('endpoint_ids', []))
    if aligned_ids - {e.endpoint_id for e in endpoints}:
        raise ValueError('Unknown endpoint in attachment_alignment.endpoint_ids')
    for endpoint in endpoints:
        if endpoint.endpoint_id in aligned_ids:
            cap = caps[cap_matches[endpoint.endpoint_id]]
            align_known_cap(endpoint, mesh, cap, cfg)
    return dict(base=base, paths=paths, manifest=manifest, transform=transform, graph=graph,
                branches=branches, endpoints=endpoints, mesh=mesh)


def polygon_wire(points):
    return cq.Wire.makePolygon([cq.Vector(*p) for p in points], close=True)


def align_known_cap(endpoint, source, cap_record, cfg):
    # Compatibility entry: a single implementation shared by every port.
    from .port_attachment_alignment import align_known_cap as align
    try:
        return align(endpoint, source, cap_record, cfg)
    except ValueError as exc:
        if 'KNOWN_CAP_' in str(exc):
            raise ValueError('CAP_PROVENANCE_UNRESOLVED_' + endpoint.endpoint_id + ': ' + str(exc)) from exc
        raise


def attachment_lead(endpoint, cfg):
    return cfg['attachment_alignment']['straight_lead_mm'] if endpoint.cap_ring is not None else 0.0


def attachment_section_qc(route, cfg):
    """Check the generated triangle mesh against the real cap contour.

    This catches swept-profile translation/rotation errors even when the result
    remains closed, connected and positive-volume.
    """
    from shapely.geometry import LineString, Polygon
    e = route.endpoint
    if e.cap_ring is None:
        return {'evaluated': False}
    offset = cfg['attachment_alignment']['section_offset_mm']
    origin = e.position + offset * e.tangent
    section = route.mesh.section(plane_origin=origin, plane_normal=e.tangent)
    paths = [] if section is None else section.discrete
    if len(paths) != 1 or np.linalg.norm(paths[0][0] - paths[0][-1]) > cfg['geometry']['coordinate_tolerance_mm']:
        raise ValueError('ATTACHMENT_SECTION_NOT_ONE_CLOSED_CONTOUR')
    u = unit(e.cap_ring[0] - e.position)
    basis = np.array([u, np.cross(e.tangent, u)]).T
    expected = (e.cap_ring - e.position) @ basis
    actual = (paths[0] - origin) @ basis
    contour_error = LineString(actual).hausdorff_distance(LineString(np.vstack([expected, expected[0]])))
    expected_center = np.array(Polygon(expected).centroid.coords[0])
    actual_center = np.array(Polygon(actual).centroid.coords[0])
    center_error = float(np.linalg.norm(actual_center - expected_center))
    axis_error = angle(route.wire.tangentAt(0).toTuple(), e.tangent)
    tolerance = cfg['attachment_alignment']['section_match_tolerance_mm']
    passed = contour_error <= tolerance and center_error <= tolerance and axis_error < 1e-4
    return dict(evaluated=True, passed=bool(passed), section_offset_mm=offset,
        center_offset_mm=center_error, contour_hausdorff_mm=float(contour_error),
        axis_error_deg=axis_error, contour_count=len(paths), expected_area_mm2=float(Polygon(expected).area),
        measured_area_mm2=float(Polygon(actual).area), tolerance_mm=tolerance)


def compute_box(bounds, cfg):
    b = cfg['box']
    inner = np.asarray(bounds, dtype=float).copy()
    inner[0] -= [b['margin_x_mm'], b['margin_y_mm'], b['margin_bottom_mm']]
    inner[1] += [b['margin_x_mm'], b['margin_y_mm'], b['margin_top_mm']]
    outer = inner.copy()
    outer[0] -= [b['wall_thickness_mm'], b['wall_thickness_mm'], b['bottom_thickness_mm']]
    outer[1] += [b['wall_thickness_mm'], b['wall_thickness_mm'], 0]
    return BoxBounds(inner, outer)


def face_axis(face):
    if face not in FACES:
        raise ValueError('Unknown box boundary: ' + str(face))
    return 'XYZ'.index(face[-1]), (0 if face[0] == '-' else 1)


def wall_edge_clearance(target, face, box):
    axis, _ = face_axis(face)
    others = [k for k in range(3) if k != axis]
    return float(min(*(target[k] - box.inner[0, k] for k in others),
                     *(box.inner[1, k] - target[k] for k in others)))


def required_spacing(a, b, cfg):
    return max(cfg['ports']['min_port_spacing_mm'], a.radius + b.radius + 2 * cfg['ports']['assembly_clearance_mm'])


def shape_mesh(shape, cfg):
    vs, fs = shape.tessellate(cfg['geometry']['mesh_tolerance_mm'], cfg['geometry']['mesh_angular_tolerance_rad'])
    mesh = trimesh.Trimesh([v.toTuple() for v in vs], fs, process=True)
    if not mesh.is_volume:
        raise ValueError('CAD_TESSELLATION_NOT_CLOSED_POSITIVE_VOLUME')
    return mesh


def sample_edge(edge, cfg):
    n = max(3, math.ceil(edge.Length() / cfg['routing']['curvature_step_mm']) + 1)
    points, params = edge.sample(n)
    curvature = np.array([edge.curvatureAt(p, mode='parameter') for p in params])
    maximum = float(np.max(curvature))
    return np.array([p.toTuple() for p in points]), math.inf if maximum < 1e-12 else 1 / maximum


def natural_bend_points(endpoint, face, radius, cfg):
    """Circular guide points for an OCC interpolating spline, not a spline kernel.

    This adds a gentle single turn followed by a straight run, avoiding an
    artificial S-bend caused by forcing the wall target to the endpoint height.
    """
    t, n = endpoint.tangent, NORMALS[face]
    theta = math.acos(float(np.clip(t @ n, -1, 1)))
    if abs(math.sin(theta)) < 1e-8:
        raise ValueError('NATURAL_BEND_COLLINEAR')
    lateral = unit(n - t * math.cos(theta))
    parameters = np.linspace(0, theta, cfg['routing']['natural_bend_interpolation_points'])
    return endpoint.position + attachment_lead(endpoint, cfg) * t + radius * (np.sin(parameters)[:, None] * t + (1 - np.cos(parameters))[:, None] * lateral)


def make_route(endpoint, face, target, box, cfg, kind, tangent_scale=1.0, bend_radius=None, *, local_first_bend=None):
    """Constant radius; local overlap is additive and never edits the source STL."""
    face_axis(face)  # Explicitly reject prohibited Z boundaries at the CAD API.
    p, t, r = endpoint.position, endpoint.tangent, endpoint.radius
    n = NORMALS[face]
    overlap = cfg['routing']['attachment_overlap_mm']
    wall = boundary_thickness(face, cfg)
    stub = cfg['ports']['outside_stub_mm']
    aligned = endpoint.cap_ring is not None
    lead = attachment_lead(endpoint, cfg)
    start = p if aligned else p - overlap * t
    bend_start = p + lead * t
    if kind == 'ROUTE_A_STRAIGHT':
        if angle(t, n) > cfg['ports']['preferred_max_turn_angle_deg']:
            raise ValueError('STRAIGHT_TURN_EXCEEDS_LIMIT')
        if angle(target - p, t) > cfg['geometry']['coordinate_tolerance_mm']:
            raise ValueError('STRAIGHT_MUST_FOLLOW_TERMINAL_TANGENT')
        end = target + (wall / (t @ n) + stub) * t
        edge = cq.Edge.makeLine(cq.Vector(*start), cq.Vector(*end))
        wire = cq.Wire.assembleEdges([edge])
        interior = np.array([p, target])
        points = np.array([p, target, end])
        bend = math.inf
        direction = t
        length = float(np.linalg.norm(end - p))
    else:
        prewall = target - cfg['routing']['wall_approach_mm'] * n
        distance = float(np.linalg.norm(prewall - bend_start))
        if local_first_bend is not None:
            spline = local_first_bend
            np.testing.assert_allclose(spline.startPoint().toTuple(), bend_start, atol=1e-7)
            np.testing.assert_allclose(spline.endPoint().toTuple(), prewall, atol=1e-7)
            if angle(spline.tangentAt(0).toTuple(), t) > 1e-4 or angle(spline.tangentAt(1).toTuple(), n) > 1e-4:
                raise ValueError('LOCAL_FIRST_BEND_AXIS_MISMATCH')
        elif bend_radius is None:
            spline = cq.Edge.makeSpline([cq.Vector(*bend_start), cq.Vector(*prewall)],
                tangents=[cq.Vector(*(t * distance * tangent_scale)), cq.Vector(*(n * distance * tangent_scale))],
                parameters=[0.0, 1.0], scale=False)
        else:
            guide = natural_bend_points(endpoint, face, bend_radius, cfg)
            prewall = guide[-1]
            if (target - prewall) @ n < cfg['routing']['wall_approach_mm']:
                raise ValueError('NATURAL_BEND_CANNOT_FINISH_BEFORE_WALL')
            spline = cq.Edge.makeSpline([cq.Vector(*v) for v in guide],
                # Chord-length parameterization needs unit derivatives. OCC's
                # automatic endpoint scaling can introduce terminal wiggles.
                tangents=[cq.Vector(*t), cq.Vector(*n)], scale=False)
        interior, bend = sample_edge(spline, cfg)
        if bend < cfg['ports']['minimum_bend_radius_mm'] * cfg['routing']['bend_radius_safety_factor']:
            raise ValueError('BEND_RADIUS_TOO_SMALL')
        end = target + (wall + stub) * n
        edges = [cq.Edge.makeLine(cq.Vector(*start), cq.Vector(*bend_start)), spline,
                 cq.Edge.makeLine(cq.Vector(*prewall), cq.Vector(*end))]
        wire = cq.Wire.assembleEdges(edges)
        if aligned:
            interior = np.vstack([p, interior])
        points = np.vstack([interior, target, end])
        direction = n
        length = lead + spline.Length() + float(np.linalg.norm(target - prewall)) + wall + stub
    # Centerline plus its radius must stay within the box until the assigned wall.
    # The intended wall crossing is checked separately by the actual box Boolean.
    axis, side = face_axis(face)
    for k in range(3):
        if k == axis:
            lo = box.inner[0, k] + (r if side else 0)
            hi = box.inner[1, k] - (r if not side else 0)
        else:
            lo, hi = box.inner[0, k] + r, box.inner[1, k] - r
        if np.min(interior[:, k]) < lo or np.max(interior[:, k]) > hi:
            raise ValueError('ROUTE_LEAVES_BOX_INTERIOR')
    if aligned:
        shape = cq.Workplane().add(polygon_wire(endpoint.cap_ring)).toPending().sweep(wire, isFrenet=False).val()
        shape = shape.fuse(endpoint.attachment_collar)
    else:
        shape = cq.Workplane(cq.Plane(origin=tuple(start), normal=tuple(t))).circle(r).sweep(wire, isFrenet=True).val()
    if not shape.isValid() or len(shape.Solids()) != 1:
        raise ValueError('INVALID_PORT_SWEEP')
    return Route(endpoint, face, kind, np.asarray(target), end, direction, wire, shape, shape_mesh(shape, cfg),
                 points, float(length), angle(t, n), bend, wall_edge_clearance(target, face, box))


def surface_distance(first, second):
    """FCL triangle-surface distance, not a nearest-vertex approximation."""
    manager = trimesh.collision.CollisionManager()
    manager.add_object('surface', first)
    if manager.in_collision_single(second):
        return 0.0
    return float(manager.min_distance_single(second))


class VesselCollision:
    """Exclude only the local parent attachment collar using companion topology.

    Faces are exempt ONLY when all three vertices map to the same terminal
    branch and lie within the configured terminal arc distance. All other faces,
    including sibling vessels and distant portions of the parent, remain active.
    This is a provenance-based face partition, not centerline extraction.
    """
    def __init__(self, mesh, graph, branches, transform, endpoints, cfg):
        samples, labels, from_node = [], [], []
        step = cfg['geometry']['centerline_label_sample_step_mm']
        for a, b in graph.edges:
            pa, pb = graph.nodes[a]['coords'][:3], graph.nodes[b]['coords'][:3]
            length = float(np.linalg.norm(pb - pa))
            for fraction in np.linspace(0, 1, max(2, math.ceil(length / step) + 1)):
                native = pa * (1 - fraction) + pb * fraction
                samples.append(transform[:3, :3] @ native + transform[:3, 3])
                labels.append(branches[a, b])
                from_node.append((a, b, fraction, length))
        nearest = cKDTree(samples).query(mesh.vertices)[1]
        labels = np.array(labels)
        self.managers, self.exempt_face_counts, self.start_checks = {}, {}, {}
        undirected = graph.to_undirected()
        for a, b in undirected.edges:
            undirected[a][b]['distance'] = float(np.linalg.norm(graph.nodes[a]['coords'][:3] - graph.nodes[b]['coords'][:3]))
        for endpoint in endpoints:
            distance = nx.single_source_dijkstra_path_length(undirected, endpoint.swc_id, weight='distance')
            sample_distance = np.array([min(distance[a] + f * length, distance[b] + (1 - f) * length)
                                        for a, b, f, length in from_node])
            local = ((labels[nearest] == endpoint.branch_id) &
                     (sample_distance[nearest] <= cfg['routing']['attachment_exemption_mm']))
            exempt = np.all(local[mesh.faces], axis=1)
            obstacle = mesh.submesh([np.flatnonzero(~exempt)], append=True, repair=False)
            manager = trimesh.collision.CollisionManager()
            manager.add_object('non_attachment_vessel', obstacle)
            self.managers[endpoint.endpoint_id] = manager
            self.exempt_face_counts[endpoint.endpoint_id] = int(exempt.sum())
            # Every route must contain this same starting cross-section. Its
            # distance is a routing-independent necessary feasibility condition.
            disk = cq.Face.makeFromWires(cq.Wire.makeCircle(endpoint.radius,
                cq.Vector(*endpoint.position), cq.Vector(*endpoint.tangent)))
            vertices, faces = disk.tessellate(cfg['geometry']['mesh_tolerance_mm'],
                                             cfg['geometry']['mesh_angular_tolerance_rad'])
            disk_mesh = trimesh.Trimesh([v.toTuple() for v in vertices], faces, process=True)
            gap, detail = manager.min_distance_single(disk_mesh, return_data=True)
            _, center_distance, _ = trimesh.proximity.closest_point(obstacle, endpoint.position[None, :])
            vessel_point = detail.point('non_attachment_vessel')
            nearest_sample = cKDTree(samples).query(vessel_point)[1]
            self.start_checks[endpoint.endpoint_id] = dict(
                distance_mm=float(gap), center_distance_mm=float(center_distance[0]), vessel_point_mm=vessel_point.tolist(),
                endpoint_disk_point_mm=detail.point('__external').tolist(),
                nearest_vessel_branch=str(labels[nearest_sample]), endpoint_branch=endpoint.branch_id,
                required_clearance_mm=cfg['ports']['route_clearance_mm'],
                passed=float(gap) >= cfg['ports']['route_clearance_mm'] + 2 * cfg['geometry']['mesh_tolerance_mm'],
                meaning='Necessary condition: each constant-radius route contains the starting disk; a failing gap cannot be fixed by choosing another wall')

    def distance(self, route):
        manager = self.managers[route.endpoint.endpoint_id]
        if manager.in_collision_single(route.mesh):
            return 0.0
        return float(manager.min_distance_single(route.mesh))


def candidate_score(distance, turn, collision_penalty, edge_penalty, box, cfg):
    values = dict(normalized_distance=distance / np.linalg.norm(box.inner_size),
                  normalized_turn_angle=turn / 180.0, collision_penalty=collision_penalty,
                  wall_edge_penalty=edge_penalty)
    return float(sum(values[key] * weight for key, weight in cfg['routing']['score_weights'].items()))


def generate_options(endpoint, box, cfg, collision=None, *, clearance_override=None, stage='STRICT'):
    """Finite, deterministic and auditable target/tangent enumeration."""
    options, rows = [], []
    p, t = endpoint.position, endpoint.tangent
    clearance = cfg['ports']['min_wall_edge_clearance_mm']
    tolerance = cfg['geometry']['mesh_tolerance_mm']
    required = cfg['ports']['route_clearance_mm'] if clearance_override is None else float(clearance_override)
    for face in enabled_faces(cfg):
        axis, side = face_axis(face)
        u, v = [k for k in range(3) if k != axis]
        nominal = p.copy()
        nominal[axis] = box.inner[side, axis]
        targets = []
        if t @ NORMALS[face] > 0 and angle(t, NORMALS[face]) <= cfg['ports']['preferred_max_turn_angle_deg']:
            target = p + t * ((nominal[axis] - p[axis]) / t[axis])
            targets.append(('ROUTE_A_STRAIGHT', target, 1.0, None))
        for radius in cfg['routing']['natural_bend_radii_mm']:
            try:
                guide = natural_bend_points(endpoint, face, radius, cfg)
                target = guide[-1].copy()
                target[axis] = box.inner[side, axis]
                targets.append(('ROUTE_B_CURVED', target, 1.0, radius))
            except ValueError:
                pass
        unique = set()
        offsets_u = cfg['routing']['target_lateral_offsets_mm']
        offsets_v = cfg['routing']['target_vertical_offsets_mm']
        for lateral, vertical in itertools.product(offsets_u, offsets_v):
            target = nominal.copy()
            target[u] = np.clip(p[u] + lateral, box.inner[0, u] + clearance, box.inner[1, u] - clearance)
            target[v] = np.clip(p[v] + vertical, box.inner[0, v] + clearance, box.inner[1, v] - clearance)
            identity = tuple(target)
            if identity in unique:
                continue
            unique.add(identity)
            targets.extend(('ROUTE_B_CURVED', target, scale, None) for scale in cfg['routing']['spline_tangent_scales'])
        face_options = []
        for index, (kind, target, scale, radius) in enumerate(targets):
            option_id = f'{stage}_{endpoint.endpoint_id}_{face}_{index:03d}'
            edge_clear = wall_edge_clearance(target, face, box)
            turn = angle(t, NORMALS[face])
            row = dict(option_id=option_id, endpoint_id=endpoint.endpoint_id, face=face, route_type=kind,
                       search_stage=stage, required_vessel_clearance_mm=required,
                       boundary_kind=boundary_kind(face), boundary_thickness_mm=boundary_thickness(face, cfg),
                       target_x=float(target[0]), target_y=float(target[1]), target_z=float(target[2]), tangent_scale=scale,
                       guide_bend_radius_mm=radius, spline_family='natural_turn' if radius else 'two_point',
                       distance_to_face_mm=abs(float(nominal[axis] - p[axis])), turn_angle_deg=turn,
                       distance_to_target_mm=float(np.linalg.norm(target - p)),
                       normalized_distance=float(np.linalg.norm(target - p) / np.linalg.norm(box.inner_size)),
                       normalized_turn_angle=turn / 180.0,
                       wall_edge_clearance_mm=edge_clear, vessel_clearance_mm=None, minimum_bend_radius_mm=None,
                       collision_penalty=1.0, wall_edge_penalty=float(edge_clear < clearance),
                       score=None, status='', selected=False, nearest_other_port_mm=None, same_wall_spacing_mm=None)
            route = None
            try:
                if collision and collision.start_checks[endpoint.endpoint_id]['distance_mm'] < required + 2 * tolerance:
                    row['vessel_clearance_mm'] = collision.start_checks[endpoint.endpoint_id]['distance_mm']
                    raise ValueError('ENDPOINT_CLEARANCE_CONFLICT')
                if edge_clear < clearance - cfg['geometry']['coordinate_tolerance_mm']:
                    raise ValueError('WALL_EDGE_CLEARANCE')
                route = make_route(endpoint, face, target, box, cfg, kind, scale, radius)
                route.required_vessel_clearance = required
                route.preferred_vessel_clearance = cfg['ports']['route_clearance_mm']
                route.vessel_clearance = collision.distance(route) if collision else math.inf
                row['vessel_clearance_mm'] = route.vessel_clearance
                row['minimum_bend_radius_mm'] = route.minimum_bend_radius
                if route.vessel_clearance < required + 2 * tolerance:
                    raise ValueError('PORT_COLLISION_RISK')
                row['status'] = 'INDIVIDUALLY_FEASIBLE'
                row['collision_penalty'] = 0.0
            except (ValueError, RuntimeError) as exc:
                row['status'] = str(exc)
            row['score'] = candidate_score(np.linalg.norm(target - p), turn, row['collision_penalty'], row['wall_edge_penalty'], box, cfg)
            if row['status'] == 'INDIVIDUALLY_FEASIBLE':
                route.score, route.option_id = row['score'], option_id
                face_options.append(route)
            rows.append(row)
        # Preserve route diversity, not just many spline scales at one target.
        best_by_target = {}
        for route in sorted(face_options, key=lambda x: (x.kind != 'ROUTE_A_STRAIGHT', x.score, x.length, x.option_id)):
            best_by_target.setdefault(tuple(route.target), route)
        options.extend(list(best_by_target.values())[:cfg['routing']['options_per_face']])
    return options, rows


def allowed_layout(routes, cfg):
    outlet_faces = {r.face for r in routes if r.endpoint.role == 'outlet'}
    if len(outlet_faces) > cfg['routing']['max_outlet_faces']:
        return False
    if len(outlet_faces) == 2:
        a, b = sorted(outlet_faces)
        if float(NORMALS[a] @ NORMALS[b]) != 0:
            return False
    return True


def select_layout(option_lists, cfg):
    """Choose an opposite-wall layout first when feasible, then adjacent fallback.

    Exhaustive small-tree search with cached library distances and lower score
    bounds. Each phase returns the least weighted-score feasible combination.
    """
    if not option_lists or any(not options for options in option_lists):
        return [], {'status': 'PORT_LAYOUT_NEEDS_MANUAL_REVIEW', 'reason': 'An endpoint has no feasible route'}
    pair_cache = {}
    tolerance = 2 * cfg['geometry']['mesh_tolerance_mm']

    def compatible(a, b):
        if a.face == b.face and np.linalg.norm(a.target - b.target) < required_spacing(a.endpoint, b.endpoint, cfg):
            return False
        key = tuple(sorted((a.option_id, b.option_id)))
        if key not in pair_cache:
            pair_cache[key] = surface_distance(a.mesh, b.mesh)
        return pair_cache[key] >= cfg['ports']['route_clearance_mm'] + tolerance

    phases = ('opposite', 'adjacent_fallback') if cfg['routing']['prefer_opposite_faces'] else ('unpreferred',)
    ordered = [sorted(options, key=lambda r: (r.score, r.option_id)) for options in option_lists]
    suffix = np.cumsum([min(r.score for r in options) for options in ordered][::-1])[::-1]
    for phase in phases:
        best, best_score = [], math.inf

        def visit(chosen, score):
            nonlocal best, best_score
            index = len(chosen)
            if index == len(ordered):
                if score < best_score:
                    best, best_score = chosen[:], score
                return
            if score + suffix[index] >= best_score:
                return
            inlet = next((r for r in chosen if r.endpoint.role == 'inlet'), None)
            for route in ordered[index]:
                if score + route.score >= best_score:
                    continue
                if phase == 'opposite' and inlet and route.endpoint.role == 'outlet' and not np.isclose(NORMALS[inlet.face] @ NORMALS[route.face], -1):
                    continue
                if phase == 'opposite' and route.kind == 'ROUTE_B_CURVED' and route.turn > cfg['routing']['preferred_max_curved_turn_angle_deg']:
                    continue
                if not allowed_layout(chosen + [route], cfg) or any(not compatible(route, old) for old in chosen):
                    continue
                visit(chosen + [route], score + route.score)
        visit([], 0.0)
        if best:
            for a in best:
                others = [b for b in best if b is not a]
                a.nearest_other_port = min(surface_distance(a.mesh, b.mesh) for b in others) if others else math.inf
                a.wall_spacing = min((np.linalg.norm(a.target - b.target) for b in others if a.face == b.face), default=math.inf)
            return best, dict(status='LAYOUT_GEOMETRY_PASS', preference_phase=phase, total_score=best_score,
                              evaluated_port_pairs=len(pair_cache))
    return [], dict(status='PORT_LAYOUT_NEEDS_MANUAL_REVIEW', reason='No combination meets spacing, collision and face constraints',
                    evaluated_port_pairs=len(pair_cache))


def inspect_native_continuation(graph, endpoint, retained_ids):
    """Trace real directed source edges, never spatially nearby vessels.

    The length is dictated by the available source topology, not a fixed graft
    length. A bifurcation is reported with every daughter; none is silently
    chosen or deleted. This audit does not extrapolate or construct anatomy.
    """
    node = endpoint.original_swc_id
    if node not in graph:
        raise ValueError(f'Original endpoint ID {node} is missing')
    ids, length = [node], 0.0
    while True:
        neighbors = sorted(graph.successors(node) if endpoint.role == 'outlet' else graph.predecessors(node))
        if not neighbors:
            reason = 'ORIGINAL_TERMINAL_NO_CONTINUATION' if len(ids) == 1 else 'ORIGINAL_CONTINUATION_AVAILABLE'
            break
        if len(neighbors) > 1:
            reason = 'ORIGINAL_BIFURCATION_AVAILABLE'
            break
        nxt = neighbors[0]
        if nxt in retained_ids:
            reason = 'CONTINUATION_ALREADY_IN_FROZEN_CORE'
            break
        length += float(np.linalg.norm(graph.nodes[nxt]['coords'][:3] - graph.nodes[node]['coords'][:3]))
        ids.append(nxt)
        node = nxt
    return dict(endpoint_id=endpoint.endpoint_id, original_swc_id=endpoint.original_swc_id,
                raw_parent_ids=sorted(graph.predecessors(endpoint.original_swc_id)),
                raw_child_ids=sorted(graph.successors(endpoint.original_swc_id)),
                status=reason, available_outward_node_ids=ids[1:], available_arc_length_mm=length,
                stopping_node_id=node, stopping_outward_neighbors=neighbors,
                restored_node_count=0, restored_length_mm=0.0,
                adaptation_basis='Only source parent-child edges beyond the current endpoint; stop at a real terminal or branching event')


def audit_native_restoration(inputs, cfg, start_checks):
    path = (ROOT / cfg['native_restoration']['source_swc']).resolve()
    from .project_paths import frozen_hash
    expected = frozen_hash(inputs['manifest']['source_hashes'], path)
    actual = sha256(path)
    if expected is None or expected != actual:
        raise ValueError('Original SWC hash does not match the frozen compact provenance')
    graph = read_source(path, allow_forest=True).graph
    compact = read_source(inputs['paths']['compact_swc']).graph
    mapping = read_csv(inputs['paths']['node_mapping'])
    retained = {int(row['original_swc_id']) for row in mapping}
    # Check ALL retained coordinates and directed edges, not just a nearby ID.
    ids = {int(row['swc_id']): int(row['original_swc_id']) for row in mapping}
    deviations = []
    for n, original in ids.items():
        deviations.append(float(np.linalg.norm(compact.nodes[n]['coords'][:3] - graph.nodes[original]['coords'][:3])))
    if max(deviations) > cfg['geometry']['coordinate_tolerance_mm'] or any(not graph.has_edge(ids[a], ids[b]) for a, b in compact.edges):
        raise ValueError('Original topology/coordinate correspondence is not verified')
    rows = []
    for endpoint in inputs['endpoints']:
        row = inspect_native_continuation(graph, endpoint, retained)
        row['restoration_needed'] = not start_checks[endpoint.endpoint_id]['passed']
        raw = graph.nodes[endpoint.original_swc_id]['coords'][:3]
        raw_print = inputs['transform'][:3, :3] @ raw + inputs['transform'][:3, 3]
        row['original_position_print_mm'] = raw_print.tolist()
        row['fitted_endpoint_displacement_mm'] = float(np.linalg.norm(raw_print - endpoint.position))
        rows.append(row)
    return dict(source_swc=str(path), source_sha256=actual, hash_matches_frozen_provenance=True,
                all_retained_directed_edges_verified=True, maximum_retained_coordinate_difference_mm=max(deviations),
                priority='STRICT_ROUTING_THEN_ORIGINAL_RESTORATION_THEN_LAST_RESORT_CLEARANCE',
                endpoints=rows, restored_total_length_mm=0.0,
                note='No downstream vessels are invented. A true original leaf cannot be extended by restoring absent source data.')


def try_clearance_fallback(endpoints, option_lists, box, cfg, collision, restoration):
    """After original-data audit, relax only endpoints without strict options.

    All CAD candidates are generated once, then tested at decreasing clearance
    thresholds. Port-to-port spacing, wall margins and actual nonintersection
    remain hard gates. A nonzero mesh safety allowance remains even at 0 mm.
    """
    audit = dict(used=False, attempts=[], affected_endpoints=[], authorization=cfg.get('clearance_fallback', {}).get('authorization'))
    if not cfg.get('clearance_fallback', {}).get('enabled'):
        return [], audit, []
    blocked = [i for i, options in enumerate(option_lists) if not options]
    if not blocked:
        audit['reason'] = 'No individually blocked endpoint; global layout constraints are retained'
        return [], audit, []
    evidence = {row['endpoint_id']: row for row in restoration['endpoints']}
    if not restoration.get('hash_matches_frozen_provenance') or not restoration.get('all_retained_directed_edges_verified'):
        raise ValueError('Last resort requires verified original anatomy')
    for index in blocked:
        row = evidence[endpoints[index].endpoint_id]
        if row['status'] != 'ORIGINAL_TERMINAL_NO_CONTINUATION':
            # Never silently bypass genuine available anatomy with a relaxed
            # artificial stem. This real dataset has no such continuation.
            audit['reason'] = 'Original vessel continuation exists and must be restored before relaxation: ' + row['endpoint_id']
            return [], audit, []
    audit['affected_endpoints'] = [endpoints[i].endpoint_id for i in blocked]
    strict = cfg['ports']['route_clearance_mm']
    lower = cfg['clearance_fallback']['minimum_mm']
    delta = cfg['clearance_fallback']['decrement_mm']
    clearance_levels = sorted({round(max(lower, strict - delta * n), 10)
        for n in range(1, math.ceil((strict - lower) / delta) + 1)}, reverse=True)
    tolerance = 2 * cfg['geometry']['mesh_tolerance_mm']
    rows, relaxed = [], {}
    for index in blocked:
        relaxed[index], records = generate_options(endpoints[index], box, cfg, collision,
            clearance_override=lower, stage='LAST_RESORT_SEARCH')
        rows.extend(records)
    for level in clearance_levels:
        pools = list(option_lists)
        for index in blocked:
            pools[index] = [route for route in relaxed[index] if route.vessel_clearance >= level + tolerance]
        routes, layout = select_layout(pools, cfg)
        audit['attempts'].append(dict(required_vessel_clearance_mm=level,
            feasible_option_counts={endpoints[i].endpoint_id: len(pools[i]) for i in blocked}, layout=layout))
        if not routes:
            continue
        for route in routes:
            if route.endpoint.endpoint_id in audit['affected_endpoints']:
                route.required_vessel_clearance = level
        audit.update(used=True, selected_clearance_mm=level, strict_clearance_mm=strict,
                     minimum_actual_vessel_clearance_mm=min(r.vessel_clearance for r in routes),
                     port_to_port_clearance_mm=strict, layout=layout,
                     note='Largest feasible threshold on the configured descending grid; only blocked endpoints relaxed')
        return routes, audit, rows
    audit['reason'] = 'No nonintersecting layout at any allowed fallback threshold'
    return [], audit, rows


def clearance_profile(routes, collision, cfg):
    """Diagnostic cross-section samples; acceptance uses the complete mesh."""
    rows = []
    step = cfg['routing']['curvature_step_mm']
    for route in routes:
        lengths = np.r_[0, np.cumsum(np.linalg.norm(np.diff(route.points, axis=0), axis=1))]
        distances = np.linspace(0, lengths[-1], max(2, math.ceil(lengths[-1] / step) + 1))
        points = np.column_stack([np.interp(distances, lengths, route.points[:, k]) for k in range(3)])
        tangents = np.gradient(points, axis=0)
        tangents[0] = route.endpoint.tangent
        manager = collision.managers[route.endpoint.endpoint_id]
        for s, p, t in zip(distances, points, tangents):
            disk = cq.Face.makeFromWires(cq.Wire.makeCircle(route.endpoint.radius, cq.Vector(*p), cq.Vector(*unit(t))))
            vertices, faces = disk.tessellate(cfg['geometry']['mesh_tolerance_mm'], cfg['geometry']['mesh_angular_tolerance_rad'])
            mesh = trimesh.Trimesh([v.toTuple() for v in vertices], faces, process=True)
            gap = 0.0 if manager.in_collision_single(mesh) else float(manager.min_distance_single(mesh))
            rows.append(dict(endpoint_id=route.endpoint.endpoint_id, arc_from_original_endpoint_mm=float(s),
                x=p[0], y=p[1], z=p[2], cross_section_gap_mm=gap,
                below_preferred_clearance=gap < cfg['ports']['route_clearance_mm']))
    return rows


def create_box(box, routes, cfg):
    outer = cq.Solid.makeBox(*box.outer_size, pnt=cq.Vector(*box.outer[0]))
    inner_size = box.inner_size.copy()
    inner_size[2] += cfg['box']['wall_thickness_mm']  # cavity cutter projects above the rim
    cavity = cq.Solid.makeBox(*inner_size, pnt=cq.Vector(*box.inner[0]))
    body = outer.cut(cavity)
    hole_checks = []
    for route in routes:
        radius = route.endpoint.radius + cfg['ports']['assembly_clearance_mm']
        n = NORMALS[route.face]
        along = boundary_thickness(route.face, cfg) / (route.direction @ n)
        extra = cfg['routing']['wall_approach_mm']
        cutter = cq.Solid.makeCylinder(radius, along + 2 * extra,
            pnt=cq.Vector(*(route.target - extra * route.direction)), dir=cq.Vector(*route.direction))
        before = body.Volume()
        body = body.cut(cutter)
        remainder = body.intersect(cutter).Volume()
        removed = before - body.Volume()
        hole_checks.append(dict(port_id=route.endpoint.endpoint_id, radius_mm=radius, removed_volume_mm3=removed,
                                remaining_cutter_intersection_mm3=remainder,
                                passed=removed > cfg['geometry']['volume_tolerance_mm3'] and remainder < cfg['geometry']['volume_tolerance_mm3']))
    if not body.isValid() or len(body.Solids()) != 1:
        raise ValueError('BOX_BOOLEAN_FAILED')
    return body, hole_checks


def mesh_qc(mesh, cfg):
    edges = np.bincount(mesh.edges_unique_inverse, minlength=len(mesh.edges_unique))
    values = dict(watertight=bool(mesh.is_watertight), manifold=bool(np.all(edges == 2)),
                  connected_components=len(mesh.split(only_watertight=False)), positive_volume=bool(mesh.volume > 0),
                  volume_mm3=float(mesh.volume), finite_vertices=bool(np.isfinite(mesh.vertices).all()),
                  degenerate_faces=int(np.sum(mesh.area_faces <= cfg['geometry']['minimum_triangle_area_mm2'])),
                  winding_consistent=bool(mesh.is_winding_consistent), vertices=len(mesh.vertices), faces=len(mesh.faces))
    values['passed'] = all(values[key] for key in ('watertight', 'manifold', 'positive_volume', 'finite_vertices', 'winding_consistent')) and values['connected_components'] == 1 and values['degenerate_faces'] == 0
    return values


def box_dimension_checks(body, box, cfg):
    """Measure material thickness with OCC solid/line intersections."""
    measured = {}
    offset = cfg['box']['wall_thickness_mm']
    for face in SIDE_FACES:
        axis, side = face_axis(face)
        start = box.inner.mean(axis=0)
        start[2] = box.inner[1, 2] - offset / 2  # safely above all hole centers
        end = start.copy()
        sign = -1 if side == 0 else 1
        start[axis] = box.outer[side, axis] + sign * offset
        end[axis] = box.inner[side, axis] - sign * offset
        section = body.intersect(cq.Edge.makeLine(cq.Vector(*start), cq.Vector(*end)))
        measured[face] = float(sum(edge.Length() for edge in section.Edges()))
    # Probe the floor near the inner corner, independently of vessel position.
    start = box.inner.mean(axis=0)
    start[:2] = box.inner[0, :2] + offset / 2
    end = start.copy()
    start[2], end[2] = box.outer[0, 2] - offset, box.inner[0, 2] + offset
    section = body.intersect(cq.Edge.makeLine(cq.Vector(*start), cq.Vector(*end)))
    bottom = float(sum(edge.Length() for edge in section.Edges()))
    # CadQuery's default AddOptimal uses cached triangulation. After tessellation
    # its deflection padding can inflate a dimension by several 1e-4 mm. Measure
    # the underlying CAD geometry, independently of the display mesh cache.
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    bounds = Bnd_Box()
    BRepBndLib.AddOptimal_s(body.wrapped, bounds, False, False)
    limits = np.asarray(bounds.Get())
    measured_size = limits[3:] - limits[:3]
    tol = cfg['geometry']['coordinate_tolerance_mm']
    passed = all(abs(v - cfg['box']['wall_thickness_mm']) < tol for v in measured.values())
    passed &= abs(bottom - cfg['box']['bottom_thickness_mm']) < tol
    passed &= bool(np.max(np.abs(measured_size - box.outer_size)) < tol)
    return dict(side_wall_thickness_mm=measured, bottom_thickness_mm=bottom,
                outer_dimensions_mm=measured_size.tolist(), passed=bool(passed),
                method='OCC solid/line intersections away from openings')


def export_checked_stl(mesh, path, cfg):
    mesh.export(path)
    exported = trimesh.load_mesh(path)
    qc = mesh_qc(exported, cfg)
    if not qc['passed']:
        Path(path).unlink()
        raise ValueError('EXPORTED_STL_QC_FAILED: ' + json.dumps(qc))
    return qc


def union_core(mesh, routes, cfg):
    result = trimesh.boolean.union([mesh] + [r.mesh for r in routes], engine='manifold', check_volume=True)
    qc = mesh_qc(result, cfg)
    qc['boolean_cleanup_used'] = False
    # A float32 STL Boolean can leave very thin triangles at coincident sweep
    # seams. Use the official kernel simplifier ONLY for this numerical case.
    # Never repair disconnected/nonmanifold geometry or the frozen input mesh.
    if qc['degenerate_faces'] and qc['watertight'] and qc['manifold'] and qc['connected_components'] == 1:
        import manifold3d as manifold
        solid = manifold.Manifold(manifold.Mesh64(np.asarray(result.vertices, dtype=np.float64),
                                                   np.asarray(result.faces, dtype=np.uint64)))
        solid = solid.simplify(cfg['geometry']['boolean_cleanup_tolerance_mm'])
        simplified = solid.to_mesh()
        result = trimesh.Trimesh(simplified.vert_properties, simplified.tri_verts, process=False)
        qc = mesh_qc(result, cfg)
        qc.update(boolean_cleanup_used=True, boolean_cleanup_tolerance_mm=cfg['geometry']['boolean_cleanup_tolerance_mm'])
    if not qc['passed']:
        raise ValueError('CORE_PORT_UNION_FAILED: ' + json.dumps(qc))
    missing = trimesh.boolean.difference([mesh, result], engine='manifold', check_volume=True)
    qc['source_subtracted_volume_mm3'] = float(abs(missing.volume)) if len(missing.faces) else 0.0
    qc['source_solid_preserved'] = qc['source_subtracted_volume_mm3'] <= cfg['geometry']['volume_tolerance_mm3']
    if not qc['source_solid_preserved']:
        raise ValueError('CORE_PORT_UNION_FAILED: original vascular volume was removed')
    return result, qc


def protected_snapshot(inputs, out):
    base = inputs['base']
    files = set(inputs['paths'].values())
    # Include current hashes, not stale pre-UI-development hashes from older runs.
    prior = base / 'protected_sources.json'
    if prior.exists():
        files.update(Path(p) for p in json.loads(prior.read_text()) if Path(p).is_file())
    files.update(p for p in base.rglob('*') if p.is_file() and not p.is_relative_to(out))
    files.update(ROOT / name for name in ('s1-2_swc_roi_generate_human.py', 's1-3_swc_roi_generate_MeVO.py'))
    return {str(p.resolve()): sha256(p) for p in sorted(files)}


def verify_snapshot(snapshot):
    current = {p: sha256(Path(p)) if Path(p).is_file() else None for p in snapshot}
    changed = [p for p in snapshot if current[p] != snapshot[p]]
    return dict(all_unchanged=not changed, protected_file_count=len(snapshot), changed=changed,
                before_sha256=snapshot, after_sha256=current)
