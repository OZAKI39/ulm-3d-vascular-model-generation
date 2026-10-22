"""Shared known-cap attachment; extracted unchanged from the verified O3 method.

All endpoints use this one implementation. The only collar fallback shortens the
buried frustum; neither its radius fraction nor the containment tolerance changes.
"""
import cadquery as cq
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
import trimesh


def collar_outside_volume(collar, source, cfg):
    from .sacrificial_fixture import shape_mesh
    outside = trimesh.boolean.difference([shape_mesh(collar, cfg), source], engine='manifold')
    return abs(float(outside.volume)) if len(outside.faces) else 0.0


def align_known_cap(endpoint, source, cap_record, cfg):
    """Match an already identified SWC/manifest cap; never discover endpoints.

    The frozen print exporter creates a fan with a known center and vertex
    count. Its boundary, rather than the 5 mm mean centerline direction, defines
    the attachment. Float32 STL noise is projected onto the measured cap plane.
    """
    from .sacrificial_fixture import unit, angle, polygon_wire
    settings = cfg['attachment_alignment']
    tolerance = settings['cap_match_tolerance_mm']
    centers = cKDTree(source.vertices).query_ball_point(endpoint.position, tolerance)
    gap, center_id = cKDTree(source.vertices).query(endpoint.position)
    if len(centers) != 1 or gap > tolerance:
        raise ValueError('KNOWN_CAP_CENTER_NOT_FOUND')
    face_ids = np.flatnonzero(np.any(source.faces == center_id, axis=1))
    boundary = nx.Graph()
    for triangle in source.faces[face_ids]:
        rim = [int(n) for n in triangle if n != center_id]
        boundary.add_edge(*rim)
    if len(face_ids) != cap_record['vertices'] or len(boundary) != cap_record['vertices'] or not nx.is_connected(boundary) or any(boundary.degree(n) != 2 for n in boundary):
        raise ValueError('KNOWN_CAP_FAN_TOPOLOGY_MISMATCH')
    order = [min(boundary)]
    while len(order) < len(boundary):
        order.append(next(n for n in sorted(boundary.neighbors(order[-1])) if n not in order))
    ring = source.vertices[order].copy()
    center = ring.mean(axis=0)
    normal = unit(np.sum(source.face_normals[face_ids] * source.area_faces[face_ids, None], axis=0))
    # Normally retain the verified O3 ring-centroid plane verbatim. A frozen
    # float32 cap fan can have its known center outside the rim's axial slab.
    # In that case its centroid plane is not a plane through the known cap:
    # anchor the SAME measured normal at the provenance-matched center vertex.
    # This is a bounded plane projection, performed BEFORE collar feasibility,
    # never a containment tolerance change or a fitted/replacement circle.
    rim_offsets = (ring - center) @ normal
    fan_offset = float((source.vertices[center_id] - center) @ normal)
    plane_anchor = 'RING_CENTROID'
    if fan_offset < rim_offsets.min() or fan_offset > rim_offsets.max():
        center = center + fan_offset * normal
        plane_anchor = 'KNOWN_CAP_CENTER_VERTEX_AXIAL_PLANE'
    if np.linalg.norm(center - endpoint.position) > tolerance:
        raise ValueError('KNOWN_CAP_RING_CENTER_MISMATCH')
    plane_errors = (ring - center) @ normal
    radial = np.linalg.norm(ring - center, axis=1)
    if max(abs(plane_errors)) > tolerance or max(abs(radial - endpoint.radius)) > tolerance:
        raise ValueError('KNOWN_CAP_NOT_PLANAR_CIRCULAR_WITH_EXPECTED_RADIUS')
    if min(source.face_normals[face_ids] @ normal) < 1 - 1e-5 or normal @ endpoint.tangent <= 0:
        raise ValueError('KNOWN_CAP_NORMAL_OR_ORIENTATION_MISMATCH')
    ring -= plane_errors[:, None] * normal
    if np.cross(ring[1] - center, ring[2] - center) @ normal < 0:
        ring = ring[::-1]
    qc = dict(method='SWC_IDENTIFIED_FROZEN_CAP_RING_AND_OUTWARD_NORMAL',
        plane_anchor=plane_anchor, known_fan_center_axial_offset_mm=fan_offset,
        reference_swc_mean_tangent=endpoint.tangent.tolist(), cap_normal=normal.tolist(),
        before_axis_error_deg=angle(endpoint.tangent, normal), after_axis_error_deg=0.0,
        known_center_vertex_distance_mm=float(gap), center_adjustment_mm=float(np.linalg.norm(center - endpoint.position)),
        cap_vertex_count=len(ring), original_ring_print_mm=source.vertices[order].tolist(),
        aligned_ring_print_mm=ring.tolist(), maximum_stl_plane_projection_mm=float(max(abs(plane_errors))),
        maximum_ring_radius_error_mm=float(max(abs(radial - endpoint.radius))),
        radius_mm=endpoint.radius, radius_changed=False, frozen_source_modified=False)
    # The overlap is tapered and buried INSIDE the original lumen solid. A
    # backward full-radius tube can protrude through a curved parent sidewall.
    attempts = []
    lengths = settings.get('collar_length_fallback_mm', [settings['collar_length_mm']])
    if lengths != sorted(set(lengths), reverse=True) or lengths[0] != settings['collar_length_mm']:
        raise ValueError('INVALID_COLLAR_FALLBACK_ORDER')
    for length in lengths:
        back = center - length * normal
        inset = back + settings['collar_inner_radius_fraction'] * (ring - center)
        collar = cq.Solid.makeLoft([polygon_wire(inset), polygon_wire(ring)], ruled=True)
        outside_volume = collar_outside_volume(collar, source, cfg)
        contained = outside_volume <= settings['collar_outside_volume_tolerance_mm3']
        attempts.append(dict(length_mm=length, outside_volume_mm3=outside_volume, contained=contained))
        if contained:
            break
    else:
        raise ValueError('CAP_COLLAR_NOT_FEASIBLE: ATTACHMENT_COLLAR_PROTRUDES_' + endpoint.endpoint_id
                         + ': ' + repr(attempts))
    qc.update(collar_length_mm=length, collar_inner_radius_fraction=settings['collar_inner_radius_fraction'],
        collar_outside_source_volume_mm3=outside_volume, collar_contained=True, collar_attempts=attempts,
        maximum_plane_projection_mm=qc['maximum_stl_plane_projection_mm'],
        provenance=dict(endpoint_id=endpoint.endpoint_id, fitted_swc_id=endpoint.swc_id,
            original_swc_id=endpoint.original_swc_id, branch_id=endpoint.branch_id,
            stl_center_vertex_id=int(center_id), stl_fan_face_ids=face_ids.tolist(),
            stl_ring_vertex_ids=order, manifest_cap=cap_record, unique=True),
        production_profile='ACTUAL_CAP_POLYGON', full_radius_backward_overlap_used=False)
    endpoint.position, endpoint.tangent, endpoint.cap_ring = center, normal, ring
    endpoint.attachment_collar, endpoint.attachment_qc = collar, qc
    return qc



def fixed_target_first_bend(endpoint, face, target, cfg):
    """Refit only the curved span, with both endpoint positions/axes fixed.

    The old mean-axis circular guide cannot reach the unchanged target after
    changing the start axis. A single deterministic local quintic solve adjusts
    two endpoint handle lengths and two interior controls. There is no face,
    target, radius, or global route enumeration. OCC builds the actual edge and
    the existing route code independently measures curvature and box clearance.
    """
    from scipy.optimize import minimize
    from scipy.special import comb
    from .sacrificial_fixture import NORMALS, attachment_lead
    normal = NORMALS[face]
    start = endpoint.position + attachment_lead(endpoint, cfg) * endpoint.tangent
    end = np.asarray(target) - cfg['routing']['wall_approach_mm'] * normal
    distance = float(np.linalg.norm(end - start))
    u = np.linspace(0, 1, 501)[:, None]
    def basis(degree):
        return np.hstack([comb(degree, j) * u**j * (1-u)**(degree-j) for j in range(degree+1)])
    derivative_basis, second_basis = basis(4), basis(3)
    def controls(v):
        return np.vstack([start, start+v[0]*endpoint.tangent, v[2:5], v[5:8], end-v[1]*normal, end])
    def curvature(v):
        points = controls(v)
        d = derivative_basis @ (5*np.diff(points, axis=0))
        dd = second_basis @ (20*np.diff(points, n=2, axis=0))
        return np.linalg.norm(np.cross(d, dd), axis=1) / np.maximum(np.linalg.norm(d, axis=1), 1e-12)**3
    initial = np.r_[distance/5, distance/5, start+2*distance/5*endpoint.tangent,
                    end-2*distance/5*normal, .2]
    result = minimize(lambda v: v[-1]+1e-7*np.sum((v[:-1]-initial[:-1])**2), initial,
        method='SLSQP', constraints=[dict(type='ineq', fun=lambda v: v[-1]-curvature(v))],
        bounds=[(.1, distance)]*2+[(None, None)]*6+[(.01, 1.)],
        options=dict(maxiter=600, ftol=1e-10))
    if not result.success:
        raise ValueError('LOCAL_FIRST_BEND_NOT_FEASIBLE_' + endpoint.endpoint_id + ': ' + result.message)
    edge = cq.Edge.makeBezier([cq.Vector(*point) for point in controls(result.x)])
    return edge, dict(method='FIXED_TARGET_LOCAL_QUINTIC_BEND', optimizer_success=bool(result.success),
        control_points_mm=controls(result.x).tolist(), estimated_minimum_bend_radius_mm=float(1/curvature(result.x).max()),
        start_mm=start.tolist(), end_mm=end.tolist(), wall_target_mm=np.asarray(target).tolist(),
        wall_target_changed=False, wall_assignment_changed=False, radius_changed=False,
        start_axis=endpoint.tangent.tolist(), end_axis=normal.tolist(), global_route_search=False)


def clean_boolean_attachment_slivers(mesh, source, endpoints, cfg):
    """Weld only sub-tolerance duplicate vertices of the NEW Boolean seam.

    Float32 Boolean intersections can leave folded triangles a few micrometres
    wide. Whole-mesh simplification need not remove them. This bounded weld uses
    the existing numerical cleanup tolerance, retains a frozen-source vertex
    when present, never averages coordinates, and never touches the input mesh.
    No holes are filled and no surface is remeshed. The result must independently
    pass closedness, volume preservation and the same feature-edge audit.
    """
    from .surface_continuity_qc import window
    from .sacrificial_fixture import mesh_qc
    tolerance=cfg['geometry']['boolean_cleanup_tolerance_mm']
    local=np.zeros(len(mesh.vertices),dtype=bool)
    for endpoint in endpoints:
        local |= window(mesh.vertices,endpoint.position,endpoint.tangent,endpoint.radius*1.7,1.)
    ids=np.flatnonzero(local)
    pairs=cKDTree(mesh.vertices[ids]).query_pairs(tolerance,output_type='ndarray')
    graph=nx.Graph();graph.add_edges_from(ids[pairs])
    remap=np.arange(len(mesh.vertices))
    source_distances=cKDTree(source.vertices).query(mesh.vertices)[0]
    audit=[]
    for component in nx.connected_components(graph):
        component=sorted(component)
        representative=min(component,key=lambda j:(source_distances[j],j))
        displacement=np.linalg.norm(mesh.vertices[component]-mesh.vertices[representative],axis=1)
        # Do not allow single-link chains to grow beyond the fixed tolerance,
        # or collapse two different vertices belonging to the frozen source.
        if displacement.max()>tolerance or np.count_nonzero(source_distances[component]<1e-12)>1:
            raise ValueError('NUMERICAL_SEAM_WELD_EXCEEDS_SCOPE')
        remap[component]=representative
        audit.append(dict(vertex_ids=component,representative_vertex_id=representative,
            maximum_displacement_mm=float(displacement.max()),source_distance_mm=float(source_distances[representative])))
    faces=remap[mesh.faces]
    keep=np.all(np.diff(np.sort(faces,axis=1),axis=1)>0,axis=1)
    cleaned=trimesh.Trimesh(mesh.vertices.copy(),faces[keep],process=False)
    cleaned.update_faces(cleaned.unique_faces());cleaned.remove_unreferenced_vertices()
    check=mesh_qc(cleaned,cfg)
    if not check['passed']:
        raise ValueError('NUMERICAL_SEAM_WELD_FAILED: '+repr(check))
    return cleaned,dict(method='LOCAL_BOOLEAN_VERTEX_WELD_WITHIN_EXISTING_NUMERICAL_TOLERANCE',
        tolerance_mm=tolerance,clusters=audit,merged_vertices=int(np.count_nonzero(remap!=np.arange(len(remap)))),
        removed_faces=int(len(mesh.faces)-len(cleaned.faces)),
        maximum_displacement_mm=max((r['maximum_displacement_mm'] for r in audit),default=0.),
        source_modified=False,outside_seam_vertices_modified=False,whole_model_remeshed=False,mesh_qc=check)
