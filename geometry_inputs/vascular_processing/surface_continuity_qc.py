"""Read-only local surface QC: VTK feature edges/cutters and trimesh adjacency.

The local cylindrical window and feature winding test exclude nearby branches and
small triangle cycles from the circumferential-ring measurement. No smoothing.
"""
import math

import networkx as nx
import numpy as np
import pyvista as pv
import trimesh
import vtk
from shapely.geometry import Polygon


def axes(normal):
    z = np.array(normal, dtype=float, copy=True); z /= np.linalg.norm(z)
    ref = np.eye(3)[np.argmin(abs(z))]
    x = np.cross(z, ref); x /= np.linalg.norm(x)
    return x, np.cross(z, x), z


def window(points, center, normal, radial, axial):
    d = np.asarray(points) - np.asarray(center)
    s = d @ axes(normal)[2]
    return (abs(s) <= axial) & (np.linalg.norm(d - s[:, None] * axes(normal)[2], axis=1) <= radial)


def surface(mesh):
    return pv.PolyData(mesh.vertices, np.column_stack((np.full(len(mesh.faces), 3), mesh.faces)).ravel())


def line_segments(poly):
    segments=[]; i=0
    while i < len(poly.lines):
        n=int(poly.lines[i]); ids=poly.lines[i+1:i+1+n]
        segments.extend(zip(ids[:-1], ids[1:])); i+=n+1
    return np.asarray(segments, dtype=int).reshape(-1, 2)


def feature_edges(mesh, center, normal, radial, axial, angle):
    detector=vtk.vtkFeatureEdges();detector.SetInputData(surface(mesh))
    detector.BoundaryEdgesOff();detector.NonManifoldEdgesOff();detector.ManifoldEdgesOff()
    detector.FeatureEdgesOn();detector.SetFeatureAngle(angle);detector.Update()
    result=pv.wrap(detector.GetOutput()); edges=line_segments(result)
    if len(edges):
        edges=edges[np.all(window(result.points, center, normal, radial, axial)[edges], axis=1)]
    graph=nx.Graph(); graph.add_edges_from(edges.tolist())
    lengths=np.linalg.norm(result.points[edges[:, 1]]-result.points[edges[:, 0]], axis=1) if len(edges) else np.array([])
    cycle_lengths=[]; rings=[]; x,y,_=axes(normal)
    for cycle in nx.cycle_basis(graph):
        points=result.points[cycle]; d=points-center
        theta=np.arctan2(d@y,d@x); delta=np.diff(np.r_[theta,theta[0]])
        winding=abs(float(np.sum(np.arctan2(np.sin(delta),np.cos(delta)))/(2*np.pi)))
        length=float(np.linalg.norm(np.roll(points,-1,axis=0)-points,axis=1).sum())
        cycle_lengths.append(length)
        if winding>.8: rings.append(length)
    component_lengths=[]
    for comp in nx.connected_components(graph):
        component_lengths.append(sum(float(np.linalg.norm(result.points[a]-result.points[b])) for a,b in graph.subgraph(comp).edges))
    shown=pv.PolyData(result.points.copy(),lines=np.column_stack((np.full(len(edges),2),edges)).ravel()) if len(edges) else pv.PolyData()
    return dict(angle_deg=float(angle),total_length_mm=float(lengths.sum()),edge_count=len(edges),
                closed_cycle_count=len(cycle_lengths),circumferential_ring_count=len(rings),
                longest_edge_mm=float(lengths.max()) if len(lengths) else 0.,
                longest_connected_feature_mm=max(component_lengths,default=0.),
                longest_circumferential_ring_mm=max(rings,default=0.)),shown


def normal_jumps(mesh, center, normal, radial, axial):
    # Trimesh provides manifold face adjacency; no separate edge detector.
    ids=mesh.face_adjacency
    local=window(mesh.triangles_center,center,normal,radial,axial)
    values=np.degrees(mesh.face_adjacency_angles[np.all(local[ids],axis=1)])
    if not len(values):raise ValueError('EMPTY_CONTINUITY_REGION')
    return dict(count=len(values),median_deg=float(np.median(values)),p90_deg=float(np.percentile(values,90)),
                p95_deg=float(np.percentile(values,95)),maximum_deg=float(values.max())),values


def cross_section(mesh, point, normal, maximum_center_distance=None):
    plane=vtk.vtkPlane();plane.SetOrigin(point);plane.SetNormal(normal)
    cut=vtk.vtkCutter();cut.SetInputData(surface(mesh));cut.SetCutFunction(plane);cut.Update()
    clean=vtk.vtkCleanPolyData();clean.SetInputConnection(cut.GetOutputPort());clean.Update()
    strip=vtk.vtkStripper();strip.SetInputConnection(clean.GetOutputPort());strip.Update()
    section=pv.wrap(strip.GetOutput());x,y,_=axes(normal);candidates=[];i=0
    while i<len(section.lines):
        n=int(section.lines[i]); ids=section.lines[i+1:i+1+n];i+=n+1
        points=section.points[ids]
        if len(points)<4 or np.linalg.norm(points[0]-points[-1])>1e-4:continue
        d=points-point; polygon=Polygon(np.column_stack((d@x,d@y)))
        if not polygon.is_valid or polygon.area<=0:continue
        centroid=np.asarray(point)+polygon.centroid.x*x+polygon.centroid.y*y
        distance=np.linalg.norm(centroid-point)
        if maximum_center_distance is not None and distance>maximum_center_distance:continue
        candidates.append((distance,polygon.area,centroid))
    if not candidates:raise ValueError('NO_VALID_LOCAL_CLOSED_SECTION')
    _,area,centroid=min(candidates,key=lambda v:v[0])
    return dict(area_mm2=float(area),equivalent_radius_mm=float(np.sqrt(area/math.pi)),
                centroid_x=float(centroid[0]),centroid_y=float(centroid[1]),centroid_z=float(centroid[2]))


def area_jump(rows):
    areas=np.array([v['area_mm2'] for v in rows])
    if len(areas)<2:return 0.
    return float(np.max(np.abs(np.diff(areas))/areas[:-1]))


def continuity_pass(features20, baseline_p95, current_p95, maximum_area_jump):
    return (features20['circumferential_ring_count']==0
            and (current_p95<=15 or current_p95<=.3*baseline_p95) and maximum_area_jump<=.1)


def circumferential_fraction(poly, center, normal, radius, half_width=.75, min_alignment=2**-.5):
    """Only near-plane edges directed within 45 degrees of the local azimuth.

    Longitudinal 24-gon facets are excluded. Length may exceed one circumference
    when there are two shoulders; the fraction is intentionally not clamped.
    """
    edges = line_segments(poly)
    if not len(edges):
        return dict(relevant_circumferential_length_mm=0., sharp_feature_circumference_fraction=0., relevant_edge_count=0)
    delta = poly.points[edges[:, 1]] - poly.points[edges[:, 0]]
    length = np.linalg.norm(delta, axis=1)
    mid = poly.points[edges].mean(axis=1) - np.asarray(center)
    normal = axes(normal)[2]
    axial = mid @ normal
    azimuth = np.cross(normal, mid)
    alignment = np.abs(np.einsum('ij,ij->i', delta, azimuth)) / np.maximum(length*np.linalg.norm(azimuth,axis=1),1e-15)
    valid = (np.abs(axial) <= half_width) & (alignment >= min_alignment) & (length > 1e-12)
    relevant = float(length[valid].sum())
    return dict(relevant_circumferential_length_mm=relevant,
                sharp_feature_circumference_fraction=relevant/(2*math.pi*radius), relevant_edge_count=int(valid.sum()))


def port_continuity(mesh, endpoint, cfg):
    """Identical seam localization for all ports, inside the configured 4 mm ROI.

    The outer context radius is NOT used to mix a neighboring vessel or a distant
    bend into seam statistics. A +/-1 mm axial band and 1.7*r radial envelope
    select the attachment, as in the previous real-data diagnosis.
    """
    settings = cfg['surface_continuity']
    radial = min(endpoint.radius*settings['radial_radius_fraction'], settings['local_radius_mm']/np.sqrt(2))
    axial = min(settings['seam_half_width_mm'], settings['local_radius_mm']/np.sqrt(2))
    normal, center = endpoint.tangent, endpoint.position
    features = []
    for angle in (10,20,30):
        metrics, poly = feature_edges(mesh, center, normal, radial, axial, angle)
        metrics.update(circumferential_fraction(poly, center, normal, endpoint.radius,
            settings['circumferential_plane_half_width_mm'], settings['circumferential_min_alignment']))
        features.append(metrics)
    normals, _ = normal_jumps(mesh, center, normal, radial, axial)
    rows = []
    step = settings['cross_section_step_mm']
    if not .1 <= step <= .2:
        raise ValueError('Cross-section step must be 0.1 to 0.2 mm')
    for s in np.linspace(-.5,1.,round(1.5/step)+1):
        rows.append(dict(arc_mm=float(s), **cross_section(mesh,center+s*normal,normal,endpoint.radius)))
    return dict(features=features, normal_jumps=normals, cross_sections=rows,
        maximum_area_jump_fraction=area_jump(rows),
        cross_section_coordinate='Signed distance along cap normal; not curved-centerline arclength',
        window=dict(context_radius_mm=settings['local_radius_mm'], radial_mm=radial, axial_half_length_mm=axial,
            center_mm=center.tolist(), normal=normal.tolist(),
            reason='Uniform seam band inside context ROI; excludes neighboring branches and remote bends'))


def all_port_gate(port_id, before, after, clearance):
    old20, new20 = before['features'][1], after['features'][1]
    oldp, newp = before['normal_jumps']['p95_deg'], after['normal_jumps']['p95_deg']
    errors=[]; warnings=[]
    if new20['circumferential_ring_count']:
        errors.append('PORT_SHARP_RING_' + port_id)
    if port_id == 'O3':
        if newp > oldp+1 or new20['total_length_mm'] > 1e-8:
            errors.append('O3_REGRESSION')
    elif not (newp <= .5*oldp or newp <= 20):
        errors.append('PORT_P95_IMPROVEMENT_FAILED_' + port_id)
    if clearance < (1.6 if port_id == 'O3' else 3.):
        errors.append('PORT_CLEARANCE_FAILED_' + port_id)
    if after['maximum_area_jump_fraction'] > .1:
        warnings.append('CROSS_SECTION_VARIATION_WARNING')
    return dict(passed=not errors, failures=errors, warnings=warnings,
        p95_reduction_fraction=1-newp/oldp if oldp else 0.,
        p95_absolute_threshold_met=newp <= 20,
        requires_visual_escalation_check=bool(new20['circumferential_ring_count'] or newp > 30))
