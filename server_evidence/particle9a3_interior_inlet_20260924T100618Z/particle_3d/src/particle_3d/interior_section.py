"""P9-A.3 interior-section geometry gate; no particle or outlet-basin sampling.

All coordinates and velocities are the canonical tetrahedral P1 representation.
The infinite cutting plane may meet unrelated distant vessels: the root section
is the connected component containing its authoritative root-centerline point.
Both global and root component counts are recorded, never summed for root flux.
"""
from __future__ import annotations

import numpy as np
import pyvista as pv
import vtk

from .inlet_flux import positive_pieces

ALGORITHM_VERSION = 'P9A3_TETRA_P1_ROOT_SECTION_V1'


def root_topology(swc, cap_center_m):
    swc = np.asarray(swc, float)
    rows = {int(r[0]): r for r in swc}
    children = {i: [] for i in rows}
    roots = []
    for i, r in rows.items():
        parent = int(r[-1])
        if parent < 0:
            roots.append(i)
        else:
            children[parent].append(i)
    if len(roots) != 1:
        raise ValueError('Authoritative ROI must have exactly one root')
    path, current = [], roots[0]
    while current not in path:
        path.append(current)
        if len(children[current]) != 1:
            break
        current = children[current][0]
    if len(children[current]) < 2:
        raise ValueError('First junction is unresolved')
    points = np.vstack([np.asarray(cap_center_m, float), [rows[i][2:5]*1e-6 for i in path]])
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if np.any(lengths <= 0):
        raise ValueError('Degenerate root edge')
    root_edges = set(zip(path[:-1], path[1:]))
    daughters = [(rows[int(r[-1])][2:5]*1e-6, r[2:5]*1e-6)
                 for r in swc if int(r[-1]) >= 0 and (int(r[-1]), int(r[0])) not in root_edges]
    return dict(points_m=points, arclength_m=np.r_[0., np.cumsum(lengths)],
                root_node_ids=path, junction_node_id=current,
                junction_m=points[-1], daughter_segments_m=np.asarray(daughters))


def candidate_sequence(topology, spacing_m, refinement=1):
    """Cell-midpoint stations based on frozen mesh h10, independent of diameters.

    A refinement of two halves the station spacing; it is a search-resolution
    audit, not root-finding for a fortuitous flux crossing or a larger bubble.
    """
    if spacing_m <= 0 or refinement not in (1, 2):
        raise ValueError('Positive mesh spacing and refinement 1 or 2 required')
    arc = topology['arclength_m']
    step = spacing_m/refinement
    for index, s in enumerate(np.arange(step/2, arc[-1], step)):
        j = int(np.searchsorted(arc, s, side='right')-1)
        a, b = topology['points_m'][j:j+2]
        normal = (b-a)/(arc[j+1]-arc[j])
        yield dict(candidate_id=index, refinement=refinement, arclength_m=float(s),
                   center_m=a+(s-arc[j])*normal, normal=normal,
                   distance_to_first_junction_m=float(arc[-1]-s))


def cut_tetrahedra(grid, center, normal):
    plane = vtk.vtkPlane()
    plane.SetOrigin(*center)
    plane.SetNormal(*normal)
    cutter = vtk.vtkCutter()
    cutter.SetInputData(grid)
    cutter.SetCutFunction(plane)
    cutter.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)
    cutter.GenerateTrianglesOn()
    cutter.Update()
    section = pv.wrap(cutter.GetOutput()).copy()
    if section.n_cells:
        # Absolute, roundoff-scale point merge; never a geometric smoothing.
        scale = max(float(np.ptp(grid.points, axis=0).max()), 1e-30)
        section = section.clean(tolerance=64*np.finfo(float).eps*scale, absolute=True)
    return section


def locator_for(surface):
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()
    return locator


def distance_to_surface(locator, point):
    closest, cell, sub, distance2 = [0., 0., 0.], vtk.reference(0), vtk.reference(0), vtk.reference(0.)
    locator.FindClosestPoint(point, closest, cell, sub, distance2)
    return float(distance2)**.5, int(cell)


def root_component(section, center, tolerance):
    if not section.n_cells:
        return section, 0, False
    connected = section.connectivity()
    distance, cell = distance_to_surface(locator_for(connected), center)
    labels = np.asarray(connected.cell_data['RegionId'])
    selected = connected.extract_cells(labels == labels[cell]).extract_surface(algorithm='dataset_surface')
    return selected, len(np.unique(labels)), distance <= tolerance


def p1_flux(section, normal):
    if not section.n_cells:
        return dict(area_m2=0., positive_Q_m3_s=0., signed_Q_m3_s=0., finite_velocity=True)
    ids = section.faces.reshape(-1, 4)[:, 1:]
    xyz = np.asarray(section.points, float)[ids]
    values = np.asarray(section.point_data['Velocity'], float)[ids] @ normal
    area = np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1)/2
    finite = bool(np.isfinite(values).all())
    positive = sum(w for x, q in zip(xyz, values) for _, _, w in positive_pieces(x, q)) if finite else float('nan')
    return dict(area_m2=float(area.sum()), positive_Q_m3_s=float(positive),
                signed_Q_m3_s=float(np.sum(area*values.mean(axis=1))), finite_velocity=finite)


def boundary_points(section):
    if not section.n_cells:
        return np.empty((0, 3))
    edges = section.extract_feature_edges(boundary_edges=True, non_manifold_edges=False,
                                         feature_edges=False, manifold_edges=False)
    if not edges.n_cells:
        return np.empty((0, 3))
    pairs = edges.lines.reshape(-1, 3)[:, 1:]
    return np.vstack([edges.points, edges.points[pairs].mean(axis=1)])


def evaluate_candidate(grid, candidate, topology, wall_locator, cap_locators, reference_q, flux_tolerance):
    center, normal = candidate['center_m'], candidate['normal']
    roundoff = 512*np.finfo(float).eps*max(np.max(np.abs(grid.points)), 1e-30)
    whole = cut_tetrahedra(grid, center, normal)
    section, global_components, contains_root = root_component(whole, center, roundoff)
    data = dict(candidate)
    data.update(p1_flux(section, normal))
    perimeter = boundary_points(section)
    wall_error = max((distance_to_surface(wall_locator, p)[0] for p in perimeter), default=float('inf'))
    # Section/cap contact can only occur at its boundary on this conforming mesh.
    cap_distances = {k: min((distance_to_surface(loc, p)[0] for p in perimeter), default=float('inf'))
                     for k, loc in cap_locators.items()}
    before = bool(candidate['arclength_m'] < topology['arclength_m'][-1])
    # Also require every section vertex to be upstream of the junction in its local tangent frame.
    if section.n_points:
        before &= bool(np.min((topology['junction_m']-section.points)@normal) > roundoff)
    daughter_intersection = False
    if section.n_cells:
        locator = locator_for(section)
        for a, b in topology['daughter_segments_m']:
            da, db = (a-center)@normal, (b-center)@normal
            if da*db <= 0 and abs(da-db) > roundoff:
                p = a+da/(da-db)*(b-a)
                if distance_to_surface(locator, p)[0] <= roundoff:
                    daughter_intersection = True
    error = abs(data['positive_Q_m3_s']/reference_q-1)
    data.update(global_plane_components=global_components,
                root_component_count=int(section.n_cells > 0), root_center_contained=contains_root,
                wall_boundary_max_distance_m=wall_error, boundary_entirely_wall=bool(len(perimeter) and wall_error <= roundoff),
                open_boundary_intersection=any(d <= roundoff for d in cap_distances.values()),
                minimum_cap_distance_m=min(cap_distances.values()), before_first_junction=before,
                daughter_intersection=daughter_intersection, flux_relative_error=error,
                flux_tolerance=flux_tolerance, triangle_count=section.n_cells,
                geometric_roundoff_m=roundoff)
    reasons = []
    checks = [('ROOT_COMPONENT_UNRESOLVED', contains_root and data['root_component_count'] == 1),
              ('BOUNDARY_NOT_ENTIRELY_WALL', data['boundary_entirely_wall']),
              ('OPEN_CAP_INTERSECTION', not data['open_boundary_intersection']),
              ('NOT_BEFORE_FIRST_JUNCTION', before), ('DAUGHTER_INTERSECTION', not daughter_intersection),
              ('NONPOSITIVE_AREA', data['area_m2'] > 0), ('NONPOSITIVE_FLUX', data['positive_Q_m3_s'] > 0),
              ('SECTION_FLUX_MISMATCH', error <= flux_tolerance), ('NONFINITE_VELOCITY', data['finite_velocity'])]
    reasons.extend(reason for reason, passed in checks if not passed)
    data.update(accepted=not reasons, rejected_reason=';'.join(reasons),
                geometry_pass=all(passed for reason, passed in checks if reason != 'SECTION_FLUX_MISMATCH'))
    return data, section
