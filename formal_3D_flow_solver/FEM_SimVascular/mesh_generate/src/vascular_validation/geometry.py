"""Small surface adapters and direct geometric measurements for the vascular case."""
import numpy as np
import pyvista as pv
import vtk
from .validation import require

ROLE_IDS = {1: 'WALL', 2: 'OUTLET_03', 3: 'OUTLET_01', 4: 'INLET', 5: 'OUTLET_02'}


def polydata(points, triangles):
    triangles = np.asarray(triangles, np.int64)
    return pv.PolyData(np.asarray(points, float), np.c_[np.full(len(triangles), 3), triangles].ravel())


def triangles_of(surface):
    faces = np.asarray(surface.faces).reshape(-1, 4)
    require(np.all(faces[:, 0] == 3), 'All surface cells must be triangles')
    return faces[:, 1:]


def surface_stats(points, triangles):
    xyz = np.asarray(points)[triangles]
    vector = .5 * np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    areas = np.linalg.norm(vector, axis=1)
    require(np.all(areas > 0), 'Zero-area surface triangle')
    area = float(areas.sum())
    vector_area = vector.sum(axis=0)
    normal = vector_area / np.linalg.norm(vector_area) if np.linalg.norm(vector_area) > 0 else np.zeros(3)
    centroid = np.sum(xyz.mean(axis=1) * areas[:, None], axis=0) / area
    edges = np.sort(np.concatenate([triangles[:, [0, 1]], triangles[:, [0, 2]], triangles[:, [1, 2]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    rim = unique[counts == 1]
    perimeter = float(np.linalg.norm(points[rim[:, 1]] - points[rim[:, 0]], axis=1).sum())
    return {'area_m2': area, 'centroid_m': centroid.tolist(), 'normal': normal.tolist(),
            'projected_area_m2': float(np.linalg.norm(vector_area)), 'perimeter_m': perimeter,
            'equivalent_radius_m': float(np.sqrt(area / np.pi)),
            'hydraulic_diameter_m': 4. * area / perimeter if perimeter > 0 else None,
            'plane_deviation_m': float(np.max(np.abs((xyz.reshape(-1, 3) - centroid) @ normal)))}


def enclosed_volume(points, triangles):
    # Translate before summing to limit cancellation at the absolute model origin.
    points = np.asarray(points, np.float64)
    shifted = points - np.mean(points, axis=0)
    xyz = shifted[triangles]
    return float(np.sum(np.einsum('ij,ij->i', xyz[:, 0], np.cross(xyz[:, 1], xyz[:, 2]))) / 6.)


def closest_surface_distances(samples, target):
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(target)
    locator.BuildLocator()
    closest = [0., 0., 0.]
    cell_id = vtk.reference(0)
    sub_id = vtk.reference(0)
    distance_squared = vtk.reference(0.)
    result = np.empty(len(samples))
    for i, point in enumerate(samples):
        locator.FindClosestPoint(point, closest, cell_id, sub_id, distance_squared)
        result[i] = np.sqrt(float(distance_squared))
    return result


def surface_samples(surface):
    triangles = triangles_of(surface)
    return np.r_[surface.points, surface.points[triangles].mean(axis=1)]
