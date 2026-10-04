"""P1 wall traction recovery in SI units; no solver/log/case imports.

Extracted unchanged algebra from production diagnostics SHA256
efa12e218bcaf99ba07db7194d3b8da79251da709c0a54a97ad6f6179be6e438.
Surface I/O imports PyVista only when called.
"""
import numpy as np


def p1_gradients(points, tetra, velocity):
    x = points[tetra]
    u = velocity[tetra]
    return np.swapaxes(np.linalg.solve(x[:, 1:] - x[:, :1], u[:, 1:] - u[:, :1]), 1, 2)


def face_keys(faces):
    a = np.ascontiguousarray(np.sort(faces, axis=1), dtype=np.int64)
    return a.view(np.dtype([('a', '<i8'), ('b', '<i8'), ('c', '<i8')])).ravel()


def boundary_owners(tetra, triangles):
    faces = tetra[:, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]].reshape(-1, 3)
    keys, first, counts = np.unique(face_keys(faces), return_index=True, return_counts=True)
    query = face_keys(triangles)
    index = np.searchsorted(keys, query)
    assert np.all(index < len(keys)), 'Boundary face absent from tetrahedra'
    assert np.array_equal(keys[index], query)
    assert np.all(counts[index] == 1), 'Boundary triangles must have exactly one adjacent element'
    return first[index] // 4


def wall_geometry(points, tetra, triangles, owners):
    xyz = points[triangles]
    cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    norm = np.linalg.norm(cross, axis=1)
    assert np.all(norm > 0)
    normal = cross / norm[:, None]
    center = xyz.mean(axis=1)
    direction = center - points[tetra[owners]].mean(axis=1)
    normal[np.einsum('ij,ij->i', normal, direction) < 0] *= -1
    assert np.all(np.einsum('ij,ij->i', normal, direction) > 0)
    return center, norm / 2, normal


def tangential_traction(gradient, normal, viscosity):
    """Tangential part of viscous stress times outward fluid normal [Pa]."""
    traction = viscosity * np.einsum('nij,nj->ni', gradient + gradient.swapaxes(1, 2), normal)
    return traction - np.einsum('ij,ij->i', traction, normal)[:, None] * normal


def nodal_average(triangles, values, areas, npoints):
    weight = np.bincount(triangles.ravel(), weights=np.repeat(areas, 3), minlength=npoints)
    total = np.bincount(triangles.ravel(), weights=np.repeat(areas * values, 3), minlength=npoints)
    result = np.zeros(npoints)
    np.divide(total, weight, out=result, where=weight > 0)
    return result, weight


def surface(points, triangles):
    import pyvista as pv
    ids, inv = np.unique(triangles, return_inverse=True)
    faces = np.column_stack((np.full(len(triangles), 3), inv.reshape(-1, 3)))
    mesh = pv.PolyData(points[ids], faces)
    mesh.point_data['GlobalNodeID_zero_based'] = ids
    return mesh, ids
