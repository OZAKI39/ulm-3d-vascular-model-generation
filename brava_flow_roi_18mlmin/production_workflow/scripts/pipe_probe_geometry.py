"""Strict fixed-position queries; distinguish actual outside points from locator tolerance."""
import numpy as np
import vtk


def frozen_probes(protocol, radius, length):
    levels = protocol['wall_distances_m'] + [v * radius for v in protocol['relative_wall_distances_R']]
    angles = 2 * np.pi * (np.arange(protocol['azimuth_count']) + protocol['azimuth_offset']) / protocol['azimuth_count']
    probes, groups, planes, azimuths = [], [], [], []
    for group, distance in enumerate(levels):
        for plane, z in enumerate(protocol['analysis_z_over_L']):
            for angle in angles:
                probes.append([(radius - distance) * np.cos(angle), (radius - distance) * np.sin(angle), z * length])
                groups.append(group); planes.append(plane); azimuths.append(angle)
    return np.array(levels), np.array(probes), np.array(groups), np.array(planes), np.array(azimuths)


def strict_cells(grid, points, tetra, probes):
    dx = points[tetra[:, 1:]] - points[tetra[:, :1]]
    cells = grid.find_containing_cell(probes)
    original = cells.copy()

    def contained(candidates, queries):
        bary = np.linalg.solve(dx[candidates].transpose(0, 2, 1), (queries - points[tetra[candidates, 0]])[..., None])[..., 0]
        weights = np.c_[1 - bary.sum(axis=1), bary]
        eps = 64 * np.finfo(float).eps * np.maximum(1, np.linalg.cond(dx[candidates]))
        valid = np.all(weights >= -eps[:, None], axis=1) & np.all(weights <= 1 + eps[:, None], axis=1)
        return valid, weights

    candidate = cells >= 0
    inside, _ = contained(cells[candidate], probes[candidate])
    cells[np.flatnonzero(candidate)[~inside]] = -1
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(grid)
    locator.BuildLocator()
    for index in np.flatnonzero(cells < 0):
        point = probes[index]
        ids = vtk.vtkIdList()
        locator.FindCellsWithinBounds(np.repeat(point, 2).tolist(), ids)
        candidates = np.array([ids.GetId(i) for i in range(ids.GetNumberOfIds())], dtype=int)
        if not len(candidates):
            continue
        inside, weights = contained(candidates, point)
        if np.any(inside):
            eligible = np.flatnonzero(inside)
            best = eligible[np.argmax(weights[eligible].min(axis=1))]
            cells[index] = candidates[best]
    return cells, original
