"""Portable evidence writers and explicit unknown terminal semantics."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, payload):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=convert, allow_nan=False)+'\n')


def write_csv(path, rows, fields=None):
    rows = list(rows)
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def terminal_audit(g, image_shape_xyz, spacing_um, analysis_component):
    """Image proximity is evidence of geometry, not certified inlet/outlet type.

    Use sample-center extent [0,(N-1)*spacing], matching the historical crop
    convention. Distance to actual voxel faces differs by half a voxel.
    """
    maximum = (np.asarray(image_shape_xyz)-1)*spacing_um
    rows = []
    for i in np.flatnonzero(g.degree == 1):
        axis_distance = np.minimum(g.xyz_um[i], maximum-g.xyz_um[i])
        rows.append(dict(original_id=int(g.ids[i]), component_id=int(g.component[i]),
                         in_analysis_component=bool(g.component[i] == analysis_component),
                         structural_root=bool(g.parent[i] == -1), terminal_class='UNKNOWN_TERMINAL',
                         x_um=g.xyz_um[i, 0], y_um=g.xyz_um[i, 1], z_um=g.xyz_um[i, 2],
                         radius_um=g.radius_um[i], image_boundary_distance_um=float(axis_distance.min()),
                         radius_ball_reaches_image_extent=bool(np.any(axis_distance <= g.radius_um[i])),
                         within_one_voxel_of_image_extent=bool(np.any(axis_distance <= spacing_um)),
                         deep_interior_candidate=bool(np.all(axis_distance > 2*g.radius_um[i]+2*np.asarray(spacing_um))),
                         reason='No endpoint boundary annotation; SWC parent/root and geometric proximity do not certify hydraulic role'))
    return rows
