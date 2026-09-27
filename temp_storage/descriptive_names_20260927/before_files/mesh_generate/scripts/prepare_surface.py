#!/usr/bin/env python3
"""Convert exactly the tagged reference surface, without changing its geometry."""
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'mesh_generate/src'))
from sv_validation.geometry import polydata, ROLE_IDS
from sv_validation.provenance import write_json, sha256
from sv_validation.visuals import geometry_figure

assert json.loads((ROOT / 'reports/sv1/official_smoke.json').read_text())['status'] == 'PASS'
source = ROOT / 'inputs/fem_reference/exterior_surface.npz'
with np.load(source) as data:
    points, triangles, tags = data['points_m'], data['triangles'], data['facet_tags']
assert set(np.unique(tags)) == set(ROLE_IDS)
surface = polydata(points, triangles)
surface.cell_data['ModelFaceID'] = tags.astype(np.int32)
surface.cell_data['CellEntityIds'] = tags.astype(np.int32)
destination = ROOT / 'outputs/sv1/model/source.vtp'
destination.parent.mkdir(parents=True, exist_ok=True)
surface.save(destination)
write_json(ROOT / 'reports/sv1/surface_adapter.json', {
    'status': 'PASS', 'source_sha256': sha256(source), 'result_sha256': sha256(destination),
    'points': len(points), 'triangles': len(triangles), 'coordinate_transform': 'NONE',
    'coordinates_changed': False, 'connectivity_changed': False, 'units': 'm',
    'label_arrays': ['ModelFaceID', 'CellEntityIds']})
geometry_figure(points, triangles, tags, 'source_geometry.png', '输入到 SimVascular 的血管模型')
print('Source surface adapter and figure complete:', destination)
