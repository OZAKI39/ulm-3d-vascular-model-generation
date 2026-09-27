#!/usr/bin/env python3
"""Freeze the one-candidate meshing scale before generation."""
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'mesh_generate/src'))
from sv_validation.mesh_diagnostics import unique_edges
from sv_validation.provenance import sha256, write_json, now

assert json.loads((ROOT / 'reports/sv1/geometry_import.json').read_text())['status'] == 'PASS'
target = ROOT / 'configs/mesh_policy.json'
assert not target.exists(), 'A frozen mesh policy must not be changed after generation'
volume = ROOT / 'inputs/fem_reference/volume_mesh.npz'
with np.load(volume) as data:
    points, tetra = data['points_m'], data['tetra']
edges = unique_edges(tetra)
h_ref = float(np.median(np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)))
with np.load(ROOT / 'inputs/fem_reference/exterior_surface.npz') as data:
    wall = data['triangles'][data['facet_tags'] == 1]
    we = np.sort(np.concatenate([wall[:, [0, 1]], wall[:, [0, 2]], wall[:, [1, 2]]]), axis=1)
    we = np.unique(we, axis=0)
    h_wall = float(np.median(np.linalg.norm(data['points_m'][we[:, 1]] - data['points_m'][we[:, 0]], axis=1)))
policy = {'timestamp': now(), 'source_mesh_sha256': sha256(volume), 'h_ref_m': h_ref,
          'h_ref_definition': 'median unique edge length of the selected Stage 1.7 tetrahedra',
          'h_wall_m': h_wall, 'primary_edge_factor': 1., 'one_fallback_factor': .8,
          'fallback_allowed_only_if': 'Primary mesh generation produces no usable volume',
          'formal_mesh_name': 'SV_MESH', 'surface_mesh_flag': True, 'volume_mesh_flag': True,
          'boundary_layers': False, 'anisotropy': False, 'adaptive_refinement': False,
          'geometry_distances': 'both directions; all vertices and triangle centroids to target triangles',
          'geometry_error_review': 'Quantitative distances, volume and port differences reported for human review',
          'hard_geometry_checks': ['complete original roles', 'no split or merged port', 'no normal flip'],
          'quality_metric': 'signed inverse Frobenius condition number for a linear tetrahedron',
          'quality_use': 'reported distribution for human review; no ranking or parameter sweep'}
write_json(target, policy)
write_json(ROOT / 'reports/sv1/mesh_policy_provenance.json', {'path': 'configs/mesh_policy.json', 'sha256': sha256(target)})
print(json.dumps(policy, indent=2))
