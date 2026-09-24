#!/usr/bin/env python3
"""Render baseline diagnosis only, then freeze the pre-repair diagnostic summary."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from fem3d.audit import sha256, timestamp, write_json
from fem3d.cap_remesh import triangle_quality

R = ROOT / 'reports/stage01_8'
O = ROOT / 'outputs/stage01_8'
I = ROOT / 'inputs/stage01_8'
if list(O.glob('repair_*/metadata/attempt.json')):
    raise RuntimeError('Cannot change pre-repair diagnosis after a repair attempt')
data = np.load(I / 'stage017_selected/mesh/volume_mesh.npz')
diagnosis = json.loads((O / 'diagnostics/residual_cells.json').read_text())
cells = diagnosis['cells']; count = len(cells)
points = data['points_m'] * 1e6
tri, tetra = data['boundary_triangles'], data['tetra']
surface = pv.PolyData(points, np.c_[np.full(len(tri), 3), tri].ravel())
grid = pv.UnstructuredGrid(np.c_[np.full(len(tetra), 4), tetra].ravel(), np.full(len(tetra), pv.CellType.TETRA, np.uint8), points)
wall = surface.extract_cells(data['facet_tags'] == 1).extract_surface()
ink, gray, red = '#233f52', '#87949c', '#ad3545'
images = {}


def plotter(shape, size):
    p = pv.Plotter(off_screen=True, shape=shape, window_size=size, border=False)
    for renderer in p.renderers:
        renderer.set_background('white')
    return p


def save3d(p, name, purpose):
    p.screenshot(str(R / name)); p.close()
    images[name] = {'sha256': sha256(R / name), 'purpose': purpose}


def camera(p, anchor, span):
    p.camera_position = [anchor + span * np.array([4., 4., 4.]), anchor, [0., 0., 1.]]
    p.enable_parallel_projection(); p.camera.parallel_scale = span * .70


def context(p, c):
    i = c['cell_index']; neighbors = c['one_ring_neighbor_tetra_ids']
    local = grid.extract_cells([i] + neighbors)
    p.add_mesh(surface.extract_cells(c['wall_diagnostics']['patch_triangle_ids']), color=gray, opacity=.14, show_edges=True, line_width=.7)
    if neighbors:
        p.add_mesh(grid.extract_cells(neighbors), color='#ccd2d6', opacity=.25, show_edges=True, edge_color=gray)
    p.add_mesh(grid.extract_cells([i]), color=red, opacity=.9, show_edges=True, line_width=2)
    p.add_text(f'Cell {i} | q = {c["min_sicn"]:.5f}', font_size=15, color=ink)
    p.add_text(f'Red: target; light gray: {len(neighbors)} face neighbors\nTransparent wall; actual shape, equal xyz scale', position='lower_left', font_size=10, color=ink)
    camera(p, local.points.mean(axis=0), float(np.linalg.norm(np.ptp(local.points, axis=0))))


p = plotter((1, 1), (1500, 1300))
p.add_mesh(wall, color='#bcc6cd', opacity=.42)
centers = np.array([c['centroid_m'] for c in cells]) * 1e6
p.add_points(centers, color=red, point_size=16, render_points_as_spheres=True)
# Nearby cells share a location label to avoid falsely implying separated positions.
remaining = {c['cell_index'] for c in cells}; by_id = {c['cell_index']: c for c in cells}
label_centers, labels = [], []
while remaining:
    cluster = {min(remaining)}; frontier = set(cluster)
    while frontier:
        following = {n for i in frontier for n in by_id[i]['one_ring_neighbor_tetra_ids'] if n in remaining} - cluster
        cluster |= following; frontier = following
    remaining -= cluster
    label_centers.append(np.mean([by_id[i]['centroid_m'] for i in cluster], axis=0) * 1e6)
    labels.append(' + '.join(map(str, sorted(cluster))) + (' (adjacent)' if len(cluster) > 1 else ''))
p.add_point_labels(np.array(label_centers), labels, font_size=15, text_color=ink, always_visible=True, shape_opacity=.8, show_points=False)
p.add_text(f'Where are the {count} residual cells?', font_size=20, color=ink)
p.add_text('Marker enlarged for visibility; markers are not tetra size\nCell IDs and exact coordinates are in residual_cells.json', position='lower_left', font_size=13, color=ink)
p.view_isometric(); p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(.86)
save3d(p, 'residual_locations.png', 'Locations on the actual vascular geometry; markers explicitly enlarged')

p = plotter((1, max(count, 1)), (max(1600, 650 * count), 850))
span = max(np.linalg.norm(np.ptp(points[tetra[c['cell_index']]], axis=0)) for c in cells)
for j, c in enumerate(cells):
    p.subplot(0, j); cell = grid.extract_cells([c['cell_index']])
    p.add_mesh(cell, color=red, show_edges=True, opacity=.85, line_width=2)
    p.add_points(cell.points, color=ink, point_size=8, render_points_as_spheres=True)
    p.add_text(f'Cell {c["cell_index"]} | actual tetra shape', font_size=14, color=ink)
    p.add_text(f'q = {c["min_sicn"]:.5f}; edge ratio = {c["edge_ratio"]:.2f}\nAngles {c["minimum_dihedral_degrees"]:.1f} to {c["maximum_dihedral_degrees"]:.1f} deg\nSame physical scale in every panel', position='lower_left', font_size=12, color=ink)
    camera(p, cell.points.mean(axis=0), span)
save3d(p, 'residual_cell_shapes.png', 'Actual tetra shapes with a common camera scale and equal xyz scaling')

p = plotter((1, max(count, 1)), (max(1800, 700 * count), 1000))
for j, c in enumerate(cells):
    p.subplot(0, j); context(p, c)
save3d(p, 'local_context_cells.png', 'Each target with its actual face-neighbor tetra and nearby fixed wall')
for c in cells:
    p = plotter((1, 1), (1250, 1100)); context(p, c)
    save3d(p, f'cell_{c["cell_index"]}_local_context.png', f'Local neighborhood of cell {c["cell_index"]}')

bq, _ = triangle_quality(data['points_m'], tri)
wq = bq[data['facet_tags'] == 1]
fig, axes = plt.subplots(1, 2, figsize=(13, 7))
ordered = np.sort(wq); axes[0].plot(ordered, np.arange(1, len(ordered) + 1) / len(ordered), color=gray, label='All Stage 1.7 wall triangles')
for c, color in zip(cells, ['#ad3545', '#d08a37', '#268579']):
    nearest = c['wall_diagnostics']['nearest_wall_triangle']; q = nearest['q_tri']
    axes[0].scatter(q, nearest['wall_quality_percentile'] / 100, color=color, s=65, label=f'Nearest wall: cell {c["cell_index"]}')
axes[0].set(xlabel='Wall triangle quality (larger is better)', ylabel='Cumulative fraction', xlim=(0, 1), ylim=(0, 1)); axes[0].legend(fontsize=9)
axes[1].boxplot([wq] + [bq[c['wall_diagnostics']['patch_triangle_ids']] for c in cells], labels=['All wall'] + [str(c['cell_index']) for c in cells], showfliers=False)
axes[1].set(ylabel='Wall triangle quality', xlabel='Two wall-edge rings around nearest triangle', ylim=(0, 1))
fig.suptitle('Are the nearby wall triangles unusual?', fontsize=18); fig.subplots_adjust(bottom=.18, top=.88, wspace=.25)
fig.text(.04, .035, 'The wall is unchanged. Local shape constraints also depend on fixed vertices and face orientations.', color=ink)
name = 'wall_surface_quality_near_residuals.png'; fig.savefig(R / name, dpi=150); plt.close(fig)
images[name] = {'sha256': sha256(R / name), 'purpose': 'Nearest/local wall triangle quality compared with the full frozen wall'}

ordinary = json.loads((O / 'diagnostics/ordinary_tetra_sample.json').read_text())['cells']
angles0 = [a['degrees'] for c in ordinary for a in c['dihedral_angles']]
fig, axes = plt.subplots(1, 2, figsize=(13, 7))
axes[0].hist(angles0, bins=np.linspace(0, 180, 37), density=True, color='#bdc8cf', label=f'{len(ordinary)} ordinary tetra')
for c, color in zip(cells, ['#ad3545', '#d08a37', '#268579']):
    angles = [a['degrees'] for a in c['dihedral_angles']]
    axes[0].plot(angles, np.full(len(angles), -.001 * (cells.index(c) + 1)), '|', ms=15, color=color, label=f'Cell {c["cell_index"]}')
    axes[1].plot(range(1, 7), sorted(angles), 'o-', color=color, label=str(c['cell_index']))
axes[0].set(xlabel='Internal dihedral angle (degrees)', ylabel='Ordinary sample density', xlim=(0, 180)); axes[0].legend(fontsize=9)
axes[1].axhline(np.degrees(np.arccos(1 / 3)), color=gray, ls='--', label='Regular tetra')
axes[1].set(xlabel='Six sorted angles per residual cell', ylabel='Internal dihedral angle (degrees)', ylim=(0, 180)); axes[1].legend(fontsize=9)
fig.suptitle('How flat are the residual tetrahedra?', fontsize=18); fig.subplots_adjust(bottom=.18, top=.88, wspace=.25)
fig.text(.04, .035, 'Ordinary sample: deterministic indices, quality at or above the whole-mesh median; not a random sample.', color=ink)
name = 'dihedral_angle_summary.png'; fig.savefig(R / name, dpi=150); plt.close(fig)
images[name] = {'sha256': sha256(R / name), 'purpose': 'All six residual internal angles and a reproducible ordinary-cell comparison'}

write_json(R / 'diagnosis_visualization_manifest.json', {'timestamp': timestamp(), 'images': images, 'mesh_sha256': diagnosis['mesh_sha256']})
summary = {'status': 'DIAGNOSIS_COMPLETE', 'timestamp': timestamp(), 'residual_count': count,
           'worst_quality': min(c['min_sicn'] for c in cells), 'mesh_sha256': diagnosis['mesh_sha256'],
           'policy_sha256': diagnosis['policy_sha256'], 'residual_cells_sha256': sha256(O / 'diagnostics/residual_cells.json'),
           'diagnostic_figures': images, 'repairs_started': False,
           'cells': [{'cell_index': c['cell_index'], 'reason': c['reason'], 'local_wall_quality': c['wall_diagnostics'],
                      'neighbor_quality': c['neighborhoods'], 'boundary_constraints': c['boundary_constraints']} for c in cells]}
write_json(R / 'diagnostic_summary.json', summary)
write_json(I / 'diagnostic_summary.json', summary)
write_json(I / 'residual_cells.json', diagnosis)
print('Diagnosis complete:', count, 'cells;', len(images), 'WSL figures. No repair has started.')
