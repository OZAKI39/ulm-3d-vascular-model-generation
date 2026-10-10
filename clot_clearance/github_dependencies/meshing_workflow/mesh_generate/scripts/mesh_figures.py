#!/usr/bin/env python3
"""Render the generated surface, actual tetrahedral cutaway and measured quality."""
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'mesh_generate/src'))
from vascular_validation.visuals import geometry_figure, save_scene, save_plot, FONT, plt, pv, np

assert json.loads((ROOT / 'reports/mesh_and_flow/mesh_validity.json').read_text())['status'] == 'PASS'
folder = ROOT / 'outputs/mesh_and_flow/solver_mesh'
with np.load(folder / 'mesh_arrays.npz') as data:
    points, triangles, tags, q = data['points_m'], data['boundary_triangles'], data['facet_tags'], data['min_sicn']
geometry_figure(points, triangles, tags, 'simvascular_surface.png', 'SimVascular 重新生成后的血管表面')
grid = pv.read(folder / 'mesh-complete.mesh.vtu')
surface = pv.read(folder / 'mesh-complete.exterior.vtp')
centers = grid.cell_centers().points
cut = grid.extract_cells(centers[:, 0] <= np.median(centers[:, 0]))
plotter = pv.Plotter(shape=(1,2), off_screen=True, window_size=(1800,1000))
plotter.subplot(0,0)
plotter.add_mesh(surface, color='#a6b8c4', opacity=.10)
plotter.add_mesh(cut.extract_surface(algorithm='dataset_surface'), color='#83bad1', show_edges=True, edge_color='#294957', line_width=.5)
focus = np.asarray(grid.points)[np.argmin(np.linalg.norm(np.asarray(grid.points)-np.asarray(grid.points).mean(axis=0),axis=1))]
plotter.add_point_labels([focus], ['Detail'], always_visible=True, font_size=16)
plotter.camera_position='iso'
plotter.reset_camera()
plotter.add_axes()
plotter.subplot(0,1)
h = json.loads((ROOT/'reports/mesh_and_flow/mesh_validity.json').read_text())['generation']['global_edge_size_m']
local = (np.linalg.norm(centers-focus,axis=1) <= 8*h) & ((centers-focus)@np.ones(3) <= 0)
local_grid=grid.extract_cells(local)
plotter.add_mesh(local_grid.extract_surface(algorithm='dataset_surface'), color='#83bad1', show_edges=True, edge_color='#294957', line_width=1.)
plotter.add_axes()
save_scene(plotter, 'mesh_cutaway.png', '血管内部四面体网格',
           '左：整体位置；右：局部放大。按单元中心移开一侧，保留真实四面体边界与内部剖面')
quality = json.loads((ROOT / 'reports/mesh_and_flow/mesh_quality.json').read_text())
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(q, bins=np.linspace(0, 1, 51), color='#397f9f', edgecolor='white')
ax.axvline(.1, color='#bb3c40', lw=1.5, linestyle='--', label='minSICN = 0.1')
ax.set_xlabel('minSICN (1 = ideal tetrahedron)')
ax.set_ylabel('四面体数量', fontproperties=FONT)
ax.legend(frameon=False)
ax.spines[['top', 'right']].set_visible(False)
note = 'q_min={:.5f}，P1={:.4f}，P5={:.4f}，median={:.4f}；低于 0.1：{} 个'.format(
    quality['q_min'], quality['P1'], quality['P5'], quality['median'], quality['N_low'])
fig.tight_layout(rect=(0, .07, 1, .93))
save_plot(fig, 'mesh_quality.png', '新网格的质量如何？', note)
print('Mesh figures generated from actual SV output')
