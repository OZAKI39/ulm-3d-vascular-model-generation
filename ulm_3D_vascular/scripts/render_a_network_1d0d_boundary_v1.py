"""Render only audited geometry; never invent pressure/flow maps for a blocked model."""
from pathlib import Path
import argparse
import csv
import json
import time
import resource
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.lines import Line2D


def render(output):
    started = time.perf_counter(); data = output/'data'; figures = output/'figures'; figures.mkdir(exist_ok=True)
    prov = json.loads((data/'a_network_provenance.json').read_text())
    mapping = json.loads((data/'roi_ports_in_a.json').read_text())
    with open(data/'a_nodes.csv') as f: nodes = list(csv.DictReader(f))
    with open(data/'a_edges.csv') as f: edges = list(csv.DictReader(f))
    xyz = np.array([[float(row[k]) for k in ['x_um', 'y_um', 'z_um']] for row in nodes])
    idx = {int(row['original_id']): i for i, row in enumerate(nodes)}
    segments = xyz[np.array([[idx[int(e['u'])], idx[int(e['v'])]] for e in edges])]
    selected = np.array([int(e['component_id']) == prov['analysis_component']['component_id'] for e in edges])
    roi = np.load(prov['lineage']['roi_archive']); roi_segments = roi['local_edge_points_um']
    plt.rcParams.update({'font.size': 10, 'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(15, 8.5)); left = fig.add_subplot(121, projection='3d'); right = fig.add_subplot(122, projection='3d')
    for ax in [left, right]:
        ax.set_facecolor('white'); ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); ax.set_zlabel('z (um)')
        ax.view_init(elev=23, azim=-62)
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]: axis.pane.fill = False
    left.add_collection3d(Line3DCollection(segments[~selected], colors='#bcc1c7', linewidths=.55, alpha=.4))
    left.add_collection3d(Line3DCollection(segments[selected], colors='#47698b', linewidths=.65, alpha=.65))
    left.add_collection3d(Line3DCollection(roi_segments, colors='#ed8b23', linewidths=2.8))
    terminals = np.array([int(n['degree']) == 1 for n in nodes]); ends = xyz[terminals]
    left.scatter(*ends.T, s=8, c='#75818b', alpha=.6)
    root_id = prov['source_root']['structural_root_ids'][0]; root = xyz[idx[root_id]]
    left.scatter(*root, s=95, marker='*', c='#ab357c', depthshade=False)
    left.text(*(root+np.array([3., 1., 7.])), '2410: structural root', fontsize=8, color='#8c2363')
    lo, hi = np.array(prov['roi_bbox_um']); corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
    box = [corners[[i, j]] for i in range(8) for j in range(i+1, 8) if np.count_nonzero(corners[i] != corners[j]) == 1]
    left.add_collection3d(Line3DCollection(box, colors='#dd7b20', linewidths=.9, linestyles='dashed'))
    left.set(xlim=(0, 195), ylim=(0, 195), zlim=(0, 385)); left.set_box_aspect([195, 195, 385])
    left.set_title('Full original SWC: 9,828 nodes / 43 components\nBlue: analysis component 42; orange: saved ROI', pad=22)
    right.add_collection3d(Line3DCollection(roi_segments, colors='#687c8a', linewidths=2.4))
    colors = ['#be303b', '#087f8c', '#4477bb', '#8855aa']
    for port, color in zip(mapping['ports'], colors):
        p = np.array(port['real_cut_xyz_um']); cap = np.array(port['fem_cap']['centroid_um'])
        for ax, size in [(left, 30), (right, 55)]: ax.scatter(*p, c=color, s=size, depthshade=False)
        right.plot(*np.array([p, cap]).T, color=color, linestyle='--', linewidth=1.6)
        right.scatter(*cap, c=color, s=36, marker='s', depthshade=False)
        offset = np.array([-7., 0., 6.]) if port['name'] == 'O1' else np.array([1.8, -2., 3.])
        right.text(*(p+offset), port['name'], color=color, fontsize=11, weight='bold')
    right.set(xlim=(75, 182), ylim=(38, 138), zlim=(75, 164)); right.set_box_aspect([107, 100, 89])
    right.set_title('Exact ROI: 3 edge cuts + node 4484 (O3)\nDashed lines: artificial FEM extensions', pad=22)
    handles = [Line2D([], [], color='#47698b', label='Analysis A (7,419 nodes)'), Line2D([], [], color='#bcc1c7', label='Other original components'),
               Line2D([], [], color='#ed8b23', linewidth=3, label='Saved connected ROI'),
               Line2D([], [], color='#75818b', marker='o', linestyle='', label='Graph endpoints (roles unknown)'),
               Line2D([], [], color='#ab357c', marker='*', linestyle='', markersize=10, label='Structural root; source unverified'),
               Line2D([], [], color='#687c8a', marker='s', linestyle='', label='Current FEM cap')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(.5, .045))
    fig.suptitle('A-Network v1 | Provenance recovered; hydraulic source unresolved', fontsize=17, y=.96)
    fig.text(.5, .018, 'Coordinates are physical um after raw voxel scaling (1, 1, 2). No blood-flow direction is inferred from the drawing.', ha='center', fontsize=10)
    fig.subplots_adjust(left=.02, right=.96, bottom=.17, top=.86, wspace=.06)
    for suffix in ['png', 'pdf']: fig.savefig(figures/f'01_full_A_roi_provenance.{suffix}', dpi=260)
    plt.close(fig)
    stats = dict(render_seconds=time.perf_counter()-started, peak_memory_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, rendered=['01_full_A_roi_provenance.png', '01_full_A_roi_provenance.pdf'], hydraulic_figures_generated=False)
    (output/'logs/render_metrics.json').write_text(json.dumps(stats, indent=2)+'\n')
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1]/'reports/a_network_1d0d_boundary_v1')
    render(p.parse_args().output)
