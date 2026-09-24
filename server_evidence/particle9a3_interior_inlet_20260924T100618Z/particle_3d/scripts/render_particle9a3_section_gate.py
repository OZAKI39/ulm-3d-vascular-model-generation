"""Replot the actual blocked-section audit; never synthesize downstream results."""
from pathlib import Path
import csv
import json
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/'particle_3d/reports/particle9a3_interior_inlet'
FROZEN = ROOT/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference'


def surface(ax, mesh, color, alpha=1., mask=None):
    faces = mesh.faces.reshape(-1, 4)[:, 1:]
    triangles = np.asarray(mesh.points, float)[faces]*1e6
    if mask is not None: triangles = triangles[mask]
    ax.add_collection3d(Poly3DCollection(triangles, facecolors=color, edgecolors='none', alpha=alpha,
                                        rasterized=True))


def axes3d(ax, points, padding=.08):
    low, high = np.min(points, axis=0), np.max(points, axis=0)
    span = high-low
    low -= padding*span; high += padding*span
    ax.set(xlim=(low[0], high[0]), ylim=(low[1], high[1]), zlim=(low[2], high[2]),
           xlabel='x (µm)', ylabel='y (µm)', zlabel='z (µm)')
    ax.set_box_aspect(high-low)
    ax.view_init(elev=20, azim=-62)
    ax.tick_params(labelsize=8, pad=1)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_locator(plt.MaxNLocator(4))
        axis.pane.fill = False
        axis._axinfo['grid'].update(color='#dddddd', linewidth=.5)


def save(fig, name):
    for suffix in ('png', 'pdf'):
        fig.savefig(REPORT/f'figures/{name}.{suffix}', dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(fig)


def main():
    (REPORT/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titlesize': 12,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    data = json.loads((REPORT/'data/interior_injection_section.json').read_text())
    if data['selected_section'] is not None:
        raise RuntimeError('This renderer is explicitly for a blocked section search')
    topo = json.loads((REPORT/'data/root_topology.json').read_text())
    root = np.array(topo['points_m'])*1e6
    wall = pv.read(FROZEN/'SV_MESH/mesh-surfaces/WALL.vtp')
    cap = pv.read(FROZEN/'SV_MESH/mesh-surfaces/INLET.vtp')
    witness = pv.read(REPORT/'data/rejected_section_witness_00.vtp')
    fig = plt.figure(figsize=(11.5, 6))
    ax = fig.add_subplot(121, projection='3d')
    surface(ax, wall, '#7b8fa3', .16)
    surface(ax, cap, '#245fa7', .95)
    ax.plot(*root.T, color='#008a91', lw=2)
    ax.scatter(*root[-1], color='#ac2244', marker='D', s=32, depthshade=False)
    ax.text(*root[0], '  Open inlet', color='#245fa7', fontsize=9)
    ax.text(*root[-1], '  First junction\n  SWC node 15', color='#92213c', fontsize=9)
    axes3d(ax, np.asarray(wall.points)*1e6)
    ax.set_title('(a) Frozen vessel and authoritative root', pad=18)
    ax = fig.add_subplot(122, projection='3d')
    center = root[0]
    tri = wall.points[wall.faces.reshape(-1, 4)[:, 1:]]*1e6
    keep = np.linalg.norm(tri.mean(axis=1)-center, axis=1) < 4
    surface(ax, wall, '#7b8fa3', .17, keep)
    surface(ax, cap, '#245fa7', .55)
    surface(ax, witness, '#ed8926', .9)
    tangent = (root[1]-root[0])/np.linalg.norm(root[1]-root[0])
    ax.plot(*np.vstack([root[0], root[0]+3*tangent]).T, color='#008a91', lw=2)
    points = np.vstack([tri[keep].reshape(-1, 3), np.asarray(cap.points)*1e6])
    axes3d(ax, points, .04)
    ax.set_title('(b) First candidate — rejected', pad=18)
    ax.text2D(.06, .91, '0.138 µm from inlet\nFlux error: 0.8874% > 0.0001%', transform=ax.transAxes,
              color='#a0440a', fontsize=10)
    fig.legend(handles=[Patch(color='#245fa7', label='Open inlet cap'),
                        Line2D([0], [0], color='#008a91', lw=2, label='Source root centerline'),
                        Patch(color='#ed8926', label='Rejected candidate (not an injection section)')],
               loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(.5, .005), fontsize=9)
    fig.suptitle('Interior injection section search | No valid section', fontsize=15, y=.98)
    fig.subplots_adjust(left=.02, right=.95, top=.86, bottom=.12, wspace=.09)
    save(fig, '01_interior_injection_section')

    rows = list(csv.DictReader((REPORT/'data/interior_section_search.csv').open()))
    fig, axs = plt.subplots(2, 1, figsize=(9, 6.8), sharex=True, gridspec_kw={'hspace': .16})
    qref = data['reference_inlet_Q_m3_s']
    for ref, color, label in [(2, '#adb8c5', 'Refined stations (408)'), (1, '#17679b', 'Primary stations (204)')]:
        subset = [r for r in rows if int(r['refinement']) == ref]
        s = np.array([float(r['arclength_m'])*1e6 for r in subset])
        q = np.array([float(r['positive_Q_m3_s'])/qref for r in subset])
        err = np.array([float(r['flux_relative_error']) for r in subset])
        axs[0].plot(s, q, '.', ms=4 if ref == 1 else 3, color=color, label=label)
        axs[1].semilogy(s, err, '.', ms=4 if ref == 1 else 3, color=color)
    axs[0].axhline(1, color='#2b2b2b', lw=1, ls='--', label='Saved inlet flux')
    axs[0].set_ylabel(r'$Q_{\mathrm{section}}/Q_{\mathrm{inlet}}$')
    axs[0].set_title('(a) Actual P1 area flux at each candidate', loc='left')
    axs[0].legend(frameon=False, fontsize=9, loc='lower left')
    axs[1].axhline(data['flux_relative_tolerance'], color='#ac2244', ls='--', lw=1.5)
    axs[1].text(2, 1.4e-6, 'Existing tolerance = 10⁻⁶ (0.0001%)', color='#ac2244', fontsize=10)
    axs[1].set(xlabel='Root arclength from open inlet (µm)', ylabel='Relative flux error', ylim=(5e-7, .9))
    axs[1].set_title('(b) All 612 candidates fail the existing flux tolerance', loc='left')
    for ax in axs:
        ax.grid(True, color='#e4e4e4', lw=.5)
        ax.axvline(data['root_length_to_first_junction_m']*1e6, color='#008a91', ls=':', lw=1)
    axs[0].text(.99, .95, 'First junction →', transform=axs[0].transAxes, ha='right', va='top', color='#008a91')
    fig.suptitle('Interior-section flux validation | Frozen 2.0 mm/s FEM', fontsize=14, y=.97)
    fig.subplots_adjust(left=.12, right=.97, top=.91, bottom=.09)
    save(fig, '07_interior_section_flux_gate')


if __name__ == '__main__': main()
