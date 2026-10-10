"""Loading figures from the existing native C1 run; no new solver states.

Arrows use the load force stored by LAMMPS, not total nodal force.  The full
saved-force history is checked against the actual current-position provider.
"""
from pathlib import Path
from pd_clot.portable_paths import relocate
import argparse
import hashlib
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image
import pyvista as pv
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pd_clot.lammps_fragmentation.current_position_provider import pulse
from pd_clot.lammps_fragmentation.fragment_analysis import read_dump
from analysis.clot_discretization_three_views_v1.render_three_views import silhouette_axes

OUT = ROOT / 'results/clot_loading_publication_v1'
V = ROOT / 'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1'
CASE = V / 'cases/C1_GAUSSIAN_A200'
FONT = Path(str(relocate('/mnt/c/Windows/Fonts/arial.ttf')))
C = dict(clot='#b85c64', anchor='#8aa7b9', debris='#586d7b', bond='#718593',
         grid='#93a3b0', axis='#647887', text='#263746', force='#b57920',
         peak='#ae4652', mean='#346a7e', ref='#788491', track='#315e61')
FIG_MM = (360, 220)
SCENE_BOX = [.018, .265, .515, .60]
ARROW_MM_PER_UN = .12


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def load_data():
    g = dict(np.load(V / 'V3_MATERIAL_POINTS_AND_BONDS.npz'))
    tr = dict(np.load(CASE / 'accepted_trajectory.npz'))
    ev = dict(np.load(CASE / 'broken_bonds.npz'))
    req = json.loads((CASE / 'request.json').read_text())
    provider = json.loads((CASE / 'LOAD_PROVIDER.json').read_text())
    history = np.loadtxt(CASE / 'force_history.dat')
    steps = tr['step']; time_s = steps * req['dt_s']; X = tr['position_m']
    center = np.asarray(req['load_centers_m'][0]); sigma = req['gaussian_sigma_m']
    distance = np.linalg.norm(X - center, axis=2)
    envelope = np.exp(-distance**2 / (2 * sigma**2))
    tangent = np.tile([1., 0., 0.], (len(g['X']), 1)) - g['normals'][:, 0, None] * g['normals']
    tangent[~g['surface']] = 0
    reconstructed = (g['areas'][None, :, None] * provider['scale_Pa'] *
                     pulse(time_s, req['duration_s'])[:, None, None] * envelope[:, :, None] * tangent)
    native = np.empty_like(reconstructed)
    max_position_error = 0.
    for frame, (step, cols, a) in enumerate(read_dump(CASE / 'trajectory.lammpstrj')):
        assert step == steps[frame]
        a = a[np.argsort(a[:, cols.index('id')])]
        assert np.array_equal(a[:, cols.index('id')], np.arange(1, len(g['X']) + 1))
        xx = a[:, [cols.index(k) for k in ['x', 'y', 'z']]]
        max_position_error = max(max_position_error, float(abs(xx - X[frame]).max()))
        native[frame] = a[:, [cols.index('v_load_f' + k) for k in 'xyz']]
    assert frame + 1 == len(steps)
    assert max_position_error < 1e-16
    force_error = float(abs(native - reconstructed).max())
    assert force_error < 1e-16
    force_mag = np.linalg.norm(native, axis=2)
    tau = force_mag / (g['areas'][None, :] + 1e-100)
    budget = force_mag.sum(axis=1) / g['areas'].sum()
    native_rows = history[steps.astype(int)]
    assert np.array_equal(native_rows[:, 0], steps)
    np.testing.assert_allclose(budget, native_rows[:, 1], rtol=5e-6, atol=1e-8)
    np.testing.assert_allclose(tau.max(axis=1), native_rows[:, 2], rtol=5e-6, atol=1e-8)
    assert np.max(abs(native[time_s >= req['duration_s']])) == 0.
    selected = int(tau.max(axis=1).argmax())
    assert steps[selected] == 540
    alive = np.ones(len(g['bonds']), bool)
    alive[ev['initial_bond_indices'][ev['events'][:, 1] <= steps[selected]]] = False
    # Native graph dumps exist only at step 0 and 4660. Validate those against
    # the same event reconstruction used for the intermediate selected frame.
    np.testing.assert_array_equal(np.sort(ev['events'][:, 2:4].astype(int)-1, axis=1),
                                  g['bonds'][ev['initial_bond_indices']])
    graph_dump_steps = []
    for step, cols, a in read_dump(CASE / 'surviving_bonds.local'):
        at_step = np.ones(len(g['bonds']), bool)
        at_step[ev['initial_bond_indices'][ev['events'][:, 1] <= step]] = False
        actual_bonds = np.sort(a[:, 1:3].astype(int) - 1, axis=1)
        order = np.lexsort((actual_bonds[:, 1], actual_bonds[:, 0]))
        np.testing.assert_array_equal(actual_bonds[order], g['bonds'][at_step])
        np.testing.assert_array_equal(a[order, 3].astype(int), g['bond_types'][at_step])
        graph_dump_steps.append(step)
    assert graph_dump_steps == [0, 4660]
    # Select a real initially strongly loaded node with greatest field decay
    # during the constant pulse plateau; not representative of every node.
    start = int(np.flatnonzero(steps == 600)[0]); end = int(np.flatnonzero(steps == 1740)[0])
    candidates = np.flatnonzero(g['surface'] & (tau[start] > 100.))
    track = int(candidates[np.argmin(envelope[end, candidates] / envelope[start, candidates])])
    assert track == 277
    chosen_track_frames = np.array([int(np.flatnonzero(steps == s)[0]) for s in [600, 800, 1200]])
    np.testing.assert_array_equal(pulse(time_s[chosen_track_frames], req['duration_s']), 1.)
    metrics = dict(case_id=req['case_id'], frame_step=int(steps[selected]),
                   frame_time_ms=float(time_s[selected]*1e3),
                   saved_frame_peak_traction_Pa=float(tau[selected].max()),
                   saved_frame_mean_traction_Pa=float(budget[selected]),
                   native_history_peak_traction_Pa=float(history[:, 2].max()),
                   native_history_peak_time_ms=float(history[history[:, 2].argmax(), 0]*req['dt_s']*1e3),
                   native_history_maximum_mean_traction_Pa=float(history[:, 1].max()),
                   reference_budget_Pa=req['tau_budget_Pa'], pulse_duration_ms=req['duration_s']*1e3,
                   reference_static_peak_traction_Pa=provider['reference_peak_traction_Pa'],
                   gaussian_sigma_mm=sigma*1e3, fixed_center_m=center.tolist(),
                   initial_exposed_surface_nodes=int(g['surface'].sum()),
                   clot_nodes=int(g['clot_count']), anchor_nodes=int(g['anchors'].sum()),
                   detached_nodes_in_selected_frame=int((~tr['attached'][selected]).sum()),
                   native_frame_surviving_bonds=int(alive.sum()),
                   graph_source='Initial bonds minus native break events with event step <=540; native endpoint dumps cross-checked at steps 0 and 4660',
                   maximum_position_crosscheck_error_m=max_position_error,
                   maximum_force_provider_crosscheck_error_N=force_error,
                   maximum_history_peak_roundoff_Pa=float(abs(tau.max(axis=1)-native_rows[:, 2]).max()),
                   native_force_history_rows=len(history), native_saved_frames=len(steps),
                   native_force_history_dt_us=req['dt_s']*1e6,
                   saved_trajectory_dt_us=req['dt_s']*req['dump_every']*1e6,
                   spatial_hard_cutoff=False, time_pulse_exactly_zero_after_duration=True,
                   surface_set_area_and_direction='Frozen initial exposed surface, area weights and projected axial tangents; no new fracture surface loading',
                   track_atom_id=track+1,
                   track_selection='Among surface nodes with traction >100 Pa at 0.30 ms, largest envelope reduction by 0.87 ms; illustrative selected node, not population average',
                   new_solver_runs=0)
    return dict(g=g, tr=tr, req=req, provider=provider, history=history,
                native=native, tau=tau, budget=budget, distance=distance, envelope=envelope,
                tangent=tangent, selected=selected, alive=alive, track=track,
                track_frames=chosen_track_frames, metrics=metrics)


def render_scene(d, dpi):
    size = np.rint(np.array(FIG_MM)/25.4 * dpi * np.array(SCENE_BOX[2:])).astype(int)
    X = d['tr']['position_m'][d['selected']]*1e3
    g = d['g']; forces = d['native'][d['selected']]; n = int(g['clot_count'])
    bounds = np.array([[-6.8, -5.2], [-.75, .85], [.55, 1.5]])
    focal = bounds.mean(axis=1); outward = np.array([.7, -2.2, -3.35]); outward /= np.linalg.norm(outward)
    right = np.cross(-outward, [0., 0., 1.]); right /= np.linalg.norm(right)
    up = np.cross(right, -outward); scale = 1.13
    in_box = np.all((X >= bounds[:, 0]) & (X <= bounds[:, 1]), axis=1)
    assert in_box.all(), 'Peak frame must show every node'
    segments = g['bonds'][d['alive']]; types = g['bond_types'][d['alive']]
    p = pv.Plotter(off_screen=True, window_size=size.tolist(), lighting='three lights')
    p.set_background('white'); p.render_window.SetMultiSamples(0); p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=100, occlusion_ratio=0.)
    grid = []
    for plane in range(3):
        fixed = bounds[plane, 0 if outward[plane] > 0 else 1]
        other = [k for k in range(3) if k != plane]
        for varying in other:
            constant = next(k for k in other if k != varying)
            for value in np.linspace(*bounds[constant], 4 if constant == 2 else 5):
                edge = np.zeros((2, 3)); edge[:, plane] = fixed; edge[:, constant] = value
                edge[:, varying] = bounds[varying]; grid.append(edge)
    grid_mesh = pv.PolyData(np.array(grid).reshape(-1, 3), lines=np.c_[np.full(len(grid), 2), np.arange(2*len(grid)).reshape(-1, 2)].ravel())
    p.add_mesh(grid_mesh, color=C['grid'], lighting=False, line_width=max(.8, dpi/180*.65), opacity=.90)
    axes = silhouette_axes(bounds, right, up)
    for axis in axes:
        p.add_mesh(pv.Line(*axis['edge']), color=C['axis'], lighting=False, line_width=max(.8, dpi/180*.65))
    for typ, opacity in [(1, .021), (2, .019), (3, .010)]:
        bb = segments[types == typ]
        mesh = pv.PolyData(X, lines=np.c_[np.full(len(bb), 2), bb].ravel())
        p.add_mesh(mesh, color='#9cabb5', opacity=opacity, line_width=max(.5, dpi/180*.42), lighting=False)
    projected = np.c_[(X-focal)@right, (X-focal)@up]; depth = (X-focal)@outward
    neighbors = cKDTree(projected).query_ball_point(projected, .055)
    front_distance = np.array([depth[js].max()-depth[i] for i, js in enumerate(neighbors)])
    front = front_distance < .045
    length = np.linalg.norm(X[segments[:, 1]]-X[segments[:, 0]], axis=1)
    emph = front[segments].all(axis=1) & (types <= 2) & (length <= .085)
    bb = segments[emph]; mesh = pv.PolyData(X, lines=np.c_[np.full(len(bb), 2), bb].ravel())
    p.add_mesh(mesh, color=C['bond'], opacity=.42, line_width=max(.85, dpi/180*.70), lighting=False)
    for color, indices, pt, alpha in [(C['anchor'], np.arange(n, len(X)), 1.45, .57),
                                     (C['clot'], np.arange(n), 1.70, .94)]:
        p.add_points(pv.PolyData(X[indices]), color=color, point_size=dpi*pt/72,
                     opacity=alpha, render_points_as_spheres=True, ambient=.55,
                     diffuse=.45, specular=.10, specular_power=20, show_scalar_bar=False)
    force_un = np.linalg.norm(forces, axis=1)*1e6
    candidates = np.flatnonzero(g['surface'] & (front_distance < .07) & (force_un >= .50))
    arrow_ids = []
    for i in candidates[np.argsort(-force_un[candidates])]:
        if arrow_ids and np.min(np.linalg.norm(projected[arrow_ids]-projected[i], axis=1)) < .12:
            continue
        arrow_ids.append(int(i))
        if len(arrow_ids) == 28:
            break
    assert len(arrow_ids) >= 12
    for i in arrow_ids:
        length = force_un[i]*ARROW_MM_PER_UN
        arrow = pv.Arrow(start=X[i], direction=forces[i], scale=length,
                         tip_length=min(.060/length, .36), tip_radius=.019/length,
                         shaft_radius=.0045/length, shaft_resolution=12, tip_resolution=20)
        p.add_mesh(arrow, color=C['force'], smooth_shading=True, ambient=.55, diffuse=.45)
    p.camera_position = [(focal+outward*20).tolist(), focal.tolist(), up.tolist()]
    p.enable_parallel_projection(); p.camera.parallel_scale = scale
    p.show(auto_close=False); p.render(); im = p.screenshot(return_img=True)
    peeling = bool(p.renderer.GetLastRenderingUsedDepthPeeling()); p.close()
    audit = dict(size=size.tolist(), box=SCENE_BOX, camera_focal_mm=focal.tolist(),
                 camera_outward=outward.tolist(), camera_right=right.tolist(), camera_up=up.tolist(),
                 parallel_scale_mm=scale, local_bounds_mm=bounds.tolist(), axis_edges=axes,
                 rendered_all_nodes=True, rendered_nodes=len(X), rendered_bonds=len(segments),
                 front_emphasized_real_bonds=int(emph.sum()), depth_peeling_used=peeling,
                 arrow_atom_ids=(np.array(arrow_ids)+1).tolist(), arrow_count=len(arrow_ids),
                 arrow_length_mm_per_microN=ARROW_MM_PER_UN,
                 arrow_semantics='Native nodal loading force; true 3D direction, common proportional length scale, starts at actual node center; sparse visible sample, not all loaded nodes',
                 glyph_semantics='Display symbols, not physical particle radii',
                 field_center_semantics='Fixed Gaussian maximum; no finite support boundary',
                 coordinate_conversion='m to mm only; no coordinate, shape or topology repair')
    d['arrow_ids'] = np.array(arrow_ids)
    return im, audit


def project(point, a, box):
    delta = np.asarray(point) - a['camera_focal_mm']; aspect = a['size'][0]/a['size'][1]
    u = .5 + delta@a['camera_right']/(2*a['parallel_scale_mm']*aspect)
    v = .5 + delta@a['camera_up']/(2*a['parallel_scale_mm'])
    return np.array([box[0]+box[2]*u, box[1]+box[3]*v])


def scene_labels(fig, ax, d, a, box, standalone=False):
    texts = []
    ticks = [[-6.8, -6., -5.2], [-.5, 0., .5], [.6, 1., 1.4]]
    for item in a['axis_edges']:
        k = item['axis']; normal = np.array(item['outward_screen']); base = np.mean(item['edge'], axis=0)
        for value in ticks[k]:
            point = base.copy(); point[k] = value
            label = f'{value:.1f}' if k == 0 else f'{value:g}'
            texts.append(ax.annotate(label, xy=project(point, a, box), xycoords=fig.transFigure,
                                     xytext=normal*18, textcoords='offset points', fontsize=16,
                                     color=C['text'], ha='center', va='center', annotation_clip=False))
        if k == 2:
            point = max(item['edge'], key=lambda q: project(q, a, box)[1]); offset = (3, 23)
        else:
            point = base; offset = normal*50
        texts.append(ax.annotate('XYZ'[k]+' (mm)', xy=project(point, a, box), xycoords=fig.transFigure,
                                 xytext=offset, textcoords='offset points', fontsize=18,
                                 color=C['text'], ha='center', va='center', annotation_clip=False))
    center = np.array(d['req']['load_centers_m'][0])*1e3
    xy = project(center, a, box)
    fig.add_artist(Line2D([xy[0]], [xy[1]], marker='+', markersize=11, markeredgewidth=1.8,
                          color=C['text'], linestyle='none', transform=fig.transFigure))
    target = (box[0]+box[2]*.035, box[1]+box[3]*.68)
    texts.append(ax.annotate('Fixed load\ncenter', xy=xy, xycoords=fig.transFigure,
                             xytext=target, textcoords=fig.transFigure, fontsize=15,
                             ha='left', va='center', color=C['text'], annotation_clip=False,
                             arrowprops=dict(arrowstyle='-', color=C['axis'], lw=.9,
                                             connectionstyle='arc3,rad=0')))
    # A true-length 2 microN key in the screen plane. Projected 3D arrows can
    # appear shorter because of perspective foreshortening (camera is parallel).
    key_y = box[1]+box[3]*.075; key_x = box[0]+box[2]*.12
    delta = box[3]*(2*ARROW_MM_PER_UN)/(2*a['parallel_scale_mm'])
    delta *= fig.get_figheight()/fig.get_figwidth()
    ax.annotate('', xy=(key_x+delta, key_y), xytext=(key_x, key_y),
                xycoords=fig.transFigure, textcoords=fig.transFigure,
                arrowprops=dict(arrowstyle='-|>', mutation_scale=15, lw=2.2, color=C['force'], shrinkA=0, shrinkB=0),
                annotation_clip=False)
    texts.append(fig.text(key_x+delta/2, key_y-.027, '2 µN', ha='center', fontsize=15, color=C['force']))
    return texts


def node_legend(fig, y, fontsize=15):
    handles = [Line2D([], [], marker='o', linestyle='none', markersize=7, color=C['clot'], label='Clot nodes'),
               Line2D([], [], marker='o', linestyle='none', markersize=7, color=C['anchor'], label='Fixed anchors'),
               Line2D([], [], color=C['bond'], lw=1.6, label='Surviving bonds'),
               Line2D([], [], marker='>', markevery=[1], color=C['force'], lw=2, markersize=7, label='Load force')]
    return fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, y),
                      ncol=4, frameon=False, fontsize=fontsize, handlelength=1.15,
                      columnspacing=1.4, handletextpad=.55)


def chart_style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    for k in ['left', 'bottom']:
        ax.spines[k].set_color(C['axis']); ax.spines[k].set_linewidth(.9)
    ax.tick_params(labelsize=15, colors=C['text'], width=.8, length=4)
    ax.grid(axis='y', color=C['grid'], alpha=.45, linewidth=.7)
    ax.set_axisbelow(True)


def compose_overview(im, a, d, dpi):
    fig = plt.figure(figsize=np.array(FIG_MM)/25.4, dpi=dpi, facecolor='white')
    texts = [fig.text(.025, .974, 'Loading at the current position', fontsize=26, va='top'),
             fig.text(.975, .974, 'C1 · prescribed Gaussian shear', fontsize=17, va='top', ha='right', color=C['axis']),
             fig.text(.045, .895, '(a) Nodal load · t = 0.270 ms', fontsize=20, color=C['text']),
             fig.text(.607, .895, '(b) Loading history', fontsize=20, color=C['text'])]
    ax = fig.add_axes(SCENE_BOX); ax.imshow(im, interpolation='none'); ax.axis('off')
    texts += scene_labels(fig, ax, d, a, SCENE_BOX)
    texts.append(fig.text(.295, .176, r'$\mathbf{F}_i^{\mathrm{load}}=A_i\,\tau(\mathbf{x}_i(t),t)\,\mathbf{t}_i^{0}$',
                          fontsize=24, ha='center', color=C['text']))
    texts.append(fig.text(.295, .131, 'Reference area and tangent; current-position field',
                          fontsize=15, ha='center', color=C['axis']))
    h = d['history']; times = h[:, 0]*d['req']['dt_s']*1e3; T = d['req']['duration_s']*1e3
    axh = fig.add_axes([.665, .60, .31, .235]); chart_style(axh)
    axh.plot(times, h[:, 2], color=C['peak'], lw=2.3, label='Actual local peak')
    axh.plot(times, h[:, 1], color=C['mean'], lw=2., label='Actual mean')
    axh.plot(times, 200*pulse(times/1e3, T/1e3), color=C['ref'], lw=1.8, ls='--', label='200 Pa × pulse')
    axh.axvline(T, color=C['axis'], lw=1., ls=':')
    axh.set(xlim=(0, 1.32), ylim=(0, 835), xticks=[0, .4, .8, 1.2], yticks=[0, 200, 400, 600, 800])
    axh.set_xlabel('Time (ms)', fontsize=17, color=C['text'], labelpad=8)
    axh.set_ylabel('Traction (Pa)', fontsize=17, color=C['text'], labelpad=7)
    axh.legend(loc='center', bbox_to_anchor=(.49, .58), frameon=False, fontsize=13.5, handlelength=1.5, labelspacing=.45)
    axh.scatter([.27], [d['metrics']['saved_frame_peak_traction_Pa']], s=28, color=C['peak'], zorder=5)
    axh.annotate('≈750 Pa', (.27, 750.258), xytext=(.035, 705), fontsize=14, color=C['peak'],
                 arrowprops=dict(arrowstyle='-', lw=.8, color=C['peak']))
    axh.text(1.21, 390, 'T = 1.1605 ms', fontsize=13, rotation=90, ha='center', va='center', color=C['axis'])
    texts.append(fig.text(.607, .462, '(c) A node leaving the loaded region', fontsize=20, color=C['text']))
    axd = fig.add_axes([.665, .174, .31, .222]); chart_style(axd)
    j = d['track']; rr = np.linspace(0, 1.45, 400)
    tau0 = d['provider']['scale_Pa']*np.linalg.norm(d['tangent'][j])
    axd.plot(rr, tau0*np.exp(-rr**2/(2*(d['req']['gaussian_sigma_m']*1e3)**2)), color=C['track'], lw=2.)
    frames = d['track_frames']; distances = d['distance'][frames, j]*1e3; yy = d['tau'][frames, j]
    axd.scatter(distances, yy, s=50, color=C['track'], edgecolors='white', linewidths=.8, zorder=4)
    for x, y, text, xytext in zip(distances, yy,
                                  ['0.30 ms', '0.40 ms', '0.60 ms\n0.017 Pa'],
                                  [(20, 12), (15, 25), (-5, 45)]):
        axd.annotate(text, (x, y), xytext=xytext, textcoords='offset points',
                     ha='right' if x > 1 else 'left', fontsize=14, color=C['track'],
                     arrowprops=dict(arrowstyle='-', lw=.8, color=C['track']))
    axd.set(xlim=(0, 1.45), ylim=(0, 760), xticks=[0, .4, .8, 1.2], yticks=[0, 200, 400, 600])
    axd.set_xlabel('Distance from fixed center (mm)', fontsize=17, color=C['text'], labelpad=8)
    axd.set_ylabel('Local traction (Pa)', fontsize=17, color=C['text'], labelpad=7)
    axd.text(.98, .95, 'Node 278 · pulse = 1\nσ = 0.28 mm', transform=axd.transAxes,
             fontsize=14, color=C['axis'], ha='right', va='top', linespacing=1.5)
    texts.extend([axh.xaxis.label, axh.yaxis.label, axd.xaxis.label, axd.yaxis.label])
    leg = node_legend(fig, .040, 16); texts.extend(leg.get_texts())
    texts.append(fig.text(.5, .020, 'Gaussian spatial decay; no hard cutoff. Arrows show a sparse sample of actual nodal load forces.',
                          fontsize=14, ha='center', color=C['axis']))
    return fig, texts


def compose_spatial(im, a, d, dpi):
    fig = plt.figure(figsize=np.array([220, 200])/25.4, dpi=dpi, facecolor='white')
    box = [.035, .25, .9, .9*220/(im.shape[1]/im.shape[0])/200]
    ax = fig.add_axes(box); ax.imshow(im, interpolation='none'); ax.axis('off')
    texts = [fig.text(.035, .975, 'Localized loading', fontsize=25, va='top'),
             fig.text(.955, .975, 't = 0.270 ms', fontsize=18, ha='right', va='top', color=C['axis'])]
    texts += scene_labels(fig, ax, d, a, box, True)
    texts.append(fig.text(.50, .163, '200 Pa reference budget · local peak 750.26 Pa',
                          fontsize=17, ha='center', color=C['text']))
    texts.append(fig.text(.50, .120, 'Pulse duration 1.1605 ms · Gaussian σ = 0.28 mm',
                          fontsize=16, ha='center', color=C['axis']))
    leg = node_legend(fig, .058, 14); texts.extend(leg.get_texts())
    texts.append(fig.text(.50, .027, 'Current-position load; fixed reference area and tangent',
                          fontsize=14, ha='center', color=C['axis']))
    return fig, texts


def save_figure(fig, texts, stem, dpi, preview):
    fig.canvas.draw(); renderer = fig.canvas.get_renderer(); w, h = fig.canvas.get_width_height()
    overflow = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        # Tick machinery may create hidden, out-of-range labels; only inspect
        # labels actually inside their axes range plus explicit text objects.
        bb = text.get_window_extent(renderer)
        if text in texts and (bb.x0 < 0 or bb.y0 < 0 or bb.x1 > w+1 or bb.y1 > h+1):
            overflow.append(text.get_text())
    assert not overflow, overflow
    boxes = [(t.get_text(), matplotlib.text.Text.get_window_extent(t, renderer)) for t in texts]
    overlaps = []
    for i, (label, b) in enumerate(boxes):
        for other, c in boxes[i+1:]:
            if min(b.x1, c.x1)-max(b.x0, c.x0) > 2 and min(b.y1, c.y1)-max(b.y0, c.y0) > 2:
                overlaps.append([label, other])
    assert not overlaps, overlaps
    name = stem + ('_preview' if preview else '')
    fig.savefig(OUT/(name+'.png'), dpi=dpi)
    if not preview:
        fig.savefig(OUT/(stem+'.pdf'), dpi=dpi, metadata={'Title':stem.replace('_', ' ')})
        fig.savefig(OUT/(stem+'.svg'), dpi=dpi)
        im = Image.open(OUT/(stem+'.png')).convert('RGB'); im.thumbnail((2400, 1600))
        im.save(OUT/(stem+'_preview.png'))
    plt.close(fig)
    return dict(filename=stem, size_pixels=[w, h], dpi=dpi, explicit_text_overflow=overflow,
                explicit_text_bbox_overlaps=overlaps)


def tables(d):
    k = d['selected']; g = d['g']; ids = np.arange(1, len(g['X'])+1); X = d['tr']['position_m'][k]
    shown = np.isin(np.arange(len(ids)), d['arrow_ids'])
    a = np.c_[ids, X, g['areas'], g['surface'].astype(int), d['native'][k], d['tau'][k], shown.astype(int)]
    np.savetxt(OUT/'selected_frame_nodal_loading.csv', a, delimiter=',', comments='',
               header='atom_id,x_m,y_m,z_m,reference_area_m2,initial_exposed_surface,load_fx_N,load_fy_N,load_fz_N,local_traction_Pa,arrow_shown')
    h = d['history']; t = h[:, 0]*d['req']['dt_s']; a = np.c_[h[:, 0], t, 200*pulse(t, d['req']['duration_s']), h[:, 1:]]
    np.savetxt(OUT/'native_loading_history.csv', a, delimiter=',', comments='',
               header='step,time_s,reference_budget_times_pulse_Pa,actual_mean_traction_Pa,actual_peak_traction_Pa,sum_load_force_magnitudes_N,net_load_fx_N,net_load_fy_N,net_load_fz_N')
    j = d['track']; tr = d['tr']; t = tr['step']*d['req']['dt_s']
    a = np.c_[tr['step'], t, tr['position_m'][:, j], d['distance'][:, j], d['envelope'][:, j],
              pulse(t, d['req']['duration_s']), d['tau'][:, j], d['native'][:, j], tr['attached'][:, j].astype(int)]
    np.savetxt(OUT/'tracked_node_278.csv', a, delimiter=',', comments='',
               header='step,time_s,x_m,y_m,z_m,distance_from_center_m,gaussian_weight,pulse,local_traction_Pa,load_fx_N,load_fy_N,load_fz_N,attached')


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--preview', action='store_true'); args = parser.parse_args()
    dpi = 140 if args.preview else 600; OUT.mkdir(parents=True, exist_ok=True)
    font_manager.fontManager.addfont(str(FONT))
    plt.rcParams.update({'font.family':'Arial', 'font.size':16, 'text.color':C['text'],
                         'mathtext.fontset':'stix', 'pdf.fonttype':42, 'svg.fonttype':'none'})
    sources = [V/'V3_MATERIAL_POINTS_AND_BONDS.npz', *[CASE/n for n in
               ['accepted_trajectory.npz', 'broken_bonds.npz', 'trajectory.lammpstrj',
                'surviving_bonds.local', 'request.json', 'LOAD_PROVIDER.json', 'in.lammps', 'force_history.dat']],
               ROOT/'analysis/clot_discretization_three_views_v1/render_three_views.py',
               ROOT/'results/clot_discretization_combined_v2/Fig_clot_discretization_three_views.png']
    before = {str(p):sha(p) for p in sources}
    d = load_data(); print('DATA_VERIFIED', json.dumps(d['metrics']), flush=True)
    # Cover the larger physical scene width of the standalone layout at the
    # requested export resolution as well as the overview's smaller panel.
    scene_dpi = int(np.ceil(dpi * 198/(FIG_MM[0]*SCENE_BOX[2])))
    im, a = render_scene(d, scene_dpi); a['scene_render_dpi'] = scene_dpi
    print('SCENE_COMPLETE', a['arrow_count'], 'force arrows', flush=True)
    figures = []
    for stem, compose in [('Fig_clot_loading', compose_overview), ('Fig_clot_loading_spatial', compose_spatial)]:
        fig, texts = compose(im, a, d, dpi)
        figures.append(save_figure(fig, texts, stem, dpi, args.preview))
        print('FIGURE_COMPLETE', stem, flush=True)
    assert before == {p:sha(p) for p in before}
    audit = dict(status='PASS', metrics=d['metrics'], rendering=a, figures=figures,
                 protected_source_sha256=before, colors=C, new_solver_runs=0,
                 formula='F_i = A_i * s * p(t) * exp(-||x_i(t)-c||^2/(2*sigma^2)) * (I-n_i^0 n_i^0) e_x',
                 tangent_is_unit_vector=False,
                 formula_scalar_tau='s*p(t)*Gaussian(x_i(t)); plotted traction is ||F_i||/A_i, including norm of projected tangent',
                 graph_c_semantics='Analytic prescribed spatial profile for the selected real node at pulse=1; markers are native force records at its actual current distances',
                 validation='All 234 native saved load-force arrays match reconstruction; every-step history agrees to its stored six-digit precision; native bond dump agrees with event reconstruction; source hashes unchanged')
    write_json(OUT/('PREVIEW_PROVENANCE.json' if args.preview else 'PROVENANCE.json'), audit)
    if not args.preview:
        tables(d)
        (OUT/'README.txt').write_text('''加载可视化：真实 C1_GAUSSIAN_A200 数据，沿用既有离散节点、存活键、锚层、白底与灰蓝网格风格。

Fig_clot_loading：主图加载箭头 + 原始加载时程 + 实际节点离开载荷中心后的衰减。
Fig_clot_loading_spatial：独立三维加载图，适合直接放入加载说明页。
PNG 为 600 dpi；PDF/SVG 保留可编辑文字与曲线，三维场景为高分辨率栅格。

主图采用真实 step 540（0.270 ms）的全部 1927 个血栓节点、1561 个锚节点和该时刻的存活键，此时尚无完全脱离的节点。形状来自加载期间的真实数据，因此与既有 2.330 ms 的最终形状不同。未挪动节点、重建光滑外壳或修改形状。箭头从真实节点出发，使用 LAMMPS 保存的 v_load_fx/y/z，而非总力或接触力；仅稀疏展示可见载荷箭头，长度按同一比例显示（1 µN 对应显示空间 0.12 mm）。箭头长度是力的显示比例，不是位移；三维投影存在缩短。节点符号不是物理球半径。

200 Pa 是初始参考几何上的面积加权平均力模预算；实际运动中的平均值并不固定为 200 Pa。保存帧最大局部牵引为 750.2578956 Pa（0.270 ms）。原始逐步时程的最大局部牵引为 750.473 Pa（0.2725 ms），最大实际平均值 183.605 Pa；脉宽 1.1605 ms，σ=0.28 mm。

载荷实现为 F_i=A_i*s*p(t)*g(x_i(t))*(I-n_i^0 n_i^0)e_x，s=780.232224778 Pa；参考面积、初始暴露面集合及投影方向冻结，当前坐标每步进入 Gaussian。公式中的 t_i^0 是投影后的轴向基向量，不是单位向量。曲线绘制的局部牵引为 ||F_i||/A_i，已包含投影向量的模。

重要文字校正：空间高斯场没有有限支持边界，不能写成“离开受力区立即归零”。建议写为“按当前节点位置计算；远离固定加载中心时迅速衰减”。时间脉冲在 1.1605 ms 后严格为零。

图(c)从真实节点中选取 0.30 ms 牵引>100 Pa、到0.87 ms高斯权重下降比例最大的节点278，用于演示该机制，不代表所有节点。0.30 / 0.40 / 0.60 ms 均处于恒定 pulse=1 的平台期，局部牵引分别约629.73 / 183.32 / 0.0170 Pa，因此衰减来自节点运动而非时间脉冲下降。曲线是该节点冻结方向下的解析空间场，标记为原始实际位置与力数据。

仅后处理已有计算；未新增求解，未改动加载规则。该载荷仍为工程解析规定剪切输入，不是声场或微泡求解反演的流体牵引。数据与校验见 PROVENANCE.json 及三个 CSV。
''')
        write_json(OUT/'DELIVERY_SHA256.json', {p.name:dict(bytes=p.stat().st_size, sha256=sha(p)) for p in OUT.iterdir() if p.suffix in ['.png', '.pdf', '.svg', '.csv']})
    print('COMPLETE', OUT, flush=True)


if __name__ == '__main__':
    main()
