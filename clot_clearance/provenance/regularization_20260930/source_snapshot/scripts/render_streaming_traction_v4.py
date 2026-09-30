"""V4, derived from render_damage_source_v3.py: local surface traction only.

Read-only decomposition at saved solver states; all interpolation is display-only.
The previous renderer, movies, configuration and simulation files stay untouched.
"""
import argparse
import copy
import csv
import hashlib
import json
import os
import shutil
import time
from pathlib import Path

import imageio.v2 as imageio
import imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
from matplotlib.legend_handler import HandlerPatch
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from PIL import Image
import pyvista as pv

from pd_clot.streaming import AnalyticTestStreaming, fluid_traction, pipe_field
from scripts.render_damage_forces_v1 import spaced_ids


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ArrowLegend(HandlerPatch):
    def create_artists(self, legend, orig_handle, xdescent, ydescent,
                       width, height, fontsize, trans):
        arrow = FancyArrowPatch((xdescent, height / 2 - ydescent),
                                (xdescent + width, height / 2 - ydescent),
                                arrowstyle='-|>', mutation_scale=fontsize,
                                linewidth=1.8, color=orig_handle.get_edgecolor())
        arrow.set_transform(trans)
        return [arrow]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--font', type=Path, default=Path('/mnt/c/Windows/Fonts/arial.ttf'))
    ap.add_argument('--preview-only', action='store_true')
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    run, out = a.run.resolve(), a.output.resolve()
    if out.exists():
        raise FileExistsError(f'Use a new directory: {out}')
    font_manager.fontManager.addfont(str(a.font))
    prop = font_manager.FontProperties(family='Arial', weight='normal')
    font = Path(font_manager.findfont(prop, fallback_to_default=False))
    assert sha(font) == sha(a.font)
    protect = [*root.glob('scripts/*.py'), *root.glob('pd_clot/*.py'), *root.glob('configs/*.json')]
    protect.extend(p for p in run.rglob('*') if p.is_file())
    protect.extend(p for p in (root / 'visualization/streaming_fragmentation_demo').rglob('*') if p.is_file())
    before = {str(p): sha(p) for p in sorted(set(protect))}
    out.mkdir(parents=True)
    (out / 'keyframes').mkdir()
    (out / 'provenance').mkdir()
    (out / 'provenance/INPUT_FILES_SHA256.json').write_text(json.dumps(before, indent=2) + '\n')
    shutil.copy2(root / 'scripts/render_damage_source_v3.py', out / 'provenance/render_damage_source_v3_original.py')
    shutil.copy2(Path(__file__), out / 'provenance/render_streaming_traction_v4.py')

    z = np.load(run / 'states.npz')
    c = json.loads((run / 'CONFIG.json').read_text())
    history = json.loads((run / 'history.json').read_text())
    identity = json.loads((run / 'IDENTITY.json').read_text())
    assert all(sha(root / p) == h for p, h in identity['source_sha256'].items())
    x, damage, fixed = z['x'], z['damage'], z['fixed']
    n = len(x)
    times = np.array([h['mechanical_time_s'] for h in history])
    local_config = copy.deepcopy(c['streaming'])
    local_config['include_pipe_traction'] = False  # In memory, for component extraction only.
    assert local_config['normal_traction_scale_Pa'] == 0
    local = AnalyticTestStreaming(local_config, c['pipe'], c['simulation']['representative_frequency_Hz'])
    normals, eligible, local_traction, pipe_traction, total_traction = [], [], [], [], []
    for k in range(n):
        mesh = pv.read(run / f'vtk/particles_{k:04d}.vtp')
        assert np.array_equal(mesh.points, x[k])
        assert np.array_equal(mesh['damage'], damage[k])
        normal = np.array(mesh['surface_normal'])
        mask = z['surface'][k] & z['attached'][k] & ~fixed & (np.linalg.norm(normal, axis=1) > 0)
        traction = np.zeros_like(x[k])
        traction[mask] = local.traction(x[k, mask], normal[mask], times[k])
        _, pressure, gradient = pipe_field(x[k], c['pipe'])
        background = fluid_traction(gradient, pressure, normal, c['pipe']['viscosity_Pa_s'])
        background[~mask] = 0
        normals.append(normal)
        eligible.append(mask)
        local_traction.append(traction)
        pipe_traction.append(background)
        total_traction.append(np.array(mesh['surface_traction_Pa']))
    normals, eligible = np.array(normals), np.array(eligible)
    local_traction, pipe_traction, total_traction = map(np.array, (local_traction, pipe_traction, total_traction))
    decomposition_error = float(np.abs(local_traction + pipe_traction - total_traction).max())
    normal_error = float(np.abs(np.sum(local_traction * normals, axis=2)).max())
    assert decomposition_error < 1e-10 and normal_error < 1e-10
    ids = spaced_ids(z['X'], np.flatnonzero(eligible.any(axis=0)), 32)
    scale_mm_per_Pa = .005
    offset_mm = c['clot']['particle_spacing_m'] * 1000 / 2
    np.savez_compressed(out / 'traction_samples.npz', local_streaming_traction_Pa=local_traction,
                        background_pipe_traction_Pa=pipe_traction, saved_total_surface_traction_Pa=total_traction,
                        surface_normals=normals, load_eligible=eligible, mechanical_time_s=times,
                        displayed_arrow_particle_ids=ids, positions_m=x, original_cycles=z['cycles'])
    fps, stride = 60, 30
    nf = (n - 1) * stride + 1
    index = np.arange(nf)
    left = np.minimum(index // stride, n - 1)
    right = np.minimum(left + 1, n - 1)
    alpha = (index % stride) / stride
    alpha[-1] = 0
    display_x = (1 - alpha[:, None, None]) * x[left] + alpha[:, None, None] * x[right]
    display_x[::stride] = x
    display_x[:, fixed] = z['X'][fixed]
    display_damage = damage[left]
    # Keep discrete load eligibility at the preceding saved state. Interpolate
    # vectors only if the particle is eligible at BOTH bracketing states.
    common = eligible[left] & eligible[right]
    interpolated = (1 - alpha[:, None, None]) * local_traction[left] + alpha[:, None, None] * local_traction[right]
    display_t = np.where(common[:, :, None], interpolated, local_traction[left])
    display_t[::stride] = local_traction
    display_mask = eligible[left][:, ids]
    # Offset glyphs by h/2 so they are legible at the model's particle surface.
    # This does not move the evaluation point or alter any solver coordinates.
    display_normals = (1 - alpha[:, None, None]) * normals[left] + alpha[:, None, None] * normals[right]
    display_normals /= np.maximum(np.linalg.norm(display_normals, axis=2, keepdims=True), 1e-30)
    display_normals = np.where(common[:, :, None], display_normals, normals[left])
    display_normals[::stride] = normals
    tails = display_x[:, ids] * 1000 + offset_mm * display_normals[:, ids]
    vectors = display_t[:, ids] * scale_mm_per_Pa
    assert np.array_equal(display_x[::stride], x)
    np.savez_compressed(out / 'display_samples.npz', display_positions_m=display_x,
                        display_damage=display_damage, display_streaming_traction_Pa=display_t,
                        source_state_left=left, source_state_right=right, interpolation_alpha=alpha,
                        video_time_s=index / fps, display_proxy_time_s=(1-alpha)*times[left]+alpha*times[right],
                        displayed_arrow_particle_ids=ids, arrow_visible=display_mask,
                        arrow_tails_mm=tails, arrow_vectors_mm=vectors, fixed=fixed)
    with (out / 'frame_map.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['frame', 'video_time_s', 'left_state', 'right_state', 'display_alpha', 'damage_state', 'visible_arrows', 'is_original_state'])
        for i in index:
            writer.writerow([int(i), i/fps, int(left[i]), int(right[i]), float(alpha[i]), int(left[i]), int(display_mask[i].sum()), bool(i % stride == 0)])

    prior = json.loads((root / 'visualization/streaming_fragmentation_demo/damage_source_v3/RENDER_MANIFEST.json').read_text())
    bounds = np.array(prior['camera']['limits_mm'])
    plt.rcParams.update({'font.family': 'Arial', 'font.weight': 'normal', 'font.size': 14,
                         'axes.labelweight': 'normal', 'axes.titleweight': 'normal'})
    fig = plt.figure(figsize=(19.2, 10.8), facecolor='#000000')
    ax = fig.add_axes([.015, .16, .95, .77], projection='3d', facecolor='#000000', computed_zorder=False)
    ax.set_proj_type('ortho')
    ax.view_init(elev=24, azim=-68)
    ax.set(xlim=tuple(bounds[0]), ylim=tuple(bounds[1]), zlim=tuple(bounds[2]), xlabel='X (mm)', ylabel='Y (mm)', zlabel='Z (mm)')
    ax.set_box_aspect(bounds[:, 1] - bounds[:, 0], zoom=1.70)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.label.set_color('white')
        axis.label.set_size(17)
        axis.line.set_color('#42556c')
        axis.set_pane_color((0, 0, 0, 1))
        axis._axinfo['grid']['color'] = '#24364c'
    ax.tick_params(colors='white', labelsize=14, pad=2)
    ax.xaxis.labelpad, ax.yaxis.labelpad, ax.zaxis.labelpad = 16, 12, 12
    colors = plt.colormaps['inferno'](damage[0])
    moving = ax.scatter(*(x[0, ~fixed] * 1000).T, c=colors[~fixed], s=38, edgecolors='#acb8c7', linewidths=.35, depthshade=False, zorder=5)
    base = ax.scatter(*(x[0, fixed] * 1000).T, c=colors[fixed], s=38, marker='s', edgecolors='#e6eaf0', linewidths=.65, depthshade=False, zorder=5)
    color = '#66e0c2'

    def arrows(i):
        mask = display_mask[i] & (np.linalg.norm(vectors[i], axis=1) > 1e-12)
        q = ax.quiver(*tails[i, mask].T, *vectors[i, mask].T, color=color, linewidth=1.7,
                      length=1, normalize=False, arrow_length_ratio=.25, pivot='tail', zorder=8)
        q.set_clip_on(False)
        return q

    quiver = arrows(0)
    for artist in [moving, base]:
        artist.set_clip_on(False)
    title = fig.text(.042, .944, 'Particle damage & streaming traction', color='white', fontsize=29, fontproperties=prop, weight='normal')
    cax = fig.add_axes([.245, .06, .51, .018])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap='inferno'), cax=cax, orientation='horizontal')
    cb.ax.tick_params(colors='white', labelsize=13, pad=3)
    cb.set_label('Particle damage', color='white', fontsize=16, labelpad=4)
    # Legend arrow identifies the component; physical scale is documented below.
    handles = [FancyArrowPatch((0, 0), (1, 0), color=color, label='Local streaming traction'),
               Line2D([], [], marker='s', markersize=9, markerfacecolor='none', markeredgecolor='white', linestyle='none', label='Fixed base')]
    fig.legend(handles=handles, handler_map={FancyArrowPatch: ArrowLegend()}, loc='center', bbox_to_anchor=(.5, .136),
               ncol=2, frameon=False, labelcolor='white', fontsize=16, handletextpad=.5, columnspacing=3)
    fig.canvas.draw()
    assert title.get_fontproperties().get_weight() == 'normal' and title.get_fontproperties().get_name() == 'Arial'
    reference = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
    points = np.concatenate([x.reshape(-1, 3) * 1000, tails[display_mask], (tails + vectors)[display_mask]])
    u, v, _ = proj3d.proj_transform(*points.T, ax.get_proj())
    px = ax.transData.transform(np.column_stack((u, v)))
    screen = np.array([px.min(axis=0), px.max(axis=0)])
    assert np.all(screen[0] > [20, 200]) and np.all(screen[1] < [1870, 945]), screen
    for artist in [moving, base, quiver]:
        artist.set_visible(False)
    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(fig.bbox)
    for artist in [moving, base, quiver]:
        artist.set_visible(True)

    def render(i):
        nonlocal quiver
        pos = display_x[i] * 1000
        colors = plt.colormaps['inferno'](display_damage[i])
        moving._offsets3d = tuple(pos[~fixed].T)
        moving.set_facecolor(colors[~fixed])
        base._offsets3d = tuple(pos[fixed].T)
        base.set_facecolor(colors[fixed])
        quiver.remove()
        quiver = arrows(i)
        fig.canvas.restore_region(background)
        for artist in [moving, base, quiver]:
            artist.do_3d_projection()
            ax.draw_artist(artist)
        return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()

    blit_error = float(np.abs(reference.astype(float) - render(0)).mean())
    assert blit_error < .1, blit_error
    if a.preview_only:
        for k in [0, 11, n - 1]:
            imageio.imwrite(out / f'preview_state_{k:04d}.png', render(k * stride))
        print(json.dumps(dict(preview_only=True, blit_pixel_MAE=blit_error, screen_bounds=screen.tolist())))
        plt.close(fig)
        return

    os.environ['IMAGEIO_FFMPEG_EXE'] = imageio_ffmpeg.get_ffmpeg_exe()
    started = time.perf_counter()
    gif_images = []
    with imageio.get_writer(out / 'particle_damage_streaming_traction_60fps.mp4', fps=fps, codec='libx264', quality=8,
                            macro_block_size=1, pixelformat='yuv420p', ffmpeg_params=['-movflags', '+faststart']) as writer:
        for i in index:
            im = render(i)
            writer.append_data(im)
            if i % stride == 0:
                imageio.imwrite(out / f'keyframes/state_{i // stride:04d}.png', im)
            if i % 3 == 0:
                gif_images.append(Image.fromarray(im).resize((1280, 720), Image.Resampling.LANCZOS))
            if i % 60 == 0 or i == nf - 1:
                print(json.dumps(dict(frame=int(i), total=nf, elapsed_s=time.perf_counter()-started)), flush=True)
    plt.close(fig)
    gif_images[0].save(out / 'particle_damage_streaming_traction_20fps.gif', save_all=True, append_images=gif_images[1:], duration=50, loop=0)
    for im in gif_images:
        im.close()
    shutil.copy2(out / 'keyframes/state_0000.png', out / 'preview.png')
    changed = [p for p, h in before.items() if not Path(p).is_file() or sha(Path(p)) != h]
    assert not changed
    manifest = dict(run=str(run), renderer=str(Path(__file__).resolve()), renderer_sha256=sha(Path(__file__)),
                    input_states_sha256=sha(run / 'states.npz'), original_state_count=n, particle_count=x.shape[1],
                    fps=fps, frame_count=nf, duration_s=nf/fps, original_endpoint_indices=(np.arange(n)*stride).tolist(),
                    positions='Piecewise linear display interpolation; exact saved positions every 30 frames; no extrapolation',
                    damage='Original recorded value held until the next macro endpoint; no interpolation',
                    traction='Local analytic_test component, evaluated at original saved particle centers and surface normals; background pipe stress excluded',
                    traction_units='Pa', traction_interpolation='Linear only for particles load-eligible at both endpoints; otherwise preceding vector held; masks held at preceding state',
                    interpolation_is_new_simulation=False, reconstructs_subcycle_oscillation=False,
                    traction_scale_mm_per_Pa=scale_mm_per_Pa, glyph_outward_offset_mm=offset_mm,
                    glyph_offset_is_display_only=True, arrow_material_ids=ids.tolist(),
                    max_component_reconstruction_error_Pa=decomposition_error, max_normal_component_Pa=normal_error,
                    peak_saved_local_traction_Pa=float(np.linalg.norm(local_traction, axis=2).max()),
                    surface_traction_arrows=True, total_particle_force_arrows=False, force_chart=False, force_numbers=False,
                    cycle_label=False, bubble_sphere=False, spatial_center_marker=False,
                    title='Particle damage & streaming traction', font=dict(family='Arial', weight='normal', path=str(font), sha256=sha(font)),
                    background='#000000', arrow_color=color,
                    camera=dict(elevation_deg=24, azimuth_deg=-68, projection='orthographic', limits_mm=bounds.tolist(), zoom=1.70),
                    displacement_scale=1, all_particles_retained=True, screen_bounds_pixels_bottom_origin=screen.tolist(),
                    collection_clipping=False, cached_background_pixel_MAE=blit_error,
                    elapsed_render_s=time.perf_counter()-started, protected_files=len(before), changed_preexisting_files=changed)
    (out / 'RENDER_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
