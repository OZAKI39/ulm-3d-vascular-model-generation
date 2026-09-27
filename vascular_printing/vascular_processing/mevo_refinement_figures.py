"""Same-camera centerline comparison figures with explicit raw cut diameters."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

from .boundary_review_figures import bounds, projected_labels

COLORS = {'Original candidate': '#bac1c7', 'Kept': '#008f9f', 'Removed': '#c62c3a',
          'M1 / upstream context': '#d98200', 'Proximal cut': '#2b9c38', 'Diameter cut': '#b529b7'}


def draw(ax, edges, coordinates, color, width=2.0, alpha=1.):
    segments = [[coordinates[a][:3], coordinates[b][:3]] for a, b in edges]
    if segments:
        ax.add_collection3d(Line3DCollection(segments, colors=color, linewidths=width, alpha=alpha))


def figure(path, title, coordinates, original, kept, context, proximal, diameter,
           limits, old_roots=(), *, show_removed=True):
    fig = plt.figure(figsize=(20, 13.34))
    # Explicit overlay order keeps coincident gray originals from hiding cyan.
    ax = fig.add_axes([.05, .12, .9, .8], projection='3d', computed_zorder=False)
    ax.view_init(elev=22, azim=-65)
    bounds(ax, [coordinates[n][:3] for n in limits])
    draw(ax, context, coordinates, COLORS['M1 / upstream context'], 1.7, .8)
    draw(ax, original, coordinates, COLORS['Original candidate'], 4, .7)
    if show_removed:
        draw(ax, set(original) - set(kept), coordinates, COLORS['Removed'], 2.0)
    draw(ax, kept, coordinates, COLORS['Kept'], 2.6)
    entries = []
    for n in sorted(set(old_roots)):
        ax.scatter(*coordinates[n][:3], marker='s', s=120, facecolors='none', edgecolors='#444444', linewidths=2, zorder=6)
        entries.append((coordinates[n][:3], f'Old root {n}'))
    for nodes, color, marker in [(proximal, COLORS['Proximal cut'], 'o'), (diameter, COLORS['Diameter cut'], 'X')]:
        for n in sorted(set(nodes)):
            radius = coordinates[n][3]
            ax.scatter(*coordinates[n][:3], color=color, marker=marker, s=100, depthshade=False, edgecolors='white', zorder=7)
            entries.append((coordinates[n][:3], f'ID {n}\nD = {2*radius:.3f} mm\nr = {radius:.3f} mm'))
    fig.canvas.draw()
    projected_labels(ax, entries, size=11)
    fig.suptitle(title, fontsize=21, y=.97)
    fig.legend(handles=[Line2D([0], [0], color=color, lw=3, label=name) for name, color in COLORS.items()],
               loc='lower center', ncol=3, fontsize=12, bbox_to_anchor=(.5, .04))
    fig.text(.5, .02, 'Original SWC centerlines in mm; marker labels use RAW radius. Fixed camera; line width is illustrative.\n'
             'Diameter filter is operational, not an anatomical M3/M4 endpoint. Isolated dips and cut endpoints may be below 0.75 mm.',
             ha='center', fontsize=11)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def export_figures(cache, originals, semantic, filtered, proximal_rows, diameter_rows, output):
    global_original, global_semantic, global_filtered = set(), set(), set()
    for component, old in originals.items():
        before = {e for graph in semantic[component] for e in graph.edges}
        after = {e for graph in filtered[component] for e in graph.edges}
        cuts = {r['refined_cut_original_swc_id'] for r in proximal_rows if r['component'] == component and r['status'] == 'SUSTAINED_MEVO_START'}
        dcuts = {r['cut_original_swc_id'] for r in diameter_rows if r['component'] == component and r['reason'] != 'DIAMETER_REBOUND'}
        root_keys = [k for k, c in cache.branch_component.items() if c == component
                     and cache.branch_component.get(cache.branches[k]['parent_branch']) != component]
        context_keys = {cache.branches[k]['parent_branch'] for k in root_keys}
        context = {e for k in context_keys if k in cache.branches
                   for e in zip(cache.branches[k]['node_ids'], cache.branches[k]['node_ids'][1:])}
        limits = set(old) | {n for e in context | before for n in e}
        old_roots = {n for n in old if old.in_degree(n) == 0}
        figure(output / f'part{component:02d}' / 'before_after_proximal.png',
               f'BG001 RMCA part{component:02d}: branch boundary to sustained point-level boundary',
               cache.coordinates, set(old.edges), before, context, cuts, [], limits, old_roots)
        figure(output / f'part{component:02d}' / 'before_after_diameter.png',
               f'BG001 RMCA part{component:02d}: semantic refinement to sustained D >= 0.75 mm filter',
               cache.coordinates, before, after, context, cuts, dcuts, limits)
        global_original.update(old.edges)
        global_semantic.update(before)
        global_filtered.update(after)
    m1 = {e for k, b in cache.branches.items() if b['label'] == 'M1'
          for e in zip(b['node_ids'], b['node_ids'][1:])}
    for name, title, old, kept in [
        ('original_candidate', 'Frozen original candidate', global_original, global_original),
        ('proximal_refined', 'Point-level proximal refinement', global_original, global_semantic),
        ('d075_filtered', 'Sustained D >= 0.75 mm filter', global_semantic, global_filtered)]:
        figure(output / f'BG001_RMCA_{name}.png', f'BG001 RMCA: {title}', cache.coordinates,
               old, set() if name == 'original_candidate' else kept, m1, [], [], set(cache.node_graph),
               show_removed=name != 'original_candidate')
