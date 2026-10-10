"""Scientific figures from actual geometry and result arrays, rendered on the CPU."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
os.environ['LIBGL_ALWAYS_SOFTWARE'] = '1'
os.environ['MPLCONFIGDIR'] = str(ROOT / 'outputs/mesh_and_flow/plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import pyvista as pv
from .geometry import ROLE_IDS, polydata, triangles_of

REPORT = ROOT / 'reports/mesh_and_flow'
FONT_PATH = Path('/mnt/c/Windows/Fonts/msyh.ttc')
FONT = FontProperties(fname=str(FONT_PATH)) if FONT_PATH.exists() else FontProperties()
COLORS = {'WALL': '#b6c6d0', 'INLET': '#146fcb', 'OUTLET_01': '#e98624',
          'OUTLET_02': '#21966c', 'OUTLET_03': '#9761c7'}


def save_plot(fig, name, title, note=None):
    fig.suptitle(title, fontproperties=FONT, fontsize=19, y=.985)
    if note:
        fig.text(.5, .018, note, ha='center', fontproperties=FONT, fontsize=10, color='#334155')
    fig.savefig(REPORT / name, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close(fig)


def save_scene(plotter, name, title, note=None):
    plotter.background_color = 'white'
    plotter.camera_position = 'iso'
    plotter.reset_camera()
    frame = plotter.screenshot(return_img=True)
    capabilities = plotter.render_window.ReportCapabilities()
    plotter.close()
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.imshow(frame)
    ax.axis('off')
    fig.subplots_adjust(top=.92, bottom=.06, left=0, right=1)
    save_plot(fig, name, title, note)
    return capabilities


def geometry_figure(points, triangles, tags, name, title, bc=False):
    plotter = pv.Plotter(off_screen=True, window_size=(1440, 1080))
    labels, centers = [], []
    for face_id, role in ROLE_IDS.items():
        selected = triangles[np.asarray(tags) == face_id]
        surface = polydata(points, selected)
        plotter.add_mesh(surface, color=COLORS[role], smooth_shading=False,
                         label=role, opacity=1. if role != 'WALL' else .8)
        if role != 'WALL':
            centers.append(points[np.unique(selected)].mean(axis=0))
            labels.append(role)
    plotter.add_point_labels(np.array(centers), labels, font_size=16, point_size=7,
                            text_color='#162b43', shape_color='white', always_visible=True)
    plotter.add_legend(bcolor='white', face='rectangle', size=(.18, .19))
    plotter.add_axes()
    note = 'SI 几何；端口身份来自原始标签' if not bc else '入口：指定总流量的恒定剖面；壁面：无滑移；三个出口：零自然牵引参考'
    save_scene(plotter, name, title, note)
