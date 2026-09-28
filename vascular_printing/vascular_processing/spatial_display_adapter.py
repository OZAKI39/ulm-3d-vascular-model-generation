"""Human spatial ROI display using the compact viewer's shared annotations.

The saved branch geometry and ROIRecord objects are read unchanged. Sampling,
radius interpolation, structural markers and cluster membership stay with the
original human workflow; this adapter only owns the window and its annotations.
"""
from datetime import datetime, timezone
import json
from pathlib import Path

from utils.rodent_vasculature import interactive as ui
from .compact_annotations import DualPanelAnnotations, display_font
from .topbrain_qc import write_json


Y_UP = (0., 1., 0.)
# Keep every existing structural marker in the single-row legend. These short
# names fit the same font/viewport as compact without implying MeVO anatomy.
SPATIAL_LEGEND_NAMES = {
    'SWC centerline': 'Centerline',
    'structural root': 'Root',
    'structural leaf': 'Terminal',
    'divergence junction': 'Branch point',
    'convergence junction': 'Merge point',
    'candidate ROI': 'Spatial ROI',
    'active ROI': 'Active ROI',
}


class SpatialROIViewer(DualPanelAnnotations):
    """Original human scene, with shared compact typography and camera rules."""

    def __init__(self, sample_root, sampling_run_root, *, max_arrows=600,
                 volume_opacity=.32, window_size=(1800, 900), show=True):
        import pyvista as pv

        self.dataset_dir = Path(sample_root).resolve()
        self.sampling_run_root = Path(sampling_run_root).resolve()
        self.run = self.sampling_run_root / 'figures'
        self.run.mkdir(parents=True, exist_ok=True)
        manifest = json.loads((self.dataset_dir / 'preprocess_manifest.json').read_text())
        self.sample_id = str(manifest['record']['sample_id'])
        if manifest.get('normalized_volume_path'):
            volume, _ = ui.load_normalized_volume(Path(manifest['normalized_volume_path']))
        elif manifest['record'].get('image_path'):
            volume = ui.load_tiff_volume(Path(manifest['record']['image_path']))
        elif manifest['record'].get('mask_path'):
            volume = ui.load_tiff_volume(Path(manifest['record']['mask_path']))
        else:
            volume = None
        self.geometry = ui._geometry_from_exports(self.dataset_dir / 'graphs', max_arrows)
        self.records = tuple(ui.load_sampling_display_rois(self.sampling_run_root))
        if not self.records:
            raise ValueError(f'No spatial ROI geometry in {self.sampling_run_root}')
        self.selected_indices = [i for i, roi in enumerate(self.records) if roi.is_representative]
        self.cluster_ids = sorted({roi.cluster_id for roi in self.records if roi.cluster_id >= 0})
        self.index = min(self.selected_indices or range(len(self.records)),
                         key=lambda i: self.records[i].selection_rank
                         if self.records[i].selection_rank > 0 else i + 100000)
        self.initial_roi = self.records[self.index].roi_id
        self.view = 'spatial'
        self.static_filename_prefix = f'{self.sample_id}_{self.initial_roi}_{self.view}'
        self.static_provenance = dict(sample_root=str(self.dataset_dir),
                                      sampling_run_root=str(self.sampling_run_root),
                                      roi_count=len(self.records),
                                      selected_count=len(self.selected_indices))
        self.group_index = 0
        self.events = []
        self.controllers = [None, None]
        self.orientation_labels = [None, None]
        self.legend_layout_sizes = [None, None]
        self.overlay_actors = []
        self.actor_to_index = {}
        self.diameter_clim_um = (ui._sampling_diameter_limits_um(self.records)
                                if ui.DIAMETER_COLORBAR_USE_GLOBAL_RANGE else None)
        self.plotter = pv.Plotter(shape=(1, 2), border=True, border_color='#606060',
                                  off_screen=not show, window_size=window_size)
        self.plotter.theme.font.family = ui.UI_FONT_FAMILY
        self.scale_alignment_observer = self.plotter.renderers[1].AddObserver(
            'EndEvent', self.align_bottom_scale)
        self.plotter.subplot(0, 0)
        self.scene_metadata, self.controllers[0] = ui._add_full_scene(
            self.plotter, volume, self.geometry,
            spacing_xyz_um=tuple(manifest['spacing_xyz_um']), volume_opacity=volume_opacity,
            sample_id=self.sample_id, sampling_available=True, left_view_up=Y_UP)
        self.plotter.remove_actor('full_scene_title', reset_camera=False, render=False)
        for name, actor in self.plotter.renderer.actors.items():
            if name.startswith('full_scene_legend_text_'):
                actor.SetInput(SPATIAL_LEGEND_NAMES.get(actor.GetInput(), actor.GetInput()))
        self.plotter.camera.Zoom(.90)
        self.style_annotations(0)
        self.redraw_boxes(self.selected_indices or list(range(len(self.records))),
                          'selected representatives')
        self.select(self.index)
        self.actions = {
            'a': lambda: self.redraw_boxes(list(range(len(self.records))), 'all candidates'),
            'r': self.show_selected, 's': self.show_selected, 'c': self.next_cluster,
            'Left': lambda: self.select((self.index - 1) % len(self.records)),
            'Right': lambda: self.select((self.index + 1) % len(self.records)),
            'PageUp': lambda: self.select((self.index - 1) % len(self.records)),
            'PageDown': lambda: self.select((self.index + 1) % len(self.records)),
            'F12': self.screenshot,
        }
        for key, callback in self.actions.items():
            for variant in ((key, key.upper()) if len(key) == 1 else (key,)):
                self.plotter.clear_events_for_key(variant)
                self.plotter.add_key_event(variant, callback)
        self.plotter.enable_mesh_picking(callback=self.pick, show=False, show_message=False,
                                        left_clicking=True, use_actor=True)
        self.rotation = None
        if show:
            self.rotation = ui._install_synchronized_horizontal_rotation(
                self.plotter, tuple(r.camera for r in self.plotter.renderers),
                left_view_up=Y_UP, right_view_up=Y_UP)

    def redraw_boxes(self, indices, mode):
        self.plotter.subplot(0, 0)
        for actor in self.overlay_actors:
            self.plotter.remove_actor(actor, reset_camera=False, render=False)
        self.overlay_actors, self.actor_to_index = ui._add_sampling_boxes(
            self.plotter, self.records, indices, mode_label=mode)
        for name in ('sampling_layer_mode', 'sampling_roi_labels-points', 'sampling_roi_labels-labels'):
            self.plotter.remove_actor(name, reset_camera=False, render=False)
        self.visible_indices = list(indices)
        self.events.append(f'BOXES:{mode}:{len(indices)}')
        self.plotter.render()

    def show_selected(self):
        self.redraw_boxes(self.selected_indices or list(range(len(self.records))),
                          'selected representatives')

    def next_cluster(self):
        if self.cluster_ids:
            cluster_id = self.cluster_ids[self.group_index % len(self.cluster_ids)]
            self.group_index += 1
            self.redraw_boxes([i for i, r in enumerate(self.records) if r.cluster_id == cluster_id],
                              f'cluster {cluster_id}')

    def pick(self, actor):
        index = self.actor_to_index.get(getattr(actor, 'memory_address', ''))
        if index is not None:
            self.select(index)

    def select(self, index):
        suppressed = self.plotter.suppress_rendering
        self.plotter.suppress_rendering = True
        try:
            self.index = index
            self.plotter.subplot(0, 0)
            ui._add_sampling_active_outline(self.plotter, self.records[index])
            self.plotter.subplot(0, 1)
            first = self.controllers[1] is None
            if not first:
                self.controllers[1].dispose()
            ui._remove_diameter_colorbar(self.plotter)
            self.plotter.renderer.clear_actors()
            self.controllers[1] = ui._add_sampling_roi_scene(
                self.plotter, self.records[index], add_orientation_axes=first,
                diameter_clim_um=self.diameter_clim_um, view_up=Y_UP)
            # The spatial 30x30x40-mm box fills more height than the compact
            # subtree. Leave room for the same larger ticks and bottom scale.
            self.plotter.camera.Zoom(.70)
            self.plotter.remove_actor('sampling_roi_information', reset_camera=False, render=False)
            self.style_annotations(1)
            self.events.append('ROI:' + self.records[index].roi_id)
        finally:
            self.plotter.subplot(0, 0)
            self.plotter.suppress_rendering = suppressed
        self.plotter.render()

    def screenshot(self):
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
        path = self.run / f'spatial_{self.records[self.index].roi_id}_{stamp}.png'
        self.plotter.screenshot(str(path))
        print(f'SPATIAL_SCREENSHOT: {path}', flush=True)
        return path

    def close(self):
        self.plotter.renderers[1].RemoveObserver(self.scale_alignment_observer)
        for controller in (*self.controllers, *self.orientation_labels):
            if controller is not None:
                controller.dispose()
        self.plotter.close()

    def run_window(self, *, show=True, screenshot_path=None, smoke_seconds=0):
        preview = Path(screenshot_path) if screenshot_path else self.run / 'interactive_sampling_layer_preview.png'
        preview.parent.mkdir(parents=True, exist_ok=True)
        report = dict(source='human-spatial', sample_id=self.sample_id,
                      **self.static_provenance, initial_roi=self.initial_roi,
                      global_branch_count=self.geometry.branch_count,
                      global_arrow_count=len(self.geometry.arrow_points_um),
                      preview=str(preview), font=display_font()[0],
                      display_units='mm', decimal_places=2, geometry_units='um',
                      both_view_up=list(Y_UP), coordinate_tick_groups_per_axis=1,
                      candidate_legend='Spatial ROI', core_recomputed_in_callback=False,
                      legend_layout='centred_single_row_below_top', legend_background_visible=False,
                      top_information_text_visible=False, roi_cluster_labels_visible=False,
                      diameter_colorbar_layout='bottom_horizontal',
                      status='GUI_STARTED' if show else 'OFFSCREEN_RENDERED')
        if show and smoke_seconds:
            def smoke(step):
                if step == 4:
                    for key in ('a', 'c', 'r', 's', 'Right', 'Left'):
                        for callback in self.plotter.iren._key_press_event_callbacks[key]:
                            callback()
                        self.events.append('KEY:' + key)
                if step >= max(8, int(smoke_seconds * 10)):
                    self.events.append('GUI_SMOKE_AUTO_CLOSE')
                    self.plotter.iren.terminate_app()
            self.plotter.add_timer_event(max_steps=max(10, int(smoke_seconds * 10) + 2),
                                         duration=100, callback=smoke)
        try:
            for _ in range(3):
                self.plotter.render_window.Render()
            self.plotter.show(title='BraVa human spatial ROIs', screenshot=str(preview),
                              interactive=show, auto_close=False)
            report['events'] = self.events
            report['gui_smoke_passed'] = bool(show and smoke_seconds and
                                             'GUI_SMOKE_AUTO_CLOSE' in self.events)
            write_json(self.run / 'spatial_ui_compatibility.json', report)
            return report
        finally:
            self.close()


def show_spatial_saved_run(run_root, *, sample_id=None, max_arrows=600, volume_opacity=.32,
                           window_size=(1800, 900), sampling_run_root=None, screenshot_path=None,
                           show=True, left_view_up=None, smoke_seconds=0):
    """Same entry contract as show_saved_run; export six views before opening UI."""
    if sampling_run_root is None or not ui.load_sampling_display_rois(Path(sampling_run_root)):
        return ui.show_saved_run(Path(run_root), sample_id=sample_id, max_arrows=max_arrows,
                                 volume_opacity=volume_opacity, window_size=window_size,
                                 sampling_run_root=sampling_run_root, screenshot_path=screenshot_path,
                                 show=show, left_view_up=left_view_up)
    samples = sorted(p for p in (Path(run_root) / 'samples').iterdir() if p.is_dir())
    if sample_id:
        samples = [p for p in samples if sample_id in {p.name, p.name.split('__', 1)[-1]}]
    if not samples:
        raise FileNotFoundError(f'No processed sample found in {run_root}')
    from .compact_static_views import _export_static_views

    options = dict(max_arrows=max_arrows, volume_opacity=volume_opacity, window_size=window_size)
    static = _export_static_views(SpatialROIViewer(samples[0], sampling_run_root, show=False, **options))
    viewer = SpatialROIViewer(samples[0], sampling_run_root, show=show, **options)
    report = viewer.run_window(show=show, screenshot_path=screenshot_path, smoke_seconds=smoke_seconds)
    report['static_views'] = static
    write_json(viewer.run / 'spatial_ui_compatibility.json', report)
    return samples[0]
