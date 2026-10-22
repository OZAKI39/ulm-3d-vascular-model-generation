"""Read-only compact ROI data for the BraVa dual-viewport UI."""
from pathlib import Path
from types import SimpleNamespace
import csv
import json

import numpy as np

from utils.rodent_vasculature import interactive as ui
from utils.sampling.sampling_types import CutPort, ROIRecord
from .brava_display_adapter import BraVaViewer
from .compact_annotations import (DualPanelAnnotations, display_font, ORIENTATION_LEGEND_LABEL, SCALE_TEXT_GAP_PX,
    CANDIDATE_LEGEND_LABEL, LEGEND_CENTRE_Y, SCALE_TITLE_FONT_SIZE, SCALE_LABEL_FONT_SIZE,
    COMPACT_LEGEND_FONT_SIZE, COMPACT_AXIS_TITLE_FONT_SIZE)
from .swc_export import read_source
from .topbrain_qc import sha256, write_json

MODES = ('MINI', 'BALANCED', 'RICH')
CONTEXT_COLORS = {0: '#B0B0B0', 1: '#FF9F1C', 2: '#35D6E3'}
GLOBAL_ARROW_SCALE = 1.6  # 20% smaller than the preceding 2x display scale.


def compact_record(mode, item, original, rank):
    """Keep compensated radii, native coordinates and exact full-tree edge IDs."""
    export = item['exports']['compensated']; path = Path(export['path'])
    if sha256(path) != export['sha256']:
        raise ValueError('COMPACT_SWC_CHANGED: ' + str(path))
    graph = read_source(path).graph
    with path.with_name(path.stem + '_mapping.csv').open() as stream:
        mapping = {int(r['swc_id']): int(r['original_swc_id']) for r in csv.DictReader(stream)}
    if set(mapping) != set(graph) or set(mapping.values()) != set(item['selected_original_ids']):
        raise ValueError('COMPACT_MAPPING_MISMATCH')
    ids = list(graph); index = {n:i for i,n in enumerate(ids)}
    source_ids = [mapping[n] for n in ids]
    xyzr = np.array([graph.nodes[n]['coords'] for n in ids])
    np.testing.assert_array_equal(xyzr[:, :3], np.array([original.nodes[n]['coords'][:3] for n in source_ids]))
    if any(graph.nodes[n]['swc_type'] != original.nodes[mapping[n]]['swc_type'] for n in ids):
        raise ValueError('COMPACT_TYPE_MISMATCH')
    global_edges = {edge:i for i,edge in enumerate(original.edges)}
    retained = {(mapping[a], mapping[b]) for a,b in graph.edges}
    if not retained <= set(global_edges):
        raise ValueError('COMPACT_NON_NATIVE_EDGE')
    points = xyzr[:, :3] * 1000.; radii = xyzr[:, 3] * 1000.
    edges = np.array([[index[a], index[b]] for a,b in graph.edges], dtype=np.int64)
    edge_ids = np.array([global_edges[mapping[a], mapping[b]] for a,b in graph.edges], dtype=np.int64)
    padding = max(float(radii.max()) * 1.2, 1.)
    low, high = points.min(axis=0)-padding, points.max(axis=0)+padding
    root = next(n for n in graph if graph.in_degree(n)==0)
    terminals = tuple(index[n] for n in ids if original.out_degree(mapping[n])==0)
    ports = []
    for n in ids:
        old = mapping[n]
        for edge in list(original.in_edges(old)) + list(original.out_edges(old)):
            if edge not in retained:
                ports.append(CutPort(f'{mode}_{len(ports)}', index[n], global_edges[edge], tuple(points[index[n]]),
                    float(radii[index[n]]), 'native_sample_boundary', 'MANUFACTURING_SUBTREE_BOUNDARY'))
    length = float(np.linalg.norm(points[edges[:, 1]]-points[edges[:, 0]], axis=1).sum())
    return ROIRecord(roi_id=mode, source_model_id='BG001', source_mouse_id='BG001', anchor_id=mapping[root],
        anchor_position_um=tuple(points[index[root]]), bbox_min_um=tuple(low), bbox_max_um=tuple(high),
        bbox_center_um=tuple((low+high)/2), bbox_size_um=tuple(high-low), global_node_ids=tuple(source_ids),
        global_edge_ids=tuple(map(int, edge_ids)), local_node_ids=np.arange(len(ids)),
        local_node_global_ids=np.array(source_ids), local_node_positions_um=points, local_node_radius_um=radii,
        local_edges=edges, local_edge_ids=np.arange(len(edges)), local_edge_global_ids=edge_ids,
        local_edge_points_um=points[edges], local_edge_radius_um=radii[edges], true_terminal_local_ids=terminals,
        true_terminal_global_ids=tuple(source_ids[i] for i in terminals), cut_ports=tuple(ports),
        raw_component_count=1, raw_total_vessel_length_um=length, retained_component_length_um=length,
        structural_features={'branch_count':item['stats']['branch_count'], 'bifurcation_count':item['stats']['bifurcation_count']},
        radius_features={f'r{q}':float(np.percentile(radii,q)) for q in (10,25,50,75,90)},
        cluster_id=rank, is_representative=True, selection_rank=rank+1)


def load_compact(directory):
    """Validate cached inputs; no manufacturing, fitting or semantic computation."""
    directory = Path(directory).resolve(); path = directory/'compact_manifest.json'
    if not path.is_file():
        raise ValueError('COMPACT_RESULTS_MISSING: ' + str(path))
    manifest = json.loads(path.read_text())
    artifact = json.loads((directory/'artifact_manifest.json').read_text())
    needed = ['compact_manifest.json']
    for mode in MODES:
        needed += [f'candidates/{mode}/compensated.swc', f'candidates/{mode}/compensated_mapping.csv',
                   f'candidates/{mode}/surface.vtk']
    for relative in needed:
        if sha256(directory/relative) != artifact['files'].get(relative):
            raise ValueError('COMPACT_ARTIFACT_CHANGED: ' + relative)
    side_dir = directory.parent
    report = json.loads((side_dir/'roi_manifest.json').read_text())
    raw_path = Path(report['source']['raw_source'])
    if sha256(raw_path) != report['source']['raw_sha256']:
        raise ValueError('BRAVA_SOURCE_CHANGED')
    original = read_source(raw_path).graph
    records = [compact_record(mode, manifest['candidates'][mode], original, rank) for rank,mode in enumerate(MODES)]
    return manifest, original, records


def full_context(original, side_dir):
    """Every raw BG001 node/edge remains present, with cached colours where known."""
    import pyvista as pv
    ids = list(original); index = {n:i for i,n in enumerate(ids)}; edges = list(original.edges)
    segments = np.array([[index[a],index[b]] for a,b in edges], dtype=np.int64)
    mesh = pv.PolyData(np.array([original.nodes[n]['coords'][:3] for n in ids])*1000.,
                       lines=np.c_[np.full(len(segments),2), segments].ravel())
    mesh.point_data['original_swc_id'] = np.array(ids)
    mesh.cell_data['global_edge_id'] = np.arange(len(edges))
    labels = np.zeros(len(edges), dtype=np.int8); edge_index = {e:i for i,e in enumerate(edges)}
    # Import saved semantic colours, while retaining all arteries outside the side.
    for path in sorted(Path(side_dir).parent.glob('*/labeled_tree.vtp')):
        source = pv.read(path); mapping = source.point_data['original_swc_id']
        np.testing.assert_allclose(source.points, np.array([original.nodes[int(n)]['coords'][:3] for n in mapping]), rtol=0, atol=1e-10)
        for line,label in zip(source.lines.reshape(-1,3)[:,1:], source.cell_data['label']):
            edge = (int(mapping[line[0]]), int(mapping[line[1]]))
            if edge not in edge_index:
                raise ValueError('CACHED_CONTEXT_NON_NATIVE_EDGE')
            labels[edge_index[edge]] = label
    mesh.cell_data['label'] = labels
    rgb = np.array([pv.Color(CONTEXT_COLORS[int(label)]).int_rgb for label in labels], dtype=np.uint8)
    mesh.cell_data['display_rgb'] = rgb.copy()
    return mesh, rgb


class CompactBraVaViewer(DualPanelAnnotations, BraVaViewer):
    """Shared geometry/controls with compact data and display-only mm annotations."""
    output_subdirectory = 'compact_visualization'
    roi_view_up = (0., 1., 0.)

    def load_dataset(self):
        return load_compact(self.dataset_dir)

    def __init__(self, directory, *, initial_roi=None, view='strict', show=True, output_dir=None):
        import pyvista as pv
        self.dataset_dir = Path(directory).resolve()
        self.manifest, self.original, self.records = self.load_dataset()
        self.run = Path(output_dir) if output_dir else self.dataset_dir.parent/self.output_subdirectory
        if self.run.resolve().is_relative_to(self.dataset_dir):
            raise ValueError('UI output must be outside frozen compact manufacturing artifacts')
        self.run.mkdir(parents=True, exist_ok=True)
        self.view=view; self.side=self.dataset_dir.parent.name
        self.full_mesh, self.base_rgb = full_context(self.original, self.dataset_dir.parent)
        self.rois = []
        modes = [record.roi_id for record in self.records]
        for mode in modes:
            item = self.manifest['candidates'][mode]
            component = dict(strict=dict(node_count=item['stats']['node_count'], original_branch_count=item['stats']['branch_count']),
                VascularMD_status=item['model']['status'], model_outputs=dict(surface_vtk=item['surface_vtk']))
            self.rois.append(SimpleNamespace(name=mode, side=self.side, component=component))
        initial_roi = (initial_roi or self.manifest['default_print_candidate']).upper()
        if initial_roi not in modes: raise ValueError('ROI_NOT_AVAILABLE: '+initial_roi)
        self.index=modes.index(initial_roi); self.initial_roi=initial_roi
        self.events=[]; self.case_delta=0; self.group_index=0
        self.controllers=[None,None]; self.overlay_actors=[]; self.actor_to_index={}
        self.orientation_labels=[None,None]
        self.legend_layout_sizes=[None,None]
        self.context_side=None; self.context_actors=[]; self.highlight_history=[]
        self.plotter=pv.Plotter(shape=(1,2),border=True,border_color='#606060',off_screen=not show,window_size=(1800,900))
        self.plotter.theme.font.family=ui.UI_FONT_FAMILY
        self.scale_alignment_observer=self.plotter.renderers[1].AddObserver('EndEvent',self.align_bottom_scale)
        self.select(self.index)
        actions={'a':lambda:self.redraw_boxes(list(range(len(self.rois))),'all candidates'),
            'r':lambda:self.redraw_boxes(list(range(len(self.rois))),'selected representatives'),
            's':lambda:self.redraw_boxes(list(range(len(self.rois))),'selected representatives'),
            'c':self.next_group,'Left':lambda:self.select((self.index-1)%len(self.rois)),
            'Right':lambda:self.select((self.index+1)%len(self.rois)),
            'PageUp':lambda:self.select((self.index-1)%len(self.rois)),
            'PageDown':lambda:self.select((self.index+1)%len(self.rois)),'F12':self.screenshot}
        for key,callback in actions.items():
            for variant in ((key,key.upper()) if len(key)==1 else (key,)):
                self.plotter.clear_events_for_key(variant);self.plotter.add_key_event(variant,callback)
        self.plotter.enable_mesh_picking(callback=self.pick,show=False,show_message=False,left_clicking=True,use_actor=True)
        self.rotation=None
        if show:
            self.rotation=ui._install_synchronized_horizontal_rotation(self.plotter,
                tuple(r.camera for r in self.plotter.renderers),left_view_up=(0.,1.,0.),
                right_view_up=self.roi_view_up)

    def draw_context(self, side):
        self.plotter.subplot(0,0)
        if self.controllers[0] is not None:self.controllers[0].dispose()
        if self.orientation_labels[0] is not None:
            self.orientation_labels[0].dispose();self.orientation_labels[0]=None
        self.plotter.renderer.clear_actors();self.overlay_actors=[]
        self.context_side=side; mesh=self.full_mesh
        self.full_actor=self.plotter.add_mesh(mesh,scalars='display_rgb',rgb=True,preference='cell',
            line_width=2.4,opacity=.86,pickable=False,show_scalar_bar=False)
        segments=mesh.lines.reshape(-1,3)[:,1:]
        starts,ends=mesh.points[segments[:,0]],mesh.points[segments[:,1]]
        vectors=ends-starts;lengths=np.linalg.norm(vectors,axis=1);keep=lengths>1e-12
        candidates=(starts[keep]+ends[keep])/2
        take=ui._spatially_spread_indices(candidates,maximum_count=ui.RIGHT_PARENT_ARROW_MAX_COUNT)
        geometry=ui.InteractiveSceneGeometry((),candidates[take],vectors[keep][take]/lengths[keep][take,None],
            np.asarray(mesh.bounds).reshape(3,2).T,('context_bounds','context_bounds'))
        _,self.controllers[0]=ui._add_full_scene(self.plotter,None,geometry,spacing_xyz_um=(1.,1.,1.),
            volume_opacity=.32,sample_id='BG001',sampling_available=True,left_view_up=(0.,1.,0.))
        # Enlarge only the existing arrow glyphs about their original anchors,
        # after camera setup, so neither vessels nor camera framing changes.
        arrow_actors=[a for a in self.plotter.renderer.actors.values() if a.IsA('vtkActor')
            and a.GetMapper() and a.GetMapper().GetInput()
            and a.GetMapper().GetInput().GetPointData().GetArray('GlyphVector') is not None]
        if len(arrow_actors)!=1:raise ValueError('EXPECTED_ONE_GLOBAL_ARROW_ACTOR')
        self.context_arrow_origins=geometry.arrow_points_um.copy()
        mesh=arrow_actors[0].mapper.dataset
        points=mesh.points.reshape(len(self.context_arrow_origins),-1,3)
        points[:]=self.context_arrow_origins[:,None,:]+GLOBAL_ARROW_SCALE*(points-self.context_arrow_origins[:,None,:])
        mesh.GetPoints().Modified();mesh.Modified()
        self.plotter.remove_actor('full_scene_title',reset_camera=False,render=False)
        # Reserve space for the requested larger, more distant tick labels.
        # Preserve the focal point, viewing direction and all world coordinates.
        self.plotter.camera.Zoom(.90)
        self.style_annotations(0)
        self.redraw_boxes(list(range(len(self.rois))),'selected representatives')

    def redraw_boxes(self,indices,mode):
        suppressed=self.plotter.suppress_rendering
        self.plotter.suppress_rendering=True
        try:super().redraw_boxes(indices,mode)
        finally:self.plotter.suppress_rendering=suppressed
        # The shared box helper recreates this annotation on every A/R/S/C action.
        self.plotter.remove_actor('sampling_layer_mode',reset_camera=False,render=False)
        for name in ('sampling_roi_labels-points','sampling_roi_labels-labels'):
            self.plotter.remove_actor(name,reset_camera=False,render=False)
        self.plotter.render()

    def select(self,index):
        # Parent owns the existing ROI box, right renderer, axes, camera and colours.
        # Updating cell RGB in the SAME full actor prevents accidental ROI-only context.
        suppressed=self.plotter.suppress_rendering
        self.plotter.suppress_rendering=True
        try:super().select(index)
        finally:self.plotter.suppress_rendering=suppressed
        record=self.records[index]; mode=self.rois[index].name
        rgb=self.base_rgb.copy(); rgb[record.local_edge_global_ids]=[255,0,0]
        self.full_mesh.cell_data['display_rgb'][:]=rgb
        self.full_mesh.GetCellData().GetArray('display_rgb').Modified()
        self.full_mesh.Modified()
        self.highlight_history.append(dict(roi=mode, highlighted_edges=record.edge_count,
            full_context_edges=self.full_mesh.n_cells, full_context_nodes=self.full_mesh.n_points))
        self.plotter.subplot(0,1)
        # Y-up makes this elongated ROI wider on screen; retain room for its
        # physical axis labels; preserve the existing camera framing.
        self.plotter.camera.Zoom(.85)
        self.plotter.remove_actor('sampling_roi_information',reset_camera=False,render=False)
        self.style_annotations(1)
        self.plotter.subplot(0,0); self.plotter.render()

    def run_window(self, *, show=True, smoke_seconds=0):
        try:
            # VTK computes the strip's pixel rectangle on the first render.
            # Settle overlay placement before the initial screenshot/window.
            for _ in range(3):self.plotter.render_window.Render()
            report=super().run_window(show=show,smoke_seconds=smoke_seconds)
        finally:
            self.plotter.renderers[1].RemoveObserver(self.scale_alignment_observer)
            for labels in self.orientation_labels:
                if labels is not None:labels.dispose()
        report.update(source='compact-brava',data_source=str(self.dataset_dir),initial_roi=self.initial_roi,
            full_global_node_count=len(self.original),full_global_edge_count=self.original.number_of_edges(),
            highlight_history=self.highlight_history,selected_colour='#FF0000',
            exact_original_edge_mapping=True,full_global_vasculature_always_retained=True,
            print_transform_applied=False,manufacturing_data_modified=False,
            top_information_text_visible=False,legend_layout='centred_single_row_below_top',
            arrow_legend=ORIENTATION_LEGEND_LABEL,arrow_legend_geometry='same vtkArrowSource as scene glyphs',
            left_arrow_scale_relative_to_previous=.8,left_arrow_scale_from_original=GLOBAL_ARROW_SCALE,
            candidate_legend=CANDIDATE_LEGEND_LABEL,legend_background_visible=False,
            legend_centre_y=LEGEND_CENTRE_Y,diameter_text_gap_px=SCALE_TEXT_GAP_PX,
            legend_font_size=COMPACT_LEGEND_FONT_SIZE,coordinate_title_font_size=COMPACT_AXIS_TITLE_FONT_SIZE,
            diameter_title_font_size=SCALE_TITLE_FONT_SIZE,diameter_label_font_size=SCALE_LABEL_FONT_SIZE,
            diameter_colorbar_layout='bottom_horizontal',
            diameter_colorbar_alignment='colour-strip centre aligned with orientation-triad origin',
            font=display_font()[0],font_file=display_font()[1],display_units='mm',decimal_places=2,
            units='Internal geometry/scalars in um unchanged; displayed values converted to mm',
            coordinate_tick_font_size=14,coordinate_tick_label_offset_px=20.,
            coordinate_tick_groups_per_axis=1,
            orientation_letters='fixed-size continuously projected text',roi_cluster_labels_visible=False,
            right_view_up=list(self.roi_view_up),renderer_modified=True,
            differences=['ROI data now reads MINI/BALANCED/RICH; BALANCED default',
                'Complete BG001 context retained; selected exact ROI edges coloured red',
                'Lower centred legend row without background; MeVO candidate ROI label; scene-matching arrow glyph',
                'Global arrows reduced 20% from previous display; right arrows unchanged',
                'Horizontal diameter bar retained; larger title and numeric labels five pixels from strip',
                'Helvetica or Arial; physical labels in mm with two decimals; larger offset axis numbers',
                'Stable orientation letters; C/R labels hidden',
                'Both subplots use Y-up, including ROI changes and automatic horizontal rotation'])
        write_json(self.run/f'brava_{self.view}_ui_compatibility.json',report)
        return report
