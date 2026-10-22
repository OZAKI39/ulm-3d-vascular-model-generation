"""Read precomputed BraVa candidates with the existing dual-viewport renderer."""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import csv
import json

import numpy as np

from utils.rodent_vasculature import interactive as ui
from utils.sampling.sampling_types import ROIRecord, CutPort
from .swc_export import read_source
from .topbrain_display_adapter import TopBrainViewer, compatibility_report, KEY_BINDINGS
from .topbrain_qc import sha256, write_json


def strict_record(component, original, side, rank):
    """Convert the saved exact subset to the renderer's micrometre data contract."""
    path=Path(component['strict']['path'])
    if sha256(path)!=component['strict']['sha256']:
        raise ValueError('ROI_SOURCE_CHANGED: '+str(path))
    graph=read_source(path).graph
    with path.with_name(path.stem+'_node_mapping.csv').open() as stream:
        mapping={int(r['roi_node_id']):int(r['original_swc_id']) for r in csv.DictReader(stream)}
    ids=list(graph);index={n:i for i,n in enumerate(ids)}
    source_ids=[mapping[n] for n in ids]
    coords=np.array([graph.nodes[n]['coords'] for n in ids])
    np.testing.assert_array_equal(coords,np.array([original.nodes[n]['coords'] for n in source_ids]))
    points,radii=coords[:,:3]*1000.,coords[:,3]*1000.
    edges=np.array([[index[a],index[b]] for a,b in graph.edges],dtype=np.int64)
    global_edges={edge:i for i,edge in enumerate(original.edges)}
    edge_ids=np.array([global_edges[(mapping[a],mapping[b])] for a,b in graph.edges])
    padding=max(float(radii.max())*1.2,1.)
    low,high=points.min(axis=0)-padding,points.max(axis=0)+padding
    root=next(n for n in graph if graph.in_degree(n)==0)
    terminal=tuple(index[n] for n in ids if original.out_degree(mapping[n])==0)
    ports=[]
    retained={(mapping[a],mapping[b]) for a,b in graph.edges}
    for node in ids:
        old=mapping[node]
        excluded=[e for e in list(original.in_edges(old))+list(original.out_edges(old)) if e not in retained]
        for e in excluded:
            ports.append(CutPort(f'{side}_{rank}_{len(ports)}',index[node],global_edges[e],
                tuple(points[index[node]]),float(radii[index[node]]),'branch_label_boundary','SEMANTIC_BRANCH_BOUNDARY'))
    length=float(np.linalg.norm(points[edges[:,1]]-points[edges[:,0]],axis=1).sum())
    record=ROIRecord(roi_id=f'{side}_MeVO_part{component["component"]:02d}',source_model_id='BG001',source_mouse_id='BG001',
        anchor_id=mapping[root],anchor_position_um=tuple(points[index[root]]),bbox_min_um=tuple(low),bbox_max_um=tuple(high),
        bbox_center_um=tuple((low+high)/2),bbox_size_um=tuple(high-low),global_node_ids=tuple(source_ids),
        global_edge_ids=tuple(map(int,edge_ids)),local_node_ids=np.arange(len(ids)),local_node_global_ids=np.array(source_ids),
        local_node_positions_um=points,local_node_radius_um=radii,local_edges=edges,local_edge_ids=np.arange(len(edges)),
        local_edge_global_ids=edge_ids,local_edge_points_um=points[edges],local_edge_radius_um=radii[edges],
        true_terminal_local_ids=terminal,true_terminal_global_ids=tuple(source_ids[i] for i in terminal),cut_ports=tuple(ports),
        raw_component_count=1,raw_total_vessel_length_um=length,retained_component_length_um=length,
        structural_features={'branch_count':component['strict']['branch_count']},cluster_id=rank,
        is_representative=True,selection_rank=rank+1)
    return record


class BraVaViewer(TopBrainViewer):
    """Reuse box selection, picking and navigation; callbacks only read results."""
    def __init__(self,result_root,*,side=None,initial_roi=None,view='strict',show=True):
        import pyvista as pv
        self.run=Path(result_root);self.view=view
        self.rois=[];self.records=[];self.manifests={}
        self.events=[];self.index=0;self.case_delta=0;self.group_index=0
        self.controllers=[None,None];self.overlay_actors=[];self.actor_to_index={}
        self.context_side=None;self.context_actors=[]
        for name in ([side] if side else ['LMCA','RMCA']):
            path=self.run/name/'roi_manifest.json'
            if not path.exists():continue
            report=json.loads(path.read_text())
            if not report.get('roi_components'):continue
            original=read_source(Path(report['source']['raw_source'])).graph
            if sha256(Path(report['source']['raw_source']))!=report['source']['raw_sha256']:
                raise ValueError('BRAVA_SOURCE_CHANGED')
            self.manifests[name]=report
            for component in report['roi_components']:
                if view=='surface' and component['VascularMD_status']!='VASCULARMD_MODELED':continue
                record=strict_record(component,original,name,len(self.records))
                self.records.append(record)
                self.rois.append(SimpleNamespace(name=record.roi_id,side=name,component=component))
        if not self.rois:raise ValueError('NO_PRECOMPUTED_BRAVA_ROI: run tools/transfer_topbrain_to_brava.py first')
        if initial_roi:
            self.index=next((i for i,r in enumerate(self.rois) if r.name==initial_roi),-1)
            if self.index<0:raise ValueError('ROI_NOT_AVAILABLE: '+initial_roi)
        self.plotter=pv.Plotter(shape=(1,2),border=True,border_color='#606060',off_screen=not show,window_size=(1800,900))
        self.plotter.theme.font.family=ui.UI_FONT_FAMILY
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
                tuple(r.camera for r in self.plotter.renderers),left_view_up=(0.,1.,0.))

    def draw_context(self,side):
        import pyvista as pv
        self.plotter.subplot(0,0)
        if self.controllers[0] is not None:self.controllers[0].dispose()
        self.plotter.renderer.clear_actors();self.overlay_actors=[]
        mesh=pv.read(self.run/side/'labeled_tree.vtp');mesh.points*=1000.
        self.context_side=side
        colors={0:('UNKNOWN','#B0B0B0'),1:('M1','#FF9F1C'),2:('MeVO','#35D6E3')}
        for code,(label,color) in colors.items():
            cells=np.flatnonzero(mesh.cell_data['label']==code)
            if len(cells):
                self.plotter.add_mesh(mesh.extract_cells(cells),color=color,line_width=2.4,opacity=.86,pickable=False)
        segments=mesh.lines.reshape(-1,3)[:,1:]
        starts,ends=mesh.points[segments[:,0]],mesh.points[segments[:,1]]
        vectors=ends-starts;lengths=np.linalg.norm(vectors,axis=1);keep=lengths>1e-12
        candidates=(starts[keep]+ends[keep])/2
        take=ui._spatially_spread_indices(candidates,maximum_count=ui.RIGHT_PARENT_ARROW_MAX_COUNT)
        geometry=ui.InteractiveSceneGeometry((),candidates[take],vectors[keep][take]/lengths[keep][take,None],
            np.asarray(mesh.bounds).reshape(3,2).T,('context_bounds','context_bounds'))
        _,self.controllers[0]=ui._add_full_scene(self.plotter,None,geometry,spacing_xyz_um=(1.,1.,1.),
            volume_opacity=.32,sample_id=f'BG001 {side}',sampling_available=True,left_view_up=(0.,1.,0.))
        self.plotter.add_text(f'BraVa BG001 | {side}\nM1: orange | MeVO: cyan | UNKNOWN: grey\nTopBrain-informed candidate labels',
            position='upper_left',font_size=11,font=ui.UI_FONT_FAMILY,color='white',name='full_scene_title')
        self.redraw_boxes([i for i,r in enumerate(self.rois) if r.side==side],'selected representatives')

    def select(self,index):
        import pyvista as pv
        self.index=index;roi=self.rois[index];record=self.records[index]
        if self.context_side!=roi.side:self.draw_context(roi.side)
        self.plotter.subplot(0,0);ui._add_sampling_active_outline(self.plotter,record)
        self.plotter.subplot(0,1)
        first=self.controllers[1] is None
        if not first:self.controllers[1].dispose()
        ui._remove_diameter_colorbar(self.plotter);self.plotter.renderer.clear_actors()
        surface=None
        if self.view=='surface':
            surface=pv.read(roi.component['model_outputs']['surface_vtk']);surface.points*=1000.
            low,high=np.asarray(surface.bounds).reshape(3,2).T
            record=replace(record,bbox_min_um=tuple(low),bbox_max_um=tuple(high),bbox_center_um=tuple((low+high)/2),
                bbox_size_um=tuple(high-low),local_edges=np.empty((0,2),dtype=int),
                local_edge_points_um=np.empty((0,2,3)),local_edge_radius_um=np.empty((0,2)),cut_ports=(),
                true_terminal_local_ids=())
        self.controllers[1]=ui._add_sampling_roi_scene(self.plotter,record,add_orientation_axes=first,
            view_up=getattr(self,'roi_view_up',None))
        if surface is not None:self.plotter.add_mesh(surface,color='#35D6E3',smooth_shading=True,opacity=.88,pickable=False)
        detail='Strict original SWC ROI' if surface is None else 'VascularMD surface | includes marked upstream M1 context'
        self.plotter.add_text(f'BG001 | {roi.name}\nANATOMICAL_TRANSFER_CANDIDATE\n{detail}\n'
            f'Nodes {roi.component["strict"]["node_count"]} | original branches {roi.component["strict"]["original_branch_count"]}\n'
            f'{roi.component["VascularMD_status"]} | MANUFACTURING_NOT_VALIDATED',
            position='upper_left',font_size=10,font=ui.UI_FONT_FAMILY,color='white',name='sampling_roi_information')
        self.plotter.subplot(0,0);self.plotter.render();self.events.append('ROI:'+roi.name)

    def run_window(self,*,show=True,smoke_seconds=0):
        preview=self.run/f'brava_{self.view}_preview.png'
        actual=dict(viewport_count=len(self.plotter.renderers),shape=list(self.plotter.shape),
            window_size=list(self.plotter.window_size),backgrounds=[list(r.GetBackground()) for r in self.plotter.renderers],
            left_camera_view_up=list(self.plotter.renderers[0].camera.GetViewUp()),
            interaction_style=self.plotter.iren.interactor.GetInteractorStyle().GetClassName(),
            registered_keys=sorted(self.plotter.iren._key_press_event_callbacks))
        if show and smoke_seconds:
            def smoke(step):
                if step==4:
                    # Invoke the installed key handlers, including the existing box modes.
                    for key in ['Right','c','a','r','s','Left','F12']:
                        for callback in self.plotter.iren._key_press_event_callbacks[key]:callback()
                        self.events.append('KEY:'+key)
                if step>=max(8,int(smoke_seconds*10)):
                    self.plotter.screenshot(preview);self.events.append('GUI_SMOKE_AUTO_CLOSE');self.plotter.iren.terminate_app()
            self.plotter.add_timer_event(max_steps=max(10,int(smoke_seconds*10)+2),duration=100,callback=smoke)
        try:self.plotter.show(title='BraVa BG001 MeVO candidates',screenshot=str(preview),auto_close=True)
        finally:
            for controller in self.controllers:
                if controller is not None:controller.dispose()
        report=compatibility_report()
        report.update(source='topbrain-brava',status='GUI_STARTED' if show else 'OFFSCREEN_RENDERED',
            preview=str(preview),view=self.view,roi_count=len(self.rois),events=self.events,actual_runtime=actual,
            core_recomputed_in_callback=False,units='source mm to display um once; no coordinate flips',
            key_bindings={**KEY_BINDINGS,'PageUp/PageDown':'previous/next cached ROI (BG001 only)'},
            differences=['left context coloured by saved branch labels','right reads exact BraVa ROI or native VascularMD surface'],
            gui_smoke_passed=bool(show and smoke_seconds and 'GUI_SMOKE_AUTO_CLOSE' in self.events and 'KEY:Right' in self.events))
        write_json(self.run/f'brava_{self.view}_ui_compatibility.json',report)
        return report
