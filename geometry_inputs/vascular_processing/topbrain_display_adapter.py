"""TopBrain scene data and minimal UI glue; the human renderer is read-only.

Reuse its axes, legends, ROI boxes, camera helpers, font and rotation. Native
segmentation surfaces replace SWC tubes by default, independently of SWC QC.
The only mm -> um conversion is in this adapter.
"""
from __future__ import annotations

import ast
from pathlib import Path
import time

import numpy as np

from utils.rodent_vasculature import interactive as ui
from utils.sampling.sampling_types import ROIRecord
from .topbrain_qc import TopBrainError,sha256,write_json

KEY_BINDINGS = {"R/S":"selected ROI boxes","A":"all ROI boxes","C":"next anatomical group of boxes",
                "Left/Right":"previous/next ROI","PageUp/PageDown":"previous/next case (precompute outside callback)",
                "F12":"save screenshot","left_click":"select ROI box","mouse":"existing PyVista trackball rotation/zoom/pan"}


def display_record(roi,mesh,case,rank):
    # mesh has already been converted to um by the viewer. Do not scale it twice.
    bounds = np.asarray(mesh.bounds).reshape(3,2)
    low,high = bounds[:,0],bounds[:,1]
    empty = np.empty((0,3),dtype=float)
    graph = roi.graph
    points = np.array([d["coords"][:3] for _,d in graph.nodes(data=True)])*1000. if graph is not None else empty
    radii = np.array([d["coords"][3] for _,d in graph.nodes(data=True)])*1000. if graph is not None else np.empty(0)
    # Geometry here is an envelope contract for the shared box/axis helpers.
    # Native mask surface is added separately, so no synthetic SWC is required.
    return ROIRecord(roi.name,case.paths.case_id,case.paths.case_id,-1,tuple((low+high)/2),tuple(low),tuple(high),
        tuple((low+high)/2),tuple(high-low),(),(),np.arange(len(points)),np.arange(len(points)),points,radii,
        np.empty((0,2),dtype=int),np.empty(0,dtype=int),np.empty(0,dtype=int),np.empty((0,2,3)),np.empty((0,2)),
        (),(),(),roi.qc.get("component_count",1),roi.qc.get("centerline_length_mm",0)*1000.,
        roi.qc.get("centerline_length_mm",0)*1000.,cluster_id=rank,is_representative=True,selection_rank=rank+1)


def compatibility_report():
    root = Path(__file__).resolve().parents[1]
    source = root/"s1-2_swc_roi_generate_human.py"
    tree = ast.parse(source.read_text())
    calls = [n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=="show_saved_run"]
    return {"reference":str(source),"reference_sha256":sha256(source),"human_display_calls":len(calls),
        "viewport_layout":[1,2],"window_size":[1800,900],"border_color":"#606060","background":"black",
        "font":ui.UI_FONT_FAMILY,"camera_helper":"view_isometric / reset_camera; original axis and horizontal rotation helpers",
        "left_view_up":[0.,1.,0.],"right_view_up":list(ui.HORIZONTAL_ROTATION_VIEW_UP),
        "axis_transform":"NIfTI affine in data layer; mm to um once in topbrain_display_adapter; no flips",
        "projection":"same PyVista default perspective as human renderer","key_bindings":KEY_BINDINGS,
        "shared_render_functions":["_add_full_scene","_add_sampling_roi_scene","_add_sampling_boxes",
            "_add_sampling_active_outline","_add_diameter_colorbar","_remove_diameter_colorbar",
            "_install_synchronized_horizontal_rotation"],
        "differences":["native unsmoothed segmentation surface replaces SWC tube by default",
            "status text identifies case, modality, exact/derived/coarse and QC",
            "case and ROI navigation added with unused keys; screenshot added on F12"],
        "renderer_sha256":sha256(Path(ui.__file__)),"renderer_modified":False}


class TopBrainViewer:
    def __init__(self,case,rois,run,manifest,*,initial_roi=None,show=True,show_native_labels=False):
        import pyvista as pv
        self.case,self.run,self.manifest = case,Path(run),manifest
        self.rois = [r for r in rois if r.mask is not None and r.mask.any()]
        if not self.rois:
            raise TopBrainError("EMPTY_ROI","No native or derived ROI mask is available for display")
        if initial_roi and initial_roi not in [r.name for r in self.rois]:
            raise TopBrainError("ROI_NOT_AVAILABLE",f"Requested {initial_roi}; see per-ROI QC for reasons")
        self.surfaces = []
        for roi in self.rois:
            mesh = pv.read(manifest["rois"][roi.name]["outputs"]["mask_surface"])
            mesh.points *= 1000.
            mesh.point_data[ui.DIAMETER_SCALAR_NAME] = mesh.point_data["display_diameter_estimate_mm"]*1000.
            self.surfaces.append(mesh)
        self.records = tuple(display_record(r,m,case,i) for i,(r,m) in enumerate(zip(self.rois,self.surfaces)))
        self.plotter = pv.Plotter(shape=(1,2),border=True,border_color="#606060",off_screen=not show,window_size=(1800,900))
        self.plotter.theme.font.family = ui.UI_FONT_FAMILY
        self.controllers = [None,None]
        self.overlay_actors = []
        self.actor_to_index = {}
        self.group_index = 0
        self.index = next((i for i,r in enumerate(self.rois) if r.name==initial_roi),
                          next((i for i,r in enumerate(self.rois) if r.name=="RMCA_M2M3"),
                               next((i for i,r in enumerate(self.rois) if r.territory=="MCA"),0)))
        self.events = []
        self.case_delta = 0
        self.show_native_labels = show_native_labels
        self.native_surfaces = {}
        if show_native_labels:
            from .topbrain_geometry import mask_surface
            for roi in self.rois:
                self.native_surfaces[roi.name] = []
                for name in roi.native_labels:
                    mask = case.mask(name) & roi.mask
                    if mask.any():
                        mesh = mask_surface(mask,case)
                        mesh.points *= 1000.
                        self.native_surfaces[roi.name].append(mesh)
        self.plotter.subplot(0,0)
        context = pv.read(manifest["context_surface"])
        context.points *= 1000.
        self.plotter.add_mesh(context,color="white",opacity=.32,smooth_shading=True,pickable=False)
        if show_native_labels:
            from .topbrain_geometry import mask_surface
            for side in ("R","L"):
                if case.label_map.has(side+"-M1") and case.mask(side+"-M1").any():
                    mesh = mask_surface(case.mask(side+"-M1"),case)
                    mesh.points *= 1000.
                    self.plotter.add_mesh(mesh,style="wireframe",color="#51E5FF",opacity=.35,pickable=False)
        corners = np.asarray(context.bounds).reshape(3,2).T
        geometry = ui.InteractiveSceneGeometry((),np.empty((0,3)),np.empty((0,3)),corners,("context_bounds","context_bounds"))
        _,self.controllers[0] = ui._add_full_scene(self.plotter,None,geometry,spacing_xyz_um=tuple(case.spacing_mm*1000.),
             volume_opacity=.32,sample_id=f"TopBrain {case.paths.case_id} {case.paths.modality.upper()}A",sampling_available=True,left_view_up=(0.,1.,0.))
        self.plotter.add_text(f"TopBrain | Case {case.paths.case_id} | {case.paths.modality.upper()}A\nMeVO anatomical ROI",
             position="upper_left",font_size=11,font=ui.UI_FONT_FAMILY,color="white",name="full_scene_title")
        self.redraw_boxes(list(range(len(self.rois))),"selected representatives")
        self.select(self.index)
        self.plotter.subplot(0,0)
        actions = {"a":lambda:self.redraw_boxes(list(range(len(self.rois))),"all candidates"),
                   "r":lambda:self.redraw_boxes(list(range(len(self.rois))),"selected representatives"),
                   "s":lambda:self.redraw_boxes(list(range(len(self.rois))),"selected representatives"),
                   "c":self.next_group,"Left":lambda:self.select((self.index-1)%len(self.rois)),
                   "Right":lambda:self.select((self.index+1)%len(self.rois)),
                   "PageUp":lambda:self.request_case(-1),"PageDown":lambda:self.request_case(1),"F12":self.screenshot}
        for key,callback in actions.items():
            for variant in ((key,key.upper()) if len(key)==1 else (key,)):
                self.plotter.clear_events_for_key(variant)
                self.plotter.add_key_event(variant,callback)
        self.plotter.enable_mesh_picking(callback=self.pick,show=False,show_message=False,left_clicking=True,use_actor=True)
        self.rotation = None
        if show:
            cameras = tuple(renderer.camera for renderer in self.plotter.renderers)
            self.rotation = ui._install_synchronized_horizontal_rotation(self.plotter,cameras,left_view_up=(0.,1.,0.))

    def redraw_boxes(self,indices,mode):
        self.plotter.subplot(0,0)
        for actor in self.overlay_actors:
            self.plotter.remove_actor(actor,reset_camera=False,render=False)
        self.overlay_actors,self.actor_to_index = ui._add_sampling_boxes(self.plotter,self.records,indices,mode_label=mode)
        self.plotter.render()

    def next_group(self):
        index = self.group_index % len(self.rois)
        self.group_index += 1
        self.redraw_boxes([index],f"cluster {index}: {self.rois[index].name}")
        self.events.append("C:changed_box_group")

    def pick(self,actor):
        index = self.actor_to_index.get(getattr(actor,"memory_address",""))
        if index is not None:
            self.select(index)

    def select(self,index):
        self.index = index
        roi,mesh,record = self.rois[index],self.surfaces[index],self.records[index]
        self.plotter.subplot(0,0)
        ui._add_sampling_active_outline(self.plotter,record)
        self.plotter.subplot(0,1)
        if self.controllers[1] is not None:
            self.controllers[1].dispose()
        ui._remove_diameter_colorbar(self.plotter)
        self.plotter.renderer.clear_actors()
        self.controllers[1] = ui._add_sampling_roi_scene(self.plotter,record,add_orientation_axes=self.controllers[1] is None)
        actor = self.plotter.add_mesh(mesh,scalars=ui.DIAMETER_SCALAR_NAME,cmap=ui.DIAMETER_COLORMAP,
             opacity=.88,smooth_shading=True,show_scalar_bar=False,pickable=False)
        ui._add_diameter_colorbar(self.plotter,actor.mapper)
        text = f"Case: {self.case.paths.case_id} | {self.case.paths.modality.upper()}A | {roi.name}\n{roi.anatomical_status}\nQC: {roi.status} | voxels {int(roi.mask.sum())} | components {roi.qc.get('component_count')}"
        text += f"\nNative mask surface | nodes {roi.qc.get('node_count','not computed')} | branches {roi.qc.get('branch_count','not computed')}"
        if roi.territory=="PCA":
            text += "\nP3/P4 unresolved"
        if self.show_native_labels:
            text += "\nNative labels: "+", ".join(roi.native_labels)
            # Per-label boundaries overlaid as fine wires, preserving the primary colour scale.
            for native_mesh in getattr(self,"native_surfaces",{}).get(roi.name,[]):
                self.plotter.add_mesh(native_mesh,style="wireframe",color="#51E5FF",opacity=.35,pickable=False)
        self.plotter.add_text(text,position="upper_left",font_size=10,font=ui.UI_FONT_FAMILY,color="white",name="sampling_roi_information")
        self.plotter.subplot(0,0)
        self.plotter.render()
        self.events.append(f"ROI:{roi.name}")

    def request_case(self,delta):
        self.case_delta = delta
        self.events.append(f"case_request:{delta}")
        # Keep the render window alive until VTK returns from its callback. Closing
        # it here can leave another active timer referencing a deleted C++ object.
        self.plotter.iren.terminate_app()

    def screenshot(self):
        path = self.run/f"topbrain_screenshot_{time.time_ns()}.png"
        self.plotter.screenshot(path)
        self.events.append(f"screenshot:{path}")

    def run_window(self,*,show=True,smoke_seconds=0):
        preview = self.run/"topbrain_preview.png"
        actual = {"viewport_count":len(self.plotter.renderers),"shape":list(self.plotter.shape),
                  "window_size":list(self.plotter.window_size),
                  "backgrounds":[list(r.GetBackground()) for r in self.plotter.renderers],
                  "left_camera_view_up":list(self.plotter.renderers[0].camera.GetViewUp()),
                  "interaction_style":self.plotter.iren.interactor.GetInteractorStyle().GetClassName(),
                  "registered_keys":sorted(self.plotter.iren._key_press_event_callbacks),
                  "surface_units":"um","surface_sources":"unsmoothed native masks"}
        if show and smoke_seconds:
            def smoke(step):
                if step==4:
                    self.select((self.index+1)%len(self.rois))
                    self.next_group()
                if step >= max(8,int(smoke_seconds*10)):
                    self.plotter.screenshot(preview)
                    self.events.append("GUI_SMOKE_AUTO_CLOSE")
                    self.plotter.iren.terminate_app()
            self.plotter.add_timer_event(max_steps=max(10,int(smoke_seconds*10)+2),duration=100,callback=smoke)
        try:
            self.plotter.show(title=f"TopBrain MeVO anatomy — {self.case.paths.case_id}",screenshot=str(preview),auto_close=True)
        finally:
            for controller in self.controllers:
                if controller is not None:
                    controller.dispose()
        report = {**compatibility_report(),"case_id":self.case.paths.case_id,"modality":self.case.paths.modality,
                  "status":"GUI_STARTED" if show else "OFFSCREEN_RENDERED","preview":str(preview),
                  "events":self.events,"case_delta":self.case_delta,"roi_count":len(self.rois),
                  "actual_runtime":actual,
                  "core_recomputed_in_callback":False}
        write_json(self.run/"topbrain_ui_compatibility.json",report)
        return report
