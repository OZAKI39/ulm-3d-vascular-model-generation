"""Compact upper caps for the frozen s1-10 mold. No vascular/mold edits."""
from __future__ import annotations
import json,hashlib
from pathlib import Path
import numpy as np
import trimesh,cadquery as cq
from .support_inputs import ROOT,discover,digest,save,verify
from .all_abs_pipeline import default_output as baseline_output
from .all_abs_geometry import build
from .casting_fixture_geometry import shape_mesh,bounds_box
from .casting_fixture_v2_geometry import world,port_transform,fuse


def load_frozen():
    base=baseline_output();s=json.loads((base/'all_abs_sacrificial_casting_summary.json').read_text())
    if s['status']!='READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED':raise ValueError('ACCEPTED_S1_10_BASELINE_REQUIRED')
    cfg=json.loads((base/'reports/config.json').read_text());inputs=discover(source=cfg.get('source'));verify(json.loads((base/'source_protection_start.json').read_text()))
    l=json.loads((base/'assembly/layout.json').read_text())
    for key in ['low','high','transform','bifurcations']:l[key]=np.asarray(l[key])
    for p in l['ports']:
        for key in ['axis','straight_start','end','attachment','trim_origin','trim_normal','approved_contact_start','approved_contact_end']:p[key]=np.asarray(p[key])
    l['core']=trimesh.load_mesh(inputs['core']);l['core'].apply_transform(l['transform']);l['central']=trimesh.load_mesh(inputs['vascular']);l['central'].apply_transform(l['transform']);l['source_ports']={}
    for p in l['ports']:
        k=p['port_id'];m=trimesh.load_mesh(Path(inputs['candidate'])/f'manufacturing/{k}.stl');m.apply_transform(l['transform']);l['source_ports'][k]=m
    line=float(s['bambu']['mold_auto_support_off_slice']['actual_profile']['line_width']);g=build(l,cfg,line)
    # Read the actual accepted CAD and validate the reconstructed read-only
    # fixture kernel against its already-exported STL, not guessed dimensions.
    frozen=cq.importers.importStep(str(base/'outer_mold/ABS_sacrificial_outer_mold.step')).val();oldmesh=trimesh.load_mesh(base/'outer_mold/ABS_sacrificial_outer_mold.stl');calc=shape_mesh(g['mold']);error=abs(calc.volume-oldmesh.volume)
    if error>.001 or not np.allclose(calc.bounds,oldmesh.bounds,atol=5e-6,rtol=0):raise ValueError('FROZEN_MOLD_RECONSTRUCTION_MISMATCH')
    g['mold']=frozen;g['mold_mesh']=oldmesh
    for k in g['retainers']:g['retainers'][k]=cq.importers.importStep(str(base/f'retainers/port_{k}_ABS_retainer.step')).val()
    return base,s,inputs,cfg,l,g,line


def validate_parameters(angle,clearance,thickness,axial):
    if not 0<angle<=180:raise ValueError('C_RETAINER_CLOSED_OR_EXCESSIVE_WRAP')
    if clearance<=0:raise ValueError('C_RETAINER_PRESS_FIT_FORBIDDEN')
    if thickness<=0 or axial<=0:raise ValueError('C_RETAINER_NONPOSITIVE_DIMENSION')


def local_cap(p,oldg,cfg,line,angle):
    k=p['port_id'];old=next(d for d in oldg['design'] if d['port_id']==k);clearance=p['anti_float_gap_mm'];guide=old['retainer_guide_clearance_mm']
    if cfg.get('port_clearance_override_mm') is not None:clearance=float(cfg['port_clearance_override_mm'])
    if cfg.get('guide_clearance_override_mm') is not None:guide=float(cfg['guide_clearance_override_mm'])
    t=line*cfg['retainer_wall_thickness_lines'];L=cfg['retainer_axial_length_mm'];validate_parameters(angle,clearance,t,L)
    if guide<=0:raise ValueError('C_RETAINER_GUIDE_PRESS_FIT_FORBIDDEN')
    Ri=p['radius_mm']+clearance;Ro=Ri+t;theta=np.radians((180-angle)/2);xo,zo=Ro*np.cos(theta),Ro*np.sin(theta);xi,zi=Ri*np.cos(theta),Ri*np.sin(theta);front=cfg['retainer_cap_front_from_trim_mm'];front=oldg['local_parts'][k]['retainer'].BoundingBox().ymin if front is None else float(front);back=front+L
    arc=cq.Workplane('XZ').moveTo(xo,zo).threePointArc((0,Ro),(-xo,zo)).lineTo(-xi,zi).threePointArc((0,Ri),(xi,zi)).close().extrude(L).val().translate((0,back,0))
    # A narrow neck initially lies inside the top of the arc, then widens and
    # rises at <=45 degrees when printed along the extension axis. The rear
    # tongue fits the EXISTING guide and stops at the existing inner wall.
    inner=oldg['local_parts'][k]['inner_guide_x'];half=inner-guide;ys=back-cfg['retainer_neck_overlap_mm'];ye=ys+cfg['retainer_neck_ramp_length_mm'];w=p['wall_distance_from_trim_mm'];small=cfg['retainer_neck_half_width_mm'];old_roof=p['radius_mm']+p['anti_float_gap_mm']+.7+3*line
    if ye>=w:raise ValueError('C_RETAINER_GUIDE_ZONE_TOO_SHORT_'+k)
    def rectangle(y,x,z0,z1):return cq.Wire.makePolygon([cq.Vector(*q) for q in [(-x,y,z0),(x,y,z0),(x,y,z1),(-x,y,z1)]],close=True)
    neck=cq.Solid.makeLoft([rectangle(ys,small,Ri+.05,Ro-.05),rectangle(ye,half,old_roof-3*line,old_roof)],ruled=True)
    tongue=bounds_box([-half,ye,old_roof-3*line],[half,w,old_roof]);ret=fuse([arc,neck,tongue]);bb=ret.BoundingBox();overhang=max((half-small)/(ye-ys),abs(old_roof-(Ro-.05))/(ye-ys),abs(old_roof-3*line-(Ri+.05))/(ye-ys))
    parameters=dict(port_id=k,role=p['role'],c_retainer_wrap_angle_deg=float(angle),retainer_opening_angle_deg=360.-angle,upper_clearance_mm=clearance,retainer_wall_thickness_lines=cfg['retainer_wall_thickness_lines'],retainer_wall_thickness_mm=t,retainer_axial_length_mm=L,total_axial_extent_mm=w-front,retainer_guide_length_mm=w-ye,retainer_max_height_above_port_mm=bb.zmax,retainer_max_height_above_port_centerline_mm=bb.zmax,retainer_volume_mm3=ret.Volume(),guide_clearance_mm=guide,cap_front_y_mm=front,cap_back_y_mm=back,neck_start_y_mm=ys,neck_end_y_mm=ye,guide_half_width_mm=half,guide_bottom_z_mm=old_roof-3*line,guide_top_z_mm=old_roof,neck_max_growth_per_print_layer=overhang,retainer_wallward_slide_mm=cfg['retainer_wallward_slide_mm'],wall_datum_contact_area_mm2=2*half*3*line,no_exposed_vertical_stem=True,no_elastic_snap_required=True,upper_cap_only=bb.zmin>0,below_nominal_fill=True)
    return ret,parameters


def candidate(l,g,cfg,line,angle):
    new=dict(g);new['retainers']={};new['compact_parameters']=[];new['compact_local']={};new['compact_ports']=l['ports']
    for p in l['ports']:
        local,record=local_cap(p,g,cfg,line,angle);k=p['port_id'];new['compact_local'][k]=local;new['retainers'][k]=world(local,p);record['below_nominal_fill']=record['retainer_max_height_above_port_mm']+p['trim_origin'][2]<l['high'][2];new['compact_parameters'].append(record)
    return new
