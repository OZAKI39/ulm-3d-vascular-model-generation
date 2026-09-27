"""Temporary peripheral print frames around a strictly read-only accepted core.

Boxes/rails/anchors are new geometry. The old central mesh supplies its bounding
box only; all production Booleans use the accepted all-port STL without editing.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import numpy as np
import trimesh
import yaml
from . import sacrificial_fixture as legacy

ROOT=Path(__file__).resolve().parents[1]
PORTS=('I1','O1','O2','O3')


def load_config(path):
    cfg=yaml.safe_load(Path(path).read_text(encoding='utf-8'))
    for group in ('frame','support_removal_probe','geometry'):
        for key,value in cfg[group].items():
            if isinstance(value,(int,float)) and (not np.isfinite(value) or value<=0):
                raise ValueError(group+'.'+key+' must be positive')
    if cfg['pdms_keep_zone']['basis']!='frozen_central_vascular_bbox_excluding_artificial_ports':
        raise ValueError('PDMS_KEEP_ZONE_BASIS_REQUIRES_EXPLICIT_CONFIGURATION')
    if not cfg['pdms_keep_zone'].get('authorization'):
        raise ValueError('PDMS_KEEP_ZONE_BASIS_NOT_CONFIRMED')
    return cfg


def load_accepted(cfg):
    folder=(ROOT/cfg['input']['alignment_root']).resolve()
    summary_path=folder/cfg['input']['summary'];report_path=folder/cfg['input']['report']
    summary=json.loads(summary_path.read_text());report=report_path.read_text()
    core_path=Path(summary['final_stl']).resolve()
    if not core_path.is_file() or core_path.name not in report or legacy.sha256(core_path)!=summary['final_sha256']:
        raise ValueError('ACCEPTED_CORE_PROVENANCE_UNRESOLVED')
    if summary['status']!='READY_FOR_HUMAN_REVIEW' or set(summary['per_port'])!=set(PORTS):
        raise ValueError('ACCEPTED_ALL_PORT_RESULTS_INCOMPLETE')
    core=trimesh.load_mesh(core_path,process=True)
    central_path=Path(summary['source_stl']).resolve()
    if legacy.sha256(central_path)!=summary['source_sha256']:
        raise ValueError('CENTRAL_REFERENCE_HASH_MISMATCH')
    central=trimesh.load_mesh(central_path,process=True)
    old_qc_path=(folder/cfg['input']['old_box_qc']).resolve()
    old_qc=json.loads(old_qc_path.read_text());old_box=Path(old_qc['box']['stl'])
    ports=[]
    for name in PORTS:
        data=summary['per_port'][name];row=data['original_selected_route']
        ports.append(dict(port_id=name,radius_mm=data['radius_mm'],face=row['face'],
            target=np.array(data['wall_target_mm']),end=np.array(data['outside_end_mm']),
            direction=legacy.NORMALS[row['face']].copy(),accepted_attachment=data['cap']))
    return dict(folder=folder,base=folder.parent,summary=summary,summary_path=summary_path,report_path=report_path,
        core_path=core_path,core=core,central_path=central_path,central=central,ports=ports,
        old_qc_path=old_qc_path,old_box_path=old_box,old_box=old_qc['box'],
        accepted_status='ALL_PORT_ATTACHMENTS_ACCEPTED')


def keep_zone_bounds(central_bounds,cfg):
    s=cfg['pdms_keep_zone'];margin=np.array([s['margin_'+a+'_mm'] for a in 'xyz'])
    if np.any(margin<0):raise ValueError('Negative PDMS keep margin')
    return np.asarray(central_bounds)+np.array([-margin,margin])


def box_mesh(bounds):
    bounds=np.asarray(bounds,dtype=float)
    if np.any(bounds[1]<=bounds[0]):raise ValueError('Nonpositive box extent')
    transform=np.eye(4);transform[:3,3]=bounds.mean(axis=0)
    return trimesh.creation.box(bounds[1]-bounds[0],transform=transform)


def member(name,role,lo,hi,port=None):
    return dict(name=name,role=role,bounds=np.array([lo,hi],dtype=float),port_id=port)


def union_members(members,cfg):
    if not members:raise ValueError('EMPTY_FRAME')
    import cadquery as cq
    solids=[cq.Solid.makeBox(*(m['bounds'][1]-m['bounds'][0]),pnt=cq.Vector(*m['bounds'][0])) for m in members]
    solid=solids[0].fuse(*solids[1:]).clean()
    vertices,faces=solid.tessellate(cfg['geometry']['mesh_tolerance_mm'],cfg['geometry']['mesh_angular_tolerance_rad'])
    return trimesh.Trimesh([v.toTuple() for v in vertices],faces,process=True)


def port_stem_bounds(port,tolerance=1e-4):
    axis='XYZ'.index(port['face'][-1]);lo=np.minimum(port['target'],port['end']);hi=np.maximum(port['target'],port['end'])
    pad=np.full(3,port['radius_mm']+tolerance);pad[axis]=tolerance
    return np.array([lo-pad,hi+pad])


def build_frame(kind,open_axis,inputs,cfg):
    if kind not in ('OPEN_PANEL','SPARSE_FRAME') or open_axis not in (0,1):raise ValueError('Only two candidate families and opposite X/Y faces are allowed')
    s=cfg['frame'];w=s['rail_width_mm'];d=s['rail_depth_mm'];members=[]
    if kind=='OPEN_PANEL':
        outer=np.array(inputs['old_box']['outer_bounds_mm']);inner=np.array(inputs['old_box']['inner_bounds_mm'])
        lo,hi=outer
        members.append(member('bottom_panel','comparison_bottom_panel',lo,[hi[0],hi[1],inner[0,2]]))
        other=1-open_axis
        for side in (0,1):
            low=lo.copy();high=hi.copy()
            if side==0:high[other]=inner[0,other]
            else:low[other]=inner[1,other]
            members.append(member(('minus' if side==0 else 'plus')+'_'+str(other)+'_panel','comparison_side_panel',low,high))
    else:
        margin=max(s['outer_margin_mm'],s['support_access_clearance_mm'],s['vessel_clearance_mm'])
        inner=np.asarray(inputs['central'].bounds)+np.array([[-margin]*3,[margin]*3])
        outer=inner+np.array([[-w,-w,-d],[w,w,d]])
        lo,hi=outer
        for zside in (0,1):
            zlo=lo[2] if zside==0 else inner[1,2];zhi=inner[0,2] if zside==0 else hi[2]
            for axis in (0,1):
                other=1-axis
                for side in (0,1):
                    low=lo.copy();high=hi.copy();low[2]=zlo;high[2]=zhi
                    if side==0:high[other]=inner[0,other]
                    else:low[other]=inner[1,other]
                    members.append(member(f'rail_{axis}_{side}_{zside}','perimeter_rail',low,high))
        post=s['corner_post_width_mm']
        for xside in (0,1):
            for yside in (0,1):
                low=lo.copy();high=hi.copy()
                for axis,side in ((0,xside),(1,yside)):
                    if side==0:high[axis]=lo[axis]+post
                    else:low[axis]=hi[axis]-post
                members.append(member(f'corner_{xside}_{yside}','corner_post',low,high))
    if kind=='OPEN_PANEL':
        for side in (0,1):
            low=outer[0].copy();high=outer[1].copy();low[2]=outer[1,2]-d
            if side==0:high[open_axis]=inner[0,open_axis]
            else:low[open_axis]=inner[1,open_axis]
            members.append(member('open_top_edge_'+str(side),'perimeter_rail',low,high))
    # Only short perimeter fingers to EXISTING external stems. They terminate
    # at the stem and never connect to middle vascular branches.
    for port in inputs['ports']:
        axis,side=legacy.face_axis(port['face']);other=1-axis
        if kind=='OPEN_PANEL' and axis!=open_axis:continue  # Existing panel already intersects this stem.
        center=port['target'];aw=s['anchor_width_mm']
        low=center.copy();high=center.copy()
        if side==0:low[axis]=outer[0,axis];high[axis]=inner[0,axis]
        else:low[axis]=inner[1,axis];high[axis]=outer[1,axis]
        low[other]-=aw/2;high[other]+=aw/2
        bottom=center[2]-inner[0,2];top=inner[1,2]-center[2]
        if bottom<=top:low[2]=outer[0,2];high[2]=center[2]+aw/2;edge='BOTTOM'
        else:low[2]=center[2]-aw/2;high[2]=outer[1,2];edge='TOP'
        members.append(member(port['port_id']+'_anchor_'+edge,'external_port_anchor',low,high,port['port_id']))
    frame=union_members(members,cfg)
    return dict(candidate=kind,layout='OPEN_'+('X' if open_axis==0 else 'Y'),open_axis=open_axis,
        open_faces=['-X','+X'] if open_axis==0 else ['-Y','+Y'],outer_bounds=outer,inner_bounds=inner,
        members=members,mesh=frame)


def combine(core,frame,cfg):
    from types import SimpleNamespace
    # Existing verified kernel and tiny-degenerate numerical cleanup only.
    return legacy.union_core(core,[SimpleNamespace(mesh=frame)],cfg)


def protected_snapshot(inputs,out):
    files=set()
    for path in inputs['base'].rglob('*'):
        if path.is_file() and not path.is_relative_to(out):files.add(path)
    previous=inputs['folder']/'protected_before.json'
    if previous.exists():files.update(Path(p) for p in json.loads(previous.read_text()) if Path(p).is_file())
    for name in ('s1-2_swc_roi_generate_human.py','s1-3_swc_roi_generate_MeVO.py','s1-4_sacrificial_box_and_ports.py'):
        files.add(ROOT/name)
    for path in (ROOT/'third_party/vascularmd').rglob('*.py'):files.add(path)
    return {str(p.resolve()):legacy.sha256(p) for p in sorted(files)}
