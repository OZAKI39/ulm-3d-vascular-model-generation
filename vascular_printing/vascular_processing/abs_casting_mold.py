"""One-piece casting mold using an immutable accepted vascular mesh.

All new solids are five intact walls. No port geometry, routing, cap alignment,
vascular smoothing or VascularMD reconstruction is performed here.
"""
import copy
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import trimesh
import yaml
import cadquery as cq
from . import sacrificial_fixture as old
from . import sacrificial_print_frame as previous

ROOT=Path(__file__).resolve().parents[1]
NEW_MODULES={'abs_casting_mold.py','casting_mold_qc.py','support_removal_qc.py','bambu_casting_mold.py','casting_mold_review.py'}


def load_config(path):
    cfg=yaml.safe_load(Path(path).read_text(encoding='utf-8'))
    if cfg['box']['policy']!='BOX_DIMENSIONS_REUSED':raise ValueError('BOX_RESIZING_NOT_AUTHORIZED')
    if cfg['support_removal_probe']['diameter_mm']<=0:raise ValueError('INVALID_PROBE_DIAMETER')
    return cfg


def load_inputs(cfg):
    inputs=previous.load_accepted(cfg)
    historical=old.load_config(inputs['folder']/'resolved_config.yaml')
    historical=copy.deepcopy(historical);historical['attachment_alignment']['endpoint_ids']=[]
    # This loader reads frozen companion topology only. Disabling alignment is
    # explicit; no saved routes are rebuilt or searched.
    companion=old.load_inputs(historical)
    inputs.update(companion=companion,historical_config=historical)
    for port in inputs['ports']:
        record=inputs['summary']['per_port'][port['port_id']]
        path=Path(record['port_stl'])
        if old.sha256(path)!=record['port_sha256']:raise ValueError('FROZEN_PORT_HASH_MISMATCH')
        port['reference_mesh']=trimesh.load_mesh(path);port['reference_path']=path
        endpoint=next(e for e in companion['endpoints'] if e.endpoint_id==port['port_id'])
        port['parent_branch']=endpoint.branch_id;port['swc_id']=endpoint.swc_id
    return inputs


def box_members(inner,outer):
    members={}
    for axis,name in ((0,'X'),(1,'Y')):
        for side,sign in ((0,'-'),(1,'+')):
            lo=outer[0].copy();hi=outer[1].copy()
            if side==0:hi[axis]=inner[0,axis]
            else:lo[axis]=inner[1,axis]
            members[sign+name]=np.array([lo,hi])
    lo=outer[0].copy();hi=outer[1].copy();hi[2]=inner[0,2];members['BOTTOM']=np.array([lo,hi])
    return members


def build_intact_box(inputs,cfg):
    record=inputs['old_box'];inner=np.array(record['inner_bounds_mm']);outer=np.array(record['outer_bounds_mm'])
    if not np.isclose(inner[1,2],outer[1,2]):raise ValueError('HISTORICAL_BOX_TOP_NOT_OPEN')
    def solid(bounds):return cq.Solid.makeBox(*(bounds[1]-bounds[0]),pnt=cq.Vector(*bounds[0]))
    cutter=inner.copy();cutter[1,2]+=1.
    shape=solid(outer).cut(solid(cutter)).clean()
    vertices,faces=shape.tessellate(cfg['geometry']['mesh_tolerance_mm'],cfg['geometry']['mesh_angular_tolerance_rad'])
    mesh=trimesh.Trimesh([v.toTuple() for v in vertices],faces,process=True)
    return dict(mesh=mesh,solid=shape,inner=inner,outer=outer,members=box_members(inner,outer),
        wall_thickness_mm=record['wall_thickness_mm'],bottom_thickness_mm=record['bottom_thickness_mm'],
        dimensions_policy='BOX_DIMENSIONS_REUSED',port_holes_cut=False,top_braces_added=False)


def combine(core,box,cfg):
    return old.union_core(core,[SimpleNamespace(mesh=box)],cfg)


def snapshot(inputs,out):
    result=previous.protected_snapshot(inputs,out)
    paths=[ROOT/'s1-5_sacrificial_print_frame.py']
    for folder in ('vascular_processing','config','tests'):
        paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in ('.py','.yaml','.json')
            and p.name not in NEW_MODULES and p.name not in ('abs_casting_mold_BG001.yaml','test_s1_6_abs_casting_mold.py','simple_five_wall_mold.json'))
    for name in ('point_predictions.csv','branch_predictions.csv'):
        paths.extend(inputs['base'].parent.rglob(name))
    result.update({str(p.resolve()):old.sha256(p) for p in paths if p.is_file()})
    return result
