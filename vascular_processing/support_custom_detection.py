"""Fuse mandatory evidence in named branch-arc coordinates for explicit support."""
import copy
import numpy as np
from scipy.spatial import cKDTree
from .support_inputs import save,table
from .support_detection import detect,contiguous_runs

def span_sensitive_raw(raw,span):
    """Change only the span limit; proximity to an FF point is not mandatory."""
    points=copy.deepcopy(raw)
    for p in points:
        mandatory=p.get('freefloating_mandatory',False) or 'LOCAL_MINIMUM' in p['detectors_triggered']
        p['risk_level']='CRITICAL' if p['unsupported_span']>span or mandatory else 'ATTENTION'
    return points

def intervals(raw,cfg):
    rows=[]
    for branch in sorted({p['branch_id'] for p in raw}):
        points=sorted([p for p in raw if p['branch_id']==branch and p['risk_level']=='CRITICAL'],key=lambda p:p['arc_position'])
        groups=[];current=[]
        for p in points:
            if current and p['arc_position']-current[-1]['arc_position']>cfg['sample_step_mm']*1.6:groups.append(current);current=[]
            current.append(p)
        if current:groups.append(current)
        for group in groups:
            s=[p['arc_position'] for p in group]
            rows.append(dict(interval_id=f'CRIT_{len(rows):03}',branch_id=branch,arc_start_mm=min(s),arc_end_mm=max(s),length_mm=max(s)-min(s),
                mean_radius_mm=float(np.mean([p['local_radius'] for p in group])),mean_angle=float(np.mean([p.get('gravity_angle_deg',0) for p in group])),
                curvature=float(np.mean([p['local_curvature'] for p in group])),detectors_triggered=sorted({v for p in group for v in p['detectors_triggered']}),severity='CRITICAL',point_ids=[p['point_id'] for p in group]))
    return rows

def fused_detect(mesh,box,branches,topology,inner,profile,cfg,out):
    detection,surface,layer,floating=detect(mesh,branches,topology,profile,cfg,out,box)
    raw=detection['raw'];tree=cKDTree([p['xyz'] for p in raw]);dist,idx=tree.query(floating) if len(floating) else ([],[]);ffrows=[]
    for p in raw:p['freefloating_mandatory']=False
    for j,(xyz,d,i) in enumerate(zip(floating,dist,idx)):
        p=raw[int(i)];accepted=d<=max(.6,p['local_radius']*1.5) and 'LAYER_UNSUPPORTED' in p['detectors_triggered']
        if accepted:
            p['risk_level']='CRITICAL'
            p['freefloating_mandatory']=True
            if 'FREEFLOATING' not in p['detectors_triggered']:p['detectors_triggered'].append('FREEFLOATING')
        ffrows.append(dict(raw_id=j,xyz=xyz.tolist(),nearest_point_id=p['point_id'],distance_mm=float(d),mapped_to_mandatory_interval=bool(accepted)))
    lookup={b['branch_id']:b for b in branches}
    starts=np.array([p['xyz'] for p in raw])+[0,0,.01];ends=starts.copy();ends[:,2]=box.bounds[1,2]+3
    hits,ray,_=mesh.ray.intersects_location(ends,(starts-ends)/np.linalg.norm(starts-ends,axis=1)[:,None],multiple_hits=False)
    # It is expected that the vessel above its own support blocks direct top
    # visibility. This candidate attribute is separate from whole-tree removal.
    visible=np.ones(len(raw),bool)
    for hit,i in zip(hits,ray):visible[i]=hit[2]<=starts[i,2]+.2
    for i,p in enumerate(raw):
        branch=lookup[p['branch_id']];p['tip_id']=p['point_id'];p['vertical_height_mm']=p['xyz'][2]-inner[0,2]
        p['removal_accessibility']=int(visible[i]);p['gravity_angle_deg']=float(np.degrees(np.arcsin(abs(branch['tangent'][p['sample_index'],2]))))
        p['distance_to_previous_supported_region_mm']=p['unsupported_span'] # Replaced below using actual supported samples.
        safe=np.flatnonzero(~np.asarray(branch['needs_support']))
        p['distance_to_previous_supported_region_mm']=float(np.min(abs(branch['arc'][safe]-p['arc_position']))) if len(safe) else None
    critical=intervals(raw,cfg);detection.update(critical_intervals=critical,critical_interval_count=len(critical),freefloating_raw_count=len(floating),freefloating_mandatory_mapped=sum(x['mapped_to_mandatory_interval'] for x in ffrows))
    # Attach anatomy and radius to real triangle records; no vertex-count proxy.
    centerlines=np.vstack([b['xyz'] for b in branches]);radii=np.concatenate([b['radius'] for b in branches]);codes=np.concatenate([np.full(len(b['xyz']),j) for j,b in enumerate(branches)])
    mapped=cKDTree(centerlines).query(surface['centroids'])[1]
    np.savez_compressed(out/'support_requirement_surface.npz',**surface,branch_index=codes[mapped],local_radius=radii[mapped])
    table(out/'tables/freefloating_raw_points.csv',ffrows);table(out/'tables/raw_support_points.csv',raw);table(out/'tables/critical_support_intervals.csv',critical);save(out/'detection.json',detection)
    return detection
