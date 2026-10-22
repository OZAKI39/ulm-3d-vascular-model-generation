"""Read-only discovery, coordinates, hashes and actual manufacturing profiles."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def save(path,value):
    def encode(v):
        if isinstance(v,Path):return str(v)
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,set):return sorted(v)
        raise TypeError(type(v).__name__)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,default=encode,allow_nan=False)+'\n')

def table(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()})

def transform(points,matrix):
    m=np.asarray(matrix);return np.asarray(points)@m[:3,:3].T+m[:3,3]

def discover(project=ROOT,source=None):
    project=Path(project).resolve()
    roots=[Path(source)] if source else [p.parent for p in (project/'outputs').rglob('roi_60x50_summary.json') if p.parent.name=='support_aware_roi_60x50']
    if not roots:raise FileNotFoundError('SUPPORT_AWARE_SOURCE_NOT_FOUND')
    root=max(roots,key=lambda p:(p/'roi_60x50_summary.json').stat().st_mtime_ns)
    summary=json.loads((root/'roi_60x50_summary.json').read_text());candidate=root/'candidates'/summary['selected_candidate']
    choices=[]
    for p in root.glob('casting_mold*/geometry_summary.json'):
        j=json.loads(p.read_text())
        if j.get('status')=='GEOMETRY_STL_CHECKS_PASSED' and j.get('selected_candidate')==summary['selected_candidate']:
            choices.append((p.stat().st_mtime_ns,p,j))
    if not choices:raise FileNotFoundError('CURRENT_MOLD_MANIFEST_NOT_FOUND')
    _,manifest,geo=max(choices,key=lambda v:v[0]);folder=manifest.parent
    names=geo['files'];oriented=[n for n in names if n.endswith('_oriented.stl')];box=[n for n in names if n.startswith('intact_open_top_box')]
    final=[n for n in names if n.startswith('final_') and not n.endswith('_oriented.stl')]
    if not (len(oriented)==len(box)==len(final)==1):raise ValueError('AMBIGUOUS_MANIFEST_STL')
    for n,h in names.items():
        if digest(folder/n)!=h:raise ValueError('SOURCE_MANIFEST_HASH_MISMATCH: '+n)
    model=json.loads((candidate/'model_result.json').read_text());swc=Path(model['VascularMD']['outputs']['swc'])
    if not swc.exists():swc=candidate/'VascularMD'/swc.name
    result=dict(root=root,candidate=candidate,manifest=manifest,vascular=candidate/'vascular_core.stl',core=candidate/'manufacturing/core_with_ports.stl',
        mold=folder/final[0],oriented_mold=folder/oriented[0],box=folder/box[0],smoothed_swc=swc,
        compensated_swc=candidate/'compensated.swc',branch_mapping=candidate/'compensated_edges.csv',
        node_mapping=candidate/'compensated_mapping.csv',geometry=candidate/'manufacturing/geometry.json',casting_transform=root/'roi_casting_transform.json',
        print_transform=geo['print_transform_4x4'],inner_bounds=geo['inner_bounds_mm'],outer_bounds=geo['outer_bounds_mm'],selected_candidate=summary['selected_candidate'])
    required=[v for v in result.values() if isinstance(v,Path)]
    for p in required:
        if not p.exists():raise FileNotFoundError(p)
    protected=set(p for p in required if p.is_file())
    protected.update(candidate.glob('manufacturing/*.stl'))
    protected.update(project.glob('s1-[2-7]*.py'))
    result['hashes']={str(p):digest(p) for p in sorted(protected)}
    return result

def verify(hashes):
    changed=[p for p,h in hashes.items() if not Path(p).is_file() or digest(p)!=h]
    result=dict(all_unchanged=not changed,checked_files=len(hashes),changed=changed)
    if changed:raise RuntimeError('UPSTREAM_GEOMETRY_MODIFIED: '+repr(changed))
    return result

def numeric(value,base=None):
    if isinstance(value,list):value=value[0]
    if isinstance(value,str) and value.endswith('%'):
        if base is None:raise ValueError('PERCENT_BASE_UNRESOLVED')
        return float(value[:-1])*base/100
    return float(value)

def profile_values(profiles):
    p=profiles.get('process',{});m=profiles.get('machine',{});fallback=[]
    result={}
    for key,src,name,default in [('nozzle_diameter',m,'nozzle_diameter',.2),('layer_height',p,'layer_height',.1),('first_layer_height',p,'initial_layer_print_height',.1),('line_width',p,'line_width',.22),('support_top_gap',p,'support_top_z_distance',.2),('support_xy_gap',p,'support_object_xy_distance',.4)]:
        try:result[key]=numeric(src[name],result.get('nozzle_diameter',.2))
        except (KeyError,ValueError,TypeError,IndexError):result[key]=default;fallback.append(key)
    result['status']='PROFILE_VALUE_FALLBACK' if fallback else 'ACTUAL_PROFILE_PARSED';result['fallback_fields']=fallback
    return result

def resolve_profile(inputs,out):
    from . import bambu_manufacturing as b
    frozen=inputs['root']/'bambu/evaluation_profile/frozen_discovery.json'
    d=json.loads(frozen.read_text());root=Path(d['active_preset_source']).parent
    active=d['active_gui_presets'].copy();active['process']='0.10mm Standard Vessel'
    try:d.update(b.active_profiles(root,active,'ABS'));source='CURRENT_SAVED_VESSEL_PROFILE'
    except (OSError,KeyError,ValueError):source='ARCHIVED_RESOLVED_VESSEL_PROFILE'
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for k,v in d['flattened'].items():save(out/(k+'.json'),v)
    result=profile_values(d['flattened']);result.update(source=source,process_name=active['process'],printer=d['model'],material='ABS',profile_sources=d['profile_sources'])
    save(out/'profile.json',result);return result,d

def swc_branches(path,matrix,edge_mapping=None,step=.25):
    """Topological branches, resampled only for analysis; source SWC is untouched."""
    a=np.loadtxt(path,comments='#',ndmin=2);nodes={int(r[0]):r for r in a};children={i:[] for i in nodes}
    for i,r in nodes.items():
        if int(r[6]) in nodes:children[int(r[6])].append(i)
    roots=[i for i,r in nodes.items() if int(r[6]) not in nodes];forks=[i for i,v in children.items() if len(v)>1];ends=[i for i,v in children.items() if not v]
    names={}
    if edge_mapping:
        with Path(edge_mapping).open() as f:
            names={(int(r['parent_id']),int(r['child_id'])):r['original_branch'] for r in csv.DictReader(f)}
    branches=[]
    for start in sorted(roots+forks):
        for child in sorted(children[start]):
            ids=[start,child]
            while len(children[ids[-1]])==1:ids.append(children[ids[-1]][0])
            coords=transform([nodes[i][2:5] for i in ids],matrix);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(coords,axis=0),axis=1))]
            if arc[-1]<1e-8:continue
            s=np.linspace(0,arc[-1],max(3,int(np.ceil(arc[-1]/step))+1));xyz=np.column_stack([np.interp(s,arc,coords[:,k]) for k in range(3)])
            radius=np.interp(s,arc,[nodes[i][5] for i in ids]);t=np.gradient(xyz,s,axis=0);t/=np.maximum(np.linalg.norm(t,axis=1)[:,None],1e-12)
            curvature=np.linalg.norm(np.gradient(t,s,axis=0),axis=1)
            branches.append(dict(branch_id=names.get((start,child),f'{start}_{ids[-1]}'),swc_ids=ids,xyz=xyz,arc=s,radius=radius,tangent=t,curvature=curvature,is_port=False))
    return branches,dict(bifurcations=transform([nodes[i][2:5] for i in forks],matrix).reshape(-1,3) if forks else np.empty((0,3)),bifurcation_ids=forks,inlets=roots,outlets=ends,
        endpoint_xyz=transform([nodes[i][2:5] for i in roots+ends],matrix))

def port_branches(inputs,step=.25):
    """Recover ONLY saved routing curves, validating length against frozen metadata.

    No port solid is created or written. Ambiguous curve families are resolved
    using the saved physical route length, not by a newly optimized route.
    """
    import cadquery as cq
    from types import SimpleNamespace
    from . import sacrificial_fixture as f
    from .support_aware_port_design import configuration
    cfg=configuration();geometry=json.loads(inputs['geometry'].read_text())
    rows=list(csv.DictReader((inputs['candidate']/'manufacturing/port_route_candidates.csv').open(encoding='utf-8-sig')))
    result=[]
    for port in geometry['ports']:
        a=port['accepted_attachment'];e=SimpleNamespace(position=np.array(a['calculation_origin_mm']),tangent=np.array(a['cap_normal']),radius=port['radius_mm'],cap_ring=np.array(a['aligned_ring_print_mm']))
        target=np.array(port['target']);end=np.array(port['end']);n=np.array(port['direction']);lead=f.attachment_lead(e,cfg);begin=e.position+lead*e.tangent
        choices={}
        for r in rows:
            if r['endpoint_id']!=port['port_id'] or r['face']!=port['face'] or any(abs(float(r['target_'+k])-v)>1e-5 for k,v in zip('xyz',target)):continue
            key=(r['route_type'],r['tangent_scale'],r['guide_bend_radius_mm'])
            if key in choices:continue
            if r['guide_bend_radius_mm']:
                guide=f.natural_bend_points(e,port['face'],float(r['guide_bend_radius_mm']),cfg);prewall=guide[-1]
                edge=cq.Edge.makeSpline([cq.Vector(*v) for v in guide],tangents=[cq.Vector(*e.tangent),cq.Vector(*n)],scale=False)
            else:
                prewall=target-cfg['routing']['wall_approach_mm']*n;length=np.linalg.norm(prewall-begin);scale=float(r['tangent_scale'])
                edge=cq.Edge.makeSpline([cq.Vector(*begin),cq.Vector(*prewall)],tangents=[cq.Vector(*(e.tangent*length*scale)),cq.Vector(*(n*length*scale))],parameters=[0.,1.],scale=False)
            length=lead+edge.Length()+np.linalg.norm(end-prewall)
            pts,_=edge.sample(max(4,int(np.ceil(edge.Length()/step))+1));pts=np.vstack([e.position,[p.toTuple() for p in pts],target,end])
            choices[key]=(abs(length-port['length_mm']),pts)
        if not choices:raise ValueError('FROZEN_PORT_CURVE_NOT_FOUND: '+port['port_id'])
        error,coords=min(choices.values(),key=lambda v:v[0])
        if error>1e-4:raise ValueError('FROZEN_PORT_CURVE_LENGTH_MISMATCH: '+port['port_id'])
        coords=transform(coords,inputs['print_transform']);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(coords,axis=0),axis=1))]
        s=np.linspace(0,arc[-1],max(3,int(np.ceil(arc[-1]/step))+1));xyz=np.column_stack([np.interp(s,arc,coords[:,k]) for k in range(3)])
        t=np.gradient(xyz,s,axis=0);t/=np.linalg.norm(t,axis=1)[:,None]
        result.append(dict(branch_id='PORT_'+port['port_id'],swc_ids=[port['swc_id']],parent_branch=port['parent_branch'],xyz=xyz,arc=s,radius=np.full(len(s),port['radius_mm']),
            tangent=t,curvature=np.linalg.norm(np.gradient(t,s,axis=0),axis=1),is_port=True,source='RECOVERED_FROZEN_ROUTE_METADATA',length_error_mm=error))
    return result
