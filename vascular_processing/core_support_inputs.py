"""Read s1-10's frozen vascular inputs; never generate or load casting solids."""
from __future__ import annotations
import copy,json,shutil
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import ROOT,save,digest,verify,resolve_profile,numeric,swc_branches,port_branches,transform
from .support_slicer_adapter import discover_engine


def default_output():
    return ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_core_support_production'


def discover_source():
    bases=list((ROOT/'outputs').glob('**/all_abs_sacrificial_casting_system/all_abs_sacrificial_casting_summary.json'))
    if len(bases)!=1:raise ValueError('AMBIGUOUS_S1_10_SOURCE')
    base=bases[0].parent;s=json.loads(bases[0].read_text())
    revisions=list(base.glob('compact_c_retainer_v2/mold_*/mold_revision_summary.json'))
    latest=max(revisions,key=lambda p:p.stat().st_mtime_ns).parent if revisions else base
    core=latest/'vascular_core/BG001_RMCA_frozen_ABS_core.stl'
    if digest(core)!=s['vascular']['source_hash']:raise ValueError('S1_10_FROZEN_CORE_HASH_MISMATCH')
    saved=json.loads((base/'reports/source_discovery.json').read_text())
    # Keep only vascular provenance. Casting geometry never enters this stage.
    keys=['root','candidate','smoothed_swc','compensated_swc','branch_mapping','node_mapping','geometry','casting_transform']
    inp={k:Path(saved[k]) for k in keys};inp.update(core=core,print_transform=saved['print_transform'],s1_10_root=base,latest_root=latest)
    required=[core,*[inp[k] for k in keys if k not in ['root','candidate']],inp['candidate']/'manufacturing/port_route_candidates.csv']
    inp['hashes']={str(p):digest(p) for p in required};return inp


def freeze(inputs,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);path=out/'validation/protected_before.json'
    if path.exists():hashes=json.loads(path.read_text());verify(hashes)
    else:
        hashes=dict(inputs['hashes'])
        # Protect every s1-10 artifact as bytes, without using its geometry.
        for p in inputs['s1_10_root'].rglob('*'):
            if p.is_file():hashes[str(p)]=digest(p)
        for p in [ROOT/'s1-10_all_abs_sacrificial_casting_system.py',*list((ROOT/'vascular_processing').glob('all_abs*.py')),*list((ROOT/'vascular_processing').glob('compact_c*.py')),*list((ROOT/'vascular_processing').glob('compact_mold*.py'))]:
            hashes[str(p)]=digest(p)
        save(path,hashes)
    (out/'source').mkdir(exist_ok=True)
    for key in ['core','compensated_swc','smoothed_swc','geometry','casting_transform']:
        name='frozen_vascular_core.stl' if key=='core' else inputs[key].name
        shutil.copy2(inputs[key],out/'source'/name)
    save(out/'source/discovery.json',{k:v for k,v in inputs.items() if k!='hashes'})
    return hashes


def current_profiles(inputs,out):
    out=Path(out);values,d=resolve_profile(inputs,out/'discovery')
    if values['source']!='CURRENT_SAVED_VESSEL_PROFILE':raise ValueError('CURRENT_SAVED_PROFILE_REQUIRED_NO_ARCHIVE_FALLBACK')
    actual=copy.deepcopy(d['flattened'])
    for kind,settings in actual.items():settings.update(type=kind,instantiation='true')
    before=copy.deepcopy(actual);actual['process']['enable_support']='0'
    if actual['filament'].get('filament_type')!=['ABS'] or numeric(actual['machine']['nozzle_diameter'])!=.2 or numeric(actual['process']['layer_height'])!=.1:raise ValueError('X1C_02_ABS_010_REQUIRED')
    p=actual['process'];w=numeric(p['line_width'],.2);sw=numeric(p.get('support_line_width',p['line_width']),.2)
    values.update(line_width=w,support_line_width=sw,first_layer_height=numeric(p['initial_layer_print_height']),layer_height=numeric(p['layer_height']),roof_gap_mm=numeric(p['support_top_z_distance']),source='CURRENT_SAVED_VESSEL_PROFILE',gap_status='PHYSICAL_REMOVAL_GAP_UNCALIBRATED',bridge_limit_mm=2.,bridge_calibration_status='UNCALIBRATED',initial_layer_line_width=numeric(p['initial_layer_line_width'],.2))
    for kind,settings in actual.items():save(out/(kind+'.json'),settings)
    engine=discover_engine('Bambu',d,out/'engine')
    if engine['status']!='ENGINE_FOUND':raise RuntimeError('ACTUAL_BAMBU_BINARY_REQUIRED')
    save(out/'profiles.json',dict(values=values,actual=actual,engine=engine,source=d,overrides={'process.enable_support':dict(before=before['process'].get('enable_support'),after='0',reason='Explicit user requirement; saved user profiles unchanged')},hashes={k:digest(out/(k+'.json')) for k in actual}))
    return values,actual,engine


def branches(inputs,step):
    matrix=np.asarray(json.loads(inputs['casting_transform'].read_text())['transform_4x4'])
    result,top=swc_branches(inputs['smoothed_swc'],matrix,inputs['branch_mapping'],step)
    inv=np.linalg.inv(inputs['print_transform'])
    for b in port_branches(inputs,step):
        b['xyz']=transform(b['xyz'],inv);b['tangent']=np.asarray(b['tangent'])@inv[:3,:3].T;result.append(b)
    return result,top


def place_branches(branches,matrix):
    out=copy.deepcopy(branches)
    for b in out:
        b['xyz']=transform(b['xyz'],matrix);b['tangent']=b['tangent']@matrix[:3,:3].T
    return out


def rotations(count=48):
    """Principal axes + Fibonacci sphere; two rolls per build direction."""
    if count<32 or count>64 or count%2:raise ValueError('32_TO_64_EVEN_ORIENTATIONS_REQUIRED')
    from scipy.spatial.transform import Rotation
    directions=[(0,0,1),(0,0,-1),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0)]
    n=count//2-len(directions);gold=np.pi*(3-np.sqrt(5))
    for i in range(n):
        z=1-2*(i+.5)/n;r=np.sqrt(1-z*z);directions.append((r*np.cos(i*gold),r*np.sin(i*gold),z))
    return [Rotation.from_euler('z',roll,degrees=True).as_matrix()@trimesh.geometry.align_vectors(d,[0,0,1])[:3,:3] for d in directions for roll in [0,90]]


def search(mesh,bs,profile,out,count=48):
    from .support_rfus_geometry import lower_points
    rows=[];w=profile['support_line_width'];xyz=np.vstack([b['xyz'][::4] for b in bs]);arc=np.concatenate([np.r_[np.linalg.norm(np.diff(b['xyz'][::4],axis=0),axis=1),0] for b in bs])
    for i,R in enumerate(rotations(count)):
        M=np.eye(4);M[:3,:3]=R;v=mesh.vertices@R.T;lo=v.min(0);hi=v.max(0);M[:3,3]=np.r_[128-(lo[:2]+hi[:2])/2,3.-lo[2]]
        m=mesh.copy();m.apply_transform(M);q=transform(xyz,M);low=lower_points(m,q);valid=np.isfinite(low).all(1)
        origins=low[valid]-[0,0,profile['roof_gap_mm']+.05];loc,idx,_=m.ray.intersects_location(origins,np.tile([0,0,-1.],(len(origins),1)),multiple_hits=False)
        heights=low[valid,2];horizontal=np.linalg.norm(np.diff(q[:,:2],axis=0),axis=1)
        rows.append(dict(orientation_id=f'O{i:03}',transform=M,rotation=R,build_fit=bool((hi-lo<[230,230,240]).all()),estimated_curtain_volume_mm3=float(np.sum(arc[valid]*heights)*2*w),maximum_curtain_height_mm=float(heights.max()) if len(heights) else 1e9,projection_collisions=int(len(idx)),invalid_rays=int((~valid).sum()),maximum_slenderness=float(heights.max()/(2*w)) if len(heights) else 1e9,open_direction_proxy=int(len(origins)-len(idx))))
    rows.sort(key=lambda r:(not r['build_fit'],r['invalid_rays'],r['estimated_curtain_volume_mm3'],r['maximum_curtain_height_mm'],r['projection_collisions'],r['maximum_slenderness'],-r['open_direction_proxy']))
    save(Path(out)/'orientations/ranking.json',rows);return rows
