"""Actual X1C extrusion evidence, per-layer warning localization and roots."""
import json,zipfile,shutil,io
from pathlib import Path
import numpy as np
import trimesh
import pyvista as pv
from scipy.spatial import cKDTree
from shapely.geometry import LineString,GeometryCollection,Polygon
from shapely.ops import unary_union
from .support_inputs import save,table
from .support_slicer_adapter import discover_engine,engine_profiles
from .vascular_support_pipeline import sliced
from .support_qc import parse_actual_paths,contact_face_metrics
from .bambu_casting_mold import sample_segments
from .support_custom_qc import feature_preservation
from .sacrificial_fixture_review import polydata
from .support_v2_layers import analyse,coverage

def graph_features(graph,profile):
    rows=[];h=profile['layer_height']
    for p in graph['primitives']:
        lo,hi=np.array(p['bbox_mm']);center=(lo+hi)/2;width=p.get('width_mm',min((hi-lo)[:2]))
        zs=[.25,.5,.75] if hi[2]-lo[2]>4 else [.5]
        if p['kind']=='ROOF_PAD':zs=[1-h/2/max(hi[2]-lo[2],h)]
        for t in zs:
            xyz=center.copy();xyz[2]=lo[2]+t*(hi[2]-lo[2])
            if 'start' in p and 'end' in p:xyz=np.array(p['start'])+(np.array(p['end'])-p['start'])*t
            rows.append(dict(feature_id=p['primitive_id']+f'_{t:.3f}',primitive_id=p['primitive_id'],kind=p['kind'],roof_id=p.get('roof_id'),xyz=xyz.tolist(),width_mm=float(width)))
    return rows

def bridge_paths(archive):
    """Reuse the audited G-code parser, retaining only bridge FEATURE blocks.

    Filter an in-memory copy for analysis; never edit delivered machine code.
    Position and extrusion modes are still parsed for every command.
    """
    from .bambu_casting_mold import parse_paths
    payload=io.BytesIO()
    with zipfile.ZipFile(archive) as src,zipfile.ZipFile(payload,'w') as dst:
        dst.writestr('Metadata/project_settings.config',src.read('Metadata/project_settings.config'))
        for name in src.namelist():
            if name.endswith('.gcode'):
                lines=src.read(name).decode(errors='replace').splitlines()
                dst.writestr(name,'\n'.join('; FEATURE:' if line.startswith('; FEATURE:') and 'bridge' not in line.lower() else line for line in lines))
    payload.seek(0);return parse_paths(payload)['segments']

def bridge_preservation(segments,graph,all_segments=None,profile=None):
    rows=[]
    for p in graph['primitives']:
        if p['kind'] not in ['SHORT_BRIDGE','LONGITUDINAL_RAIL']:continue
        if p.get('span_mm',0)<=0:
            rows.append(dict(primitive_id=p['primitive_id'],span_mm=0.,preserved=True,status='NO_FREE_SPAN_SUPPORT_FOOTPRINTS_OVERLAP',actual_bridge_extrusions=0));continue
        lo,hi=np.array(p['bbox_mm']);mid=segments.mean(axis=1);delta=segments[:,1,:2]-segments[:,0,:2]
        axis=np.array(p.get('bridge_direction_xy',[1,0] if hi[0]-lo[0]>=hi[1]-lo[1] else [0,1]),float);axis/=np.linalg.norm(axis)
        length=np.linalg.norm(delta,axis=1);alignment=np.abs(delta@axis)/np.maximum(length,1e-9)
        mask=np.all(mid>=lo-[.3,.3,.05],axis=1)&np.all(mid<=hi+[.3,.3,.15],axis=1)&(alignment>.85)&(length>max(.2,p.get('span_mm',0)*.6))
        spanning=0;already_supported=0
        if all_segments is not None and profile is not None:
            h=profile['layer_height'];w=profile['line_width'];amid=all_segments.mean(axis=1);delta=all_segments[:,1,:2]-all_segments[:,0,:2];length=np.linalg.norm(delta,axis=1)
            aligned=np.abs(delta@axis)/np.maximum(length,1e-9)>.85
            eligible=np.all(amid>=lo-[w,w,.03],axis=1)&np.all(amid<=hi+[w,w,.03],axis=1)&(amid[:,2]<=lo[2]+1.7*h)&aligned&(length>max(.2,p.get('span_mm',0)*.6))
            for z in np.unique(amid[eligible,2]):
                previous=all_segments[(abs(amid[:,2]-(z-h))<.2*h)&np.all(amid[:,:2]>=lo[:2]-3*w,axis=1)&np.all(amid[:,:2]<=hi[:2]+3*w,axis=1)]
                previous_points=sample_segments(previous,w/3)
                if not len(previous_points):continue
                previous_tree=cKDTree(previous_points[:,:2])
                for segment in all_segments[eligible&(abs(amid[:,2]-z)<.001)]:
                    probes=sample_segments(segment[None],w/3);distance=previous_tree.query(probes[:,:2])[0]
                    if distance[0]<=w*.75 and distance[-1]<=w*.75:
                        if (distance>w*.75).any():spanning+=1
                        else:already_supported+=1
        rows.append(dict(primitive_id=p['primitive_id'],span_mm=p.get('span_mm'),bridge_direction_xy=axis.tolist(),actual_bridge_extrusions=int(mask.sum()),
            actual_two_end_anchored_span_extrusions=spanning,actual_fully_supported_member_extrusions=already_supported,
            preserved=bool(mask.any() or spanning or already_supported),evidence='Native Bridge FEATURE, or first-member-layer extrusion with both ends on actual previous-layer paths; fully supported members do not need a bridge tag'))
    return dict(passed=all(r['preserved'] for r in rows),rows=rows)

def audit(record,result,profile):
    folder=Path(record['folder']);target=folder/'bambu';graph=json.loads((folder/'support_graph_v2.json').read_text())
    if not result.get('slicer_success'):return dict(slicer_success=False,error=result.get('status'))
    with zipfile.ZipFile(result['archive']) as z:settings=json.loads(z.read('Metadata/project_settings.config'))
    paths=parse_actual_paths(result['archive']);segments=paths['segments'][paths['kind']==0];points=sample_segments(segments,.12)
    support=trimesh.load_mesh(folder/'custom_support.stl');vascular=trimesh.load_mesh(folder/'vascular_with_ports.stl');box=trimesh.load_mesh(folder/'box.stl')
    lo,hi=support.bounds;mask=np.all(points>=lo-.12,axis=1)&np.all(points<=hi+.12,axis=1);subset=points[mask]
    distance=pv.PolyData(subset).compute_implicit_distance(polydata(support))['implicit_distance'];custom=subset[distance<.08];np.savez_compressed(target/'actual_custom_paths.npz',points=custom)
    features=feature_preservation(custom,graph_features(graph,profile),profile);table(target/'feature_preservation.csv',features['rows'])
    roots=[];tree=cKDTree(custom)
    for r in graph['roots']:
        p=np.array(r['xyz']);p[2]=profile['first_layer_height'] if r['domain']=='BUILD_PLATE' else p[2]+profile['layer_height']*.6
        indices=tree.query_ball_point(p,np.hypot(r['foot_width_mm']/2,profile['layer_height']*2));near=custom[indices]
        present=bool(len(near) and np.any(np.abs(near[:,2]-p[2])<profile['layer_height']*1.5))
        roots.append(dict(root_id=r['root_id'],domain=r['domain'],target_xyz=p.tolist(),actual_start_path_present=present))
    table(target/'root_start_paths.csv',roots)
    roofs=json.loads((folder/'roofs.json').read_text());preserved={r['roof_id'] for r in features['rows'] if r['kind']=='ROOF_PAD' and r['preserved']}
    mid=segments.mean(axis=1)
    mid_distance=pv.PolyData(mid).compute_implicit_distance(polydata(support))['implicit_distance']
    support_segment_mask=mid_distance<.06
    for r in roofs:
        center=np.array(r['xyz'][:2]);radius=np.hypot(r['length_mm'],r['width_mm'])/2+profile['line_width'];mid=segments.mean(axis=1)
        eligible=support_segment_mask&(np.linalg.norm(mid[:,:2]-center,axis=1)<radius)&(mid[:,2]<=r['top_z_mm']+.025)&(mid[:,2]>=r['top_z_mm']-profile['layer_height']*1.6)
        lines=[LineString(s[:,:2]).buffer(profile['line_width']/2,cap_style=2) for s in segments[eligible]]
        from shapely import affinity
        from shapely.geometry import box as rect
        nominal=affinity.translate(affinity.rotate(rect(-r['length_mm']/2,-r['width_mm']/2,r['length_mm']/2,r['width_mm']/2),r['angle_rad'],use_radians=True),*center)
        r['_footprint']=unary_union(lines).intersection(nominal) if r['roof_id'] in preserved and lines else GeometryCollection()
    actual_bridges=bridge_paths(result['archive'])
    field=analyse(vascular,profile,'vascular',box);boxfield=analyse(box,profile,'box');vcov=coverage(field,roofs,profile,preserved,segments);bcov=coverage(boxfield,roofs,profile,preserved,segments)
    table(target/'actual_roof_coverage.csv',vcov['rows']);table(target/'actual_box_roof_coverage.csv',bcov['rows'])
    det=json.loads((folder/'detection.json').read_text());raw=det['raw'];rawtree=cKDTree([p['xyz'] for p in raw]);levels=np.unique(segments[:,:,2].mean(axis=1));locations=[]
    for row in vcov['rows']+bcov['rows']:
        xyz=[*row['representative_xy'],row['z_mm']];nearest=raw[rawtree.query(xyz)[1]];kind='box' if row['geometry_class']=='box' else 'port' if nearest['is_port'] else 'vascular'
        layer=int(np.argmin(abs(levels-row['z_mm'])))
        locations.append(dict(region_id=row['region_id'],gcode_layer_index=layer,gcode_layer_z_mm=float(levels[layer]),geometry_class=kind,xyz=xyz,
            attribution='INTENTIONAL_GAP_WITH_ACTUAL_ROOF' if row['covered'] else 'NO_SUFFICIENT_ACTUAL_ROOF',uncovered_area_mm2=row['uncovered_area_mm2'],
            location_source='Independent geometry-layer/G-code correlation; native warning gives no coordinates'))
    for row in record.get('self_buildability',{}).get('rows',[]):locations.append(dict(region_id=row['region_id'],gcode_layer_index=int(np.argmin(abs(levels-row['z_mm']))),geometry_class='support',attribution='SUPPORT_SELF_UNSUPPORTED',xyz=[*row['representative_xy'],row['z_mm']]))
    table(target/'floating_warning_localization.csv',locations)
    # Anatomical branch/port preservation is tested independently of support.
    vlo,vhi=vascular.bounds;eligible=np.all(points>=vlo-.12,axis=1)&np.all(points<=vhi+.12,axis=1);vpoints=points[eligible]
    vpoints=vpoints[pv.PolyData(vpoints).compute_implicit_distance(polydata(vascular))['implicit_distance']<.06]
    modeltree=cKDTree(vpoints);samples=[]
    for p in raw:
        xyz=np.array(p['xyz'])+[0,0,min(p['local_radius']*.5,.2)];distance=float(modeltree.query(xyz)[0]);limit=max(profile['line_width'],p['local_radius']*.75)
        samples.append(dict(point_id=p['point_id'],branch=p['branch_id'],geometry_class='port' if p['is_port'] else 'vascular',nearest_toolpath_mm=distance,passed=distance<=limit))
    table(target/'vascular_toolpath_preservation.csv',samples)
    metric,face_mask=contact_face_metrics(vascular,cKDTree(custom),profile);np.save(target/'contact_face_mask.npy',face_mask)
    native_warnings=result.get('warning_audit',{}).get('critical',[]);unresolved=sum(r['attribution']!='INTENTIONAL_GAP_WITH_ACTUAL_ROOF' for r in locations) if native_warnings else 0
    bridge=bridge_preservation(actual_bridges,graph,segments,profile);table(target/'bridge_toolpath_preservation.csv',bridge['rows'])
    audit_result=dict(slicer_success=True,auto_support_off=str(settings.get('enable_support')) in ['0','false','False'],native_support_toolpaths=int((paths['kind']==1).sum()),bridge_toolpaths=bridge,
        all_features_preserved=features['passed'],feature_preservation=features,all_roots_have_start_paths=all(r['actual_start_path_present'] for r in roots),
        vascular_toolpaths_preserved=all(r['passed'] for r in samples),vascular_failed_samples=sum(not r['passed'] for r in samples),
        actual_vascular_coverage={k:v for k,v in vcov.items() if k!='rows'},actual_box_coverage={k:v for k,v in bcov.items() if k!='rows'},
        contact_area_mm2=metric['contact_area_mm2'],warnings=native_warnings,unresolved_warning_count=unresolved,warning_locations=len(locations),
        estimated_print_time_seconds=result.get('estimated_print_time_seconds'),production_geometry=result['production_geometry'])
    save(target/'bambu_validation.json',audit_result);return audit_result

def slice_candidate(record,profile,discovery,secondary=False):
    folder=Path(record['folder']);engine=discover_engine('Bambu',discovery,folder);profiles=engine_profiles(discovery,'Bambu',{'enable_support':'0'})
    source=folder/'mold_with_custom_support.stl';result=sliced(engine,profiles,source,source,folder/'bambu',900);metrics=audit(record,result,profile)
    if result.get('slicer_success'):shutil.copy2(result['archive'],folder/f"candidate_{record['candidate_id']}_X1C.3mf")
    # Support alone retains legitimate box foundations; no vascular material
    # may be used to rescue a recursively unsupported support.
    support=trimesh.load_mesh(folder/'custom_support.stl');box=trimesh.load_mesh(folder/'box.stl');base=trimesh.boolean.union([support,box],engine='manifold')
    from .support_custom_qc import canonical_stl
    canonical_stl(base).export(folder/'support_with_approved_box_foundation.stl');path=folder/'support_with_approved_box_foundation.stl'
    solo=sliced(engine,profiles,path,path,folder/'bambu_support_only',900);metrics['support_only_slicer_success']=solo.get('slicer_success',False)
    if secondary:
        engine2=discover_engine('Orca',discovery,folder);p2=engine_profiles(discovery,'Orca',{'enable_support':'0'});other=sliced(engine2,p2,source,source,folder/'orca_secondary',900)
        metrics['orca_secondary']=dict(slicer_success=other.get('slicer_success',False),status=other.get('status'),warnings=other.get('warning_audit',{}),production_geometry=other.get('production_geometry'))
    save(folder/'bambu_validation.json',metrics);print('V2 SLICE',record['candidate_id'],metrics.get('slicer_success'),'feature preservation',metrics.get('all_features_preserved'),flush=True);return metrics
