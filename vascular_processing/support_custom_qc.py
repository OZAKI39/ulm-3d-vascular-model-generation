"""Geometry, root, self-printability and real auto-support-OFF extrusion audits."""
from __future__ import annotations
import io,json,zipfile
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from .support_inputs import save,table
from .support_detection import layers
from .support_contact_planner import coverage
from .support_breakaway import bottom_root_allowed
from .support_tree_planner import strut
from .support_custom_graph import minimum_distance
from .support_qc import parse_actual_paths,contact_face_metrics
from .bambu_casting_mold import sample_segments

def canonical_stl(mesh):
    """Validate the exact binary-STL precision before any production audit.

    Boolean seams can contain zero-area faces after float32 serialization.
    Remove only those degenerate faces, then require a closed positive volume.
    This does not smooth, shrink, remesh or change the frozen source geometry.
    """
    result=trimesh.load_mesh(io.BytesIO(mesh.export(file_type='stl')),file_type='stl')
    result.update_faces(result.nondegenerate_faces());result.remove_unreferenced_vertices()
    if not result.is_volume:raise ValueError('EXPORTED_SUPPORT_STL_NOT_CLOSED_VOLUME')
    return result

def intersection_volume(a,b):
    if not len(a.faces) or not len(b.faces):return 0.
    q=trimesh.boolean.intersection([a,b],engine='manifold')
    return abs(float(q.volume)) if len(q.faces) else 0.

def root_union_audit(support,vascular,box,feet,roots,inner):
    approved=trimesh.boolean.union(feet,engine='manifold') if feet else trimesh.Trimesh()
    box_intersection=trimesh.boolean.intersection([support,box],engine='manifold')
    outside=trimesh.boolean.difference([box_intersection,approved],engine='manifold') if len(box_intersection.faces) and len(approved.faces) else box_intersection
    vessel_overlap=intersection_volume(support,vascular);illegal_box=abs(float(outside.volume)) if len(outside.faces) else 0.
    rows=[]
    for i,part in enumerate(support.split()):
        bv=intersection_volume(part,box);vv=intersection_volume(part,vascular)
        rows.append(dict(component_id=i,bottom_overlap_volume_mm3=bv,vascular_and_ports_overlap_volume_mm3=vv))
    return dict(passed=vessel_overlap<1e-6 and illegal_box<1e-6 and all(bottom_root_allowed(r,inner) for r in roots),
        vascular_and_ports_overlap_volume_mm3=vessel_overlap,illegal_wall_or_bottom_overlap_volume_mm3=illegal_box,
        approved_bottom_overlap_volume_mm3=float(box_intersection.volume),roots_only_inner_bottom=all(bottom_root_allowed(r,inner) for r in roots),components=rows)

def self_printability(support,bottom,profile,out=None):
    result=layers(support,profile,.005,bottom);bad=[]
    # Complete support deficit polygons, not their centroids, are used. Every
    # component larger than a tiny tessellation sliver is a conservative failure.
    for z,shape in zip(result['levels'],result['unsupported_polygons']):
        pieces=list(shape.geoms) if hasattr(shape,'geoms') else [shape]
        for p in pieces:
            if p.area>.005 and max(p.bounds[2]-p.bounds[0],p.bounds[3]-p.bounds[1])>profile['line_width']:
                bad.append(dict(z_mm=float(z),area_mm2=float(p.area),bounds=list(p.bounds)))
    record=dict(passed=not bad,status='SELF_PRINTABILITY_SCREEN_PASSED' if not bad else 'SUPPORT_REQUIRES_SUPPORT',unsupported_regions=len(bad),
        unsupported_area_sum_mm2=sum(r['area_mm2'] for r in bad),rows=bad,legal_foundation_only='Frozen inner bottom; no vascular or side-wall material is used as support')
    if out:save(out/'self_printability.json',record);table(out/'self_printability.csv',bad,fields=['z_mm','area_mm2','bounds'])
    return record

def graph_features(graph,profile):
    rows=[]
    for t in graph['tips']:
        if t['root_connected']:rows.append(dict(feature_id=t['tip_id'],kind='tip',xyz=(np.array(t['xyz'])-[0,0,profile['layer_height']/2]).tolist(),width_mm=t['width_mm'],contact_ids=t['contact_ids']))
    for r in graph['roots']:rows.append(dict(feature_id=r['root_id'],kind='neck',xyz=r['estimated_break_location'],width_mm=r['neck_width_mm']))
    for e in graph['edges']:
        if e['kind'] in ['pillar','branch','bridge']:
            a=np.array(graph['nodes'][e['start']]['xyz']);b=np.array(graph['nodes'][e['end']]['xyz'])
            for fraction in [.25,.5,.75]:rows.append(dict(feature_id=f"{e['edge_id']}_{fraction}",kind=e['kind'],xyz=(a+(b-a)*fraction).tolist(),width_mm=e['diameter_mm']))
    return rows

def feature_preservation(points,features,profile):
    tree=cKDTree(points) if len(points) else None;rows=[]
    for f in features:
        p=np.array(f['xyz']);width=f['width_mm'];radius=width/2+profile['line_width']*.75
        indices=tree.query_ball_point(p,np.hypot(radius,profile['layer_height']*2)) if tree else []
        nearby=points[indices] if indices else np.empty((0,3))
        present=bool(len(nearby) and np.any((np.linalg.norm(nearby[:,:2]-p[:2],axis=1)<=radius)&(np.abs(nearby[:,2]-p[2])<=profile['layer_height']*1.5)))
        rows.append(dict(**f,preserved=present))
    return dict(passed=all(r['preserved'] for r in rows),missing=[r['feature_id'] for r in rows if not r['preserved']],rows=rows,
        status='FEATURES_PRESERVED' if all(r['preserved'] for r in rows) else 'SUPPORT_TIP_DROPPED')

def audit_slice(record,support,vascular,raw,contacts,graph,profile,cfg,out):
    if not record.get('slicer_success'):return dict(slicer_success=False,status='SLICER_VALIDATION_PENDING',coverage=None,features=None,auto_support_off=False)
    with zipfile.ZipFile(record['archive']) as z:settings=json.loads(z.read('Metadata/project_settings.config'))
    enabled=settings.get('enable_support');off=str(enabled) in ['0','false','False'] or enabled==['0']
    paths=parse_actual_paths(record['archive']);segments=paths['segments'][paths['kind']!=1];points=sample_segments(segments,.10)
    # Restrict ordinary model toolpaths to the explicit support volume; the
    # vessel and mold are ordinary model paths as well and must not count.
    lo,hi=support.bounds;keep=np.all(points>=lo-.12,axis=1)&np.all(points<=hi+.12,axis=1);subset=points[keep]
    import pyvista as pv
    from .sacrificial_fixture_review import polydata
    signed=pv.PolyData(subset).compute_implicit_distance(polydata(support))['implicit_distance'];custom_points=subset[signed<=.08]
    np.savez_compressed(out/'actual_custom_extrusion.npz',points=custom_points)
    features=feature_preservation(custom_points,graph_features(graph,profile),profile);preserved={r['feature_id'] for r in features['rows'] if r['kind']=='tip' and r['preserved']}
    ids={c for t in graph['tips'] if t['tip_id'] in preserved for c in t['contact_ids']};actual=[c for c in contacts if c['contact_id'] in ids];cov=coverage(raw,actual,cfg)
    metrics,mask=contact_face_metrics(vascular,cKDTree(custom_points),profile) if len(custom_points) else ({'contact_area_mm2':None},None)
    result=dict(slicer_success=True,auto_support_off=off,native_support_toolpath_count=int((paths['kind']==1).sum()),features=features,coverage=cov,contact_area_mm2=metrics['contact_area_mm2'],
        actual_preserved_tip_count=len(preserved),estimated_print_time_seconds=record.get('estimated_print_time_seconds'),critical_warnings=record.get('warning_audit',{}).get('critical',[]),
        actual_profile={k:settings.get(k) for k in ['layer_height','line_width','nozzle_diameter','enable_support','filament_type']},extrusion_point_count=len(custom_points))
    if mask is not None:np.save(out/'actual_contact_face_mask.npy',mask)
    save(out/'custom_slice_metrics.json',result);table(out/'feature_preservation.csv',features['rows']);return result

def hard_gates(record,cfg):
    failures=[]
    if not record.get('source_unchanged'):failures.append('UPSTREAM_GEOMETRY_MODIFIED')
    cov=record.get('actual_coverage') or record.get('planned_coverage') or {}
    if cov.get('percent',0)<100:failures.append('CRITICAL_INTERVAL_COVERAGE_INCOMPLETE')
    if cov.get('max_critical_unsupported_span',1e9)>cfg['critical_max_span_mm']:failures.append('CRITICAL_UNSUPPORTED_SPAN_EXCEEDED')
    if not record.get('self_printability',{}).get('passed'):failures.append('SUPPORT_REQUIRES_SUPPORT')
    if not record.get('union_audit',{}).get('passed'):failures.append('ILLEGAL_SUPPORT_INTERSECTION_OR_ROOT')
    if record.get('gap_audit',{}).get('passed') is not True:failures.append('SUPPORT_GAP_NOT_MET')
    if any(e['status']=='UNSTABLE_SUPPORT_EDGE' for e in record.get('edges',[])):failures.append('UNSTABLE_SUPPORT_EDGE')
    if record.get('removal',{}).get('trapped_support_count')!=0:failures.append('TRAPPED_SUPPORT_COMPONENT')
    actual=record.get('slice',{})
    if not actual.get('slicer_success'):failures.append('SLICER_VALIDATION_PENDING')
    else:
        if not actual.get('auto_support_off') or actual.get('native_support_toolpath_count',1)>0:failures.append('NATIVE_AUTO_SUPPORT_NOT_OFF')
        if (actual.get('features') or {}).get('passed') is not True:failures.append('SUPPORT_TIP_DROPPED')
        if actual.get('critical_warnings'):failures.append('SLICER_CRITICAL_WARNING')
        if record.get('contact_area_reduction_percent',-1)<25:failures.append('CONTACT_REDUCTION_TARGET_NOT_MET')
        support_only=record.get('support_only_slice',{})
        if not support_only.get('slicer_success') or not support_only.get('auto_support_off') or not (support_only.get('features') or {}).get('passed'):
            failures.append('SUPPORT_ONLY_SLICE_NOT_PRESERVED')
    return failures
