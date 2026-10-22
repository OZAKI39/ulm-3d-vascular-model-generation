"""Actual extrusion evidence, independently measured contact and critical coverage."""
import json,re,zipfile
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .support_inputs import save,table
from .support_contact_planner import coverage
from .bambu_casting_mold import parse_paths,sample_segments
from .support_accessibility import component_access

def nozzle_gate(diameter,profile):
    if float(diameter)<max(profile['line_width'],profile['nozzle_diameter']):
        raise ValueError('SUPPORT_TIP_BELOW_RELIABLE_LINE_WIDTH')
    return True

def parse_actual_paths(path):
    # Both current engines emit FEATURE tags. Fail explicitly if a future
    # engine uses an unknown dialect; never report missing toolpaths as safe.
    with zipfile.ZipFile(path) as z:
        contents=[z.read(p) for p in z.namelist() if p.endswith('.gcode')]
        if not any(b'; FEATURE:' in p for p in contents):raise ValueError('GCODE_FEATURE_DIALECT_UNSUPPORTED')
    return parse_paths(path)

def candidate_gates(record,cfg):
    failed=[]
    if not record.get('source_unchanged',False):failed.append('UPSTREAM_GEOMETRY_MODIFIED')
    if record.get('coverage',{}).get('percent',0)<100:failed.append('CRITICAL_COVERAGE_INCOMPLETE')
    if record.get('coverage',{}).get('max_critical_unsupported_span',1e6)>cfg['critical_max_span_mm']:failed.append('CRITICAL_UNSUPPORTED_SPAN_EXCEEDED')
    if not record.get('slicer_success',False):failed.append('SLICER_PENDING_OR_FAILED')
    if record.get('critical_warnings'):failed.append('SLICER_CRITICAL_WARNING')
    if record.get('trapped_support_count')!=0:failed.append('TRAPPED_SUPPORT_RISK')
    if not record.get('no_illegal_collision',False):failed.append('SUPPORT_COLLISION_NOT_CLEARED')
    if record.get('missing_planned_tips'):failed.append('SLICER_DROPPED_SUPPORT_TIP')
    return failed

def contact_face_metrics(mesh,tree,profile):
    center=mesh.triangles_center;n=mesh.face_normals;w=profile['line_width'];gap=profile['support_top_gap'];h=profile['layer_height']
    probes=center-[0,0,gap+h*.5]
    dist,idx=tree.query(probes)
    mask=(n[:,2]<-.1)&(dist<np.sqrt((w*.7)**2+(h*1.5)**2))
    ids=np.flatnonzero(mask);parent=np.full(len(mask),-1);parent[ids]=np.arange(len(ids));adj=mesh.face_adjacency
    keep=mask[adj].all(axis=1);ij=parent[adj[keep]]
    graph=coo_matrix((np.ones(len(ij)),(ij[:,0],ij[:,1])),shape=(len(ids),len(ids))).tocsr()
    count,labels=connected_components(graph,directed=False) if len(ids) else (0,np.array([],int))
    areas=np.bincount(labels,weights=mesh.area_faces[ids],minlength=count)
    return dict(contact_area_mm2=float(mesh.area_faces[mask].sum()),tip_count=int((areas>=.015).sum()),all_face_patch_count=count,
        patch_area_threshold_mm2=.015,method='Actual support extrusion samples shifted by top gap; downward triangle area and mesh-connected patches',contact_face_ids=ids.tolist(),patch_areas_mm2=areas.tolist()),mask

def audit_slice(record,mesh,box,outer,raw,planned,profile,cfg,out):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    if not record.get('slicer_success'):
        result=dict(status='SLICE_FAILED',slicer_success=False,source_unchanged=True,contact_area_mm2=None,tip_count=None,coverage=None,
            trapped_support_count=None,no_illegal_collision=False,critical_warnings=record.get('warning_audit',{}).get('critical',[]))
        save(out/'support_metrics.json',result);return result
    paths=parse_actual_paths(record['archive']);segments=paths['segments'][paths['kind']==1];e=paths['extrusion_filament_mm'][paths['kind']==1]
    points=sample_segments(segments,.12);np.savez_compressed(out/'actual_support_paths.npz',segments=segments)
    if not len(points):raise ValueError('ACTUAL_SLICER_SUPPORT_PATHS_MISSING')
    tree=cKDTree(points);h=profile['layer_height'];w=profile['line_width'];gap=profile['support_top_gap'];touch=[];contact_rows=[]
    # XY and Z checks separately prevent an extrusion beside a vessel from
    # falsely becoming a contact underneath it.
    for p in raw:
        target=np.array(p['xyz'])-[0,0,gap+h*.5];indices=tree.query_ball_point(target,.45)
        pts=points[indices];xy=np.linalg.norm(pts[:,:2]-target[:2],axis=1) if len(pts) else np.array([])
        dz=np.abs(pts[:,2]-target[2]) if len(pts) else np.array([])
        actual=bool(np.any((xy<=w*.75+.04)&(dz<=h*1.5+.03)))
        if actual:touch.append(p)
        contact_rows.append(dict(point_id=p['point_id'],branch_id=p['branch_id'],arc_position=p['arc_position'],actual_support_underneath=actual))
    cov=coverage(raw,touch,cfg);metrics,face_mask=contact_face_metrics(mesh,tree,profile)
    # Contact count along arc is more stable on a fine vascular tessellation;
    # expose it alongside surface patch count rather than silently substituting.
    clusters=[]
    for branch in sorted({p['branch_id'] for p in touch}):
        ps=sorted([p for p in touch if p['branch_id']==branch],key=lambda p:p['arc_position']);current=[]
        for p in ps:
            if current and p['arc_position']-current[-1]['arc_position']>cfg['sample_step_mm']*1.6:clusters.append(current);current=[]
            current.append(p)
        if current:clusters.append(current)
    contact_tips=[g[len(g)//2] for g in clusters]
    # Approximate support component connectivity directly from occupied
    # extrusion cells; the voxel resolution and uncertainty are explicit.
    pitch=max(w,.25);cells=np.unique(np.floor(points/pitch).astype(int),axis=0);lookup={tuple(c):i for i,c in enumerate(cells)}
    edges=[]
    neighbors=[np.array([x,y,z]) for x in [-1,0,1] for y in [-1,0,1] for z in [-1,0,1] if (x,y,z)>(0,0,0)]
    for i,c in enumerate(cells):
        for delta in neighbors:
            j=lookup.get(tuple(c+delta))
            if j is not None:edges.append((i,j))
    ij=np.array(edges).reshape(-1,2);graph=coo_matrix((np.ones(len(ij)),(ij[:,0],ij[:,1])),shape=(len(cells),len(cells))).tocsr();count,labels=connected_components(graph,directed=False)
    sizes=np.bincount(labels);components=[]
    for label in range(count):
        ps=(cells[labels==label]+.5)*pitch;lo=ps.min(0)-pitch/2;hi=ps.max(0)+pitch/2
        shape=trimesh.creation.box(hi-lo);shape.apply_translation((hi+lo)/2);shape.metadata['support_occupied_points']=ps;components.append(shape)
    access=component_access(components,mesh,box,outer,cfg['probe_diameter_mm'])
    # Collision of actual support paths against the unmodified solid. Sampled
    # centerline containment is conservative evidence, not an exact support STL.
    inside=mesh.contains(points);box_inside=box.contains(points)
    # Extrusion paths lie on legal top-of-bottom boundaries; exclude a 0.02 mm
    # numerical skin before reporting a material penetration.
    bad=[]
    for name,obj,mask in [('vascular',mesh,inside),('box',box,box_inside)]:
        if mask.any():
            dist=trimesh.proximity.closest_point(obj,points[mask])[1];bad.append(dict(object=name,inside_sample_count=int(mask.sum()),deep_penetrations=int((dist>.03).sum()),maximum_depth_mm=float(dist.max())))
    # Bound the extrusion envelope as well as testing its center. G-code Z is
    # the top of the layer, so envelope centers are half a layer lower. A sphere
    # of radius max(width,height)/2 encloses its cross-section. The extra .08
    # covers half the sample spacing plus arc chord approximation error.
    import pyvista as pv
    from .sacrificial_fixture_review import polydata
    midpoints=points-[0,0,h/2]
    signed=pv.PolyData(midpoints).compute_implicit_distance(polydata(mesh))['implicit_distance']
    clearance=float(np.min(signed))-max(w,h)/2-.08
    envelope=dict(minimum_vascular_clearance_lower_bound_mm=clearance,passed=clearance>=-.02,
        sample_step_mm=.12,between_sample_and_arc_bound_mm=.08,extrusion_radius_bound_mm=max(w,h)/2)
    missing=[]
    for p in planned:
        nearby=[q for q in touch if q['branch_id']==p['branch_id'] and abs(q['arc_position']-p['arc_position'])<max(.5,p.get('tip_diameter',.6))]
        if not nearby:missing.append(p['contact_id'])
    material=float(e.sum()*np.pi*(1.75/2)**2);bif=sum(any(p['distance_to_bifurcation']<cfg['bifurcation_keepout_mm'] for p in group) for group in clusters)
    metrics['surface_contact_patch_count']=metrics['tip_count'];metrics['tip_count']=len(clusters)
    metrics['tip_count_method']='Approximate disconnected contact intervals on sampled vascular/port arcs; continuous wraparound may have multiple micro-gaps. Surface patch count also reported.'
    result=dict(status='ACTUAL_TOOLPATH_AUDITED',slicer_success=True,source_unchanged=True,coverage=cov,**metrics,
        arc_contact_patch_count=len(clusters),support_near_bifurcation_count=bif,estimated_support_volume_mm3=material,estimated_support_mass_g=material*.00104,
        support_path_length_mm=float(np.linalg.norm(segments[:,1]-segments[:,0],axis=1).sum()),trapped_support_count=access['trapped_support_count'],top_accessible_fraction=access['top_accessible_fraction'],
        support_component_sampling_mm=pitch,sampled_component_count=len(components),tiny_voxel_components_excluded=0,voxel_connectivity=26,
        no_illegal_collision=not any(r['deep_penetrations'] for r in bad) and envelope['passed'],collision=bad,extrusion_envelope=envelope,missing_planned_tips=missing,
        estimated_print_time_seconds=record.get('estimated_print_time_seconds'),critical_warnings=record.get('warning_audit',{}).get('critical',[]),slicer_warnings=record.get('warning_audit',{}).get('all',[]),
        evidence_limit='Extrusion sampling and contact/3 mm probe proxies; not physical strength or removal validation')
    result['hard_gate_failures']=candidate_gates(result,cfg)
    save(out/'support_metrics.json',result);save(out/'support_accessibility.json',access);table(out/'actual_contact_coverage.csv',contact_rows);np.save(out/'contact_face_mask.npy',face_mask)
    return result

def compare(records):
    baseline=records.get('BAMBU_NATIVE_BASELINE',{}).get('contact_area_mm2');rows=[]
    for name,r in records.items():
        area=r.get('contact_area_mm2');rows.append(dict(name=name,contact_area_mm2=area,support_tip_count=r.get('tip_count'),arc_contact_count=r.get('arc_contact_patch_count'),
            coverage_percent=(r.get('coverage') or {}).get('percent'),maximum_unsupported_mm=(r.get('coverage') or {}).get('max_critical_unsupported_span'),
            contact_area_reduction_percent=100*(1-area/baseline) if area is not None and baseline else None,bifurcation_contacts=r.get('support_near_bifurcation_count'),
            trapped_support_count=r.get('trapped_support_count'),top_accessible_fraction=r.get('top_accessible_fraction'),estimated_support_mass_g=r.get('estimated_support_mass_g'),
            estimated_print_time_seconds=r.get('estimated_print_time_seconds'),hard_gate_failures=r.get('hard_gate_failures',['SLICER_PENDING_OR_FAILED'])))
    feasible=[r for r in rows if not r['hard_gate_failures']]
    feasible.sort(key=lambda r:(r['contact_area_mm2'],r['bifurcation_contacts'],r['support_tip_count'],-r['top_accessible_fraction'],r['estimated_support_mass_g'],r['estimated_print_time_seconds'] or float('inf')))
    return rows,feasible[0]['name'] if feasible else None
