"""Cut-first load-path and order-aware removal proxies for frozen s1-11 cores.

This module never predicts ABS fracture forces. Coverage is an independent
hard prerequisite; a failed coverage candidate has no production removal plan.
"""
from __future__ import annotations
import json,pickle,shutil
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from .support_inputs import save,table,digest,verify
from .core_support_geometry import boolean,union,overlap
from .core_support_removal import directions,wrap_audit
from .support_removal import swept_path


def lever_arm(force_point,contact_point,direction):
    """Perpendicular distance to the proposed applied force's line of action."""
    d=np.asarray(direction,dtype=float);d/=np.linalg.norm(d)
    return float(np.linalg.norm(np.cross(np.asarray(force_point)-contact_point,d)))


def direction_check(direction,normal,tangent,near_bifurcation=False):
    d=np.asarray(direction,dtype=float);d/=np.linalg.norm(d);n=np.asarray(normal,dtype=float);n/=np.linalg.norm(n);t=np.asarray(tangent,dtype=float);t/=np.linalg.norm(t)
    angle=float(np.degrees(np.arccos(np.clip(d@n,-1,1))))
    axial=abs(float(d@t))>.85
    return dict(passed=not axial and d@n>0 and (not near_bifurcation or angle<=60),removal_direction_vs_local_surface_normal_deg=angle,axial_drag=axial)


def wrap_status(angle,closed=False):
    if closed:return 'CLOSED_SUPPORT_CAGE'
    if angle>150:return 'WRAP_OVER_150_REJECTED'
    if angle>120:return 'HIGH_REMOVAL_RISK'
    return 'ACCEPTABLE' if angle>=90 else 'PREFERRED'


def cross_section_area(mesh,point,normal):
    """Measure the local lumen section, excluding other branches in the plane."""
    from shapely.geometry import Point
    path=mesh.section(plane_origin=np.asarray(point),plane_normal=np.asarray(normal))
    if path is None:return None
    flat,to3d=path.to_2D();q=np.linalg.inv(to3d)@np.r_[point,1.]
    areas=[p.area for p in flat.polygons_full if p.covers(Point(q[:2]))]
    # A missed centreline cannot silently borrow a distant vessel's section.
    return float(min(areas)) if areas else None


def actual_surface_normal(mesh,point):
    """Frozen STL facet normal; the downward ray direction is not a normal."""
    _,_,faces=mesh.nearest.on_surface([point]);face=int(faces[0])
    return np.asarray(mesh.face_normals[face]).copy(),face


def weak_ratio(width,thickness,vascular_area):
    if vascular_area is None or vascular_area<=0:return None
    return float(width*thickness/vascular_area)


def closed_cage(mesh,point,tangent):
    """A closed section hole enclosing the vascular axis is a topology trap."""
    from shapely.geometry import Point,Polygon
    path=mesh.section(plane_origin=point,plane_normal=tangent)
    if path is None:return False
    planar,to3d=path.to_2D();q=np.linalg.inv(to3d)@np.r_[point,1.]
    return any(Polygon(hole).covers(Point(q[:2])) for p in planar.polygons_full for hole in p.interiors)


def preserved_feature(feature,footprints,profile):
    """All witness sections of a designed link must retain actual extrusion."""
    from shapely.geometry import Point
    from .core_support_slicing import coordinate_tolerance
    rows=[]
    for q in np.asarray(feature['witness_points']):
        z=min(footprints,key=lambda z:abs(z-q[2]),default=None)
        valid=z is not None and abs(z-q[2])<=profile['layer_height']/2
        distance=float(footprints[z].distance(Point(q[:2]))) if valid else None
        rows.append(dict(xyz=q.tolist(),layer_top=z,distance_mm=distance,preserved=valid and distance<=coordinate_tolerance(q,profile)))
    return dict(passed=bool(rows) and all(r['preserved'] for r in rows),rows=rows)


def tool_access(point,vascular,obstacles,diameter=3.,allowed=None):
    """Finite cylindrical approach stopping at a material surface, no box."""
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular)
    for name,part in obstacles.items():manager.add_object(name,part)
    candidates=directions() if allowed is None else allowed
    attempts=[]
    for d in candidates:
        d=np.asarray(d,dtype=float);d/=np.linalg.norm(d)
        # Stop at the actual material-side lip, rather than aiming a solid
        # cylinder at the interior of a planned (not-yet-cut) kerf.
        nearest=None
        for part in obstacles.values():
            hits,_,_=part.ray.intersects_location([np.asarray(point)+d*100.],[-d],multiple_hits=True)
            if len(hits):
                q=hits[np.argmin(np.linalg.norm(hits-point,axis=1))]
                if np.linalg.norm(q-point)<=.8 and (nearest is None or np.linalg.norm(q-point)<np.linalg.norm(nearest-point)):nearest=q
        target=np.asarray(point) if nearest is None else nearest
        end=target+d*.05;start=end+d*100.
        M=trimesh.geometry.align_vectors([0,0,1],d);M[:3,3]=(start+end)/2
        probe=trimesh.creation.cylinder(radius=diameter/2,height=100.,sections=32,transform=M)
        hit=bool(manager.in_collision_single(probe));attempts.append(dict(direction=d.tolist(),collision=hit))
        if not hit:return dict(passed=True,direction=d.tolist(),probe_diameter_mm=diameter,attempts=attempts)
    return dict(passed=False,probe_diameter_mm=diameter,attempts=attempts)


def release_path(component,vascular,obstacles,normal=None,tangent=None,near_bifurcation=False):
    pivot=component.centroid;local=component.copy();local.apply_translation(-pivot);start=np.eye(4);start[:3,3]=pivot
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular)
    for name,part in obstacles.items():manager.add_object(name,part)
    extent=float(np.linalg.norm(vascular.extents)+np.linalg.norm(component.extents)+5);attempts=[]
    normals=[np.array([0.,0.,-1.]),*[np.array(d,dtype=float) for d in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0)]]]
    for d in normals:
        policy=direction_check(d,normal,tangent,near_bifurcation) if normal is not None and tangent is not None else dict(passed=True)
        if not policy['passed']:continue
        end=start.copy();end[:3,3]+=extent*d
        audit=swept_path(local,[],[start,end],static_manager=manager);attempts.append(dict(direction=d.tolist(),**audit))
        if audit['passed']:return dict(passed=True,direction=d.tolist(),poses=[start.tolist(),end.tolist()],attempts=attempts,**{k:v for k,v in policy.items() if k!='passed'})
    for axis in ['x','y']:
        for degrees in [-10,10,-20,20]:
            first=start.copy();first[2,3]-=.1
            tilted=first.copy();tilted[:3,:3]=Rotation.from_euler(axis,degrees,degrees=True).as_matrix()
            end=tilted.copy();end[2,3]-=extent
            policy=direction_check([0,0,-1],normal,tangent,near_bifurcation) if normal is not None and tangent is not None else dict(passed=True)
            if not policy['passed']:continue
            audit=swept_path(local,[],[start,first,tilted,end],static_manager=manager);attempts.append(dict(axis=axis,degrees=degrees,**audit))
            if audit['passed']:return dict(passed=True,direction=[0,0,-1],rotation_deg=degrees,axis=axis,poses=[x.tolist() for x in [start,first,tilted,end]],attempts=attempts,**{k:v for k,v in policy.items() if k!='passed'})
    return dict(passed=False,status='NO_RIGID_RELEASE_PATH_CERTIFIED',requires_vessel_deflection=False,vessel_deflection_permitted=False,attempts=attempts)


def sequential_plan(parts,vascular,contacts,cut_access,release,cut_safe=None):
    """Recheck access after each removal; failed pieces remain obstacles.

    Callbacks are geometric queries in production and permit isolated order
    regressions. No inaccessible cut or colliding motion is silently accepted.
    """
    remaining=dict(parts);steps=[];rounds=[];initial=[];round_id=0
    while remaining:
        round_id+=1;new=[];moved=[]
        for name in sorted(remaining):
            obstacles={k:v for k,v in remaining.items() if k!=name}
            if cut_safe is not None and not cut_safe.get(name,False):continue
            access=cut_access(name,remaining)
            if not access['passed']:continue
            if round_id==1:initial.append(name)
            new.append(name);motion=release(name,obstacles)
            if not motion['passed']:continue
            steps.append(dict(step=len(steps)+1,round=round_id,segment_id=name,stage_order=['A_CUT_BULK','B_REMOVE_CURTAIN','C_EXPOSE_ROOF','D_REMOVE_SHORT_CONTACT'],cut_points=contacts[name]['planned_cut_points'],access=access,removal=motion,obstacles_before=sorted(remaining)))
            remaining.pop(name);moved.append(name)
        rounds.append(dict(round=round_id,accessible=new,removed=moved,remaining=sorted(remaining)))
        if not moved:break
    return dict(passed=not remaining and bool(parts),complete=not remaining and bool(parts),initial_accessible=initial,rounds=rounds,steps=steps,remaining=sorted(remaining),all_notches_accessible_at_start_required=False,physical_validation=False)


def virtual_groups(g,out):
    """Finite virtual post-print cuts through existing seams, source untouched.

    These are inspection cuts, not changes to printing geometry. A cut is
    executable only after its real material-side approach has been certified.
    """
    out=Path(out);out.mkdir(parents=True,exist_ok=True);v=g['vascular'];s=g['support'];w=g['profile']['support_line_width'];h=g['profile']['layer_height'];cutters=[];rows=[];kerfs={}
    for n in g['notches']:
        tangent=np.asarray(n['tangent']);lateral=np.asarray(n['lateral']);xy=np.asarray(n['point'])[:2]
        # Cutting beyond the old notch apex must include the continuous ridge.
        near=np.linalg.norm(s.vertices[:,:2]-xy,axis=1)<3*w
        if not near.any():continue
        top=float(s.vertices[near,2].max()+h)
        M=np.eye(4);M[:3,0]=tangent;M[:3,1]=lateral;M[:3,2]=[0,0,1];M[:3,3]=[*xy,(top-h)/2]
        cutter=trimesh.creation.box([w,6*w,top+h],transform=M)
        intersection=overlap(cutter,v)
        manager=trimesh.collision.CollisionManager();manager.add_object('vascular',v)
        distance=float(manager.min_distance_single(cutter))
        passed=intersection<1e-7 and distance>0
        row=dict(cut_id=n['notch_id'],point=n['point'],normal=tangent,nearest_vascular_surface_mm=distance,vascular_intersection_mm3=intersection,passed=passed,stage='A',kerf_mm=w,initial_material_cut=True)
        rows.append(row)
        if passed:
            cutters.append(cutter);kerfs[n['notch_id']]=boolean([s,cutter],'intersection')
    after=boolean([s,union(cutters)],'difference') if cutters else s.copy()
    after.export(out/'virtual_after_cut.stl')
    parts={f'G{i:03}':part for i,part in enumerate(after.split(only_watertight=False)) if part.volume>1e-7}
    for name,part in parts.items():part.export(out/(name+'.stl'))
    return parts,rows,kerfs


def describe_components(g,parts,wrap_rows,cutrows):
    v=g['vascular'];bs=g['branches'];samples=np.vstack([b['lower'] for b in bs]);xyz=np.vstack([b['xyz'] for b in bs]);tangents=np.vstack([b['tangent'] for b in bs]);branch_names=np.concatenate([[b['branch_id']]*len(b['lower']) for b in bs]);arcs=np.concatenate([b['arc'] for b in bs]);tree=cKDTree(samples)
    forks=np.array([b['xyz'][0] for b in bs if not b['is_port']][1:]).reshape(-1,3);rows=[]
    for name,part in parts.items():
        distances,ids=tree.query(part.vertices);near=distances<=g['profile']['roof_gap_mm']+3*g['profile']['layer_height']
        ids=np.unique(ids[near]);ids=ids if len(ids) else np.array([tree.query(part.centroid)[1]])
        contacts=samples[ids]
        # A branched support's centroid can be in empty space. Proposed grip
        # points must lie on actual support material, even in a failed audit.
        force=part.nearest.on_surface([part.centroid])[0][0];direction=np.array([0.,0.,-1.])
        # The mean of contacts on several branches can lie in empty space.
        # Use a real lower-surface sample with the greatest downward-load
        # lever arm, rather than underestimating the group's worst arm.
        j=int(ids[np.argmax(np.linalg.norm(np.cross(contacts-force,direction),axis=1))])
        point=samples[j];normal,face_id=actual_surface_normal(v,point);tangent=tangents[j]
        spans={str(b):float(np.ptp(arcs[ids][branch_names[ids]==b])) for b in np.unique(branch_names[ids])}
        span=sum(spans.values());distance_bif=float(np.min(cKDTree(forks).query(contacts)[0])) if len(forks) else None
        ranges={str(b):(float(arcs[ids][branch_names[ids]==b].min()),float(arcs[ids][branch_names[ids]==b].max())) for b in np.unique(branch_names[ids])}
        branch_lookup={b['branch_id']:b for b in bs}
        local_angles=[r['angular_coverage_deg'] for r in wrap_rows if r['branch_id'] in ranges and ranges[r['branch_id']][0]-g['profile']['support_line_width']<=branch_lookup[r['branch_id']]['arc'][r['sample_index']]<=ranges[r['branch_id']][1]+g['profile']['support_line_width']]
        angle=max(local_angles,default=0.)
        # A large group may contact several branches. Its representative point
        # alone cannot rule out a closed trap around another contacted branch.
        cage_samples=[int(k) for k in ids if closed_cage(part,xyz[k],tangents[k])]
        closed=bool(cage_samples)
        lever=lever_arm(force,point,direction)
        necessary=[r for r in cutrows if np.linalg.norm(np.asarray(r['point'])[:2]-part.centroid[:2])<np.linalg.norm(part.extents[:2])/2+2]
        q=part.vertices[np.argmin(np.linalg.norm(part.vertices-point,axis=1))]
        area=cross_section_area(v,xyz[j],tangent)
        rows.append(dict(segment_id=name,contact_branch=';'.join(spans),contact_arc=spans,contact_length_mm=span,support_height_mm=float(part.extents[2]),segment_length_mm=float(np.linalg.norm(part.extents)),lever_arm_to_vessel_mm=lever,distance_to_bifurcation_mm=distance_bif,wrap_angle_deg=angle,wrap_status=wrap_status(angle,closed),closed_cage=closed,last_contact_piece_length_mm=span,last_contact_piece_volume_mm3=float(part.volume),last_contact_piece_status='UNSEPARATED_CONTACT_GROUP' if span>3 else 'SHORT_CONTACT_GROUP_UNCERTIFIED',planned_cut_points=[np.asarray(r['point']).tolist() for r in necessary],planned_cut_ids=[r['cut_id'] for r in necessary],all_cuts_geometrically_safe=all(r['passed'] for r in necessary),force_point=force.tolist(),contact_point=point.tolist(),normal=normal.tolist(),tangent=tangent.tolist(),tool_target=q.tolist(),removal_direction=direction.tolist(),removal_direction_vs_local_surface_normal_deg=direction_check(direction,normal,tangent)['removal_direction_vs_local_surface_normal_deg'],bbox_mm=part.bounds.tolist(),bbox_size_mm=part.extents.tolist(),nearest_vascular_cross_section_mm2=area,weak_link_width_mm=None,weak_link_thickness_mm=None,weak_link_area_mm2=None,geometric_weak_link_ratio=None,weak_link_status='NOT_OPTIMIZED_BEFORE_COVERAGE_CLOSURE',direct_long_lever_removal_prohibited=span>3 or lever>3))
        rows[-1].update(contact_reference_method='MAX_DOWNWARD_LOAD_LEVER_ON_ACTUAL_LOWER_SAMPLES',contact_sample_index=j,normal_source_face_index=face_id,closed_cage_checked_samples=len(ids),closed_cage_sample_indices=cage_samples)
    for row in rows:
        links=[x for x in g.get('weak_links',[]) if x['link_id'] in row['planned_cut_ids']]
        measured=[x for x in links if x.get('geometric_weak_link_ratio') is not None]
        row['weak_links_preserved']=bool(links) and all(x.get('actual_preservation',{}).get('passed',False) for x in links)
        if measured:
            link=min(measured,key=lambda x:x['geometric_weak_link_ratio'])
            row.update(weak_link_width_mm=link['weak_link_width_mm'],weak_link_thickness_mm=link['weak_link_thickness_mm'],weak_link_area_mm2=link['weak_link_area_mm2'],geometric_weak_link_ratio=link['geometric_weak_link_ratio'],weak_link_status='ACTUAL_FEATURE_PRESERVED' if row['weak_links_preserved'] else 'ACTUAL_FEATURE_DROPPED')
    return rows


def audit(g,out,coverage_closed=False):
    out=Path(out);parts,cuts,kerfs=virtual_groups(g,out/'virtual_components');wrap=wrap_audit(g['vascular'],g['support'],g['branches'],g['profile']);rows=describe_components(g,parts,wrap['rows'],cuts);byid={r['segment_id']:r for r in rows}
    # Inspections below are not an approved removal sequence. In particular a
    # single tool position cannot certify all cuts of a large connected group.
    # Kerf material remains an obstacle until that particular cut is reached.
    # This avoids the unsafe shortcut of checking access on a model in which
    # all future cuts have already been applied.
    remaining=dict(parts);uncut=dict(kerfs);events=[];rounds=[];initial=[];round_no=0
    while remaining:
        round_no+=1;new_access=[];removed=[];performed=[]
        for cut in cuts:
            key=cut['cut_id']
            if key not in uncut or not cut['passed']:continue
            obstacles={**remaining,**{'KERF_'+k:m for k,m in uncut.items()}}
            access=tool_access(cut['point'],g['vascular'],obstacles)
            if not access['passed']:continue
            new_access.append(key)
            if round_no==1:initial.append(key)
            events.append(dict(step=len(events)+1,round=round_no,action='CUT',cut_id=key,stage='A',access=access,nearest_vascular_surface_mm=cut['nearest_vascular_surface_mm'],obstacles_before=sorted(obstacles)))
            uncut.pop(key);performed.append(key)
        for name in sorted(remaining):
            r=byid[name]
            if any(k in uncut for k in r['planned_cut_ids']):continue
            if not r['all_cuts_geometrically_safe']:continue
            if r['direct_long_lever_removal_prohibited']:continue
            obstacles={**{k:v for k,v in remaining.items() if k!=name},**{'KERF_'+k:m for k,m in uncut.items()}}
            motion=release_path(remaining[name],g['vascular'],obstacles,r['normal'],r['tangent'],r['distance_to_bifurcation_mm'] is not None and r['distance_to_bifurcation_mm']<2)
            if not motion['passed']:continue
            events.append(dict(step=len(events)+1,round=round_no,action='REMOVE',stage='D',segment_id=name,removal=motion,obstacles_before=sorted(obstacles)))
            remaining.pop(name);removed.append(name)
        rounds.append(dict(round=round_no,new_accessible_cuts=new_access,cuts_performed=performed,removed=removed,remaining=sorted(remaining),uncut=sorted(uncut)))
        if not performed and not removed:break
    sequence=dict(passed=not remaining and not uncut,complete=not remaining and not uncut,initial_accessible=initial,rounds=rounds,steps=events,remaining=sorted(remaining),uncut=sorted(uncut),all_notches_accessible_at_start_required=False,kerf_material_retained_until_cut=True,physical_validation=False)
    for r in rows:
        step=next((s for s in sequence['steps'] if s.get('segment_id')==r['segment_id']),None)
        r['sequence_step']=None if step is None else step['step'];r['rigid_removal_pass']=step is not None
    report=dict(coverage_geometry_frozen=coverage_closed,removal_geometry_optimized=False,analysis_scope='READ_ONLY_VIRTUAL_DIAGNOSTIC' if not coverage_closed else 'FROZEN_COVERAGE_REMOVAL_INSPECTION',planned_segment_count=len(g['segments']),actual_component_count=len(parts),max_removal_component_arc_span_mm=max((r['contact_length_mm'] for r in rows),default=0),max_removal_component_bbox=[max((r['bbox_size_mm'][i] for r in rows),default=0) for i in range(3)],cut_count=len(cuts),safe_cut_count=sum(r['passed'] for r in cuts),max_wrap_angle=max((r['wrap_angle_deg'] for r in rows),default=0),wrap_over_120=sum(r['wrap_angle_deg']>120 for r in rows),wrap_over_150=sum(r['wrap_angle_deg']>150 for r in rows),closed_cages=sum(r['closed_cage'] for r in rows),max_lever_arm_mm=max((r['lever_arm_to_vessel_mm'] for r in rows),default=0),complete_sequential_plan=sequence['complete'],physical_validation=False,rows=rows,cut_checks=cuts,sequence=sequence)
    save(out/'safe_removal_diagnostic.json',report);table(out/'safe_removal_segments.csv',rows);table(out/'removal_load_paths.csv',rows)
    table(out/'weak_links.csv',[{k:r[k] for k in ['segment_id','weak_link_width_mm','weak_link_thickness_mm','weak_link_area_mm2','nearest_vascular_cross_section_mm2','geometric_weak_link_ratio','weak_link_status']} for r in rows])
    table(out/'final_contact_pieces.csv',[{k:r[k] for k in ['segment_id','last_contact_piece_length_mm','last_contact_piece_volume_mm3','direct_long_lever_removal_prohibited']} for r in rows])
    table(out/'bifurcation_removal_audit.csv',[r for r in rows if r['distance_to_bifurcation_mm'] is not None and r['distance_to_bifurcation_mm']<2],fields=list(rows[0]) if rows else ['segment_id'])
    save(out/'safe_removal_sequence.json',dict(production_plan=None,diagnostic_sequence=sequence,status='BLOCKED_BY_COVERAGE' if not coverage_closed else 'PENDING_ALL_SAFETY_GATES'))
    return report


def safety_gates(e,removal,record):
    from .core_support_closure import coverage_pass
    rotation=np.asarray(record['orientation']['rotation'])
    return dict(COVERAGE_100=coverage_pass(e),SUPPORT_SELF_ZERO=e['support_self']['self_unsupported_regions']==0,NO_ILLEGAL_INTERSECTION=record['intersection_mm3']<1e-7,NO_CLOSED_CAGE=removal['closed_cages']==0,NO_VESSEL_DEFLECTION_REQUIRED=removal['complete_sequential_plan'],COMPLETE_SEQUENTIAL_REMOVAL=removal['complete_sequential_plan'],NO_LONG_LEVER_CONTACT=all(not r['direct_long_lever_removal_prohibited'] for r in removal['rows']),WRAP_AT_MOST_150=removal['wrap_over_150']==0,HIGH_WRAP_EXCEPTION_JUSTIFIED=removal['wrap_over_120']==0 or bool(record.get('all_lower_wrap_trials_failed_at_100pct',False)),WEAK_LINKS_VALIDATED=bool(removal['rows']) and all(r['geometric_weak_link_ratio'] is not None and r.get('weak_links_preserved',False) for r in removal['rows']),REAL_SLICE=e['slice']['slicer_success'],AUTO_SUPPORT_OFF=e['auto_support_off'],FEATURES_PRESERVED=e['feature_preservation']['passed'],NATIVE_GEOMETRY_UNCHANGED=e.get('archive_geometry',{}).get('passed',False),SEPARATE_PARTS=e.get('independent_part_count')==2,RIGID_ONLY=bool(np.allclose(rotation.T@rotation,np.eye(3),atol=1e-9) and np.isclose(np.linalg.det(rotation),1.,atol=1e-9)),WARNINGS_RESOLVED=not e['unresolved_warnings'])


def prepare_removal(out,record,e,actual,engine):
    from .core_support_closure import closure_pass
    folder=Path(record['folder']);closed=closure_pass(e)
    if closed:
        save(out/'validation/coverage_freeze'/f"{record['candidate_id']}.json",dict(support=digest(folder/'custom_support.stl'),vascular=digest(folder/'vascular_core.stl'),evaluation=e))
        from .core_support_safe_geometry import optimize_after_closure
        record,e=optimize_after_closure(record,e,out/'candidates'/(record['candidate_id']+'_seam_trial'),actual,engine)
        folder=Path(record['folder'])
    g=pickle.load((folder/'geometry.pkl').open('rb'))
    diagnostic_dir=out/'removal_diagnostics'/record['candidate_id'];diagnostic_dir.mkdir(parents=True,exist_ok=True)
    signature=dict(geometry_pickle_sha256=digest(folder/'geometry.pkl'),analysis_code_sha256=digest(__file__),coverage_closed=closed)
    cache=diagnostic_dir/'input_signature.json'
    if cache.exists() and json.loads(cache.read_text())==signature and (diagnostic_dir/'safe_removal_diagnostic.json').exists():
        removal=json.loads((diagnostic_dir/'safe_removal_diagnostic.json').read_text())
    else:
        removal=audit(g,diagnostic_dir,closed);save(cache,signature)
    return dict(record=record,evaluation=e,geometry=g,removal=removal,diagnostic_dir=diagnostic_dir,gates=safety_gates(e,removal,record),coverage_closed=closed)


def choose_safe_state(states):
    """Safety comparisons are allowed only after all real production gates."""
    valid=[x for x in states if all(x['gates'].values())]
    if not valid:return None
    def rank(x):
        r=x['removal']
        return (r['max_wrap_angle'],r['max_lever_arm_mm'],r['max_removal_component_arc_span_mm'],r['cut_count'],x['record']['support_volume_mm3'])
    return min(valid,key=rank)


def finish(base,out,record,e,iterations,stop,hashes,actual,engine,alternatives=None):
    from .core_support_closure import closure_pass,diagnostic_rank
    from .core_support_travel import audit as travel_audit
    # Current failed cases need only a read-only diagnostic. If several routes
    # really close, inspect every closed route before comparing removal safety;
    # the first 0.25 mm entry may not displace a valid 0.20 mm production route.
    finalists=alternatives or [(record,e)]
    states=[prepare_removal(out,r,ev,actual,engine) for r,ev in finalists]
    state=choose_safe_state(states) or min(states,key=lambda x:diagnostic_rank(x['evaluation']))
    record=state['record'];e=state['evaluation'];g=state['geometry'];removal=state['removal'];diagnostic_dir=state['diagnostic_dir'];gates=state['gates'];closed=state['coverage_closed'];folder=Path(record['folder'])
    for file in diagnostic_dir.iterdir():
        if file.is_file() and file.name!='input_signature.json':shutil.copy2(file,out/file.name)
    if (out/'virtual_components').exists():shutil.rmtree(out/'virtual_components')
    shutil.copytree(diagnostic_dir/'virtual_components',out/'virtual_components')
    travel_file=folder/'travel_near_miss.json';travel=json.loads(travel_file.read_text()) if travel_file.exists() else travel_audit(e['slice']['archive'],folder,g['profile']['line_width'])
    failures=[k for k,v in gates.items() if not v];status='NEEDS_ADJUSTMENT' if failures else 'READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_REMOVAL_UNCALIBRATED'
    save(out/'safe_removal_sequence.json',dict(production_plan=None if failures else removal['sequence'],diagnostic_sequence=removal['sequence'],status='BLOCKED_BY_COVERAGE' if not closed else status))
    for src,dst in [('vascular_core.stl','vascular_core.stl'),('custom_support.stl','custom_support.stl'),('vascular_core_plus_support_review.stl','vascular_plus_support_review.stl'),('vascular_print_transform.json','print_transform.json'),('inverse_transform.json','inverse_transform.json')]:shutil.copy2(folder/src,out/'selected'/dst)
    shutil.copy2(e['slice']['archive'],out/'selected/X1C_project.3mf')
    for src,dst in [('coverage.json','coverage_final.json'),('support_self.json','support_self_final.json'),('layer_coverage.json','layer_coverage_final.json')]:shutil.copy2(folder/src,out/dst)
    protection=verify(hashes);save(out/'validation/source_protection.json',protection)
    # verify() has just rehashed every source and raises on any mismatch.
    # Persist the matching end hashes as well as the explicit zero-change gate.
    save(out/'validation/protected_after.json',hashes)
    summary=dict(status=status,production_selected=None if failures else record['candidate_id'],diagnostic_selected=record['candidate_id'],coverage=e['coverage'],layer_coverage=e['layer_coverage'],support_self=e['support_self'],gates=gates,failures=failures,removal={k:v for k,v in removal.items() if k not in ['rows','sequence','cut_checks']},source_protection=protection,actual_profile_hashes={k:digest(out/'bambu/profiles'/(k+'.json')) for k in actual},physical_validation=False,stop_conditions=stop,iteration_count=len(iterations),travel={k:v for k,v in travel.items() if k!='rows'},no_printer_upload_or_start=True)
    save(out/'safe_removal_summary.json',summary)
    from .core_support_safe_review import render,report
    render(base,out,g,record,e,removal,iterations);report(base,out,summary,record,e,removal)
    (out/'selected/README.md').write_text(f'# s1-11 本轮候选\n\n状态：{status}。\n\n'+('仅为诊断模型，尚未达到生产要求。' if failures else '通过软件检查，拆除力仍未经过实物校准。')+'\n\n血管及支撑为两个独立 part；Auto Support OFF。\n')
    print('CLOSURE_FINAL',status,record['candidate_id'],failures,flush=True)
    return summary
