"""Material graph and ordered, collision-checked practical removal evidence.

All dimensions are geometry proxies. A hypothetical bulk/roof split is never
substituted for an executable cut or reported as a successful removal.
"""
from pathlib import Path
import json,pickle
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from .core_support_geometry import boolean,union,overlap
from .core_support_practical_geometry import material_components
from .core_support_practical_topology import contact_measure,cut_effect
from .core_support_safe_removal import (actual_surface_normal,closed_cage,
    preserved_feature,release_path,direction_check)
from .core_support_removal import directions,wrap_audit
from .support_inputs import save,table

REMOVAL_AUDIT_VERSION='practical-material-access-v7-construction-roof-residue'


def wrap_risk(angle,closed=False):
    if closed:return 'CLOSED_CAGE'
    if angle>180:return 'VERY_HIGH_REMOVAL_RISK'
    if angle>=150:return 'HIGH_RISK'
    if angle>=120:return 'MODERATE_HIGH'
    return 'PREFERRED'


def last_piece_check(length,ordinary_limit,near_bifurcation=False):
    # Scale with the actual local vessel diameter; bifurcation contact residue
    # has a stricter relative target, not an absolute hardcoded millimetre size.
    limit=ordinary_limit*(.5 if near_bifurcation else 1.)
    return dict(passed=length<=limit,length_mm=float(length),target_mm=float(limit),
                near_bifurcation=bool(near_bifurcation),physical_calibration=False)


def practical_gates(audit,removal,source_unchanged=True,archive_valid=True):
    s=audit['structural'];m=audit['macro']
    gates=dict(TRUE_MISSING_SUPPORT_ZERO=s['true_missing_support_count']==0,
        NO_DANGEROUS_SPAN=s['max_true_unsupported_arc_mm']==0,
        NO_FLOATING_VASCULAR=audit['floating_vascular_regions_unresolved']==0,
        NO_MACRO_SUPPORT_FAILURE=m['macro_fail_count']==0,
        NO_SUPPORT_VASCULAR_INTERSECTION=audit.get('support_vascular_intersection_mm3',0.)<1e-7,
        NO_CLOSED_CAGE=removal['closed_cage_count']==0,
        MAJOR_COMPONENTS_SEPARATED=removal['major_components_successfully_separated'],
        COMPLETE_SAFE_REMOVAL=removal['sequence']['complete'],
        SHORT_FINAL_CONTACT_PIECES=removal['short_final_contact_pieces'],
        ACTUAL_C3_PRESERVED=removal['thin_web_preserved'],
        VASCULAR_UNCHANGED=source_unchanged,
        REAL_BAMBU_SLICE=audit['native_slice_pass'] and archive_valid,
        VASCULAR_FEATURES_PRESERVED=audit['vascular_features_preserved'],
        AUTO_SUPPORT_OFF=audit['auto_support_off'],
        CRITICAL_WARNING_ZERO=audit['warnings']['critical_unresolved']==0)
    return gates


def actual_lip_access(kerf,vascular,obstacles,row,diameter=3.):
    """Probe actual material lips, never a centroid floating in a future gap."""
    if not len(kerf.faces):return dict(passed=False,reason='NO_ACTUAL_CUT_MATERIAL',attempts=[])
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular)
    for name,part in obstacles.items():manager.add_object(name,part)
    t=np.asarray(row['tangent']);side=np.array([-t[1],t[0],0.]);candidate_dirs=[side,-side,*directions()]
    # Low, real side surfaces are the intended bulk grip, with upper surfaces
    # considered only if they have a clear approach. No minus-Z tool approach.
    centers=kerf.triangles_center;normals=kerf.face_normals;attempts=[]
    for d in candidate_dirs:
        d=np.asarray(d,dtype=float)
        if np.linalg.norm(d)<1e-10:continue
        d/=np.linalg.norm(d)
        candidates=np.flatnonzero(normals@d>.8)
        if not len(candidates):continue
        candidates=sorted(candidates,key=lambda i:(centers[i,2],-kerf.area_faces[i]))
        # Select at most three material lips (bottom, middle, top), not a
        # stochastic search for a numerical collision loophole.
        ids=[candidates[0],candidates[len(candidates)//2],candidates[-1]]
        for i in dict.fromkeys(ids):
            point=centers[i];end=point+d*.05;M=trimesh.geometry.align_vectors([0,0,1],d);M[:3,3]=end+d*50
            probe=trimesh.creation.cylinder(radius=diameter/2,height=100,sections=32,transform=M)
            hit=bool(manager.in_collision_single(probe))
            attempts.append(dict(point=point.tolist(),direction=d.tolist(),collision=hit,material_face=int(i)))
            if not hit:return dict(passed=True,point=point.tolist(),direction=d.tolist(),probe_diameter_mm=diameter,attempts=attempts,claim='ACTUAL_MATERIAL_LIP_APPROACH_AND_SEPARATE_FINITE_BLADE_ENVELOPE')
    return dict(passed=False,probe_diameter_mm=diameter,attempts=attempts,reason='NO_CLEAR_ACTUAL_MATERIAL_LIP_APPROACH')


def component_metrics(parts,g,wrap_rows):
    rows=contact_measure(parts,g,exact=True);samples=np.vstack([b['lower'] for b in g['branches']]);centers=np.vstack([b['xyz'] for b in g['branches']]);tangents=np.vstack([b['tangent'] for b in g['branches']]);radii=np.concatenate([b['radius'] if 'radius' in b else b['radii'] for b in g['branches']])
    names=np.concatenate([[b['branch_id']]*len(b['arc']) for b in g['branches']]);indices=np.concatenate([np.arange(len(b['arc'])) for b in g['branches']]);forks=np.unique(np.array([b['xyz'][0] for b in g['branches'] if not b['is_port']][1:]),axis=0)
    wrapped={(r['branch_id'],r['sample_index']):r['angular_coverage_deg'] for r in wrap_rows}
    for part,row in zip(parts,rows):
        ids=np.asarray(row['sample_indices'],dtype=int);contacts=samples[ids]
        if not len(ids):
            row.update(effective_lever_arm_mm=0.,wrap_angle_deg=0.,closed_cage=False,bifurcation_contacts=[],grip_point=part.vertices[0].tolist(),ordinary_contact_limit_mm=None,short_contact=True);continue
        # A real triangle-centre grip, nearest to the contact centroid; its
        # tool accessibility is verified separately in sequence construction.
        candidates=part.triangles_center;target=contacts.mean(0);grip=candidates[np.argmin(np.linalg.norm(candidates-target,axis=1))];ds=np.linalg.norm(contacts-grip,axis=1);i=int(ids[np.argmin(ds)]);normal,facet=actual_surface_normal(g['vascular'],samples[i])
        bifdistance=float(cKDTree(contacts).query(forks)[0].min()) if len(forks) else None
        near=bifdistance is not None and bifdistance<=4*float(np.median(radii[ids]));limit=6*float(np.median(2*radii[ids]));short=last_piece_check(row['contact_arc_mm'],limit,near)
        angles=[wrapped.get((str(names[j]),int(indices[j]//2*2)),0.) for j in ids];angle=max(angles,default=0.)
        # Audit actual tangent-normal sections at each branch contact end,
        # middle, and every measured high-wrap point.
        inspect=set()
        for name in np.unique(names[ids]):
            subset=ids[names[ids]==name];inspect.update([int(subset[0]),int(subset[len(subset)//2]),int(subset[-1])]);inspect.update(int(j) for j in subset if wrapped.get((str(names[j]),int(indices[j]//2*2)),0.)>=150)
        cages=[j for j in sorted(inspect) if closed_cage(part,centers[j],tangents[j])]
        row.update(grip_point=grip.tolist(),grip_access_certified=False,nearest_vascular_contact=samples[i].tolist(),effective_lever_arm_mm=float(ds.min()),farthest_contact_distance_mm=float(ds.max()),distance_to_bifurcation_mm=bifdistance,bifurcation_contacts=[q.tolist() for q in forks if cKDTree(contacts).query(q)[0]<=4*np.median(radii[ids])],wrap_angle_deg=angle,wrap_risk=wrap_risk(angle,bool(cages)),closed_cage=bool(cages),closed_cage_sample_indices=cages,ordinary_contact_limit_mm=limit,short_contact=short['passed'],short_contact_check=short,normal=normal.tolist(),tangent=tangents[i].tolist(),normal_facet=facet,near_bifurcation=near,prohibited_directions=['AXIAL_DRAG','TOWARD_VASCULAR','WHOLE_COMPONENT_TWIST','LATERAL_BIFURCATION_PULL'])
    return rows


def bulk_cut_proposal(part,row,g,obstacles):
    """Check a finite far-from-contact bulk cut without inventing execution.

    A failed downstream roof sequence cannot be filled in by pretending a
    Boolean roof extraction was an accessible physical cut.
    """
    points=np.vstack([b['lower'] for b in g['branches']])[row['sample_indices']]
    if not len(points):return dict(component_id=row['component_id'],status='NO_VASCULAR_CONTACT')
    p=g['profile'];w=p['support_line_width'];z=float(points[:,2].min()-p['roof_gap_mm']-1.5)
    result=dict(component_id=row['component_id'],cut_id=row['component_id']+'_BULK',cut_z_mm=z,tool_diameter_mm=3.,executed=False)
    if z<=part.bounds[0,2]+w:return {**result,'status':'NO_BULK_HEIGHT_BELOW_CONTACT_WITH_TOOL_CLEARANCE'}
    cutter=trimesh.creation.box([part.extents[0]+2*w,part.extents[1]+2*w,w]);cutter.apply_translation([*part.bounds.mean(0)[:2],z])
    intersection=overlap(cutter,g['vascular']);result['finite_blade_vascular_intersection_mm3']=intersection
    if intersection>=1e-7:return {**result,'status':'FINITE_BULK_BLADE_INTERSECTS_ANOTHER_VASCULAR_BRANCH'}
    kerf=boolean([part,cutter],'intersection');split=material_components(boolean([part,cutter],'difference'))
    result['after_material_components']=len(split)
    if len(split)<2:return {**result,'status':'INEFFECTIVE_BULK_CUT'}
    access=actual_lip_access(kerf,g['vascular'],obstacles,dict(tangent=[1,0,0]));result['access']=access
    if not access['passed']:return {**result,'status':'BULK_CUT_TOOL_ACCESS_BLOCKED'}
    contacts=contact_measure(split,g,exact=True);bulk=[(m,r) for m,r in zip(split,contacts) if r['contact_arc_mm']==0];contact=[r for r in contacts if r['contact_arc_mm']>0]
    result.update(point=access['point'],approach_direction=access['direction'],detached_bulk_volume_mm3=sum(float(m.volume) for m,_ in bulk),remaining_roof_max_contact_arc_mm=max((r['contact_arc_mm'] for r in contact),default=0.))
    # The connector can be approached, but short roof pieces and their removal
    # must also be certified before this can become a complete safe sequence.
    result['status']='BULK_ENTRY_CERTIFIED_DOWNSTREAM_SHORT_ROOF_SEQUENCE_UNRESOLVED'
    return result


def audit_removal(folder):
    folder=Path(folder);out=folder/'practical_removal';out.mkdir(exist_ok=True)
    g=pickle.load((folder/'geometry.pkl').open('rb'));p=g['profile'];support=g['support'];vascular=g['vascular'];plan=json.loads((folder/'strategic_cut_plan.json').read_text());selected=pickle.load((folder/'strategic_cut_solids.pkl').open('rb'))
    actual=pickle.load((folder/'actual_footprints.pkl').open('rb'));precision=json.loads((folder/'gcode_precision.json').read_text());profile={**p,'gcode_xy_half_cell_mm':precision['half_coordinate_cell_mm']}
    links=[]
    from .core_support_closure import local_box
    for link in plan['weak_links']:
        evidence=preserved_feature(link,actual['support'],profile)
        cut=next(x['row'] for x in selected if x['row']['cut_id']==link['link_id']);mid=np.asarray(link['center']);thick=link['thickness_mm'];w=p['support_line_width']
        band=local_box(mid[:2],mid[2]+thick/2,cut['lateral_width_mm'],w/2,np.asarray(cut['tangent'])[:2],mid[2]-thick/2)
        measured=overlap(support,band)/(w/2)
        intended=link['width_mm']*link['thickness_mm'];isolated=measured<=intended+1e-4
        links.append(dict(link_id=link['link_id'],fixture=link['source_fixture'],**evidence,slicer_preserved=evidence['passed'],measured_neck_section_mm2=measured,intended_C3_section_mm2=intended,isolated_C3_geometry=isolated,production_C3_valid=evidence['passed'] and isolated))
    before=material_components(support);remaining=support.copy();cuts=[];kerfs={}
    for item in selected:
        cutter=item['cutter'];kerf=boolean([remaining,cutter],'intersection');after=boolean([remaining,cutter],'difference');effect=cut_effect(material_components(remaining),material_components(after),g);row={**item['row'],**effect};row['vascular_intersection_mm3']=overlap(cutter,vascular)
        if not effect['effective']:row['status']='INEFFECTIVE_CUT_REMOVED';cuts.append(row);continue
        row['status']='VIRTUAL_MATERIAL_GRAPH_CUT';cuts.append(row);kerfs[row['cut_id']]=kerf;remaining=after
    parts=material_components(remaining);wrap=wrap_audit(vascular,support,g['branches'],p);wrap['method']='Near-surface angular geometry proxy; 120 degrees is not an acceptance gate';rows=component_metrics(parts,g,wrap['rows']);part_map={r['component_id']:m for r,m in zip(rows,parts)}
    source_folder=folder;visited=set()
    while source_folder not in visited:
        visited.add(source_folder);record=json.loads((source_folder/'geometry.json').read_text())
        if not record.get('practical_parent'):break
        source_folder=Path(record['practical_parent'])
    else:raise ValueError('CYCLIC_SUPPORT_GEOMETRY_PROVENANCE')
    original=pickle.load((source_folder/'geometry.pkl').open('rb'))
    additions=boolean([support,original['support']],'difference')
    roof_template=union([original['roofs'],additions]) if len(additions.faces) else original['roofs']
    contact_rows=[];contact_meshes={};bulk_meshes={}
    for row,part in zip(rows,parts):
        # Extract real top-layer residue for diagnosis. This is a proposed
        # bulk/roof separation, not a claim that an unverified cut took place.
        # Use the original construction roof plus actual local additions.
        # A support-minus-downshift shell also contains vertical curtain
        # ledges and must not be mistaken for final vessel-contact residue.
        residue=boolean([part,roof_template],'intersection');bulk=boolean([part,residue],'difference')
        bulk_meshes[row['component_id']]=bulk
        for i,res in enumerate(material_components(residue)):
            if res.volume<=1e-7:continue
            measure=contact_measure([res],g,exact=True)[0]
            if not measure['sample_indices']:continue
            length=float(max(res.bounding_box_oriented.primitive.extents))
            key=f"{row['component_id']}_R{i:02}";check=last_piece_check(max(measure['contact_arc_mm'],length),row['ordinary_contact_limit_mm'],row['near_bifurcation'])
            contact_meshes[key]=res
            contact_rows.append(dict(piece_id=key,parent_component=row['component_id'],last_contact_piece_length_mm=length,length_method='ORIENTED_BOUNDING_BOX_LONGEST_DIMENSION',last_contact_piece_volume_mm3=float(res.volume),contact_arc_mm=measure['contact_arc_mm'],distance_to_bifurcation_mm=row['distance_to_bifurcation_mm'],wrap_angle_deg=row['wrap_angle_deg'],branch_contacts=measure['branch_contacts'],near_bifurcation=row['near_bifurcation'],short_piece_check=check,bulk_separation_cut_verified=False,state='DIAGNOSTIC_ROOF_RESIDUE_NOT_EXECUTABLE_CUT_PROOF'))
    # Reconstruct all uncut kerf material as obstacles. No future kerf is
    # treated as empty air before its approach has actually passed.
    obstacles=dict(part_map);pending=dict(kerfs);steps=[];attempts=[];cut_by_id={r['cut_id']:r for r in cuts};row_by_id={r['component_id']:r for r in rows};removed=[]
    manager=trimesh.collision.CollisionManager()
    dependencies={name:[] for name in part_map}
    for cutid,kerf in kerfs.items():
        manager=trimesh.collision.CollisionManager();manager.add_object(cutid,kerf)
        for name,part in part_map.items():
            if manager.min_distance_single(part)<1e-5:dependencies[name].append(cutid)
    while pending or obstacles:
        progress=False
        for cutid,kerf in list(pending.items()):
            current={**obstacles,**{'UNCUT_'+k:v for k,v in pending.items()}};check=actual_lip_access(kerf,vascular,current,cut_by_id[cutid]);attempts.append(dict(cut_id=cutid,obstacles=sorted(current),**check))
            if not check['passed']:continue
            steps.append(dict(step=len(steps)+1,action='CUT_STRATEGIC_CONNECTOR',cut_id=cutid,point=check['point'],access_direction=check['direction'],nearby_branch=cut_by_id[cutid]['branch'],prohibited='Do not pull vessel, twist tree, or enlarge cut toward vascular',obstacles_before=sorted(current),evidence=check))
            pending.pop(cutid);progress=True
        for name,part in list(obstacles.items()):
            row=row_by_id[name]
            if any(k in pending for k in dependencies[name]):continue
            if not row.get('short_contact',False) or row['closed_cage']:continue
            current={**{k:v for k,v in obstacles.items() if k!=name},**{'UNCUT_'+k:v for k,v in pending.items()}}
            grip=actual_lip_access(part,vascular,{**current,name:part},dict(tangent=row.get('tangent',[1,0,0])))
            row['grip_access']=grip
            if not grip['passed']:continue
            row['grip_access_certified']=True;row['grip_point']=grip['point']
            samples=np.vstack([b['lower'] for b in g['branches']])[row['sample_indices']]
            ds=np.linalg.norm(samples-np.asarray(grip['point']),axis=1)
            row['effective_lever_arm_mm']=float(ds.min());row['farthest_contact_distance_mm']=float(ds.max())
            motion=release_path(part,vascular,current,row.get('normal'),row.get('tangent'),row.get('near_bifurcation',False));row['rigid_release']=motion
            if not motion['passed']:continue
            steps.append(dict(step=len(steps)+1,action='REMOVE_SHORT_COMPONENT',component_id=name,grip_point=row['grip_point'],direction=motion['direction'],nearby_branches=list(row['branch_contacts']),prohibited_directions=row['prohibited_directions'],obstacles_before=sorted(obstacles),motion=motion));obstacles.pop(name);removed.append(name);progress=True
        if not progress:break
    complete=not pending and not obstacles and bool(parts)
    # For long remaining components the geometric roof diagnostic is not an
    # executable bulk-cut plan. Do not manufacture safe instructions from it.
    blocked=[dict(component_id=name,reason='LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS' if not row_by_id[name].get('short_contact') else 'NO_CERTIFIED_RIGID_RELEASE_OR_PENDING_CUT',pending_cuts=[k for k in dependencies[name] if k in pending]) for name in obstacles]
    current_obstacles={**obstacles,**{'UNCUT_'+k:v for k,v in pending.items()}}
    bulk_attempts=[]
    for name,part in obstacles.items():
        row=row_by_id[name]
        if row.get('short_contact'):continue
        proposal=bulk_cut_proposal(part,row,g,current_obstacles);bulk_attempts.append(proposal)
        row['bulk_cut_point']=proposal.get('point');row['bulk_cut_id']=proposal['cut_id'];row['bulk_cut_status']=proposal['status']
    sequence=dict(complete=complete,steps=steps,remaining_components=list(obstacles),pending_cuts=list(pending),blocked=blocked,all_cuts_initially_accessible_required=False,requires_vessel_deflection=False,vessel_deflection_permitted=False,physical_removal_validated=False,diagnostic_only=not complete)
    result=dict(initial_components=len(before),major_component_count=len(parts),major_components_successfully_separated=len(parts)>len(before) and all(r['effective'] for r in cuts if r['status']=='VIRTUAL_MATERIAL_GRAPH_CUT'),maximum_component_contact_arc_mm=max((r['contact_arc_mm'] for r in rows),default=0.),max_effective_lever_arm_mm=max((r['effective_lever_arm_mm'] for r in rows),default=0.),max_farthest_contact_distance_mm=max((r.get('farthest_contact_distance_mm',0) for r in rows),default=0.),max_last_contact_piece_length_mm=max((r['last_contact_piece_length_mm'] for r in contact_rows),default=0.),closed_cage_count=sum(r['closed_cage'] for r in rows),maximum_wrap_angle_deg=wrap['maximum_wrap_angle_deg'],very_high_wrap_count=sum(r['wrap_angle_deg']>180 for r in rows),thin_web_preserved=bool(links) and all(r['passed'] for r in links),thin_web_checks=links,short_final_contact_pieces=bool(contact_rows) and all(r['short_piece_check']['passed'] for r in contact_rows) and complete,major_components=rows,final_contact_pieces=contact_rows,cuts=cuts,sequence=sequence,primary_cut_count=len(kerfs),actual_executed_virtual_cut_count=sum(r['action']=='CUT_STRATEGIC_CONNECTOR' for r in steps),physical_break_force_validated=False)
    for name,meshes in [('major_components',part_map),('contact_residues',contact_meshes),('curtain_bulk',bulk_meshes),('cut_kerfs',kerfs)]:
        (out/name).mkdir(exist_ok=True)
        for old in (out/name).glob('*.stl'):
            if old.stem not in meshes:
                archive=out/'obsolete_diagnostic_meshes'/name;archive.mkdir(parents=True,exist_ok=True)
                old.replace(archive/old.name)
        for key,m in meshes.items():m.export(out/name/(key+'.stl'),file_type='stl_ascii')
    result['thin_web_preserved']=bool(links) and all(r['production_C3_valid'] for r in links)
    result['actual_extrusion_thin_link_witnesses_preserved']=bool(links) and all(r['slicer_preserved'] for r in links)
    result['audit_version']=REMOVAL_AUDIT_VERSION
    result['contact_residue_definition']='Intersection with original construction roofs and actual local added material; hypothetical bulk separation, not executable-cut proof'
    result['bulk_cut_attempts']=bulk_attempts
    initial_metrics=component_metrics(before,g,wrap['rows'])
    result['before_split_metrics']=initial_metrics
    result['before_max_component_contact_arc_mm']=max((r['contact_arc_mm'] for r in initial_metrics),default=0.)
    result['before_max_effective_lever_arm_mm']=max((r['effective_lever_arm_mm'] for r in initial_metrics),default=0.)
    save(out/'removal_audit.json',result);save(out/'safe_removal_sequence.json',sequence);save(out/'sequential_access.json',dict(attempts=attempts,complete=complete,pending_cuts=list(pending),obstacle_updates_verified=True));save(out/'wrap_audit.json',wrap)
    table(out/'strategic_cut_seams.csv',cuts);table(out/'major_removal_components.csv',rows);table(out/'final_contact_pieces.csv',contact_rows)
    save(out/'bulk_cut_attempts.json',bulk_attempts)
    print('PRACTICAL_REMOVAL',folder.name,'major',len(parts),'maxcontact',result['maximum_component_contact_arc_mm'],'C3',result['thin_web_preserved'],'complete',complete,flush=True)
    return result
