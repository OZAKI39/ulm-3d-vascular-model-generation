"""s1-10 assembly, nominal solvent topology and post-dissolution trim proofs."""
from __future__ import annotations
import csv
import numpy as np
import networkx as nx
import trimesh
from scipy.spatial.transform import Rotation
from shapely.geometry import Polygon
from .support_inputs import save,table,swc_branches
from .casting_fixture_geometry import shape_mesh,bounds_box
from .casting_fixture_v2_geometry import world,port_transform,fuse
from .casting_fixture_v2_audit import manifold,cad_manifold,mesh,contact_zone
from .casting_fixture_audit import sweep_translation
from .support_removal import swept_path
from .casting_mold_qc import fcl_distance


def assembly(l,g,cfg,out):
    core=l['core'];central=l['central'];mold=shape_mesh(g['mold']);H=l['high'][2];tol=cfg['intersection_volume_tolerance_mm3']
    # Zero-gap unilateral contact cannot pass a positive-distance FCL gate.
    # Prove the open path with FCL, then the last seating interval by the
    # vertical hypograph construction of lands on actual downward facets.
    print('s1-10: continuous insertion and nominal seating',flush=True)
    unlanded=shape_mesh(g['mold_without_lands']);full=sweep_translation(core,[unlanded],[0,0,H+5],cfg['sweep_tolerance_mm'])
    elevated=core.copy();elevated.apply_translation([0,0,.02]);near=sweep_translation(elevated,[mold],[0,0,H+5],cfg['sweep_tolerance_mm'])
    penetration=float((manifold(core)^cad_manifold(g['mold'])).volume())
    contactrows=[]
    for kind in ['fixtures','retainers','guides']:
        for p in l['ports']:
            k=p['port_id'];part=shape_mesh(g[kind][k]);d,cp,fp=fcl_distance(core,part);dc,_,_=fcl_distance(central,part)
            approved=contact_zone(cp,p,distal_stop=p['role']=='MASTER') if kind!='guides' else d>.05
            contactrows.append(dict(part=kind+'_'+k,port_id=k,minimum_core_gap_mm=d,central_gap_mm=dc,central_contact=dc<=1e-6,closest_core_point=cp.tolist(),approved_extension_zone=approved))
    witnesses_ok=all(contact_zone(w,p,tol=2e-4) for p in l['ports'] for r in g['access_records'] if r['port_id']==p['port_id'] for w in r['contact_witnesses'])
    monotonic=all(np.all(l['core'].face_normals[r['source_face_ids'],2]<-.5) for r in g['access_records'])
    top=dict(passed=bool(full['passed'] and near['passed'] and penetration<tol and witnesses_ok and monotonic),positive_clearance_without_lands=full,full_mold_to_0p02mm_standoff=near,final_seating_interval_mm=.02,final_seating_proof='Land is the vertical hypograph of downward-facing frozen facets over fixed XY polygons; raising core increases separation at every XY. Reverse path seats without penetration.',downward_face_condition=monotonic,final_interference_mm3=penetration)
    present=[mold,core];insertions=[];locks=[];manager=trimesh.collision.CollisionManager();manager.add_object('MOLD',mold)
    for p in l['ports']:
        k=p['port_id'];part=shape_mesh(g['retainers'][k]);slide=g['local_parts'][k]['retainer_slide_mm'];a=np.eye(4);a[:3,3]=-np.array(p['axis'])*slide+[0,0,H+5];b=np.eye(4);b[:3,3]=-np.array(p['axis'])*slide
        r=swept_path(part,present,[a,b,np.eye(4)],tolerance=cfg['sweep_tolerance_mm']);insertions.append(dict(port_id=k,sequence=['TOP_DOWN','WALLWARD_SLIDE'],slide_mm=slide,**r));present.append(part)
        m=np.eye(4);m[2,3]=.8;locks.append(dict(port_id=k,blocked=bool(manager.in_collision_single(part,transform=m))))
    cm=trimesh.collision.CollisionManager()
    for i,part in enumerate([mold,*[shape_mesh(r) for r in g['retainers'].values()]]):cm.add_object(str(i),part)
    pivot=np.mean([p['trim_origin'] for p in l['ports']],axis=0);poses=[]
    for name,ax,ang,dz in [('Z+',None,0,1.),('pitch+','y',2,0),('pitch-','y',-2,0),('roll+','x',2,0),('roll-','x',-2,0)]:
        m=np.eye(4)
        if ax:m[:3,:3]=Rotation.from_euler(ax,ang,degrees=True).as_matrix()
        m[:3,3]=pivot-m[:3,:3]@pivot+[0,0,dz];poses.append(dict(name=name,degrees=ang,dz_mm=dz,blocked=bool(cm.in_collision_single(core,transform=m))))
    walls_missing=float(g['walls'].cut(g['mold']).Volume());bottom_missing=float(g['bottom'].cut(g['mold']).Volume());centralgap=fcl_distance(central,mold)[0]
    result=dict(top_insertion_pass=top['passed'],top_insertion=top,retainer_insertions=insertions,retainer_insertion_pass=all(x['passed'] for x in insertions),contacts=contactrows,fixture_contacts=len(g['access_records']),central_tree_contacts=sum(x['central_contact'] for x in contactrows),minimum_central_auxiliary_clearance_mm=min(centralgap,*[x['central_gap_mm'] for x in contactrows]),all_contacts_on_approved_extensions=witnesses_ok and all(x['approved_extension_zone'] for x in contactrows),walls_closed=walls_missing<tol,bottom_closed=bottom_missing<tol,wall_missing_volume_mm3=walls_missing,bottom_missing_volume_mm3=bottom_missing,top_open=top['passed'],core_interference_mm3=penetration,anti_float_geometric_pass=all(x['blocked'] for x in poses+locks),anti_float_poses=poses,retainer_locks=locks,physical_force_validated=False,status='ANTI_FLOAT_FORCE_NOT_PHYSICALLY_VALIDATED')
    save(out/'assembly_audit.json',result);table(out/'fixture_contacts.csv',contactrows);return result


def graph_status(graph,ports):
    reachable=[p for p in ports if 'external_bath' in graph and p in graph and nx.has_path(graph,'external_bath',p)]
    return dict(status='SOLVENT_ACCESS_GEOMETRIC_PASS' if len(reachable)==len(ports) and ports else 'SOLVENT_ACCESS_BLOCKED',reachable_ports=reachable)


def contact_surface_evidence(core,land,witnesses,area,cfg):
    points=np.asarray(witnesses)
    _,core_dist,_=core.nearest.on_surface(points)
    _,land_dist,_=land.nearest.on_surface(points)
    gap=max(float(np.max(core_dist)),float(np.max(land_dist)))
    return dict(passed=area>cfg['contact_area_tolerance_mm2'] and gap<cfg['contact_distance_tolerance_mm'],maximum_contact_surface_error_mm=gap,core_error_mm=float(np.max(core_dist)),land_error_mm=float(np.max(land_dist)))


def network_access(swc,ports,accessible,mapping=None):
    data=np.loadtxt(swc,comments='#',ndmin=2);G=nx.Graph();G.add_nodes_from(int(r[0]) for r in data)
    for r in data:
        if int(r[6]) in G:G.add_edge(int(r[0]),int(r[6]))
        elif int(r[6])!=-1:raise ValueError('SWC_PARENT_MISSING')
    starts={p['port_id']:int(p['swc_id']) for p in ports if p['port_id'] in accessible}
    if any(v not in G for v in starts.values()):raise ValueError('PORT_SWC_ANCHOR_MISSING')
    distances={k:nx.single_source_shortest_path_length(G,v) for k,v in starts.items()}
    branches,_=swc_branches(swc,np.eye(4),mapping);rows=[];covered=set()
    for b in branches:
        ids=set(b['swc_ids']);covered|=ids;reach=[k for k,d in distances.items() if ids<=d.keys()]
        rows.append(dict(branch_id=b['branch_id'],swc_ids=b['swc_ids'],reachable_ports=reach,minimum_graph_hops=min((min(distances[k][i] for i in ids) for k in reach),default=None),status='REACHABLE' if reach else 'SOLVENT_ACCESS_BLOCKED'))
    for node in sorted(set(G)-covered):
        reach=[k for k,d in distances.items() if node in d];rows.append(dict(branch_id='UNASSIGNED_NODE_'+str(node),swc_ids=[node],reachable_ports=reach,minimum_graph_hops=min((distances[k][node] for k in reach),default=None),status='REACHABLE' if reach else 'SOLVENT_ACCESS_BLOCKED'))
    return dict(all_vascular_branches_reachable=bool(rows) and all(r['reachable_ports'] for r in rows),every_branch_multiple_ports=bool(rows) and all(len(r['reachable_ports'])>1 for r in rows),component_count=nx.number_connected_components(G),node_count=len(G),branch_count=len(branches),rows=rows)


def solvent(l,g,cfg,inputs,out):
    out.mkdir(parents=True,exist_ok=True);G=nx.DiGraph();G.add_node('external_bath',kind='BATH');rows=[];tol=cfg['intersection_volume_tolerance_mm3'];corem=manifold(l['core']);cavity=g['cavity'];core_components=len(l['core'].split());target='vascular_component_0';G.add_node(target,kind='VASCULAR_CORE',mesh_components=core_components)
    exterior_access=bool(g['bottom'].Volume()>0 and g['shell'].intersect(cavity).Volume()<tol)
    G.add_node('outer_bottom',kind='BOTTOM')
    if exterior_access:G.add_edge('external_bath','outer_bottom',evidence='Closed bottom exterior is outside PDMS casting domain and directly bath exposed')
    # Every edge is conditional on actual geometry, not material identity.
    for p in l['ports']:
        k=p['port_id'];record=next(r for r in g['access_records'] if r['port_id']==k);wall='wall_'+k;fixture='fixture_'+k;land='access_land_'+k;ret='retainer_cavity_'+k;port='port_'+k
        for name,kind in [(wall,'WALL_REGION'),(fixture,'FIXTURE_BLADE'),(land,'ACCESS_LAND'),(ret,'RETAINER_CAVITY'),(port,'PORT_EXTENSION')]:G.add_node(name,kind=kind,port_id=k)
        root=float(g['fixture_bases'][k].intersect(g['walls']).Volume());basejoin=float(g['lands'][k].intersect(g['fixture_bases'][k]).Volume());witness=np.asarray(record['contact_witnesses']);_,dist,_=l['core'].nearest.on_surface(witness);gap=float(np.max(dist));overlap=float((cad_manifold(g['lands'][k])^corem).volume());z=l['high'][2];retmesh=shape_mesh(g['retainers'][k]);exposed=float(g['retainers'][k].intersect(bounds_box(np.r_[l['low'][:2],z],l['high']+[0,0,5])).Volume())
        G.add_edge('external_bath',wall,evidence='Uncut continuous exterior wall face, directly bath exposed')
        if root>tol:G.add_edge(wall,fixture,attached_volume_mm3=root)
        if basejoin>tol:G.add_edge(fixture,land,attached_volume_mm3=basejoin)
        contact=contact_surface_evidence(l['core'],shape_mesh(g['lands'][k]),witness,record['contact_area_mm2'],cfg);gap=contact['maximum_contact_surface_error_mm']
        surface_ok=contact['passed'] and overlap<tol and all(contact_zone(w,p) for w in witness)
        if surface_ok:G.add_edge(land,port,finite_contact_area_mm2=record['contact_area_mm2'],maximum_frozen_surface_error_mm=gap,intersection_mm3=overlap,assumption='Nominal dry surface contact; no intervening PDMS membrane. Physical print/vapor/assembly validation pending.')
        if exposed>tol:G.add_edge('external_bath',ret,exposed_stem_volume_above_fill_mm3=exposed)
        if core_components==1:G.add_edge(port,target,evidence='Frozen core is one closed connected mesh; SWC topology separately checked')
        route=nx.shortest_path(G,'external_bath',port) if nx.has_path(G,'external_bath',port) else []
        rows.append(dict(port_id=k,route=route,wall_fixture_attachment_mm3=root,fixture_land_attachment_mm3=basejoin,contact_area_mm2=record['contact_area_mm2'],contact_surface_error_mm=gap,interference_mm3=overlap,retainer_exposed_volume_mm3=exposed,retainer_reachable=nx.has_path(G,'external_bath',ret),status='REACHABLE' if route else 'SOLVENT_ACCESS_BLOCKED'))
    status=graph_status(G,['port_'+p['port_id'] for p in l['ports']]);accessible=[r['port_id'] for r in rows if r['route']];net=network_access(inputs['smoothed_swc'],l['ports'],accessible,inputs['branch_mapping']);comp=network_access(inputs['compensated_swc'],l['ports'],accessible,inputs['branch_mapping'])
    independent=len(list(nx.node_disjoint_paths(G,'external_bath',target))) if nx.has_path(G,'external_bath',target) else 0
    graph=dict(nodes=[dict(id=k,**v) for k,v in G.nodes(data=True)],edges=[dict(source=a,target=b,**v) for a,b,v in G.edges(data=True)],scope='TOPOLOGICAL_AND_NOMINAL_GEOMETRIC_ONLY',dissolution_time_predicted=False,physical_contact_verified=False)
    result=dict(status=status['status'],external_mold_access=exterior_access,port_route_count=len(accessible),independent_access_routes=independent,all_ports_accessible=len(accessible)==len(l['ports']),all_retainers_accessible=all(r['retainer_reachable'] for r in rows),all_vascular_branches_reachable=net['all_vascular_branches_reachable'] and comp['all_vascular_branches_reachable'],every_branch_multiple_ports=net['every_branch_multiple_ports'],network=net,compensated_network=comp,port_paths=rows,geometry_only=True,physical_contact_verified=False)
    if not(result['all_ports_accessible'] and result['all_retainers_accessible'] and result['all_vascular_branches_reachable']):result['status']='SOLVENT_ACCESS_BLOCKED'
    save(out/'acetone_access_graph.json',graph);save(out/'solvent_access_summary.json',result);table(out/'port_access_paths.csv',rows);table(out/'vascular_branch_access.csv',net['rows']);table(out.parent/'tables/vascular_solvent_access.csv',net['rows']);return result


def trim(l,g,cfg,out):
    out.mkdir(parents=True,exist_ok=True);diag=out.parent/'diagnostic_PDMS';diag.mkdir(exist_ok=True);lo=l['low'];hi=l['high'];cav=g['cavity'];cvm=cad_manifold(cav);core=manifold(l['core']);central=manifold(l['central']);parts=cad_manifold(g['mold'])
    for part in g['retainers'].values():parts=parts+cad_manifold(part)
    cuts={};planes=[];records=[];rows=[];wrong=[]
    for p in l['ports']:
        k=p['port_id'];half=cfg['trim_notch_half_width_mm'];z0=lo[2]-p['trim_origin'][2];z1=hi[2]-p['trim_origin'][2];w=p['wall_distance_from_trim_mm']
        cuts[k]=world(bounds_box([-half,0,z0-1],[half,w+3,z1+1]),p).intersect(cav)
        patches=[('PORT_EXPOSURE',[-half,-1e-5,z0],[half,1e-5,z1],p['trim_origin'],p['axis'])]
        mat=port_transform(p)
        for sign in [-1,1]:patches.append(('NOTCH_SIDE',[sign*half-1e-5,0,z0],[sign*half+1e-5,w+.01,z1],p['trim_origin']+mat[:3,0]*sign*half,mat[:3,0]*sign))
        for kind,a,b,origin,normal in patches:
            patch=world(bounds_box(a,b),p).intersect(cav)
            if patch.Volume()>1e-7:planes.append(shape_mesh(patch));records.append(dict(port=k,kind=kind,origin=np.asarray(origin).tolist(),normal=np.asarray(normal).tolist()))
    cutm=cad_manifold(fuse(list(cuts.values())));retained=cvm-cutm;aux=(parts^cvm)-core;raw_before=cvm-core-parts
    # Sew the zero-gap contact for the DIAGNOSTIC Boolean only. Add each
    # land to the core first, with its roof 0.0001 mm inside the already
    # occupied vascular volume. This removes coincident face ambiguity.
    # No production CAD/STL is altered. Added contact material outside the
    # original occupied region must be < 1e-6 mm3. Separately assembled CAD
    # operands have slightly different tessellations: bound the whole
    # diagnostic symmetric difference to 0.01 mm3 (not a production offset).
    contact_epsilon=1e-4;occupied=core;contact_added=0.
    for land in g['lands'].values():
        shifted=cad_manifold(land.translate((0,0,contact_epsilon)))
        contact_added+=abs(float((shifted-core-cad_manifold(land)).volume()));occupied=occupied+shifted
    occupied=occupied+cad_manifold(g['mold_without_lands'])
    for part in g['retainers'].values():occupied=occupied+cad_manifold(part)
    before=cvm-occupied
    contact_delta=abs(float((raw_before-before).volume()))+abs(float((before-raw_before).volume()))
    if contact_added>1e-6 or contact_delta>.01:raise ValueError('DIAGNOSTIC_CONTACT_SEWING_CHANGED_PDMS_DOMAIN')
    after=before-cutm;ideal=retained-core
    residual=float((aux^retained).volume());missing=float((ideal-after).volume());extra=float((after-ideal).volume());central_removed=float((cutm^central).volume())
    finalmesh=mesh(after)
    for p in l['ports']:
        k=p['port_id'];pt=p['trim_origin'];n=p['axis'];pm=l['source_ports'][k];section=pm.section(plane_origin=pt,plane_normal=n)
        if section is None:raise ValueError('TRIM_PLANE_MISSES_EXTENSION_'+k)
        mat=port_transform(p);coords=(section.discrete[0]-pt)@mat[:3,:3];area=float(Polygon(coords[:,[0,2]]).area);patches=[m for m,r in zip(planes,records) if r['port']==k];distances=[fcl_distance(l['central'],m)[0] for m in patches]
        for other,om in l['source_ports'].items():
            if other!=k:
                distances.extend(fcl_distance(om,m)[0] for m in patches);wrong.append(float((cad_manifold(cuts[k])^manifold(om)).volume()))
        angles=np.linspace(0,2*np.pi,24,endpoint=False);ring=pt+(p['radius_mm']+.08)*(np.cos(angles)[:,None]*mat[:3,0]+np.sin(angles)[:,None]*mat[:3,2]);inside=finalmesh.contains(ring-n*.025);outside=finalmesh.contains(ring+n*.025)
        exposed=bool(inside.all() and not outside.any() and len(section.discrete)==1 and (cad_manifold(cuts[k])^manifold(pm)).volume()>0)
        rows.append(dict(port_id=k,origin=pt.tolist(),normal=n.tolist(),axis_plane_angle_degrees=90.,section_area_mm2=area,area_equivalent_diameter_mm=2*np.sqrt(area/np.pi),usable_straight_length_mm=p['usable_straight_length_after_trim_mm'],distance_to_curvature_mm=p['distance_to_curvature_mm'],distance_to_bifurcation_mm=p['distance_to_bifurcation_mm'],minimum_nonport_clearance_mm=min(distances),port_exposed=exposed))
    objs=dict(PDMS_as_cast=before,PDMS_after_all_ABS_removed=before,PDMS_after_trim=after,auxiliary_voids_before_trim=aux,all_ABS_voids_inside_cavity=occupied^cvm,vascular_lumen_after_trim=core^retained,sacrificial_trim_regions=cutm,retained_region=retained);qc=[];meshes={}
    for name,obj in objs.items():
        m=mesh(obj);path=diag/(name+'.stl');m.export(path);read=trimesh.load_mesh(path);removed=0
        if not read.is_watertight:
            count=len(read.faces);read.process(validate=True);removed=count-len(read.faces);read.export(path);read=trimesh.load_mesh(path)
        # A touching cavity/wall Boolean can export a closed, zero-thickness
        # double-sided sheet. It has no volume but falsely paints whole walls
        # as "void" in review. Remove only independently planar, zero-volume
        # components; never discard negative-volume lumen boundaries.
        kept=[];zero_sheets=0
        for component in read.split(only_watertight=False):
            tri=component.triangles
            signed=float(np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6)
            if np.min(component.extents)<1e-7 and abs(signed)<1e-8:zero_sheets+=1
            else:kept.append(component)
        if zero_sheets:
            read=trimesh.util.concatenate(kept);read.export(path);read=trimesh.load_mesh(path)
        err=abs(float(read.volume)-obj.volume());qc.append(dict(name=name,watertight=bool(read.is_watertight),winding_consistent=bool(read.is_winding_consistent),removed_collapsed_export_faces=removed,removed_zero_thickness_components=zero_sheets,roundtrip_volume_error_mm3=err,passed=bool(read.is_watertight and read.is_winding_consistent and err<.05)));meshes[name]=read
    from .casting_fixture_review import poly
    poly(trimesh.util.concatenate(planes)).save(out/'trim_planes.vtp');trimesh.util.concatenate([shape_mesh(x) for x in g['guides'].values()]).export(out/'trim_guide_features.stl')
    result=dict(auxiliary_void_volume_before_trim=float(aux.volume()),auxiliary_void_volume_after_trim=residual,ghost_missing_volume_mm3=missing,ghost_extra_volume_mm3=extra,central_vascular_removed_mm3=central_removed,core_outside_cavity_mm3=float((core-cvm).volume()),maximum_unrelated_port_cut_mm3=max(wrong,default=0.),all_ports_exposed=all(r['port_exposed'] for r in rows),minimum_nonport_clearance=min(r['minimum_nonport_clearance_mm'] for r in rows),usable_port_lengths={r['port_id']:r['usable_straight_length_mm'] for r in rows},pdms_before_volume_mm3=float(before.volume()),pdms_after_volume_mm3=float(after.volume()),pdms_final_material_components=len(after.decompose()),actual_cut_plane_count=len(records),actual_cut_planes=records,diagnostic_mesh_qc=qc,diagnostic_mesh_qc_pass=all(x['passed'] for x in qc),cad_snap_mm=1e-8,frozen_core_snapped=False,as_cast_and_after_dissolution_same_PDMS_domain=True,retainer_stems_above_fill_not_counted_as_PDMS_void=True,guide_depth_mm=cfg['trim_guide_depth_mm'],guide_method='Shallow inner-bottom / inner-wall mold relief transferred to PDMS; guides dissolve with mold, transferred marks remain until trim.',physical_trim_validated=False)
    result.update(diagnostic_contact_sewing_epsilon_mm=contact_epsilon,diagnostic_contact_added_outside_occupied_domain_mm3=contact_added,diagnostic_sewing_PDMS_symmetric_difference_mm3=contact_delta,diagnostic_PDMS_difference_limit_mm3=.01,production_geometry_changed_for_sewing=False)
    save(out/'post_dissolution_trim.json',result);table(out/'port_trim_sections.csv',rows);meshes['planes']=planes;return result,rows,meshes
