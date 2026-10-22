"""Actual compact-cap assembly, finite motion, fit envelope and wall routes."""
from __future__ import annotations
import json
import numpy as np
import trimesh,networkx as nx
from scipy.spatial.transform import Rotation
from .support_inputs import save,table
from .casting_fixture_geometry import shape_mesh
from .casting_fixture_v2_geometry import port_transform
from .casting_fixture_v2_audit import manifold,cad_manifold,contact_zone
from .casting_mold_qc import fcl_distance
from .support_removal import swept_path
from .all_abs_audit import solvent as original_solvent


def insertion(l,g,cfg):
    present=[g['mold_mesh'],l['core']];inside=[l['core'],shape_mesh(g['bottom']),*[shape_mesh(x) for k in ['fixtures','guides'] for x in g[k].values()]];rows=[]
    for p in l['ports']:
        k=p['port_id'];part=shape_mesh(g['retainers'][k]);n=p['axis'];slide=cfg['retainer_wallward_slide_mm'];standoff=cfg['wall_contact_standoff_mm'];a=np.eye(4);a[:3,3]=-n*slide+[0,0,l['high'][2]+5];b=np.eye(4);b[:3,3]=-n*slide;c=np.eye(4);c[:3,3]=-n*standoff
        open_path=swept_path(part,present,[a,b,c],tolerance=.004);last=swept_path(part,inside,[c,np.eye(4)],tolerance=.004)
        v=part.vertices;v0=v-n*standoff;halfspace=bool(np.all(np.minimum(v[:,:2].min(0),v0[:,:2].min(0))>=l['low'][:2]-1e-7) and np.all(np.maximum(v[:,:2].max(0),v0[:,:2].max(0))<=l['high'][:2]+1e-7))
        interference=float(g['retainers'][k].intersect(g['mold']).Volume());central=fcl_distance(part,l['central'])[0];passed=bool(open_path['passed'] and last['passed'] and halfspace and abs(interference)<1e-6 and central>0)
        rows.append(dict(port_id=k,passed=passed,sequence=['VERTICAL_FROM_TOP','SHORT_WALLWARD_SLIDE','PLANAR_WALL_DATUM_SEATING'],slide_mm=slide,positive_clearance_path=open_path,final_seating_path_without_walls=last,wall_seating_halfspace_proof=halfspace,wall_seating_distance_mm=standoff,mold_interference_mm3=interference,central_clearance_mm=central,elastic_snap_required=False))
        present.append(part);inside.append(part)
    return dict(passed=all(x['passed'] for x in rows),order=[p['port_id'] for p in l['ports']],rows=rows,obstacles=['FROZEN_CORE','FROZEN_MOLD_AND_CRADLES','PREVIOUSLY_INSTALLED_RETAINERS'],end_wall_contact_is_allowed=True,no_elastic_deformation=True)


def contacts(l,g):
    rows=[]
    for p in l['ports']:
        k=p['port_id'];m=shape_mesh(g['retainers'][k]);gap,point,_=fcl_distance(l['core'],m);tree=fcl_distance(l['central'],m)[0];dist=m.nearest.on_surface(l['bifurcations'])[1];transition=float(m.nearest.on_surface(np.asarray([p['straight_start']]))[1][0]);rec=next(d for d in g['compact_parameters'] if d['port_id']==k)
        rows.append(dict(port_id=k,upper_clearance_baseline_mm=p['anti_float_gap_mm'],upper_clearance_parameter_mm=rec['upper_clearance_mm'],measured_minimum_core_gap_mm=gap,guide_clearance_mm=rec['guide_clearance_mm'],central_tree_gap_mm=tree,central_tree_contact=tree<1e-6,closest_core_point=point.tolist(),approved_extension_zone=contact_zone(point,p),distance_to_nearest_bifurcation_mm=float(min(dist)),distance_to_transition_point_mm=transition,retainer_below_fill=rec['below_nominal_fill'],no_lower_weight_support=rec['upper_cap_only']))
    return dict(rows=rows,central_tree_contacts=sum(r['central_tree_contact'] for r in rows),extension_only_contacts=all(r['approved_extension_zone'] for r in rows),positive_clearance=all(r['measured_minimum_core_gap_mm']>0 and r['upper_clearance_parameter_mm']>=r['upper_clearance_baseline_mm']-1e-8 for r in rows))


def anti_float(l,g):
    cm=trimesh.collision.CollisionManager();cm.add_object('MOLD',g['mold_mesh']);mold=trimesh.collision.CollisionManager();mold.add_object('MOLD',g['mold_mesh']);rows=[]
    for k,s in g['retainers'].items():cm.add_object(k,shape_mesh(s))
    pivot=np.mean([p['trim_origin'] for p in l['ports']],axis=0)
    for name,axis,angle,dz in [('Z+',None,0,1.),('pitch+','y',2,0),('pitch-','y',-2,0),('roll+','x',2,0),('roll-','x',-2,0)]:
        mat=np.eye(4)
        if axis:mat[:3,:3]=Rotation.from_euler(axis,angle,degrees=True).as_matrix()
        mat[:3,3]=pivot-mat[:3,:3]@pivot+[0,0,dz];rows.append(dict(mode=name,degrees=angle,z_mm=dz,blocked=bool(cm.in_collision_single(l['core'],transform=mat))))
    locks=[]
    for p in l['ports']:
        k=p['port_id'];mat=np.eye(4);mat[2,3]=.8;locked=bool(mold.in_collision_single(shape_mesh(g['retainers'][k]),transform=mat));gauge=l['source_ports'][k];one=trimesh.collision.CollisionManager();one.add_object(k,shape_mesh(g['retainers'][k]));mat[2,3]=1.;upper=bool(one.in_collision_single(gauge,transform=mat));locks.append(dict(port_id=k,retainer_up_0p8_blocked=locked,port_up_1mm_blocked=upper))
    return dict(passed=all(r['blocked'] for r in rows) and all(r['retainer_up_0p8_blocked'] and r['port_up_1mm_blocked'] for r in locks),poses=rows,ports=locks,physical_force_validated=False)


def axial_range(l,old,new,cfg):
    results=[];limits=cfg['axial_search_limit_mm'];step=cfg['axial_search_step_mm'];base_mesh=shape_mesh(old['mold_without_lands'])
    for p in l['ports']:
        k=p['port_id'];n=p['axis'];gauge=l['source_ports'][k];values={};line=old['line_width_mm'];front=old['local_parts'][k]['retainer'].BoundingBox().ymin
        for label,g in [('old',old),('new',new)]:
            cm=trimesh.collision.CollisionManager();cm.add_object('LOWER_FIXTURE_WITHOUT_TANGENT_LAND',base_mesh);cm.add_object('RETAINER',shape_mesh(g['retainers'][k]));capfront=front if label=='old' else next(x['cap_front_y_mm'] for x in g['compact_parameters'] if x['port_id']==k);capback=3. if label=='old' else next(x['cap_back_y_mm'] for x in g['compact_parameters'] if x['port_id']==k)
            def valid(d):
                mat=np.eye(4);mat[:3,3]=n*d;gap=float(cm.min_distance_single(gauge,transform=mat));overlap=min(p['straight_length_mm']/2+d,capback)-max(-p['straight_length_mm']/2+d,capfront)
                return gap>2e-5 and overlap>=line
            bounds=[];censored=[]
            for sign in [-1,1]:
                last=0.;stop=None
                for v in np.arange(step,limits+step/2,step):
                    if not valid(sign*v):stop=v;break
                    last=v
                censored.append(stop is None)
                if stop is not None:
                    while stop-last>cfg['axial_bisection_mm']:
                        mid=(last+stop)/2
                        if valid(sign*mid):last=mid
                        else:stop=mid
                bounds.append(sign*last)
            values[label]=dict(min_mm=float(bounds[0]),max_mm=float(bounds[1]),search_censored=censored)
        target=p['axial_float_mm'];required_path={}
        for label,gg in [('old',old),('new',new)]:
            mats=[]
            for d in [-target,0.,target]:
                m=np.eye(4);m[:3,3]=n*d;mats.append(m)
            required_path[label]=swept_path(gauge,[base_mesh,shape_mesh(gg['retainers'][k])],mats,tolerance=.004)
        preserved=values['new']['min_mm']<=-target+1e-3 and values['new']['max_mm']>=target-1e-3;not_worse=values['new']['min_mm']<=values['old']['min_mm']+1e-3 and values['new']['max_mm']>=values['old']['max_mm']-1e-3
        # A +/-target placement-error envelope keeps the solvent land under
        # the original straight portion. Nominal exact contact is not treated
        # as interference by the fit-envelope collision search.
        land=next(x for x in old['access_records'] if x['port_id']==k);a=p['straight_length_mm']/2+land['center_from_trim_mm']-land['length_mm']/2-target;b=p['straight_length_mm']/2+land['center_from_trim_mm']+land['length_mm']/2+target;straight_land=a>=.8-1e-6 and b<=p['straight_length_mm']-.8+1e-6
        results.append(dict(port_id=k,role=p['role'],required_axial_float_mm=target,old_min_mm=values['old']['min_mm'],old_max_mm=values['old']['max_mm'],new_min_mm=values['new']['min_mm'],new_max_mm=values['new']['max_mm'],old_search_censored=values['old']['search_censored'],new_search_censored=values['new']['search_censored'],required_float_preserved=preserved and required_path['new']['passed'],required_float_continuous_paths=required_path,not_worse_than_old=not_worse,land_stays_on_straight_at_required_float=straight_land))
    return dict(passed=all(x['required_float_preserved'] and x['not_worse_than_old'] and x['land_stays_on_straight_at_required_float'] for x in results),rows=results,scope='Relative rigid fit envelope of each frozen port versus its fixture; not independent motion of a connected vascular branch. No production port is moved.',land_policy='Analytic tangent land contact; excludes zero-gap lands from positive-distance search, separately checks straight-span engagement at approved float. Print-specific contact is uncalibrated.')


def wall_routes(l,g,basecfg,inputs,out):
    # Rebuild actual wall/land evidence and both SWC graphs using the original
    # algorithm. A buried compact retainer has no above-fill stem edge.
    original=original_solvent(l,g,basecfg,inputs,out);data=json.loads((out/'acetone_access_graph.json').read_text());G=nx.DiGraph()
    for node in data['nodes']:G.add_node(node['id'],**{k:v for k,v in node.items() if k!='id'})
    for edge in data['edges']:G.add_edge(edge['source'],edge['target'],**{k:v for k,v in edge.items() if k not in ['source','target']})
    rows=[];removal=[]
    for p in l['ports']:
        k=p['port_id'];rec=next(x for x in g['compact_parameters'] if x['port_id']==k);shape=g['retainers'][k];w=p['wall_distance_from_trim_mm'];M=port_transform(p);witness=np.array([[0,w,(rec['guide_bottom_z_mm']+rec['guide_top_z_mm'])/2]])@M[:3,:3].T+M[:3,3]
        err=float(shape_mesh(shape).nearest.on_surface(witness)[1][0]);wallerr=float(g['mold_mesh'].nearest.on_surface(witness)[1][0]);area=float(shape.translate(tuple(p['axis']*1e-4)).intersect(g['walls']).Volume())/1e-4;penetration=float(shape.intersect(g['walls']).Volume());contact=err<1e-5 and wallerr<1e-5 and area>.01 and abs(penetration)<1e-6
        if contact:G.add_edge('wall_'+k,'retainer_cavity_'+k,role='REMOVAL_OF_RETAINER_ITSELF_ONLY',contact_area_mm2=area,surface_error_mm=max(err,wallerr),interference_mm3=penetration,nominal_wall_seating_contact=True)
        removal.append(dict(port_id=k,wall_seated_contact=contact,finite_contact_area_mm2=area,retainer_external_route=nx.has_path(G,'external_bath','retainer_cavity_'+k),no_stem=True))
    without=G.copy();without.remove_nodes_from([n for n in G if n.startswith('retainer_')])
    for p in l['ports']:
        k=p['port_id'];path=nx.shortest_path(without,'external_bath','port_'+k) if nx.has_path(without,'external_bath','port_'+k) else [];rows.append(dict(port_id=k,wall_based_solvent_access='PASS' if path else 'PORT_WALL_SOLVENT_ROUTE_FAILED_'+k,route=path,retainer_route_required=False))
    count=len(list(nx.node_disjoint_paths(without,'external_bath','vascular_component_0'))) if nx.has_path(without,'external_bath','vascular_component_0') else 0
    result=dict(retainer_stem_required=False,wall_based_routes=rows,all_ports_wall_access=all(r['route'] for r in rows),independent_port_routes=count,all_vascular_branches_reachable=original['all_vascular_branches_reachable'],every_branch_multiple_ports=original['every_branch_multiple_ports'],network=original['network'],compensated_network=original['compensated_network'],retainer_removal_routes=removal,all_buried_retainers_accessible=all(r['retainer_external_route'] for r in removal),scope='TOPOLOGICAL_AND_GEOMETRIC_ONLY',no_dissolution_time_flow_rate_refresh_prediction=True)
    data['nodes']=[dict(id=k,**v) for k,v in G.nodes(data=True)];data['edges']=[dict(source=a,target=b,**v) for a,b,v in G.edges(data=True)];data['retainer_stem_required']=False
    original.update(status='SOLVENT_ACCESS_GEOMETRIC_PASS' if result['all_ports_wall_access'] and result['all_buried_retainers_accessible'] and result['all_vascular_branches_reachable'] else 'SOLVENT_ACCESS_BLOCKED',all_retainers_accessible=result['all_buried_retainers_accessible'],retainer_stem_required=False)
    for row,rem in zip(original['port_paths'],removal):row.update(retainer_reachable=rem['retainer_external_route'],retainer_removal_route='WALL_DATUM_CONTACT')
    save(out/'solvent_access_summary.json',original);table(out/'port_access_paths.csv',original['port_paths'])
    save(out/'acetone_access_graph.json',data);save(out/'graph_without_retainer_routes.json',dict(nodes=list(without),edges=list(without.edges),port_routes=rows));save(out/'wall_based_solvent_summary.json',result);table(out.parent/'tables/wall_based_solvent_access.csv',rows);return result
