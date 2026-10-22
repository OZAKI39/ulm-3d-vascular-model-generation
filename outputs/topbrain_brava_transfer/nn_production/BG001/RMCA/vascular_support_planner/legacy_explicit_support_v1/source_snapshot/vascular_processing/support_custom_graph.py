"""Explicit Clever-style downward graph with bottom-only breakaway roots.

Algorithmic basis: Vanek et al. 2014, section 5 (highest-point cone sweep),
cross-checked with meshy js/supportGenerator.js, MIT, Anton Bovin 2017.
The sqrt(load) radius rule is adapted from meshy; its license is retained under
third_party/support_research/meshy. Collision uses FCL, never custom triangles.
Differences: no mesh/vascular anchoring, fixed casting orientation, exact gaps,
controlled feet and necks, conservative stability limits and rejected-route logs.
"""
from __future__ import annotations
import numpy as np
import trimesh
from .support_tree_planner import strut,load_diameter
from .support_breakaway import root_geometry,bottom_root_allowed

def cone_intersection(a,b,angle_deg=35.,offset=.25):
    a=np.asarray(a,float);b=np.asarray(b,float);delta=b[:2]-a[:2];d=np.linalg.norm(delta);t=np.tan(np.radians(angle_deg))
    if d<1e-8:return np.r_[a[:2],min(a[2],b[2])-offset]
    z=min(a[2],b[2],(a[2]+b[2]-d/t)/2)-offset
    along=np.clip((d+(a[2]-b[2])*t)/2,0,d)
    return np.r_[a[:2]+delta*(along/d),z]

def stability(a,b,diameter,load,cfg):
    a=np.asarray(a);b=np.asarray(b);length=float(np.linalg.norm(b-a));angle=float(np.degrees(np.arccos(np.clip(abs(b[2]-a[2])/max(length,1e-12),0,1))))
    slenderness=length/diameter;failed=slenderness>cfg['maximum_slenderness'] or angle>cfg['maximum_branch_inclination_deg']+.01
    return dict(length_mm=length,diameter_mm=diameter,slenderness=slenderness,angle_to_vertical_deg=angle,supported_load=load,height_mm=float(max(a[2],b[2])),
        status='UNSTABLE_SUPPORT_EDGE' if failed else 'GEOMETRIC_STABILITY_SCREEN_PASSED',physical_strength_validated=False)

def minimum_distance(mesh,other):
    manager=trimesh.collision.CollisionManager();manager.add_object('frozen',other)
    if manager.in_collision_single(mesh):return 0.
    return float(manager.min_distance_single(mesh))

def pillar_strut(a,b,diameter_a,diameter_b):
    """Loft horizontal rings so an offset pillar starts on its whole foot.

    Rotating a conventional frustum tilts its end caps. At the vertical foot
    and tip shoulder that creates unsupported wedges despite a legal axis
    angle. A shear keeps every cross-section parallel to the print layers.
    """
    a=np.asarray(a,float);b=np.asarray(b,float);height=b[2]-a[2]
    if height<=1e-8:raise ValueError('PILLAR_MUST_ASCEND')
    shape=strut([0,0,0],[0,0,height],diameter_a,diameter_b)
    matrix=np.eye(4);matrix[:2,2]=(b[:2]-a[:2])/height;matrix[:3,3]=a
    shape.apply_transform(matrix)
    return shape

def edge_geometry(graph,edge):
    a=graph['nodes'][edge['start']]['xyz'];b=graph['nodes'][edge['end']]['xyz']
    builder=pillar_strut if edge['kind']=='pillar' else strut
    return builder(a,b,edge['diameter_start_mm'],edge['diameter_end_mm'])

def build(contacts,vascular,box,inner,profile,cfg,tip_width,neck_width,pillar_diameter):
    inner=np.asarray(inner);floor=float(inner[0,2]);gap=profile['support_top_gap'];xygap=profile['support_xy_gap']
    taper_height=10*profile['line_width']
    vessels=trimesh.collision.CollisionManager();vessels.add_object('vascular',vascular)
    walls=trimesh.collision.CollisionManager();walls.add_object('box',box)
    def legal(shape,clearance=xygap):
        return not vessels.in_collision_single(shape) and not walls.in_collision_single(shape) and vessels.min_distance_single(shape)>=clearance-.015
    nodes=[];edges=[];pieces=[];tips=[];roots=[];feet=[];rejected=[];repairs=[];active=[]
    def node(xyz,kind,load=1,contact_ids=None):
        i=len(nodes);nodes.append(dict(id=i,xyz=np.asarray(xyz).tolist(),kind=kind,load=float(load),contact_ids=contact_ids or []));return i
    def edge(i,j,da,db,kind):
        builder=pillar_strut if kind=='pillar' else strut
        shape=builder(nodes[i]['xyz'],nodes[j]['xyz'],da,db);pieces.append(shape)
        row=stability(nodes[i]['xyz'],nodes[j]['xyz'],min(da,db),nodes[j]['load'],cfg)
        row.update(edge_id=f'E{len(edges):03}',start=i,end=j,diameter_start_mm=da,diameter_end_mm=db,kind=kind);edges.append(row)
    for c in contacts:
        top=np.array(c['xyz'])-[0,0,gap];shoulder=top-[0,0,taper_height]
        matched=next((p for p in tips if np.linalg.norm(np.array(p['contact_xyz'])-c['xyz'])<.05),None)
        if matched:
            matched['contact_ids'].append(c['contact_id'])
            for identity in [matched['node_id'],matched['shoulder_node_id']]:nodes[identity]['contact_ids'].append(c['contact_id'])
            continue
        tip=strut(shoulder,top,1.1,tip_width);distance=float(vessels.min_distance_single(tip))
        # Correct the finite-width cap against actual triangles; vertical point
        # subtraction alone does not establish the true surface gap.
        previous=0.;shift=0.;original_top=top.copy()
        for _ in range(30):
            if distance>=gap and not vessels.in_collision_single(tip):break
            previous=shift;shift+=.05;top=original_top-[0,0,shift];shoulder=top-[0,0,taper_height]
            tip=strut(shoulder,top,1.1,tip_width);distance=float(vessels.min_distance_single(tip))
        if shift and distance>=gap:
            lower,upper=previous,shift
            for _ in range(8):
                mid=(lower+upper)/2;mtop=original_top-[0,0,mid];trial=strut(mtop-[0,0,taper_height],mtop,1.1,tip_width)
                if not vessels.in_collision_single(trial) and vessels.min_distance_single(trial)>=gap:upper=mid
                else:lower=mid
            top=original_top-[0,0,upper];shoulder=top-[0,0,taper_height];tip=strut(shoulder,top,1.1,tip_width);distance=float(vessels.min_distance_single(tip))
        if not legal(tip,gap) or abs(distance-gap)>.08:
            rejected.append(dict(contact_id=c['contact_id'],point_id=c['point_id'],xyz=c['xyz'],reason='TIP_GAP_OR_COLLISION',measured_gap_mm=distance));continue
        a=node(top,'tip',contact_ids=[c['contact_id']]);b=node(shoulder,'shoulder',contact_ids=[c['contact_id']]);edge(b,a,1.1,tip_width,'taper');active.append(b)
        tips.append(dict(tip_id=f'T{len(tips):03}',node_id=a,shoulder_node_id=b,contact_ids=[c['contact_id']],point_id=c['point_id'],branch_id=c['branch_id'],arc_position=c['arc_position'],
            contact_xyz=c['xyz'],xyz=top.tolist(),width_mm=tip_width,measured_top_gap_mm=distance,distance_to_bifurcation=c['distance_to_bifurcation'],
            status='SUPPORT_NEAR_BIFURCATION_EXCEPTION' if c['distance_to_bifurcation']<cfg['bifurcation_keepout_mm'] else 'SELECTED',
            exception_reason=c.get('keepout_reason'),root_connected=False))
    # Top-down priority from the paper, with deterministic nearest legal cone
    # intersection. Rejected merges fall back to a wider vertical rooted pillar,
    # never to an additional contact on the vessel.
    terminals=[]
    while active:
        i=min(active,key=lambda n:(-nodes[n]['xyz'][2],n));active.remove(i);a=np.array(nodes[i]['xyz']);choices=[]
        for j in active:
            b=np.array(nodes[j]['xyz']);m=cone_intersection(a,b,min(35,cfg['maximum_branch_inclination_deg']))
            if m[2]<floor+4 or np.linalg.norm(a-b)>cfg.get('merge_max_distance_mm',8.):continue
            d=load_diameter(pillar_diameter,nodes[i]['load']+nodes[j]['load'])
            da=max(1.1,load_diameter(pillar_diameter,nodes[i]['load']) if nodes[i]['kind']=='merge' else 1.1)
            db=max(1.1,load_diameter(pillar_diameter,nodes[j]['load']) if nodes[j]['kind']=='merge' else 1.1)
            arms=[strut(m,a,d,da),strut(m,b,d,db)]
            if all(legal(v) for v in arms):choices.append((float(np.linalg.norm(a-m)),j,m,d,da,db))
        if choices:
            _,j,m,d,da,db=min(choices,key=lambda v:(v[0],v[1]));active.remove(j)
            k=node(m,'merge',nodes[i]['load']+nodes[j]['load'],nodes[i]['contact_ids']+nodes[j]['contact_ids']);edge(k,i,d,da,'branch');edge(k,j,d,db,'branch');active.append(k)
        else:terminals.append(i)
    for n,i in enumerate(terminals):
        a=np.array(nodes[i]['xyz']);load=nodes[i]['load'];d=max(load_diameter(pillar_diameter,load),(a[2]-floor)/cfg['maximum_slenderness'])
        root,foot,r=root_geometry(a[:2],floor,d,neck_width,profile,f'R{n:03}',load)
        # A solid side wall separates exterior ports below its top from every
        # allowed inner-bottom root. Do not silently root them on the bed.
        if np.any(a[:2]<inner[0,:2]) or np.any(a[:2]>inner[1,:2]):
            rejected.append(dict(node_id=i,contact_ids=nodes[i]['contact_ids'],xyz=a.tolist(),reason='EXTERIOR_TIP_NO_INNER_BOTTOM_ROUTE',wall_top_mm=float(inner[1,2])));continue
        radius=r['foot_diameter_mm']/2+.08;proposed=[]
        offsets=cfg.get('custom_support',{}).get('root_route_offsets_mm',[0.,1.,2.,4.,6.])
        for offset in offsets:
            for angle in ([0] if offset==0 else np.arange(0,360,45)):
                direction=np.array([np.cos(np.radians(angle)),np.sin(np.radians(angle))])
                proposed.append(np.clip(a[:2]+offset*direction,inner[0,:2]+radius,inner[1,:2]-radius))
        found=False;original_end=load_diameter(pillar_diameter,load) if nodes[i]['kind']=='merge' else 1.1
        for xy in proposed:
            root,foot,r=root_geometry(xy,floor,d,neck_width,profile,f'R{n:03}',load);start=np.r_[xy,r['trunk_start_mm']-.015]
            end_d=max(original_end,float(np.linalg.norm(a-start))/cfg['maximum_slenderness']+.005)
            if stability(start,a,min(d,end_d),load,cfg)['status']=='UNSTABLE_SUPPORT_EDGE':continue
            trunk=pillar_strut(start,a,d,end_d)
            if bottom_root_allowed(r,inner) and legal(trunk) and not vessels.in_collision_single(root):found=True;break
        if not found:
            rejected.append(dict(node_id=i,contact_ids=nodes[i]['contact_ids'],xyz=a.tolist(),reason='BOTTOM_ROOT_OR_TRUNK_COLLISION'));continue
        if end_d>original_end:
            repairs.append(dict(node_id=i,action='WIDEN_LONG_TRUNK_UPPER_END',original_diameter_mm=original_end,repaired_diameter_mm=end_d,added_vascular_tips=0))
        roots.append(r);feet.append(foot);pieces.append(root);k=node(start,'trunk_base',load,nodes[i]['contact_ids']);edge(k,i,d,end_d,'pillar')
        shifted=float(np.linalg.norm(xy-a[:2]));repairs.append(dict(node_id=i,action='OFFSET_PILLAR_COLLISION_REPAIR' if shifted>.1 else 'VERTICAL_PILLAR_FALLBACK',root_shift_mm=shifted,added_vascular_tips=0,trunk_diameter_mm=d,neck_width_mm=neck_width))
        for tip in tips:
            if set(tip['contact_ids'])&set(nodes[i]['contact_ids']):tip['root_connected']=True;tip['root_id']=r['root_id']
    # Only retain geometry that is linked to an approved root. Unrooted tips
    # remain in rejected_routes for coverage diagnostics, not floating STL parts.
    rooted_contacts={q for t in tips if t['root_connected'] for q in t['contact_ids']}
    retained=[];retained_edges=[]
    for e in edges:
        if set(nodes[e['end']]['contact_ids'])&rooted_contacts:
            retained.append(edge_geometry({'nodes':nodes},e));retained_edges.append(e)
    for root in roots:
        m,_,_=root_geometry(root['xyz'][:2],floor,root['trunk_diameter_mm'],neck_width,profile,root['root_id'],root['trunk_load']);retained.append(m)
    model=trimesh.boolean.union(retained,engine='manifold') if retained else trimesh.Trimesh()
    graph=dict(nodes=nodes,edges=retained_edges,tips=tips,roots=roots,merge_nodes=[v['id'] for v in nodes if v['kind']=='merge'],rejected_routes=rejected,stability_repairs=repairs,
        method='CUSTOM_LITERATURE_AWARE_SUPPORT',algorithm='Top-down nearest collision-free cone intersections; MIT meshy load law; bottom-only breakaway roots',
        tip_width_mm=tip_width,neck_width_mm=neck_width,pillar_diameter_mm=pillar_diameter,taper_height_mm=taper_height,physical_strength_validated=False)
    return model,graph,feet
