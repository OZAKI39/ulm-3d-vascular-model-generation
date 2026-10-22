"""Research support trees: binary merging, load radii, exact collision queries.

Independent Python geometry implementation informed by meshy/supportGenerator.js
(MIT, Anton Bovin 2017). Its radius law r + k*sqrt(load) is reused with attribution;
no WebGL/UI or slicer code is incorporated. See third_party/support_sources.json.
"""
import numpy as np
import trimesh

def load_diameter(base_diameter,load,k=.08):return float(base_diameter+2*k*np.sqrt(max(0.,load)))

def strut(a,b,diameter_a,diameter_b=None,sections=24):
    a=np.asarray(a,float);b=np.asarray(b,float);v=b-a;length=np.linalg.norm(v)
    if length<1e-8:raise ValueError('ZERO_SUPPORT_STRUT')
    rb=(diameter_a if diameter_b is None else diameter_b)/2;ra=diameter_a/2
    theta=np.arange(sections)*2*np.pi/sections
    # A capped revolved profile with exact endpoint radii.
    mesh=trimesh.creation.revolve(np.array([[0,0],[ra,0],[rb,length],[0,length]]),sections=sections)
    mat=trimesh.geometry.align_vectors([0,0,1],v/length);mat[:3,3]=a;mesh.apply_transform(mat)
    return mesh

def illegal_collision(support,vascular,box=None,tolerance=1e-6):
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular)
    if box is not None:manager.add_object('box',box)
    intersects=manager.in_collision_single(support)
    if not intersects:return dict(intersects=False,volume_mm3=0.,passed=True)
    volumes=[]
    for m in [vascular]+([box] if box is not None else []):
        hit=trimesh.boolean.intersection([support,m],engine='manifold');volumes.append(abs(float(hit.volume)))
    return dict(intersects=any(v>tolerance for v in volumes),volume_mm3=sum(volumes),passed=sum(volumes)<=tolerance)

def merge_node(a,b,angle_deg,floor):
    a=np.asarray(a);b=np.asarray(b);xy=(a[:2]+b[:2])/2
    z=min(a[2]-np.linalg.norm(xy-a[:2])/np.tan(np.radians(angle_deg)),b[2]-np.linalg.norm(xy-b[:2])/np.tan(np.radians(angle_deg)))-.3
    return np.r_[xy,z] if z>floor+1 else None

def build(contacts,vascular,box,inner,profile,cfg,base_diameter=2.5,tip_diameter=.6):
    inner=np.asarray(inner);floor=inner[0,2];gap=profile['support_top_gap'];nodes=[];edges=[];pieces=[];failed=[]
    manager=trimesh.collision.CollisionManager();manager.add_object('vessel',vascular);manager.add_object('box',box)
    def addnode(p,kind,load=1,contact=None):
        i=len(nodes);nodes.append(dict(id=i,xyz=np.asarray(p).tolist(),kind=kind,load=float(load),contact_id=contact));return i
    def addedge(i,j,da,db,kind):
        a=np.array(nodes[i]['xyz']);b=np.array(nodes[j]['xyz']);shape=strut(a,b,da,db);pieces.append(shape)
        length=float(np.linalg.norm(b-a));incl=float(np.degrees(np.arccos(abs(b[2]-a[2])/length)))
        edges.append(dict(start=i,end=j,kind=kind,length_mm=length,diameter_start_mm=da,diameter_end_mm=db,inclination_deg=incl,supported_tip_count=nodes[j]['load'],
            estimated_load_mm=nodes[j]['load']*length,slenderness=length/min(da,db),volume_mm3=float(shape.volume)))
    active=[]
    for c in contacts:
        top=np.array(c['xyz'])-[0,0,gap];base=top-[0,0,1.6]
        tip=strut(base,top,min(1.2,base_diameter),tip_diameter)
        if not illegal_collision(tip,vascular,box)['passed']:failed.append(dict(contact_id=c['contact_id'],reason='TIP_COLLISION'));continue
        i=addnode(top,'tip',contact=c['contact_id']);j=addnode(base,'shoulder');addedge(j,i,min(1.2,base_diameter),tip_diameter,'taper');active.append(j)
    # Pair highest nodes; merge only where both arms and future trunk are legal.
    while len(active)>1:
        possibilities=[]
        for n,i in enumerate(active):
            for j in active[n+1:]:
                a=np.array(nodes[i]['xyz']);b=np.array(nodes[j]['xyz']);m=merge_node(a,b,cfg['maximum_branch_inclination_deg'],floor)
                if m is None or np.linalg.norm(a-b)>9:continue
                d=load_diameter(base_diameter,nodes[i]['load']+nodes[j]['load'])
                arm1=strut(m,a,d,1.2);arm2=strut(m,b,d,1.2)
                if manager.in_collision_single(arm1) or manager.in_collision_single(arm2):continue
                possibilities.append((np.linalg.norm(a-m)+np.linalg.norm(b-m),i,j,m,d))
        if not possibilities:break
        _,i,j,m,d=min(possibilities,key=lambda q:(q[0],q[1],q[2]));k=addnode(m,'merge',nodes[i]['load']+nodes[j]['load'])
        addedge(k,i,d,1.2,'branch');addedge(k,j,d,1.2,'branch');active=[v for v in active if v not in (i,j)]+[k]
    for i in active:
        a=np.array(nodes[i]['xyz']);inside=np.all(a[:2]>inner[0,:2]) and np.all(a[:2]<inner[1,:2]);ground=floor+profile.get('support_top_gap',.2) if inside else 0.
        b=np.r_[a[:2],ground];d=load_diameter(base_diameter,nodes[i]['load'])
        pillar=strut(b,a,d,1.2)
        if manager.in_collision_single(pillar):failed.append(dict(node=i,reason='PILLAR_COLLISION'));continue
        j=addnode(b,'pillar_base',nodes[i]['load']);addedge(j,i,d,1.2,'pillar')
    # Union support pieces with EACH OTHER only; vascular/box never operands.
    support=trimesh.boolean.union(pieces,engine='manifold') if pieces else trimesh.Trimesh()
    graph=dict(nodes=nodes,edges=edges,tips=[n['id'] for n in nodes if n['kind']=='tip'],merge_nodes=[n['id'] for n in nodes if n['kind']=='merge'],
        pillars=[e for e in edges if e['kind']=='pillar'],bridges=[],failed=failed,tip_diameter_mm=tip_diameter,base_diameter_mm=base_diameter,
        material_volume_mm3=float(support.volume),total_length_mm=sum(e['length_mm'] for e in edges),production=False)
    return support,graph
