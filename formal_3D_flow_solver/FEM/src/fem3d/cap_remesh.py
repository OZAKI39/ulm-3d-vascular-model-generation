"""Stage 1.5 cap-only geometry validation; no FEM solver imports."""
import numpy as np


def rim_loop(triangles):
    """Extract exactly one simple closed boundary cycle from a triangle patch."""
    triangles=np.asarray(triangles,dtype=np.int64)
    if triangles.ndim!=2 or triangles.shape[1]!=3 or not len(triangles):
        raise ValueError('Invalid cap triangles')
    if np.any(np.diff(np.sort(triangles,axis=1),axis=1)==0):
        raise ValueError('Degenerate cap triangle')
    edges=np.concatenate([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]])
    edges,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
    if np.any(counts>2): raise ValueError('Nonmanifold cap edge')
    boundary=edges[counts==1]
    adjacency={}
    for x,y in boundary:
        adjacency.setdefault(int(x),[]).append(int(y))
        adjacency.setdefault(int(y),[]).append(int(x))
    if not adjacency or any(len(v)!=2 for v in adjacency.values()):
        raise ValueError('Open or branching rim; every vertex must have degree 2')
    first=min(adjacency); loop=[first]; previous=None; current=first
    while True:
        choices=sorted(v for v in adjacency[current] if v!=previous)
        nxt=choices[0]
        if nxt==first: break
        if nxt in loop: raise ValueError('Repeated rim vertex')
        loop.append(nxt); previous,current=current,nxt
    if len(loop)!=len(adjacency): raise ValueError('Multiple rim loops')
    return np.asarray(loop,dtype=np.int64),boundary


def plane_basis(normal):
    """Deterministic right-handed orthonormal basis; no fitted circle/plane."""
    n=np.asarray(normal,dtype=float)
    if n.shape!=(3,) or not np.isfinite(n).all() or np.linalg.norm(n)==0:
        raise ValueError('Invalid source plane normal')
    n=n/np.linalg.norm(n)
    seed=np.eye(3)[np.argmin(np.abs(n))]
    t1=np.cross(n,seed); t1/=np.linalg.norm(t1)
    t2=np.cross(n,t1)
    return t1,t2,n


def project(points,origin,basis):
    t1,t2,n=basis
    delta=np.asarray(points)-np.asarray(origin)
    return delta@np.column_stack([t1,t2]),delta@n


def lift(interior_xy,origin,basis):
    """Lift NEW interior vertices only. Original rim coordinates must be reused."""
    return np.asarray(origin)+np.asarray(interior_xy)@np.asarray(basis[:2])


def triangle_quality(points,triangles):
    p=np.asarray(points)[np.asarray(triangles)]
    edges=np.stack([p[:,1]-p[:,0],p[:,2]-p[:,1],p[:,0]-p[:,2]],axis=1)
    lengths=np.linalg.norm(edges,axis=2)
    area=np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)/2
    if not np.isfinite(area).all() or (area<=0).any(): raise ValueError('Nonpositive cap triangle')
    return 4*np.sqrt(3)*area/np.sum(lengths**2,axis=1),lengths.max(axis=1)/lengths.min(axis=1)


def quality_summary(points,triangles):
    from .mesh_qc import quantiles
    q,r=triangle_quality(points,triangles)
    return {'triangle_count':len(triangles),'q_tri':quantiles(q),'edge_ratio_max_to_min':quantiles(r),
            'low_quality_counts':{str(v):int(np.count_nonzero(q<v)) for v in (.1,.2,.3)}}


def edge_keys(edges):
    return {tuple(e) for e in np.sort(np.asarray(edges),axis=1)}


def check_rim(original_points,original_edges,points,triangles):
    """Persistent original vertex IDs; exact coordinates AND exact edge set."""
    loop,edges=rim_loop(triangles)
    ids=np.unique(original_edges)
    if ids.max()>=len(points) or not np.array_equal(np.asarray(points)[ids],np.asarray(original_points)[ids]):
        raise ValueError('Rim vertex moved or removed')
    if edge_keys(edges)!=edge_keys(original_edges): raise ValueError('Rim edge set changed (split/deleted edge)')
    return {'maximum_rim_displacement_m':0.0,'rim_edge_set_identical':True,'rim_vertex_count':len(loop)}


def check_wall(source,derived):
    old=source['triangles'][source['facet_tags']==1]
    new=derived['triangles'][derived['facet_tags']==1]
    if not np.array_equal(old,new): raise ValueError('Wall connectivity or semantic tag changed')
    ids=np.unique(old)
    if ids.max()>=len(derived['points_m']) or not np.array_equal(source['points_m'][ids],derived['points_m'][ids]):
        raise ValueError('Wall vertex moved')
    if set(derived['facet_tags'])!={1,2,3,4,5}: raise ValueError('Entity mapping changed')
    return {'wall_triangle_count':len(old),'wall_connectivity_identical':True,'maximum_wall_displacement_m':0.0}


def cross2(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def signed_area(polygon):
    return .5*float(np.sum(cross2(polygon,np.roll(polygon,-1,axis=0))))


def convex_intersection_area(subject,clip):
    """Sutherland-Hodgman clipping; both triangles counter-clockwise, scaled 2D."""
    result=list(subject)
    for a,b in zip(clip,np.roll(clip,-1,axis=0)):
        if not result: return 0.0
        output=[]; prev=result[-1]; dp=cross2(b-a,prev-a)
        for cur in result:
            dc=cross2(b-a,cur-a)
            if (dp>=0)!=(dc>=0): output.append(prev+(cur-prev)*(dp/(dp-dc)))
            if dc>=0: output.append(cur)
            prev,dp=cur,dc
        result=output
    return abs(signed_area(np.asarray(result))) if len(result)>=3 else 0.0


def inside_polygon(points,polygon,tol):
    points=np.asarray(points); inside=np.zeros(len(points),bool); on=np.zeros(len(points),bool)
    for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
        edge=b-a; delta=points-a
        on|=(np.abs(cross2(edge,delta))<=tol*np.linalg.norm(edge))&(delta@edge>=-tol)&(delta@edge<=edge@edge+tol)
        hit=(a[1]>points[:,1])!=(b[1]>points[:,1])
        if b[1]!=a[1]: inside^=hit&(points[:,0]<(b[0]-a[0])*(points[:,1]-a[1])/(b[1]-a[1])+a[0])
    return inside|on


def coverage_check(xy,triangles,rim_ids,relative_tolerance=1e-12):
    """Boundary identity + oriented disk topology + clipping/coverage in 2D.

    Coordinates scaled by polygon diameter only for geometric predicates.
    Overlap tolerance and area deficit/excess tolerance are relative to polygon
    area. Strict rim edge identity also rejects holes smaller than this tolerance.
    """
    xy=np.asarray(xy); triangles=np.asarray(triangles)
    scale=float(np.ptp(xy[rim_ids],axis=0).max())
    p=(xy-xy[rim_ids].mean(axis=0))/scale; polygon=p[rim_ids]
    if signed_area(polygon)<0: polygon=polygon[::-1]
    pa=signed_area(polygon); t=p[triangles]
    areas=.5*cross2(t[:,1]-t[:,0],t[:,2]-t[:,0])
    if np.any(areas<=0): raise ValueError('Normal flipped or degenerate projected triangle')
    _,actual_edges=rim_loop(triangles)
    expected=np.c_[rim_ids,np.roll(rim_ids,-1)]
    if edge_keys(actual_edges)!=edge_keys(expected): raise ValueError('Cap hole or changed polygon boundary')
    edges=np.unique(np.sort(np.concatenate([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]]),axis=1),axis=0)
    if len(np.unique(triangles))-len(edges)+len(triangles)!=1: raise ValueError('Cap is not a topological disk')
    samples=np.concatenate([t.mean(axis=1),(p[edges[:,0]]+p[edges[:,1]])/2,p[np.unique(triangles)]])
    tol=relative_tolerance
    if not inside_polygon(samples,polygon,tol).all(): raise ValueError('Triangle outside rim polygon')
    # Reject proper crossings of mesh edges with the polygon boundary.
    for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
        c=p[edges[:,0]]; d=p[edges[:,1]]
        s1=cross2(b-a,c-a); s2=cross2(b-a,d-a)
        s3=cross2(d-c,a-c); s4=cross2(d-c,b-c)
        if np.any((s1*s2 < -tol**2)&(s3*s4 < -tol**2)):
            raise ValueError('Triangle edge crosses outside polygon')
    lower=t.min(axis=1); upper=t.max(axis=1); max_overlap=0.0
    for i in range(len(t)):
        possible=np.flatnonzero(np.all(upper[i]>=lower[i+1:],axis=1)&np.all(upper[i+1:]>=lower[i],axis=1))+i+1
        for j in possible:
            overlap=convex_intersection_area(t[i],t[j])
            max_overlap=max(max_overlap,overlap)
            if overlap>tol*pa: raise ValueError('Cap triangles overlap')
    error=abs(float(areas.sum())-pa)/pa
    if error>tol: raise ValueError('Cap polygon coverage hole or excess')
    return {'status':'PASS','holes':0,'overlaps':0,'triangles_outside_polygon':0,'relative_area_coverage_error':error,'maximum_pair_overlap_area_normalized':max_overlap/pa,'relative_area_tolerance':tol,'method':'Exact rim-edge identity; Euler disk topology; positive orientation; point-in-polygon and proper edge crossings; all bbox-overlapping triangle pairs clipped with Sutherland-Hodgman; total signed area closure','predicate_coordinate_scale_m':scale}


def port_geometry(points,triangles,port,new_interior_ids,policy):
    from .mesh_qc import triangle_geometry
    a,c,v=triangle_geometry(points,triangles)
    total=float(a.sum()); center=np.sum(c*a[:,None],axis=0)/total
    normal=v.sum(axis=0); normal/=np.linalg.norm(normal)
    origin=np.asarray(port['plane_origin_m']); n=np.asarray(port['outward_normal'])
    plane=float(np.abs((points[np.unique(triangles)]-origin)@n).max())
    interior=float(np.abs((points[new_interior_ids]-origin)@n).max(initial=0))
    roundoff=512*np.finfo(float).eps*float(np.max(np.abs(points)))
    rel=abs(total-port['area_m2'])/port['area_m2']; shift=float(np.linalg.norm(center-origin)); dot=float(normal@n)
    g=policy['geometry']
    checks={'area':rel<=g['area_relative_error_max'],'centroid':shift<=g['centroid_displacement_m_max'],'normal':dot>=g['normal_dot_min'],'source_planarity':plane<=port['max_plane_deviation_m']+roundoff,'interior_on_original_plane':interior<=roundoff}
    checks={name:bool(value) for name,value in checks.items()}
    return {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'area_m2':total,'source_area_m2':port['area_m2'],'relative_area_error':rel,'centroid_m':center.tolist(),'centroid_displacement_m':shift,'normal_dot_source':dot,'outward_normal':normal.tolist(),'max_plane_deviation_m':plane,'source_max_plane_deviation_m':port['max_plane_deviation_m'],'new_interior_max_plane_deviation_m':interior,'roundoff_tolerance_m':roundoff}
