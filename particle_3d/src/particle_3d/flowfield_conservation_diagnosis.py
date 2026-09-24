"""Read-only P1 conservation diagnostics. Not imported by particle production.

All fluxes are signed. Slab outward balance is Q_section + Q_wall - Q_in = D.
Tetra gradients and affine polygon fluxes are exact up to floating point;
clipped tetra volumes use convex-hull determinants in translated/scaled space.
No field modification, geometric smoothing, finite differences, or quadrature
approximation of divergence is used.
"""
from itertools import combinations
import numpy as np
from scipy.spatial import ConvexHull
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

EDGES = tuple(combinations(range(4), 2))
FACES = np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])


def p1_divergence(points, tetra, velocity):
    x, u = np.asarray(points,float)[tetra], np.asarray(velocity,float)[tetra]
    # Solve edge equations independently of FrozenFEMField's inverse/einsum.
    dx, du = x[:,1:]-x[:,:1], u[:,1:]-u[:,:1]
    gradient_transpose = np.linalg.solve(dx, du)
    return np.trace(gradient_transpose, axis1=1, axis2=2), np.abs(np.linalg.det(dx))/6


def clip_polygon(x, u, center, normal):
    """Sutherland-Hodgman clipping to dot(x-center,normal)<=0, affine u."""
    d = (x-center)@normal
    px, pu = [], []
    for i in range(len(x)):
        j = (i+1) % len(x)
        if d[i] <= 0:
            px.append(x[i]); pu.append(u[i])
        if (d[i] < 0 < d[j]) or (d[j] < 0 < d[i]):
            t = d[i]/(d[i]-d[j])
            px.append(x[i]+t*(x[j]-x[i])); pu.append(u[i]+t*(u[j]-u[i]))
    return np.asarray(px).reshape(-1,3), np.asarray(pu).reshape(-1,3)


def polygon_flux(x, u, normal=None):
    """Signed integral for affine velocity on an oriented convex polygon."""
    if len(x) < 3:
        return 0., 0.
    av = np.cross(x[1:-1]-x[0], x[2:]-x[0])/2
    if normal is not None:
        av = np.linalg.norm(av,axis=1)[:,None]*normal
    q = np.einsum('ij,ij->', av, (u[0]+u[1:-1]+u[2:])/3)
    return float(q), float(np.linalg.norm(av,axis=1).sum())


def tetra_plane_polygon(x, u, center, normal):
    d = (x-center)@normal
    px, pu = [], []
    for i in np.flatnonzero(d == 0):
        px.append(x[i]); pu.append(u[i])
    for i,j in EDGES:
        if d[i]*d[j] < 0:
            t = d[i]/(d[i]-d[j])
            px.append(x[i]+t*(x[j]-x[i]));pu.append(u[i]+t*(u[j]-u[i]))
    if len(px)<3:
        return np.empty((0,3)),np.empty((0,3))
    px, pu = np.asarray(px), np.asarray(pu)
    e1=np.cross(normal,np.eye(3)[np.argmin(np.abs(normal))]);e1/=np.linalg.norm(e1)
    e2=np.cross(normal,e1);z=px-px.mean(axis=0)
    order=np.argsort(np.arctan2(z@e2,z@e1))
    return px[order],pu[order]


def clipped_tetra_volume(x, center, normal):
    d=(x-center)@normal
    if np.all(d<=0):return float(abs(np.linalg.det(x[1:]-x[0]))/6)
    if np.all(d>=0):return 0.
    vertices=list(x[d<=0])
    for i,j in EDGES:
        if d[i]*d[j]<0:
            t=d[i]/(d[i]-d[j]);vertices.append(x[i]+t*(x[j]-x[i]))
    v=np.asarray(vertices);scale=float(np.max(np.ptp(x,axis=0)))
    # Translation/scaling only conditions the volume calculation, not geometry.
    return float(ConvexHull((v-x[0])/scale).volume*scale**3)


def face_topology(tetra):
    faces=tetra[:,FACES].reshape(-1,3)
    keys=np.sort(faces,axis=1)
    _, first, inverse, counts=np.unique(keys,axis=0,return_index=True,return_inverse=True,return_counts=True)
    if counts.max()>2:raise ValueError('Nonmanifold tetra face')
    order=np.argsort(inverse,kind='stable');starts=np.r_[0,np.cumsum(counts)[:-1]]
    internal=counts==2
    left=order[starts[internal]]//4;right=order[starts[internal]+1]//4
    boundary=counts==1
    return dict(internal_faces=faces[first[internal]],left=left,right=right,
                boundary_faces=faces[first[boundary]],boundary_owners=first[boundary]//4)


def boundary_owners(topology, triangles):
    mapping={tuple(sorted(f)):int(t) for f,t in zip(topology['boundary_faces'],topology['boundary_owners'])}
    return np.array([mapping[tuple(sorted(f))] for f in triangles],dtype=int)


def upstream_cells(points,tetra,topology,center,normal,seeds):
    d=(points-center)@normal
    active=np.min(d[tetra],axis=1)<0
    valid=np.min(d[topology['internal_faces']],axis=1)<0
    left,right=topology['left'][valid],topology['right'][valid]
    adj=coo_matrix((np.ones(2*len(left)),(np.r_[left,right],np.r_[right,left])),shape=(len(tetra),len(tetra))).tocsr()
    _,labels=connected_components(adj,directed=False)
    seeds=np.asarray(seeds)[active[seeds]]
    if not len(seeds) or len(np.unique(labels[seeds]))!=1:
        raise ValueError('Inlet does not seed exactly one connected upstream slab')
    mask=active & (labels==labels[seeds[0]])
    return mask,d


def clipped_boundary_flux(points,velocity,triangles,owners,centroids,selected,center,normal):
    q,area,abs_q,max_vn,count=0.,0.,0.,0.,0
    for tri,owner in zip(triangles,owners):
        if not selected[owner]:continue
        x=points[tri];u=velocity[tri]
        av=np.cross(x[1]-x[0],x[2]-x[0])
        if av@(x.mean(axis=0)-centroids[owner])<0:
            x=x[[0,2,1]];u=u[[0,2,1]];av=-av
        n=av/np.linalg.norm(av)
        x,u=clip_polygon(x,u,center,normal)
        if len(x)<3:continue
        value,a=polygon_flux(x,u)
        q+=value;area+=a;abs_q+=abs(value);count+=1
        max_vn=max(max_vn,float(np.max(np.abs(u@n))))
    return dict(Q_m3_s=q,area_m2=area,sum_absolute_triangle_Q_m3_s=abs_q,max_abs_normal_velocity_m_s=max_vn,triangle_count=count)


def slab_balance(points,tetra,velocity,divergence,volumes,topology,boundaries,center,normal):
    """Actual inlet-connected tetra intersection with upstream halfspace.

    Other disconnected vessels below the same infinite plane are excluded by
    face connectivity. All original boundary roles are integrated separately,
    so unwanted outlet intersections cannot silently be called WALL.
    """
    mask,d=upstream_cells(points,tetra,topology,center,normal,boundaries['INLET'][1])
    maxd=np.max(d[tetra],axis=1)
    full=np.flatnonzero(mask & (maxd<=0));partial=np.flatnonzero(mask & (maxd>0))
    weights=np.zeros(len(tetra));weights[full]=volumes[full]
    q_plane,area=0.,0.
    for t in partial:
        x=points[tetra[t]];u=velocity[tetra[t]]
        weights[t]=clipped_tetra_volume(x,center,normal)
        polygon,values=tetra_plane_polygon(x,u,center,normal)
        q,a=polygon_flux(polygon,values,normal);q_plane+=q;area+=a
    centroids=points[tetra].mean(axis=1)
    surface={role:clipped_boundary_flux(points,velocity,tri,owner,centroids,mask,center,normal)
             for role,(tri,owner) in boundaries.items()}
    D=float(weights@divergence)
    total=q_plane+sum(s['Q_m3_s'] for s in surface.values())
    return dict(Q_section_numpy_m3_s=q_plane,section_area_numpy_m2=area,
                volume_integral_divergence_m3_s=D,closure_residual_m3_s=total-D,
                clipped_volume_m3=float(weights.sum()),full_tetra_count=len(full),partial_tetra_count=len(partial),
                surface_fluxes=surface),weights


def nearest_root_arclength(centroids,root_points):
    best=np.full(len(centroids),np.inf);arc=np.zeros(len(centroids));acc=0.
    for a,b in zip(root_points[:-1],root_points[1:]):
        edge=b-a;length=np.linalg.norm(edge)
        t=np.clip((centroids-a)@edge/length**2,0,1)
        distance=np.linalg.norm(centroids-a-t[:,None]*edge,axis=1)
        pick=distance<best;best[pick]=distance[pick];arc[pick]=acc+t[pick]*length;acc+=length
    return arc,best


def stats(values,weights):
    values,weights=np.asarray(values),np.asarray(weights)
    if not len(values) or weights.sum()<=0:return dict(count=0)
    q=np.quantile(values,[.01,.05,.5,.95,.99])
    return dict(count=len(values),volume_m3=float(weights.sum()),min=float(values.min()),max=float(values.max()),
                mean=float(values.mean()),volume_weighted_mean=float(weights@values/weights.sum()),
                median=float(q[2]),P01=float(q[0]),P05=float(q[1]),P95=float(q[3]),P99=float(q[4]),
                RMS=float(np.sqrt(np.mean(values**2))),volume_weighted_RMS=float(np.sqrt(weights@(values**2)/weights.sum())))
