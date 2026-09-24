"""Signed Euclidean convex/triangle distance by exact support-feature reduction.

Triangle faces, edges and vertices are all included. Ellipsoid edge queries
reduce to a projected ellipse; vertex queries to the ellipsoid. Every real
stationary root of the ellipsoid secular equation is enumerated (including
singular-axis solutions), so overlap uses minimum translation depth rather
than a fabricated negative center distance. Capsule queries use the exact
segment/triangle Minkowski prism dilated by a sphere. No point-cloud proxy.
"""
from dataclasses import dataclass
from itertools import product
import numpy as np
from scipy.optimize import brentq
from .particle_shapes import Sphere, Ellipsoid, Capsule, EPS, unit, roundoff_length


class ConvexDistanceError(RuntimeError): pass


@dataclass(frozen=True)
class TriangleGap:
    gap_m: float
    wall_point_m: np.ndarray
    particle_point_m: np.ndarray
    normal_inward: np.ndarray
    wall_feature: str
    roundoff_m: float


def triangle_normal(triangle):
    return unit(np.cross(triangle[1]-triangle[0],triangle[2]-triangle[0]))


def triangle_closest_many(point, triangles):
    """Independent face projection plus three clamped finite edge segments."""
    t=np.asarray(triangles,dtype=float); p=np.asarray(point,dtype=float)
    a=t[:,0];e=t[:,1]-a;f=t[:,2]-a;v=p-a
    ee=np.einsum('ij,ij->i',e,e);ef=np.einsum('ij,ij->i',e,f);ff=np.einsum('ij,ij->i',f,f)
    ev=np.einsum('ij,ij->i',e,v);fv=np.einsum('ij,ij->i',f,v);den=ee*ff-ef*ef
    good=den>0
    u=np.divide(ff*ev-ef*fv,den,out=np.zeros_like(den),where=good)
    w=np.divide(ee*fv-ef*ev,den,out=np.zeros_like(den),where=good)
    face=a+u[:,None]*e+w[:,None]*f
    weights=np.column_stack((1-u-w,u,w))
    candidates=[face];bary=[weights];distance=[np.where(good & np.all(weights>=0,axis=1),np.sum((p-face)**2,axis=1),np.inf)]
    for i,j in [(0,1),(1,2),(2,0)]:
        edge=t[:,j]-t[:,i];length2=np.sum(edge*edge,axis=1)
        s=np.clip(np.divide(np.sum((p-t[:,i])*edge,axis=1),length2,out=np.zeros_like(length2),where=length2>0),0,1)
        cp=t[:,i]+s[:,None]*edge; bw=np.zeros((len(t),3));bw[:,i]=1-s;bw[:,j]=s
        candidates.append(cp);bary.append(bw);distance.append(np.sum((p-cp)**2,axis=1))
    which=np.argmin(np.stack(distance,axis=1),axis=1);indices=np.arange(len(t))
    return np.stack(candidates,axis=1)[indices,which],np.stack(bary,axis=1)[indices,which]


def feature(weights,tolerance=128*EPS):
    count=np.count_nonzero(np.asarray(weights)>tolerance)
    return 'VERTEX' if count<=1 else 'EDGE' if count==2 else 'FACE'


def _barycentric(point,tri):
    e=np.column_stack((tri[1]-tri[0],tri[2]-tri[0]))
    tail=np.linalg.lstsq(e,point-tri[0],rcond=None)[0]
    return np.r_[1-tail.sum(),tail]


def _finite_triangle_witness(point,tri,tolerance):
    """Validate a point in physical length units, including skinny triangles.

    A fixed barycentric threshold is not scale invariant: a tiny negative
    weight along a long edge can mean a different length error. Nonnegative
    normalized weights provide an actual finite-triangle witness; accept only
    when its Euclidean discrepancy is within the unchanged roundoff budget.
    This rounds an auxiliary witness, never the particle center.
    """
    weights=_barycentric(point,tri)
    nonnegative=np.maximum(weights,0.)
    nonnegative/=nonnegative.sum()
    return np.linalg.norm(nonnegative@tri-point)<=tolerance,nonnegative


def _stationary_ellipsoid_points(eigenvalues,y,extra_directions=()):
    """All stationary surface points relative to a point y, in principal axes.

Scale-free secular equation F(l)=sum(alpha*y²/(alpha+l)²)-1. Between
active poles F is strictly convex; bracketing its unique derivative root
therefore isolates both possible roots. Exterior intervals each have one.
No optimizer start guess or discretized surface is used.
"""
    alpha=np.asarray(eigenvalues,dtype=float);y=np.array(y,dtype=float,copy=True)
    y[np.abs(y)<=32*EPS*max(1.,np.linalg.norm(y))]=0.
    groups=[]
    for i,a in enumerate(alpha):
        found=next((g for g in groups if abs(a-alpha[g[0]])<=32*EPS*max(a,alpha[g[0]])),None)
        if found is None:groups.append([i])
        else:found.append(i)
    for group in groups:alpha[group]=np.mean(alpha[group])
    active=[group for group in groups if np.any(y[group]!=0)]
    poles=sorted(-alpha[g[0]] for g in active)
    coeff=np.array([alpha[g[0]]*np.dot(y[g],y[g]) for g in active],dtype=np.longdouble)
    av=np.array([alpha[g[0]] for g in active],dtype=np.longdouble)
    def f(l):
        z=av+np.longdouble(l)
        if np.any(z==0):return float('inf')
        return float(np.sum(coeff/z**2)-1)
    def df(l):
        z=av+np.longdouble(l)
        return float(np.sum(-2*coeff/z**3))
    roots=[]
    def root(fn,lo,hi):
        return brentq(fn,lo,hi,xtol=np.finfo(float).tiny,rtol=4*EPS,maxiter=150)
    if poles:
        extent=max(2.,2*np.linalg.norm(y)*np.sqrt(np.max(alpha)),2*np.max(alpha))
        lo=poles[0]-extent;hi=poles[-1]+extent
        while f(lo)>0:lo-=extent;extent*=2
        while f(hi)>0:hi+=extent;extent*=2
        roots.append(root(f,lo,np.nextafter(poles[0],-np.inf)))
        roots.append(root(f,np.nextafter(poles[-1],np.inf),hi))
        for left,right in zip(poles[:-1],poles[1:]):
            a=np.nextafter(left,np.inf);b=np.nextafter(right,-np.inf)
            middle=root(df,a,b);minimum=f(middle)
            if minimum<0:
                roots.extend([root(f,a,middle),root(f,middle,b)])
            elif abs(minimum)<128*EPS:roots.append(middle)
    points=[]
    for l in roots:
        denominator=alpha+l
        if np.any(denominator==0):continue
        z=alpha*y/denominator
        residual=abs(np.sum(z*z/alpha)-1)
        if residual<=2048*EPS:
            points.append(z)
    # At an inactive pole, the coordinates in its eigenspace are free.
    for group in groups:
        if np.any(y[group]!=0):continue
        other=[i for i in range(len(alpha)) if i not in group]
        z=np.zeros(len(alpha));z[other]=alpha[other]*y[other]/(alpha[other]-alpha[group[0]])
        remainder=1-np.sum(z[other]**2/alpha[other])
        if remainder < -128*EPS:continue
        radius=np.sqrt(max(0.,remainder)*alpha[group[0]])
        directions=[np.eye(len(alpha))[i] for i in group]
        directions.extend(np.asarray(d,dtype=float) for d in extra_directions)
        if len(group)>1:
            directions.extend(np.array(s,dtype=float) for s in product([-1,0,1],repeat=len(alpha)))
        for d in directions:
            restricted=np.zeros(len(alpha));restricted[group]=d[group]
            length=np.linalg.norm(restricted)
            if length:
                for sign in [-1,1]:points.append(z+sign*radius*restricted/length)
    return points


def _ellipsoid_triangle(shape,tri,tolerance):
    rho=shape.bounding_radius_m
    local=(tri-shape.center_m)/rho
    q=shape.quadratic/rho**2
    normal=triangle_normal(local);candidates=[]
    # A valid positive face support gives an exact global separation certificate.
    for n in [-normal,normal]:
        particle=-(q@n)/np.sqrt(n@q@n)
        gap=float(n@(particle-local[0]));wall=particle-gap*n
        valid,weights=_finite_triangle_witness(wall,local,tolerance/rho)
        if valid:
            candidates.append((gap,wall,particle,n,feature(weights)))
            if gap>tolerance/rho:
                return candidates[-1]
    def add(n,wall_feature,edge=None,vertex=None):
        n=unit(n);particle=-(q@n)/np.sqrt(n@q@n)
        supports=local@n;maximum=float(np.max(supports));gap=float(n@particle-maximum)
        wall=particle-gap*n
        if edge is not None:
            i,j=edge;d=local[j]-local[i];s=float((wall-local[i])@d/(d@d))
            projected=local[i]+np.clip(s,0.,1.)*d
            valid=np.linalg.norm(projected-wall)<=tolerance/rho
        else:valid=np.linalg.norm(wall-local[vertex])<=tolerance/rho
        if valid:
            valid,weights=_finite_triangle_witness(wall,local,tolerance/rho)
            if valid:
                candidates.append((gap,wall,particle,n,feature(weights)))
    for i,j in [(0,1),(1,2),(2,0)]:
        tangent=unit(local[j]-local[i]);basis=np.column_stack((normal,unit(np.cross(tangent,normal))))
        alpha,ev=np.linalg.eigh(basis.T@q@basis);basis=basis@ev;y=basis.T@local[i]
        extra=[basis.T@n for n in [normal,-normal]]
        for z in _stationary_ellipsoid_points(alpha,y,extra):
            add(-basis@(z/alpha),'EDGE',edge=(i,j))
    alpha,ev=np.linalg.eigh(q)
    extra=[ev.T@n for n in [normal,-normal]]+[ev.T@v for v in local]
    for i in range(3):
        y=ev.T@local[i]
        for z in _stationary_ellipsoid_points(alpha,y,extra):
            add(-ev@(z/alpha),'VERTEX',vertex=i)
    if not candidates:raise ConvexDistanceError('No certified ellipsoid/triangle stationary feature; do not fabricate a gap')
    return max(candidates,key=lambda c:c[0])


def _capsule_triangle(shape,tri,tolerance):
    h=shape.cylindrical_length_m/2*shape.axis_world
    if np.linalg.norm(h)<=tolerance:
        return triangle_gap(Sphere(shape.center_m,shape.radius_m),tri)
    # Minkowski prism triangle +/- half shaft, with exact original-vertex provenance.
    vertices=np.vstack((tri-h,tri+h))
    face_ids=np.array([[0,1,2],[3,4,5],[0,1,4],[0,4,3],[1,2,5],[1,5,4],[2,0,3],[2,3,5]])
    triangles=vertices[face_ids];cross=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    valid=np.linalg.norm(cross,axis=1)>tolerance*max(shape.bounding_radius_m,np.max(np.linalg.norm(tri-tri[0],axis=1)))
    triangles=triangles[valid];face_ids=face_ids[valid];cross=cross[valid]
    closest,weights=triangle_closest_many(shape.center_m,triangles)
    distances=np.linalg.norm(closest-shape.center_m,axis=1);k=int(np.argmin(distances));distance=float(distances[k])
    wall=weights[k]@tri[face_ids[k]%3]
    full_rank=abs(np.dot(h,triangle_normal(tri)))>tolerance
    inside=False
    if full_rank:
        normals=cross/np.linalg.norm(cross,axis=1)[:,None]
        prism_center=vertices.mean(axis=0)
        inward=np.einsum('ij,ij->i',normals,prism_center-triangles[:,0])>0
        normals[inward]*=-1
        inside=bool(np.all(np.einsum('ij,ij->i',normals,shape.center_m-triangles[:,0])<=0))
    if distance>tolerance:
        n=(closest[k]-shape.center_m)/distance if inside else (shape.center_m-closest[k])/distance
    else:n=-triangle_normal(tri)
    signed=-distance if inside else distance
    gap=signed-shape.radius_m
    particle=wall+gap*n
    wall_feature=feature(_barycentric(wall,tri))
    nin=n if n@(-triangle_normal(tri))>=0 else -n
    if wall_feature=='FACE':nin=-triangle_normal(tri)
    return TriangleGap(float(gap),wall,particle,nin,wall_feature,tolerance)


def capsule_triangle_many(shape,triangles):
    """Vectorized exact Minkowski-prism queries, same geometry as scalar kernel."""
    tri=np.asarray(triangles);count=len(tri);center=shape.center_m
    tol=roundoff_length(center,tri,shape.bounding_radius_m)
    normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
    if shape.cylindrical_length_m/2<=tol:
        wp,bary=triangle_closest_many(center,tri);d=center-wp;distance=np.linalg.norm(d,axis=1)
        ns=np.divide(d,distance[:,None],out=-normals.copy(),where=distance[:,None]>0)
        on_face=np.all(bary>128*EPS,axis=1)
        ns[on_face]=np.where((np.einsum('ij,ij->i',d,-normals)[on_face]>=0)[:,None],-normals[on_face],normals[on_face])
        gaps=distance-shape.radius_m;pp=wp+gaps[:,None]*ns
        pp[on_face]=center-shape.radius_m*ns[on_face]
    else:
        h=shape.cylindrical_length_m/2*shape.axis_world
        vertices=np.concatenate((tri-h,tri+h),axis=1)
        ids=np.array([[0,1,2],[3,4,5],[0,1,4],[0,4,3],[1,2,5],[1,5,4],[2,0,3],[2,3,5]])
        faces=vertices[:,ids].reshape(-1,3,3)
        cross=np.cross(faces[:,1]-faces[:,0],faces[:,2]-faces[:,0]);size=np.linalg.norm(cross,axis=1)
        valid=size>tol*max(shape.bounding_radius_m,float(np.max(np.linalg.norm(tri-tri[:,:1],axis=2))))
        closest,weights=triangle_closest_many(center,faces)
        distance=np.linalg.norm(closest-center,axis=1).reshape(count,8)
        distance[~valid.reshape(count,8)]=np.inf
        which=np.argmin(distance,axis=1);index=np.arange(count)*8+which
        point=closest[index];dist=distance[np.arange(count),which];wt=weights[index]
        wp=np.einsum('ni,nij->nj',wt,tri[np.arange(count)[:,None],ids[which]%3])
        bary=np.zeros((count,3))
        for j in range(3):
            for a in range(3):bary[:,a]+=np.where(ids[which,j]%3==a,wt[:,j],0.)
        fn=np.divide(cross,size[:,None],out=np.zeros_like(cross),where=size[:,None]>0)
        prism_center=np.repeat(tri.mean(axis=1),8,axis=0)
        fn[np.einsum('ij,ij->i',fn,prism_center-faces[:,0])>0]*=-1
        side=np.einsum('ij,ij->i',fn,center-faces[:,0]).reshape(count,8)
        side[~valid.reshape(count,8)]=-np.inf
        inside=(np.abs(normals@h)>tol)&np.all(side<=0,axis=1)
        delta=center-point;delta[inside]*=-1
        ns=np.divide(delta,dist[:,None],out=-normals.copy(),where=dist[:,None]>tol)
        gaps=np.where(inside,-dist,dist)-shape.radius_m
        pp=wp+gaps[:,None]*ns
    ns[np.einsum('ij,ij->i',ns,-normals)<0]*=-1
    ns[np.all(bary>128*EPS,axis=1)]=-normals[np.all(bary>128*EPS,axis=1)]
    return gaps,wp,pp,ns,bary


def triangle_gap(shape,triangle):
    tri=np.asarray(triangle,dtype=float)
    if tri.shape!=(3,3) or not np.isfinite(tri).all():raise ValueError('finite triangle (3,3) required')
    inward=-triangle_normal(tri)
    tolerance=roundoff_length(tri,shape.center_m,shape.bounding_radius_m)
    if isinstance(shape,Sphere):
        wall,bary=triangle_closest_many(shape.center_m,tri[None,:,:]);wall=wall[0]
        delta=shape.center_m-wall;distance=float(np.linalg.norm(delta));n=delta/distance if distance else inward
        if feature(bary[0])=='FACE':n=inward if delta@inward>=0 else -inward
        particle=shape.center_m-shape.radius_m*n;gap=distance-shape.radius_m
        return TriangleGap(float(gap),wall,particle,n if n@inward>=0 else -n,feature(bary[0]),tolerance)
    if isinstance(shape,Capsule):return _capsule_triangle(shape,tri,tolerance)
    if isinstance(shape,Ellipsoid):
        g,w,p,n,f=_ellipsoid_triangle(shape,tri,tolerance)
        rho=shape.bounding_radius_m
        return TriangleGap(float(g*rho),shape.center_m+w*rho,shape.center_m+p*rho,n if n@inward>=0 else -n,f,tolerance)
    raise TypeError('Sphere, Ellipsoid or Capsule required')
