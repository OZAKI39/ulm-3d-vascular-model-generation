"""Continuous straight-segment finite-radius certificates, float64/no fastmath.

Full pair minimum is analytic. WALL certificate first uses the 1-Lipschitz
distance bound, then exact segment-triangle distance for all overlapping expanded
AABBs. Caps are NOT WALL triangles. VTK only supplies conservative spatial search.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def point_segment2(p,a,b):
    d=b-a; den=np.dot(d,d)
    f=max(0.,min(1.,np.dot(p-a,d)/den)) if den>0 else 0.
    v=p-a-f*d
    return np.dot(v,v)


@njit(cache=True)
def point_triangle2(p,a,b,c):
    ab=b-a;ac=c-a;n=np.cross(ab,ac);nn=np.dot(n,n)
    if nn>0:
        q=p-n*(np.dot(p-a,n)/nn)
        v=q-a
        aa=np.dot(ab,ab);bb=np.dot(ab,ac);cc=np.dot(ac,ac)
        determinant=aa*cc-bb*bb
        if determinant>0:
            s=(np.dot(v,ab)*cc-np.dot(v,ac)*bb)/determinant
            t=(np.dot(v,ac)*aa-np.dot(v,ab)*bb)/determinant
            if s>=0 and t>=0 and s+t<=1:
                return np.dot(p-q,p-q)
    return min(point_segment2(p,a,b),point_segment2(p,b,c),point_segment2(p,c,a))


@njit(cache=True)
def segment_segment2(p,q,a,b):
    # Interior stationary point plus all four boundary optima of [0,1]^2.
    u=q-p;v=b-a;w=p-a
    aa=np.dot(u,u);bb=np.dot(u,v);cc=np.dot(v,v);dd=np.dot(u,w);ee=np.dot(v,w)
    result=min(point_segment2(p,a,b),point_segment2(q,a,b),
               point_segment2(a,p,q),point_segment2(b,p,q))
    det=aa*cc-bb*bb
    if det>1e-15*aa*cc:
        s=(bb*ee-cc*dd)/det;t=(aa*ee-bb*dd)/det
        if 0<=s<=1 and 0<=t<=1:
            d=w+s*u-t*v;result=min(result,np.dot(d,d))
    return result


@njit(cache=True)
def segment_triangle2(p,q,a,b,c):
    n=np.cross(b-a,c-a);d=q-p;den=np.dot(n,d)
    if den!=0:
        t=np.dot(n,a-p)/den
        if 0<=t<=1:
            hit=p+t*d
            # barycentric oriented-edge signs include the triangle boundaries.
            if (np.dot(np.cross(b-a,hit-a),n)>=0 and
                np.dot(np.cross(c-b,hit-b),n)>=0 and
                np.dot(np.cross(a-c,hit-c),n)>=0):
                return 0.
    return min(point_triangle2(p,a,b,c),point_triangle2(q,a,b,c),
               segment_segment2(p,q,a,b),segment_segment2(p,q,b,c),segment_segment2(p,q,c,a))


@njit(cache=True)
def triangle_minimum(p,q,triangles):
    result=np.inf
    for tri in triangles:
        result=min(result,segment_triangle2(p,q,tri[0],tri[1],tri[2]))
    return np.sqrt(result)


def pair_sweep(x,y,r):
    ii,jj=np.triu_indices(len(x),1)
    delta=x[jj]-x[ii];change=(y-x)[jj]-(y-x)[ii]
    denom=np.einsum('ij,ij->i',change,change)
    fraction=np.clip(-np.einsum('ij,ij->i',delta,change)/np.maximum(denom,np.finfo(float).tiny),0,1)
    gap=np.linalg.norm(delta+fraction[:,None]*change,axis=1)-r[ii]-r[jj]
    return float(np.min(gap,initial=np.inf))


class WallSweep:
    def __init__(self,geometry,wall_path):
        import pyvista as pv
        import vtk
        self.geometry=geometry
        self.mesh=pv.read(wall_path).triangulate()
        faces=self.mesh.faces.reshape(-1,4)
        if np.any(faces[:,0]!=3) or not np.isfinite(self.mesh.points).all():
            raise ValueError('GEOMETRY_DATA_CORRUPTED')
        self.triangles=np.asarray(self.mesh.points,dtype=np.float64)[faces[:,1:]]
        self.locator=vtk.vtkStaticCellLocator();self.locator.SetDataSet(self.mesh);self.locator.BuildLocator()
        self.ids=vtk.vtkIdList();self.path=str(wall_path)
        self.calls=self.exact_particles=0
        # Compile outside measured scientific advance.
        triangle_minimum(np.zeros(3),np.ones(3),self.triangles[:1])

    def __call__(self,x,y,r,tolerance):
        self.calls+=1
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            return dict(pass_gate=False,reason='NONFINITE_INPUT')
        initial=self.geometry.wall_gap_and_normal(x,r)[0]
        final=self.geometry.wall_gap_and_normal(y,r)[0]
        length=np.linalg.norm(y-x,axis=1)
        bound=initial-length
        queried=0;triangle_count=0
        for i in np.flatnonzero(bound < -tolerance):
            # Closed exterior distance is a conservative lower bound on distance
            # to its WALL subset; close paths require the actual open WALL mesh.
            low=np.minimum(x[i],y[i])-r[i]-tolerance
            high=np.maximum(x[i],y[i])+r[i]+tolerance
            self.locator.FindCellsWithinBounds(np.column_stack([low,high]).ravel(),self.ids)
            cells=np.array([self.ids.GetId(k) for k in range(self.ids.GetNumberOfIds())],int)
            minimum=triangle_minimum(x[i],y[i],self.triangles[cells]) if len(cells) else np.inf
            # Omitted triangles lie outside expanded AABB, distance > R+tol.
            bound[i]=min(minimum-r[i],tolerance)
            queried+=1;triangle_count+=len(cells)
        # Floating-point guard subtracted from the computed certificate, never
        # added to the physical tolerance. SI coordinates here are millimetres.
        roundoff=64*np.finfo(float).eps*max(np.max(abs(x),initial=0),np.max(abs(y),initial=0),r.max())
        pair=pair_sweep(x,y,r)-roundoff
        wall=float(bound.min(initial=np.inf))-roundoff
        self.exact_particles+=queried
        return dict(pass_gate=bool(min(pair,wall)>=-tolerance),
                    minimum_pair_swept_gap_m=pair if np.isfinite(pair) else None,
                    minimum_wall_swept_gap_bound_m=wall if np.isfinite(wall) else None,
                    minimum_wall_endpoint_gap_m=float(final.min()),
                    wall_certificate='Lipschitz lower bound or expanded-AABB segment/triangle exact minimum',
                    wall_triangle_candidates=triangle_count,exact_wall_particles=queried,
                    projection_max_m=0.,roundoff_guard_m=roundoff)
