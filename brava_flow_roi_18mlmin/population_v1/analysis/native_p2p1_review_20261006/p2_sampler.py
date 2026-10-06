"""Validation-only FP64 TET10 velocity/TET4 pressure sampling; no production import."""
import numpy as np
from scipy.spatial import cKDTree

EDGES=np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]],dtype=np.int64)
DL=np.vstack([np.eye(3),-np.ones(3)])

def basis(xi):
    xi=np.asarray(xi,dtype=np.float64)
    L=np.concatenate([xi,1-xi.sum(axis=-1,keepdims=True)],axis=-1)
    N=np.concatenate([L*(2*L-1),4*L[...,EDGES[:,0]]*L[...,EDGES[:,1]]],axis=-1)
    D=np.concatenate([(4*L-1)[...,None]*DL,
         4*(L[...,EDGES[:,0],None]*DL[EDGES[:,1]]+L[...,EDGES[:,1],None]*DL[EDGES[:,0]])],axis=-2)
    return N,D,L

def elevate_affine(points,tets):
    """One arithmetic midpoint per sorted unique edge; unchanged original corners."""
    points=np.asarray(points,dtype=np.float64);tets=np.asarray(tets,dtype=np.int64)
    unique,index=np.unique(np.sort(tets[:,EDGES].reshape(-1,2),axis=1),axis=0,return_inverse=True)
    return np.vstack([points,points[unique].mean(axis=1)]),np.column_stack([tets,len(points)+index.reshape(-1,6)]),unique

def quantities(u,G):
    return dict(u=u,G=G,curl=np.array([G[2,1]-G[1,2],G[0,2]-G[2,0],G[1,0]-G[0,1]]),strain=(G+G.T)/2)

class P2Sampler:
    """Closed-element containment with roundoff tolerance; smallest cell ID wins.

    Curved inverse mapping is supported for fixtures, not thereby certified for CFD.
    Native pressure is evaluated in the shared reference coordinate. The native
    solver's different pressure gradient/Jacobian mapping is audited separately.
    """
    def __init__(self,points,cells,velocity,pressure):
        self.points=np.asarray(points,dtype=np.float64)
        self.cells=np.asarray(cells,dtype=np.int64)
        self.velocity=np.asarray(velocity,dtype=np.float64)
        self.pressure=np.asarray(pressure,dtype=np.float64)
        if self.cells.ndim!=2 or self.cells.shape[1]!=10:raise ValueError('TET10 connectivity required')
        if self.velocity.shape!=(len(self.points),3) or self.pressure.shape!=(len(self.points),):raise ValueError('nodal arrays required')
        if not all(np.isfinite(a).all() for a in (self.points,self.velocity,self.pressure)):raise ValueError('nonfinite nodal input')
        X=self.points[self.cells]
        # Bernstein control points bound curved quadratic geometry, unlike nodal AABBs.
        B=np.concatenate([X[:,:4],2*X[:,4:]-.5*(X[:,EDGES[:,0]]+X[:,EDGES[:,1]])],axis=1)
        self.low=B.min(axis=1);self.high=B.max(axis=1)
        self.center=(self.low+self.high)/2
        self.radius=np.linalg.norm((self.high-self.low)/2,axis=1)
        self.tree=cKDTree(self.center)
        self.scale=np.max(self.high-self.low,axis=1)
        if np.any(self.scale<=0):raise ValueError('degenerate element')
        self.maxradius=float(self.radius.max())

    def in_cell(self,x,cell):
        X=self.points[self.cells[cell]];origin=X[3];Y=X-origin
        J=Y[:3].T
        try:xi=np.linalg.solve(J,x-origin)
        except np.linalg.LinAlgError:return None
        eps=np.finfo(float).eps
        # Coordinate error budget, not physical position projection or relaxed wall clearance.
        tolx=128*eps*max(self.scale[cell],np.max(np.abs(X)),np.max(np.abs(x)))
        converged=False
        for _ in range(24):
            N,D,L=basis(xi);J=Y.T@D
            residual=N@Y-(x-origin)
            if not np.isfinite(J).all() or np.linalg.det(J)<=0:return None
            if np.max(np.abs(residual))<=tolx:
                converged=True;break
            try:xi-=np.linalg.solve(J,residual)
            except np.linalg.LinAlgError:return None
        if not converged:return None
        # Same tolerance accounts for inversion conditioning, with a hard ceiling.
        toll=min(1e-10,max(128*eps,4*tolx*np.linalg.norm(np.linalg.inv(J),ord=np.inf)))
        if np.min(L)<-toll or np.max(L)>1+toll:return None
        return N,D,L,J

    def sample(self,positions):
        positions=np.atleast_2d(np.asarray(positions,dtype=np.float64))
        n=len(positions)
        result={k:np.full((n,)+shape,np.nan) for k,shape in [('u',(3,)),('G',(3,3)),('curl',(3,)),('strain',(3,3)),('pressure',())]}
        result.update(inside_lumen=np.zeros(n,dtype=bool),element_id=np.full(n,-1,dtype=np.int64))
        for i,x in enumerate(positions):
            if not np.isfinite(x).all():continue
            pad=128*np.finfo(float).eps*max(np.max(np.abs(x)),self.maxradius)
            candidates=self.tree.query_ball_point(x,self.maxradius+pad)
            for c in sorted(candidates):
                if np.any(x<self.low[c]-pad) or np.any(x>self.high[c]+pad):continue
                found=self.in_cell(x,c)
                if found is None:continue
                N,D,L,J=found;ids=self.cells[c]
                u=N@self.velocity[ids];G=self.velocity[ids].T@D@np.linalg.inv(J)
                for k,v in quantities(u,G).items():result[k][i]=v
                result['pressure'][i]=L@self.pressure[ids[:4]]
                result['inside_lumen'][i]=True;result['element_id'][i]=c
                break
        return result
