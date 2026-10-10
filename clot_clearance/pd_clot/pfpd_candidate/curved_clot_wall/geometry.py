"""Physical P1 tetrahedron volumes and outward boundary triangles, in SI."""
import numpy as np
from scipy.spatial import cKDTree


def tetra_geometry(X, cells):
    X=np.asarray(X,float);cells=np.asarray(cells,np.int64).copy()
    edge=X[cells[:,1:]]-X[cells[:,0,None]]
    det=np.linalg.det(edge)
    bad=det<0
    cells[bad,1],cells[bad,2]=cells[bad,2].copy(),cells[bad,1].copy()
    vol=np.abs(det)/6
    if np.any(vol<=0):raise ValueError('Zero-volume physical tetrahedron')
    lump=np.bincount(cells.ravel(),weights=np.repeat(vol/4,4),minlength=len(X))
    return cells,vol,lump,int(bad.sum())


def surface_geometry(X, triangles, wall):
    triangles=np.asarray(triangles,np.int64).copy()
    a,b,c=X[triangles[:,0]],X[triangles[:,1]],X[triangles[:,2]]
    vec=np.cross(b-a,c-a);center=(a+b+c)/3
    # A point strictly inside the convex cloud also handles faces through origin.
    interior=np.mean(X,axis=0)
    reverse=np.where(wall,vec[:,2]>0,np.einsum('ij,ij->i',vec,center-interior)<0)
    triangles[reverse,1],triangles[reverse,2]=triangles[reverse,2].copy(),triangles[reverse,1].copy()
    vec[reverse]*=-1;twice=np.linalg.norm(vec,axis=1)
    if np.any(twice<=0):raise ValueError('Zero-area physical face')
    return triangles,twice/2,vec/twice[:,None],center


def ring_distance(X,R):
    return np.hypot(np.linalg.norm(np.asarray(X)[...,:2],axis=-1)-R,np.asarray(X)[...,2])


def cloud_topology(X,delta):
    pairs=cKDTree(X).query_pairs(delta,output_type='ndarray')
    pairs=np.sort(pairs,axis=1);pairs=pairs[np.lexsort((pairs[:,1],pairs[:,0]))]
    # Exactly-on-horizon zero weights do not supply rank/support.
    pairs=pairs[np.linalg.norm(X[pairs[:,1]]-X[pairs[:,0]],axis=1)<delta]
    return pairs


def nearest_statistics(X):
    d=cKDTree(X).query(X,k=2)[0][:,1]
    return d,dict(zip(['min','P05','median','P95','max'],np.quantile(d,[0,.05,.5,.95,1]).tolist()))


def prescribed_wall_positions(X, fixed, displacement):
    """Accept a constant SI translation or an explicit spatial wall displacement."""
    u=displacement(X[fixed]) if callable(displacement) else np.asarray(displacement,float)
    result=np.array(X,copy=True);result[fixed]+=u
    if not np.isfinite(result).all():raise ValueError('Invalid prescribed wall displacement')
    return result
