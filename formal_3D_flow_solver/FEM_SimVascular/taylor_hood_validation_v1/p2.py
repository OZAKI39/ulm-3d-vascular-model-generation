"""Straight-sided native TET10/TRI6 basis and exact polynomial quadrature.

Ordering is derived from pinned nn.cpp solver_to_basis_node_map and
FE/Basis/NodeOrderingConventions.cpp, not inferred from an assumed VTK order.
"""
import numpy as np
from itertools import permutations

EDGES=np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])
TRI_EDGES=np.array([[0,1],[1,2],[2,0]])
TRI_Q=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])
TRI_W=np.full(3,1/3)

def basis(lam):
    lam=np.asarray(lam)
    return np.concatenate([lam*(2*lam-1),4*lam[...,EDGES[:,0]]*lam[...,EDGES[:,1]]],axis=-1)

def tri_basis(lam):
    lam=np.asarray(lam)
    return np.concatenate([lam*(2*lam-1),4*lam[...,TRI_EDGES[:,0]]*lam[...,TRI_EDGES[:,1]]],axis=-1)

def reference_gradient(lam):
    # Solver reference corners are (1,0,0),(0,1,0),(0,0,1),(0,0,0).
    gl=np.array([[1.,0,0],[0,1.,0],[0,0,1.],[-1.,-1.,-1.]])
    lam=np.asarray(lam)
    corners=(4*lam-1)[...,None]*gl
    edges=4*(lam[...,EDGES[:,0],None]*gl[EDGES[:,1]]+lam[...,EDGES[:,1],None]*gl[EDGES[:,0]])
    return np.concatenate([corners,edges],axis=-2)

def corner_gradients(x):
    # x has (...,4,3); returns physical gradients of the four barycentrics.
    inv=np.linalg.inv(np.swapaxes(x[...,1:,:]-x[...,:1,:],-1,-2))
    return np.concatenate([-inv.sum(axis=-2,keepdims=True),inv],axis=-2)

def gradient_at(lam,gl):
    lam=np.asarray(lam)
    return np.concatenate([(4*lam-1)[...,None]*gl,
        4*(lam[...,EDGES[:,0],None]*gl[...,EDGES[:,1],:]+lam[...,EDGES[:,1],None]*gl[...,EDGES[:,0],:])],axis=-2)

def tetra_rule15():
    # Pinned nn_elem_gip.h:493-543, normalized to sum=1 here.
    pts=[[.25]*4];weights=[.0302836780970890]
    for a,b,w in [(0.,1/3,.0060267857142860),(.727272727272727,.090909090909091,.0116452490860290)]:
        for i in range(4):
            p=[b]*4;p[i]=a;pts.append(p);weights.append(w)
    for p in sorted(set(permutations([.066550153573664]*2+[.433449846426336]*2))):
        pts.append(p);weights.append(.0109491415613860)
    return np.array(pts),np.array(weights)*6

def triangle_rule7():
    # Exact literals/order from pinned nn_elem_gip.h:740-766, area-normalized.
    xy=np.array([[.333333333333333,.333333333333333],
        [.797426985353087,.101286507323456],[.101286507323456,.797426985353087],
        [.101286507323456,.101286507323456],
        [.059715871789770,.470142064105115],[.470142064105115,.059715871789770],
        [.470142064105115,.470142064105115]])
    return np.column_stack([xy,1-xy.sum(axis=1)]),np.array([.225]+[.125939180544827]*3+[.132394152788506]*3)

def unique_elevation(points,tetra):
    edge_pairs=np.sort(tetra[:,EDGES].reshape(-1,2),axis=1)
    unique,inverse=np.unique(edge_pairs,axis=0,return_inverse=True)
    newpoints=np.vstack([points,points[unique].mean(axis=1)])
    return newpoints,np.column_stack([tetra,inverse.reshape(-1,6)+len(points)]),unique

def face_elevation(tri,edges,ncorner):
    # Integer packed keys provide deterministic exact global edge lookup.
    keys=edges[:,0].astype(np.int64)*ncorner+edges[:,1]
    e=np.sort(tri[:,TRI_EDGES],axis=-1)
    requested=e[...,0].astype(np.int64)*ncorner+e[...,1]
    idx=np.searchsorted(keys,requested)
    assert np.array_equal(keys[idx],requested)
    return np.column_stack([tri,idx+ncorner])

def barycentric(x,corner):
    tail=np.linalg.solve(np.swapaxes(corner[...,1:,:]-corner[...,:1,:],-1,-2),(x-corner[...,0,:])[...,None])[...,0]
    return np.concatenate([1-tail.sum(axis=-1,keepdims=True),tail],axis=-1)
