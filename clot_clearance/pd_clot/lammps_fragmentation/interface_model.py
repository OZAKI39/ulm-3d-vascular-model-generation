"""Cubic clot + fixed one-layer anchors; bulk/interface/wall real bond types."""
import numpy as np
from scipy.spatial import cKDTree
from .geometry_adapter import prepare_geometry

BULK,INTERFACE,WALL=1,2,3

def attach_wall(g,r):
    X=g['X'];n=len(X);base=np.flatnonzero(g['base'])
    layers=int(r.get('anchor_layers',1))
    if layers<1:raise ValueError('At least one fixed anchor layer is required')
    anchor=np.vstack([X[base]-[0,0,k*r['dx_m']] for k in range(1,layers+1)])
    links=[];tree=cKDTree(anchor)
    for i in base:
        links.extend((int(i),n+j) for j in tree.query_ball_point(X[i],r['delta_m']*(1+1e-12)))
    wall=cKDTree(anchor).query_pairs(r['delta_m']*(1+1e-12),output_type='ndarray')+n
    bonds=np.vstack([g['bonds'],np.asarray(links,int),wall]);types=np.r_[np.full(len(g['bonds']),BULK),np.full(len(links),INTERFACE),np.full(len(wall),WALL)]
    order=np.lexsort((bonds[:,1],bonds[:,0]));bonds=bonds[order];types=types[order];N=n+len(anchor)
    out=dict(g);out.update(X=np.vstack([X,anchor]),bonds=bonds,bond_types=types,clot_count=n,
        anchors=np.arange(N)>=n,volume=np.full(N,r['dx_m']**3),mass=np.full(N,r['rho_kg_m3']*r['dx_m']**3),
        surface=np.r_[g['surface'],np.zeros(len(anchor),bool)],base=np.r_[g['base'],np.zeros(len(anchor),bool)],
        areas=np.r_[g['areas'],np.zeros(len(anchor))],normals=np.vstack([g['normals'],np.zeros_like(anchor)]))
    out['info']={**g['info'],'node_count':N,'clot_nodes':n,'anchor_nodes':len(anchor),'bond_count':len(bonds),
        'bulk_bonds':int((types==BULK).sum()),'interface_bonds':int((types==INTERFACE).sum()),'wall_bonds':int((types==WALL).sum()),
        'average_bonds_per_node':2*len(bonds)/N,'wall_interface_model':'BREAKABLE_INTERFACE_LIU2024_V1',
        'permanently_fixed_clot_nodes':0,'anchors_in_fragment_statistics':False,'anchor_layers':layers}
    assert np.all(bonds[types==INTERFACE,0]<n) and np.all(bonds[types==INTERFACE,1]>=n)
    assert np.all(g['base'][bonds[types==INTERFACE,0]])
    return out

def prepare_interface_geometry(r):return attach_wall(prepare_geometry(r),r)

def calibration_geometry(r,kind):
    dx=r['dx_m'];nz=8 if kind=='bulk' else 4
    axes=[(np.arange(8)-3.5)*dx,(np.arange(8)-3.5)*dx,(np.arange(nz)+.5)*dx]
    X=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3);N=len(X)
    b=cKDTree(X).query_pairs(r['delta_m']*(1+1e-12),output_type='ndarray');b=b[np.lexsort((b[:,1],b[:,0]))]
    g=dict(X=X,bonds=b,bond_types=np.ones(len(b),int),clot_count=N,anchors=np.zeros(N,bool),
           volume=np.full(N,dx**3),mass=np.full(N,r['rho_kg_m3']*dx**3),base=X[:,2]<dx,
           surface=np.zeros(N,bool),areas=np.zeros(N),normals=np.zeros_like(X),
           info=dict(node_count=N,bond_count=len(b),average_bonds_per_node=2*len(b)/N,validation='CALIBRATION_CELL'))
    return attach_wall(g,r) if kind=='interface' else g

def write_typed_data(p,g,r):
    X=g['X'];pad=4*r['delta_m'];lo=X.min(0)-pad;hi=X.max(0)+pad;mass=r['rho_kg_m3']*r['dx_m']**3
    with p.open('w') as f:
        f.write(f'Liu strength-mapped LPS; anchors excluded from clot statistics\n\n{len(X)} atoms\n{len(g["bonds"])} bonds\n\n2 atom types\n3 bond types\n\n')
        for j,a in enumerate('xyz'):f.write(f'{lo[j]:.17g} {hi[j]:.17g} {a}lo {a}hi\n')
        f.write(f'\nMasses\n\n1 {mass:.17g}\n2 {mass:.17g}\n\nAtoms # bond\n\n')
        for i,x in enumerate(X):f.write(f'{i+1} 1 {2 if g["anchors"][i] else 1} {x[0]:.17g} {x[1]:.17g} {x[2]:.17g}\n')
        f.write('\nBonds\n\n')
        for i,((a,b),t) in enumerate(zip(g['bonds'],g['bond_types']),1):f.write(f'{i} {t} {a+1} {b+1}\n')
