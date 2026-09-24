from itertools import permutations, product
from pathlib import Path
import os
import numpy as np
import pyvista as pv
import pytest
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.flowfield_conservation_diagnosis import face_topology, p1_divergence, slab_balance


@pytest.fixture(scope='session')
def cube():
    unit=np.array(list(product([0.,1.],repeat=3)))
    ids={tuple(x):i for i,x in enumerate(unit)};tetra=[]
    for order in permutations(range(3)):
        x=np.zeros(3);row=[ids[tuple(x)]]
        for j in order:
            x=x.copy();x[j]=1;row.append(ids[tuple(x)])
        tetra.append(row)
    t=np.array(tetra);scale=2e-6;origin=np.array([100e-6,57e-6,150e-6]);p=origin+scale*unit
    top=face_topology(t);bf=top['boundary_faces'];owner=top['boundary_owners']
    left=np.all(unit[bf,0]==0,axis=1);right=np.all(unit[bf,0]==1,axis=1)
    boundaries={'INLET':(bf[left],owner[left]),'OUTLET':(bf[right],owner[right]),'WALL':(bf[~(left|right)],owner[~(left|right)])}
    return dict(points=p,tetra=t,unit=unit,scale=scale,origin=origin,topology=top,boundaries=boundaries,
                center=origin+scale*np.array([.37,.5,.5]),normal=np.array([1.,0,0]))


def balance(cube,u):
    p,t=cube['points'],cube['tetra'];d,v=p1_divergence(p,t,u)
    return slab_balance(p,t,u,d,v,cube['topology'],cube['boundaries'],cube['center'],cube['normal'])[0]


def grid(cube,u):
    g=pv.UnstructuredGrid(np.c_[np.full(len(cube['tetra']),4),cube['tetra']].ravel(),
                          np.full(len(cube['tetra']),10,np.uint8),cube['points'])
    g['Velocity']=u;return g


@pytest.fixture(scope='session')
def real():
    repo=Path(__file__).resolve().parents[3]
    summary,mesh,flow,boundaries=read_frozen(repo/'formal_3D_flow_solver/FEM_SimVascular')
    field=FrozenFEMField.from_grids(mesh,flow)
    default='/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps'
    case=Path(os.environ.get('P9A3B_CASE',default))
    return repo,case,mesh,flow,boundaries,field
