"""Explain the existing P1 no-slip wall stencil on actual solved fields.

This is a discrete identity check, not an independent physical validation.
"""
from pathlib import Path
import argparse, csv
import numpy as np
from case_common import *
from flow_solver_support.wss_case import material
from analyze_vessel import csvout

p = argparse.ArgumentParser()
p.add_argument('--case', type=Path, required=True)
a = p.parse_args()
case = a.case.resolve()
mesh = np.load(case/'SV_MESH/mesh_arrays.npz')
flow = np.load(case/'frozen_flow/flow_arrays_si.npz')
x,t,b,tags = [mesh[k] for k in ['points_m','tetra','boundary_triangles','facet_tags']]
u = flow['velocity_m_s']
wallids = np.flatnonzero(tags == 1)
wall = b[wallids]
own = wss.boundary_owners(t,wall)
c,area,n = wss.wall_geometry(x,t,wall,own)
assert np.max(np.abs(u[wall])) == 0
is_wall_vertex = np.any(t[own][:,:,None] == wall[:,None,:], axis=2)
assert np.all(is_wall_vertex.sum(axis=1) == 3)
fourth = t[own][~is_wall_vertex]
h = -np.einsum('ij,ij->i',x[fourth]-c,n)
assert np.all(h > 0)
ut = u[fourth] - np.einsum('ij,ij->i',u[fourth],n)[:,None]*n
mu = material(case)['mu_Pa_s']
simple_vector = -mu*ut/h[:,None]
production_vector = wss.tangential_traction(wss.p1_gradients(x,t[own],u),n,mu)
delta = np.abs(simple_vector-production_vector)
dump(case/'reports/wall_stencil_identity.json',dict(
    case=case.name,formula='tau=-mu*u_fourth_tangent/h_normal for three exact no-slip face vertices',
    wall_faces=len(wall),max_vector_component_difference_Pa=float(delta.max()),
    max_magnitude_difference_Pa=float(np.max(np.abs(np.linalg.norm(simple_vector,axis=1)-np.linalg.norm(production_vector,axis=1)))),
    meaning='Identity of the present discrete stencil only; not a physical accuracy test'))
lookup={int(f):i for i,f in enumerate(wallids)}
rows=[]
for row in csv.DictReader((case/'reports/known_jump_pairs.csv').open()):
    result=dict(row)
    for side in ['a','b']:
        i=lookup[int(row['facet_'+side])]
        result['normal_height_'+side+'_um']=float(h[i]*1e6)
        result['fourth_node_'+side+'_zero_based']=int(fourth[i])
        result['fourth_tangent_speed_'+side+'_mm_s']=float(np.linalg.norm(ut[i])*1e3)
        for j,axis in enumerate('xyz'):
            result['fourth_'+side+'_'+axis+'_um']=float(x[fourth[i],j]*1e6)
    rows.append(result)
csvout(case/'reports/known_jump_stencils.csv',rows)
allrows=[]
for p in sorted(V.glob('stage[34]/*/reports/known_jump_stencils.csv')):
    allrows.extend(list(csv.DictReader(p.open())))
csvout(V/'data/known_jump_stencils.csv',allrows)
print(case.name,'max discrete vector identity difference',delta.max(),'Pa')
