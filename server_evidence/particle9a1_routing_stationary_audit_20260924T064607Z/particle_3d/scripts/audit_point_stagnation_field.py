#!/usr/bin/env python3
"""Analytic audit of saved point paths in rank-one, no-slip P1 wall tetrahedra.
No integration and no field reconstruction: exactly the original nodal P1 field.
"""
import numpy as np
from particle_3d.routing_stationary_audit import *
from particle_3d.particle8_replay import REPO

def main():
 R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';D=R/'data';env=environment();g=env.field.geometry
 wall_faces={tuple(sorted(row)):i for i,row in enumerate(env.wall.global_node_ids)};records=[]
 for r in read(D/'point_resolution.json')['results']:
  if r['baseline_outlet']!='NO_EXIT':continue
  pid=r['particle_id'];p=np.load(R/'diagnostic_outputs/point_resolution'/f'point_{pid:06d}_level2.npz')['path'];sample=env.field.sample(p[-1,1:]);tet=sample.tetra_id;nodes=g.tetra[tet];v=env.field.velocity_nodes_m_s[nodes];G=sample.velocity_gradient_s_inv
  zero=np.all(v==0,axis=1);nz=np.flatnonzero(~zero);sigma=float(np.trace(G));tol=float(g.weight_tolerance[tet]);w123=(p[:,1:]-g.origins[tet])@g.inverse[tet].T;weights=np.column_stack([1-w123.sum(1),w123]);indices=np.flatnonzero(np.min(weights,axis=1)>=tol*4)
  proof=False;details={}
  if zero.sum()==3 and len(indices) and sigma<0:
   j=int(nz[0]);candidates=indices[weights[indices,j]>max(100*tol,1e-8)];k=int(candidates[0]) if len(candidates) else int(indices[0]);x=p[k,1:];u=weights[k]@v;xinf=x-u/sigma
   b123=g.inverse[tet]@(xinf-g.origins[tet]);binf=np.r_[1-b123.sum(),b123];wall_key=tuple(sorted(nodes[zero]));gradb=-g.inverse[tet].sum(0) if j==0 else g.inverse[tet][j-1]
   expected=np.outer(v[j],gradb);relative_error=float(np.linalg.norm(G-expected)/max(np.linalg.norm(G),np.finfo(float).tiny));roundoff=256*np.finfo(float).eps*max(1,np.linalg.cond(g.inverse[tet]))
   proof=bool(wall_key in wall_faces and np.min(binf)>=-4*tol and abs(binf[j])<=4*tol and relative_error<=roundoff and abs(sigma)>roundoff*np.linalg.norm(G))
   details=dict(interior_node_local_index=j,wall_triangle_id=wall_faces.get(wall_key),witness_sample_index=k,witness_time_s=float(p[k,0]),witness_position_m=x.tolist(),witness_barycentric=weights[k].tolist(),limit_position_m=xinf.tolist(),limit_barycentric=binf.tolist(),gradient_rank_one_relative_error=relative_error,rank_one_verification_budget=roundoff,barycentric_tolerance=tol,
    analytic_ode='u(x)=b_j(x)*v_j; db_j/dt=sigma*b_j; x(t)=x0+u0*(exp(sigma*t)-1)/sigma',
    explanation='sigma<0, three zero-velocity WALL nodes, and limit barycentric coordinates on the same WALL face prove asymptotic approach without a finite-time cap crossing while the exact P1 path stays in this convex tetrahedron.')
  records.append(dict(particle_id=pid,accepted=r['accepted'],tetra_id=tet,tetra_global_node_ids_zero_based=nodes.tolist(),nodal_velocities_m_s=v.tolist(),velocity_gradient_s_inv=G.tolist(),divergence_s_inv=sigma,zero_velocity_nodes=int(zero.sum()),analytic_wall_stagnation_confirmed=proof,**details));print(pid,'P1 wall stagnation proof',proof,'divergence',sigma)
 dump(D/'point_stagnation_field.json',dict(role='FROZEN_P1_FIELD_LOCAL_NON_SOLENOIDAL_STAGNATION_AUDIT_NO_PHYSICS_CHANGE',records=records,all_confirmed=all(r['analytic_wall_stagnation_confirmed'] for r in records),count=len(records),not_a_physiological_trapping_claim=True))
 csvwrite(D/'point_stagnation_field.csv',[dict(particle_id=r['particle_id'],accepted=r['accepted'],tetra_id=r['tetra_id'],divergence_s_inv=r['divergence_s_inv'],zero_velocity_nodes=r['zero_velocity_nodes'],analytic_wall_stagnation_confirmed=r['analytic_wall_stagnation_confirmed'],wall_triangle_id=r.get('wall_triangle_id')) for r in records])
if __name__=='__main__':main()
