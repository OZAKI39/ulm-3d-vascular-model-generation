"""Read-only algebra/geometry checks, not a CFD or material simulation."""
from analysis.dual_track_v1.common import *
import sys
level=sys.argv[1];out=B/level;a=np.load(out/'SV_MESH/mesh_arrays.npz');t=np.load(out/'CLOT_SURFACE_TRACTIONS.npz');P=a['points_m'];T=a['tetra'];ids=t['adjacent_mesh_tetra_ids']-1;xyz=P[T[ids]];edges=xyz[:,1:]-xyz[:,:1]
A0=np.array([[.3,.5,1.2],[-.2,.7,-.6],[.9,-.4,-.8]])
velocity=np.einsum('ij,ekj->eki',A0,xyz)+np.array([.002,-.001,.0003]);gradient=np.linalg.solve(edges,velocity[:,1:]-velocity[:,:1]).transpose(0,2,1);error=float(np.max(np.linalg.norm(gradient-A0,axis=(1,2)))/np.linalg.norm(A0));assert error<1e-7
# Exact linear-pressure moment via barycentric monomial integrals, independent of quadrature.
r=t['triangles_xyz_m']-t['reference_point_m'];n=t['n_clot'];area=t['triangle_area_m2'];q=t['quadrature_coordinates_m'];qp=t['quadrature_pressure_Pa'];bary=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]]);pv=np.linalg.solve(bary,qp.T).T;analytic=np.zeros(3)
for i in range(3):
 for j in range(3):analytic-=np.einsum('f,fj->j',area*pv[:,j]*(1/6 if i==j else 1/12),np.cross(r[:,i],n))
analytic+=np.einsum('f,fj->j',area,np.cross(r.mean(1),t['viscous_traction_Pa']));moment_error=float(np.linalg.norm(analytic-t['integrated_total_torque_Nm'])/np.linalg.norm(analytic));assert moment_error<1e-10
normal_error=float(abs(t['n_clot']+t['n_fluid']).max());tangent_error=float(abs(np.einsum('ij,ij->i',t['tangential_traction_Pa'],t['n_clot'])).max());assert normal_error==0 and tangent_error<1e-12
save(out/'TRACTION_MATH_AUDIT.json',dict(status='PASS',scope='Read-only manufactured affine gradient on actual adjacent elements + independent exact linear-pressure moment integral; no new flow/mechanical solve',affine_gradient_relative_max=error,gradient_gate=1e-7,linear_pressure_moment_relative_error=moment_error,moment_gate=1e-10,normal_sign_error=normal_error,tangent_dot_normal_max_Pa=tangent_error,solver_fields_unchanged=True))
