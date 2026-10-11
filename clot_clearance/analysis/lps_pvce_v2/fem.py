"""Conforming slit, genuine independent crack faces, mixed 3D linear elasticity."""
import sys
sys.path.insert(0,'/workspace/y_clot/envs/pvce_v2_packages')
from skfem import MeshTet,Basis,ElementVector,ElementTetP2,ElementTetP1,BilinearForm,asm
from skfem.helpers import sym_grad,ddot,div
from scipy.sparse import bmat
from scipy.sparse.linalg import spsolve
from analysis.lps_pvce_v2.common import *
import resource,traceback
L=1e-3
@BilinearForm
def shear(u,v,w):return 2*ddot(sym_grad(u),sym_grad(v))
@BilinearForm
def coupling(u,q,w):return div(u)*q
@BilinearForm
def mass(p,q,w):return p*q

def mesh(level,case):
 f=json.loads((BASE/'V2_DESIGN_FREEZE.json').read_text())['A']['mesh'];xx=np.array(f['coarse_x_um'],float);yy=np.array(f['y_um'],float);zz=np.array(f['coarse_z_um'],float)
 if level=='refined':
  mids=(xx[1:]+xx[:-1])/2;xx=np.sort(np.r_[xx,mids[(mids>-150)&(mids<50)]]);mids=(zz[1:]+zz[:-1])/2;zz=np.sort(np.r_[zz,mids[(mids>-100)&(mids<100)]])
 tip=[-100,-50,0][int(case[-1])]*1e-3
 orig=MeshTet.init_tensor(xx*1e-3,yy*1e-3,zz*1e-3);p=orig.p.copy();t=orig.t.copy();ids=np.flatnonzero((abs(p[2])<1e-12)&(p[0]<tip-1e-12));dup=np.arange(len(ids))+p.shape[1];mapping=np.arange(p.shape[1]);mapping[ids]=dup;upper=p[2,t].mean(0)>0;t[:,upper]=mapping[t[:,upper]];p=np.c_[p,p[:,ids]];m=MeshTet(p,t,sort_t=False)
 faces=m.boundary_facets();fc=m.p[:,m.facets[:,faces]].mean(1);cr=faces[(abs(fc[2])<1e-12)&(fc[0]<tip-1e-12)];coords=m.p[:,m.facets[:,cr]];areas=np.linalg.norm(np.cross((coords[:,1]-coords[:,0]).T,(coords[:,2]-coords[:,0]).T),axis=1)/2
 vol=abs(np.linalg.det(np.moveaxis(m.p[:,m.t[1:]]-m.p[:,m.t[0]][:,None,:],2,0))).sum()/6
 assert len(cr)>0 and abs(areas.sum()-2*(tip+.3)*.6)<1e-12
 assert abs(vol-.36)<1e-12
 # Each split vertex belongs only to its respective half; front remains shared.
 for k,d in zip(ids,dup):
  assert np.all(m.p[2,m.t[:,np.any(m.t==k,axis=0)]].mean(0)<0)
  assert np.all(m.p[2,m.t[:,np.any(m.t==d,axis=0)]].mean(0)>0)
 return m,ids,dup,cr,dict(vertices=m.nvertices,tetrahedra=m.nelements,split_vertices=len(ids),crack_face_triangles=len(cr),two_face_area_m2=float(areas.sum()*L**2),volume_m3=float(vol*L**3),tip_m=tip*L)

def run(level,case):
 out=BASE/'continuum'/level/case;out.mkdir(parents=True,exist_ok=True);assert not (out/'ATTEMPT.json').exists();save(out/'ATTEMPT.json',dict(utc=time.strftime('%FT%TZ',time.gmtime()),independent_equilibrium_attempt=1));start=time.perf_counter();cpu=time.process_time();cfg=json.loads((BASE/'V2_DESIGN_FREEZE.json').read_text())['A'];m,ids,dup,cr,ma=mesh(level,case)
 bu=Basis(m,ElementVector(ElementTetP2()),intorder=4);bp=Basis(m,ElementTetP1(),intorder=4);A=asm(shear,bu).tocsr();B=asm(coupling,bu,bp).tocsr();M=asm(mass,bp).tocsr();lam=K/G-2/3;C=M/lam;H=bmat([[A,B.T],[B,-C]],format='csr')
 # Vector P2 dofs are interleaved by component; locate ALL nodal+edge control DOFs.
 idx=np.arange(bu.N);zdof=idx[idx%3==2];lo=zdof[np.isclose(bu.doflocs[2,zdof],-.475,atol=1e-12,rtol=0)];hi=zdof[np.isclose(bu.doflocs[2,zdof],.475,atol=1e-12,rtol=0)]
 na=np.flatnonzero(np.all(np.isclose(m.p.T,[-.3,-.3,-.475],atol=1e-12,rtol=0),axis=1))[0];nb=np.flatnonzero(np.all(np.isclose(m.p.T,[.3,-.3,-.475],atol=1e-12,rtol=0),axis=1))[0];gauge=np.array([bu.nodal_dofs[0,na],bu.nodal_dofs[1,na],bu.nodal_dofs[1,nb]]);D=np.unique(np.r_[lo,hi,gauge]);free=np.setdiff1d(np.arange(H.shape[0]),D);q=np.zeros(H.shape[0]);q[lo]=-.001;q[hi]=.001
 assembly_seconds=time.perf_counter()-start;ts=time.perf_counter();q[free]=spsolve(H[free][:,free].tocsc(),-H[free][:,D]@q[D],permc_spec='COLAMD');solve_seconds=time.perf_counter()-ts
 assert np.all(np.isfinite(q));u=q[:bu.N];p=q[bu.N:];res=H@q;Fscale=G*L**2;Uscale=G*L**3;U=.5*(u@(A@u)+p@(C@p))*Uscale;top=res[hi].sum()*Fscale;bottom=res[lo].sum()*Fscale;Q=(top-bottom)/2;ud=bu.nodal_dofs;gap=(u[ud[2,dup]]-u[ud[2,ids]])*L
 ctrl=np.r_[lo,hi];ufree=np.setdiff1d(np.arange(bu.N),D);frel=np.linalg.norm(res[ufree])/np.linalg.norm(res[ctrl]);mixed=np.linalg.norm(res[free])/np.linalg.norm(res[D]);work=q[D]@res[D]*Uscale;volerr=abs(ma['volume_m3']/3.6e-10-1)
 # Physical unprojected energy is an additional locking/projection diagnostic, not the mixed discrete variational energy used for reaction consistency.
 field=bu.interpolate(u);dv=div(field);pv=bp.interpolate(p);vol_mismatch=float(np.sqrt(np.sum((dv-pv/lam)**2*bu.dx)/np.sum(dv**2*bu.dx)));U_primal=U+.5*lam*np.sum((dv-pv/lam)**2*bu.dx)*Uscale
 c=cfg['criteria'];gates=dict(free_force=frel<c['relative_free_force'],mixed_residual=mixed<c['relative_mixed_residual'],work_energy=abs(work/(2*U)-1)<c['relative_work_energy'],reaction_balance=abs(top+bottom)/abs(Q)<c['relative_reaction_balance'],minimal_gauge_force=np.linalg.norm(res[gauge])*Fscale/abs(Q)<c['gauge_reaction_over_Q'],positive_crack_opening=gap.min()>=c['minimum_crack_gap_m'],physical_volume=volerr<c['relative_volume_error'],independent_faces=len(ids)>0,pressure_coercivity=lam>0)
 result=dict(level=level,case=case,accepted=all(gates.values()),gates=gates,U_J=float(U),top_reaction_N=float(top),bottom_reaction_N=float(bottom),relative_free_force=float(frel),max_free_force_N=float(abs(res[ufree]).max()*Fscale),relative_mixed_residual=float(mixed),work_energy_relative=float(abs(work/(2*U)-1)),gauge_force_over_reaction=float(np.linalg.norm(res[gauge])*Fscale/abs(Q)),crack_vertex_gap_range_m=[float(gap.min()),float(gap.max())],mesh=ma,displacement_dofs=bu.N,pressure_dofs=bp.N,axial_control_dofs_each=len(lo),transverse_gauge_dofs=gauge.tolist(),pressure_projection_relative_L2=vol_mismatch,primal_energy_J=float(U_primal),primal_vs_mixed_energy_relative=float(U_primal/U-1),energy_scope='Mixed variational stored energy .5uAu+.5pCp, consistent with reaction; unprojected lambda*div(u)^2 energy also reported to expose approximation.',crack_traction='Zero applied natural traction on distinct boundary facets. Discrete weak equilibrium verified; pointwise stress traction at singular front is not claimed zero.',wall_seconds=time.perf_counter()-start,assembly_seconds=assembly_seconds,sparse_solve_seconds=solve_seconds,cpu_seconds=time.process_time()-cpu,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,threads=1,GPU=False)
 np.savez_compressed(out/'equilibrium.npz',vertices_m=m.p.T*L,tetrahedra=m.t.T,facets=m.facets.T,crack_facets=cr,split_original=ids,split_duplicate=dup,displacement_dof_coordinates_m=bu.doflocs.T*L,displacement_dofs_m=u*L,vertex_displacement_m=u[ud].T*L,pressure_Pa=p*G,nodal_dofs=ud,low_dofs=lo,high_dofs=hi,gauge_dofs=gauge,residual_N=res[:bu.N]*Fscale,crack_vertex_gap_m=gap)
 save(out/'RESULT.json',result);print((out/'RESULT.json').read_text(),flush=True);return result
if __name__=='__main__':
 try:run(sys.argv[1],sys.argv[2])
 except BaseException:
  p=BASE/'continuum'/sys.argv[1]/sys.argv[2];p.mkdir(parents=True,exist_ok=True);(p/'ERROR.txt').write_text(traceback.format_exc());raise
