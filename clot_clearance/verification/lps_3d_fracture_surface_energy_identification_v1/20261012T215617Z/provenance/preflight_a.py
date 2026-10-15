from analysis.dual_track_v1.core_a import *
from analysis.dual_track_v1.native_a import Native
import time,resource,subprocess,traceback

def main():
 start=time.monotonic();src=ROOT/'analysis/dual_track_v1/energy_y.cpp';lib=src.with_suffix('.so');cmd=['g++','-O3','-shared','-fPIC',str(src),'-o',str(lib)];r=subprocess.run(cmd,capture_output=True,text=True);save(A/'EVALUATOR_BUILD.json',dict(command=cmd,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr,source_sha=sha(src)));assert r.returncode==0
 g,e,a0,mask,low,high,xx,cross=model();X=g['X'];b=e.b;n=len(X)
 # Independent integer-lattice enumeration with only Y wrapped.
 ijk=np.stack(np.meshgrid(np.arange(100),np.arange(8),np.arange(20),indexing='ij'),axis=-1).reshape(-1,3);keys=[]
 for di in range(-3,4):
  for dj in range(-3,4):
   for dk in range(-3,4):
    if not 0<di*di+dj*dj+dk*dk<=(DELTA/DX)**2:continue
    target=ijk+[di,dj,dk];valid=(target[:,0]>=0)&(target[:,0]<100)&(target[:,2]>=0)&(target[:,2]<20);target[:,1]%=8;j=(target[:,0]*8+target[:,1])*20+target[:,2];i=np.arange(n);valid &=i<j;keys.append(i[valid]*n+j[valid])
 lattice=np.unique(np.concatenate(keys));kd=b[:,0].astype(np.int64)*n+b[:,1];assert np.array_equal(kd,lattice)
 image=-np.rint((X[b[:,1],1]-X[b[:,0],1])/PERIOD).astype(np.int8);assert np.all(b[:,0]<b[:,1]);assert e.r.max()<=DELTA*(1+1e-12);assert np.all(e.m>0)
 tips=[-.001,-.00095,-.00090];alive=np.array([~(cross&(xx<t-1e-12)) for t in tips],np.uint8);removed=(alive==0).sum(1);assert np.all(alive[1:]<=alive[:-1]);assert removed[1]-removed[0]==removed[2]-removed[1]>0
 # Translation by one Y lattice layer leaves each virtual cut and m invariant.
 perm=((ijk[:,0]*8+(ijk[:,1]+1)%8)*20+ijk[:,2]);mapped=np.sort(perm[b],axis=1);mk=mapped[:,0]*n+mapped[:,1];idx=np.searchsorted(kd,mk);assert np.array_equal(kd[idx],mk);assert np.all(alive[:,idx]==alive);assert np.allclose(e.m[perm],e.m,rtol=1e-14,atol=0)
 np.savez_compressed(A/'reference_topology.npz',X=X,bonds=b,bond_reference_id=np.arange(len(b))+1,periodic_image_y=image,r0_m=e.r,h=e.h,m=e.m,alive=alive,tip_x_m=tips,mask=mask,low=low,high=high,reference_intersection_x_m=xx,cross_z0=cross)
 report=dict(status='PENDING_NATIVE',nodes=n,bonds=len(b),cross_periodic_bonds=int(np.count_nonzero(image)),total_volume_m3=n*DX**3,total_mass_kg=n*DX**3*RHO,weighted_m_range_m4=[float(e.m.min()),float(e.m.max())],cut_counts=removed.tolist(),cut_increments=np.diff(removed).tolist(),controlled_nodes=int((low|high).sum()),controlled_DOFs=int(mask.sum()),free_DOFs=int((~mask).sum()),integer_stencil_matches_KD=True,Y_shift_invariant_cuts_and_m=True,minimum_image_rmax_m=float(e.r.max()),Y_period_gt_twice_delta=PERIOD>2*DELTA,evaluator_sha=sha(lib),checks=[])
 save(A/'PERIODIC_TOPOLOGY_AUDIT.json',report)
 for case,al in [('full',np.ones(len(b),np.uint8)),('A0',alive[0])]:
  native=Native(A/'preflight'/case,g,al)
  for state in ['reference','periodic_nonaffine']:
   x=X.copy()
   if state!='reference':
    x[:,2]*=1+2e-6/.00095;x[:,0]+=.35e-6*np.sin(2*np.pi*X[:,1]/PERIOD)*np.cos(np.pi*X[:,2]/.001);x[:,1]+=.30e-6*np.sin(2*np.pi*X[:,1]/PERIOD)
   fn,diag,tn=native.evaluate(x);U,f=e.eval(x,al,True);Uc,fc,cp=contact(x,g,al);th=e.theta(x,al);norm=np.linalg.norm(fn);error=float(np.linalg.norm(fn-f-fc)/norm) if norm else float(np.linalg.norm(f+fc));terr=float(np.linalg.norm(tn-th)/np.linalg.norm(th)) if np.linalg.norm(th) else float(np.linalg.norm(tn));rr=dict(case=case,state=state,U_J=U,force_relative=error,theta_relative=terr,native=diag,contact_independent_J=Uc,contact_pairs=len(cp),contact_force_N=float(np.linalg.norm(fc)),max_free_N=float(abs(fn[~mask]).max()))
   report['checks'].append(rr);save(A/'PERIODIC_TOPOLOGY_AUDIT.json',report);np.savez_compressed(A/'preflight'/case/(state+'.npz'),x=x,native_force=fn,independent_force=f,native_theta=tn,independent_theta=th,contact_force=fc)
   assert error<1e-9 and terr<1e-9,(case,state,error,terr);assert Uc<1e-24 and diag['contact_J']<1e-24,'Unexpected periodic contact'
  native.close()
 # Independent finite differences in a periodic smooth displacement, fixed active A0.
 rng=np.random.default_rng(917);checks=[];U,f=e.eval(x,alive[0],True)
 for _ in range(4):
  d=rng.normal(size=x.shape);d[mask]=0;d/=np.linalg.norm(d);h=1e-9;fd=(e.eval(np.ascontiguousarray(x+h*d),alive[0])-e.eval(np.ascontiguousarray(x-h*d),alive[0]))/(2*h);analytic=-float(np.sum(f*d));err=abs(fd-analytic)/np.linalg.norm(f);checks.append(dict(FD_N=fd,analytic_N=analytic,scaled_error=err))
 assert max(v['scaled_error'] for v in checks)<1e-6
 report.update(status='PASS',gradient_checks=checks,wall_s=time.monotonic()-start,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,native_scope='Single MPI rank, public atom virtual-cut arrays, full initial native m, static force calls without integrator, unchanged production lib; explicit force/theta comparison for full and A0 across seam')
 save(A/'PERIODIC_TOPOLOGY_AUDIT.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
 try:main()
 except BaseException as e:
  p=A/'PERIODIC_TOPOLOGY_AUDIT.json';d=json.loads(p.read_text()) if p.exists() else {};d.update(status='PERIODIC_TOPOLOGY_NOT_VALIDATED',error=repr(e),traceback=traceback.format_exc());save(p,d);raise
