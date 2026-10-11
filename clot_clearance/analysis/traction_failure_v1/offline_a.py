"""Geometry weights and one small-specimen energy diagnostic; never native failure."""
from analysis.traction_failure_v1.common import *
import sys,ctypes
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import bicgstab
from scipy.spatial import cKDTree

class Potential:
 """Locked uniform LPS expanded potential, optionally symmetric bond weights.

 U=V/2*(K-5G/3)*sum(theta^2)+15G V^2/2*sum_b h*w*e^2/r*(1/m_i+1/m_j).
 theta_i=3 V/m_i sum_alive h*w*e. m is the ORIGINAL pre-damage reference
 sum(h*r*V), never altered mid-path. w=1 reproduces the independent evaluator.
 Weighted version is a variational diagnostic, NOT exact Shojaei reproduction.
 """
 def __init__(self,X,b,w=None):
  self.X=X;self.b=b;self.i,self.j=b.T;self.r=np.linalg.norm(X[self.j]-X[self.i],axis=1);self.h=np.where(abs(self.r-DELTA)<=DX/2,(DELTA+DX/2-self.r)/DX,1.);self.m=np.bincount(b.ravel(),weights=np.repeat(self.h*self.r*V,2),minlength=len(X));self.w=np.ones(len(b)) if w is None else (w[self.i]+w[self.j])/2
 def evaluate(self,x):
  y=x[self.j]-x[self.i];length=np.linalg.norm(y,axis=1);e=length-self.r;e[np.abs(e)<2.2204e-16]=0
  theta=3*V*np.bincount(self.b.ravel(),weights=np.repeat(self.h*self.w*e,2),minlength=len(x))/self.m
  U=.5*(K-5*G/3)*V*np.dot(theta,theta)+7.5*G*V*V*np.sum(self.h*self.w/self.r*e**2*(1/self.m[self.i]+1/self.m[self.j]))
  magnitude=V*V*self.h*self.w*((3*K-5*G)*(theta[self.i]/self.m[self.i]+theta[self.j]/self.m[self.j])+15*G/self.r*(1/self.m[self.i]+1/self.m[self.j])*e)
  pair=magnitude[:,None]*y/length[:,None];force=np.stack([np.bincount(self.i,weights=pair[:,d],minlength=len(x))-np.bincount(self.j,weights=pair[:,d],minlength=len(x)) for d in range(3)],axis=1)
  return U,force,pair,theta

def families(b,n):
 f=[[] for _ in range(n)]
 for i,j in b:f[i].append(int(j));f[j].append(int(i))
 return f

def weight_audit(name,X,b):
 sys.path.insert(0,str(A/'sources/perifit'))
 import perifit.core.local_system as ls
 from perifit.core.targets import build_targets_osb
 original=ls._solve_local;conditions=[];residuals=[];regularizers=[];rowcounts=[];systems=[]
 def instrument(ML,targets,nF,rel,ab):
  norm=np.linalg.norm(ML,axis=1);valid=norm>1e-20;MN=ML.copy();MN[valid]/=norm[valid,None];gram=MN.T@MN;lam=max(ab,rel*np.max(abs(np.diag(gram))));conditions.append(np.linalg.cond(gram+lam*np.eye(10)));regularizers.append(lam);rowcounts.append(int(valid[nF:].sum()));systems.append((ML[nF:].copy(),targets.copy(),norm[nF:].copy()));return original(ML,targets,nF,rel,ab)
 ls._solve_local=instrument
 coords=(X-X.mean(0))/DELTA;vol=np.full(len(X),V/DELTA**3);fam=families(b,len(X));ds,es=build_targets_osb(1.);mat=lil_matrix((len(X),len(X)));rhs=np.empty(len(X));degenerate=[]
 for i,ids in enumerate(fam):
  ret=ls.build_local_system_osb(coords[ids]-coords[i],vol[ids],1.,ds,es)
  mat[i,i]=1
  if ret is None:rhs[i]=1;degenerate.append(i);continue
  beta,coupling=ret;mat[i,ids]=-coupling;rhs[i]=beta
 matrix=mat.tocsr();weights,info=bicgstab(matrix,rhs,x0=np.ones(len(X)),rtol=1e-10,maxiter=5000);ls._solve_local=original
 # Evaluate the actual 54 frozen upstream OSB operator rows, not README's 109 claim.
 before=[];after=[]
 from perifit.core.local_system import _dphi_all
 for i,ids in enumerate(fam):
  xi=coords[ids]-coords[i];r2=np.sum(xi**2,axis=1);dphi=_dphi_all(xi);op=[]
  for axis in range(3):
   for k in range(9):op.append(xi[:,axis]/r2*dphi[:,k]*vol[ids])
  for axis in range(3):
   for k in range(9):op.append(dphi[:,k]**2*xi[:,axis]**2/r2*vol[ids])
  op=np.array(op);target=np.r_[ds,es];scale=np.maximum(abs(target),np.sum(abs(op),axis=1));valid=scale>0
  before.append(np.linalg.norm(((op.sum(1)-target)/scale)[valid]));after.append(np.linalg.norm(((op@((weights[i]+weights[ids])/2)-target)/scale)[valid]))
 residual=np.linalg.norm(matrix@weights-rhs)/np.linalg.norm(rhs);wb=(weights[b[:,0]]+weights[b[:,1]])/2
 result=dict(source_commit='66b2192a9f95db29d16adefa162acb579759b925',model='osb; frozen source has54 operator rows, README says109',coordinate_units='nondimensional x/delta,V/delta^3, horizon1 to avoid source absolute1e-20 row cutoff silently dropping SI rows',family='Exact supplied reference topology, never geometric neighbor reconstruction across initial crack',weights=stats(weights),negative_nodes=int((weights<0).sum()),nonpositive_bonds=int((wb<=0).sum()),bond_weights=stats(wb),global_solver_info=int(info),global_relative_residual=residual,local_regularized_gram_condition=stats(conditions),local_operator_rows=stats(rowcounts),degenerate_nodes=degenerate,operator_residual_before=stats(before),operator_residual_after=stats(after),positivity_enforced_by_upstream=False,weight_clipping=False,physical_response_certified=False)
 np.savez_compressed(A/(name+'_WEIGHTS.npz'),X=X,bonds=b,weights=weights,bond_weights=wb,operator_residual_before=before,operator_residual_after=after,local_condition=conditions)
 return result,weights

def neighborhood(name,X,b,selected=None,removed=None,mask=None):
 e=Potential(X,b);n=len(X);xi=X[b[:,1]]-X[b[:,0]];unit=xi/e.r[:,None];second=np.zeros((n,3,3));odd=np.zeros((n,3));count=np.bincount(b.ravel(),minlength=n)
 for a in range(3):
  odd[:,a]=np.bincount(b[:,0],weights=e.h*unit[:,a],minlength=n)-np.bincount(b[:,1],weights=e.h*unit[:,a],minlength=n)
  for d in range(3):second[:,a,d]=3*np.bincount(b.ravel(),weights=np.repeat(e.h*e.r*V*unit[:,a]*unit[:,d],2),minlength=n)/e.m
 odd/=np.maximum(count,1)[:,None];anis=np.linalg.norm(second-np.eye(3),axis=(1,2));eig=np.linalg.eigvalsh(second)
 if selected is None:selected=np.arange(n)
 # Actual same-lattice full-horizon numerical stencil, not continuum m or fabricated reference nodes.
 q=np.stack(np.meshgrid(np.arange(-4,5),np.arange(-4,5),np.arange(-4,5),indexing='ij'),axis=-1).reshape(-1,3)*DX;r=np.linalg.norm(q,axis=1);r=r[(r>0)&(r<=DELTA)];h=np.where(abs(r-DELTA)<=DX/2,(DELTA+DX/2-r)/DX,1.);fullm=np.sum(h*r*V);fullcount=len(r)
 labels=np.full(n,'free_surface',dtype='U25');labels[count==fullcount]='interior'
 if mask is not None:
  grip=mask.any(1);labels[grip]='grip';corners=grip&((abs(X[:,0])>X[:,0].max()-DELTA)|(abs(X[:,1])>X[:,1].max()-DELTA));labels[corners]='grip_corner'
 if removed is not None:
  labels[np.unique(removed)]='initial_notch';tip=np.sqrt((X[:,0]+.0001)**2+X[:,2]**2)<DELTA;labels[tip]='notch_tip'
 if name=='V3':
  labels[np.setdiff1d(np.arange(n),selected)]='anchor'
 result=dict(nodes=len(selected),reference_full_stencil_neighbors=fullcount,full_stencil_weighted_volume_m4=fullm,neighbor_count=stats(count[selected]),weighted_volume_ratio=stats(e.m[selected]/fullm),direction_matrix_anisotropy=stats(anis[selected]),direction_matrix_eigenvalues=stats(eig[selected].ravel()),odd_direction_bias=stats(np.linalg.norm(odd[selected],axis=1)),regions={str(k):dict(count=int(np.sum(labels[selected]==k)),neighbors=stats(count[selected][labels[selected]==k]),m_ratio=stats((e.m/fullm)[selected][labels[selected]==k]),anisotropy=stats(anis[selected][labels[selected]==k])) for k in np.unique(labels[selected])},reference_m='Actual initial topology; m fixed before loading. V3 includes real typed anchor/interface neighbors. M5C excludes saved2570 initial notch bonds.',interpretation='Truncated-neighborhood operator diagnostic; nonzero affine boundary force is not by itself a boundary-condition error.')
 np.savez_compressed(A/(name+'_NEIGHBORHOOD.npz'),X=X,bonds=b,selected=selected,neighbors=count,m=e.m,full_m=fullm,second_moment=second,anisotropy=anis,odd=odd,labels=labels)
 return result

def mechanical(X,b,w):
 baseline=Potential(X,b);corr=Potential(X,b,w);u=np.c_[.003*X[:,0]+.001*X[:,2],-.001*X[:,1],.004*X[:,2]+.0005*(X[:,0]**2/DX)];x=X+u;result={};saved={}
 for name,model in [('original',baseline),('weighted_diagnostic',corr)]:
  U,f,pair,theta=model.evaluate(x);rng=np.random.default_rng(1241);d=rng.normal(size=x.shape);d/=np.linalg.norm(d);h=1e-9;fd=-(model.evaluate(x+h*d)[0]-model.evaluate(x-h*d)[0])/(2*h);project=np.sum(f*d);error=abs(fd-project)/max(abs(fd),abs(project))
  chosen=np.unravel_index(np.argmax(abs(f)),f.shape);xp=x.copy();xm=x.copy();xp[chosen]+=h;xm[chosen]-=h;fdcomp=-(model.evaluate(xp)[0]-model.evaluate(xm)[0])/(2*h);comp_error=abs(fdcomp-f[chosen])/max(abs(fdcomp),abs(f[chosen]));angle=.4;rot=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]]);ur,fr,_,_=model.evaluate(x@rot.T+np.array([.0001,-.0002,.0003]));u0,f0,_,_=model.evaluate(X@rot.T+np.array([.0001,-.0002,.0003]));
  result[name]=dict(U_J=U,force_L2_N=np.linalg.norm(f),max_node_force_N=np.linalg.norm(f,axis=1).max(),directional_FD_relative_error=error,component_FD_relative_error=comp_error,difference_h_m=h,force_sum_N=f.sum(0),current_torque_Nm=np.cross(x,f).sum(0),rigid_reference_energy_J=u0,rigid_reference_max_force_N=np.linalg.norm(f0,axis=1).max(),objectivity_energy_relative_error=abs(ur-U)/U,rotated_force_relative_error=np.linalg.norm(fr-f@rot.T)/np.linalg.norm(f));saved[name+'_force_N']=f;saved[name+'_theta']=theta
 # Baseline is independently compared to locked existing C++ evaluator (no compile/native change).
 from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy
 native=Energy({'X':X,'bonds':b},DX,ROOT/'analysis/material_macro_audit_v1/energy.so');un,fn=native.eval(np.ascontiguousarray(x),np.ones(len(b),np.uint8),True);ub,fb,_,_=baseline.evaluate(x);result['original_cpp_reference']=dict(relative_energy_error=abs(un-ub)/ub,relative_force_error=np.linalg.norm(fn-fb)/np.linalg.norm(fb),source_sha256=sha(ROOT/'analysis/material_macro_audit_v1/energy.cpp'))
 result['definition']='Frozen-m kernel-preserving expanded variational diagnostic with upstream geometry weights; NOT paper-exact operator or validated mechanical surface correction. theta uses w, all e^2 terms use w, m original fixed. Ordinary deviatoric-square rewrite acquires extra (5G/6)(1-m_weighted/m)*theta^2 term; omitting it changes the potential.'
 result['passes_FD']=all(result[k]['directional_FD_relative_error']<1e-5 and result[k]['component_FD_relative_error']<1e-5 for k in ['original','weighted_diagnostic']);np.savez_compressed(A/'FORCE_ENERGY_DIAGNOSTIC.npz',X=X,x=x,weights=w,**saved);return result,baseline

def paths(model):
 X=model.X;aa=np.c_[.003*X[:,0],np.zeros((len(X),2))];bb=np.c_[np.zeros((len(X),2)),.004*X[:,2]];corners=[X,X+aa,X+aa+bb,X+bb,X];gp,gw=np.polynomial.legendre.leggauss(8);work=np.zeros(len(model.b));hist=[work.copy()];total=[]
 for a,b in zip(corners[:-1],corners[1:]):
  dx=b-a;dy=dx[model.j]-dx[model.i];increment=np.zeros(len(work))
  for q,w in zip((gp+1)/2,gw/2):increment+=w*np.einsum('ij,ij->i',model.evaluate(a+q*dx)[2],dy)
  work+=increment;hist.append(work.copy());total.append(work.sum())
 energy=np.array([model.evaluate(x)[0] for x in corners]);hist=np.array(hist);scale=energy.max();best=int(np.argmax(abs(work)));straight=np.zeros(len(work));dx=aa+bb;dy=dx[model.j]-dx[model.i]
 for q,w in zip((gp+1)/2,gw/2):straight+=w*np.einsum('ij,ij->i',model.evaluate(X+q*dx)[2],dy)
 # A direct loading and its exact reverse: quadrature reversal gives opposite pair work.
 np.savez_compressed(A/'INTERACTION_PATH_WORK.npz',corner_energy_J=energy,all_pair_cumulative_work_J=hist,selected_bond_index=best,selected_bond_ids=model.b[best]+1,straight_pair_work_J=straight,selected_pair_cumulative_work_J=hist[:,best],corner_coordinates_m=np.array(corners))
 return dict(actual_paths='Two monotone orderings to same endpoint, their four-edge closed rectangle, plus straight load and exact reverse; four prescribed shape corners, 8-point path quadrature per segment. No equilibrium solves or failure.',total_work_at_corners_J=hist.sum(1),energy_at_corners_J=energy,max_total_work_energy_error_scaled=float(np.max(abs(hist.sum(1)-energy))/scale),closed_total_work_J=work.sum(),closed_sum_abs_pair_work_J=np.abs(work).sum(),closed_max_abs_pair_work_J=abs(work).max(),selected_bond_ids=model.b[best]+1,two_orderings_same_endpoint_total_error_J=float(hist[2].sum()+(hist[4]-hist[2]).sum()),straight_total_work_J=straight.sum(),straight_energy_error_relative=float(abs(straight.sum()-energy[2])/scale),double_counting='One unordered pair with Vi*Vj times directed-state difference; directed-all-pairs sum requires factor1/2.',conclusion='Total hyperelastic path work is recoverable and agrees with global U. Individual pair loop works need not vanish: path-defined interaction work is computable, independent unique pair elastic potential is NOT_ESTABLISHED. No Foster threshold authorized.')

def main():
 t=time.perf_counter();cpu=time.process_time();assert not (A/'A_DESIGN_FREEZE.json').exists();save(A/'A_DESIGN_FREEZE.json',dict(time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),geometries=['M5C2880','V3material1927'],weight_method='Pinned upstream OSB with explicit actual families, nondimensional units, no clipping/scan',weight_gates='Converged relative residual<1e-8; all bond weights positive for energy application; no bonds restored; large values and operator residual disclosed.',mechanical='At most M5C single deterministic field and its rigid transforms; finite differences1nm; frozen original m expanded variational diagnostic only if positive weights. This is an adaptation, not paper exact.',force_energy_relative_gate=1e-5,pathwork='Original undamaged M5C, four analytic corners and one straight path; Gauss8 integration, no dynamics; total work error<1e-8',production_changes=0))
 gm=dict(np.load(M5C/'geometry.npz'));gv=dict(np.load(V3/'V3_MATERIAL_POINTS_AND_BONDS.npz'));X=gm['X'];b=gm['bonds'];n=int(gv['clot_count']);bv=gv['bonds'][gv['bond_types']==1]
 diag={'M5C':neighborhood('M5C',X,b,removed=gm['removed_notch_bonds'],mask=gm['prescribed_mask']),'V3':neighborhood('V3',gv['X'],gv['bonds'],np.arange(n))};save(A/'SURFACE_NEIGHBORHOOD_AUDIT.json',diag)
 wm,w=weight_audit('M5C',X,b);wv,v=weight_audit('V3',gv['X'][:n],bv);save(A/'SURFACE_CORRECTION_FEASIBILITY.json',dict(M5C=wm,V3=wv,scope='Paper unit-kernel geometry weights only. Current native uses omega=1/r with horizon taper and frozen discrete m. Direct substitution not certified; V3 optimization omits anchors, so does not model wall attachment mechanics.',broken_bonds_restored=0,no_native_weight_reader=True))
 valid=wm['global_solver_info']==0 and wm['global_relative_residual']<1e-8 and wm['nonpositive_bonds']==0
 if valid:mech,model=mechanical(X,b,w)
 else:mech={'status':'NOT_RUN_INVALID_WEIGHTS','reason':wm};model=Potential(X,b)
 save(A/'FORCE_ENERGY_COMPATIBILITY.json',mech);pw=paths(model);save(A/'LPS_INTERACTION_WORK_AUDIT.json',pw)
 old={};base=ROOT/'verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z'
 for name in ['M1','M2']:
  h=np.genfromtxt(base/name/'history.csv',names=True,delimiter=',');U=h['U_bulk_J'];old[name]=dict(history_sha256=sha(base/name/'history.csv'),samples=len(h),peak_U_J=U.max(),final_U_J=U[-1],final_external_work_J=h['W_ext_J'][-1],final_kinetic_J=h['K_free_J'][-1],final_numerical_dissipation_J=h['D_num_J'][-1],final_energy_residual_J=h['R_energy_J'][-1],broken_bulk_max=int(h['broken_bulk'].max()),scope='Dynamic damped load-unload; W differs from elastic energy due to KE and recorded numerical damping; not material dissipation.')
 save(A/'EXISTING_LOAD_UNLOAD_AUDIT.json',old)
 states=dict(ORIGINAL_LPS_BASELINE_PROTECTED='PENDING_FINAL_SHA',SURFACE_NEIGHBORHOOD_DIAGNOSTIC='COMPLETED_ACTUAL_TOPOLOGIES',SURFACE_CORRECTION_WEIGHT_VALIDITY='GEOMETRIC_PREPROCESSING_ONLY_SEE_NEGATIVITY_AND_RESIDUALS',SURFACE_CORRECTION_FORCE_ENERGY_CONSISTENCY='PASS_VARIATIONAL_DIAGNOSTIC_ONLY' if mech.get('passes_FD') else 'NOT_ESTABLISHED',FOSTER_ENERGY_DEFINITION='FULL_TEXT_READ_DIMENSION_AND_PATH_AUDITED',LPS_INTERACTION_WORK_IDENTIFIABILITY='PATH_WORK_COMPUTABLE_UNIQUE_PAIR_POTENTIAL_NOT_ESTABLISHED',ENERGY_FAILURE_PORTABILITY='NOT_ESTABLISHED',STRENGTH_GC_INDEPENDENCE='ADDITIONAL_INDEPENDENT_PARAMETERS_REQUIRED',MATERIAL_GC_UNAVAILABLE=True,NATIVE_FAILURE_LAW_UNCHANGED=True)
 save(A/'FINAL_STATUS.json',dict(run_id=RUN,states=states,controlled_small_specimen_diagnostics=int(valid),new_native_solves=0,new_fracture=0));save(A/'RESOURCES.json',cost(t,cpu));print(json.dumps({'weight_M5C':wm['weights'],'negative_M5C':wm['negative_nodes'],'weight_V3':wv['weights'],'negative_V3':wv['negative_nodes'],'mechanical':valid,'work_error':pw['max_total_work_energy_error_scaled']},default=convert,indent=2))
if __name__=='__main__':main()
