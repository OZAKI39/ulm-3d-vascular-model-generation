from common import *
from model import SofteningLPS
import time,resource
start=time.perf_counter();freeze=json.loads((O/'DESIGN_FREEZE.json').read_text());m=SofteningLPS(np.load(O/'provenance/reference.npz'));X=m.X
assert not (O/'PATH_RESULTS.json').exists()
def path(knots,counts,events=()):
 ss=[knots[0]];ph=[0]
 for n,(a,b,count) in enumerate(zip(knots[:-1],knots[1:],counts)):
  q=np.linspace(a,b,count+1)[1:]
  ev=[v for v in events if min(a,b)<v<max(a,b)]
  q=np.unique(np.r_[q,ev]);q=q if b>a else q[::-1];ss.extend(q);ph.extend([n+1]*len(q))
 return np.array(ss),np.array(ph)
primary,phase=path(freeze['path']['knots'],freeze['path']['intervals'],freeze['path']['event_points'])
aux,aphase=path(freeze['auxiliary']['knots_over_dx'],freeze['auxiliary']['intervals'])
sign=np.where(X[:,0]<-30e-6,-.5,.5);auxv=np.zeros_like(X);auxv[:,0]=sign
rng=np.random.default_rng(20261013);summaries={}
for name,sc,ss,ph,Vvec in [('instant_reference',None,primary,phase,X/DX),('soft_sc018',.18,primary,phase,X/DX),('soft_sc024',.24,primary,phase,X/DX),('planar_partial_sc024',.24,aux,aphase,auxv)]:
 t0=time.perf_counter();n=len(ss);kappa=np.zeros(m.nb);d=np.zeros(m.nb);uprev=0.;prevx=X.copy();W=0.;Diss=0.;records=[];stored={k:[] for k in ['x','force','d','kappa','theta']};maxhealing=0.;mindiss=0.;maxfailed=0.;maxfd=0.;maxmoment=0.
 evaluator=lambda xx,dd:m.baseline(xx,1-dd) if sc is None else m.evaluate(xx,dd)
 for it,s in enumerate(ss):
  x=X+s*DX*Vvec;old=evaluator(x,d);dW=0.;dx=x-prevx
  for g in [.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)]:
   fg=evaluator(prevx+g*dx,d)['force'];dW+=-.5*float(np.sum(fg*dx))
  kapnew=np.maximum(kappa,np.maximum(old['e']/m.r,0))
  if sc is None:dnew=(kapnew>S0+1e-12).astype(float)
  else:kapnew,dnew=m.update(x,kappa,S0,sc)
  state=evaluator(x,dnew);drop=old['U']-state['U'];W+=dW;Diss+=drop
  maxhealing=max(maxhealing,float(np.max(d-dnew)));mindiss=min(mindiss,drop)
  f=state['force'];Q=-np.sum(f*Vvec);Rwork=W-(state['U']+Diss)
  failed=dnew==1.;maxfailed=max(maxfailed,float(np.max(abs(state['t'][failed]))) if failed.any() else 0.)
  maxfd=max(maxfd,float(np.linalg.norm(f.sum(0))));maxmoment=max(maxmoment,float(np.linalg.norm(np.cross(x,f).sum(0))))
  records.append([s,s*DX,Q,state['U'],W,Diss,Rwork,drop,dW,float(dnew.mean()),float(dnew.max()),float(kapnew.max()),int(failed.sum()),int(((dnew>1e-12)&(dnew<1)).sum())])
  for key,val in [('x',x),('force',f),('d',dnew),('kappa',kapnew),('theta',state['theta'])]:stored[key].append(val.copy())
  d=dnew;kappa=kapnew;prevx=x;uprev=state['U']
 arr=np.array(records);cols=['s','controlled_displacement_m','reaction_N','stored_J','external_work_J','damage_dissipation_J','balance_residual_J','signed_damage_increment_J','work_increment_J','mean_damage','max_damage','max_history_stretch','fully_failed_bonds','partial_bonds']
 np.savez_compressed(O/f'{name}.npz',**{k:np.array(v) for k,v in stored.items()},series=arr,columns=cols,phase=ph,reference_X=X,bonds=m.b)
 force_scale=float(abs(arr[:,2]).max());energy_scale=float(max(arr[:,4].max(),1e-30));at0=np.flatnonzero((arr[:,0]==0)&(ph>1));peakidx=int(np.argmax(arr[:,2]));onsetidx=np.flatnonzero(arr[:,10]>1e-8);completeidx=np.flatnonzero(arr[:,12]==m.nb)
 # Compare unload and reload states at the same s=.165, excluding later damage extension.
 turn=np.flatnonzero(np.isclose(ss,.165,rtol=0,atol=1e-14)) if name!='planar_partial_sc024' else np.flatnonzero(np.isclose(ss,.18,rtol=0,atol=1e-14))
 turn=turn[(ph[turn]>=2)&(ph[turn]<=5)] if name!='planar_partial_sc024' else turn
 reload_error=float(np.max(abs(stored['d'][turn[0]]-stored['d'][turn[-1]]))) if len(turn)>1 else None
 expected=None
 if name!='planar_partial_sc024':expected=m.baseline(X*(1+S0))['U']*(1 if sc is None else sc/S0)
 summary=dict(name=name,s_complete=sc,states=n,wall_s=time.perf_counter()-t0,peak_reaction_N=float(arr[peakidx,2]),peak_strain=float(ss[peakidx]),peak_local_edge_force_N=None,first_observed_damage_strain=float(ss[onsetidx[0]]) if len(onsetidx) else None,first_complete_strain=float(ss[completeidx[0]]) if len(completeidx) else None,total_dissipation_J=float(arr[-1,5]),final_energy_J=float(arr[-1,3]),final_reaction_N=float(arr[-1,2]),max_balance_relative=float(abs(arr[:,6]).max()/energy_scale),max_damage_decrease=maxhealing,min_signed_damage_increment_J=mindiss,max_force_failed_edges_N=maxfailed,max_net_force_N=maxfd,max_torque_Nm=maxmoment,reload_damage_error=reload_error,uniform_analytic_dissipation_J=expected,uniform_splitting_relative_error=abs(arr[-1,5]-expected)/expected if expected else None,zero_displacement_energies_J=arr[at0,3],no_physical_time=True,all_nodes_controlled=True,equilibrium_solved=False)
 if name=='planar_partial_sc024':summary.update(protocol_area_m2=freeze['auxiliary']['area_m2'],PROTOCOL_DISSIPATION_PER_AREA_J_m2=float(arr[-1,5]/freeze['auxiliary']['area_m2']),fully_separated=False,final_failed_bonds=int(arr[-1,12]),final_partial_bonds=int(arr[-1,13]))
 # Local edge peak, exact on first loading at onset for the uniform protocol.
 if name!='planar_partial_sc024':summary['peak_local_edge_force_N']=float(np.max(abs(m.baseline(X*(1+S0))['t'])))
 summaries[name]=summary;dump(O/(name+'_SUMMARY.json'),summary);print(name,summary['total_dissipation_J'],summary['max_balance_relative'],flush=True)
 if mindiss < -1e-20 or maxhealing>1e-14 or maxfailed>1e-20:raise RuntimeError('Constitutive stop gate; all failed evidence retained')
 if summary['max_balance_relative']>2e-8:raise RuntimeError('Work balance failed')
 if expected and summary['uniform_splitting_relative_error']>.012:raise RuntimeError('Uniform splitting accuracy failed')
dump(O/'PATH_RESULTS.json',summaries)
# One deterministic analytic work budget, on the actual same reference geometry.
# U_d=1/2 e^T H_d e with 0<=H_d<=H: |Q| <= sqrt(e^T H e)*sqrt(eprime^T H eprime).
a,b=m.b.T;dy0=X[b]-X[a];dv=auxv[b]-auxv[a]
quad,w=np.polynomial.legendre.leggauss(48);limit=.245*DX;bound=0.
for delta,weight in zip((quad+1)*limit/2,w*limit/2):
 y=dy0+delta*dv;rr=np.linalg.norm(y,axis=1);e=rr-m.r;ep=np.sum(y*dv,axis=1)/rr
 en=np.dot(m.D,e*e)+m.A*np.dot(m.L@e,m.L@e)
 ed=np.dot(m.D,ep*ep)+m.A*np.dot(m.L@ep,m.L@ep)
 bound+=weight*np.sqrt(en*ed)
area=freeze['auxiliary']['area_m2']
dump(O/'ENERGY_TARGET_FEASIBILITY.json',dict(status='HIGH_ENERGY_TARGET_NOT_YET_ADMISSIBLE',scope='TARGET_ENERGY_OUTSIDE_CURRENT_MODEL_SCOPE',protocol='monotone actual 341-node rigid-half opening, max bond stretch .245',area_m2=area,opening_limit_m=limit,controlled_work_upper_bound_J=bound,work_upper_bound_J_m2=bound/area,target_lower_J_m2=6,target_over_bound=6/(bound/area),proof='Cauchy-Schwarz in positive H_d metric and 0 <= H_d <= H, integrated over prescribed path. No Gc or extrapolated sc formula.',physical_validity='No experimentally certified strain ceiling; .245 is a generous numerical diagnostic ceiling, not proven small-strain validity.',Gc_J_m2=None,G0_THEORY_J_m2=.07170070833333,G0_DISCRETE_J_m2=.06969623734234))
dump(O/'PATH_RESOURCES.json',dict(wall_s=time.perf_counter()-start,cpu_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,formal_path_runs=4,candidate_mechanisms=1,energy_parameters=[.18,.24],threads=1,GPU_used=False))
print('PATHS_COMPLETE',flush=True)
