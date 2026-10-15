from common import *
import sys,time,resource,csv,ctypes,traceback
sys.path.insert(0,str(R))
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator,cg
from analysis.dual_track_v1.core_a import Energy,delta_y,contact
T0=time.perf_counter();main=json.loads((O/'MAIN_SOLVER_EXIT.json').read_text())
assert main['error'] and 'MAIN_WALL_BUDGET' in main['error']
oldrec=json.loads((O/'RECOVERY_SOLVER_EXIT.json').read_text());assert oldrec['invocations']==0
assert not (O/'RECOVERY_STARTED.json').exists()
import shutil
for name in ['STATIC_SOLVER_TRACE.csv','function_evaluations.csv','INDEPENDENT_EQUILIBRIUM_CANDIDATE.json','final_candidate.npz','RECOVERY_SOLVER_EXIT.json']:
 shutil.copy2(O/name,O/('before_recovery_'+name))
freeze=json.loads((O/'LOW_LOAD_DESIGN_FREEZE.json').read_text());plan=json.loads((O/'STATIC_SOLVER_RECOVERY_PLAN.json').read_text());gates=freeze['gates'];T=T0-main['wall_s']
dump(O/'RECOVERY_STARTED.json',dict(epoch=time.time(),frozen_plan_sha256=sha(O/'STATIC_SOLVER_RECOVERY_PLAN.json'),source_sha256=sha(__file__),main_stopped='MAIN_WALL_BUDGET',reserved_recovery_invocation=1,scope='Use the already-frozen 120s recovery after the one bounded main solve. Fix an orchestration-only omission that unnecessarily excluded a main wall-budget exit. No change to physics, force, plateau, tangent, line search or recovery budget. Preserve all original stop records.',original_driver_sha256=sha(Path(__file__).with_name('solve.py'))))
rt=np.load(O/'reference_topology.npz');X=rt['X'];b=rt['bonds'];alive=np.ascontiguousarray(rt['alive'][0],dtype=np.uint8);am=alive.astype(bool);g=dict(X=X,bonds=b);e=Energy(g);fext=np.load(O/'load_identity.npz')['external_force_N'];Fnorm=np.linalg.norm(fext)
xs=plan['main']['coordinate_scale_m'];es=plan['main']['energy_scale_J'];calls=0;rows=[];hist=[];phase='recovery';state_dir=O/'accepted_optimizer_states';state_dir.mkdir(exist_ok=True);last=None;recovery_start=T0;Hcalls=0;current=None
def project(v):
 a=v.reshape(X.shape);return (a-a.mean(0)).reshape(v.shape)
def budget():
 elapsed=time.perf_counter()-T
 if elapsed>plan['budgets']['total_solver_wall_s']:raise RuntimeError('TOTAL_WALL_BUDGET')
 if phase=='main' and elapsed>plan['main']['wall_s']:raise RuntimeError('MAIN_WALL_BUDGET')
 if phase=='recovery' and time.perf_counter()-recovery_start>plan['recovery']['wall_s']:raise RuntimeError('RECOVERY_WALL_BUDGET')
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024>plan['budgets']['peak_RSS_MiB']:raise RuntimeError('MEMORY_BUDGET')
 if len(rows)*X.nbytes*2>plan['budgets']['output_MiB']*1024**2:raise RuntimeError('OUTPUT_BUDGET_CONSERVATIVE')
eval_stream=(O/'function_evaluations.csv').open('a',buffering=1);ew=csv.writer(eval_stream);None # Existing function log header retained
def evaluate(x):
 global calls,current
 budget();U,fi=e.eval(x,alive,True);Uc,fc,cp=contact(x,g,alive);fr=fi+fc+fext;work=float(np.sum(fext*(x-X)));calls+=1
 current=dict(x=x.copy(),U=float(U),Uc=float(Uc),Fdotu=work,Pi=float(U+Uc-work),fi=fi,fc=fc,res=fr,max_force=float(np.linalg.norm(fr,axis=1).max()),relative_force=float(np.linalg.norm(fr)/Fnorm),contact_pairs=len(cp))
 ew.writerow([calls,phase,time.perf_counter()-T,U,Uc,work,current['Pi'],current['max_force'],current['relative_force']]);return current
def fun(q):
 d=evaluate(X+project(q.reshape(X.shape))*xs);return d['Pi']/es,-project(d['res']).ravel()*xs/es
trace=(O/'STATIC_SOLVER_TRACE.csv').open('a',buffering=1);tw=None
def accept(x,stage,alpha=None,CG_info=None,HVP_FD=None):
 global last,tw
 d=current if current is not None and np.array_equal(current['x'],x) else evaluate(x)
 if last is not None and np.array_equal(last['x'],x):return False
 s=np.linalg.norm(delta_y(x,b),axis=1)/e.r-1
 hist.append((x.copy(),d['U']));hist[:]=hist[-10:];h=hist[-gates['plateau_accepted_distinct_states']:]
 span=(max(v[1] for v in h)-min(v[1] for v in h))/max(abs(d['U']),1e-30) if len(h)>=3 else None
 pos=float(max(np.linalg.norm(x-v[0],axis=1).max() for v in h)) if len(h)>=3 else None
 row=dict(accepted_state=len(rows),phase=stage,evaluations=calls,evaluations_since_last=calls-rows[-1]['evaluations'] if rows else calls,wall_s=time.perf_counter()-T,U_LPS_J=d['U'],contact_U_J=d['Uc'],Fdotu_J=d['Fdotu'],Pi_J=d['Pi'],max_force_N=d['max_force'],relative_force=d['relative_force'],max_positive_stretch=float(max(0,s[am].max())),step_max_m=float(np.linalg.norm(x-last['x'],axis=1).max()) if last is not None else 0.,Pi_decrease_J=last['Pi']-d['Pi'] if last is not None else 0.,energy_span_last3_relative=span,position_span_last3_m=pos,energy_span_last10_relative=(max(v[1] for v in hist)-min(v[1] for v in hist))/max(abs(d['U']),1e-30) if len(hist)==10 else None,position_span_last10_m=float(max(np.linalg.norm(x-v[0],axis=1).max() for v in hist)) if len(hist)==10 else None,global_force_relative=float(np.linalg.norm(d['res'].sum(0))/freeze['F_face_N']),global_torque_relative=float(np.linalg.norm(np.cross(x,d['res']).sum(0))/(freeze['F_face_N']*.005)),external_torque_relative=float(np.linalg.norm(np.cross(x,fext).sum(0))/(freeze['F_face_N']*.005)),alpha=alpha,CG_info=CG_info,HVP_FD_relative=HVP_FD,HVP_calls=Hcalls,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 if tw is None:tw=csv.DictWriter(trace,fieldnames=list(row));tw.writeheader()
 tw.writerow(row);rows.append(row);last=d
 np.savez_compressed(state_dir/f'state_{len(rows)-1:05d}.npz',x=x,residual_force_N=d['res'])
 if len(rows)%100==0 or stage=='recovery':print('accepted',len(rows)-1,stage,'evaluations',calls,'maxforce',row['max_force_N'],'rel',row['relative_force'],'wall',round(row['wall_s'],2),flush=True)
 return True
def mechanical_pass(row):
 return row['global_force_relative']<=gates['global_force_relative_to_face'] and row['global_torque_relative']<=gates['global_torque_relative_to_face_times_W'] and row['external_torque_relative']<=gates['global_torque_relative_to_face_times_W'] and row['max_force_N']<=gates['force_max_abs_N'] and row['relative_force']<=gates['force_L2_relative'] and row['max_positive_stretch']<=gates['max_positive_stretch'] and row['contact_U_J']<=gates['contact_U_max_J'] and row['energy_span_last3_relative'] is not None and row['energy_span_last3_relative']<=gates['energy_span_relative'] and row['position_span_last3_m']<=gates['position_span_max_m']
def callback(q):
 accept(X+project(q.reshape(X.shape))*xs,'main')
 if mechanical_pass(rows[-1]):raise StopIteration

with (O/'before_recovery_STATIC_SOLVER_TRACE.csv').open() as f:
 reader=csv.DictReader(f);fields=reader.fieldnames
 for rr in reader:
  rows.append({k:(v if k=='phase' else None if v=='' else float(v)) for k,v in rr.items()})
tw=csv.DictWriter(trace,fieldnames=fields)
calls=int(rows[-1]['evaluations'])
for row in rows[-10:]:
 saved=np.load(state_dir/f"state_{int(row['accepted_state']):05d}.npz");hist.append((saved['x'],row['U_LPS_J']))
d=np.load(O/'before_recovery_final_candidate.npz');row=rows[-1]
last=dict(x=d['x'],U=row['U_LPS_J'],Uc=row['contact_U_J'],Fdotu=row['Fdotu_J'],Pi=row['Pi_J'],fi=d['internal_LPS_force_N'],fc=d['contact_force_N'],res=d['residual_force_N'],max_force=row['max_force_N'],relative_force=row['relative_force'])
current=last;main_error=None
recovery=dict(invocations=0,status='NOT_NEEDED',steps=[]);rec_error=None
lib=ctypes.CDLL(str(R/'analysis/dual_track_v1/hvp_y.so'));hv=lib.hessian_vector;hv.restype=None;hv.argtypes=[ctypes.c_int]*2+[ctypes.c_void_p]*6+[ctypes.c_double]*3+[ctypes.c_void_p]*2
def H(x,v):
 global Hcalls
 budget();Hcalls+=1
 if Hcalls>plan['recovery']['maximum_total_HVP']:raise RuntimeError('HVP_TOTAL_BUDGET')
 v=np.ascontiguousarray(v.reshape(X.shape));out=np.zeros_like(x);hv(len(x),len(b),*[a.ctypes.data for a in [x,e.b,alive,e.r,e.h,e.m]],V,K,G,v.ctypes.data,out.ctypes.data);return out
if not mechanical_pass(rows[-1]) and main_error is None:
 phase='recovery';recovery_start=time.perf_counter();recovery.update(invocations=1,status='RUNNING',before_accepted_state=len(rows)-1)
 try:
  for k in range(plan['recovery']['maximum_Newton_iterations']):
   x=np.ascontiguousarray(last['x']);d=evaluate(x)
   if d['Uc']>gates['contact_U_max_J']:raise RuntimeError('NONNEGLIGIBLE_CONTACT_TANGENT_NOT_IMPLEMENTED')
   v=np.random.default_rng(928).normal(size=X.shape);v=project(v);v/=np.linalg.norm(v);h=plan['recovery']['FD_step_m'];fd=-(e.eval(x+h*v,alive,True)[1]-e.eval(x-h*v,alive,True)[1])/(2*h);hv0=H(x,v);err=float(np.linalg.norm(fd-hv0)/np.linalg.norm(fd))
   if err>plan['recovery']['FD_relative_gate']:raise RuntimeError('HVP_FD_FAILED')
   def matvec(z):
    vv=z.reshape(X.shape);return (project(H(x,project(vv)))+np.broadcast_to(vv.mean(0),X.shape)).ravel()
   rhs=project(d['res']).ravel();dq,code=cg(LinearOperator((X.size,X.size),matvec=matvec,dtype=float),rhs,rtol=plan['recovery']['CG_rtol'],atol=plan['recovery']['CG_atol_N'],maxiter=plan['recovery']['CG_maxiter_each'])
   if code!=0:raise RuntimeError('CG_STOP_'+str(code))
   accepted=False
   for alpha in plan['recovery']['line_search_alphas']:
    xt=X+project(x-X+alpha*project(dq.reshape(X.shape)));trial=evaluate(xt)
    if np.linalg.norm(trial['res'])<np.linalg.norm(d['res']) and trial['Pi']<=d['Pi']+100*np.finfo(float).eps*max(abs(d['U']),abs(d['Pi']),1e-30):
     accepted=accept(xt,'recovery',alpha,code,err)
     if accepted:break
   recovery['steps'].append(dict(iteration=k,CG_info=code,HVP_FD_relative=err,HVP_calls=Hcalls,accepted=accepted,alpha=alpha))
   if not accepted:raise RuntimeError('NO_DISTINCT_RESIDUAL_REDUCING_CORRECTION')
   if mechanical_pass(rows[-1]):recovery['status']='PHYSICAL_GATES_MET';break
  else:recovery['status']='NEWTON_ITERATION_LIMIT'
 except Exception:rec_error=traceback.format_exc();recovery['status']='STOPPED_WITH_FAILURE'
 recovery.update(error=rec_error,wall_s=time.perf_counter()-recovery_start,after_accepted_state=len(rows)-1,HVP_calls=Hcalls)
elif main_error is not None:recovery.update(status='NOT_RUN_MAIN_RESOURCE_OR_IMPLEMENTATION_FAILURE')
dump(O/'RECOVERY_SOLVER_EXIT.json',recovery)
# Save the last ACCEPTED optimizer state, never a rejected line-search trial.
phase='finalization';x=last['x'];U,fi=e.eval(x,alive,True);Uc,fc,cp=contact(x,g,alive);res=fi+fc+fext;s=np.linalg.norm(delta_y(x,b),axis=1)/e.r-1
np.savez_compressed(O/'final_candidate.npz',X=X,x=x,atom_id=np.arange(1,len(X)+1),bonds=b,alive=alive,reference_m=e.m,theta=e.theta(x,alive),stretch=s,external_force_N=fext,internal_LPS_force_N=fi,contact_force_N=fc,residual_force_N=res,displacement_m=x-X)
dump(O/'INDEPENDENT_EQUILIBRIUM_CANDIDATE.json',dict(status='INDEPENDENT_GATES_PASS_NATIVE_PENDING' if mechanical_pass(rows[-1]) else 'STATIC_STATE_NOT_ACCEPTED',last=rows[-1],main=main,recovery=recovery,new_bond_breaks=0,global_external_force_N=fext.sum(0),global_internal_force_N=(fi+fc).sum(0),global_residual_N=res.sum(0),current_external_torque_N_m=np.cross(x,fext).sum(0),current_residual_torque_N_m=np.cross(x,res).sum(0),U0_LPS_J=e.eval(X,alive),wall_s=time.perf_counter()-T,CPU_user_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime,CPU_system_s=resource.getrusage(resource.RUSAGE_SELF).ru_stime,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,evaluations=calls,accepted_states=len(rows),HVP_calls=Hcalls))
trace.close();eval_stream.close();print('SOLVER_FINISHED',rows[-1],flush=True)
