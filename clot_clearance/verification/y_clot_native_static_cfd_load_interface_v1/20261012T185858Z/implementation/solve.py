"""One pre-registered full CFD dead-load equilibrium with independent energy/native verification."""
from common import *
import numpy as np,time,resource,csv,sys
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from native import Native
from mpi4py import MPI
from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy
assert MPI.COMM_WORLD.size==1
assert json.loads((O/'ATOM_ID_LOAD_APPLICATION_AUDIT.json').read_text())['status'].startswith('PASS')
assert json.loads((O/'UNLOADED_REFERENCE_AUDIT.json').read_text())['status'].startswith('PASS')
assert not (O/'STATIC_SOLVE_STARTED.json').exists();t=time.perf_counter();cpu=time.process_time();freeze=json.loads((O/'STATIC_LOAD_DESIGN_FREEZE.json').read_text());dump(O/'STATIC_SOLVE_STARTED.json',{'epoch':time.time(),'formal_equilibrium_attempt':1,'freeze_sha256':sha(O/'STATIC_LOAD_DESIGN_FREEZE.json')})
r=np.load(O/'REFERENCE_AND_LOAD.npz');X=r['X'];b=r['bonds'];types=r['bond_types'];fl=r['force_N'];x0=r['reference_point_m'];alive=np.ones(len(b),np.uint8);e=Energy({'X':X,'bonds':b},DX,R/'analysis/material_macro_audit_v1/energy.so');uscale=freeze['scales']['displacement_m'];escale=freeze['scales']['energy_J'];fscale=freeze['scales']['nodal_force_L2_N'];gates=freeze['gates'];lim=freeze['budget'];contact0=freeze['contact_no_activation_bounds'];n=Native(O/'native_static');native0,diag0=n.evaluate(X,'unloaded',False,initial=True)
U0,f0=e.eval(np.ascontiguousarray(X),alive,True);assert np.linalg.norm(f0,axis=1).max()<gates['unloaded_max_force_N'];dump(O/'INDEPENDENT_REFERENCE_ENERGY.json',{'U0_J':U0,'max_independent_force_N':float(np.linalg.norm(f0,axis=1).max()),'native_m_error':diag0['m_max_relative']})
numcalls=0;hist=[];states=[];last=None;maxtrial=0.;maxu=0.;end_reason=None
class BudgetStop(Exception):pass

def unpack(q):
 x=X.copy();x[:NC]+=q.reshape(NC,3)*uscale;return np.ascontiguousarray(x)
def evaluate(q):
 global numcalls,last,maxtrial,maxu
 numcalls+=1
 if numcalls>lim['max_function_evaluations'] or time.perf_counter()-t>lim['wall_seconds']:raise BudgetStop('FROZEN_COMPUTE_BUDGET')
 x=unpack(q);un=np.linalg.norm(x-X,axis=1);um=float(un.max());maxu=max(maxu,um)
 if um>lim['max_trial_displacement_m']:raise BudgetStop('TRIAL_DISPLACEMENT_STOP')
 stretch=(np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)-r['reference_lengths_m'])/r['reference_lengths_m'];util=stretch/np.choose(types-1,[.15,.140624859375,1.]);maxtrial=max(maxtrial,float(util.max()))
 if util.max()>=1:raise BudgetStop('NATURAL_BREAK_CRITERION_REACHED')
 # Conservative distance lower bounds; all trials remain in a proven inactive contact region.
 if np.min(r['wall_clearance_m']-un[:NC])<=20e-6 or contact0['nonbonded_min_reference_distance_m']-2*um<=DX:raise BudgetStop('CONTACT_REGION_ENTERED_REQUIRES_SEPARATE_MODEL_AUDIT')
 U,fi=e.eval(x,alive,True);work=float(np.sum(fl*(x-X)));Pi=U-work;res=fi+fl
 if not np.isfinite(Pi) or not np.isfinite(res).all():raise BudgetStop('NONFINITE')
 last=(q.copy(),x,U,fi,work,Pi,res)
 return Pi/escale,(-res[:NC]*uscale/escale).ravel()
def record(q):
 if last is None or not np.array_equal(q,last[0]):evaluate(q)
 _,x,U,fi,work,Pi,res=last;states.append((x.copy(),U,Pi));states[:]=states[-10:];reaction=-fi[NC:];fb=fl.sum(0)+reaction.sum(0);mc=np.cross(x-x0,fl).sum(0)+np.cross(X[NC:]-x0,reaction).sum(0);mr=np.cross(X-x0,fl).sum(0)+np.cross(X[NC:]-x0,reaction).sum(0)
 row={'iteration':len(hist),'function_evaluations':numcalls,'U_LPS_J':U,'dead_load_work_F_dot_u_J':work,'Pi_J':Pi,'free_force_max_N':float(np.linalg.norm(res[:NC],axis=1).max()),'free_force_L2_relative':float(np.linalg.norm(res[:NC])/fscale),'global_force_residual_N':float(np.linalg.norm(fb)),'current_moment_residual_Nm':float(np.linalg.norm(mc)),'reference_moment_residual_Nm':float(np.linalg.norm(mr)),'last10_max_position_change_m':float(np.max(np.linalg.norm(x-states[0][0],axis=1))) if len(states)==10 else None,'last10_U_span_J':float(max(z[1] for z in states)-min(z[1] for z in states)) if len(states)==10 else None,'max_displacement_m':float(np.linalg.norm(x[:NC]-X[:NC],axis=1).max())}
 hist.append(row)
 with (O/'STATIC_CONVERGENCE.csv').open('w') as fp:w=csv.DictWriter(fp,fieldnames=list(row));w.writeheader();w.writerows(hist)
 if len(hist)%50==0:print(json.dumps(row),flush=True)
 ok=row['free_force_max_N']<=gates['free_max_force_N'] and row['free_force_L2_relative']<=gates['free_force_L2_relative'] and row['global_force_residual_N']<=gates['global_force_abs_N'] and row['current_moment_residual_Nm']<=gates['global_moment_abs_Nm'] and row['last10_max_position_change_m'] is not None and row['last10_max_position_change_m']<=gates['last10_displacement_m'] and row['last10_U_span_J']<=gates['last10_U_span_J']
 return ok
q0=np.zeros(NC*3);evaluate(q0);record(q0)
def callback(intermediate_result):
 if record(intermediate_result.x):raise StopIteration
result=None
try:
 result=minimize(evaluate,q0,jac=True,method='L-BFGS-B',callback=callback,options=freeze['optimizer']);q=result.x;end_reason=str(result.message)
except BudgetStop as ex:
 end_reason=str(ex);q=last[0] if last is not None else q0
# The last valid evaluated candidate is retained on budget stop; this does not certify equilibrium.
x=unpack(q)
if last is None or not np.array_equal(last[0],q):evaluate(q)
if not np.array_equal(states[-1][0],x):record(q)
_,x,U,fi,work,Pi,res=last;native,diag=n.evaluate(x,'loaded_candidate',True);rn=native['force'];reaction=-rn[NC:];Fext=fl.sum(0);Freac=reaction.sum(0);Fb=Fext+Freac;Mext_ref=np.cross(X-x0,fl).sum(0);Mext_cur=np.cross(x-x0,fl).sum(0);Mr=np.cross(X[NC:]-x0,reaction).sum(0);Mbr=Mext_ref+Mr;Mbc=Mext_cur+Mr;native_error=rn-(fi+fl);force_cons=float(np.linalg.norm(native_error)/fscale);stretch=(np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)-r['reference_lengths_m'])/r['reference_lengths_m'];critical=np.choose(types-1,[.15,.140624859375,1.]);util=stretch/critical
Fs=freeze['scales']['face_absolute_force_sum_N'];Ms=Fs*freeze['scales']['length_m'];row=hist[-1]
checks={'free_force_absolute':float(np.linalg.norm(rn[:NC],axis=1).max())<=gates['free_max_force_N'],'free_force_relative':float(np.linalg.norm(rn[:NC])/fscale)<=gates['free_force_L2_relative'],'global_force_absolute':float(np.linalg.norm(Fb))<=gates['global_force_abs_N'],'global_force_relative':float(np.linalg.norm(Fb)/Fs)<=gates['global_force_scaled'],'global_current_moment_absolute':float(np.linalg.norm(Mbc))<=gates['global_moment_abs_Nm'],'global_current_moment_relative':float(np.linalg.norm(Mbc)/Ms)<=gates['global_moment_scaled'],'global_reference_moment_absolute':float(np.linalg.norm(Mbr))<=gates['global_moment_abs_Nm'],'global_reference_moment_relative':float(np.linalg.norm(Mbr)/Ms)<=gates['global_moment_scaled'],'position_plateau':row['last10_max_position_change_m'] is not None and row['last10_max_position_change_m']<=gates['last10_displacement_m'],'energy_plateau':row['last10_U_span_J'] is not None and row['last10_U_span_J']<=gates['last10_U_span_J'],'native_energy_gradient':force_cons<=gates['native_independent_L2_over_load'],'subcritical':float(util.max())<1,'topology':diag['new_breaks']==0 and diag['bonds']==len(b),'anchors_fixed':np.array_equal(x[NC:],X[NC:]),'contact_inactive':diag['contact_pairs']==0 and diag['wall_active_nodes']==0,'finite':bool(np.isfinite(x).all() and np.isfinite(rn).all())}
accepted=all(checks.values());disp=x-X;norm=np.linalg.norm(disp[:NC],axis=1);imax=int(np.argmax(norm));dump(O/'STATIC_EQUILIBRIUM_AUDIT.json',{'status':'PASS_STATIC_EQUILIBRIUM' if accepted else 'STATIC_EQUILIBRIUM_NOT_ESTABLISHED','accepted':accepted,'gates':checks,'criteria':gates,'iterations':len(hist)-1,'function_evaluations':numcalls,'optimizer_status':int(result.status) if result is not None else None,'optimizer_exit_reason':end_reason,'last':row,'max_native_free_force_N':float(np.linalg.norm(rn[:NC],axis=1).max()),'native_free_force_L2_over_load_L2':float(np.linalg.norm(rn[:NC])/fscale),'native_independent_error_L2_over_load_L2':force_cons,'native_max_force_difference_N':float(np.linalg.norm(native_error,axis=1).max()),'max_displacement_m':float(norm[imax]),'max_displacement_atom_id':imax+1,'max_displacement_reference_position_m':X[imax].tolist(),'mean_displacement_m':float(norm.mean()),'mean_displacement_vector_m':disp[:NC].mean(0).tolist(),'max_all_trial_utilization':maxtrial,'max_all_trial_displacement_m':maxu,'formal_equilibrium_attempts':1,'no_dynamic_integration':True,'candidate_file':'CLOT_STATIC_DISPLACEMENTS.npz','native_state':diag})
np.savez_compressed(O/'CLOT_STATIC_DISPLACEMENTS.npz',atom_id=np.arange(1,N+1),reference_coordinates_m=X,position_m=x,displacement_m=disp,force_CFD_N=fl,native_force_total_N=rn,independent_internal_force_N=fi,reference_m=r['reference_m'],native_m=native['m'],native_theta=native['theta'],accepted=np.array(accepted),anchors=np.arange(N)>=NC,support_node=np.r_[r['support_node'],np.zeros(N-NC,bool)],bonds=b,bond_types=types,stretch=stretch,critical_utilization=util,last10_positions_m=np.array([z[0] for z in states]),last10_energy_J=np.array([z[1] for z in states]))
np.savez_compressed(O/'CLOT_STATIC_REACTIONS.npz',anchor_atom_id=np.arange(NC+1,N+1),anchor_position_m=x[NC:],native_anchor_internal_force_N=rn[NC:],constraint_reaction_N=reaction,reference_point_m=x0,external_force_N=fl,wall_force_N=native['wall_force_N'],pair_contact_force_N=native['contact_force_N'],accepted=np.array(accepted))
dump(O/'FORCE_MOMENT_BALANCE.json',{'reference_point_m':x0.tolist(),'external_force_N':Fext.tolist(),'anchor_reaction_N':Freac.tolist(),'force_residual_N':Fb.tolist(),'force_residual_norm_N':float(np.linalg.norm(Fb)),'force_scaled_error':float(np.linalg.norm(Fb)/Fs),'external_reference_moment_Nm':Mext_ref.tolist(),'external_current_moment_Nm':Mext_cur.tolist(),'anchor_reaction_moment_Nm':Mr.tolist(),'reference_moment_residual_Nm':Mbr.tolist(),'current_moment_residual_Nm':Mbc.tolist(),'reference_moment_scaled_error':float(np.linalg.norm(Mbr)/Ms),'current_moment_scaled_error':float(np.linalg.norm(Mbc)/Ms),'configuration_moment_shift_Nm':(Mext_cur-Mext_ref).tolist(),'force_scale_N':Fs,'moment_scale_Nm':Ms,'current_configuration_is_physical_torque_balance':True,'reference_configuration_is_original_mapping_convention':True,'reaction_source':'Minus independently read native anchor internal force, not minus total external force','other_external_wall_force_N':native['wall_force_N'].sum(0).tolist()})
dump(O/'ANCHOR_REACTION_AUDIT.json',{'status':'PASS_NATIVE_INTERNAL_FORCE_CONSTRAINT_RECOVERY' if accepted else 'CANDIDATE_REACTION_ONLY','anchor_count':N-NC,'anchor_max_displacement_m':float(np.linalg.norm(disp[NC:],axis=1).max()),'reaction_sum_N':Freac.tolist(),'max_anchor_reaction_N':float(np.linalg.norm(reaction,axis=1).max()),'max_anchor_reaction_atom_id':int(np.argmax(np.linalg.norm(reaction,axis=1))+NC+1),'method':'Fixed all anchor coordinates in optimizer, no nve and no setforce clearing; reaction is negative of actual constrained-DOF native internal force','CFD_load_on_anchors_N':float(np.linalg.norm(fl[NC:])),'inferred_from_external_sum':False})
# Few independent energy-gradient/path checks at the candidate; zero additional equilibria.
direction=disp.copy();direction/=np.linalg.norm(direction);fd=[]
for eps in [1e-10,5e-11]:
 du=(e.eval(np.ascontiguousarray(x+eps*direction),alive)-e.eval(np.ascontiguousarray(x-eps*direction),alive))/(2*eps);pred=-float(np.sum(fi*direction));fd.append({'step_m':eps,'directional_gradient_N':du,'minus_force_dot_direction_N':pred,'error_over_load_L2':abs(du-pred)/fscale})
points,weights=np.polynomial.legendre.leggauss(5);pathwork=0.
for z,w in zip((points+1)/2,weights/2):
 _,fp=e.eval(np.ascontiguousarray(X+z*disp),alive,True);pathwork-=w*float(np.sum(fp*disp))
dump(O/'INDEPENDENT_LPS_ENERGY_AUDIT.json',{'status':'PASS_GRADIENT_CHECK' if force_cons<=gates['native_independent_L2_over_load'] and max(a['error_over_load_L2'] for a in fd)<=gates['FD_over_load_L2'] else 'FORCE_ENERGY_MISMATCH','U_reference_J':U0,'U_loaded_J':U,'DeltaU_J':U-U0,'constant_dead_load_F_dot_u_J':work,'external_potential_J':-work,'total_potential_Pi_J':Pi,'independent_internal_path_work_J':pathwork,'path_work_minus_DeltaU_J':pathwork-(U-U0),'half_F_dot_u_J':.5*work,'DeltaU_minus_half_F_dot_u_J':U-U0-.5*work,'half_work_scope':'Small-displacement linearity diagnostic only; no separately simulated ramped-load path','finite_difference':fd,'native_independent_L2_over_load':force_cons,'native_PE_times_volume_J_excluded_from_objective':diag['native_bond_PE_times_V_J'],'pair_contact_energy_J':diag['pair_contact_PE_J'],'wall_contact_energy_J':0.,'numerical_damping_work_J':0.,'Gc':'NOT_CALIBRATED_NO_FRACTURE','kinetic_energy':'No time integrator; KE=0 is not used as evidence of static equilibrium'})
n.close();qres=resource.getrusage(resource.RUSAGE_SELF);dump(O/'SOLVE_PERFORMANCE.json',{'formal_equilibrium_attempts':1,'optimizer_iterations':len(hist)-1,'function_evaluations':numcalls,'native_instances':1,'native_force_states':2,'dynamic_steps':0,'wall_s':time.perf_counter()-t,'CPU_s':time.process_time()-cpu,'peak_RSS_KiB':qres.ru_maxrss,'MPI_ranks':1,'OMP_BLAS_threads':1,'GPU_used':False,'exit_reason':end_reason});print(json.dumps({'accepted':accepted,'checks':checks,'last':row,'native_force_consistency':force_cons,'maxu_m':float(norm[imax]),'wall_s':time.perf_counter()-t}),flush=True)
