from common import *
import sys,time,resource,csv,traceback
sys.path.insert(0,str(R))
from scipy.optimize import minimize
from analysis.dual_track_v1.core_a import Energy,delta_y,contact
start=time.perf_counter();freeze=json.loads((O/'PRECHECK_FREEZE.json').read_text());assert not (O/'PRECHECK_STARTED.json').exists();dump(O/'PRECHECK_STARTED.json',{'epoch':time.time(),'freeze_sha256':sha(O/'PRECHECK_FREEZE.json')})
rt=np.load(O/'reference_topology.npz');X=rt['X'];b=rt['bonds'];alive=np.ascontiguousarray(rt['alive'][0],dtype=np.uint8);g={'X':X,'bonds':b};e=Energy(g);assert np.max(abs(e.m-rt['m'])/rt['m'])<1e-12
fext=np.zeros_like(X);fext[rt['high'],2]=freeze['force_per_loaded_node_N'];fext[rt['low'],2]=-freeze['force_per_loaded_node_N'];xs=1e-6;es=1e-12;evals=0;last=np.zeros(X.size);hist=[]
stream=(O/'precheck_iterations.csv').open('w',buffering=1);writer=csv.writer(stream);writer.writerow(['evaluation','wall_s','U_LPS_J','external_potential_negative_J','total_potential_J','residual_force_max_N','residual_L2_relative','mean_opening_m'])
U0,f0=e.eval(X,alive,True);uc0,fc0,cp0=contact(X,g,alive)
dump(O/'INDEPENDENT_INITIAL_STATE.json',{'U_LPS_J':U0,'contact_U_J':uc0,'max_internal_force_N':float(np.linalg.norm(f0+fc0,axis=1).max()),'velocity_m_s':0.,'zero_load_reaction_N':0.,'source':'Independent reference equilibrium, pending native confirmation'})
def fun(q):
 global evals,last
 if time.perf_counter()-start>freeze['budget']['wall_s']:raise RuntimeError('PRECHECK_WALL_BUDGET')
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024>freeze['budget']['peak_RSS_MiB']:raise RuntimeError('PRECHECK_MEMORY_BUDGET')
 x=X+q.reshape(X.shape)*xs;u,f=e.eval(x,alive,True);res=f+fext;potext=-np.sum(fext*(x-X));evals+=1;last=q.copy()
 row=[evals,time.perf_counter()-start,u,potext,u+potext,float(np.linalg.norm(res,axis=1).max()),float(np.linalg.norm(res)/np.linalg.norm(fext)),float((x[rt['high'],2]-X[rt['high'],2]).mean()-(x[rt['low'],2]-X[rt['low'],2]).mean())];writer.writerow(row)
 if evals%100==0:print('eval',evals,'wall',round(row[1],2),'maxres',row[5],flush=True)
 return (u+potext)/es,-res.ravel()*xs/es
err=None;opt=None
try:opt=minimize(fun,last,jac=True,method='L-BFGS-B',options={k:v for k,v in freeze['optimizer'].items() if k in ['maxiter','maxfun','ftol','gtol','maxcor','maxls']})
except Exception:err=traceback.format_exc()
x=X+(opt.x if opt is not None else last).reshape(X.shape)*xs;u,f=e.eval(x,alive,True);uc,fc,cp=contact(x,g,alive);res=f+fc+fext;s=np.linalg.norm(delta_y(x,b),axis=1)/e.r-1;theta=e.theta(x,alive)
d=delta_y(X,b);mid=X[b[:,0]]+.5*d;labels=np.full(len(b),4,np.int8);tip=np.hypot(mid[:,0]+.001,mid[:,2])<=DELTA;mouth=(mid[:,0]-X[:,0].min()<=DELTA)&(abs(mid[:,2])<=DELTA);grip=X[:,2].max()-abs(mid[:,2])<=DELTA;edge=np.minimum(mid[:,0]-X[:,0].min(),X[:,0].max()-mid[:,0])<=DELTA
labels[edge]=3;labels[grip]=2;labels[mouth]=1;labels[tip]=0;names=['notch_front','notch_mouth','loading_ends','other_X_free_surface','interior'];regions={}
for k,name in enumerate(names):
 idx=np.flatnonzero((labels==k)&alive.astype(bool));j=idx[np.argmax(s[idx])];regions[name]={'bonds':len(idx),'max_stretch':float(s[j]),'max_utilization':float(s[j]/S00),'reference_bond_id':int(j+1),'atom_ids':b[j]+1,'reference_midpoint_m':mid[j]}
gates={'residual_abs':float(np.linalg.norm(res,axis=1).max())<=freeze['gates']['force_max_abs_N'],'residual_relative':float(np.linalg.norm(res)/np.linalg.norm(fext))<=freeze['gates']['force_L2_relative_to_external'],'small_strain':float(s[alive.astype(bool)].max())<freeze['gates']['max_bond_stretch'],'contact_negligible':uc<freeze['gates']['noncontact_equilibrium_requires_final_contact_U_J_below'],'front_leads':regions['notch_front']['max_stretch']>max(v['max_stretch'] for n,v in regions.items() if n!='notch_front'),'resource_and_solver_completed':err is None}
np.savez_compressed(O/'elastic_precheck.npz',X=X,x=x,bonds=b,alive=alive,theta=theta,reference_m=e.m,external_force_N=fext,internal_LPS_force_N=f,contact_force_N=fc,stretch=s,bond_regions=labels,precut_reference_bond_ids=np.flatnonzero(alive==0)+1)
dump(O/'ELASTIC_PRECHECK.json',{'status':'PASS' if all(gates.values()) else 'DESIGN_NOT_ADMITTED','gates':gates,'error':err,'optimizer_message':None if opt is None else str(opt.message),'optimizer_success':None if opt is None else bool(opt.success),'evaluations':evals,'regions':regions,'U_LPS_J':u,'contact_U_J':uc,'max_residual_force_N':float(np.linalg.norm(res,axis=1).max()),'residual_L2_relative_to_external':float(np.linalg.norm(res)/np.linalg.norm(fext)),'mean_opening_m':float((x[rt['high'],2]-X[rt['high'],2]).mean()-(x[rt['low'],2]-X[rt['low'],2]).mean()),'external_force_sum_N':fext.sum(0),'reference_external_torque_N_m':np.cross(X,fext).sum(0),'current_external_torque_N_m':np.cross(x,fext).sum(0),'dead_load_potential_magnitude_J':float(np.sum(fext*(x-X))),'warning':'Single elastic endpoint. Optimizer iterations are not physical loading time; dead-load potential is not measured path work. Regional trend extrapolation is not a measured onset.','wall_s':time.perf_counter()-start,'CPU_user_s':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'CPU_system_s':resource.getrusage(resource.RUSAGE_SELF).ru_stime,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
print(json.dumps({'gates':gates,'regions':{k:v['max_utilization'] for k,v in regions.items()}},indent=2),flush=True);stream.close()
