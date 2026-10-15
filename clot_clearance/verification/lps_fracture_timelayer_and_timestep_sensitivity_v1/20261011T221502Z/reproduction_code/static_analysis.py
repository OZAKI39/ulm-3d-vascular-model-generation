from common import *
from model import LPS
import numpy as np,csv,sys,time
sys.path.insert(0,str(R))
from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy
t0=time.perf_counter();d=O/'historical_static';g=np.load(d/'reference.npz');model=LPS(g)
pre=np.load(d/'T0.npz');br=np.load(d/'break.npz');post=np.load(d/'postbreak_same_coordinates.npz')
x=br['x'];before=pre['alive'];after=br['alive'];assert np.array_equal(x,post['x'])
thpre,fpre,upre=model.eval(x,before);thpost,fpost,upost=model.eval(x,after);_,mixed,_=model.eval(x,after,thpre)
df_top=fpost-fpre;df_lag=mixed-fpost;broken=before&~after;scale=float(np.linalg.norm(fpost,axis=1).max())
eng=Energy({'X':g['X'],'bonds':g['bonds']},50e-6,R/'analysis/material_macro_audit_v1/energy.so')
u_cpp,f_cpp=eng.eval(np.ascontiguousarray(x),after.astype(np.uint8),True)
fd=[]
for seed in range(3):
 v=np.random.default_rng(seed+520).normal(size=x.shape);v/=np.linalg.norm(v);eps=5e-10
 actual=(eng.eval(x+eps*v,after.astype(np.uint8))-eng.eval(x-eps*v,after.astype(np.uint8)))/(2*eps);pred=-float(np.sum(fpost*v))
 fd.append(dict(derivative_N=actual,negative_force_dot_direction_N=pred,absolute_N=abs(actual-pred),relative=abs(actual-pred)/abs(pred)))
err=float(np.linalg.norm(mixed-br['force_N'],axis=1).max());err_post=float(np.linalg.norm(fpost-f_cpp,axis=1).max())
assert err/scale<2e-9 and err_post/scale<2e-9 and max(q['relative'] for q in fd)<2e-6
assert np.max(abs(br['m']-g['m'])/g['m'])<2e-12
event_nodes=np.unique(g['bonds'][broken]);adj=np.zeros(len(x),bool);adj[event_nodes]=True
np.savez_compressed(O/'HISTORICAL_SAME_COORDINATE_DECOMPOSITION.npz',id=g['id'],X=g['X'],x=x,bonds=g['bonds'],alive_pre=before,alive_post=after,broken_bond_id=np.flatnonzero(broken)+1,broken_atom_id_pairs=g['bonds'][broken]+1,theta_pre=thpre,theta_post=thpost,F_pre_N=fpre,F_mixed_N=mixed,F_post_N=fpost,DeltaF_topology_N=df_top,DeltaF_timelayer_N=df_lag,m=g['m'],native_F_mixed_N=br['force_N'],event_related_node=adj)
fields=['atom_id','reference_x_m','reference_y_m','reference_z_m','x_m','y_m','z_m','event_endpoint','region','theta_pre','theta_post']+[f'{name}_{axis}_N' for name in ['F_pre','F_mixed','F_post','DeltaF_topology','DeltaF_timelayer'] for axis in 'xyz']
with (O/'HISTORICAL_SAME_COORDINATE_NODES.csv').open('w') as f:
 w=csv.writer(f);w.writerow(fields)
 for i in range(len(x)):w.writerow([i+1,*g['X'][i],*x[i],int(adj[i]),'loaded_node' if i==89 else ('break_endpoint' if adj[i] else 'other'),thpre[i],thpost[i],*fpre[i],*mixed[i],*fpost[i],*df_top[i],*df_lag[i]])
rec=dict(status='PASS_STATIC_DECOMPOSITION',source=str(PRIOR/'final_fracture'),native_calculations_added=0,physical_time_s=None,physical_time_reason='Historical coordinate-driven test has no integrator; step labels are not a loading history',broken_bond_ids=(np.flatnonzero(broken)+1).tolist(),broken_atom_id_pairs=(g['bonds'][broken]+1).tolist(),force_reference_scale_N=scale,max_DeltaF_topology_N=float(np.linalg.norm(df_top,axis=1).max()),max_DeltaF_timelayer_N=float(np.linalg.norm(df_lag,axis=1).max()),max_theta_change=float(abs(thpre-thpost).max()),U_pre_J=upre,U_post_J=upost,DeltaU_topology_J=upost-upre,F_mixed_vs_native_max_abs_N=err,F_mixed_vs_native_relative=err/scale,F_post_vs_independent_CPP_max_abs_N=err_post,F_post_vs_independent_CPP_relative=err_post/scale,independent_CPP_energy_J=u_cpp,finite_difference_checks=fd,m_relative_error_by_actual_ID=float(np.max(abs(br['m']-g['m'])/g['m'])),interpretation='Signed topology energy change is not calibrated irreversible fracture energy. Historical DeltaF*dt is only a dimensional proxy: no historical integrated impulse exists.',wall_s=time.perf_counter()-t0)
dump(O/'HISTORICAL_STATIC_DECOMPOSITION_AUDIT.json',rec);print(json.dumps(rec))
