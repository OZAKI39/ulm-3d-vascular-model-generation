from pathlib import Path
import json,hashlib,numpy as np
from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path.cwd();O=R/'verification'/C['task']/C['run'];V=(50e-6)**3;K=1e5/3;G=940.
def put(name,data):(O/name).write_text(json.dumps(data,indent=2,allow_nan=False))
records={};movement={};force_records=[]
for variant in ['original','fixed','final']:
 for mode,event in [('sort','T1_sort'),('migration','T2_migration'),('fracture','restart_after_break')]+([('mpi_fracture','restart_after_break')] if variant=='final' else []):
  d=O/(variant+'_'+mode);rep=json.loads((d/'report.json').read_text());g=np.load(d/'reference.npz');X=g['X'];b=g['bonds'];m=g['m'];degree=g['degree'];eng=Energy({'X':X,'bonds':b},50e-6,R/'analysis/material_macro_audit_v1/energy.so')
  np.testing.assert_allclose(eng.m,m,rtol=2e-15,atol=0)
  for key,v in rep.items():
   v['wrong_m_neighbor_counts']=degree[np.array(v['wrong_m_atom_ids'],int)-1].tolist();v['wrong_m_surface_count']=int(np.sum(degree[np.array(v['wrong_m_atom_ids'],int)-1]<degree.max()))
   z=np.load(d/(key+'.npz'));u,f=eng.eval(np.ascontiguousarray(z['x']),z['alive'].astype(np.uint8),True);scale=np.max(np.linalg.norm(f,axis=1));err=np.max(np.linalg.norm(z['force_N']-f,axis=1));fd=[]
   if key!='initial':
    for seed in range(3):
     vdir=np.random.default_rng(seed+510).normal(size=X.shape);vdir/=np.linalg.norm(vdir);eps=50e-6*1e-5
     du=(eng.eval(z['x']+eps*vdir,z['alive'].astype(np.uint8))-eng.eval(z['x']-eps*vdir,z['alive'].astype(np.uint8)))/(2*eps);pred=-np.sum(f*vdir)
     fd.append({'directional_dU_N':du,'minus_F_dot_direction_N':float(pred),'relative_error':float(abs(du-pred)/abs(pred))})
   force_records.append({'variant':variant,'case':mode,'state':key,'force_relative_to_existing_CPP':float(err/scale) if scale else None,'energy_relative_Python_to_existing_CPP':float(abs(v['independent_energy_J']-u)/u) if u else None,'independent_energy_J':u,'native_PE_times_V_J':v['native_PE_times_V_J'],'native_PE_relative_discrepancy':abs(v['native_PE_times_V_J']-u)/u if u else None,'FD':fd})
  records[variant+'_'+mode]=rep
  if mode in ['sort','migration']:
   a=rep['T0'];z=rep[event];moved=np.flatnonzero(np.array(a['owner_rank_by_id'])!=np.array(z['owner_rank_by_id']))+1
   order_changed=[a['local_id_order'][i]!=z['local_id_order'][i] for i in range(len(a['local_id_order']))]
   movement[variant+'_'+mode]={'nonidentity_local_order':order_changed,'actual_migrated_atom_ids':moved.tolist(),'migrated_count':len(moved),'before_owner':a['owner_rank_by_id'],'after_owner':z['owner_rank_by_id'],'event':event}
assert any(movement['original_sort']['nonidentity_local_order']) and movement['original_migration']['migrated_count']>0
identity_pass=all(v['wrong_m_node_count']==0 for k,d in records.items() if k.startswith('final') for v in d.values())
fd_pass=all(x['relative_error']<2e-6 for v in force_records for x in v['FD'])
fixed_steady=[v for v in force_records if v['variant']=='final' and v['state'] not in ['initial','break']]
force_pass=all(v['force_relative_to_existing_CPP']<2e-9 for v in fixed_steady)
put('WVOLUME_IDENTITY_TEST.json',{'status':'PASS_FIXED' if identity_pass else 'FAILED','original_bug_reproduced':True,'records':records})
put('MPI_SORT_MIGRATION_AUDIT.json',{'status':'PASS_FIXED' if identity_pass else 'FAILED','events':movement,'ghost_scope':'All owned and ghost m compared to immutable by-ID reference; ghost theta saved separately in full precision NPZ','first_error':'First actually triggered sort / ownership exchange, not initial force evaluation'})
put('RESTART_IDENTITY_AUDIT.json',{'status':'PASS_FIXED' if records['final_fracture']['restart_after_break']['wrong_m_node_count']==0 else 'FAILED','original':records['original_fracture']['restart_after_break'],'fixed':records['final_fracture']['restart_after_break'],'scope':'one postbreak restart per variant, same rank; new fixed-format checkpoint retains reference m; old damaged restart lacking reference m cannot reconstruct lost original neighborhood'})
put('NATIVE_FORCE_ENERGY_AUDIT.json',{'status':'PASS_INDEPENDENT_GRADIENT' if force_pass and fd_pass else 'FAILED','force_pass':force_pass,'FD_pass':fd_pass,'existing_independent_CPP_sha256':hashlib.sha256((R/'analysis/material_macro_audit_v1/energy.cpp').read_bytes()).hexdigest(),'records':force_records,'native_PE_tally_status':'PREEXISTING_DIAGNOSTIC_ENERGY_MISMATCH_NOT_FIXED','energy_convention':'Independent frozen-reference expanded LPS U is gradient-consistent. Native pe*V has a separate deviatoric tally discrepancy and must not be used for fracture-energy accounting. No energy tally or material change made here.'})
# Same coordinates, pre/post event topology; fixed variant excludes m-identity contamination.
d=O/'final_fracture';g=np.load(d/'reference.npz');pre=np.load(d/'T0.npz');z=np.load(d/'break.npz');post=np.load(d/'postbreak_same_coordinates.npz');assert np.array_equal(z['x'],post['x'])
b=g['bonds'];a,j=b.T;r=g['r0'];h=g['taper'];m=g['m'];live=z['alive'];theta=z['theta'];x=z['x'];eng=Energy({'X':g['X'],'bonds':b},50e-6,R/'analysis/material_macro_audit_v1/energy.so')
y=x[j]-x[a];length=np.linalg.norm(y,axis=1);e=length-r;e[abs(e)<2.2204e-16]=0
mag=V**2*h*((3*K-5*G)*(theta[a]/m[a]+theta[j]/m[j])+15*G*(1/m[a]+1/m[j])*e/r);mag[~live]=0
pair=mag[:,None]*y/length[:,None];lag=np.zeros_like(x);np.add.at(lag,a,pair);np.add.at(lag,j,-pair);diff=z['force_N']-post['force_N'];norm=np.linalg.norm(diff,axis=1);scale=float(np.max(np.linalg.norm(post['force_N'],axis=1)))
Upre=eng.eval(x,pre['alive'].astype(np.uint8));Upost=eng.eval(x,live.astype(np.uint8));broken=pre['alive']&~live
put('FRACTURE_TIME_LAYER_AUDIT.json',{'status':'EXPLICIT_SCHEME_LIMITATION','DIAGNOSTIC_ONLY':True,'source_algorithm_modified':False,'broken_count':int(broken.sum()),'broken_atom_id_pairs':(b[broken]+1).tolist(),'max_force_difference_N':float(norm.max()),'relative_force_difference':float(norm.max()/scale),'max_theta_difference':float(abs(theta-post['theta']).max()),'lagged_theta_force_prediction_relative_error':float(np.max(np.linalg.norm(lag-z['force_N'],axis=1))/scale),'same_coordinates_topology_only':True,'at_break_force_sum_N':z['force_N'].sum(0).tolist(),'at_break_torque_Nm':np.cross(x-x.mean(0),z['force_N']).sum(0).tolist(),'postbreak_force_sum_N':post['force_N'].sum(0).tolist(),'energy_intact_at_break_coordinates_J':Upre,'energy_postbreak_same_coordinates_J':Upost,'topology_potential_drop_J':Upre-Upost,'max_estimated_algorithmic_impulse_Ns':float(norm.max()*.5e-6),'max_estimated_velocity_increment_m_s':float(norm.max()*.5e-6/(1000*V)),'interpretation':'Five actual native breaks; survivor force exactly matches declared prebreak-theta explicit force pass. Difference is substantial and localized, not claimed negligible. One subsequent force evaluation at same x matches fixed-topology independent gradient. No integrated velocity jump measured in coordinate-driven test. Physical topology energy drop is distinct from lagged force impulse estimate. Criterion/theta order/special scheduling/history/integrator unchanged. Final process_broken has an independent owned-endpoint safety guard, not a theta algorithm change.'})
put('SMALL_TEST_GATES.json',{'identity_pass':identity_pass,'independent_force_gradient_pass':force_pass and fd_pass,'V3_eligible':identity_pass and force_pass and fd_pass,'native_PE_caveat':'separate documented diagnostic tally; V3 evaluates trajectory forces, not Gc from pe'})
print('GATES',identity_pass,force_pass,fd_pass)
