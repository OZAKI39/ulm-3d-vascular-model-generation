"""Independent finalizer from raw response matrices and frozen contracts.
No import of solver code, no reuse of solver PASS labels, no fitted wall offset.
"""
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np,h5py
from analyze_matrices import C,R,Q,D,reference,compare,props,extract,groups,phase_statistics

def dump(name,data):
 (R/name).write_text(json.dumps(data,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def csvout(name,rows):
 if not rows:rows=[{'status':'NOT_RUN_PLANAR_GATE_FAILED'}]
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (R/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
paths=sorted((R/'remote_raw').glob('RAW_*.h5'))
for file in paths:
 mode,beta_text=file.stem.removeprefix('RAW_').split('_beta');receipt=R/'remote_raw'/f"{mode}_{beta_text.replace('.', 'p')}_RECEIPT.json"
 assert receipt.exists(),f'INCOMPLETE_STAGING {file}'
 record=json.loads(receipt.read_text())
 with h5py.File(file) as check:
  if record['exit_code']==0:assert check.attrs.get('terminal_state')=='COMPLETED_WITH_RECORDED_GATES'
rows,geos,leaks,perf,failures=extract(paths)
phase=phase_statistics(rows);anis=[];rough=[];contact=[];patch=[];raw_lookup={}
# One self-contained HDF5 archive retains remote raw operators and adds independent evidence.
out=R/'FIXED_MULTIBLOB_WALL_AUDIT.h5'
with h5py.File(out,'w') as dst:
 dst.attrs['contract_sha256']=hashlib.sha256((R/'contracts/FIXED_MULTIBLOB_WALL_FEASIBILITY_CONTRACT.json').read_bytes()).hexdigest();dst.attrs['wall_representation']='FIXED_HARD_SPHERES';dst.attrs['geometry_length_unit']='target radius a';dst.attrs['force_velocity_scaling']='work-conjugate sqrt(bulk drag); SI reconstructions saved for frozen d10/d50/d90';dst.create_dataset('matrix_action_test_vectors',data=Q);dst.create_dataset('dissipation_test_vectors',data=D)
 for fn in paths:
  root=dst.create_group(fn.stem)
  with h5py.File(fn) as src:
   for k,v in src.attrs.items():root.attrs[k]=v
   for gn,g in src.items():
    src.copy(g,root,name=gn);gg=root[gn]
    if g.attrs['category']=='PLANAR' and g.attrs.get('status')=='ASSEMBLED':
     x=g['coords'][:];a=g['radii'][:];b=float(g.attrs['beta']);meta=json.loads(g.attrs['lattice_metadata']);v1=np.array(meta['b1']);v2=np.array(meta['b2']);uv=np.array([(i/100,j/100) for i in range(101) for j in range(101)]);xy=uv[:,0,None]*v1+uv[:,1,None]*v2;dist2=np.sum((xy[:,None,:]-x[None,:,:2])**2,axis=2)
     surface=np.max(np.where(dist2<=a[None,:]**2,x[None,:,2]+np.sqrt(np.maximum(0,a[None,:]**2-dist2)),-np.inf),axis=1);valid=np.isfinite(surface);z=surface[valid]
     rr=dict(beta=b,layout=g.attrs['layout'],layers=int(g.attrs['layers']),surface_covered_fraction=float(valid.mean()),surface_peak_valley_over_a=float(np.ptp(z)),surface_RMS_about_mean_over_a=float(np.sqrt(np.mean((z-z.mean())**2))),surface_RMS_about_nominal_over_a=float(np.sqrt(np.mean(z*z))),uncovered_vertical_columns=int((~valid).sum()),geometry=gn);rough.append(rr)
     gg.create_dataset('roughness_xy',data=xy);gg.create_dataset('wall_surface_envelope_z',data=surface);gg.attrs['uncovered_surface_columns_are_not_filled']=True
    for qn,q in g.items():
     if not isinstance(q,h5py.Group) or 'R_scaled' not in q:continue
     A=q['R_scaled'][:];qq=gg[qn];m=A.shape[0];qq.create_dataset('test_generalized_velocities_scaled',data=np.eye(m));qq.create_dataset('test_generalized_force_torque_scaled',data=A);qq.create_dataset('inverse_residual_independent',data=np.linalg.norm(A@q['M_scaled'][:]-np.eye(m),ord=np.inf))
     if m!=6:
      ev12=np.linalg.eigvalsh((A+A.T)/2);rec12=float(np.linalg.norm(A-A.T)/np.linalg.norm(A));inv12=float(np.linalg.norm(A@q['M_scaled'][:]-np.eye(m),ord=np.inf));assert ev12.min()>0 and rec12<=1e-8 and inv12<=1e-10;qq.attrs['independent_min_eigenvalue']=float(ev12.min());qq.attrs['independent_reciprocity_error']=rec12
      continue
     B=reference(q.attrs['epsilon']);qq.create_dataset('RMBW_reference_total_scaled',data=B);qq.create_dataset('R_excess_scaled',data=A-np.eye(6));qq.create_dataset('RMBW_reference_excess_scaled',data=B-np.eye(6));qq.create_dataset('action_error_per_vector',data=np.linalg.norm(Q@(A-B).T,axis=1)/np.linalg.norm(Q@B.T,axis=1));qq.create_dataset('force_responses_256',data=Q@A.T);qq.create_dataset('power_10000',data=np.einsum('ni,ij,nj->n',D,A,D));si=qq.create_group('SI_scaled_from_certified_similarity')
     for label,ar in C['radii_m'].items():
      d=np.array([6*np.pi*C['mu_Pa_s']*ar]*3+[8*np.pi*C['mu_Pa_s']*ar**3]*3);si.create_dataset(label+'_R',data=A*np.sqrt(np.outer(d,d)))
     key=(float(g.attrs['beta']),str(g.attrs['layout']),int(g.attrs['layers']),float(g.attrs['spacing_ratio']),float(g.attrs['extent']),float(q.attrs['epsilon']),q.attrs.get('phase','NONE'));raw_lookup[key]=A
     if g.attrs['category']=='PLANAR':
      dirs=np.array([[1.,0,0],[2**-.5,2**-.5,0],[0,1.,0]]);tangent=np.einsum('ni,ij,nj->n',dirs,A[:3,:3],dirs);rot=np.einsum('ni,ij,nj->n',dirs,A[3:,3:],dirs);cross=np.cross([0,0,1],dirs);tr=np.einsum('ni,ij,nj->n',dirs,A[:3,3:],cross)
      ar=dict(beta=float(g.attrs['beta']),layout=g.attrs['layout'],layers=int(g.attrs['layers']),epsilon=float(q.attrs['epsilon']),phase=q.attrs['phase'],tangent_spread=float(np.ptp(tangent)/abs(tangent.mean())),RR_spread=float(np.ptp(rot)/abs(rot.mean())),TR_spread=float(np.ptp(tr)/max(abs(tr.mean()),1e-12)))
      for field,val in [('TT',tangent),('RR',rot),('TR',tr)]:
       for angle,v in zip([0,45,90],val):ar[f'{field}_{angle}']=float(v)
      anis.append(ar)
      if q.attrs['epsilon']==.01:
       xx=q['targets'][0,:2];dd=np.sum((g['coords'][:,:2]-xx)**2,axis=1);aw=g['radii'][:];z=np.where(dd<=(1+aw)**2,g['coords'][:,2]+np.sqrt(np.maximum(0,(1+aw)**2-dd)),-np.inf);contact.append(dict(beta=float(g.attrs['beta']),layout=g.attrs['layout'],layers=int(g.attrs['layers']),phase=q.attrs['phase'],x=float(xx[0]),y=float(xx[1]),first_contact_nominal_h_over_a=float(np.max(z)-1)))
 free=np.load(R/'raw/FREE_SPACE_FIRST.npz')['matrices'];fg=dst.create_group('free_space');fg.create_dataset('native_SI_resistance',data=free)
 free_rows=[];idx=0
 for size,a in C['radii_m'].items():
  for mu in [.0005,.001,.002]:
   d=np.array([6*np.pi*mu*a]*3+[8*np.pi*mu*a**3]*3);A=free[idx]/np.sqrt(np.outer(d,d));idx+=1;free_rows.append(dict(size=size,radius=a,mu=mu,TT_error=float(np.max(abs(A[:3,:3]-np.eye(3)))),RR_error=float(np.max(abs(A[3:,3:]-np.eye(3)))),TR_absolute=float(max(np.max(abs(A[:3,3:])),np.max(abs(A[3:,:3]))))))
 for key,A in raw_lookup.items():
  b,layout,layers,sp,extent,e,ph=key
  if extent not in [6,8,12] or ph not in ['BEAD','BRIDGE','PORE']:continue
  previous={6:4,8:6,12:8}[extent];B=raw_lookup.get((b,layout,layers,sp,previous,e,ph))
  if B is None:continue
  def rel(ids):return float(np.linalg.norm((A-B)[ids])/max(np.linalg.norm(A[ids]),1e-12))
  rr=dict(beta=b,layout=layout,layers=layers,spacing_ratio=sp,epsilon=e,phase=ph,patch_radius=extent,previous_radius=previous,TT_change=rel(np.s_[:3,:3]),RR_change=rel(np.s_[3:,3:]),TR_change=rel(np.s_[:3,3:]),normal_change=float(abs(A[2,2]-B[2,2])/abs(A[2,2])),tangent_change=float(max(abs(A[0,0]-B[0,0])/abs(A[0,0]),abs(A[1,1]-B[1,1])/abs(A[1,1]))));rr['pass']=rr['TT_change']<=.02 and rr['normal_change']<=.02 and rr['tangent_change']<=.02 and rr['RR_change']<=.05 and rr['TR_change']<=.05;patch.append(rr)
 for groupname,status in [('cylinder','NOT_RUN_PLANAR_GATE_FAILED'),('real_STL','NOT_RUN_PLANAR_GATE_FAILED'),('dynamic_smoke','NOT_RUN_PLANAR_GATE_FAILED')]:dst.create_group(groupname).attrs['status']=status
 cache_errors=[]
 for p in sorted((R/'raw').glob('NATIVE_CACHE_*_FIRST.npz')):
  z=np.load(p);cg=dst.create_group(p.stem)
  for k in z.files:cg.create_dataset(k,data=z[k])
  matrix_keys=[k for k in z.files if z[k].shape==(6,6)]
  assert len(matrix_keys)==2
  er=float(np.linalg.norm(z[matrix_keys[0]]-z[matrix_keys[1]])/np.linalg.norm(z[matrix_keys[0]]));cache_errors.append(er);assert er<=1e-10

wall_wall_checks=[]
for p in sorted((R/'raw').glob('WALL_WALL_LUBRICATION_S*.npz')):
 z=np.load(p);A=z['full_native_scaled'];T=z['target_native'];K=z['target_cached'];err=float(np.linalg.norm(T-K)/np.linalg.norm(T));ev=np.linalg.eigvalsh((A+A.T)/2);assert err<=1e-10 and ev.min()>0
 xx=z['coords'][1:];a=z['radii'][1:];dd=np.linalg.norm(xx[:,None]-xx[None,:],axis=2)-a[:,None]-a[None,:];np.fill_diagonal(dd,np.inf)
 wall_wall_checks.append(dict(file=p.name,full_R_condition=float(ev[-1]/ev[0]),wall_M_condition=float(np.linalg.cond(z['wall_M_scaled'])),target_R_condition=float(np.linalg.cond(T)),native_cached_error=err,min_wall_gap=float(dd.min())))
 with h5py.File(out,'a') as dst:
  g=dst.create_group(p.stem)
  for k in z.files:g.create_dataset(k,data=z[k])
csvout('WALL_WALL_LUBRICATION_DIAGNOSTIC.csv',wall_wall_checks)

remote_native_errors=[]
for p in sorted((R/'remote_raw').glob('REMOTE_NATIVE_CAPABILITY_*.npz')):
 z=np.load(p);er=float(np.linalg.norm(z['native']-z['cached'])/np.linalg.norm(z['native']));assert er<=1e-10;remote_native_errors.append(dict(file=p.name,error=er))
 with h5py.File(out,'a') as dst:
  g=dst.create_group(p.stem)
  for k in z.files:g.create_dataset(k,data=z[k])

core=[r for r in rows if r['category']=='PLANAR' and .01<=r['epsilon']<=.2];ultra=[r for r in rows if r['category']=='PLANAR' and .001<=r['epsilon']<=.01]
corephase=[r for r in phase if .01<=r['epsilon']<=.2];coreanis=[r for r in anis if .01<=r['epsilon']<=.2 and r['layout']=='HEX']
resolutions=[]
for b in sorted(set(r['beta'] for r in core),reverse=True):
 rr=[r for r in core if r['beta']==b];pr=[r for r in corephase if r['beta']==b and r['mode'] in ['normal','tangent']];resolutions.append(dict(beta=b,core_queries=len(rr),worst_action=max(r['matrix_action_error'] for r in rr),mean_action=float(np.mean([r['matrix_action_error'] for r in rr])),worst_TT_phase=max(r['spread'] for r in pr)))
matched=[]
bs=sorted({r['beta'] for r in core},reverse=True)
for before,after in zip(bs[:-1],bs[1:]):
 key=lambda r:(r['layout'],r['layers'],r['epsilon'],r['phase'])
 aa={key(r):r for r in core if r['beta']==before};bb={key(r):r for r in core if r['beta']==after};keys=sorted(set(aa)&set(bb));va=np.array([aa[k]['matrix_action_error'] for k in keys]);vb=np.array([bb[k]['matrix_action_error'] for k in keys]);matched.append(dict(beta_before=before,beta_after=after,matched_query_count=len(keys),mean_action_before=float(va.mean()),mean_action_after=float(vb.mean()),mean_absolute_improvement=float((va-vb).mean()),worst_action_before=float(va.max()),worst_action_after=float(vb.max()),unmatched_configurations_excluded=True))
csvout('MATCHED_REFINEMENT.csv',matched)

stop=False
if len(resolutions)>=2:
 b0=next(r for r in resolutions if r['beta']==.5);b1=next(r for r in resolutions if r['beta']==.25);improvement=b0['mean_action']-b1['mean_action'];stop=b1['worst_action']>.5 and b1['worst_TT_phase']>.1 and improvement<.05
else:improvement=None
if any(r['min_eigenvalue']<=0 for r in rows):stop=True
# No selecting a passing phase or averaging away bad phase. Each configuration must pass every required core point.
configs=[]
for key,rr in groups(core,['beta','layout','layers']).items():
 b,la,ly=key;ph=[r for r in corephase if (r['beta'],r['layout'],r['layers'])==key];an=[r for r in anis if (r['beta'],r['layout'],r['layers'])==key and .01<=r['epsilon']<=.2]
 vals=dict(max_normal=max(r['TT_normal_error'] for r in rr),max_tangent=max(r['TT_parallel_error'] for r in rr),max_action=max(r['matrix_action_error'] for r in rr),max_RR=max(r['RR_error'] for r in rr),max_TR=max(r['TR_error'] for r in rr));acc=vals['max_normal']<=.05 and vals['max_tangent']<=.05 and vals['max_action']<=.05 and vals['max_RR']<=.1 and vals['max_TR']<=.1
 phasepass=all(r['spread']<=(.05 if r['mode'] in ['normal','tangent'] else .1) for r in ph);anpass=all(r['tangent_spread']<=.05 and r['RR_spread']<=.1 for r in an)
 pp=[r for r in patch if (r['beta'],r['layout'],r['layers'])==key and r['patch_radius']==12 and .01<=r['epsilon']<=.2];patchpass=len(pp)==15 and all(r['pass'] for r in pp)
 ll=[r for r in leaks if r['beta']==b and r['layout']==la and r['layers']==ly and r['spacing_ratio']==1.01];one=next((r for r in leaks if r['beta']==b and r['layout']==la and r['layers']==1 and r['spacing_ratio']==1.01),None);three=next((r for r in leaks if r['beta']==b and r['layout']==la and r['layers']==3 and r['spacing_ratio']==1.01),None)
 leakpass=bool(ll) and all(abs(r['transmission'])<=.05 for r in ll) and one is not None and three is not None and abs(three['transmission'])/max(abs(one['transmission']),1e-12)<=.8
 configs.append(dict(beta=b,layout=la,layers=ly,**vals,accuracy_pass=acc,phase_pass=phasepass,anisotropy_pass=anpass,patch_pass=patchpass,leakage_pass=leakpass,configuration_pass=len(rr)==95 and acc and phasepass and anpass and patchpass and leakpass,core_query_count=len(rr),core_grid_complete=len(rr)==95))
planar=any(r['configuration_pass'] for r in configs)
if planar:raise RuntimeError('PLANAR_PASS_REQUIRES_AUTHORIZED_CONDITIONAL_CURVED_WORK_BEFORE_FINAL_REPORT')
summary=dict(status='FAIL_WALL_REPRESENTATION',scope='tested fixed hard-sphere Pecnut implementation; no universal exclusion of regularized blobs or finer untested limits',independent_finalizer='PASS',native_cache_errors=cache_errors,wall_wall_checks=wall_wall_checks,remote_native_errors=remote_native_errors,resolutions=resolutions,matched_refinement=matched,coarse_to_fine_stop=stop,mean_action_improvement=improvement,configurations=configs,selected_configuration=None,core_max_normal=max(r['TT_normal_error'] for r in core),core_max_tangent=max(r['TT_parallel_error'] for r in core),core_max_action=max(r['matrix_action_error'] for r in core),core_max_RR=max(r['RR_error'] for r in core),core_max_TR=max(r['TR_error'] for r in core),ultra_max_action=max(r['matrix_action_error'] for r in ultra),max_core_TT_phase=max(r['spread'] for r in corephase if r['mode'] in ['normal','tangent']),max_core_TR_phase=max(r['spread'] for r in corephase if r['mode']=='TR'),max_core_HEX_TT_anisotropy=max(r['tangent_spread'] for r in coreanis),max_core_HEX_RR_anisotropy=max(r['RR_spread'] for r in coreanis),max_reciprocity=max(r['reciprocity'] for r in rows),min_eigenvalue=min(r['min_eigenvalue'] for r in rows),max_condition=max(r['condition'] for r in rows),minimum_dissipation=min(r['minimum_dissipation'] for r in rows),max_linear_residual=max(r['linear_residual'] for r in rows),completed_response_queries=len(rows),failed_queries=len(failures),max_wall_count_measured=max(r['N_wall'] for r in rows),max_leakage=max((abs(r['transmission']) for r in leaks),default=None),min_leakage=min((abs(r['transmission']) for r in leaks),default=None),max_surface_roughness_RMS=max(r['surface_RMS_about_mean_over_a'] for r in rough),free_space=free_rows,free_space_pass=all(r['TT_error']<=1e-4 and r['RR_error']<=1e-4 and r['TR_absolute']<=1e-8 for r in free_rows),patch_largest_pass_count=sum(r['pass'] for r in patch if r['patch_radius']==12),patch_largest_count=sum(r['patch_radius']==12 for r in patch),performance=perf,curved='NOT_RUN_PLANAR_GATE_FAILED',real_STL='NOT_RUN_PLANAR_GATE_FAILED',online='UNVERIFIED_NO_ACCEPTED_CONFIGURATION_OR_PARTICLE_TIMESTEP_BUDGET',offline='UNVERIFIED_CURVED_NOT_RUN',human_visual_review='PENDING')
assert summary['free_space_pass'] and summary['max_reciprocity']<=1e-8 and summary['min_eigenvalue']>0 and summary['minimum_dissipation']>0 and summary['max_linear_residual']<=1e-10
# Recompute leakage layer trend from actual response, not labels.
for r in leaks:
 one=next((a for a in leaks if a['beta']==r['beta'] and a['layout']==r['layout'] and a['spacing_ratio']==r['spacing_ratio'] and a['layers']==1),None);r['relative_to_one_layer']=r['transmission']/one['transmission'] if one else None;r['magnitude_relative_to_one_layer']=abs(r['transmission'])/max(abs(one['transmission']),1e-12) if one else None
layer_checks=[dict(beta=r['beta'],layout=r['layout'],spacing_ratio=r['spacing_ratio'],layer3_to_layer1_magnitude_ratio=r['magnitude_relative_to_one_layer'],pass_gate=r['magnitude_relative_to_one_layer']<=.8) for r in leaks if r['layers']==3]
summary['leakage_absolute_gate_pass']=all(abs(r['transmission'])<=.05 for r in leaks);summary['leakage_layer_checks']=layer_checks;summary['leakage_layer_gate_pass']=all(r['pass_gate'] for r in layer_checks) and bool(layer_checks)
csvout('LEAKAGE_LAYER_TREND.csv',layer_checks)
csvout('PLANAR_CONVERGENCE.csv',rows);csvout('LATERAL_PHASE_SENSITIVITY.csv',phase);csvout('ANISOTROPY.csv',anis);csvout('LEAKAGE_SCREENING.csv',leaks);csvout('MATRIX_PROPERTIES.csv',[{k:r[k] for k in ['raw_file','geometry','query','reciprocity','min_eigenvalue','condition','minimum_dissipation','linear_residual']} for r in rows]);csvout('CYLINDER_CURVATURE.csv',[]);csvout('REAL_STL_PATCH_AUDIT.csv',[]);csvout('PERFORMANCE_SCALING.csv',perf);csvout('PATCH_CONVERGENCE.csv',patch);csvout('WALL_ROUGHNESS.csv',rough);csvout('FIRST_CONTACT.csv',contact);csvout('EXECUTION_FAILURES.csv',failures or [{'status':'NO_FAILED_QUERIES'}]);csvout('GEOMETRY_EXECUTION.csv',geos);dump('validation/INDEPENDENT_FINALIZER.json',summary)
print(json.dumps({k:v for k,v in summary.items() if k not in ['configurations','free_space','performance']},indent=2))
