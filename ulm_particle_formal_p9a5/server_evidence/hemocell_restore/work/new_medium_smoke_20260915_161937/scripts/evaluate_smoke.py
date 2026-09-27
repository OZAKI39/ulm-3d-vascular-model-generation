#!/usr/bin/env python3
import pathlib,json,hashlib,sys,csv,math,traceback
import numpy as np
sys.dont_write_bytecode=True
R=pathlib.Path(__file__).resolve().parents[1]
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def table(p):
 a=np.genfromtxt(p,delimiter=',',names=True,dtype=float);return np.atleast_1d(a)
def snapshot(p):
 with p.open('rb') as f:
  n=int(np.fromfile(f,dtype='<u8',count=1)[0]);a=np.fromfile(f,dtype=[('id','<u8'),('rho','<f8'),('u','<f8',(3,))])
 assert len(a)==n and np.all(np.diff(a['id'].astype(np.int64))>0),'Snapshot size/order'
 assert np.all(np.isfinite(a['rho'])) and np.all(np.isfinite(a['u'])),'Snapshot nonfinite'
 return a
def evaluate(n):
 D=R/f'run_{n}';G=D/'diagnostics';c=read(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json');k=read(R/'contracts/SMOKE_SAFETY_CONTRACT.json');checks={};out=dict(horizon=n,checks=checks)
 def check(name,v):
  checks[name]='PASS' if bool(v) else 'FAIL'
 def near(a,b,tol):return abs(float(a)-float(b))<=tol
 try:
  t=read(D/'RUN_TERMINAL.json');s=read(G/'solver_status.json');a=read(G/'ACTUAL_NUMERICS_READBACK.json')
  check('process',t['status']=='PASS' and t['returncode']==0 and t['completed_steps']==n and t['binding_status']=='PASS')
  check('runtime_status',s['runtime_safety']=='PASS' and s['timesteps']==n and s['safety_checks']==n+1 and s['cell_count']==0 and s['mpi_ranks']==1 and s['mpi_runtime_correctness']=='PASS')
  check('static',read(R/'verification/NUMERICS_STATIC_CHECK.json')['status']=='PASS')
  check('contract_sha256',sha(R/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json')==t['numerics_contract_sha256']==a['numerics_contract_sha256'] and t['contract_unchanged_after_run'])
  check('actual_dt',a['dt_s']==c['dt_s'])
  check('actual_outlet_rho',all(a[f'outlet_{i:02}_rho_lu']==c[f'outlet_{i:02}_rho_lu'] for i in (1,2,3)))
  check('actual_units_and_inputs',a['pressure_unit_pa_read']==c['pressure_unit_pa'] and a['Qtarget']==c['physical_Qtarget_m3_s'] and a['multiplier']==c['inlet_numerical_multiplier'] and all(a[x]==c[x] for x in ['dx_m','tau','rho_kg_m3','nu_m2_s']) and not a['dt_or_outlet_density_recomputed_in_solver'])
  dt=c['dt_s'];dx=c['dx_m'];rho=c['rho_kg_m3'];qt=c['physical_Qtarget_m3_s'];mu=rho*dx**3
  check('independent_dt_formula',near(dt,(1/3)*(c['tau']-.5)*dx*dx/c['nu_m2_s'],8*math.ulp(dt)))
  check('independent_pressure_unit',near(c['pressure_unit_pa'],rho*(dx/dt)**2,8*math.ulp(c['pressure_unit_pa'])))
  pressures=[14.544978101274268,132.20454922317552,-13.700626673311461]
  check('independent_pressure_conversion',all(near(c[f'outlet_{i+1:02}_rho_lu'],1+3*p/c['pressure_unit_pa'],8*math.ulp(1.)) for i,p in enumerate(pressures)))
  sf=table(G/'safety_counts.csv');rh=table(G/'NEW_MEDIUM_RUNTIME_HISTORY.csv');fl=table(G/'NEW_MEDIUM_FLUX_HISTORY.csv');ma=table(G/'NEW_MEDIUM_MASS_HISTORY.csv')
  for name,x,spacing in [('safety',sf,1),('runtime',rh,50),('flux',fl,100),('mass',ma,500)]:
   check(name+'_cadence',np.array_equal(x['iteration'],np.arange(0,n+1,spacing)))
   check(name+'_finite',all(np.all(np.isfinite(x[col])) for col in x.dtype.names))
   if name!='safety':check(name+'_physical_time',np.all(np.abs(x['time_s']-x['iteration']*dt)<=8*np.spacing(np.maximum(x['time_s'],dt))))
  nonfinite=int(sum(np.sum(sf[x]) for x in sf.dtype.names if x!='iteration'));out['nonfinite_count']=nonfinite;check('nonfinite',nonfinite==0)
  rmin=s['all_step_fluid_rho_min'];rmax=s['all_step_fluid_rho_max'];mach=s['all_step_max_mach_including_ghosts']
  check('density',rmin>=k['rho_min'] and rmax<=k['rho_max'] and s['all_step_ghost_rho_min']>0)
  check('Mach_safety_margin',mach<=k['mach_smoke_pass_limit'] and max(np.max(rh['Mach_max']),np.max(rh['ghost_Mach_max']))<=mach)
  m0=rh['total_mass_kg'][0];cv0=ma['CV_mass_kg'][0]
  drift=(rh['total_mass_kg']-m0)/m0;massjump=float(np.max(np.abs(np.diff(rh['total_mass_kg']))/m0))
  check('mass_drift_identity',np.max(np.abs(drift-rh['total_mass_drift']))<1e-14)
  check('mass_stability',np.min(rh['total_mass_kg'])>0 and s['max_abs_mass_drift']<=k['max_total_mass_drift_fraction'] and massjump<=k['max_mass_jump_per50_fraction'])
  increments=np.diff(np.abs(ma['total_mass_drift']));limit=k['fast_monotonic_mass_growth_per500_fraction'];window=k['fast_monotonic_consecutive_intervals']
  rapid=any(np.all(increments[j:j+window]>limit) for j in range(len(increments)-window+1));check('no_fast_monotonic_mass_growth',not rapid)
  net=sum(fl[f'mass_outward_g{j}'] for j in (2,8,14,20));integ=[0.]
  for j in range(1,len(ma)):
   ix=(fl['iteration']>=ma['iteration'][j-1])&(fl['iteration']<=ma['iteration'][j]);tt=fl['time_s'][ix];yy=net[ix];integ.append(float(np.sum(.5*(yy[1:]+yy[:-1])*np.diff(tt))))
  integ=np.array(integ);delta=np.r_[0.,np.diff(ma['CV_mass_kg'])];err=delta+integ;cum=np.cumsum(err)
  check('CV_balance_identity',np.max(np.abs(integ-ma['integrated_net_outward_mass_kg']))/cv0<1e-12 and np.max(np.abs(err-ma['CV_balance_error_kg']))/cv0<1e-12)
  check('CV_balance',np.min(ma['CV_mass_kg'])>0 and np.max(np.abs(err))/cv0<=k['max_abs_CV_interval_balance_over_initial_mass'] and np.max(np.abs(cum))/cv0<=k['max_abs_CV_cumulative_balance_over_initial_mass'])
  qnames=[f'{p}_g{g}' for p in ('Qin','Qout01','Qout02','Qout03') for g in range(6)];q=np.stack([fl[x] for x in qnames],axis=1)
  check('flux_sanity',np.max(np.abs(q))/qt<=k['max_abs_flux_over_Qtarget'] and q[-1,2]/qt>=k['final_inlet_Q_min_over_Qtarget'])
  back={p:dict(count=int(np.sum(fl[p+'_g2']<0)),iterations=fl['iteration'][fl[p+'_g2']<0].astype(int).tolist()) for p in ('Qout01','Qout02','Qout03')}
  check('backflow_record_identity',bool(s['transient_outlet_backflow'])==any(v['count'] for v in back.values()))
  cvown=read(G/'control_volume_ownership_check.json');check('ownership',cvown['status']=='PASS' and cvown['CV_CELL_GLOBAL_COUNT']==180543 and cvown['physical_fluid_global_count']==182694 and s['quadrature_ownership_checks']>=n//100+1)
  ev=list(csv.DictReader((G/'snapshot_events.csv').open()));expected=k['full_field_snapshots'][str(n)]
  check('snapshot_schedule',[int(x['iteration']) for x in ev]==expected)
  samplesdir=G/'field_samples';files=list(samplesdir.glob('*.bin'));out['snapshot_files']=[str(p.relative_to(D)) for p in files]
  check('snapshot_files_schedule',all(int(p.stem.rsplit('_',1)[1]) in expected for p in files))
  # Full field and quadrature snapshots use the original Stage4 binary format.
  initial=snapshot(samplesdir/'fields_0.bin');check('initial_condition',len(initial)==182694 and np.all(initial['rho']==1) and np.all(initial['u']==0))
  check('initial_mass_identity',near(m0,len(initial)*mu,k['independent_snapshot_mass_relative_tolerance']*m0))
  initial_samples=snapshot(samplesdir/'samples_0.bin');check('initial_sample_velocity',np.all(initial_samples['u']==0))
  if n==5000:
   field=snapshot(samplesdir/'fields_5000.bin');sp=snapshot(samplesdir/'samples_5000.bin');cvIds=np.loadtxt(R/'contracts/control_volume_indices.txt',dtype=np.uint64);cvslots=np.searchsorted(field['id'],cvIds)
   check('snapshot_geometry_and_CV',np.array_equal(field['id'],initial['id']) and np.array_equal(field['id'][cvslots],cvIds))
   vr=float(np.max(np.linalg.norm(field['u'],axis=1)));recomputed=dict(rho_min=float(np.min(field['rho'])),rho_max=float(np.max(field['rho'])),rho_mean=float(np.mean(field['rho'])),Mach_max=vr/math.sqrt(1/3),total_mass_kg=float(np.sum(field['rho'])*mu),CV_mass_kg=float(np.sum(field['rho'][cvslots])*mu))
   check('snapshot_density_identity',all(near(recomputed[x],rh[x][-1],k['independent_snapshot_rho_abs_tolerance']) for x in ('rho_min','rho_max','rho_mean')))
   check('snapshot_velocity_identity',near(recomputed['Mach_max']*math.sqrt(1/3),rh['Mach_max'][-1]*math.sqrt(1/3),k['independent_snapshot_velocity_abs_tolerance']))
   check('snapshot_mass_identity',near(recomputed['total_mass_kg'],rh['total_mass_kg'][-1],k['independent_snapshot_mass_relative_tolerance']*m0) and near(recomputed['CV_mass_kg'],ma['CV_mass_kg'][-1],k['independent_snapshot_mass_relative_tolerance']*cv0))
   quad=np.loadtxt(R/'contracts/multiplane_quadrature.tsv');label=quad[:,0].astype(int);xyz=quad[:,1:4];base=np.floor(xyz).astype(np.int64);frac=xyz-base;vel=np.zeros((len(quad),3));rr=np.zeros(len(quad))
   header=(R/'contracts/solver_parameters.txt').read_text().splitlines()[2].split();nx,ny=int(header[0]),int(header[1])
   for ox in range(2):
    for oy in range(2):
     for oz in range(2):
      ix=base+np.array([ox,oy,oz]);ids=((ix[:,2]*ny)+ix[:,1])*nx+ix[:,0];slots=np.searchsorted(sp['id'],ids);assert np.array_equal(sp['id'][slots],ids),'Quadrature snapshot cells missing'
      w=np.prod(np.where(np.array([ox,oy,oz]),frac,1-frac),axis=1);vel+=sp['u'][slots]*w[:,None];rr+=sp['rho'][slots]*w
   dq=quad[:,4]*np.sum(vel*quad[:,5:8],axis=1)*dx**3/dt;fq=np.bincount(label,weights=dq,minlength=24);fm=np.bincount(label,weights=dq*rr*rho,minlength=24);fq[:6]*=-1
   qerr=float(np.max(np.abs(fq-q[-1]))/qt);merr=float(np.max(np.abs(fm-np.array([fl[f'mass_outward_g{j}'][-1] for j in range(24)])))/(qt*rho))
   check('snapshot_all24_flux_identity',qerr<=k['independent_snapshot_flux_over_Qtarget_tolerance'] and merr<=k['independent_snapshot_flux_over_Qtarget_tolerance']);recomputed.update(all24flux_max_error_over_Qtarget=qerr,all24massflux_max_error_over_rho_Qtarget=merr);out['independent_final_snapshot']=recomputed
  resources=list(csv.DictReader((D/'provenance/resource_trace.csv').open()));steady=[x for x in resources if n//10<=int(x['last_completed_step'])<n and int(x['solver_RSS_bytes'])>0];thirds=[]
  for j in range(3):
   lo=n//10+j*(n-n//10)/3;hi=n//10+(j+1)*(n-n//10)/3;z=[x for x in steady if lo<=int(x['last_completed_step'])<hi];thirds.append(dict(samples=len(z),vram_median_MiB=float(np.median([float(x['gpu_memory_used_MiB']) for x in z])) if z else None,RSS_median_MiB=float(np.median([float(x['solver_RSS_bytes'])/2**20 for x in z])) if z else None))
  sufficient=all(x['samples']>=2 for x in thirds);growth={}
  if sufficient:
   for key,lim in [('vram_median_MiB',k['VRAM_growth_limit_MiB']),('RSS_median_MiB',k['RSS_growth_limit_MiB'])]:
    v=[x[key] for x in thirds];growth[key]=dict(delta=v[-1]-v[0],limit=lim,growing=v[0]<=v[1]<=v[2] and v[-1]-v[0]>lim)
   check('GPU_memory_stability',not any(x['growing'] for x in growth.values()))
  elif n==5000:check('GPU_memory_stability',False)
  else:checks['GPU_memory_stability']='NOT_APPLICABLE_SHORT_WINDOW';out['quick_memory_limit']='5000-step gate supplies sustained-thirds evidence'
  out.update(rho_min=rmin,rho_max=rmax,max_mach=mach,total_mass_drift=float(drift[-1]),max_abs_mass_drift=s['max_abs_mass_drift'],max_mass_jump_per50_fraction=massjump,CV_interval_balance_max_fraction=float(np.max(np.abs(err))/cv0),CV_cumulative_balance_max_fraction=float(np.max(np.abs(cum))/cv0),final_flux_m3_s={p:float(fl[p+'_g2'][-1]) for p in ('Qin','Qout01','Qout02','Qout03')},max_abs_flux_over_Qtarget=float(np.max(np.abs(q))/qt),backflow=back,transient_backflow=any(v['count'] for v in back.values()),resource_thirds=thirds,memory_growth=growth,gpu_steps_per_second=t['solver_timing']['steps_per_second'],end_to_end_seconds=t['end_to_end_seconds'],initialization_seconds=t['initialization_seconds'],gpu_vram_peak_MiB=t['gpu_vram_peak_MiB'],gpu_utilization_timed_median=float(np.median([float(x['gpu_utilization_percent']) for x in steady])) if steady else None)
 except Exception as exc:
  checks['audit_execution']='FAIL';out['error']=str(exc);out['traceback']=traceback.format_exc()
 out['status']='FAIL' if 'FAIL' in checks.values() else 'PASS';out['SOLVER_FINALIZER_IDENTITY']=out['status']
 save(R/f'verification/SMOKE_{n}_AUDIT.json',out);save(R/('verification/QUICK_SMOKE_500_GATE.json' if n==500 else 'verification/SUSTAINED_SMOKE_5000_GATE.json'),dict(status=out['status'],source=f'SMOKE_{n}_AUDIT.json',failed_checks=[key for key,v in checks.items() if v=='FAIL']))
 print(json.dumps(out,indent=2));return out['status']=='PASS'
if __name__=='__main__':sys.exit(0 if evaluate(int(sys.argv[1])) else 2)
