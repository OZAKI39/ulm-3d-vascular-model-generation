"""Independent audit of frozen field, original inputs, raw trajectories and dumps.
Does not call the C++ sampler, adapter, or a LAMMPS solver.
"""
from pathlib import Path
import argparse,json,hashlib,math,re
import numpy as np,h5py
from reference_flow_sampler import ReferenceField
from reference_drag import drag,tau,uniform_trajectory,linear_trajectory

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def maximum(a):return float(np.max(np.abs(a)))
def table(p):return np.atleast_1d(np.genfromtxt(p,delimiter=',',names=True))
def cols(t,names):return np.column_stack([t[k] for k in names])
def dump(p):
 with Path(p).open() as f:
  assert f.readline().strip()=='ITEM: TIMESTEP';step=int(f.readline());assert f.readline().strip()=='ITEM: NUMBER OF ATOMS';n=int(f.readline());assert f.readline().startswith('ITEM: BOX BOUNDS')
  box=np.array([[float(x) for x in f.readline().split()] for _ in range(3)]);names=f.readline().split()[2:];a=np.atleast_2d(np.loadtxt(f))
 assert len(a)==n;a=a[np.argsort(a[:,names.index('id')])]
 return step,{name:a[:,i] for i,name in enumerate(names)}
def integrity(root):
 for line in (root/'FROZEN_INPUT_SHA256SUMS').read_text().splitlines():
  expected,name=line.split('  ',1)
  if sha(root/name)!=expected:return False
 return True

def unit_audit(root):
 root=Path(root);ct=read(root/'contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json');g=ct['gates'];syn=read(root/'contracts/SYNTHETIC_FIELDS_CONTRACT.json');checks={};metrics={};out=root/'validation/cpp';q=root/'validation/queries'
 def check(k,b):checks[k]=bool(b)
 check('frozen_identity',integrity(root));fc=read(root/'contracts/FROZEN_FLOW_FIELD_CONTRACT.json');check('field_sha',sha(root/fc['field_file'])==fc['field_sha256']==ct['field_sha256'])
 sources=[Path(fc['source_file']),Path('/workspace/hemocell_restore/results/new_medium_smoke_20260915_161937/run_5000/diagnostics/field_samples/fields_5000.bin')]
 source=next((p for p in sources if p.exists()),None)
 check('original_snapshot_available',source is not None)
 if source is not None:
  check('original_snapshot_sha',sha(source)==fc['source_sha256'])
  with source.open('rb') as f:
   n=int(np.fromfile(f,dtype='<u8',count=1)[0]);raw=np.fromfile(f,dtype=[('id','<u8'),('rho','<f8'),('u','<f8',(3,))])
  real=ReferenceField(root/fc['field_file']);num=read(root/'provenance/PBS_BSA_NUMERICS.json')
  check('original_vector_export',n==len(raw)==len(real.ids) and np.array_equal(raw['id'],real.ids) and np.array_equal(real.ids,real.fluid) and np.array_equal(raw['u']*(num['dx_m']/num['dt_s']),real.u))
  check('geometry_identity',sha(root/'provenance/closed_geometry_m.stl')==fc['geometry_sha256'])
 mapping=table(root/'validation/COORDINATE_NODE_CHECK.csv');xyz=cols(mapping,['mapped_x','mapped_y','mapped_z']);original=cols(mapping,['vtk_x','vtk_y','vtk_z']);derived=np.array(fc['origin_m'])+cols(mapping,['ix','iy','iz'])*fc['dx_m'];metrics['max_coordinate_error_m']=max(maximum(xyz-original),maximum(derived-original));check('coordinate_mapping',metrics['max_coordinate_error_m']<=g['coordinate_m'] and len(mapping)>=100)
 for label,file,kind in [('uniform','UNIFORM_FIELD.h5','uniform'),('linear','LINEAR_FIELD.h5','linear'),('real_nodes','FROZEN_FLOW_FIELD_V0.h5','real_nodes'),('real_interior','FROZEN_FLOW_FIELD_V0.h5','real')]:
  inp=table(q/(label+'.csv'));r=table(out/(label+'.csv'));positions=cols(inp,['x','y','z']);cpp=cols(r,['ux','uy','uz']);field=ReferenceField(root/'fields'/file);py,status=field.sample(positions)
  check(label+'_valid',np.all(r['status']==0) and np.all(status==0) and np.array_equal(inp['id'],r['id']) and np.all(np.isfinite(cpp)))
  metrics[label+'_cpp_python_m_s']=maximum(cpp-py);check(label+'_cpp_python',metrics[label+'_cpp_python_m_s']<=g['cpp_python_sampler_m_s'])
  if kind=='uniform':ref=np.tile(syn['uniform_velocity_m_s'],(len(inp),1));tol=g['uniform_sampler_m_s']
  elif kind=='linear':ref=np.array(syn['linear_intercept_m_s'])+positions@np.array(syn['linear_matrix_per_s']).T;tol=g['linear_sampler_m_s']
  elif kind=='real_nodes':ids=np.load(q/'real_node_ids.npy');loc=np.searchsorted(field.ids,ids);assert np.array_equal(field.ids[loc],ids);ref=field.u[loc];tol=g['real_node_identity_m_s']
  else:ref=py;tol=g['cpp_python_sampler_m_s']
  metrics[label+'_reference_error_m_s']=maximum(cpp-ref);check(label+'_reference',metrics[label+'_reference_error_m_s']<=tol)
  check(label+'_sample_count',len(inp)>=(1000 if label=='real_nodes' else 10000))
 for label,expected in [('invalid_domain',[1,1,1,4]),('invalid_missing',[3]),('invalid_solid',[2])]:
  r=table(out/(label+'.csv'));check(label,np.array_equal(r['status'],expected))
 fin=table(q/'force.csv');fout=table(out/'force.csv');reference=drag(.001,fin['d'],cols(fin,['ufx','ufy','ufz']),cols(fin,['vx','vy','vz']));computed=cols(fout,['fx','fy','fz']);norm=np.linalg.norm(reference,axis=1);nz=norm>0
 metrics['max_force_relative_error']=float(np.max(np.linalg.norm(computed[nz]-reference[nz],axis=1)/norm[nz]));metrics['max_zero_force_absolute_N']=maximum(computed[~nz]);check('force',metrics['max_force_relative_error']<=g['force_relative'] and metrics['max_zero_force_absolute_N']<=g['force_zero_absolute_N'] and len(fin)>=10000)
 return {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'metrics':metrics}

def case_audit(root,name):
 root=Path(root);ct=read(root/'contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json');g=ct['gates'];sp=next(s for s in ct['cases'] if s['name']==name);p=root/'cases'/name;checks={};m={}
 def check(k,b):checks[k]=bool(b)
 try:
  check('frozen_identity',integrity(root));check('case_field_sha',sha(root/sp['field'])==sp['field_sha256'] and sp['canonical_real_field_sha256']==ct['field_sha256'])
  rows=np.concatenate([table(f) for f in sorted(p.glob('trajectory_rank*.csv'))]);rows.sort(order=['step','id']);check('trace_not_empty',len(rows)>0)
  m['nonfinite_count']=sum(int(np.sum(~np.isfinite(rows[n]))) for n in rows.dtype.names);check('trace_finite',m['nonfinite_count']==0);m['invalid_query_count']=int(np.sum(rows['query_status']!=0));check('all_queries_valid',m['invalid_query_count']==0)
  first,di=dump(p/'initial.dump');last,df=dump(p/'final.dump');check('steps',first==0 and last==sp['steps']);ids=np.arange(1,sp['N']+1)
  check('ids',np.array_equal(di['id'],ids) and np.array_equal(df['id'],ids));check('fixed_diameter',np.array_equal(di['diameter'],df['diameter']))
  pop=table(root/'populations'/(sp['population']+'.csv'));src_ids=np.array(sp['source_bubble_ids']);diam=pop['diameter_um'][src_ids]*1e-6;check('diameter_source',np.array_equal(pop['bubble_id'][src_ids],src_ids) and maximum(df['diameter']/1e-6-pop['diameter_um'][src_ids])<=g['diameter_roundtrip_um'] and np.array_equal(diam,np.array(sp['diameters_m'])))
  meta=read(root/'populations'/(sp['population']+'.metadata.json'));check('population_sha',sha(root/'populations'/(sp['population']+'.csv'))==meta['population_sha256'])
  check('mass',maximum((df['mass']-1000*np.pi*diam**3/6)/(1000*np.pi*diam**3/6))<1e-12)
  counts=np.unique(rows['step'],return_counts=True);check('every_record_ids',np.all(counts[1]==sp['N']) and len(set(zip(rows['step'],rows['id'])))==len(rows))
  expected=np.arange(0,sp['steps']+1,sp['trajectory_stride']);expected=np.unique(np.r_[expected,sp['steps']]);check('cadence',np.array_equal(counts[0],expected))
  check('physical_time',maximum(rows['time_s']-rows['step']*sp['dt_s'])<1e-18)
  position=cols(rows,['x','y','z']);velocity=cols(rows,['vx','vy','vz']);field=ReferenceField(root/sp['field']);uf,status=field.sample(position);check('independent_query_valid',np.all(status==0));m['max_cpp_python_sampled_velocity_error_m_s']=maximum(cols(rows,['ufx','ufy','ufz'])-uf);check('sampler',m['max_cpp_python_sampled_velocity_error_m_s']<=g['cpp_python_sampler_m_s'])
  ref=drag(.001,rows['diameter_m'],uf,velocity);actual=cols(rows,['Fx','Fy','Fz']);applied_ref=drag(.001,rows['diameter_m'],uf,cols(rows,['force_vx','force_vy','force_vz']));applied=cols(rows,['applied_Fx','applied_Fy','applied_Fz'])
  # A field sample can contain roundoff. Normalize force differences by the drag scale,
  # independently check the pure function at 1e-12, and retain exact velocity phases.
  denom=3*np.pi*.001*rows['diameter_m']*np.maximum(np.linalg.norm(uf,axis=1),1e-12)
  m['max_force_scale_error']=float(np.max(np.linalg.norm(actual-ref,axis=1)/denom));m['max_applied_force_scale_error']=float(np.max(np.linalg.norm(applied-applied_ref,axis=1)/denom));check('stokes_force',max(m['max_force_scale_error'],m['max_applied_force_scale_error'])<=g['force_relative'])
  finrows=rows[rows['step']==sp['steps']];check('dump_vs_trace',maximum(cols(finrows,['x','y','z'])-cols(df,['xu','yu','zu']))<=1e-18 and maximum(cols(finrows,['vx','vy','vz'])-cols(df,['vx','vy','vz']))<=1e-18);check('force_reaches_LAMMPS_array',maximum(cols(finrows,['applied_Fx','applied_Fy','applied_Fz'])-cols(df,['fx','fy','fz']))<=g['force_zero_absolute_N'])
  log=(p/'log.lammps').read_text();runtime=read(p/'RUN_METRICS.json');check('runtime',runtime['returncode']==0 and 'COUPLING_CASE_COMPLETED' in log and 'ERROR' not in log);check('log_finite',re.search(r'\b(?:nan|inf)\b',log,re.I) is None)
  m['lost_atoms']=sp['N']-len(df['id']);check('lost_atoms',m['lost_atoms']==0 and 'Lost atoms:' not in log)
  counters=[read(f) for f in sorted(p.glob('trajectory_counters_rank*.json'))];check('every_step_counters',len(counters)==sp['mpi_ranks'] and sum(x['force_queries'] for x in counters)==sp['N']*(sp['steps']+2) and all(x['invalid_query_count']==x['nonfinite_count']==0 and x['completed_step']==sp['steps'] for x in counters))
  m['force_queries']=sum(x['force_queries'] for x in counters);m['wall_seconds']=runtime['wall_seconds'];m['peak_rss_kib']=runtime['process_tree_peak_rss_kib']
  if sp['N']==1:
   x0=cols(di,['xu','yu','zu'])[0];v0=cols(di,['vx','vy','vz'])[0];d=float(diam[0]);tp=float(tau(d));m['diameter_um']=d/1e-6;m['tau_p_s']=tp
   if sp['field_kind'] in ['uniform','linear']:
    syn=read(root/'contracts/SYNTHETIC_FIELDS_CONTRACT.json');u0=np.array(syn['uniform_velocity_m_s']);vscale=float(np.linalg.norm(u0));xscale=vscale*tp
    if sp['field_kind']=='uniform':xref,vref=uniform_trajectory(rows['time_s'],x0,v0,u0,d)
    else:
     u0=np.array(syn['linear_intercept_m_s'])+np.array(syn['linear_matrix_per_s'])@x0;vscale=float(np.linalg.norm(u0));xscale=vscale*tp;xref,vref,odemeta=linear_trajectory(rows['time_s'],x0,v0,d,syn['linear_intercept_m_s'],syn['linear_matrix_per_s']);m['ODE_reference']=odemeta
     analytic=np.array(syn['linear_intercept_m_s'])+position@np.array(syn['linear_matrix_per_s']).T;m['sampled_vs_affine_error_m_s']=maximum(cols(rows,['ufx','ufy','ufz'])-analytic);check('linear_sampler_during_trajectory',m['sampled_vs_affine_error_m_s']<=g['linear_sampler_m_s'])
    m['max_velocity_normalized_error']=float(np.max(np.linalg.norm(velocity-vref,axis=1))/vscale);m['max_position_normalized_error']=float(np.max(np.linalg.norm(position-xref,axis=1))/xscale);label='analytic' if sp['field_kind']=='uniform' else 'linear'
    check('reference_trajectory',m['max_velocity_normalized_error']<=g[label+'_velocity_normalized'] and m['max_position_normalized_error']<=g[label+'_position_normalized'])
    if sp['field_kind']=='uniform':
     slip=np.linalg.norm(velocity-u0,axis=1)/vscale;m['final_normalized_slip']=float(slip[-1]);check('relaxation',slip[-1]<=g['final_slip_over_initial_max'] and rows['time_s'][-1]>=5*tp)
     # Interpolate the actual decay crossing at e^-1, no fitted relaxation parameter.
     cross=np.flatnonzero(slip<=math.exp(-1))[0];j=max(cross-1,0);m['measured_e_folding_time_s']=float(np.interp(math.exp(-1),slip[[cross,j]],rows['time_s'][[cross,j]]))
   else:
    safe=sp['safety'];dist=np.linalg.norm(position-np.array(safe['spawn_m']),axis=1);lower=safe['center_wall_distance_m']-dist;m['min_center_wall_lower_bound_m']=float(lower.min());m['min_surface_clearance_lower_bound_m']=float(lower.min()-d/2);m['max_displacement_m']=float(dist.max());m['steps']=sp['steps'];m['physical_time_s']=sp['steps']*sp['dt_s'];check('safe_region',np.all(lower>d/2+3*field.dx) and maximum(lower-rows['center_wall_lower_bound_m'])<=1e-18)
    check('real_field_drives_particle',np.linalg.norm(position[-1]-x0)>0 and np.linalg.norm(velocity[-1])>0)
  if sp['N']>1:
   from scipy.spatial.distance import pdist
   initial=cols(di,['xu','yu','zu']);delta=np.abs(initial[:,None,:]-initial[None,:,:]);box=80e-6;delta=np.minimum(delta,box-delta);pairdist=np.linalg.norm(delta,axis=2);np.fill_diagonal(pairdist,np.inf)
   mincenter=float(pairdist.min());maxtravel=float(np.linalg.norm(read(root/'contracts/SYNTHETIC_FIELDS_CONTRACT.json')['uniform_velocity_m_s']))*sp['steps']*sp['dt_s'];m['noncontact_surface_gap_lower_bound_m']=mincenter-float(diam.max())-2*maxtravel
   check('nonoverlap_entire_short_horizon',m['noncontact_surface_gap_lower_bound_m']>0)
  if name=='case_e_mpi4':
   _,ref=dump(root/'cases/case_e_mpi1/final.dump');m['mpi_max_position_difference_m']=maximum(cols(df,['xu','yu','zu'])-cols(ref,['xu','yu','zu']));m['mpi_max_velocity_difference_m_s']=maximum(cols(df,['vx','vy','vz'])-cols(ref,['vx','vy','vz']));check('mpi_compare',np.array_equal(df['id'],ref['id']) and np.array_equal(df['diameter'],ref['diameter']) and m['mpi_max_position_difference_m']<=g['mpi_position_m'] and m['mpi_max_velocity_difference_m_s']<=g['mpi_velocity_m_s'])
  if sp['engine']=='gpu':
   check('kokkos_enabled','KOKKOS mode with Kokkos version' in log);m['host_fix_synchronization']={'host_fix':True,'kokkosable':counters[0]['kokkosable'],'datamask_read':counters[0]['datamask_read'],'datamask_modify':counters[0]['datamask_modify'],'mechanism':'unchanged ModifyKokkos sync and modified around POST_FORCE and END_OF_STEP host callbacks','GPU_PERFORMANCE_READY':'NO'}
   r=table(p/'GPU_RESOURCE_TRACE.csv');m['gpu_vram_peak_mib']=float(np.max(r['memory_mib']));m['gpu_utilization_peak_percent']=float(np.max(r['utilization_percent']))
  return {'status':'PASS' if all(checks.values()) else 'FAIL','case':name,'checks':checks,'metrics':m}
 except Exception as exc:
  import traceback
  return {'status':'FAIL','case':name,'checks':checks,'metrics':m,'exception':repr(exc),'traceback':traceback.format_exc()}

def audit_all(root):
 root=Path(root);unit=unit_audit(root);ct=read(root/'contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json');results=[case_audit(root,s['name']) for s in ct['cases']];lookup={x['case']:x for x in results};ds=[lookup['case_d_'+q]['metrics'] for q in ['d10','d50','d90']]
 scaling=all(d.get('measured_e_folding_time_s',0)>0 for d in ds) and all(ds[j]['measured_e_folding_time_s']<ds[j+1]['measured_e_folding_time_s'] for j in [0,1])
 return {'status':'PASS' if unit['status']=='PASS' and all(x['status']=='PASS' for x in results) and scaling else 'FAIL','unit_audit':unit,'cases':results,'diameter_relaxation_ordering':'PASS' if scaling else 'FAIL','nonfinite_count':sum(x['metrics'].get('nonfinite_count',0) for x in results),'lost_atoms':sum(x['metrics'].get('lost_atoms',0) for x in results)}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--unit',action='store_true');ap.add_argument('--case');ap.add_argument('--output',required=True);a=ap.parse_args()
 try:r=unit_audit(a.root) if a.unit else (case_audit(a.root,a.case) if a.case else audit_all(a.root))
 except Exception as e:
  import traceback
  r={'status':'FAIL','exception':repr(e),'traceback':traceback.format_exc()}
 Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps({'status':r['status'],'output':a.output}));raise SystemExit(0 if r['status']=='PASS' else 1)
