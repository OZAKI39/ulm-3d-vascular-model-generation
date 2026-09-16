"""Independent audit: frozen inputs + raw CSV/dumps + Python reference, no solver calls."""
from pathlib import Path
import json,hashlib,re,sys,argparse
import numpy as np
from reference_passive_transport import ReferenceField,velocity,rk2,dop853

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(p):return np.atleast_1d(np.genfromtxt(p,delimiter=',',names=True))
def cols(t,names):return np.column_stack([t[n] for n in names])
def maxnorm(a):return float(np.max(np.linalg.norm(a,axis=-1)))
def xyz(t):return cols(t,['x_m','y_m','z_m'])
def vel(t):return cols(t,['vx_m_s','vy_m_s','vz_m_s'])
def trajectory(p):
 a=np.concatenate([table(f) for f in sorted(Path(p).glob('trajectory_rank*.csv'))]);a.sort(order=['step','particle_id']);return a

def integrity(root):
 return all(sha(root/n)==h for h,n in (l.split('  ',1) for l in (root/'FROZEN_INPUT_SHA256SUMS').read_text().splitlines()))
def dump(p):
 with Path(p).open() as f:
  assert f.readline().strip()=='ITEM: TIMESTEP';step=int(f.readline());assert f.readline().strip()=='ITEM: NUMBER OF ATOMS';n=int(f.readline());assert f.readline().startswith('ITEM: BOX BOUNDS');box=[f.readline() for _ in range(3)];names=f.readline().split()[2:];a=np.atleast_2d(np.loadtxt(f))
 assert len(a)==n;a=a[np.argsort(a[:,names.index('id')])];return step,{name:a[:,i] for i,name in enumerate(names)}
def spec(root,name):return read(root/'cases'/name/'CASE_CONTRACT.json')
def audit_case(root,name):
 root=Path(root);p=root/'cases'/name;sp=spec(root,name);ct=read(root/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json');g=ct['gates'];checks={};metrics={}
 def ck(k,b):checks[k]=bool(b)
 try:
  ck('frozen_identity',integrity(root));ck('field_sha',sha(root/sp['field'])==sp['field_sha256']);t=trajectory(p);states=[read(f) for f in sorted(p.glob('trajectory_state_rank*.json'))]
  start,di=dump(p/'initial.dump');end,df=dump(p/'final.dump');metrics['completed_steps']=end;metrics['physical_time_s']=end*sp['dt_s'];metrics['termination_reason']=states[0]['termination_reason'];ck('state_steps',len(states)==sp['mpi_ranks'] and all(s['completed_step']==s['accepted_steps']==end and s['termination_reason']==states[0]['termination_reason'] for s in states))
  metrics['nonfinite_count']=sum(int(np.sum(~np.isfinite(t[n]))) for n in t.dtype.names)+sum(s['nonfinite_count'] for s in states);metrics['invalid_query_count']=sum(s['invalid_query_count'] for s in states)+int(np.sum(t['query_status']!=0));metrics['lost_atoms']=sp['N']-len(df['id']);ck('finite_and_valid',metrics['nonfinite_count']==metrics['invalid_query_count']==metrics['lost_atoms']==0)
  ids=np.arange(1,sp['N']+1);ck('ids_exact',np.array_equal(di['id'],ids) and np.array_equal(df['id'],ids));ck('diameter_exact',np.array_equal(di['diameter'],df['diameter']))
  pop=table(root/'populations'/(sp['population']+'.csv'));diam=pop['diameter_um'][sp['source_bubble_ids']];ck('diameter_source',np.max(np.abs(df['diameter']*1e6-diam))<=1e-10 and np.max(np.abs(diam-np.array(sp['diameters_um'])))<=1e-10);ck('initial_positions',maxnorm(cols(di,['xu','yu','zu'])-sp['initial_positions_m'])<=1e-18)
  unique,counts=np.unique(t['step'],return_counts=True);ck('record_identity',np.all(counts==sp['N']) and len(set(zip(t['step'],t['particle_id'])))==len(t));expected=np.unique(np.r_[np.arange(0,end+1,sp['trajectory_stride']),end]);ck('output_cadence',np.array_equal(unique,expected));ck('time',np.max(np.abs(t['time_s']-t['step']*sp['dt_s']))<=1e-14)
  ck('initial_final_steps',start==0 and unique[0]==0 and unique[-1]==end);last=t[t['step']==end];ck('dump_vs_csv',maxnorm(xyz(last)-cols(df,['xu','yu','zu']))<=1e-18 and maxnorm(vel(last)-cols(df,['vx','vy','vz']))<=1e-16)
  field=ReferenceField(root/sp['field']);u,status=field.sample(xyz(t));ck('independent_query_valid',np.all(status==0));metrics['sampled_velocity_error_m_s']=maxnorm(cols(t,['fluid_ux_m_s','fluid_uy_m_s','fluid_uz_m_s'])-u);metrics['particle_vs_fluid_velocity_error_m_s']=maxnorm(vel(t)-u);ck('velocity_equal_fluid',max(metrics['sampled_velocity_error_m_s'],metrics['particle_vs_fluid_velocity_error_m_s'])<=g['cpp_python_velocity_m_s']);ck('speed_fields',np.max(np.abs(t['speed_m_s']-np.linalg.norm(vel(t),axis=1)))<=1e-13 and np.max(np.abs(t['fluid_speed_m_s']-np.linalg.norm(u,axis=1)))<=1e-13)
  runtime=read(p/'RUN_METRICS.json');log=(p/'log.lammps').read_text();inp=(p/'in.lammps').read_text();ck('process',runtime['returncode']==0 and 'PASSIVE_CASE_COMPLETED' in log and 'ERROR' not in log);ck('no_inertial_fix','nve/sphere' not in inp and 'sonovue/passive/flow' in inp);ck('log_finite',re.search(r'\b(?:nan|inf)\b',log,re.I) is None)
  metrics.update(wall_seconds=runtime['wall_seconds'],peak_rss_kib=runtime['process_tree_peak_rss_kib'],particle_steps_per_s=sp['N']*end/runtime['wall_seconds'],mpi_ranks=sp['mpi_ranks'],N=sp['N'])
  if not sp['real_geometry']:ck('full_steps',end==sp['steps'] and metrics['termination_reason']=='MAX_PHYSICAL_TIME')
  else:ck('termination',metrics['termination_reason'] in ['MAX_PHYSICAL_TIME','WALL_SAFETY_STOP','SPHERE_OVERLAP_SAFETY_STOP'] and end>0)
  if sp['field_kind']=='uniform':
   syn=read(root/'contracts/SYNTHETIC_FIELDS_CONTRACT.json');uf=np.array(syn['uniform_velocity_m_s']);x0=np.array(sp['initial_positions_m'])[t['particle_id'].astype(int)-1];ref=x0+t['time_s'][:,None]*uf;metrics['max_position_error_m']=maxnorm(xyz(t)-ref);metrics['max_velocity_error_m_s']=maxnorm(vel(t)-uf);ck('uniform_analytic',metrics['max_position_error_m']<=g['uniform_position_m'] and metrics['max_velocity_error_m_s']<=g['uniform_velocity_m_s'])
  elif sp['field_kind']=='linear':
   syn=read(root/'contracts/SYNTHETIC_FIELDS_CONTRACT.json');ref,ode=dop853(field,np.array(sp['initial_positions_m'][0]),t['time_s'],affine=(syn['linear_intercept_m_s'],syn['linear_matrix_per_s']));metrics['ODE_reference']=ode;metrics['max_position_error_m']=maxnorm(xyz(t)-ref);metrics['final_position_error_m']=float(np.linalg.norm(xyz(t)[-1]-ref[-1]));ck('linear_DOP853',metrics['max_position_error_m']<=g['linear_max_position_m'] and metrics['final_position_error_m']<=g['linear_final_position_m'])
  if sp['case']==2:
   metrics['net_displacement_m']=float(np.linalg.norm(xyz(t)[-1]-xyz(t)[0]));metrics['min_surface_wall_gap_m']=float(np.min(t['surface_wall_gap_um']))*1e-6;ck('minimum_transport',metrics['net_displacement_m']>=g['case2_min_displacement_m']);ck('recorded_clearance',metrics['min_surface_wall_gap_m']>=3*ct['dx_m'])
  if sp['mpi_ranks']==4:
   counterpart=root/'cases'/name.replace('mpi4','mpi1');_,ref=dump(counterpart/'final.dump');_,ini=dump(counterpart/'initial.dump');metrics['mpi_position_difference_m']=maxnorm(cols(df,['xu','yu','zu'])-cols(ref,['xu','yu','zu']));metrics['mpi_velocity_difference_m_s']=maxnorm(cols(df,['vx','vy','vz'])-cols(ref,['vx','vy','vz']));ck('mpi_comparison',np.array_equal(df['id'],ref['id']) and np.array_equal(df['diameter'],ref['diameter']) and np.array_equal(cols(di,['xu','yu','zu']),cols(ini,['xu','yu','zu'])) and metrics['mpi_position_difference_m']<=g['mpi_position_m'] and metrics['mpi_velocity_difference_m_s']<=g['mpi_velocity_m_s'])
   metrics['particles_with_observed_MPI_migration']=sum(len(np.unique(t['owner_rank'][t['particle_id']==i]))>1 for i in ids);ck('migration_exercised',metrics['particles_with_observed_MPI_migration']>0)
  if sp['engine']=='gpu':
   ck('kokkos_enabled','KOKKOS mode with Kokkos version' in log);r=table(p/'GPU_RESOURCE_TRACE.csv');metrics['gpu_utilization_peak_percent']=float(r['utilization_percent'].max());metrics['gpu_vram_peak_mib']=float(r['memory_mib'].max());metrics['GPU_PERFORMANCE_READY']='NO'
  return {'status':'PASS' if all(checks.values()) else 'FAIL','name':name,'checks':checks,'metrics':metrics}
 except Exception as e:
  import traceback
  return {'status':'FAIL','name':name,'checks':checks,'metrics':metrics,'exception':repr(e),'traceback':traceback.format_exc()}

def convergence(root):
 root=Path(root);ct=read(root/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json');g=ct['gates'];details={}
 for case,kind in [(0,'uniform'),(1,'linear'),(2,'real')]:
  series={};base=.5*ct['field_speed_statistics'][kind]['dt_base_s']
  for suffix in ['c050','c025','c0125']:
   t=trajectory(root/'cases'/f'case{case}_{suffix}');k=np.rint(t['time_s']/base).astype(int);ok=np.abs(t['time_s']/base-k)<1e-6;series[suffix]={int(i):x for i,x in zip(k[ok],xyz(t)[ok])}
  result={}
  for left,right in [('c050','c025'),('c025','c0125')]:
   keys=sorted(set(series[left])&set(series[right]));assert len(keys)>=2
   diff=np.array([series[left][i]-series[right][i] for i in keys]);mx=maxnorm(diff);final=float(np.linalg.norm(diff[-1]));result[left+'_vs_'+right]={'max_trajectory_difference_m':mx,'final_position_difference_m':final,'comparison_final_time_s':keys[-1]*base,'common_records':len(keys),'status':'PASS' if mx<=g['dt_max_trajectory_difference_m'] and final<=g['dt_final_difference_m'] else 'FAIL'}
  details[kind]=result
 selected=.25 if all(d['c025_vs_c0125']['status']=='PASS' for d in details.values()) else .125
 stable=all(audit_case(root,f'case{i}_c0125')['status']=='PASS' for i in range(3))
 return {'status':'PASS' if stable else 'FAIL','selected_c_adv':selected,'particle_dt_s':selected*ct['field_speed_statistics']['real']['dt_base_s'],'details':details,'comparison_policy':'same recorded physical times only, ending at last common safe record; differing safety-stop tail timestamps reported separately, never compare different physical times','selection_policy':ct['selection_policy'],'C_adv_tested':[.5,.25,.125]}

def audit_all(root):
 root=Path(root);selection=read(root/'contracts/PASSIVE_TRANSPORT_TIMESTEP_CONTRACT.json');label='c025' if selection['selected_c_adv']==.25 else 'c0125';names=[f'case{i}_{s}' for i in range(3) for s in ['c050','c025','c0125']]+['case3_'+q+'_'+label for q in ['d10','d50','d90']]+['case4_mpi1_'+label,'case4_mpi4_'+label,'kokkos_'+label]
 if (root/'cases'/('case5_'+label)/'RUN_METRICS.json').exists():names.append('case5_'+label)
 results=[audit_case(root,n) for n in names];conv=convergence(root);ds=[trajectory(root/'cases'/('case3_'+q+'_'+label)) for q in ['d10','d50','d90']];dx=max(maxnorm(xyz(ds[i])-xyz(ds[0])) for i in [1,2]);dv=max(maxnorm(vel(ds[i])-vel(ds[0])) for i in [1,2]);g=read(root/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json')['gates'];diam={'status':'PASS' if dx<=g['diameter_position_m'] and dv<=g['diameter_velocity_m_s'] else 'FAIL','max_position_difference_m':dx,'max_velocity_difference_m_s':dv};
 return {'status':'PASS' if all(a['status']=='PASS' for a in results) and conv['status']==diam['status']=='PASS' else 'FAIL','cases':results,'timestep_convergence':conv,'diameter_independence':diam,'selected_label':label,'geometry_audit':'PENDING_LOCAL_EXACT_STL_AND_RK2_REPLAY','nonfinite_count':sum(a['metrics'].get('nonfinite_count',0) for a in results),'invalid_query_count':sum(a['metrics'].get('invalid_query_count',0) for a in results),'lost_atoms':sum(a['metrics'].get('lost_atoms',0) for a in results)}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--case');ap.add_argument('--output',required=True);args=ap.parse_args();r=audit_case(args.root,args.case) if args.case else audit_all(args.root);Path(args.output).write_text(json.dumps(r,indent=2)+'\n');print(r['status']);raise SystemExit(0 if r['status']=='PASS' else 1)
