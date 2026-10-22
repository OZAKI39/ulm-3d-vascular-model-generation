"""One independent CFD case, actual frozen solver; resumable via saved checkpoints.
Run under the validation-only supervisor, never writes production case paths.
"""
from pathlib import Path
import sys,os,json,time,subprocess,hashlib,argparse,re,traceback,signal
import numpy as np
import pyvista as pv
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'production'))
from flow_solver_support import wss
from flow_solver_support.wss_case import coordinate_identity
from log_acceptance import classify

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,a):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.pending');tmp.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n');tmp.replace(p)

def terminate_case_processes(run):
 """MPI ranks may have their own process groups; select only this fresh case cwd."""
 stopped=[]
 for proc in Path('/proc').iterdir():
  if not proc.name.isdigit():continue
  try:
   cwd=(proc/'cwd').resolve();cmd=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
  except OSError:continue
  if cwd==run and any(name in cmd for name in ['svmultiphysics','mpiexec','mpirun']):
   try:os.kill(int(proc.name),signal.SIGTERM);stopped.append(dict(pid=int(proc.name),command=cmd))
   except ProcessLookupError:pass
 return stopped

def solver_resident_memory(run):
 total=0
 for proc in Path('/proc').iterdir():
  if not proc.name.isdigit():continue
  try:
   if (proc/'cwd').resolve()!=run:continue
   cmd=(proc/'cmdline').read_bytes()
   if b'/svmultiphysics\x00' not in cmd:continue
   total+=int((proc/'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
  except (OSError,IndexError,ValueError):continue
 return total

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--case',type=Path,required=True);args=parser.parse_args();case=args.case.resolve();run=case/'run';reports=case/'reports';reports.mkdir(exist_ok=True)
 def interrupted(signum,frame):
  stopped=terminate_case_processes(run);write(reports/'external_interruption.json',dict(signal=signum,stopped_case_processes=stopped,unix=time.time()));raise SystemExit(128+signum)
 signal.signal(signal.SIGTERM,interrupted)
 signal.signal(signal.SIGINT,interrupted)
 assert 'wss_validation_v2' in str(case) and not (run/'solver.log').exists(),'Independent fresh case required'
 for rel,h in json.loads((case/'input_hashes.json').read_text()).items():assert sha(case/rel)==h,rel
 policy=json.loads((case/'policy.json').read_text());dt=policy['dt_s'];ranks=int(policy.get('MPI_ranks',1));result_dir=run/(str(ranks)+'-procs');mesh=np.load(case/'SV_MESH/mesh_arrays.npz');x=mesh['points_m'];t=mesh['tetra'];b=mesh['boundary_triangles'];tags=mesh['facet_tags'];vol=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6
 wall=b[tags==1];own=wss.boundary_owners(t,wall);cent,area,norm=wss.wall_geometry(x,t,wall,own)
 bx=x[b];av=.5*np.cross(bx[:,1]-bx[:,0],bx[:,2]-bx[:,0]);wallids=np.unique(wall)
 def l2(u):
  a=u[t];return np.sqrt(np.dot(vol,(np.sum(a*a,axis=tuple(range(1,a.ndim)))+np.sum(a.sum(axis=1)**2,axis=-1))/20)) if u.ndim==2 else np.sqrt(np.dot(vol,(np.sum(a*a,axis=1)+a.sum(axis=1)**2)/20))
 def load(p):
  g=pv.read(p);coordinate_identity(g.points,x);return np.asarray(g['Velocity']),np.asarray(g['Pressure']).reshape(-1)
 def wvalues(u):return np.linalg.norm(wss.tangential_traction(wss.p1_gradients(x,t[own],u),norm,policy.get('mu_Pa_s',.00345312)),axis=1)
 buildpath=Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/reports/svmp_reuse_build.json');build=json.loads(buildpath.read_text());ps=build['PETSc_build'];exe=build['executable']
 assert sha(exe)=='0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7'
 assert sha(Path(ps['prefix'])/'lib/libpetsc.so')=='b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef'
 env={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
 env.update(PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',LD_LIBRARY_PATH=build['runtime_library_path'],PETSC_OPTIONS=(run/'PETSC_OPTIONS.txt').read_text().strip())
 cmd=[ps['candidate_wrapper'],ps['launcher'],'-n',str(ranks),ps['candidate_wrapper'],exe,'solver.xml']
 write(reports/'launch.json',dict(command=cmd,solver_sha256=sha(exe),build_sha256=sha(buildpath),PETSC_OPTIONS=env['PETSC_OPTIONS'],runner_sha256=sha(__file__),log_classifier_sha256=sha(HERE/'log_acceptance.py'),flow_parser_sha256=sha(HERE/'production/flow_solver_support/flow_parser.py'),wss_core_sha256=sha(wss.__file__),case_policy_sha256=sha(case/'policy.json'),config_sha256=sha(run/'solver.xml'),input_manifest_sha256=sha(case/'input_hashes.json'),initial_state='zero',MPI_ranks=ranks,dt_s=dt,started_unix=time.time()))
 previous=None;states=[];seen=set();steady=0;stop=None;read_errors=[];start=time.time();early_failure=None;log_offset=0;log_tail='';peak_solver_RSS=0
 def observe(p):
  nonlocal previous,steady,stop
  step=int(p.stem.rsplit('_',1)[1]);u,pr=load(p);ww=wvalues(u);fl=np.einsum('ij,ij->i',av,u[b].mean(axis=1));flows={str(int(tag)):float(fl[tags==tag].sum()) for tag in np.unique(tags)}
  row=dict(step=step,time_s=step*dt,epsilon_Q=abs(-flows['4']-policy['Q_target_m3_s'])/policy['Q_target_m3_s'],epsilon_mass=abs(sum(flows.values()))/policy['Q_target_m3_s'],outward_flows_m3_s=flows,wall_max_speed_m_s=float(np.linalg.norm(u[wallids],axis=1).max()),wss_area_mean_Pa=float(np.average(ww,weights=area)),max_speed_m_s=float(np.linalg.norm(u,axis=1).max()))
  if previous is not None:
   old_u,old_p,old_w=previous
   row.update(velocity_relative_change=float(l2(u-old_u)/max(l2(u),1e-50)),pressure_relative_change=float(l2(pr-old_p)/max(l2(pr),np.sqrt(vol.sum()))),wss_area_L2_relative_change=float(np.sqrt(np.average((ww-old_w)**2,weights=area)/np.average(ww**2,weights=area))))
   passed=max(row[k] for k in ['velocity_relative_change','pressure_relative_change','wss_area_L2_relative_change'])<policy['steady_change_limit'] and row['epsilon_mass']<=policy['mass_limit'] and row['epsilon_Q']<=1e-6 and row['wall_max_speed_m_s']<=1e-12
   steady=steady+1 if passed else 0
   row['steady_interval_pass']=bool(passed)
  previous=(u.copy(),pr.copy(),ww.copy());states.append(row);seen.add(step)
  if step>=policy['minimum_steps'] and steady>=policy['steady_consecutive_intervals'] and stop is None:
   stop=dict(reason='consecutive velocity, pressure, WSS and mass checks',step=step,consecutive=steady);(run/'STOP_SIM').write_text('0\n')
 with (run/'solver.log').open('w') as log:
  proc=subprocess.Popen(cmd,cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  while proc.poll() is None:
   rss=solver_resident_memory(run);peak_solver_RSS=max(peak_solver_RSS,rss)
   if rss>policy.get('maximum_solver_RSS_bytes',float('inf')):
    early_failure=dict(reason='Registered resident memory budget exceeded; incomplete solve, not a numerical precision failure',RSS_bytes=rss,limit_bytes=policy['maximum_solver_RSS_bytes'],terminated_case_processes=terminate_case_processes(run));break
   with (run/'solver.log').open(errors='replace') as live:
    live.seek(log_offset);log_tail+=live.read();log_offset=live.tell()
   bad=re.search(r'FRESH_RETRY_FAILED|^\s*NS\s+\d+-[^\n]*WARNING: The linear system solution has not converged',log_tail,re.M)
   if bad:
    early_failure=dict(reason='Unrecovered failed linear correction: no qualified nonlinear state',marker=bad[0],observed_elapsed_s=time.time()-start)
    early_failure['terminated_case_processes']=terminate_case_processes(run);break
   log_tail=log_tail[-2000:]
   for p in sorted(result_dir.glob('result_*.vtu'),key=lambda p:int(p.stem.rsplit('_',1)[1])):
    step=int(p.stem.rsplit('_',1)[1])
    if step==0 or step in seen or step%5 or time.time()-p.stat().st_mtime<1:continue
    try:observe(p)
    except Exception as e:read_errors.append(dict(step=step,error=str(e)));continue
   if time.time()-start>policy['maximum_wall_time_s'] and stop is None:stop=dict(reason='wall_time_limit_without_convergence');(run/'STOP_SIM').write_text('0\n')
   write(reports/'progress.json',dict(pid=proc.pid,peak_solver_RSS_bytes=peak_solver_RSS,elapsed_s=time.time()-start,states=states,stop=stop,read_errors=read_errors))
   time.sleep(2)
  code=proc.wait()
 outputs=sorted(result_dir.glob('result_*.vtu'),key=lambda p:int(p.stem.rsplit('_',1)[1]))
 if not outputs:write(reports/'execution.json',dict(status='FAIL',exit_code=code,no_results=True,early_failure=early_failure,solver_log_sha256=sha(run/'solver.log'),elapsed_s=time.time()-start));raise RuntimeError('No result VTU; inspect solver.log')
 final=outputs[-1];step=int(final.stem.rsplit('_',1)[1])
 if step not in seen:observe(final)
 u,pr=load(final);frozen=case/'frozen_flow';frozen.mkdir(exist_ok=True)
 import shutil
 shutil.copy2(final,frozen/'steady_flow.vtu');np.savez_compressed(frozen/'flow_arrays_si.npz',points_m=x,tetra=t,boundary_triangles=b,facet_tags=tags,velocity_m_s=u,pressure_pa=pr)
 write(frozen/'manifest.json',dict(case=case.name,step=step,input_configuration_sha256=sha(run/'solver.xml'),mesh_arrays_sha256=sha(case/'SV_MESH/mesh_arrays.npz'),policy_sha256=sha(case/'policy.json'),input_manifest_sha256=sha(case/'input_hashes.json'),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in frozen.iterdir() if p.is_file() and p.name!='manifest.json'}))
 text=(run/'solver.log').read_text(errors='replace');diverged=re.findall(r'Linear solve did not converge[^\n]*|DIVERGED_[A-Z_]+|NaN|nan\b|Segmentation fault',text)
 acceptance=classify(text,dt);write(reports/'linear_attempt_acceptance.json',acceptance)
 reasons=re.findall(r'CONVERGED_[A-Z_]+',text);record=dict(status='PASS' if not early_failure and code==0 and stop and 'step' in stop and acceptance['accepted'] else 'FAIL',exit_code=code,early_failure=early_failure,peak_solver_RSS_bytes=peak_solver_RSS,final_step=step,final_vtu=str(final),elapsed_s=time.time()-start,stop=stop,states=states,read_errors=read_errors,converged_linear_reasons=len(reasons),divergence_markers=diverged,linear_attempt_acceptance=acceptance,solver_log_sha256=sha(run/'solver.log'),input_hashes_unchanged=all(sha(case/k)==h for k,h in json.loads((case/'input_hashes.json').read_text()).items()))
 write(reports/'execution.json',record);print(json.dumps({k:v for k,v in record.items() if k not in ['states','read_errors']},indent=2),flush=True)
 if record['status']!='PASS':raise SystemExit(2)

if __name__=='__main__':main()
