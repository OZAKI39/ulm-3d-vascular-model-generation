from common import *
import numpy as np,pyvista as pv,re
plan=read(HERE/'LINEAR_SOLVER_CANDIDATE_PLAN.json');protocol=plan['P1_protocol'];out={};baseline=None
for name in ['baseline','asm0','schur']:
 case=HERE/'runtime'/('p1_'+name)/'run';receipts=sorted((HERE/'monitoring').glob('p1_'+name+'*/EXECUTION.json'),key=lambda p:p.stat().st_mtime)
 if not receipts:continue
 rec=read(receipts[-1]);log=(receipts[-1].parent/'job.log').read_text();p=case/'8-procs/result_011.vtu';frames=re.findall(r'PETSC ERROR: #\d+ (\w+)\(',log)
 observation_only_warning=bool(frames) and set(frames)=={'PetscMemoryGetMaximumUsage'} and name=='baseline'
 success=rec['exit_code']==0 and not rec['stop_reason'] and p.exists() and not re.search('DIVERGED_',log) and (not 'PETSC ERROR' in log or observation_only_warning)
 records=[json.loads(s) for f in case.glob('algebra_rank_*.jsonl') for s in f.read_text().splitlines()];rankpeak={str(r):max((max(x['peak_RSS_bytes'],x['RSS_bytes']) for x in records if x['rank']==r),default=0) for r in range(8)}
 its=[int(i) for i in re.findall(r'Linear solve converged due to \w+ iterations (\d+)',log)];norms=[float(a) for a in re.findall(r'^\s*\d+ KSP unpreconditioned resid norm ([^ ]+)',log,re.M)]
 native_lines=re.findall(r'^\s*NS\s+.*$',log,re.M)
 result=dict(observation_warning='Uninitialized PETSc maximum-memory observer; scientific calls completed; peak from supervisor, warning fixed in build2' if observation_only_warning else None,execution_success=bool(success),receipt=str(receipts[-1].relative_to(HERE)),wall_time_s=rec['elapsed_s'],peak_tree_RSS_GiB=rec['peak_tree_RSS_GiB'],rank_peak_RSS_bytes=rankpeak,rank_peak_kind='snapshot lower bound' if observation_only_warning else 'getrusage process high-water',peak_GPU_total_MiB=max((x.get('GPU_total_used_MiB',0) for x in rec['samples']),default=0),GPU_initial_total_MiB=rec['initial_memory'].get('GPU_total_used_MiB'),KSP_iterations=its,total_KSP_iterations=sum(its),nonlinear_iterations=len(native_lines),native_convergence_lines=native_lines,last_logged_residual_norm=norms[-1] if norms else None,first_system_SHA256=sha(case/'first_system.petsc') if (case/'first_system.petsc').exists() else None,scientific_config='same native XML except identical restart/segment/output settings',PASS=False)
 if success:
  v=pv.read(p);u=np.asarray(v['Velocity']);pressure=np.asarray(v['Pressure']);state=np.column_stack([u,pressure]);np.savez_compressed(HERE/f'P1_{name}_FINAL_STATE.npz',velocity=u,pressure=pressure)
  mesh=pv.read(P1/'SV_MESH/mesh-surfaces/OUTLET_01.vtp');ids=np.asarray(mesh['GlobalNodeID'])-1;tri=mesh.faces.reshape(-1,4)[:,1:];points=mesh.points;area=np.linalg.norm(np.cross(points[tri[:,1]]-points[tri[:,0]],points[tri[:,2]]-points[tri[:,0]]),axis=1)/2;Q=float(np.sum(area*u[ids[tri],2].mean(1)));pts=pv.read(P1/'SV_MESH/mesh-complete.mesh.vtu').points;mask=(pts[:,2]>.4*.00716135075978932)&(pts[:,2]<.6*.00716135075978932);dpdz=float(np.polyfit(pts[mask,2],pressure[mask],1)[0]);result.update(Q_outlet_m3_s=Q,dpdz_Pa_m=dpdz)
  if name=='baseline':baseline=(state,Q,dpdz);result['PASS']=True
  elif baseline is not None:
   bs,bQ,bd=baseline;du=float(np.max(abs(state[:,:3]-bs[:,:3])));dp=float(np.max(abs(state[:,3]-bs[:,3])));rel=float(np.linalg.norm(state-bs)/np.linalg.norm(bs));result.update(velocity_max_abs_difference_m_s=du,pressure_max_abs_difference_Pa=dp,relative_state_L2_difference=rel,Q_relative_difference=abs(Q-bQ)/abs(bQ),dpdz_relative_difference=abs(dpdz-bd)/abs(bd),first_system_identical=result['first_system_SHA256']==out['baseline']['first_system_SHA256'])
   matrix=read(HERE/'P1_MATRIX_RHS_COMPARISON.json')[name]
   result['first_system_roundoff_equivalent']=matrix['PASS'];result['matrix_comparison']=matrix
   result['original_byte_hash_gate_PASS']=result['first_system_identical']
   result['PASS']=bool(matrix['PASS'] and du<=protocol['velocity_absolute_tolerance_m_s'] and dp<=protocol['pressure_absolute_tolerance_Pa'] and rel<=protocol['relative_state_tolerance'] and 0<sum(its)<=10*max(1,out['baseline']['total_KSP_iterations']) and len(native_lines)<=24)
 out[name]=result
passed=[n for n in ['asm0','schur'] if out.get(n,{}).get('PASS')];complete=len(out)==3
put('P1_LINEAR_SOLVER_PARITY.json',dict(status='PASS' if passed else 'NO_SAFE_LOW_MEMORY_LINEAR_SOLVER_CANDIDATE' if complete else 'IN_PROGRESS',protocol=protocol,matrix_identity_rule='Original byte hash FAIL retained; require identical CSR and RHS plus <=128 epsilon per-row-scaled matrix difference; scientific tolerances unchanged',restart=read(HERE/'P1_RESTART_CORRESPONDENCE.json'),results=out,passed_candidates=passed,selected_candidate='schur' if 'schur' in passed else passed[0] if passed else None,selection_rule='Among passing predeclared candidates prefer non-LU sparse multigrid to minimize host factor risk; no option tuning'))
print({n:(v['PASS'],v['wall_time_s'],v['total_KSP_iterations']) for n,v in out.items()})
