from common import *
import numpy as np,pyvista as pv,re
case=HERE/'runtime/p1_candidate_fixed_observer/run';protocol=read(HERE/'LINEAR_SOLVER_CANDIDATE_PLAN.json')['P1_protocol']
receipts=sorted((HERE/'monitoring').glob('p1_candidate*/EXECUTION.json'),key=lambda p:p.stat().st_mtime);rp=receipts[-1];rec=read(rp);log=(rp.parent/'job.log').read_text()
base=read(PRIOR/'P1_LINEAR_SOLVER_PARITY.json')['results']['baseline']
def matrix(p):
 with p.open('rb') as f:
  tag,m,n,nz=np.fromfile(f,dtype='>i4',count=4);assert tag==1211216
  row=np.fromfile(f,dtype='>i4',count=m);col=np.fromfile(f,dtype='>i4',count=nz);a=np.fromfile(f,dtype='>f8',count=nz)
  tag,nb=np.fromfile(f,dtype='>i4',count=2);assert tag==1211214 and nb==m
  rhs=np.fromfile(f,dtype='>f8',count=nb);assert not f.read()
 return row,col,a,rhs
result=dict(execution_success=rec['exit_code']==0 and rec['stop_reason'] is None and 'PETSC ERROR' not in log and 'DIVERGED_' not in log,receipt=str(rp.relative_to(HERE)),wall_time_s=rec['elapsed_s'],peak_tree_RSS_GiB=rec['peak_tree_RSS_GiB'],peak_GPU_total_MiB=max(x.get('GPU_total_used_MiB',0) for x in rec['samples']),KSP_iterations=[int(x) for x in re.findall(r'Linear solve converged due to \w+ iterations (\d+)',log)],nonlinear_iterations=len(re.findall(r'^\s*NS\s+.*$',log,re.M)),CPU_Pmat_identity_checks=log.count('CPU_PMAT_IDENTICAL 1'),PASS=False)
result['total_KSP_iterations']=sum(result['KSP_iterations'])
if result['execution_success']:
 ref=matrix(PRIOR/'runtime/p1_baseline/run/first_system.petsc');got=matrix(case/'first_system.petsc');structure=all(np.array_equal(ref[i],got[i]) for i in [0,1]);rhs=ref[3].tobytes()==got[3].tobytes();ptr=np.r_[0,np.cumsum(ref[0])];scale=np.maximum.reduceat(abs(ref[2]),ptr[:-1]);delta=np.maximum.reduceat(abs(got[2]-ref[2]),ptr[:-1]);rel=float(np.max(delta/np.maximum(scale,np.finfo(float).tiny)))
 ident=dict(CSR_bitwise_equal=structure,RHS_bitwise_equal=rhs,byte_identical=sha(case/'first_system.petsc')==sha(PRIOR/'runtime/p1_baseline/run/first_system.petsc'),max_A_abs_difference=float(np.max(abs(got[2]-ref[2]))),max_row_scaled_difference=rel,row_roundoff_units=rel/np.finfo(float).eps,PASS=bool(structure and rhs and rel<=128*np.finfo(float).eps))
 put('P1_MATRIX_RHS_COMPARISON.json',ident)
 old=np.load(PRIOR/'P1_baseline_FINAL_STATE.npz');v=pv.read(case/'8-procs/result_011.vtu');u=np.asarray(v['Velocity']);p=np.asarray(v['Pressure']);du=float(np.max(abs(u-old['velocity'])));dp=float(np.max(abs(p-old['pressure'])));state=np.column_stack([u,p]);oldstate=np.column_stack([old['velocity'],old['pressure']]);dr=float(np.linalg.norm(state-oldstate)/np.linalg.norm(oldstate));np.savez_compressed(HERE/'P1_candidate_FINAL_STATE.npz',velocity=u,pressure=p)
 mesh=pv.read(P1/'SV_MESH/mesh-surfaces/OUTLET_01.vtp');ids=np.asarray(mesh['GlobalNodeID'])-1;tri=mesh.faces.reshape(-1,4)[:,1:];x=mesh.points;area=np.linalg.norm(np.cross(x[tri[:,1]]-x[tri[:,0]],x[tri[:,2]]-x[tri[:,0]]),axis=1)/2;Q=float(np.sum(area*u[ids[tri],2].mean(1)));pts=pv.read(P1/'SV_MESH/mesh-complete.mesh.vtu').points;mask=(pts[:,2]>.4*.00716135075978932)&(pts[:,2]<.6*.00716135075978932);dpdz=float(np.polyfit(pts[mask,2],p[mask],1)[0])
 records=[json.loads(s) for f in case.glob('algebra_rank_*.jsonl') for s in f.read_text().splitlines()];cpu=[r for r in records if r['stage'].startswith('CPU_') or r['stage']=='after_PC_setup_Pmat'];cpuok=bool(cpu) and all(r['matrix_type'] in {'mpiaij','seqaij'} for r in cpu)
 norm=re.findall(r'^\s*\d+ KSP unpreconditioned resid norm ([^ ]+)',log,re.M)
 result.update(matrix_identity=ident,velocity_max_abs_difference_m_s=du,pressure_max_abs_difference_Pa=dp,relative_state_L2_difference=dr,Q_outlet_m3_s=Q,dpdz_Pa_m=dpdz,Q_relative_difference=abs(Q-base['Q_outlet_m3_s'])/abs(base['Q_outlet_m3_s']),dpdz_relative_difference=abs(dpdz-base['dpdz_Pa_m'])/abs(base['dpdz_Pa_m']),last_logged_residual_norm=float(norm[-1]),CPU_PC_hierarchy_verified=cpuok,rank_peak_RSS_bytes={str(i):max(r['peak_RSS_bytes'] for r in records if r['rank']==i) for i in range(8)})
 result['PASS']=bool(ident['PASS'] and cpuok and du<=protocol['velocity_absolute_tolerance_m_s'] and dp<=protocol['pressure_absolute_tolerance_Pa'] and dr<=protocol['relative_state_tolerance'] and result['CPU_Pmat_identity_checks']==result['nonlinear_iterations'] and 0<result['total_KSP_iterations']<=10*base['total_KSP_iterations'] and result['nonlinear_iterations']<=24)
put('P1_LINEAR_SOLVER_PARITY.json',dict(status='PASS' if result['PASS'] else 'NO_SAFE_LOW_MEMORY_LINEAR_SOLVER_CANDIDATE',baseline_inherited_from=str((PRIOR/'P1_LINEAR_SOLVER_PARITY.json').relative_to(BASE)),baseline=base,candidate=result,selected_candidate='schur_cpu_pc' if result['PASS'] else None,protocol=protocol,attempt1='Observer incorrectly expected Schur matrix-free finest operator to be AIJ; fixed to inspect actual preconditioning matrix. No tolerance or algorithm change; failed logs retained.'))
print(result,flush=True)
