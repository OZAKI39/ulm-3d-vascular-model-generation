from common import *
import sys,os,time,subprocess,shutil,resource
sys.path.insert(0,str(R))
import scipy
from analysis.dual_track_v1.core_a import Energy,delta_y
t=time.perf_counter();O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
assert not (O/'LOW_LOAD_DESIGN_FREEZE.json').exists()
protect=json.loads((OLD/'HISTORICAL_PROTECTION_BEFORE.json').read_text())
for p,v in protect.items():assert sha(p)==v['sha256'],p
for base in [OLD,R/'results'/C['prior_mode'].split('/')[1]/C['prior_mode'].split('/')[2],R/'analysis/mode_i_v1',R/'analysis/dual_track_v1']:
 for p in base.rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:protect[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
dump(O/'HISTORICAL_PROTECTION_BEFORE.json',protect)
cmds={}
for k,cmd in [('branch',['git','branch','--show-current']),('HEAD',['git','rev-parse','HEAD']),('status',['git','status','--short']),('cpu',['lscpu']),('RAM',['free','-b']),('disk',['df','-B1',str(O)]),('jobs',['ps','-eo','pid,ppid,etimes,pcpu,rss,comm']),('GPU',['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu','--format=csv'])]:
 q=subprocess.run(cmd,capture_output=True,text=True);cmds[k]={'command':cmd,'returncode':q.returncode,'stdout':q.stdout,'stderr':q.stderr}
paths=['/sys/fs/cgroup/cpu/cpu.cfs_quota_us','/sys/fs/cgroup/cpu/cpu.cfs_period_us','/sys/fs/cgroup/memory/memory.limit_in_bytes']
dump(O/'PREFLIGHT.json',dict(cwd=os.getcwd(),python=sys.executable,python_version=sys.version,scipy_version=scipy.__version__,affinity=sorted(os.sched_getaffinity(0)),cgroup={p:Path(p).read_text() for p in paths},commands=cmds,threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},GPU_used=False))
ident=json.loads((OLD/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text())
for p,v in ident['files'].items():assert sha(p)==v['sha256'],p
for p in [R/'analysis/dual_track_v1/hvp_y.cpp',R/'analysis/dual_track_v1/hvp_y.so']:
 ident['files'][str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
ident.update(actual_LAMMPS_instance_loaded=False,scope='Current files checked; runtime mapping pending this task',previous_task_reference=str(OLD))
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',ident)
assert sha(OLD/'reference_topology.npz')=='ef145b6d6db6c17584dce6cbea1fa1b9d95847d71a69265a73df6334d1695c03'
shutil.copy2(OLD/'reference_topology.npz',O/'reference_topology.npz')
rt=np.load(O/'reference_topology.npz');X=rt['X'];b=rt['bonds'];alive=np.ascontiguousarray(rt['alive'][0],dtype=np.uint8);e=Energy(dict(X=X,bonds=b))
z0=np.load(OLD/'native_initial_A0.npz');prior=np.load(OLD/'elastic_precheck.npz')
assert np.array_equal(prior['X'],X) and np.array_equal(prior['bonds'],b) and np.array_equal(prior['alive'],alive)
# Native zero-load coordinates are checked against the exact original reference.
assert np.array_equal(z0['x'],X)
assert np.max(abs(e.m-rt['m'])/rt['m'])<1e-12
assert len(X)==16000 and len(b)==909392 and alive.sum()==892976
assert np.array_equal(alive==0,rt['cross_z0']&(rt['reference_intersection_x_m']<-.001-1e-12))
lo=np.isclose(X[:,2],X[:,2].min(),atol=1e-12,rtol=0);hi=np.isclose(X[:,2],X[:,2].max(),atol=1e-12,rtol=0)
assert np.array_equal(lo,rt['low']) and np.array_equal(hi,rt['high']) and lo.sum()==hi.sum()==800
spec=json.loads((OLD/'SPECIMEN_AND_PRECRACK_AUDIT.json').read_text());spec.update(original_atom_IDs_unchanged=True,accepted_zero_load_source=str(OLD/'native_initial_A0.npz'),accepted_zero_load_sha256=sha(OLD/'native_initial_A0.npz'),nodal_volume_m3=V,nodal_mass_kg=1000*V,reference_m_recomputed_relative=float(np.max(abs(e.m-rt['m'])/rt['m'])))
dump(O/'SPECIMEN_IDENTITY_AUDIT.json',spec)
fext=np.zeros_like(X);fext[hi,2]=2*DX**2;fext[lo,2]=-2*DX**2
np.savez_compressed(O/'load_identity.npz',atom_id=np.arange(1,len(X)+1),external_force_N=fext,top_atom_ids=np.flatnonzero(hi)+1,bottom_atom_ids=np.flatnonzero(lo)+1,precut_reference_bond_ids=np.flatnonzero(alive==0)+1)
dump(O/'LOW_LOAD_FORCE_APPLICATION_AUDIT.json',dict(status='FROZEN_NATIVE_APPLICATION_PENDING',traction_Pa=2.,per_node_force_N=2*DX**2,face_area_m2=int(hi.sum())*DX**2,top_force_N=fext[hi].sum(0),bottom_force_N=fext[lo].sum(0),global_force_N=fext.sum(0),reference_torque_N_m=np.cross(X,fext).sum(0),actual_native_application_pending=True))
# Read-only 20 Pa diagnostics: no minimization or new loading case.
x=prior['x'];U=e.eval(x,alive);p20=prior['external_force_N'];P=U-np.sum(p20*(x-X));s=np.linalg.norm(delta_y(x,b),axis=1)/e.r-1
shift=np.array([3e-6,-7e-6,2e-6]);xt=x+shift;Ut=e.eval(xt,alive);Pt=Ut-np.sum(p20*(xt-X));st=np.linalg.norm(delta_y(xt,b),axis=1)/e.r-1
rot=[]
for axis in range(3):
 v=np.eye(3)[axis];A=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]);a=1e-4;Q=np.eye(3)+np.sin(a)*A+(1-np.cos(a))*(A@A);xr=x@Q.T;Ur=e.eval(xr,alive);sr=np.linalg.norm(delta_y(xr,b),axis=1)/e.r-1
 rot.append(dict(axis='XYZ'[axis],angle_rad=a,periodic_vector_change_m=Q@np.array([0,PERIOD,0])-np.array([0,PERIOD,0]),U_change_J=Ur-U,Pi_change_J=Ur-np.sum(p20*(xr-X))-P,max_stretch_change=float(abs(sr-s).max())))
dump(O/'RIGID_MODE_AND_GAUGE_AUDIT.json',dict(status='TRANSLATION_GAUGE_JUSTIFIED',source_state='Saved UNACCEPTED 20 Pa endpoint; read-only transformation, not an equilibrium solve',translation_shift_m=shift,U_change_J=Ut-U,Pi_change_J=Pt-P,max_stretch_change=float(abs(st-s).max()),external_force_sum_N=p20.sum(0),external_torque_reference_N_m=np.cross(X,p20).sum(0),external_torque_saved_state_N_m=np.cross(x,p20).sum(0),rotation_probes=rot,interpretation='Only 3 translations are removed. Fixed Y lattice vector excludes rigid rotations about X/Z. Rotation about Y preserves internal bond lengths but the fixed laboratory dead-load potential is not invariant, so it is not removed. No physical support force is added.',local_fixed_nodes=0,removed_translation_modes=3,removed_rotation_modes=0))
audit='''# 20 Pa 提前停止：只读诊断\n\n直接证据：旧ELASTIC_PRECHECK.json记录L-BFGS-B因RELATIVE REDUCTION OF F <= FACTR*EPSMCH停止；最大节点残差1.954014e-12 N超过1e-13 N。相对残差虽通过，不覆盖绝对门槛。\n\n旧变量q=u/(1e-6 m)，目标Pi/(1e-12 J)，梯度=-F×1e6。gtol=1e-10对应分量1e-16 N，实际未达到；ftol=2.220446049250313e-15允许目标下降足够小时先退出。旧Pi约-3.115687e-9 J，缩放目标约-3115.687，终止下降尺度约6.92e-24 J。最终评估日志在该尺度内往返，力仍未合格。能量不再可分辨不能证明力平衡。\n\n旧日志是1132次函数评估，不是接受迭代序列。没有完整方向、步长或Wolfe诊断，不能确定哪次线搜索触发、精确条件数、L-BFGS历史矩阵是否恶化。高载几何非线性、柔顺裂纹臂、体积/剪切刚度差造成的病态性只是候选解释。未重跑20 Pa。\n\n独立扩展LPS储能包含两端体积V²、完整参考m；物理力已是N。上一轮原生内部力经一次V转换与独立力相对差1.73e-13，未见体积漏乘或重乘证据。当前沿用同一源码，不改材料或原生算法。\n\n三个全局平移是精确参考系自由度，总外力为零。数值平移与三个转动探针见RIGID_MODE_AND_GAUGE_AUDIT.json；固定Y周期盒排除X/Z转动，Y转动虽保持内部能，固定方向dead load势能不保持。仅对全部节点均匀平移作质心规约，所有48000个材料坐标仍参与力学，不钉扎节点、不投影转动。最终验收采用未投影原生残余，不隐去支承力。\n\n新主路径：q=(u-mean(u))/(1e-6 m)，Pi/(1e-11 J)，梯度=-P(F)×1e5，ftol=0，gtol=1e-11（分量1e-16 N）。callback按真实力及平台停止。目标尺度约O(1)，并不声称缩放本身改善Hessian条件数。必要时一次既有解析HVP Newton-CG残差恢复；只投影三个平移，线搜索检查真实Pi与未投影残差。\n\n平台预注册：上轮20 Pa PRECHECK_FREEZE没有显式平台阈值。本轮取既有wide静态工具的能量相对1e-8、位置1e-10 m；采用最后3个不同的接受状态，适应短Newton恢复，并额外报告最后10态供审查，不用重复评价填充平台。平台是必要条件，绝对力1e-13 N与相对1e-5同时保留。\n\n参考：[SciPy L-BFGS-B](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html)。程序/日志实际身份见SOURCE_AND_LIBRARY_IDENTITY及只读历史清单。\n'''
(O/'SOLVER_FAILURE_ROOT_CAUSE_AUDIT.md').write_text(audit)
plan=dict(main=dict(method='L-BFGS-B',coordinate_scale_m=1e-6,energy_scale_J=1e-11,maxiter=3000,maxfun=6000,ftol=0.,gtol=1e-11,maxcor=30,maxls=40,wall_s=240),recovery=dict(method='One bounded projected analytic-HVP Newton-CG continuation',maximum_invocations=1,maximum_Newton_iterations=10,maximum_total_HVP=1600,CG_maxiter_each=600,CG_rtol=1e-7,CG_atol_N=1e-19,wall_s=120,FD_step_m=1e-9,FD_relative_gate=1e-6,line_search_alphas=[1.,.5,.25,.125,.0625],acceptance='Unprojected L2 residual strictly decreases; Pi <= previous Pi + 100 eps max(|U|,|Pi|,1e-30). No duplicate coordinates accepted.',contact='Abort recovery for Ucontact >1e-24 J; HVP is existing noncontact tangent'),budgets=dict(total_solver_wall_s=360,peak_RSS_MiB=4096,output_MiB=1500,physical_load_cases=1,native_static_instances=1,native_coordinate_evaluations=4,formal_dynamic_fracture_cases=0),old_cost=dict(wall_s=24.896568,peak_RSS_MiB=233.238),scope='One main solve from zero A0, at most one same-2Pa recovery; no other optimizer/load/specimen. Numerical accepted states are not physical time.')
dump(O/'STATIC_SOLVER_RECOVERY_PLAN.json',plan)
freeze=dict(created_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),geometry_sha256=sha(O/'reference_topology.npz'),initial_state='Exact accepted A0 zero displacement; NOT saved 20Pa endpoint',p_Pa=2.,F_node_N=2*DX**2,F_face_N=800*2*DX**2,material=dict(K_Pa=K,G_Pa=G,rho_kg_m3=1000,s00=S00,alpha=0,dx_m=DX,delta_m=DELTA),Pi='Existing expanded fixed-m LPS U + existing conservative contact U - sum(Fext*u)',gates=dict(force_max_abs_N=1e-13,force_L2_relative=1e-5,max_positive_stretch=.03,plateau_accepted_distinct_states=3,energy_span_relative=1e-8,position_span_max_m=1e-10,contact_U_max_J=1e-24,native_force_gradient_relative=1e-9,reference_m_relative=1e-12,global_force_relative_to_face=1e-10,global_torque_relative_to_face_times_W=1e-10,topology_new_breaks=0),gauge='Remove only 3 uniform translations in optimization coordinates; evaluate unprojected forces',plan_sha256=sha(O/'STATIC_SOLVER_RECOVERY_PLAN.json'),prior_freeze_sha256=sha(OLD/'PRECHECK_FREEZE.json'),stopping='One bounded main + at most one recovery; then native certification/report, regardless of success. No fracture.')
dump(O/'LOW_LOAD_DESIGN_FREEZE.json',freeze)
dump(O/'PREPARATION_RESOURCES.json',dict(wall_s=time.perf_counter()-t,protected_files=len(protect),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
print('FROZEN',O,flush=True)
