from common import *
import sys,time,resource,ctypes,csv,difflib
sys.path.insert(0,str(R))
from analysis.dual_track_v1.core_a import Energy,delta_y
T=time.perf_counter();d=np.load(O/'final_candidate.npz');rt=np.load(O/'reference_topology.npz');x=d['x'];X=d['X'];b=d['bonds'];alive=np.ascontiguousarray(d['alive'],dtype=np.uint8);e=Energy(dict(X=X,bonds=b));u=x-X
lib=ctypes.CDLL(str(R/'analysis/dual_track_v1/hvp_y.so'));hv=lib.hessian_vector;hv.restype=None;hv.argtypes=[ctypes.c_int]*2+[ctypes.c_void_p]*6+[ctypes.c_double]*3+[ctypes.c_void_p]*2
def H(v):
 v=np.ascontiguousarray(v);out=np.zeros_like(x);hv(len(x),len(b),*[a.ctypes.data for a in [x,e.b,alive,e.r,e.h,e.m]],V,K,G,v.ctypes.data,out.ctypes.data);return out
v=np.random.default_rng(928).normal(size=x.shape);v-=v.mean(0);v/=np.linalg.norm(v);h=1e-9;fd=-(e.eval(x+h*v,alive,True)[1]-e.eval(x-h*v,alive,True)[1])/(2*h);err=float(np.linalg.norm(fd-H(v))/np.linalg.norm(fd));rays={}
for name,vv in [('random_projected',v),('rotation_Y_tangent',np.cross(np.broadcast_to([0.,1.,0.],x.shape),x)),('saved_displacement',u)]:
 vv=vv-vv.mean(0);vv/=np.linalg.norm(vv);hh=H(vv);rays[name]=dict(Rayleigh_N_m=float(np.sum(vv*hh)),Hnorm_N_m=float(np.linalg.norm(hh)))
trans=[]
for j in range(3):
 v=np.zeros_like(x);v[:,j]=1/np.sqrt(len(x));trans.append(float(np.linalg.norm(H(v))))
dump(O/'READ_ONLY_TANGENT_AUDIT.json',dict(source_state='Final UNACCEPTED 2Pa candidate',new_solve=False,FD_step_m=h,HVP_FD_relative=err,translation_H_norms=trans,directional_Rayleigh=rays,scope='Read-only tangent diagnostics after bounded CG failure. Directional quotients are not extreme eigenvalues, a condition number, a preconditioner benchmark or an accepted elastic mode. Rotation-Y tangent remains unprojected. No additional CG solve.'))
mag=np.linalg.norm(u,axis=1);im=int(np.argmax(mag));q=float(u[rt['high'],2].mean()-u[rt['low'],2].mean());di=dict(status='UNACCEPTED_CANDIDATE_ONLY',mean_opening_m=q,max_displacement_m=float(mag[im]),max_displacement_atom_id=im+1,max_displacement_reference_m=X[im],max_displacement_current_m=x[im],mean_displacement_vector_m=u.mean(0),mean_displacement_magnitude_m=float(mag.mean()),top_mean_displacement_m=u[rt['high']].mean(0),bottom_mean_displacement_m=u[rt['low']].mean(0),certified_secant_compliance_m_N=None)
dump(O/'UNACCEPTED_DISPLACEMENT_DIAGNOSTIC.json',di)
main=json.loads((O/'MAIN_SOLVER_EXIT.json').read_text());rec=json.loads((O/'RECOVERY_SOLVER_EXIT.json').read_text());prior=np.load(O/'before_recovery_final_candidate.npz');unchanged=all(np.array_equal(d[k],prior[k]) for k in d.files)
with (O/'STATIC_SOLVER_TRACE.csv').open() as f:rows=list(csv.DictReader(f))
with (O/'function_evaluations.csv').open() as f:ev=list(csv.DictReader(f))
with (O/'before_recovery_function_evaluations.csv').open() as f:before=list(csv.DictReader(f))
statefiles=sorted((O/'accepted_optimizer_states').glob('state_*.npz'));assert len(statefiles)==len(rows)==989
assert np.array_equal(np.load(statefiles[-1])['x'],x)
assert main['evaluations']==len(before)
impl=dict(status='IMPLEMENTATION_CHECKS_COMPLETE_WITH_SCIENTIFIC_FAILURE_RETAINED',exact_original_topology_sha256=sha(O/'reference_topology.npz'),native_library_sha256=sha(B/'liblammps.so'),HVP_FD_relative=err,HVP_FD_gate=1e-6,HVP_FD_pass=err<1e-6,all_three_translation_tangents_exact_zero=max(trans)==0.,recovery_accepted_corrections=0,recovery_candidate_unchanged=unchanged,accepted_state_files=len(statefiles),trace_rows=len(rows),last_saved_state_matches_final=True,physical_traction_cases=1,main_solves=1,reserved_recoveries=rec['invocations'],formal_dynamic_cases=0,new_natural_breaks=0,raw_energy_force_evaluations=len(ev),main_energy_force_evaluations=len(before),recovery_energy_force_evaluations=len(ev)-len(before),function_count_semantics='Recorded energy/force evaluations include initial/callback evaluations; SciPy nfev unavailable after budget exception. Raw row count is authoritative.',logging_issue='Recovery seeded display evaluation counter from last accepted state 1034, whereas main had evaluated two later rejected trials (1036). Raw log rows are retained unchanged; total rows 1037 are the correct evaluation count. No accepted state or physics lost.',recovery_dispatch_correction='Original driver excluded any main exception, inadvertently excluding reserved recovery after MAIN_WALL_BUDGET. Standalone recovery honored the pre-frozen 120s, 600-CG/1600-HVP/10-Newton limits exactly. Old stop JSON, candidate, logs and trace preserved with before_recovery_ prefix. No second recovery.',plateau_scope='Last3 frozen before any 2Pa solve; last10 energy also reported and fails. Neither overrides force failure.',torque_scope='Fixed periodic Y cell does not admit arbitrary X/Z rigid rotations; raw-cell X/Z torque may include periodic boundary moments. Frozen full-vector check retained as conservative diagnostic. Y-compatible torque also exceeds its frozen absolute scale, and nodal absolute/relative residuals fail independently.',wall_s=time.perf_counter()-T)
dump(O/'IMPLEMENTATION_CHECKS.json',impl)
# Audit-only diff of the reused native adapter, no changes to historical code.
old=(R/'analysis/dual_track_v1/native_a.py').read_text();new=Path(__file__).with_name('native_adapter.py').read_text();(O/'native_adapter.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='original_native_a.py',tofile='task_native_adapter.py')))
print('READ_ONLY_CHECKS',json.dumps({'HVP_FD':err,'states':len(rows),'evals':len(ev),'same_candidate':unchanged,'directions':rays}),flush=True)
