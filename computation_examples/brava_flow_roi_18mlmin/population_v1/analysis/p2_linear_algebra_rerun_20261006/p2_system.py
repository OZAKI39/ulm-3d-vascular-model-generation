from launch import *
if (HERE/'P2_SYSTEM_ATTEMPT.json').exists():raise SystemExit('P2_SYSTEM_ALREADY_ATTEMPTED_NO_RETRY')
parity=read(HERE/'P1_LINEAR_SOLVER_PARITY.json');assert parity['status']=='PASS' and parity['candidate']['CPU_PC_hierarchy_verified']
assert read(TARGET/'P2_CURVED_MESH_QA.json')['status']=='PASS'
assert read(TARGET/'P2_PIPE_READINESS_GATE.json')['nine_implementation_prerequisites_pass']
for n,h in read(HERE/'INPUT_SHA256.json').items():assert sha(BASE/n)==h,n
for n,h in read(HERE/'INHERITED_EVIDENCE.json').items():assert sha(BASE/n)==h,n
budget=read(HERE/'CPU_SCHUR_GRAPH_BUDGET.json');m=managed().cgroup_memory();limit=int(Path('/sys/fs/cgroup/memory.max').read_text())/2**30
free=float(subprocess.check_output(['nvidia-smi','--query-gpu=memory.free','--format=csv,noheader,nounits'],text=True).strip())
# No CUDA sparse-matrix products in the preconditioner; observed P1 CPU hierarchy
# plus exact native source establishes this backend path. Allow large contingency
# above fine operator and complete FGMRES basis; runtime allocator guard is additional.
resident_GiB=budget['full_CSR_payload_bytes']/2**30+2*1037320*8*105/2**30+5
pre=dict(candidate='schur_cpu_pc',host_resource=m,cgroup_limit_GiB=limit,required_host_headroom_GiB=.2*limit,host_job_forecast_GiB=budget['planning_host_job_GiB'],host_safe=m['nonreclaimable_upper_GiB']+budget['planning_host_job_GiB']<=.8*limit,GPU_free_MiB=free,GPU_resident_forecast_GiB=resident_GiB,GPU_reserve_MiB=2048,GPU_safe=resident_GiB*1024+2048<=free,forecast_kind='Explicit array budget plus contingency, NOT measured total peak or mathematical bound',CUDA_PC_SpGEMM_workspace_bytes=0,CPU_hierarchy_P1_verified=True,allocation_guard=read(HERE/'CUDA_ADMISSION_GUARD.json'))
put('P2_SYSTEM_PRELAUNCH.json',pre)
if not (pre['host_safe'] and pre['GPU_safe']):
 put('P2_LINEAR_SYSTEM_RESOURCE_RESULT.json',dict(status='P2_LINEAR_SYSTEM_MEMORY_UNSAFE',setup='NOT_RUN_PREFLIGHT',solve='NOT_RUN',headroom_20percent_proven=False,prelaunch=pre));raise SystemExit(0)
put('P2_SYSTEM_ATTEMPT.json',dict(attempts=1,scope='one native assembled P2 linear system; zero completed CFD timesteps; no retry',prelaunch_sha256=sha(HERE/'P2_SYSTEM_PRELAUNCH.json'),binary_sha256=sha(HERE/'runtime/svmultiphysics_audit_v4')))
run('p2_system','single')
