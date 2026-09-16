from pathlib import Path
import json,csv,hashlib,sys,datetime
R=Path(sys.argv[1]);LOCAL='/home/lzy/projects/compre_output/lammps_particle_engine/20260915_214759'
def read(p,default=None):return json.loads(p.read_text()) if p.exists() else default
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ct=read(R/'contracts/MINIMAL_LAMMPS_PARTICLE_ENGINE_CONTRACT.json');b=read(R/'build_provenance/LAMMPS_BUILD_PROVENANCE.json');state=read(R/'STAGE_STATE.json',{});final=read(R/'validation/INDEPENDENT_FINALIZER.json',{})
transfer=read(R/'validation/WSL_DOWNLOAD_AUDIT.json',{'status':'PENDING_DOWNLOAD'})
cases={}
for sp in ct['cases']:
 key=sp['case'];cases[key]=read(R/'cases'/key/'INDEPENDENT_CASE_VALIDATION.json',{'status':'NOT_RUN','checks':{},'metrics':{}})
def status(k):return cases[k]['status']
def metric(k,n):return cases[k]['metrics'].get(n,'NOT_AVAILABLE')
a='case_a_single_ballistic';c='case_c_technical_contact';e='case_e_mpi_migration/cpu_mpi4';f='case_f_kokkos_gpu';cg='cpu_gpu_comparison/gpu'
getcmp=lambda key,sub,name:cases[key]['metrics'].get(sub,{}).get(name,'NOT_AVAILABLE')
s={'MINIMAL_LAMMPS_PARTICLE_ENGINE':'PASS' if state.get('state')=='PASS' and final.get('status')=='PASS' and transfer['status']=='PASS' else ('NUMERICAL_PASS_AWAITING_DOWNLOAD' if state.get('state')=='PASS' else 'FAIL'),
'PLATFORM':'VAST RTX4090','LAMMPS_RELEASE':b['release']['release_name'],'LAMMPS_GIT_COMMIT':b['release']['git_commit'],'LAMMPS_SOURCE_SHA256':b['source']['source_sha256'], 'CPU_BINARY_SHA256':b['builds']['cpu']['sha256'],'GPU_BINARY_SHA256':b['builds']['gpu']['sha256'],'LAMMPS_UNITS':'si','ATOM_STYLE':'sphere','LATBOLTZ':'OFF','SONOVUE_SAMPLER_STATUS':'PASS','SONOVUE_HISTOGRAM_SHA256':'2c9f783c6169421b06c057e16653878e7385139a179682c0648519e72a9f5198','DIAMETER_SUPPORT_UM':'0.75–5.25','DIAMETER_ROUNDTRIP':status('diameter_roundtrip'),'MAX_DIAMETER_ROUNDTRIP_ERROR_UM':metric('diameter_roundtrip','max_diameter_roundtrip_error_um'),'CASE_A_SINGLE_BALLISTIC':status(a),'CASE_A_MAX_POSITION_ERROR_M':metric(a,'max_ballistic_position_error_m'),'CASE_B_TWO_SEPARATED':status('case_b_two_separated'),'CASE_C_TECHNICAL_CONTACT':status(c),'CASE_C_FORCE_SYMMETRY_ERROR':metric(c,'force_symmetry_relative_error'),'BUBBLE_BUBBLE_PHYSICS_VALIDATED':'NO','CASE_D_POLYDISPERSE_1000':status('case_d_polydisperse_1000'),'CASE_E_MPI_MIGRATION':status(e),'MPI1_MPI4_MAX_POSITION_DIFFERENCE_M':getcmp(e,'mpi1_mpi4','max_position_difference_m'),'CASE_F_KOKKOS_GPU_SMOKE':status(f),'GPU_KERNEL_ACTIVITY':cases[f]['metrics'].get('gpu_kernel_activity',{}).get('status','NOT_AVAILABLE'),'GPU_VRAM_PEAK_MIB':metric(f,'gpu_vram_peak_mib'),'GPU_UTILIZATION':{'median_percent':metric(f,'gpu_utilization_median_percent'),'peak_percent':metric(f,'gpu_utilization_peak_percent')},'CPU_GPU_NO_CONTACT_COMPARISON':status(cg),'CPU_GPU_MAX_POSITION_DIFFERENCE_M':getcmp(cg,'cpu_gpu','max_position_difference_m'),'CPU_GPU_MAX_VELOCITY_DIFFERENCE_M_S':getcmp(cg,'cpu_gpu','max_velocity_difference_m_s'),'NONFINITE_COUNT':sum(x['metrics'].get('nonfinite_count',0) for x in cases.values()),'LOST_ATOMS':sum(x['metrics'].get('lost_atoms',0) for x in cases.values()),'INDEPENDENT_FINALIZER':final.get('status','NOT_COMPLETED'),'TOTAL_RESULT_GIB':0,'REMOTE_TO_WSL_INTEGRITY':transfer['status'],'FORMAL_PURE_FLUID_BASELINE_MODIFIED':'NO','PALABOS_USED':'NO','FLUID_DRAG':'OFF','REAL_VASCULAR_GEOMETRY':'OFF','ADHESION':'OFF','RBC':'OFF','TECHNICAL_CONTACT_ENGINE':status(c),'BUBBLE_BUBBLE_INTERACTION_MODEL':'PENDING','BUBBLE_BUBBLE_LUBRICATION':'PENDING','SONOVUE_CONTACT_MATERIAL_PROPERTIES':'NOT_DEFINED','MICROBUBBLE_TRANSPORT_MODEL':'PENDING','PALABOS_LAMMPS_COUPLING':'PENDING','MICROBUBBLE_WALL_MODEL':'PENDING','MICROBUBBLE_ADHESION':'PENDING','REMOTE_RESULT_DIR':str(R),'LOCAL_RESULT_DIR':LOCAL,'REPORT':LOCAL+'/MINIMAL_LAMMPS_PARTICLE_ENGINE_REPORT.md','NEXT_STAGE':'USER_REVIEW_BEFORE_PALABOS_LAMMPS_COUPLING'}
rows=[]
for sp in ct['cases']:
 key=sp['case'];m=cases[key]['metrics'];runtime=read(R/'cases'/key/'RUN_METRICS.json',{})
 rows.append({'case':key,'status':status(key),'N':sp['N'],'steps':sp['steps'],'seed':sp['seed'],'dt_s':sp['dt_s'],'mpi_ranks':sp['mpi_ranks'],'engine':sp['engine'],'wall_seconds':runtime.get('wall_seconds'),'solver_loop_seconds':m.get('solver_loop_seconds'),'solver_steps_per_s':m.get('solver_steps_per_s'),'end_to_end_steps_per_s':runtime.get('end_to_end_steps_per_s'),'process_tree_peak_rss_kib':m.get('process_tree_peak_rss_kib'),'lost_atoms':m.get('lost_atoms'),'profiled':runtime.get('profiled')})
with (R/'validation/RUNTIME_SUMMARY.csv').open('w') as fobj:
 wr=csv.DictWriter(fobj,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
cc=read(R/'contracts/TECHNICAL_CONTACT_TEST_CONTRACT.json');gpu=cases[f]['metrics'];em=cases[e]['metrics'];cm=cases[c]['metrics']
lines=['# MINIMAL_LAMMPS_PARTICLE_ENGINE_V0 核查报告','',f'阶段状态：**{s["MINIMAL_LAMMPS_PARTICLE_ENGINE"]}**。本阶段只验证有限尺寸刚性球的软件工程能力。真实 SonoVue 微泡间相互作用、lubrication、流体阻力、壁面、黏附和流固耦合均未验证；RBC 关闭。','',
'运行采用 Vast RTX4090（CC 8.9），独立 work/result 目录。未运行 Palabos，未修改正式纯流体基线、Stage4、PBS/BSA、Step3C 或 RBC 资产。', '',
f'LAMMPS 版本：{s["LAMMPS_RELEASE"]}；tag `{b["release"]["git_tag"]}`；commit `{s["LAMMPS_GIT_COMMIT"]}`。官方来源：[release]({b["release"]["release_url"]})。源码压缩包 SHA256：`{s["LAMMPS_SOURCE_SHA256"]}`。两套 binary 共用同一份源码。官方压缩包中 15013 个普通文件已逐一比对 SHA256，零差异。', '',
f'CPU binary SHA256：`{s["CPU_BINARY_SHA256"]}`。GPU binary SHA256：`{s["GPU_BINARY_SHA256"]}`。CMake 构建均启用 GRANULAR；GPU 另外启用 KOKKOS CUDA，`Kokkos_ARCH_ADA89=ON` 来自该版本源码的架构声明。LATBOLTZ 显式关闭。编译使用 gcc/g++ 13.3.0、CUDA nvcc 13.2.86、OpenMPI 4.1.6、CMake 3.28.3；并行编译限制 8。完整 flag、cache、帮助和日志在 `build_provenance/`。', '',
'单位为 SI；原始直径 µm 仅乘 `1e-6` 一次。`atom_style sphere` 保存 diameter、radius、mass、位置、速度、omega。密度 1000 kg/m³ 是技术测试值，**不是 SonoVue 微泡真实密度**；质量按球体积计算，半径固定。', '',
'粒径由原冻结 sampler 在原 WSL NumPy 1.26.4 环境生成；远端只读取结果。每组输入绑定 sampler source、contract、histogram、N、seed 和 population CSV 哈希。没有新采样逻辑、重新加权、裁尾或粒径调整。Case D 的 1000 颗为 ENGINE_TEST_ONLY；Case F 的 50000 颗为 GPU_ENGINE_SMOKE_ONLY，均不是实验浓度或正式科研 population。', '',
f'首个 run 0 前完成单位、阶段和技术接触合同封存。技术参数：Kn={cc["Kn_N_per_m"]} N/m，Kt=gamma_n=gamma_t=friction=0，dampflag=0；不启用 limit_damping。初始 overlap={cc["initial_overlap_m"]:.17g} m，dt={cc["dt_s"]} s，{cc["steps"]} 步。预检 dt·ω={cc["stability_preflight"]["dt_times_omega"]:.8g} < 0.05。CPU `gran/hooke/history` 与 GPU counterpart 的语法和 half-neighbor 要求均据冻结源码核查。该接触参数不表示 SonoVue 材料特性。', '',
'| 算例 | 状态 | N | steps | MPI / engine | wall s | loop steps/s | peak RSS KiB |','|---|---|---:|---:|---|---:|---:|---:|']
for row in rows:
 def fmt(x):return f'{x:.6g}' if isinstance(x,(int,float)) else str(x)
 lines.append(f'| {row["case"]} | {row["status"]} | {row["N"]} | {row["steps"]} | {row["mpi_ranks"]} / {row["engine"]} | {fmt(row["wall_seconds"])} | {fmt(row["solver_steps_per_s"])} | {fmt(row["process_tree_peak_rss_kib"])} |')
lines+=['', 'wall time 包括 MPI 启动、初始化、稀疏输出；F 还包含 profiler 开销，不能用于 CPU/GPU 性能比较。RSS 为每 0.1 s 采样的进程树合计峰值，短时峰值可能漏采。loop 时间来自对应实际 steps/count 的 LAMMPS 日志，逐项记录在 runtime CSV。', '',
f'直径 round-trip 最大误差：{s["MAX_DIAMETER_ROUNDTRIP_ERROR_UM"]} µm（门槛 1e-10）。A 解析位置最大误差：{s["CASE_A_MAX_POSITION_ERROR_M"]} m。B 检查静止位置、速度和零接触力。C 初始作用力对称相对误差：{s["CASE_C_FORCE_SYMMETRY_ERROR"]}，初始 overlap={cm.get("initial_overlap_m")} m，最终 overlap={cm.get("final_overlap_m")} m；另外检查初始 Hooke 力、排斥方向和能量尺度速度上界。', '',
f'E 使用相同输入、固定 4×1×1 分解，初始 ownership={em.get("initial_owners")}，最终 ownership={em.get("final_owners")}，实际 owner 改变粒子数={em.get("ownership_changes")}。证据来自 dump 的 proc 字段，而非只凭 mpirun rank 数。MPI1/MPI4 最大位置差={s["MPI1_MPI4_MAX_POSITION_DIFFERENCE_M"]} m；使用 unwrapped 坐标。', '',
f'CPU/GPU 无接触比较：位置最大差={s["CPU_GPU_MAX_POSITION_DIFFERENCE_M"]} m，速度最大差={s["CPU_GPU_MAX_VELOCITY_DIFFERENCE_M_S"]} m/s，ID 与 diameter 按 ID 排序精确匹配。没有使用接触算例强行要求 CPU/GPU 逐点一致。', '',
f'F 的 CUDA kernel 证据：{gpu.get("gpu_kernel_activity")}。原始 trace 在 `cases/case_f_kokkos_gpu/gpu_trace.nsys-rep`，SQLite 导出同目录，按 kernel 汇总在 `validation/GPU_KERNEL_SUMMARY.csv`。GPU VRAM 峰值={s["GPU_VRAM_PEAK_MIB"]} MiB；GPU-resident 样本利用率中位数={gpu.get("gpu_utilization_median_percent")}%，全部采样峰值={gpu.get("gpu_utilization_peak_percent")}%。排除 context 建立/退出端点后，驻留样本中央一半共 {gpu.get("gpu_steady_samples")} 个，VRAM 范围={gpu.get("gpu_vram_steady_range_mib")} MiB（冻结门槛 256 MiB）。该稳定性只覆盖本次短 smoke，不代表长程无内存增长。', '',
'日志中的 Dangerous builds 为 not checked：因为冻结配置采用 every 1、delay 0、check no，每一步无条件重建 neighbor list；本轮不把缺失计数记作零。独立补充核查会验证该日志声明和实际 200/1000 次重建次数。冻结 finalizer 中对应布尔值只表示没有发现非零计数，不是已测得危险重建次数为零。', '',
'独立 finalizer 不导入 converter，也不以 LAMMPS 的成功退出代替数值判断。它重新读取 population、实际 data、initial/final dumps、log、MPI ownership 和 CUDA activity，逐项核验。采样输出与 thermo 日志的 nonfinite 检查、particle count/ID 和 lost error 联合使用；本阶段没有高频 full dump。', '',
f'独立数值核查={s["INDEPENDENT_FINALIZER"]}，记录到 nonfinite={s["NONFINITE_COUNT"]}，lost atoms={s["LOST_ATOMS"]}。远端到 WSL 哈希核验={s["REMOTE_TO_WSL_INTEGRITY"]}。下载只包含紧凑证据与代码，不含源码包、build tree、checkpoint 或 restart。', '',
'构建后的帮助启动曾在 hwloc GL 显示扫描中等待：strace 显示尝试 X0–X6 并在 localhost:6006 的 X11 握手阻塞。移除 DISPLAY 的排查无效（该环境本身未设 DISPLAY）；仅在任务子进程设置 HWLOC_COMPONENTS=-gl 后，MPI1 GPU 帮助在 0.95 s 完成。此设置只关闭显示拓扑探测，不关闭 CUDA。原始失败/超时排查证据保留，所有排查均为零 timestep，未重新编译 binary。依据：[hwloc 官方组件选择说明](https://www.open-mpi.org/projects/hwloc/doc/v2.0.4/a00324.php)。', '',
'预检期间修正过与冻结版本打印格式不匹配的 KOKKOS 日志识别器，发生于首次正式 run 前，保留修改前后哈希；未改数值门槛。不存在为通过结果而调整接触参数、粒径或容差的操作。', '',
'TECHNICAL_CONTACT_ENGINE = '+s['TECHNICAL_CONTACT_ENGINE']+'；BUBBLE_BUBBLE_PHYSICS_VALIDATED = NO；BUBBLE_BUBBLE_INTERACTION_MODEL = PENDING；BUBBLE_BUBBLE_LUBRICATION = PENDING；SONOVUE_CONTACT_MATERIAL_PROPERTIES = NOT_DEFINED。', '',
'MICROBUBBLE_TRANSPORT_MODEL、PALABOS_LAMMPS_COUPLING、HYDRODYNAMIC_DRAG、MICROBUBBLE_WALL_MODEL、MICROBUBBLE_ADHESION 均为 PENDING。本阶段完成后停止，等待用户审核；没有自动进入耦合开发。', '',f'远端结果：`{R}`。本地结果：`{LOCAL}`。']
if state.get('state')!='PASS':lines+=['', '**失败/中止记录：** `'+json.dumps(state,ensure_ascii=False)+'`。未运行的依赖算例保持 NOT_RUN。']
(R/'MINIMAL_LAMMPS_PARTICLE_ENGINE_REPORT.md').write_text('\n'.join(lines)+'\n')
for _ in range(3):
 s['TOTAL_RESULT_GIB']=round(sum(p.stat().st_size for p in R.rglob('*') if p.is_file())/1024**3,9)
 (R/'FINAL_SUMMARY.json').write_text(json.dumps(s,indent=2,ensure_ascii=False)+'\n')
 (R/'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(k+' = '+(json.dumps(v,ensure_ascii=False) if isinstance(v,dict) else str(v)) for k,v in s.items())+'\n')
assert s['TOTAL_RESULT_GIB']<1
print(json.dumps({'stage':s['MINIMAL_LAMMPS_PARTICLE_ENGINE'],'gib':s['TOTAL_RESULT_GIB']},ensure_ascii=False))
