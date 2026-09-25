from pathlib import Path
import json,hashlib,sys,csv,platform
import numpy as np
R=Path('/home/lzy/projects/compre_output/palabos_lammps_coupling_v0/20260915_224019')
REMOTE='/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019'
def read(p):return json.loads((R/p).read_text())
def save(p,o):(R/p).write_text(json.dumps(o,indent=2,ensure_ascii=False)+'\n')
def sha(p):return hashlib.sha256((R/p).read_bytes()).hexdigest()
a=read('LOCAL_INDEPENDENT_FINALIZER.json');remote=read('REMOTE_INDEPENDENT_FINALIZER.json');geo=read('LOCAL_EXACT_GEOMETRY_AUDIT.json');field=read('contracts/FROZEN_FLOW_FIELD_CONTRACT.json');contract=read('contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json');step=read('contracts/PARTICLE_TIMESTEP_CONTRACT.json');build=read('provenance/COUPLING_BUILD_PROVENANCE.json');sync=read('cases/kokkos_compatibility/HOST_DEVICE_SYNC_AUDIT.json');download=read('REMOTE_TO_WSL_INTEGRITY.json')
assert a['status']==remote['status']==geo['status']==download['status']=='PASS'
m={c['case']:c['metrics'] for c in a['cases']};A=m['case_a_uniform'];B=m['case_b_linear'];C=m['case_c_real'];E=m['case_e_mpi4'];G=m['kokkos_compatibility'];u=a['unit_audit']['metrics']
provenance={'Palabos_GPU_commit':field['Palabos_GPU_commit'],'Palabos_commit_evidence':'provenance/PBS_BSA_BUILD_PROVENANCE.json (original verified smoke provenance)','LAMMPS_release':'22Jul2025 Update 6','LAMMPS_commit':'9c5ab448c78a14fd534619622162ba418d6a1fb1','histogram_sha256':sha('provenance/frozen_sampler/input/FROZEN_SONOVUE_HISTOGRAM.csv'),'sampler_contract_sha256':sha('provenance/frozen_sampler/contracts/SONOVUE_SAMPLER_CONTRACT_V0.json'),'frozen_field':field,'PBS_BSA_numerics_sha256':sha('provenance/PBS_BSA_NUMERICS.json'),'coupling_source_SHA256':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((R/'src').glob('*')) if p.is_file()},'adapter_type':'EMBEDDED_CUSTOM_FIX','plugin_SHA256':'NOT_APPLICABLE: no dynamic plugin; custom fix sources and embedding binaries hashed separately','build':build,'local_python':sys.version,'local_numpy':np.__version__,'frozen_input_manifest_sha256':sha('FROZEN_INPUT_SHA256SUMS'),'numerical_evidence_manifest_sha256':sha('NUMERICAL_EVIDENCE_SHA256SUMS'),'local_baseline_audit':read('LOCAL_FROZEN_BASELINE_AUDIT.json'),'remote_baseline_audit':read('provenance/BASELINE_INTEGRITY_AFTER.json'),'remote_to_WSL':download}
save('PALABOS_LAMMPS_COUPLING_PROVENANCE.json',provenance)
rows=[]
for s in contract['cases']:
 d=m[s['name']];rows.append({'case':s['name'],'N':s['N'],'mpi_ranks':s['mpi_ranks'],'steps':s['steps'],'physical_time_s':s['steps']*s['dt_s'],'wall_s':d['wall_seconds'],'peak_rss_kib':d['peak_rss_kib'],'status':'PASS','max_velocity_normalized_error':d.get('max_velocity_normalized_error','NOT_APPLICABLE'),'max_position_normalized_error':d.get('max_position_normalized_error','NOT_APPLICABLE'),'invalid_queries':d['invalid_query_count'],'nonfinite':d['nonfinite_count'],'lost_atoms':d['lost_atoms']})
with (R/'CASE_VALIDATION_SUMMARY.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'README.md').write_text('''本目录是 PALABOS_LAMMPS_COUPLING_V0 的独立工程验证归档。

阅读 PALABOS_LAMMPS_COUPLING_V0_REPORT.md 和 FINAL_TERMINAL_SUMMARY.txt。
冻结物理与阈值：contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json。
真实场只有一份：fields/FROZEN_FLOW_FIELD_V0.h5；所有算例相对路径引用同一文件。

不运行任何求解器的独立复核：

```bash
sha256sum -c SHA256SUMS
/usr/bin/python3 -B src/finalize_palabos_lammps_coupling.py --root "$PWD" --output /tmp/coupling_reaudit.json
/usr/bin/python3 -B scripts/audit_geometry_local.py "$PWD"
```

Python reference 使用 NumPy、h5py、SciPy；精确 STL 后处理另需 VTK。
源快照身份复核需要原始 WSL PBS/BSA smoke 或其原 Vast 路径；归档没有重复复制完整原快照。
几何后处理命令会重新写本目录的两个几何审计文件；需要保持整包字节不变时请在独立副本运行。

重新编译只读 LAMMPS baseline 的外部适配器：

```bash
cmake -S . -B /path/to/new/build -G Ninja -DENGINE_ROOT=/path/to/validated/lammps/work -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/new/build -j 4
```

CMake 使用当前已验证平台的 CUDA 13.2 与 Kokkos 静态库布局；迁往其他环境需配置依赖路径。
不修改 LAMMPS core。CPU/GPU embedding 在启动时向现有 public fix_map 注册 frozen/flow/drag。
`binaries/` 是实际执行二进制的证据，不能假定跨机器直接可用。
运行使用 `HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0`；前者只排除 GL 显示拓扑探测。
`scripts/run_validation.py` 是本轮固定远端一次性执行记录：检查启动标记，拒绝覆盖或二次求解，不是可直接重复提交的模板。
任何新科研 case 应另建目录并明确授权；本阶段结束后停止。
''')
size=round(sum(p.stat().st_size for p in R.rglob('*') if p.is_file())/2**30,3)
summary={
'PALABOS_LAMMPS_COUPLING_V0':'PASS','PLATFORM':'VAST RTX4090','COUPLING_ARCHITECTURE':'FROZEN_FLOW_ONE_WAY','PALABOS_RUNTIME_DURING_PARTICLE_CASES':'NO','PALABOS_GPU_COMMIT':field['Palabos_GPU_commit'],'LAMMPS_RELEASE':'22Jul2025 Update 6','LAMMPS_COMMIT':'9c5ab448c78a14fd534619622162ba418d6a1fb1','SONOVUE_SAMPLER_STATUS':'PASS','SONOVUE_HISTOGRAM_SHA256':provenance['histogram_sha256'],
'FROZEN_FLOW_FIELD_FILE':'fields/FROZEN_FLOW_FIELD_V0.h5','FROZEN_FLOW_FIELD_SHA256':field['field_sha256'],'FLOW_FIELD_SOURCE':'PBS_BSA_5000_STEP_SMOKE','FLOW_FIELD_SCIENTIFIC_STATUS':'ENGINEERING_TRANSIENT_FIELD_ONLY','VELOCITY_VECTOR_SOURCE':'Velocity_m_s','INVALID_VELOCITY_MAGNITUDE_FIELD_USED':'NO','DX_M':field['dx_m'],'ORIGIN_M':field['origin_m'],'LATTICE_DIMS':'493 x 497 x 280','FLOW_FIELD_FORMAT':'FROZEN_FLOW_FIELD_V0','INTERPOLATION':'TRILINEAR','INVALID_QUERY_POLICY':'HARD_FAIL',
'COORDINATE_MAPPING_AUDIT':'PASS','MAX_COORDINATE_ERROR_M':u['max_coordinate_error_m'],'CASE_0A_UNIFORM_SAMPLER':'PASS','MAX_UNIFORM_VELOCITY_ERROR_M_S':u['uniform_reference_error_m_s'],'CASE_0B_LINEAR_INTERPOLATION':'PASS','MAX_LINEAR_INTERPOLATION_ERROR_M_S':u['linear_reference_error_m_s'],'CPP_PYTHON_INTERPOLATION':'PASS',
'FORCE_MODEL':'TECHNICAL_STOKES_DRAG_V0','MU_PA_S':.001,'TECHNICAL_PARTICLE_DENSITY_KG_M3':1000,'TECHNICAL_DENSITY_IS_SONOVUE_PHYSICS':'NO','MICROBUBBLE_TRANSPORT_MODEL':'PENDING','FORCE_UNIT_TEST':'PASS','MAX_FORCE_RELATIVE_ERROR':u['max_force_relative_error'],'PARTICLE_TIMESTEP_S':contract['cases'][0]['dt_s'],
'CASE_A_UNIFORM_ANALYTIC_COUPLING':'PASS','CASE_A_DIAMETER_UM':A['diameter_um'],'CASE_A_TAU_P_S':A['tau_p_s'],'CASE_A_MAX_VELOCITY_ERROR':A['max_velocity_normalized_error'],'CASE_A_MAX_POSITION_ERROR':A['max_position_normalized_error'],
'CASE_B_LINEAR_FIELD_TRAJECTORY':'PASS','CASE_B_MAX_POSITION_ERROR':B['max_position_normalized_error'],'CASE_B_MAX_VELOCITY_ERROR':B['max_velocity_normalized_error'],'REAL_FIELD_NODE_IDENTITY':'PASS','REAL_FIELD_CPP_PYTHON_INTERPOLATION':'PASS',
'CASE_C_REAL_PALABOS_FIELD':'PASS','CASE_C_DIAMETER_UM':C['diameter_um'],'CASE_C_STEPS':C['steps'],'CASE_C_PHYSICAL_TIME_S':C['physical_time_s'],'CASE_C_INVALID_QUERY_COUNT':C['invalid_query_count'],'CASE_C_MIN_WALL_CLEARANCE_UM':geo['min_exact_surface_clearance_m']*1e6,'CASE_C_CLEARANCE_DEFINITION':'sphere surface to closed STL = center-wall distance - radius',
'CASE_D_DIAMETER_SCALING':'PASS','D10_TEST_DIAMETER_UM':m['case_d_d10']['diameter_um'],'D50_TEST_DIAMETER_UM':m['case_d_d50']['diameter_um'],'D90_TEST_DIAMETER_UM':m['case_d_d90']['diameter_um'],
'CASE_E_MULTIPARTICLE_COUPLING':'PASS','CASE_E_PARTICLE_COUNT':1000,'MPI_COUPLING_COMPARISON':'PASS','MPI1_MPI4_MAX_POSITION_DIFFERENCE_M':E['mpi_max_position_difference_m'],'MPI1_MPI4_MAX_VELOCITY_DIFFERENCE_M_S':E['mpi_max_velocity_difference_m_s'],
'KOKKOS_COMPATIBILITY_SMOKE':'PASS','COUPLING_GPU_ACCELERATION':'PENDING','GPU_PERFORMANCE_READY':'NO','NONFINITE_COUNT':a['nonfinite_count'],'LOST_ATOMS':a['lost_atoms'],'INDEPENDENT_FINALIZER':'PASS','TOTAL_RESULT_GIB':size,'REMOTE_TO_WSL_INTEGRITY':'PASS',
'FORMAL_PURE_FLUID_BASELINE_MODIFIED':'NO','PBS_BSA_SMOKE_MODIFIED':'NO','LAMMPS_ENGINE_BASELINE_MODIFIED':'NO','SONOVUE_SAMPLER_MODIFIED':'NO','BUBBLE_BUBBLE_INTERACTION_MODEL':'PENDING','BUBBLE_BUBBLE_LUBRICATION':'PENDING','MICROBUBBLE_WALL_HYDRODYNAMICS':'PENDING','MICROBUBBLE_ADHESION':'PENDING','TWO_WAY_COUPLING':'OFF','RBC':'OFF','REMOTE_RESULT_DIR':REMOTE,'LOCAL_RESULT_DIR':str(R),'REPORT':str(R/'PALABOS_LAMMPS_COUPLING_V0_REPORT.md'),'NEXT_STAGE':'USER_REVIEW_BEFORE_BUBBLE_BUBBLE_INTERACTION_V0'}
save('FINAL_SUMMARY.json',summary);(R/'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(f'{k} = {v}' for k,v in summary.items())+'\n')
lines=['# PALABOS_LAMMPS_COUPLING_V0 核查报告','',
'本阶段 **PASS**：在冻结速度场与技术性 Stokes 阻力条件下，已验证单位/坐标映射、三线性插值、力注入、完整粒子推进、MPI 一致性和 Kokkos host-fix 兼容性。独立 Python 复核在 Vast 和下载后的 WSL 均通过。',
'', '这些证据只支持工程耦合链正确。真实场是 PBS/BSA 5000 步瞬态 smoke，尚未证明长期稳态；入口倍率在新介质中未重新校准，原 smoke 观察到瞬态回流。技术密度 1000 kg/m³ 和 Stokes 力仅供数值验证，不能据此声称 SonoVue 水动力学、真实输运、margination、壁面停留或分支选择已验证。',
'', '本轮没有执行 Palabos 求解，没有改变流体参数、STL、入口倍率、出口条件或原有工程源码。bubble-bubble、lubrication、wall hydrodynamics、adhesion 与最终 transport model 均为 PENDING；two-way coupling 与 RBC 为 OFF。',
'', '## 冻结输入与实现', '',
f'- Palabos GPU commit：`{field["Palabos_GPU_commit"]}`；沿用已核验 smoke 的构建 provenance。',
'- LAMMPS：22Jul2025 Update 6，commit `9c5ab448c78a14fd534619622162ba418d6a1fb1`。复用原 CPU/GPU 静态库，外部 embedding 程序注册 custom Fix，未修改 core 或 baseline binary。动态 plugin 未使用；共享库、Fix 源文件和可执行文件均有 SHA256。',
f'- SonoVue histogram SHA256：`{provenance["histogram_sha256"]}`。原 sampler 和 contract 逐文件核验，所有 population 保留原 seed、CSV、metadata 与 source IDs。',
f'- 唯一真实场：`{field["field_file"]}`，SHA256 `{field["field_sha256"]}`，{(R/field["field_file"]).stat().st_size:,} bytes，182,694 个流体节点。稀疏索引排序且唯一；完整 fluid lookup 区分 SOLID 与 MISSING。',
f'- 原始 `fields_5000.bin` SHA256：`{field["source_sha256"]}`。原快照保存 PackedNode 的三分量 LU 速度；按冻结 dx/dt 一次转换为 `Velocity_m_s`。独立 finalizer 对全部 ID 和向量重新比较。没有读取或替代使用无效 `VelocityMagnitude_m_s`。',
f'- `dx={field["dx_m"]:.17g} m`；origin `{field["origin_m"]} m`；dims 493×497×280。公式 `x_m=origin+(ix,iy,iz)*dx`；展平 `id=ix+nx*(iy+ny*iz)`。原 VTI CellData 几何中心对应 Palabos node，不另加半格偏移。',
'- synthetic uniform/linear/missing-corner 与真实场使用同一 HDF5 schema 和同一 C++ sampler 路径。无 nearest、clamp、silent extrapolation；八角必须全部合法流体节点。',
'- `F=3*pi*mu*d*(u_f-v_p)`，mu=0.001 Pa·s；SI 单位。每步在 POST_FORCE 内存注入 atom->f。独立核验 dump 中实际力，以及采样流速、直径与 force-evaluation velocity。',
'- 保留原 `fix nve/sphere`。速度依赖力在其标准半步速度阶段计算；轨迹同时保存 force_v、applied_F 和完整步状态的诊断 F，避免把两种速度时相混为一谈。',
f'- 冻结 particle dt={contract["cases"][0]["dt_s"]:.17g} s。支持下限给出 tau_min=3.125e-8 s、tau_max=1.53125e-6 s；采用 tau_min/500，满足 <=tau_min/50。该选择在任何正式数值输出产生前完成，未按结果改动。',
'- 冻结清单包含 86 个输入/代码文件，所有正式运行与后处理均重新验证。初次 supervisor 启动遇到 Unix socket 路径过长，求解器启动次数为 0；只缩短控制 socket 路径后启动。本轮 9 个算例各执行一次，没有数值重试。',
'', '## 插值、坐标和力硬门', '', '| 检查 | 最大误差 | 阈值 | 结果 |','|---|---:|---:|---|',
f'| 256 个原 VTI 节点坐标 | {u["max_coordinate_error_m"]:.9g} m | 1e-12 m | PASS |',
f'| 10,000 个 uniform 查询 | {u["uniform_reference_error_m_s"]:.9g} m/s | 1e-13 | PASS |',
f'| 10,000 个 linear 查询 vs 解析 | {u["linear_reference_error_m_s"]:.9g} m/s | 1e-12 | PASS |',
f'| linear C++ vs Python | {u["linear_cpp_python_m_s"]:.9g} m/s | 1e-12 | PASS |',
f'| 1,000 个真实节点 identity | {u["real_nodes_reference_error_m_s"]:.9g} m/s | 1e-12 | PASS |',
f'| 10,000 个真实 interior C++ vs Python | {u["real_interior_cpp_python_m_s"]:.9g} m/s | 1e-12 | PASS |',
f'| 10,000 组力函数相对误差 | {u["max_force_relative_error"]:.9g} | 1e-12 | PASS |',
f'| 其中 100 组零力绝对误差 | {u["max_zero_force_absolute_N"]:.9g} N | 1e-25 N | PASS |',
'', 'outside x/y/z、NaN、缺失角点、真实 solid 六类查询均返回明确 invalid 状态；正式 Fix 遇到 invalid 会中止，不使用 CLI invalid 行中的数值占位作为流速。',
'', '## 完整轨迹结果', '', '| 算例 | N / MPI | steps | 物理时长 (s) | 归一速度误差 | 归一位置误差 | 结果 |','|---|---:|---:|---:|---:|---:|---|']
for row in rows:
 def fmt(x):return f'{x:.9g}' if isinstance(x,(float,int)) else '—'
 lines.append(f'| {row["case"]} | {row["N"]} / {row["mpi_ranks"]} | {row["steps"]} | {row["physical_time_s"]:.9g} | {fmt(row["max_velocity_normalized_error"])} | {fmt(row["max_position_normalized_error"])} | PASS |')
lines += ['', 'A/D 使用整个 trajectory 与解析解比较：速度误差除以 |u_f-v0|，位置误差除以 |u_f-v0|·tau_p；所有初速为零。位置误差没有除以全局坐标或大 box 长度。A/D 门为 1e-3。B 使用 DOP853（rtol=1e-12，位置 atol=1e-20 m、速度 atol=1e-15 m/s），相同归一方式，门为 2e-3。',
f'',f'Case A 第一颗 seed=20260918 的直径为 {A["diameter_um"]:.15g} µm，tau={A["tau_p_s"]:.15g} s；末态 slip/initial slip={A["final_normalized_slip"]:.9g}，覆盖 >=5 tau。',
'', 'Case D 根据 seed=20260920 的 10,000 个原样本选择最近 d10/d50/d90：', '', '| 分位目标 | 实际直径 (µm) | 公式 tau (s) | 数值 1/e crossing (s) |','|---|---:|---:|---:|']
for q in ['d10','d50','d90']:
 d=m['case_d_'+q];lines.append(f'| {q} | {d["diameter_um"]:.15g} | {d["tau_p_s"]:.12g} | {d["measured_e_folding_time_s"]:.12g} |')
lines += ['', '三者解析轨迹均通过；公式 tau 与 d² 一致，测得的 1/e crossing 严格递增。crossing 使用相邻轨迹记录插值，未拟合或更改 drag 参数。',
f'',f'Case E 使用 seed=20260921 的 1000 粒子，同一初态分别执行 MPI1/MPI4；最终 ID 和直径逐项完全相同，最大位置差 {E["mpi_max_position_difference_m"]:.9g} m，最大速度差 {E["mpi_max_velocity_difference_m_s"]:.9g} m/s。初始周期中心间距和保守运动上界证明整个短程无接触；表面间隙下界 {E["noncontact_surface_gap_lower_bound_m"]*1e6:.9g} µm。没有粒子间物理作用。',
'', '## Case C 的几何与科学边界', '',
f'seed=20260919 的 100 个样本中，先按最接近冻结 d50 的规则选定 {C["diameter_um"]:.15g} µm 粒子，之后才做几何放置。初始 center-wall 距离 2.38282920197874 µm；要求 radius+3dx=1.58705430916387 µm。粒径、seed、STL、安全 margin 均未为几何适配而改变。',
f'',f'冻结短程 {C["steps"]} 步，{C["physical_time_s"]:.15g} s。总位移 {C["max_displacement_m"]:.12g} m，说明真实场链路实际驱动了粒子，但路径非常短，仅用于工程查询和推进检查。',
f'',f'每一步均检查有限值、完整八角查询和 Lipschitz 保守壁距下界。下载后又对全部 {geo["number_of_trajectory_points"]} 个记录用原闭合 STL 计算精确三角面距离及 enclosure：全部在腔内；最小 center-wall 距离 {geo["min_exact_center_wall_distance_m"]*1e6:.12g} µm；最小球表面到壁间隙 {geo["min_exact_surface_clearance_m"]*1e6:.12g} µm，大于 3dx={geo["required_surface_margin_m"]*1e6:.12g} µm。终端 WALL_CLEARANCE 指球表面间隙，不是中心到壁距离。',
'', '该安全规则是防止当前无壁面模型的工程算例进入未定义区域；它不是 wall force，不提供真实 bubble wall-interaction 结论。',
'', '## Kokkos 兼容性与资源证据', '',
f'同一 1000 粒子短例由 GPU/Kokkos binary 执行，非有限值、丢粒子、invalid query 均为 0。NVML 100 ms 采样记录的 utilization 峰值 {G["gpu_utilization_peak_percent"]:g}%，VRAM 峰值 {G["gpu_vram_peak_mib"]:g} MiB。采样可能遗漏短峰值，不据此评价性能。',
'', '保留原始 `kokkos_host_fix_trace.nsys-rep`、SQLite 和 CUDA API/memcpy 明细。trace 涵盖整次 invocation（含初始化和 run 0），不能把这些总量直接标成纯每步性能：',
'', '| 观测 | 次数 | 数据/累计时长 |','|---|---:|---:|']
for x in sync['CUDA_memcpy']:
 name={1:'H2D',2:'D2H'}.get(x['CUPTI_copyKind'],str(x['CUPTI_copyKind']));lines.append(f'| {name} | {x["calls"]} | {x["bytes"]:,} bytes；{x["seconds"]:.6g} s |')
for x in sync['CUDA_runtime_API']:
 if 'Synchronize' in x['name']:lines.append(f'| {x["name"]} | {x["calls"]} | {x["total_seconds"]:.6g} s API 时间 |')
lines += ['', '这些计数与保守 host Fix 的 ALL_MASK 同步路径一致。CUDA kernel 累计时间约 '+f'{sync["CUDA_kernel_summed_seconds"]:.6g}'+' s；API、copy、kernel 时间可能重叠，不能相加当总 wall time。没有改写 Kokkos kernel 或进行 GPU 优化；GPU_PERFORMANCE_READY=NO，COUPLING_GPU_ACCELERATION=PENDING。',
'', '## 完整性、大小与复核入口', '',
f'远端与 WSL 独立 finalizer 均为 PASS；全部数值记录中 NONFINITE_COUNT=0、LOST_ATOMS=0。原 LAMMPS 源码/二进制在前后 inventory 中一致，原 LAMMPS 与 PBS/BSA 结果清单远端再次通过；WSL 原 sampler、LAMMPS 结果、PBS/BSA 结果也逐文件通过。详情见 PALABOS_LAMMPS_COUPLING_PROVENANCE.json。',
f'',f'远端数值证据共 {download["checked_files"]} 个校验项，下载后 `sha256sum -c SHA256SUMS` 全部通过；NUMERICAL_EVIDENCE_SHA256SUMS 保留第一次数值封存，最终 SHA256SUMS 另覆盖本地审计与报告。结果总量约 {size:.3f} GiB（逻辑文件字节数，四舍五入），小于 2 GiB；真实 HDF5 仅一份。',
'', '- `contracts/`：物理、选择、时步、geometry safety 和误差门。','- `src/`：C++ library/Fix、Python reference 与独立 finalizer。','- `cases/`：原输入、初末 dump、全/稀疏轨迹、每步计数器、执行日志和资源 trace。','- `validation/`：原始查询、C++ 输出、原 VTI 对照、单元审计。','- `LOCAL_INDEPENDENT_FINALIZER.json`、`LOCAL_EXACT_GEOMETRY_AUDIT.json`：本地独立最终核查。','- `FINAL_TERMINAL_SUMMARY.txt`：完整用户要求的终端字段。','',f'远端：`{REMOTE}`。',f'本地：`{R}`。','', '本阶段到此停止。NEXT_STAGE=USER_REVIEW_BEFORE_BUBBLE_BUBBLE_INTERACTION_V0。']
(R/'PALABOS_LAMMPS_COUPLING_V0_REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({'status':'PASS','report':str(R/'PALABOS_LAMMPS_COUPLING_V0_REPORT.md'),'size_gib':size}))
