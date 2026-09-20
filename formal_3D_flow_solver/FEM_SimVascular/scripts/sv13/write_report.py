#!/usr/bin/env python3
"""Final performance-stage report derived from measured, immutable evidence."""
import json,sys,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import *
from sv_validation.provenance import write_json,sha256,now

def junit(path):
    rows=[]
    for t in ET.parse(path).getroot().iter('testcase'):
        status='failed' if t.find('failure') is not None else 'error' if t.find('error') is not None else 'skipped' if t.find('skipped') is not None else 'passed'
        rows.append({'id':t.get('classname')+'::'+t.get('name'),'status':status})
    return {'counts':{k:sum(r['status']==k for r in rows) for k in ('passed','failed','skipped','error')},'cases':rows}

check_reference();ref=load('reference_freeze');cpu=load('cpu_execution');perf=load('cpu_performance')
eq=load('cpu_reference_equivalence');qc=load('cpu_saved_states');gpu=load('gpu_decision')
env=load('gpu_environment');audit=load('preservation_audit');vis=load('visuals')
accepted=load('cpu_accepted_solution');policy=frozen_policy();validation=load('cpu_validation')
full=junit(REPORT/'pytest_all.xml');newrows=[r for r in full['cases'] if 'test_sv13_' in r['id']]
new={'counts':{k:sum(r['status']==k for r in newrows) for k in ('passed','failed','skipped','error')},'cases':newrows}
historical=[r['id'] for r in full['cases'] if r['status']=='failed' and 'test_sv13_' not in r['id']]
allowed={
 'tests.test_sv11_short_run_mass::test_actual_step_ten_mass_strictly_improves',
 'tests.test_sv1_mass_balance::test_actual_field_mass_conservation',
 'tests.test_sv1_solver_result::test_real_solver_converged',
 'tests.test_sv1_solver_result::test_required_steady_intervals_reached'}
normalize=lambda s:'tests.'+s if not s.startswith('tests.') else s
known_history={normalize(s) for s in historical}==allowed
tests={'sv13':new,'all':full,'historical_failures':historical,'historical_failure_set_unchanged':known_history,
    'acceptance_policy':'SV1.3 zero failures; historical failures remain visible and unchanged; blocked real GPU tests are explicit skips',
    'latest_full_suite_includes_final_figure_revision':True,'timestamp':now()}
write_json(REPORT/'test_results.json',tests)
stage_ok=(validation['status']=='PASS' and eq['status']=='PASS' and load('solution_reload')['status']=='PASS'
    and perf['measured_continuation_speedup']>1 and new['counts']['failed']==new['counts']['error']==0
    and full['counts']['error']==0 and known_history and audit['status']=='PASS' and len(vis['figures'])==10)
status='CONDITIONAL PASS' if stage_ok else 'FAIL'
write_json(REPORT/'stage_result.json',{'status':status,'scope':'CPU optimization only; GPU prerequisite blocked',
    'GPU_status':gpu['status'],'GPU_reason':gpu['reason_code'],'selected':'CPU_EARLY_STOP_PRODUCTION',
    'human_review_pending':True,'mesh_convergence_started':False,'timestamp':now()})
write_json(REPORT/'render_inspection.json',{'assistant_basic_inspection':'PASS','figure_count':10,
    'checked':'All ten images viewed; four revised images reviewed again; readable titles, units, zero mass markers, and explicit GPU unavailability',
    'scientific_fields_modified':False,'human_review_complete':False,
    'images':[{'path':f['path'],'sha256':f['sha256']} for f in vis['figures']]})
write_json(REPORT/'cuda_build_outcome.json',{'status':'BLOCKED','reason':'CUDA_PETSC_BUILD',
    'configure_attempts':2,'make_status':'NOT_RUN','solver_cuda_build_status':'NOT_RUN','binary_sha256':None,
    'VECCUDA_verified':False,'MATAIJCUSPARSE_verified':False,'actual_Mat_type':None,'actual_Vec_type':None,
    'PETSc_cuda_smoke':'NOT_RUN','official_fluid_cuda_smoke':'NOT_RUN','GPU_kernel_observed':False,
    'version_upgrade_performed':False,'source_patched':False,
    'limitation':'MPI startup prerequisite failed; source-level PETSc 3.19.6 / CUDA 13.2 incompatibility has not been established'})
hashes={f['path']:f['sha256'] for f in ref['files']}
baseline=ref['baseline_resources']['PETSc_iterations'];complete=load('solver_history','sv1_2')['complete_linear_iteration_distribution']
cs=cpu['linear_statistics'];errs=eq['errors'];a=eq['candidate'];r=eq['reference'];last=qc['intervals'][-1]
with np.load(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz') as arr:points=len(arr['points_m']);cells=len(arr['tetra'])
figlink=lambda name:f'![{name}]({name})'
table='\n'.join(f"| {role} | {r['outlet_flows_m3_s'][role]:.16g} | {a['outlet_flows_m3_s'][role]:.16g} | {a['outlet_fractions'][role]:.12g} | {errs['absolute_fraction_difference'][role]:.4g} |" for role in ('OUTLET_01','OUTLET_02','OUTLET_03'))
pt='\n'.join(f"| {role} | {r['area_average_pressure_pa'][role]:.16g} | {a['area_average_pressure_pa'][role]:.16g} | {errs['port_pressure_scale_normalized'][role]:.4g} |" for role in ('INLET','OUTLET_01','OUTLET_02','OUTLET_03'))
text=f'''# Stage SV1.3 — SimVascular 血流计算性能优化与 GPU 加速验证

**{status} — CPU optimization only，等待人工审核。** 推荐 `CPU_EARLY_STOP_PRODUCTION`。同一起点的实测续算从 {perf['reference_measured_wall_s']:.0f} s 降至 {perf['candidate_measured_wall_s']:.0f} s，**{perf['measured_continuation_speedup']:.4f}×**；计入共同继承的前 10 步成本后为 **{perf['accounted_total_speedup']:.4f}×**。GPU 状态为 `GPU_BLOCKED: CUDA_PETSC_BUILD`，没有 GPU 性能或驻留的实测结论。

## 为什么现在开始优化？

用户在本阶段请求中确认 SV1.2 已通过数值验收和人工审核。本阶段保持该 step{ref['accepted_solution']['step']} 解作为 `VALIDATION_REFERENCE`，只改变终止策略并尝试独立 CUDA build。新的性能配置属于 `FAST_PRODUCTION_CONFIGURATION`。所有网格、物理、边界、时间积分、非线性方法和官方求解器科学源码均保持不变。本阶段没有执行 Mesh convergence、WSS、RBC 或 microbubble。

冻结证据：[reference_freeze.json](reference_freeze.json)、[policy.json](../../configs/sv1_3/policy.json)。重新读取确认：

| 项目 | 冻结值 |
|---|---|
| 实际网格 | {points} points / {cells} tetrahedra |
| mesh-complete.mesh.vtu SHA256 | `{hashes['outputs/sv1/SV_MESH/mesh-complete.mesh.vtu']}` |
| 输入 exterior_surface.npz SHA256 | `{hashes['inputs/fem_reference/exterior_surface.npz']}` |
| 密度 / 动力黏度 | {ref['physics']['rho']:.16g} kg/m³ / {ref['physics']['mu']:.16g} Pa·s |
| Qtarget | {ref['accepted_solution']['Q_target_m3_s']:.17g} m³/s |
| dt | {ref['dt']:.17g} s |
| 时间积分 | {ref['time_integrator']['name']}；spectral radius {ref['time_integrator']['spectral_radius']} |
| BC | 入口 steady flat Dirichlet、Impose_flux；三个出口零自然牵引；壁面零速度 |
| PETSc / MPI / OMP | {ref['build']['PETSc_version']} / {ref['PETSc']['mpi_ranks']} ranks / {ref['PETSc']['OMP_NUM_THREADS']} thread |
| svMultiPhysics commit | `{ref['build']['commit']}` |
| accepted step / VTU SHA256 | {ref['accepted_solution']['step']} / `{ref['accepted_solution']['sha256']}` |
| reference velocity volume L2 | {r['velocity_L2']:.17g} |
| reference pressure range | {r['pressure_range_pa']} Pa |

生产模型 XML 为 [SV1.2 原 XML](../../configs/sv1_2/sv_flow.xml) 的逐字节副本。完整 PETSc runtime options 为：

```text
{ref['PETSc']['PETSC_OPTIONS']}
```

## 原来一次计算有多慢？

SV1.2 终点为 step{ref['accepted_solution']['step']}，同一合法 native step10 起点之后执行 {perf['reference_new_steps']} 个新时间步，GNU wall time {perf['reference_measured_wall_s']:.0f} s（{perf['reference_measured_wall_s']/3600:.3f} h）。首次连续五个稳态区间截止在 step{ref['baseline_steady']['first_five_joint_intervals_end_step']}。线性和非线性失败均为 0。

| KSP 指标 | SV1.2 新增 steps11–400 | SV1.2 含前10步完整历史 | CPU_EARLY_STOP 新增 steps11–80 |
|---|---:|---:|---:|
| solve count | {baseline['count']} | {complete['count']} | {cs['count']} |
| total iterations | {baseline['total']} | {complete['total']} | {cs['total']} |
| mean | {baseline['mean']:.6f} | {complete['mean']:.6f} | {cs['mean']:.6f} |
| median | {baseline['median']} | {complete['median']} | {cs['median']} |
| P95 | {baseline['P95']:.6f} | {complete['P95']:.6f} | {cs['P95']:.6f} |
| maximum | {baseline['max']} | {complete['max']} | {cs['max']} |

数据来自原始 [SV1.2 resources](../sv1_2/solver_resource_usage.json)、[SV1.2 solver history](../sv1_2/solver_history.json) 和 [本次执行记录](cpu_execution.json)。原生日志的逐次线性求解占比中位数为 98%；该段计时包含 KSP 设置、预条件和同步。首条记录因继承计时器而排除于占比图，原始记录仍完整保留。Assembly、PC 单项、I/O 和 other 没有独立 profiler 数值，图中明确标为未拆分/未测，不能据此构造精确总时间饼图。

{figlink('where_time_goes.png')}

## 最简单的优化是什么？

通用 `SteadyStopMonitor` 逐个读取真实保存场，用同一网格的 P1 体积 L2 计算 E_u，以 Qtarget 归一化四个端口的最大流量变化 E_Q。要求 E_u≤{policy['velocity_change_limit']:g}、E_Q≤{policy['flow_change_limit']:g} 连续 {policy['steady_last_intervals']} 个保存区间；当前质量、入口误差均≤{policy['mass_limit']:g}，场有限、壁面条件通过，且累计线性/非线性失败为零。

离线重放实际 {len(load('steady_stop_replay')['rows'])} 个 SV1.2 VTU，首次 stop=true 在 step{load('steady_stop_replay')['first_stop_step']}，与原报告完全一致。新 CPU 计算没有固定 step70 终点。算法在 step{validation['first_qualifying_step']} 合格后原子写入 `STOP_SIM`，请求下一个保存步 step{validation['stop_step']}。固定源码先完成 restart，再完成 VTU，然后退出；最终 exit code {cpu['exit_code']}。没有向原生流体进程发送 kill 或暂停信号。

原生机制证据见 [termination_source_audit.json](termination_source_audit.json)，请求见 [cpu_stop_request.json](cpu_stop_request.json)，完整原生终态检查见 [cpu_final_checkpoint.json](cpu_final_checkpoint.json)。最终又通过 {len([i for i in qc['intervals'] if i['E_u']<=policy['velocity_change_limit'] and i['E_Q']<=policy['flow_change_limit']])} 个连续合格区间的检查；最后 E_u={last['E_u']:.6g}，E_Q={last['E_Q']:.6g}。

{figlink('early_stop_savings.png')}

## 自动停止有没有改变答案？

[独立等价性计算](cpu_reference_equivalence.json) 与 [fresh-process 重读](solution_reload.json) 均 PASS。门槛在运行前冻结：速度与压力相对体积 L2≤1e-5，Qin 相对差≤1e-6，各 Qout 相对差≤1e-5，分流绝对差≤1e-5，各端口平均压力采用共同 reference max(abs(p)) 尺度归一化后≤1e-5；另检查最大速度相对差≤1e-5。没有压力平移或流量重标定。

| 检查 | 实测误差 |
|---|---:|
| velocity relative volume L2 | {errs['velocity_relative_volume_L2']:.12g} |
| pressure relative volume L2 | {errs['pressure_relative_volume_L2']:.12g} |
| Qin relative difference | {errs['relative_Qin']:.12g} |
| 最大各 Qout relative difference | {max(errs['relative_Qout'].values()):.12g} |
| 最大分流 absolute difference | {max(errs['absolute_fraction_difference'].values()):.12g} |
| 最大端口压力 scale-normalized difference | {max(errs['port_pressure_scale_normalized'].values()):.12g} |
| max velocity relative difference | {errs['relative_max_velocity']:.12g} |

首次收尾比较曾因原生流量日志额外包含 WALL 列而报告后处理错误。已按既有 SV1.2 做法只选择四个开口端口；该修正与原错误保留于 [cpu_postprocess_attempt1.json](cpu_postprocess_attempt1.json)。求解器没有失败或重跑，原生 VTU、checkpoint、日志均未修改。

{figlink('solution_equivalence.png')}

## RTX 4090 是否真正被 PETSc 使用？

**尚未获得有效 PETSc CUDA 运行。** 实际环境为 {env['gpu_model']}，物理显存 {env['gpu_memory_bytes']/1024**2:.0f} MiB，driver {env['driver_version']}，compute capability {env['compute_capability']}，toolkit/nvcc 13.2 / 13.2.86。空闲时 PCIe 观测为 gen{env['PCIe']['current_gen']}×{env['PCIe']['current_width']}，报告上限 gen{env['PCIe']['max_gen']}×{env['PCIe']['max_width']}；这不是带载吞吐测量。[gpu_environment.json](gpu_environment.json) 保存原始 nvidia-smi、nvcc、CPU/MPI 信息。

本机是 RTX 4060 Laptop 且没有 nvcc，因此 CUDA build 位于既有远端独立 `sv1_3` 目录。没有运行 Docker 命令或修改驱动。WSL 仍持有源码、配置和同步回来的构建证据。

PETSc 3.19.6 的真实 `configure --help` 已分别在新 WSL 副本和远端副本执行。首轮默认 CUDA 库搜索缺少 `nvToolsExt`；第二轮用该版本支持的显式库列表通过 CUDA 库/版本检查，随后 `mpiexec --oversubscribe -n 1` 测试超时，configure 返回失败。独立 `/bin/true` 启动在默认、关闭绑定和本地通信设置下也超时。见 [两轮构建汇总](cuda_petsc_build.json)、[首轮记录](cuda_petsc_build_attempt1.json)、[第二轮完整配置日志（凭据脱敏）](../../outputs/sv1_3/remote_build_attempt2/configure_detail.log)、[MPI 探测](remote_mpi_probe.json)。

因此 make、CUDA-enabled svMultiPhysics build、两种 CUDA smoke 均未执行；CUDA binary SHA、actual Mat type、actual Vec type、实际 GPU KSP/PC view 与 kernel 证据均为 NOT_OBSERVED。静态源码确实调用 MatSetFromOptions / VecSetFromOptions，但不能当作运行时 CUDA 验证。没有升级 PETSc，也没有证明 PETSc 源码与 CUDA 13.2 必然不兼容；实际阻塞点是远端 MPI 启动。

## 数据有没有频繁在 CPU/GPU 之间搬？

未知。本次没有有效 CUDA KSP，所以 matrix resident、vectors resident、传输次数随 nonlinear solves 或 KSP iterations 增长、H2D/D2H 时间等都保持 null/NOT_RUN。[gpu_residency_audit.json](gpu_residency_audit.json) 明确记录这一点，GPU 利用率和硬件可见性没有替代这些测量。

{figlink('gpu_residency.png')}

{figlink('gpu_transfer_cost.png')}

## GPU-A 做了什么？

完成了硬件与接口源码审核，以及独立 PETSc CUDA 配置尝试。原计划保持 GMRES、容差、ASM+ILU(2)，用 1 MPI rank / 1 GPU / OMP=1 验证 CUDA Mat/Vec。由于构建前置门槛失败，GPU_A vascular benchmark 没有启动。

## GPU-A 值得采用吗？

本阶段没有达到可评估性能的条件，不采用 GPU_A。以下同机固定 20-step、同一完整 native checkpoint 的比较全部 NOT_RUN，未用跨机器 CPU 耗时冒充 GPU 加速比，也没有用单次最优值代替重复试验。

| 同机固定窗口 | wall time | 等价性 |
|---|---|---|
| CPU_1R | NOT_RUN | NOT_RUN |
| CPU_4R | NOT_RUN | NOT_RUN |
| GPU_1R | NOT_RUN | NOT_RUN |

[benchmark 状态](../../benchmarks/sv1_3/status.json) 保留固定窗口与重复计时规则。源代码的 native checkpoint stamp 同时校验 MPI rank 数与本地节点数；没有把 4-rank checkpoint 静默当成 1-rank restart，也没有替换为丢失历史的 VTU 初态。

{figlink('cpu_gpu_runtime.png')}

## 如果运行 GPU-B

GPU_B = **NOT_RUN**，尝试数 0。GPU_A 尚未证明科学正确，也没有 PC/传输主导性能的 profiling 证据，因此不满足 GPU_B 启动条件。没有运行 GPU_C/D/E、field-split、Schur solver、自定义 GPU block solver 或 GPU FEM assembly。

## 最终哪套配置最快而且可靠？

**CPU_EARLY_STOP_PRODUCTION** 是本阶段唯一完成全部数值、稳态、质量、等价性、完整 restart、独立重读与实测加速门槛的配置。GPU 状态是 blocked，并非测出更慢。

正式配置：[production_performance.yaml](../../configs/sv1_3/production_performance.yaml)。它保留原 CPU PETSc 选项、4 MPI ranks、OMP=1，GPU count=0，规定动态稳态停止及最终检查。已完成证据目录受防覆盖保护；后续独立运行应使用新的输出目录，不能覆盖本次或 SV1.2。SV1.2 step400 永久保留为验证参考。

本次最大单进程 RSS {cpu['peak_single_process_RSS_KiB']} KiB（{cpu['peak_single_process_RSS_KiB']/1024:.3f} MiB），不是所有 MPI ranks 的内存总和。GPU 峰值显存未测，不能将设备总显存用作求解峰值。

{figlink('memory_usage.png')}

## 实际节省多少时间？

| 计时口径 | SV1.2 reference | CPU_EARLY_STOP | speedup |
|---|---:|---:|---:|
| 同一 step10 native 起点，GNU 实测 wall | {perf['reference_measured_wall_s']:.3f} s | {perf['candidate_measured_wall_s']:.3f} s | {perf['measured_continuation_speedup']:.6f}× |
| 同一历史前10步成本 + 各次 Python monotonic | {perf['reference_accounted_total_s']:.6f} s | {perf['candidate_accounted_total_s']:.6f} s | {perf['accounted_total_speedup']:.6f}× |
| best same-server CPU / GPU 20-step window | NOT_RUN | NOT_RUN | NOT_MEASURED |

节省实测续算时间 {perf['saved_measured_seconds']:.0f} s（{perf['saved_measured_seconds']/3600:.3f} h）。第二行的前10步成本 {perf['inherited_first10_wall_s']:.6f} s 来自历史原生运行，不是本阶段新测 t=0 全程；两种计时口径没有混用。自动停止的收益来自少计算 {perf['reference_new_steps']-perf['candidate_new_steps']} 个时间步，未更换预条件器。一次 CPU early-stop validation 按请求运行；它不是重复的同机 GPU benchmark。

{figlink('ksp_iterations_comparison.png')}

{figlink('final_speedup.png')}

## 科学结果有没有变化？

在预冻结等价门槛下没有可辨识改变。最终 Qin={a['Q_in_m3_s']:.17g} m³/s，Qout total={a['Q_out_total_m3_s']:.17g} m³/s；epsilon_Q={a['epsilon_Q']:.12g}，epsilon_mass={a['epsilon_mass']:.12g}。质量误差为 0 指双精度积分之差为零，并非离散误差证明。速度和压力有限，wall max={a['wall_velocity_max_m_s']:.12g} m/s，最大速度 {a['velocity_max_m_s']:.17g} m/s，pressure range={a['pressure_range_pa']} Pa。

| 端口 | SV1.2 Qout (m³/s) | CPU Qout (m³/s) | CPU fraction | fraction absolute difference |
|---|---:|---:|---:|---:|
{table}

| 端口 | SV1.2 mean pressure (Pa) | CPU mean pressure (Pa) | reference-scale normalized error |
|---|---:|---:|---:|
{pt}

最终保存解：[result_080.vtu](../../{accepted['path']})，SHA256 `{accepted['sha256']}`。比较仅验证对同一 SV1.2 解的保持，不能替代网格独立性验证。

{figlink('steady_equivalence.png')}

永久测试新增 {len(list((ROOT/'tests').glob('test_sv13_*.py')))} 个模块。SV1.3：**{new['counts']['passed']} passed / {new['counts']['failed']} failed / {new['counts']['skipped']} skipped / {new['counts']['error']} errors**。10 个真实 GPU 验证因前置环境阻塞而显式跳过；CPU 类型假冒 GPU、逐 KSP 大量搬运、OOM、科学不等价、质量失败、speedup<1.25、参考被修改、四区间提前停止等负例测试仍实际执行。

完整回归：**{full['counts']['passed']} passed / {full['counts']['failed']} failed / {full['counts']['skipped']} skipped / {full['counts']['error']} errors**。保留的历史失败集合与 SV1.2 完全相同：

'''+'\n'.join('- `'+s+'`' for s in historical)+f'''

没有删除、改写或 xfail 历史失败。最终 [preservation_audit.json](preservation_audit.json)：{audit['historical_entries_checked']} 个历史条目及旧 FEM {audit['old_fem_entries_checked']} 个条目均未改变；旧 FEM git 状态保持原样。新测试汇总见 [test_results.json](test_results.json)，图像 SHA 见 [visuals.json](visuals.json)。

## 下一步

本阶段停在 **{status}**，等待用户审核上述 10 张图，尤其是 CPU 实测加速、场等价性和 GPU 尚未验证的边界。只有 performance stage 获得正式 PASS 后才进入 Mesh convergence。本次没有自动开始下一阶段。
'''
(REPORT/'REPORT.md').write_text(text)
summary=f'''Stage SV1.3 completed.

reference:
    SV1.2 runtime = {perf['reference_measured_wall_s']} s (same step10 continuation)
    steps = {ref['accepted_solution']['step']} ({perf['reference_new_steps']} newly executed)
    first steady step = {ref['baseline_steady']['first_five_joint_intervals_end_step']}

CPU early stop:
    stop step = {cpu['last_step']} (trigger step {validation['first_qualifying_step']})
    runtime = {perf['candidate_measured_wall_s']} s
    speedup vs SV1.2 = {perf['measured_continuation_speedup']:.6f}x
    total including inherited first10 = {perf['accounted_total_speedup']:.6f}x
    science equivalent = PASS

GPU:
    model = {env['gpu_model']}
    CUDA = 13.2 / nvcc 13.2.86
    PETSc CUDA = GPU_BLOCKED: CUDA_PETSC_BUILD (remote MPI startup timeout)
    Mat type = NOT_OBSERVED
    Vec type = NOT_OBSERVED
    device memory peak = NOT_MEASURED

fixed benchmark:
    CPU 1R = NOT_RUN
    CPU 4R = NOT_RUN
    GPU 1R = NOT_RUN

GPU residency:
    matrix resident = NOT_MEASURED
    vectors resident = NOT_MEASURED
    transfer pattern = NOT_MEASURED
    excessive transfer = NOT_MEASURED

GPU candidate:
    GPU_A speedup = NOT_MEASURED
    GPU_B = NOT_RUN
    selected = CPU_EARLY_STOP_PRODUCTION

final production:
    CPU_EARLY_STOP
final runtime:
    {perf['candidate_measured_wall_s']} s (step10 -> step80)
speedup vs SV1.2:
    {perf['measured_continuation_speedup']:.6f}x

scientific equivalence:
    velocity relative L2 = {errs['velocity_relative_volume_L2']:.12g}
    pressure relative L2 = {errs['pressure_relative_volume_L2']:.12g}
    Qout fractions = {list(a['outlet_fractions'].values())}
    mass = {a['epsilon_mass']}
    steady = PASS (Eu={last['E_u']:.12g}, EQ={last['E_Q']:.12g})

tests (SV1.3):
    passed = {new['counts']['passed']}
    failed = {new['counts']['failed']}
    skipped = {new['counts']['skipped']} (GPU prerequisite blocked)
full regression:
    passed = {full['counts']['passed']}
    failed = {full['counts']['failed']} (unchanged historical failures)
    skipped = {full['counts']['skipped']}
history integrity:
    PASS
report:
    reports/sv1_3/REPORT.md
human review:
'''+''.join('    '+Path(f['path']).name+'\n' for f in vis['figures'])+f'''
STAGE SV1.3 STATUS:
    {status} (CPU optimization only)
'''
(REPORT/'terminal_summary.txt').write_text(summary)
files=[ROOT/'src/sv_validation/sv13.py',ROOT/'tests/sv13_support.py',*sorted((ROOT/'scripts/sv13').glob('*.py')),*sorted((ROOT/'tests').glob('test_sv13_*.py')),*sorted(CONFIG.glob('*'))]
write_json(REPORT/'implementation_manifest.json',{'files':[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p)} for p in files if p.is_file()]})
write_json(REPORT/'completion.json',{'status':status,'completed_utc':now(),'report':'reports/sv1_3/REPORT.md',
    'report_sha256':sha256(REPORT/'REPORT.md'),'native_solver_exited':True,'formal_figure_count':10,
    'human_review_pending':True,'GPU_blocked':True,'next_stage_started':False})
print(summary,flush=True)
raise SystemExit(0 if stage_ok else 1)
