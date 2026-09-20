"""Close the GPU development report; all science equivalence stays deferred."""
import json,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3n';sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,now
def read(n):return json.loads((R/(n+'.json')).read_text())
def write(n,d):(R/(n+'.json')).write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
s=read('petsc325_source');build=read('petsc_gpu13_build');a=read('compatibility_adapter');smoke=read('svmp_gpu_smoke');run=read('REAL_VASCULAR_GPU_acceptance');c=read('gpu_steady_candidate');h=read('gpu_steady_history')
assert c['status']=='PASS' and c['scientific_equivalence']=='DEFERRED'
assert read('preservation_audit')['status']==read('remote_preservation')['status']=='PASS'
cases=ET.parse(R/'pytest_full.xml').findall('.//testcase')
def count(rows):return dict(passed=sum(not any(c.find(t) is not None for t in ('failure','error','skipped')) for c in rows),failed=sum(c.find('failure') is not None for c in rows),errors=sum(c.find('error') is not None for c in rows),skipped=sum(c.find('skipped') is not None for c in rows))
new=[x for x in cases if 'test_sv13n_' in x.get('classname','')];old=[x for x in cases if x not in new]
tests=dict(stage=count(new),full=count(cases),historical=count(old),failures=[dict(test=x.get('classname')+'::'+x.get('name'),historical=x in old,message=x.find('failure').get('message')) for x in cases if x.find('failure') is not None]);write('test_summary',tests)
assert tests['stage']['failed']==tests['stage']['errors']==0
state=dict(completed=True,status='GPU_DEVELOPMENT_PASS',classification='GPU_STEADY_CANDIDATE',scientific_validation='NOT YET SCIENTIFICALLY VALIDATED',CPU_GPU_equivalence='DEFERRED',PETSc_version_equivalence='DEFERRED',performance_benchmark='DEFERRED',transfer_profiling='DEFERRED',GPU_native_PC_tuning='DEFERRED',supported_MPI_ranks=1,optional_MPI_CUDA_2R='FAIL: reverse/local-form stale owned values; not adopted',production='CPU_EARLY_STOP_PRODUCTION',production_changed=False,repair_iterations=a['repair_iterations'],timestamp=now());write('stage_result',state)
policy=json.loads((ROOT/'configs/sv1_3/policy.json').read_text());last=c['last_interval'];m=c['measurement'];rows=run['history']['linear_solves'];iterations=c['total_KSP_iterations']
device_samples=read('remote/vascular_gpu_device_snapshots')['samples']
device_memory=[float(line.split(',')[1].strip()) for sample in device_samples if sample['exit_code']==0 for line in sample['raw'].splitlines() if line.strip()]
memory=dict(status='OBSERVATIONAL ONLY',sampled_device_memory_max_MiB=max(device_memory),samples=len(device_memory),cadence_s=15,scope='Device-wide sampled maximum; not an exact process peak and not a residency/profile measurement',source='remote/vascular_gpu_device_snapshots.json')
write('observed_gpu_memory',memory)
oldcpu=read('OLD_PETSC_CPU_PROOF_20_acceptance')
portrows='\n'.join('| '+key+' | '+format(value,'.10g')+' |' for key,value in m['outlet_flows_m3_s'].items())
pressure_rows='\n'.join('| '+key+' | '+format(value,'.10g')+' |' for key,value in m['area_average_pressure_pa'].items())
report=f"""# Stage SV1.3N — 直接完成 SimVascular RTX 4090 GPU 血流求解开发

**GPU_DEVELOPMENT_PASS — GPU_STEADY_CANDIDATE（限定 1 MPI rank）。**

**NOT YET SCIENTIFICALLY VALIDATED。CPU/GPU 科学等价验证为 DEFERRED，正式 production 仍为 `CPU_EARLY_STOP_PRODUCTION`。**

## 为什么这次直接做 GPU？

用户在 [最新策略](USER_GPU_PRIORITY_UPDATE.txt) 中明确将目标改为 GPU 功能开发，把版本科学对照、CPU/GPU 等价、重复性能统计和重型 profiling 推迟。已通过的 source/API/PetscBool/CPU build/official smoke 直接继承。旧 CPU 作业收到原生 STOP_SIM 请求后在第 {oldcpu['steps_completed']} 步正常退出；既有 step10 VTU 和 checkpoint 留作开发记录，没有重启，也未用来生成科学等价结论。

本阶段只执行一条真实血管 GPU dynamic-to-steady run；没有 CPU short comparison、CPU4、重复 GPU timing 或独立长 profiling。WSL 是源码主目录。

## 使用了哪一版 PETSc？

版本 `{s['version']}`；tag `{s['tag']}`；完整 commit `{s['commit']}`；官方 source archive SHA256 `{s['source_archive_sha256']}`。

继承的 [来源核验](petsc325_source.json)、[API 清单](petsc_api_inventory.json)、[3.20—3.25 变更审计](release_change_audit.json)、[PetscBool 审计](petsc_bool_audit.json) 保持原样。采用官方 v3.25.5，不使用 main、浮动 release HEAD 或旧 ghost backport。PETSc 源码未修改。

## CUDA / MPI 是什么？

实际 CUDA **{build['CUDA_version']}**，RTX 4090 / sm_89，服务器默认 GCC/G++ **13.3**，Stage L 已验证 Open MPI **4.1.6**（完整 Fortran datatype 支持）。PETSc 本身未构建 Fortran binding。real double、32-bit indices、optimized shared、相同系统 BLAS/LAPACK。

CUDA13.2 的编译器兼容性由本次 configure 和实际 CUDA 自测试确认；没有使用 allow-unsupported-compiler。初次新脚本使用了错误的 configure 参数名，修正为源码实际声明的 CUDAFLAGS 后编译、安装及 CPU1/MPI2/CUDA checks 全部通过。随后修正了成功日志 CUDA 大小写解析，复用原成功运行，没有重做构建。[GPU build](petsc_gpu13_build.json)、[CUDA 路线](cuda_fallback_policy.json)。

![GPU 链路](gpu_stack_pipeline.png)

## 为新版 PETSc 修改了什么？

固定 svMultiPhysics 科学基线 `{read('svmp_gpu_build')['commit']}`。共 {a['repair_iterations']} 次 controlled repairs：

- `repair_01`：2 个 svMP PETSc 适配文件，+{a['added_lines']} / -{a['removed_lines']} 行。实现受 ownership/初始化状态保护的进程级 cleanup，检查 MPI 仍 active，销毁 PETSc 对象后恰好调用一次 PetscFinalize。未使用方程槽以零初始化的 null handles 安全清理，新增 PETSc cleanup 调用全部检查错误。
- `repair_02`：仅修正本阶段构建脚本的 CUDAFLAGS 参数和大小写观察逻辑，不修改 PETSc 源码或 CFD 数学。

实际单次官方 GPU 测试在 application MPI_Finalize 入口观察到 initialized=0、finalized=1、PetscFinalize call count=1，并正常退出。另一个最小 CUDA 对象创建、运算、销毁和 finalization 测试也以 exit0 结束。[生命周期跟踪](solver_finalize_lifecycle.json)、[最小测试](remote/minimal_lifecycle.json)、[修复明细](compatibility_adapter.json)、[补丁及行号依据](../../patches/sv1_3n/PATCH_PROVENANCE.md)。

这符合 [PETSc PetscFinalize 的官方生命周期约定](https://petsc.org/release/manualpages/Sys/PetscFinalize/)。FEM、Navier–Stokes、矩阵组装、边界、mesh、rho、mu、Q、dt、时间积分及容差均未修改。原 CPU binary 和历史工程保留。

## CUDA ghost vector 解决了吗？

**本阶段采用的 1-rank seqcuda 路径 PASS。** 实际核对了 local form、owned values、forward/reverse 和 restore；该路径与真实单 rank solver 一样没有跨 rank ghost。

**额外 2-rank mpicuda 检查 FAIL，未采用。** Forward 值正确，但 reverse 后 CPU local form 读取到旧 owned 值（最大误差分别为 111 和 101），程序按数值错误拒绝并 exit77。CPU seq/MPI 对照通过。原始失败完整保存，不能把新版 ghost 功能概括为全部通过，也不宣称支持多个 MPI rank 共享这张 GPU。本阶段用户指定 1 MPI rank，故继续验证通过的小测试、官方流体和真实血管路径。[逐元素原始结果](ghost_probe_new.json)、[采用范围](gpu_ghost_target.json)。

![Ghost 状态](cuda_ghost_status.png)

## SimVascular 是否真的用了 GPU？

真实记录 Mat=`{run['mat_type']}`，Vec=`{run['vec_type']}`，KSP=`{run['KSP']}`，PC=`{run['PC']}`。每次线性求解的 KSPView/PCView 都核对 right preconditioning、restart100、rtol1e-10、atol1e-24、max2000、ASM overlap2、ILU(2)、natural ordering 和 zero-pivot 默认值。

1 MPI rank、1 RTX4090、OMP=1。[链接审计](svmp_gpu_link.json)、[真实运行验收](REAL_VASCULAR_GPU_acceptance.json)。CUDA Mat/Vec 是实际类型证据；本阶段不把它扩展成全部预条件器操作常驻 GPU 或性能收益结论。

## 官方 GPU fluid smoke 是否通过？

**PASS**：{smoke['steps_completed']} 个完整非线性步、{smoke['linear_solves']} 次线性求解、0 线性失败、0 非线性失败，有限速度/压力、VTU 输出、fresh-process reload、正常退出全部通过。历史 MPI_ERR_TYPE、Vector is not ghosted 和 PETSc-after-MPI 错误均未出现在本次采用的 1-rank 路径。[官方测试](svmp_gpu_smoke.json)。

## 真实 vascular GPU 是否跑起来？

**PASS**：从 t=0 开始，完成 {c['steps']} 步、{c['linear_solves']} 次 KSP solve、总 {iterations} 次 KSP iteration。0 线性失败、0 非线性失败；正式检查点的有限场、边界流量和壁面检查通过。[完整求解记录](REAL_VASCULAR_GPU_acceptance.json)。

![实际求解进展](gpu_solver_progress.png)

## GPU 流场什么时候达到稳态？

首次完整通过 production 联合判据：**step {c['first_full_steady_step']}**；要求 **{c['required_consecutive_intervals']} 个连续区间**，每 {policy['save_interval_steps']} 步检查。冻结阈值 E_u≤{policy['velocity_change_limit']:.0e}、E_Q≤{policy['flow_change_limit']:.0e}，mass 和 inlet gate≤{policy['mass_limit']:.0e}，并检查有限场、无滑移及无 solver failure。

最后正式区间 E_u={last['E_u']:.10g}、E_Q={last['E_Q']:.10g}；最终 mass={m['epsilon_mass']:.10g}。计算直接复用生产 SteadyStopMonitor/stop_gate，KSP convergence 不代替流场稳态。[监测历史](gpu_steady_history.json)、[冻结 production policy](../../configs/sv1_3/policy.json)。

![正式稳态判断](gpu_steady_convergence.png)

## 为什么在那里停止？

第一次完整满足连续门槛后立即原子写入 STOP_SIM=0；只让已经运行的当前步完成并保存，**step {c['stop_step']}** 正常退出，物理时间 **{c['physical_time_s']:.12g} s**。没有额外安全余量步、重复验证步或统计运行。

本阶段唯一 I/O 调整是每步保存 VTU，使原生 STOP_SIM 能在当前步产生相同时间点的 final VTU 和 restart；正式 steady check cadence 始终为 10。原 production 的 maximum_total_steps 只作未收敛保护上限，不是运行目标。未向 solver 发送 kill 信号。[安全停止请求](gpu_steady_stop_request.json)、[配置策略](../../configs/sv1_3n/convergence_strategy.json)。

## 最终 GPU candidate 是否完整？

**PASS**：final VTU、完整 1-rank native checkpoint、fresh-process reload、finite velocity/pressure、Qin/三出口流量、wall no-slip 和 mass 均核验。原生 checkpoint 包含 Y_n、A_n、step/time 和方程初始范数，长度及所有数值均检查。

VTU：`{c['VTU']}`；SHA256 `{c['VTU_sha256']}`。

Checkpoint：`{c['checkpoint']['path']}`；SHA256 `{c['checkpoint']['sha256']}`。

Qin={m['Q_in_m3_s']:.10g} m³/s；Qout total={m['Q_out_total_m3_s']:.10g} m³/s；max velocity={m['velocity_max_m_s']:.10g} m/s；wall max={m['wall_velocity_max_m_s']:.10g} m/s。

| 出口 | 实际流量 / m³/s |
|---|---:|
{portrows}

| 端口 | 原生面积平均压力 / Pa |
|---|---:|
{pressure_rows}

这些是 GPU 本身的完整性记录，没有与 CPU pressure/flow fraction 做等价比较。[候选清单和全部哈希](gpu_steady_candidate.json)、[WSL 源码镜像](WSL_source_mirror.json)、[原生构建归档](native_artifact_mirror.json)。

![最终 GPU 原生场](gpu_final_field.png)

必要真实 GPU 运行的 wall time={c['wall_time_s']:.6f} s，steps={c['steps']}，KSP solves={c['linear_solves']}，total KSP iterations={iterations}，mean iterations/solve={iterations/c['linear_solves']:.6g}。**OBSERVATIONAL ONLY**：没有 median、speedup 或正式性能分类。每 15 秒采样的设备总显存最大观测值为 {memory['sampled_device_memory_max_MiB']:g} MiB（{memory['samples']} 个样本）；这是整卡采样值，不是精确进程峰值，也不证明 GPU residency。[原始记录](remote/vascular_gpu_device_snapshots.json)。

## 哪些验证被明确推迟？

- OLD vs NEW PETSc science comparison — **DEFERRED**。
- CPU vs GPU field equivalence — **DEFERRED**。
- Formal speed benchmark — **DEFERRED**。
- GPU transfer profiling — **DEFERRED**。
- GPU-native PC tuning — **DEFERRED**。

原因是用户明确要求当前先完成功能性 GPU steady solver；没有运行对应额外 CFD，也没有把 NOT RUN 写为 PASS。

## 当前结论

**GPU_DEVELOPMENT_PASS；GPU_STEADY_CANDIDATE；NOT YET SCIENTIFICALLY VALIDATED。** 采用范围仅为本次验证的单 rank CUDA 路径；2-rank CUDA ghost reverse 问题明确保留为限制。

正式 `CPU_EARLY_STOP_PRODUCTION` 保持不变。没有自动启动科学对照、性能调优、mesh convergence、particle coupling 或额外 timestep。

本阶段测试 {tests['stage']}；完整 pytest {tests['full']}；历史测试结果 {tests['historical']}。新增测试读取实际 artifact，不重新启动 CFD。额外 MPI CUDA 失败通过拒绝门槛与已知不支持范围检查保留，未伪称功能成功。[测试明细](test_summary.json)、[WSL 历史保留审计](preservation_audit.json)、[远端保留审计](remote_preservation.json)、[原生 mesh/geometry 输入哈希](remote/final_native_input_integrity.json)。
"""
(R/'REPORT.md').write_text(report)
write('report_manifest',dict(status='PASS',report_sha256=sha256(R/'REPORT.md'),timestamp=now(),stage_status=state['status']))
print('GPU development report written from closed native records; comparisons remain deferred.')
