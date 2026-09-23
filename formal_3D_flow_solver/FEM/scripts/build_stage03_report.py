#!/usr/bin/env python3
"""Report the measured Stage 3 failure, preserving unavailable results as unavailable."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
from fem3d.vascular import load_config,reynolds
R=ROOT/'reports/stage03';B=ROOT/'outputs/stage03/reference';read=lambda p:json.loads(p.read_text())
config,config_hash=load_config(ROOT);f=read(B/'metadata/failure.json');d=read(B/'qc/singularity_diagnosis.json')
a=read(ROOT/'outputs/stage03/preflight/assembly_verified.json');m=read(ROOT/'outputs/stage03/mesh/qc.json')
resources=read(B/'metadata/resources.json');accounting=read(B/'metadata/resource_accounting.json')
initial=read(ROOT/'outputs/stage03/resources/initial.json');plan=read(R/'stage018_cleanup_plan.json')
ref=read(R/'reference_integrity.json');history=read(R/'history_preservation.json');visual=read(R/'visualization_manifest.json')
assert f['status']=='FAIL' and ref['status']==history['status']=='PASS'
suites=ET.parse(R/'pytest_results.xml').getroot().findall('testsuite')
stats={k:sum(int(s.attrib.get(k,0)) for s in suites) for k in ('tests','failures','errors','skipped')}
stats['passed']=stats['tests']-stats['failures']-stats['errors']-stats['skipped']
stage03=[c for s in suites for c in s.findall('testcase') if 'test_stage03_' in c.attrib.get('classname','')]
stats['stage03']={'total':len(stage03),'skipped':sum(c.find('skipped') is not None for c in stage03)}
stats['stage03']['passed']=stats['stage03']['total']-stats['stage03']['skipped']
assert stats['failures']==stats['errors']==0
write_json(R/'test_summary.json',stats)
names=plan['forbidden_package_names'];tools_label=names[0].upper()+'/'+names[-1][:3].title()+names[-1][3:].title()
physics=config['physics'];re=reynolds(config);Q=physics['inlet_volume_flow_m3_s'];gib=1024**3
contract={'schema_version':1,'stage':3,'stage_status':'FAIL','availability':'NO_VALID_SOLUTION',
          'reason':f['reason'],'downstream_consumption_allowed':False,'experimental':False,
          'warning':'NOT_EXPERIMENTAL_PUMP_FLOW','mesh_source':config['mesh']['source'],
          'mesh_sha256':config['mesh']['volume_mesh_sha256'],'config_sha256':config_hash,
          'velocity':{'family':'Lagrange','order':2,'value_shape':[3],'units':'m/s','component_order':['x','y','z'],
                      'primary_representation':'Original continuous P2 nodal coefficients; no projection','available':False},
          'pressure':{'family':'Lagrange','order':1,'units':'Pa','gauge_convention':'Natural zero traction at all three outlets; no pointwise p=0 condition, no pressure pin or registered pressure nullspace',
                      'observed_discrete_failure':'Two independent localized pressure modes decouple under wall Dirichlet elimination','available':False},
          'gradient':{'sampling_representation':'Cell-center DG0 samples of the P2 gradient, intended for visualization and QC only',
                      'definition':'gradient[i,j] = du_i/dx_j','units':'s^-1','component_order':['xx','xy','xz','yx','yy','yz','zx','zy','zz'],
                      'particle_interpolation_field':False,'available':False},
          'strain_rate':{'definition':'0.5*(gradient+transpose(gradient))','units':'s^-1','available':False},
          'vorticity':{'definition':'curl(u)','units':'s^-1','component_order':['x','y','z'],'available':False},
          'lambda':{'meaning':'Single global inlet normal traction multiplier: sigma*n = -lambda*n','units':'Pa','available':False},
          'boundary_flux':{'units':'m^3/s','sign_convention':'Outward-normal flux; prescribed inlet signed flux = -Q; positive outlet net flow means outflow',
                           'tags':{'WALL':1,'OUTLET_03':2,'OUTLET_01':3,'INLET':4,'OUTLET_02':5,'FLUID':100},
                           'source_required':'Actual FEM facet integrals, never configuration flow substituted as a measurement','available':False},
          'checkpoint':{'expected_primary':'outputs/stage03/reference/checkpoints/primary.npz','exists':False,'roundtrip':'NOT_EXECUTED_NO_SOLUTION'},
          'wss_validated':False,'particle_sampler_implemented':False}
write_json(R/'flow_field_contract.json',contract)
cleanup={site:read(ROOT/f'outputs/stage03/cleanup/{site}_final_audit.json') for site in ('wsl','remote')}
assert all(v['status']=='PASS' and v['fem_environment_unchanged'] for v in cleanup.values())
write_json(R/'final_cleanup_audit.json',cleanup)
status={'stage':3,'status':'FAIL','reason':'SINGULAR_PRESSURE_SUPPORT','timestamp':timestamp(),
        'cleanup':'PASS','selected_mesh_sha':'PASS','mesh_boundary_qc':'PASS','spaces':'PASS',
        'assembly_only_preflight':'PASS','production_factorization':'FAIL','solution_available':False,
        'flux_constraint':'NOT_EVALUATED_NO_SOLUTION','mass_conservation':'NOT_EVALUATED_NO_SOLUTION',
        'field_finiteness':'NOT_EVALUATED_NO_SOLUTION','wall_solution_values':'NOT_EVALUATED_NO_SOLUTION',
        'local_flow_finiteness':'NOT_EVALUATED_NO_SOLUTION','solution_roundtrip':'NOT_EXECUTED_NO_SOLUTION',
        'reference_integrity':ref['status'],'history_preservation':history['status'],'pytest':stats,
        'formal_factorization_attempts':1,'resource_blocked':False,'experimental_run':'NOT_EXECUTED',
        'flow_plots_available':False,'human_review':'Failure evidence only; cannot approve flow results that do not exist',
        'next_stage_started':False}
write_json(R/'final_status.json',status)
zero_rows='\n'.join(f"| {r['matrix_global_row']} | {r['physical_source_match']['source_vertex_index']} | {', '.join(map(str,r['physical_source_match']['support_source_cells']))} | {r['actual_nonzero_entries']} | {r['unit_pressure_mode_matrix_product_norm']:g} |" for r in d['zero_rows'])
report=f'''# Stage 3 — 真实三维血管 FEM 流场首次求解

**STAGE 3 STATUS: FAIL — SINGULAR_PRESSURE_SUPPORT。**
正式 4-rank MUMPS 数值分解失败，未得到有效流场。独立的仅组装诊断确认：两个局部压力基函数对应完全为零的行和列，两个独立单位压力向量均满足 `A·e=0`。本阶段停止，没有改变网格或经过验证的 Stage 2 核心来绕过失败。

## 开始前清理了什么？

按用户要求，先建立 [cleanup plan](stage018_cleanup_plan.json)，通过导入图、Git 状态、文本搜索和目录清单确定范围，再删除 {tools_label} 活动实现、专用测试、配置、原始结果及独立 mesh-tools 环境。计划包含 {len(plan['remove_files'])} 个文件条目和 {len(plan['remove_directories'])} 个目录树条目。没有新建代码备份。

保留通用质量、边界和成本诊断模块。Stage 1.8 的最终 REPORT、assessment_summary、final_status、11 张最终 PNG 保持原字节；旧版归档仅保留最终科学证据，不再保留代码与原始实验。新增 [CLEANUP_NOTICE](../stage01_8/CLEANUP_NOTICE.md)，历史结论仍为 KEEP_STAGE017。

WSL 与远端最终活动代码/包检查均 PASS。远端原 FEM Python 与全部 conda 元数据文件的指纹未变；独立环境确已删除。清理后首次完整 remaining tests 为 269 passed、5 skipped、0 failed；本阶段最终回归见末节。[最终清理审计](final_cleanup_audit.json)

## 这一步在求什么？

Stage 1.7 提供真实四面体血管网格，Stage 2 提供通过解析基准的三维稳态不可压 Newtonian Stokes 方程、P2/P1 + 一个全局 Real 及完整应力形式 `σ=2με(u)−pI`。Stage 3 首次组合它们，检验真实网格上的可解性和数值守恒。

唯一正式来源为 `outputs/stage01_7/selected`，volume_mesh SHA256：`{config['mesh']['volume_mesh_sha256']}`。全部 {a['mesh']['cells']:,} 个 tetra 重新计算质量：minSICN={m['quality']['min_sicn']['minimum']:.12g}、P1={m['quality']['min_sicn']['P1']:.12g}、P5={m['quality']['min_sicn']['P5']:.12g}，低于 0.1 的仍为 3 个。连通域 1；零、负、非有限体积均为 0；wall/rim 位移为 0。几何质量通过不代表与指定混合有限元空间组合后一定满秩。

## 使用的物性来自哪里？

来自只读当前正式配置 `{config['source']['configuration_path']}`，并与 Stage 0 几何源 lineage 对齐；不是手写 Stage 2 管道参数。[物性来源与哈希](reference_condition_provenance.json)

| 量 | SI 数值 |
|---|---:|
| ρ | {physics['density_kg_m3']:.12g} kg/m³ |
| ν | {physics['kinematic_viscosity_m2_s']:.12g} m²/s |
| μ=ρν | {physics['dynamic_viscosity_pa_s']:.12g} Pa·s |
| Q_reference | {Q:.16e} m³/s |

**REFERENCE NUMERICAL CONDITION — NOT EXPERIMENTAL PUMP FLOW。** 正式配置 `{config_hash}` 冻结于第一次组装之前。源配置中的出口压力没有导入：本次三个出口按用户指定采用相同零牵引参考。因此也没有声称本阶段已匹配 LBM 全部边界条件。

远端冻结环境无 PyYAML；初次组装在读取配置前退出。随后使用与原 YAML SHA 绑定的 JSON 执行快照，源 YAML 未改、环境未安装新包。实验模板 Q 仍为 null，**EXPERIMENTAL RUN: NOT EXECUTED**。

## 边界条件是什么？

入口只约束泵入的总流量 `∫inlet u·n dS = −Q`，其自然牵引为 `σn=−λn`；三个出口处于相同 atmospheric zero-traction reference，即 `σn=0`；血管壁 no-slip `u=0`。没有入口速度剖面、预设出口比例、pressure pin 或 pressure Dirichlet；wall/cap rim 仍属于壁面约束。

标签完全保留：WALL=1、OUTLET_03=2、OUTLET_01=3、INLET=4、OUTLET_02=5、FLUID=100。恰好一个入口、三个出口，未标记外表面为 0。

![真实几何与边界](real_geometry_and_bc.png)

## 计算规模有多大？

| 实测量 | 数值 |
|---|---:|
| vertices / edges / tetra | {m['proxy']['N_vertex']:,} / {m['proxy']['N_edge']:,} / {a['mesh']['cells']:,} |
| velocity P2 DOF | {a['velocity_dofs']:,} |
| pressure P1 DOF | {a['pressure_dofs']:,} |
| global Real DOF | {a['real_dofs']}，仅 rank 0 拥有 |
| 总未知量、矩阵行列 | {a['total_dofs']:,} × {a['total_dofs']:,} |
| 存储 nnz | {a['nnz']:,} |
| 壁面速度 Dirichlet DOF | {a['boundary_conditions']['velocity_dirichlet_dofs']:,} |
| 矩阵/向量组装时间 | {a['assembly_wall_time_s']:.6f} s |
| 空间、表单、JIT 设置 | {a['spaces_forms_jit_s']:.6f} s |
| 组装最大单 rank peak RSS | {a['max_rank_peak_rss_kib']/1024**2:.6f} GiB |

[核实后的仅组装记录](../../outputs/stage03/preflight/assembly_verified.json) 明确没有 factorize 或 KSP solve。首份 `assembly.json` 的 nnz 错将 PETSc 默认全局和再次累加；已通过显式 `InfoType.LOCAL` + MPI SUM 重新组装核实，保留原错误记录用于追踪，报告只使用核实值。nnz 是存储结构数量，并不保证每个值非零或矩阵可逆。

## 求解器是否正常结束？

**没有。** 唯一正式尝试采用冻结 Stage 2 `preonly + LU + MUMPS`，4 MPI ranks、OMP_NUM_THREADS=1、GPU=false。MUMPS 报 `INFOG(1)={d['mumps_infog_1']}, INFO(2)={d['mumps_info_2']}`；PETSc 将 −10 分类为数值奇异/零 pivot。[PETSc 官方实现](https://petsc.org/release/src/mat/impls/aij/mpi/mumps/impl/imumps.c.html)

KSP 在 setup/数值分解阶段抛错，因此没有可报告的成功 converged reason、solution residual、λ 或有限字段。全部原始 stdout/stderr、正式代码 SHA 与配置快照仍在 [正式运行记录]({'../../'+f['formal_run_record']})，摘要见 [failure_summary](failure_summary.json)。没有第二次 factorization、MPI rank sweep、iterative fallback 或新增压力固定点。

随后以相同网格、表单与壁面条件**仅重新组装**做诊断：以下两行所有数值严格为 0；对相应单位向量计算的矩阵乘积也严格为 0，独立于 MUMPS 的报错分类。[实际矩阵证据](singularity_diagnosis.json)

| MPI 矩阵全局行（0-based） | 物理匹配源顶点 | P1 支撑源单元 | 实际非零项 | ‖A·e‖ |
|---:|---:|---|---:|---:|
{zero_rows}

这里的源数组编号只用于核实后的溯源；在 DOLFINx 中通过物理坐标匹配，不用旧编号定位。已经证明至少两个独立精确零模式，没有宣称完成全矩阵的完整谱分解。

## 流量是否真的满足？

未能验证。只有配置目标 Q={Q:.16e} m³/s 已知；实际入口、三个出口、总出口、入口约束误差和质量闭合误差均为 **NOT_EVALUATED_NO_SOLUTION**。没有把 Q_reference 当作实际 FEM 积分，也没有把缺失值填成 0。

预先冻结的 `epsilon_Q ≤ 1e−10` 和 `epsilon_mass ≤ 1e−10` 保持不变；它们没有通过，原因是不存在可验收的解。

![流量结果不可用](flux_balance.png)

## 三个出口各分到多少流量？

三个带符号 Q_i 和流量比例均未知；没有预设分流。正常期望正净流出仍不是自动 hard gate；若将来实际解出现负净流出，需标 MANUAL_PHYSICS_REVIEW。本次不存在可解释的出口符号。

![分流结果不可用](outlet_flow_split.png)

## 速度场是什么样？

没有有效 P2 速度系数，不能生成真实三维速度图或 2–4 个内部速度切面。实际壁面 DOF 速度、divergence L2 和相对诊断也未获得；创建了 no-slip BC 不等于验证了求解后 no-slip 数值。

![速度场不可用](velocity_global.png)

![速度切面不可用](velocity_slices.png)

![梯度切面不可用](velocity_gradient_slice.png)

## 压力场是什么样？

没有有效 P1 压力系数或 λ。入口、三个出口和内部截面面积平均压力均未计算。**零牵引出口不等于 pointwise p=0**；本次两个局部压力零模式属于壁面约束后离散支撑缺失，不是可随意用 pressure pin 消除的全局 gauge 自由度。

![压力场不可用](pressure_global.png)

![平均压力不可用](pressure_sections.png)

## Stokes 假设是否合理？

用真实入口投影几何和参考 Q 计算，而非使用不存在的求解速度：A={re['inlet_area_m2']:.12e} m²，P={re['inlet_perimeter_m']:.12e} m，Dh=4A/P={re['hydraulic_diameter_m']:.12e} m，Umean=Q/A={re['mean_inlet_velocity_m_s']:.12e} m/s。

**Re=ρ Umean Dh/μ={re['Re']:.12g}**，远低于 1，支持该参考条件下采用 Stokes 近似。它不证明网格收敛、实验物理充分成立或当前离散系统可解；也不改变本阶段 FAIL。

## 那 3 个低质量 tetra 有影响吗？

三个物理质心全部在当前 DOLFINx 网格中重新定位，匹配位移均为 0 m。两个相邻 residual tetra（minSICN≈0.062380、0.088715）的 4 个顶点和 6 条边中点全部被 wall no-slip 固定，因此整个 P2 速度多项式在这两个 tetra 上均被约束为零。

两个 P1 压力基函数的全部支撑仅落在这些 tetra 上。壁面消元后，它们与自由速度的散度耦合消失，产生上述两个严格零行、零列。这个实际秩缺陷已经超出“低质量但解是否有局部 spike”的筛查问题；是本次不可解的直接证据。第三个 residual tetra 完成位置核实，未获得流场样本。

速度、压力、梯度 Frobenius norm、应变率 norm、涡量以及一阶面邻居 min/median/max/ratio 均 **未评估**。没有使用固定 2×/5× 判据，没有平滑，也不能给 local finite hard gate 盖章。对低质量单元的完整流场影响仍需有效解和后续收敛研究；当前已经证实的是两个局部压力零模式。

![三个位置](residual_cell_locations.png)

![残余单元失败诊断](residual_cell_flow_check.png)

## 真实计算花了多少资源？

远端 `{f['hostname']}`：AMD Ryzen 7 7800X3D，8 物理核/16 线程，CPU cgroup 额度 7.68 核；物理 RAM {initial['memory_bytes']['MemTotal']/gib:.3f} GiB，容器 memory.max {int(initial['cgroup']['memory.max'])/gib:.3f} GiB。正式尝试前有效可用 RAM {resources['fresh_memory_before']['effective_available_bytes']/gib:.3f} GiB；磁盘约 {initial['disk']['free']/gib:.1f} GiB 空余。

冻结环境：DOLFINx {initial['dolfinx_version']}、PETSc {'.'.join(map(str,initial['petsc_version']))}、MUMPS {initial['packages']['mumps-mpi']}、MPICH {initial['packages']['mpich']}；无 GPU。仅组装 {a['assembly_wall_time_s']:.6f} s，正式失败尝试整体 {resources['elapsed_time_s']:.6f} s。后者包含 MPI 启动、组装和失败分解，**不能称为成功 solve time**；没有成功 factorization/solve 的独立时长。

Linux `resource.getrusage(RUSAGE_CHILDREN)` 记录最大单子进程峰值 **{accounting['valid_peak_rss_kib']/1024**2:.6f} GiB**。最初的进程组 RSS 采样未包含使用不同 process group 的 MPICH worker，约 5.5 MiB 的错误总量已明确排除，不能当成整个求解的内存。容器整体采样峰值 {accounting['peak_sampled_cgroup_memory_bytes']/gib:.3f} GiB 包含其他进程与缓存，也不是 FEM RSS。完整 MPI 同时 RSS 总峰值不可恢复，未为弥补计量重复分解。[资源计量说明](resource_accounting.json)

运行期间最小有效可用内存 {accounting['minimum_sampled_available_memory_bytes']/gib:.3f} GiB，cgroup OOM/oom_kill 均为 0，没有资源守卫终止；本次归类数学/离散系统 FAIL，**不是 DIRECT_SOLVER_RESOURCE_BLOCKED**。

![时间与内存](solver_resource_usage.png)

## 结果能否重新读取？

**NOT_EXECUTED_NO_SOLUTION。** 分解未成功，未创建 `checkpoints/primary.npz`，没有速度、压力或 λ checkpoint。不能用 mesh reload 冒充 solution roundtrip。保存并复核的是实际网格、组装、失败日志、资源记录、两个矩阵零模式与三个 residual 的物理定位证据。

[flow_field_contract.json](flow_field_contract.json) 记录 P2/P1、单位、梯度采样表示、λ 含义和流量符号，同时明确 `NO_VALID_SOLUTION`、`downstream_consumption_allowed=false`。预备的结果重读/派生场代码未能在真实解上执行，不视为已经验证；cell-center DG0 样本也不被宣称为 particle interpolation field。

## 当前还没有做什么？

没有 FEM–LBM 正式比较、mesh convergence、实验 Q_pump 正式 run、WSS validation、RBC/MB、粒子采样器或线性 scaling utility。没有更换网格、修改 Stage 2 核心、增加 stabilization 或 pressure pin 来强行获得解。

11 个要求的图路径均已建立，其中 **3 张实际几何/资源图、1 张实际失败诊断表、7 张明确标注“未生成流场图”的状态页**；实际流场图数量为 0。这些状态页不是原请求的成功流场交付。全部在 WSL 生成，详见 [可视化清单](visualization_manifest.json)。

## 是否可以进入下一阶段？

**不能。STAGE 3 STATUS: FAIL。REASON: SINGULAR_PRESSURE_SUPPORT。**
清理、网格 SHA/边界、P2/P1/单 Real 和仅组装检查通过；正式分解失败。流场有限性、实际 no-slip、入口流量、质量守恒、局部流场有限性及 solution roundtrip 均未评估，不能因为 regression tests 无失败就标 CONDITIONAL PASS。

完整 pytest：**{stats['passed']} passed、{stats['failures']} failed、{stats['errors']} errors、{stats['skipped']} skipped**。Stage 3 的 20 个永久测试模块合计 {stats['stage03']['passed']} passed、{stats['stage03']['skipped']} skipped；跳过均因正式解不存在，成功路径测试永久保留。另有 5 个历史跳过。矩阵零模式、拓扑支撑、无伪造结果和残余物理定位的失败路径测试已通过。[测试 XML](pytest_results.xml)

两个只读参考项目共 {ref['entries']:,} 个条目的全量完整性检查 PASS；{history['frozen_file_count']:,} 个历史文件和全部保留的旧 fem3d 源文件未改。用户已有 Stage 1 REPORT 修改按原字节保留，不纳入本次提交。[参考审计](reference_integrity.json) · [历史审计](history_preservation.json)

需要另行决定如何处理这个真实网格与指定混合空间的兼容性问题；本阶段按既定边界停止。没有启动下一阶段。
'''
(R/'REPORT.md').write_text(report)
summary=f'''Stage 3 completed (failed numerical preflight; no valid solution).

cleanup:
    {tools_label} active code removed = true
    mesh tools env removed = true
mesh:
    source = Stage 1.7 selected
    tetra = {a['mesh']['cells']}
    velocity DOF = {a['velocity_dofs']}
    pressure DOF = {a['pressure_dofs']}
    Real DOF = 1
    total DOF = {a['total_dofs']}
reference condition:
    Q = {Q:.16e} m^3/s
    rho = {physics['density_kg_m3']}
    nu = {physics['kinematic_viscosity_m2_s']}
    mu = {physics['dynamic_viscosity_pa_s']}
    experimental = false
solver:
    MPI ranks = 4
    converged = false
    factorization = FAILED, MUMPS INFOG(1)=-10
    solve time = UNAVAILABLE (no completed solve)
    failed attempt elapsed = {resources['elapsed_time_s']:.6f} s
    peak RSS = {accounting['valid_peak_rss_kib']} KiB (maximum individual child; MPI simultaneous sum unavailable)
flux:
    target = {Q:.16e}
    inlet = NOT_EVALUATED_NO_SOLUTION
    outlet_01 = NOT_EVALUATED_NO_SOLUTION
    outlet_02 = NOT_EVALUATED_NO_SOLUTION
    outlet_03 = NOT_EVALUATED_NO_SOLUTION
    total outlet = NOT_EVALUATED_NO_SOLUTION
    inlet error = NOT_EVALUATED_NO_SOLUTION
    mass closure = NOT_EVALUATED_NO_SOLUTION
lambda:
    UNAVAILABLE
Re:
    {re['Re']:.12g} (reference Q and physical projected inlet geometry)
residual cells:
    relocated = 3
    flow monitored = 0 (no solution)
    finite = NOT_EVALUATED
    exact independent pressure null modes = 2
    manual anomaly review = no flow to review; singularity evidence available
roundtrip:
    NOT_EXECUTED_NO_SOLUTION
tests:
    passed = {stats['passed']}
    failed = {stats['failures']}
    skipped = {stats['skipped']} (5 historical; 9 without a Stage 3 solution)
report:
    reports/stage03/REPORT.md
human review:
'''+''.join('    '+name+' ['+record['kind']+']\n' for name,record in visual['figures'].items())+'''
STAGE 3 STATUS:
    FAIL
REASON:
    SINGULAR_PRESSURE_SUPPORT
STOP. No next stage started.
'''
(R/'terminal_summary.txt').write_text(summary)
write_json(R/'delivery_manifest.json',{'timestamp':timestamp(),'status':'FAIL','files':{str(p.relative_to(ROOT)):sha256(p) for p in sorted(R.iterdir()) if p.is_file() and p.name not in ('reference_final.json','delivery_manifest.json')}})
print(summary)
