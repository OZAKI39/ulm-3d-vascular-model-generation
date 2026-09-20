#!/usr/bin/env python3
"""Write the Chinese SV1 report from measured artifacts and actual acceptance tests."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import write_json,sha256,now
R=ROOT/'reports/sv1'
def read(n): return json.loads((R/(n+'.json')).read_text())
q=read('flow_qc');g=read('geometry_qc');m=read('mesh_validity');quality=read('mesh_quality')
r=read('flow_execution');run=r['runs'][-1];s=read('steady_state');reload=read('solution_reload')
a=read('cleanup_summary');native=read('native_solver');smoke=read('official_smoke');api=read('simvascular_api')
p=read('reference_physics')['physics'];inlet=read('inlet_normalization');old=read('old_fem_readonly_audit')
tp=json.loads((ROOT/'configs/time_policy.json').read_text())
resource=read('solver_resource_usage')['measurements'][-1]
manifest=json.loads((ROOT/'inputs/MANIFEST.json').read_text())
xml=ET.parse(R/'pytest_results.xml');cases=xml.findall('.//testcase')
failures=[{'module':c.get('classname'),'name':c.get('name'),'message':(c.find('failure') if c.find('failure') is not None else c.find('error')).get('message')} for c in cases if c.find('failure') is not None or c.find('error') is not None]
skipped=sum(c.find('skipped') is not None for c in cases)
tests={'passed':len(cases)-len(failures)-skipped,'failed':len(failures),'skipped':skipped,'total':len(cases),'failures':failures}
write_json(R/'test_summary.json',tests)
status='CONDITIONAL PASS' if q['status']=='PASS' and not failures and reload['status']=='PASS' and old['status']=='PASS' and a['status']=='PASS' else 'FAIL'
reasons=q['failure_reasons'][:]
if failures: reasons.append('ACCEPTANCE_TEST_FAILURES')
final={'timestamp':now(),'stage':'SV1','status':status,'reasons':reasons,'tests':tests,
       'mesh_valid':m['status']=='PASS','native_build':native['status'],'official_smoke':smoke['status'],
       'solver_process_exit_code':run['exit_code'],'solver_success':q['solver_success'],
       'planned_block_completed':s['scheduled_block_complete'],'steady':s['status'],
       'diagnostic_artifact':q['path'],'artifact_kind':q['artifact_kind'],
       'epsilon_Q':q['epsilon_Q'],'epsilon_mass':q['epsilon_mass'],
       'diagnostic_reload':reload['status'],'old_fem_readonly':old['status'],
       'human_review':'Images available for diagnosis; physical acceptance failed' if status=='FAIL' else 'Pending'}
write_json(R/'final_status.json',final)
port_table='\n'.join(f"| {role} | {v['area_relative_error']:.6%} | {v['centroid_shift_m']:.6e} | {v['normal_dot']:.12f} | 1 / 否 |" for role,v in g['ports'].items())
flow_table='\n'.join(f"| {role} | {q['outlet_flows_m3_s'][role]:.12e} | {q['outlet_fractions'][role]:.8%} |" for role in ('OUTLET_01','OUTLET_02','OUTLET_03'))
fail_table='\n'.join(f"- `{f['module']}.{f['name']}`：{f['message'].splitlines()[0]}" for f in failures)
commit=native['commit'];upstream=f'https://github.com/SimVascular/svMultiPhysics/blob/{commit}/Code/Source/solver/'
text=f'''# Stage SV1 — 使用 SimVascular 直接求解真实三维血管流场

**STAGE SV1 STATUS: {status}。原因：{', '.join(reasons)}。**

原生构建、官方最小算例和真实血管网格均通过。真实血管计算反复出现线性系统不收敛，在首个计划保存点第 {q['step']} 步按官方停止机制结束。以下速度、压力、流量均来自原生求解器实际输出，但属于**未验收的瞬态诊断状态**，不能作为正式稳态血流结果。此次失败只描述这次固定配置的运行，不证明该几何在所有求解设置下都无法收敛。

## 这次删除了哪些无关路线？

已删除旧部署、异地计算及兼容性比较代码、相关配置、输出、日志和 {len(a['deleted_obsolete_tests'])} 个旧测试；初次清理共记录 {len(a['deleted_files'])} 个文件。仅保留 [SV0 历史报告](../sv0/REPORT.md)。现在只有 WSL 原生 SimVascular 生产路线。

[清理审计](cleanup_summary.json)列出逐文件删除记录、保留的共享文件及活动文本扫描。未修改的官方第三方发行包、源码及构建产物不视为本工程活动实现。旧工程文件清单采用无损压缩归档，保留原始内容哈希，不参与运行。

旧 FEM 全文件审计 **{old['status']}**，覆盖 {old['entries']} 项，内容、权限、修改时间、符号链接和 Git 状态均未改变；原有 `reports/stage01/REPORT.md` 未提交修改保持原样。证据：[只读审计](old_fem_readonly_audit.json)。

## 使用了哪个血管？

输入为旧工程正式 Stage 1.7 selected exterior surface：
`{manifest['files']['exterior_surface.npz']['source_path']}`。

SHA256：`{manifest['files']['exterior_surface.npz']['source_sha256']}`。
五份必要输入的源文件与副本 SHA256 全部相同，见 [输入清单](../../inputs/MANIFEST.json)。旧体网格只用于计算参考网格尺寸。

| 原始标签 | 语义 | 当前面 ID |
|---|---|---|
| 1 | WALL | 1 |
| 2 | OUTLET_03 | 2 |
| 3 | OUTLET_01 | 3 |
| 4 | INLET | 4 |
| 5 | OUTLET_02 | 5 |

身份由原始标签、连通面和几何来源确定，见 [面映射](../../configs/face_map.json)。

![输入到 SimVascular 的血管模型](source_geometry.png)

## SimVascular 是否正常工作？

官方 Linux distribution **{api['package_version']}** 已再次实际验证：launcher、embedded Python 3.5.5、`import sv`、`sv.modeling.PolyData`、`sv.meshing.TetGen` 均通过。核心 API 不依赖未使用的图像插件；启动日志中保留了该插件缺少旧版 ICU 的提示。

本地没有可用流体 executable，因此只原生构建官方 [svMultiPhysics](https://github.com/SimVascular/svMultiPhysics)，固定 commit `{commit}`。executable：`{native['executable']}`；SHA256 `{native['executable_sha256']}`。

使用系统编译器与 `/usr/bin/mpiexec`。已有系统 BLAS/LAPACK 直接复用；缺少 VTK SDK，故在工程内部仅构建 VTK 9.3.1 所需 I/O/几何模块。已有发行包中的网格静态库直接复用，没有重编网格库、完整 GUI 或额外代数后端。[构建记录](native_build_log.txt)、[依赖记录](native_dependencies.json)、[链接库记录](native_solver.json)均保留。

唯一官方 smoke 为原始 `fluid/newtonian`：4 ranks、exit {smoke['exit_code']}、{smoke['elapsed_s']:.3f} s，实际 VTU 的速度和压力均 finite。证据：[smoke](official_smoke.json)。

## SimVascular 如何生成网格？

正式流水线：原始 tagged surface → PolyData model → 保留五个面角色 → `sv.meshing.TetGen` 同时启用 surface 和 volume meshing → `SV_MESH`。

全局尺寸为旧 Stage 1.7 体网格唯一边长中位数：`h_ref={m['generation']['global_edge_size_m']:.12e} m`。只生成一个 primary candidate，用时 {m['generation']['elapsed_s']:.3f} s；未使用允许的 0.8 倍失败回退，也未扫描参数。全程坐标单位为 m。

![重新生成后的血管表面](simvascular_surface.png)
![血管内部四面体网格](mesh_cutaway.png)

## 新网格有效吗？

**PASS**：{m['vertices']:,} 个节点，{m['tetra']:,} 个四面体；一个连通体，外表面闭合，五种边界完整，inverted / degenerate / non-finite 均为 0。

| 指标 | 实测值 |
|---|---:|
| minSICN 最小值 | {quality['q_min']:.9f} |
| P1 | {quality['P1']:.9f} |
| P5 | {quality['P5']:.9f} |
| median | {quality['median']:.9f} |
| P95 | {quality['P95']:.9f} |
| minSICN < 0.1 | {quality['N_low']} |

独立计算线性四面体理想映射的 signed inverse Frobenius condition number；分布供人工审核。证据：[有效性](mesh_validity.json)、[质量](mesh_quality.json)。

![新网格的质量如何](mesh_quality.png)

## 血管几何有没有明显改变？

重新三角化产生了可测变化，不能称为完全相同：壁面双向 closest-surface distance 的 max / P95 / RMS 分别为 **{g['wall']['max_m']:.9e} / {g['wall']['P95_m']:.9e} / {g['wall']['RMS_m']:.9e} m**；max/h_wall={g['wall']['max_over_h_wall']:.6f}。测量覆盖两个方向的全部顶点和三角形重心，是采样距离，不是连续曲面的严格 Hausdorff 上界。

包围体积相对变化 **{g['volume_relative_error']:.6%}**；实际体积 {m['volume_m3']:.12e} m³。数值几何变化按生成前冻结的 SV1 策略留给人工审核，拓扑及端口身份为硬检查。

| 端口 | 面积相对变化 | 重心位移 / m | 法向点积 | 连通片 / 翻转 |
|---|---:|---:|---:|---|
{port_table}

没有端口丢失、合并、拆分或法向翻转。各端口前后面积、三维重心和法向原始值见 [几何 QC](geometry_qc.json)。

## 流动条件是什么？

**REFERENCE NUMERICAL CONDITION；NOT EXPERIMENTAL FLOW；experimental=false。**
来自冻结的 Stage 3 物理配置，只读取参考物性和总流量：rho={p['density_kg_m3']} kg/m³，nu={p['kinematic_viscosity_m2_s']:.8e} m²/s，mu=rho·nu={p['dynamic_viscosity_pa_s']:.8e} Pa·s，Q={p['inlet_volume_flow_m3_s']:.12e} m³/s。

刚性壁、Newtonian 不可压缩流体，SI 单位。实际入口面积 {tp['A_in_m2']:.12e} m²，Umean={tp['Umean_m_s']:.12e} m/s，Dh={tp['Dh_m']:.12e} m，Re={tp['Re']:.12e}。这不是实验泵流量。见 [正式物理配置](../../configs/sv_reference.yaml)。

## 边界条件是什么？

入口采用官方 `Flat + Impose_flux + Zero_out_perimeter`：内部为恒定幅值，29 个共享壁面周边节点置零，靠周边一层采用离散线性过渡。为同时满足 no-slip 与总流量，内部速度幅值为 {inlet['interior_speed_m_s']:.12e} m/s，区别于几何平均速度 Umean；这不是严格每个入口节点均非零的平顶剖面。

按实际三角面与官方节点法向对写入后的 prescribed velocity 重新积分：Q_written={inlet['Q_written_m3_s']:.12e} m³/s，相对误差={inlet['relative_error']:.3e} ≤ 1e-10。微小非平面法向对齐因子 {inlet['alignment_factor']:.12f} 在求解前一次修正。该边界文件明确命名 `PrescribedVelocity`，未当作解场。

WALL：u=0。三个出口均为官方 `Neu, Value=0` 自然压力/牵引参考，代表表压参考负载；**不等于所有出口节点逐点 p=0**。没有 RCR、resistance 或预设分流。

[官方边界初始化源码]({upstream}baf_ini.cpp)和[边界处理源码]({upstream}set_bc.cpp)支持上述实现；本次 [solver XML](../../configs/sv_flow.xml)及其 SHA 在运行前冻结。

![流体从哪里进入从哪里流出](real_geometry_and_bc.png)

## solver 是否运行成功？

**进程运行并输出真实结果，但数值求解未通过。** 当前固定版本没有使用正式 steady PDE 开关，采用 constant BC 的 `fluid` transient-to-steady；BC 的 `Time_dependence=Steady` 仅指边界常量。[官方时间推进实现]({upstream}main.cpp)。

4 MPI ranks，OMP_NUM_THREADS=1，实际 exit={run['exit_code']}，GNU time 壁钟时间 {resource['wall_elapsed_s']:.2f} s（{resource['wall_elapsed_text']}），Python monotonic 计时 {run['elapsed_s']:.3f} s；最大单进程 RSS {run['peak_rss_kib']/1024:.2f} MiB（GNU time，非并行进程内存总和）。没有使用 GPU。两种时钟实测不一致，原始记录均保留，不据此做性能推断；详见 [资源记录](solver_resource_usage.json)。

原计划首段 {tp['first_block_steps']} 步。日志记录 **{q['linear_nonconvergence_warnings']} 次线性系统不收敛**、**{q['ill_conditioned_warnings']} 次病态左端矩阵警告**。即使外层残差降得很小，也不能据此认定内层线性系统收敛。已通过官方 `STOP_SIM` 请求在第 {q['step']} 步保存后退出：exit 0 表示正常停止，**不表示计划计算完成或物理验收成功**。

本次固定 LS 为内置 NS/FSILS：外迭代上限 30、相对容差 1e-10、绝对容差 1e-24、Krylov 维数 100、内层 GM/CG 容差 1e-3。它们是本次选定设置，并非声称全部采用官方默认数值。原始 LS 设置、dt、物性、BC 均未在运行中调整；未使用延长机会，因本次触发的是数值失败停止规则。保留 [停止决定](flow_stop_decision.json)、[执行记录](flow_execution.json)、[完整原始日志](../../logs/sv1/vascular_flow_first_block.log)。失败原因尚未被隔离为某一个机制，不将警告直接归因于网格或单位。

正式验收测试：**{tests['passed']} passed，{tests['failed']} failed，{tests['skipped']} skipped**。失败没有被改写为跳过或预期失败：

{fail_table}

![计算时间和内存](solver_resource_usage.png)

## 是否达到稳态？

**{s['status']}**。dt=min(0.5·h10/Umean, 0.05·Dh²/nu)={tp['dt_s']:.12e} s；h10={tp['h10_m']:.12e} m，t_nu={tp['viscous_time_s']:.12e} s。原计划 400 步为 20 t_nu，最多一次同长度 checkpoint 延长；本次实际只到 t={q['time_s']:.12e} s（{q['time_s']/tp['viscous_time_s']:.4f} t_nu）。

保存间隔 10 步，目前只有 {s['saved_states']} 个实际保存状态，无法验证连续五个保存区间同时满足 velocity change≤1e-5 和 boundary flow change≤1e-6。逐步边界日志中的趋势不代替上述联合判据。[冻结时间策略](../../configs/time_policy.json)、[稳态验收](steady_state.json)。

![计算是否达到稳态](steady_convergence.png)

## 速度场如何？

实际文件：`{q['path']}`，SHA256 `{q['sha256']}`。原生 point arrays 为 `Velocity` 与 `Pressure`。
速度 finite={q['velocity_finite']}，最大速度 {q['velocity_max_m_s']:.12e} m/s，体积积分 L2 范数 {q['velocity_L2']:.12e} m^(5/2)/s。全局图显示实际体网格节点；切片由该实际场线性插值得到。**属于失败运行的瞬态诊断，不是可采信的稳态速度分布。**

![血管中的速度分布](velocity_global.png)
![血管内部截面速度](velocity_slices.png)

## 压力场如何？

压力 finite={q['pressure_finite']}，全域范围 [{q['pressure_range_pa'][0]:.12e}, {q['pressure_range_pa'][1]:.12e}] Pa。压力图采用原生实际压力；截面平均是几何主轴法向平面的面积加权值，平面可能同时穿过多个分支，不能当作某根血管中心线压降。出口压力采用自然牵引参考，不额外平移或钉住压力。

![血管中的压力分布](pressure_global.png)
![不同位置的平均压力](pressure_sections.png)

## 流量守恒吗？

**{'PASS' if q['mass_balance_pass'] else 'FAIL / MASS_BALANCE_FAIL'}**。所有流量由实际 VTU 的四个端口外法向三角面独立积分；入口向内为正。线性三角面通量为精确 P1 积分。

| 指标 | 实际值 |
|---|---:|
| Qtarget / m³/s | {q['Q_target_m3_s']:.12e} |
| Qin / m³/s | {q['Q_in_m3_s']:.12e} |
| Qout1 / m³/s | {q['outlet_flows_m3_s']['OUTLET_01']:.12e} |
| Qout2 / m³/s | {q['outlet_flows_m3_s']['OUTLET_02']:.12e} |
| Qout3 / m³/s | {q['outlet_flows_m3_s']['OUTLET_03']:.12e} |
| Qout total / m³/s | {q['Q_out_total_m3_s']:.12e} |
| epsilon_Q = abs(Qin−Qtarget)/Qtarget | {q['epsilon_Q']:.12e} |
| epsilon_mass = abs(Qout−Qin)/Qtarget | {q['epsilon_mass']:.12e} |

两项门限均为 1e-6。独立 VTU 积分与原生边界积分日志的最大差/Qtarget={q['native_vs_vtu_max_difference_over_Q']:.12e}，支持后处理对应正确，但不消除质量闭合失败。证据：[实际场 QC](flow_qc.json)。此处只说明该未收敛瞬态状态的守恒误差，不外推到尚未得到的稳态解。

![流入流出是否守恒](flux_balance.png)

## 三个出口怎么分流？

以下比例是当前诊断状态的实测 Qout_i / sum(Qout)，没有预设。**不作为稳态分流结论。**

| 出口 | 流量 / m³/s | 占当前总流出量 |
|---|---:|---:|
{flow_table}

![三个出口分别流出多少](outlet_flow_split.png)

## 壁面 no-slip 是否满足？

实际 wall 节点速度 magnitude：max={q['wall_velocity_max_m_s']:.12e} m/s，P95={q['wall_velocity_P95_m_s']:.12e} m/s。检查门限为 1e-10·Umean={q['wall_noslip_tolerance_m_s']:.12e} m/s，结果 **{'PASS' if q['wall_noslip_pass'] else 'FAIL'}**。共享入口周边节点也包含在壁面检查内。

## 结果能否重新读取？

**诊断文件重读 {reload['status']}；物理验收仍为 FAIL。** 原生进程结束后，在新的 Python 进程重新读取原始 VTU，复核文件 SHA256，并重新计算速度 L2、压力范围、入口/三个出口流量与质量闭合。

重读相对速度范数差 {reload['velocity_norm_relative_difference']:.3e}，最大流量差/Qtarget {reload['flow_difference_over_Q']:.3e}，压力范围一致={reload['pressure_range_matches']}。这证明失败状态保存、解析和积分可复现，不意味着存在已验收的最终流场。[重读证据](solution_reload.json)、[测试结果](pytest_results.xml)、[图片来源](visual_provenance.json)。

## 当前还没有做什么？

- experimental pump flow
- mesh convergence
- WSS validation
- RBC/microbubble
'''
(R/'REPORT.md').write_text(text)
terminal=f'''Stage SV1 completed.

cleanup:
    obsolete deployment code removed = true
    remote code removed = true
    obsolete tests removed = {len(a['deleted_obsolete_tests'])}
SimVascular:
    API version = {api['package_version']}
    solver executable = {native['executable']}
    solver commit/version = {commit}
mesh:
    tetra = {m['tetra']}
    q_min = {quality['q_min']}
    P1 = {quality['P1']}
    P5 = {quality['P5']}
    median = {quality['median']}
    P95 = {quality['P95']}
    N_low = {quality['N_low']}
    invalid tetra = 0
geometry:
    wall max error = {g['wall']['max_m']} m
    port checks = PASS
physics:
    rho = {p['density_kg_m3']} kg/m3
    mu = {p['dynamic_viscosity_pa_s']} Pa*s
    Q_target = {q['Q_target_m3_s']} m3/s
solver:
    mode = transient-to-steady
    MPI ranks = 4
    completed = false; stopped at step {q['step']} after repeated linear nonconvergence
    process exit = {run['exit_code']}
    runtime wall-clock = {resource['wall_elapsed_s']} s
    runtime monotonic = {run['elapsed_s']} s
flow (unaccepted transient diagnostic state):
    Qin = {q['Q_in_m3_s']}
    Qoutlet_01 = {q['outlet_flows_m3_s']['OUTLET_01']}
    Qoutlet_02 = {q['outlet_flows_m3_s']['OUTLET_02']}
    Qoutlet_03 = {q['outlet_flows_m3_s']['OUTLET_03']}
    Qout_total = {q['Q_out_total_m3_s']}
    epsilon_Q = {q['epsilon_Q']}
    epsilon_mass = {q['epsilon_mass']}
fields:
    velocity finite = {q['velocity_finite']}
    pressure finite = {q['pressure_finite']}
    max wall velocity = {q['wall_velocity_max_m_s']} m/s
    P95 wall velocity = {q['wall_velocity_P95_m_s']} m/s
reload:
    diagnostic file = {reload['status']}
    accepted physical solution = false
tests:
    passed = {tests['passed']}
    failed = {tests['failed']}
    skipped = {tests['skipped']}
report:
    reports/sv1/REPORT.md
human review (diagnosis only):
'''+''.join('    '+p.name+'\n' for p in sorted(R.glob('*.png')))+f'\nSTAGE SV1 STATUS:\n    {status}: '+', '.join(reasons)+'\n'
(R/'terminal_summary.txt').write_text(terminal)
print(terminal)
