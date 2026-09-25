from pathlib import Path
import json,csv,hashlib
S=Path(__file__).resolve().parents[1];P=json.loads((S/'provenance/STAGE_PATHS.json').read_text());reference=json.loads((S/'provenance/REFERENCE_2D_SOURCE.json').read_text());get=lambda n:json.loads((S/'validation'/n).read_text())
cases={n:get(n+'_mpi1.json') for n in ['LOW','MEDIUM','LONG_TRANSPORT']};states={n:json.loads((S/'runs'/f'{n}_mpi1'/'RUN_STATE.json').read_text()) for n in cases};long=states['LONG_TRANSPORT'];top=get('WALL_SURFACE_TOPOLOGY_AUDIT.json');synthetic=get('SYNTHETIC_SURFACE_QUALIFICATION.json');sens=get('REAL_LONG_TOLERANCE_SENSITIVITY.json');mpi=get('MPI_DETERMINISM.json');pres=get('ORIGINAL_INPUT_PRESERVATION.json');fatal=get('FATAL_DIAGNOSTIC_QUALIFICATION.json')
for n in ['SYNTHETIC_SURFACE_QUALIFICATION.json','PRODUCTION_SURFACE_CPP_PYTHON.json','FEATURE_NORMAL_CONTINUITY_AND_SINGLE_NORMAL_LIMIT.json','REMOTE_NATIVE_SYNTHETIC_INDEPENDENT.json','REAL_LONG_TOLERANCE_SENSITIVITY.json','MPI_DETERMINISM.json','RAW_PHYSICS_SOURCE_IDENTITY.json','ORIGINAL_INPUT_PRESERVATION.json','FATAL_DIAGNOSTIC_QUALIFICATION.json','NATIVE_INHERITED_CHECKS.json']:
 assert get(n)['status']=='PASS',n
assert all(c['status']=='PASS' for c in cases.values()) and long['time_s']>.494122765648
status=dict(REFERENCE_2D_REVIEW_STATUS='PASS',WALL_SURFACE_TOPOLOGY_STATUS='PASS',OLD_AMBIGUITY_ROOT_CAUSE_STATUS='CLASSIFIED',OLD_AMBIGUITY_ROOT_CAUSE='MESH_DISCRETIZATION_TIE',SMOOTH_SURFACE_CLUSTERING_STATUS='PASS',PSEUDONORMAL_STATUS='PASS',TRUE_MULTI_SURFACE_STATUS='PASS',TRUE_MULTI_SURFACE_SCOPE='SYNTHETIC_WEDGES_AND_THREE_SURFACE_CORNER; not triggered in formal real runs',TRIANGLE_TESSELLATION_INVARIANCE_STATUS='PASS',TOLERANCE_SENSITIVITY_STATUS='PASS',RAW_RESISTANCE_SOLVER_STATUS='PASS',RK2_STATUS='PASS',HARD_WALL_NONPENETRATION_STATUS='PASS',PAIR_NONOVERLAP_STATUS='PASS',PAIR_NONOVERLAP_SCOPE='Accepted states and swept paths satisfy inherited -1e-12 m tolerance; terminal pair rejection is retained.',PARTICLE_ACCOUNTING_STATUS='PASS',ADAPTIVE_INJECTION_STATUS='PASS',MPI_STATUS='PASS',OLD_LOW_STALL_RESOLVED='YES',OLD_MEDIUM_STALL_RESOLVED='YES',LONG_OLD_AMBIGUITY_RESOLVED='YES',LONG_FINAL_TIME=long['time_s'],LONG_COMPLETED_0P8S='NO',LONG_STOP_REASON=long['reason'],NATURAL_OUTLET_STATUS='NOT_OBSERVED',WALL_HYDRODYNAMICS_STATUS='NOT_IMPLEMENTED',FLOW_PHYSICS_STATUS='ENGINEERING_TRANSIENT_FIELD_ONLY',PRODUCTION_BUBBLE_CONCENTRATION_STATUS='UNSPECIFIED',SURFACE_CONTACT_MODULE_STATUS='PASS',WORKFLOW_TECHNICAL_STATUS='PARTIAL',REAL_SCIENTIFIC_PRODUCTION_READY='NO',RMBW_USED='NO',CF2003_USED='NO',FALADE_BRENNER_USED='NO',RALLABANDI_USED='NO',HIGDON_MULDOWNEY_USED='NO',ADHESION_USED='NO',RBC_USED='NO',GPU_OPTIMIZATION_USED='NO',STL_MODIFIED='NO',SOURCE_MODIFIED_OUTSIDE_NEW_STAGE='NO',GIT_COMMIT='NO',disclaimer='NOT EXPERIMENTAL CONCENTRATION')
(S/'WORKFLOW_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
rows=[]
for n,c in cases.items():
 st=states[n];rows.append(f"| {n} | {st['time_s']:.15g} | {st['accepted_steps']} | {st['reason']} | {st['source_drawn']}/{st['size_rejected']}/{st['admitted']}/{st['active']} | {c['min_wall_gap_m']:.12g} | {c['min_all_pair_swept_gap_m']:.12g} |")
checks=[]
for n,c in cases.items():checks.append(f"| {n} | {c['wall_constraints']['wall_stage_records']} | {c['independent_dense_solves']} | {c['max_dense_reference_scaled_error_m_s']:.6g} | {c['wall_constraints']['position_corrections_accepted']} | {c['wall_constraints']['max_independent_normal_error']:.6g} |")
statusrows='\n'.join(f'| `{k}` | `{v}` |' for k,v in status.items() if k.endswith('STATUS') or k in ['OLD_AMBIGUITY_ROOT_CAUSE','OLD_LOW_STALL_RESOLVED','LONG_OLD_AMBIGUITY_RESOLVED','LONG_FINAL_TIME','LONG_COMPLETED_0P8S','REAL_SCIENTIFIC_PRODUCTION_READY'])
report=f'''CONTINUOUS_SURFACE_CONTACT_3D_V0 — 开发与独立验证报告

1. 上一阶段 LONG 为什么停？旧算法把最近三角形的数值并列当作法向歧义，在 0.49412276564811575 s 的已接受状态之后停止。
2. 是网格造成的，还是真实多个壁面？是 **MESH_DISCRETIZATION_TIE**。两片 BIOLOGICAL_CORE 三角形共享流形边，二面角仅 2.22694°，属于同一个光滑局部壁面；不是区域接缝，也不是非流形边。
3. 2D simulator 给了什么参考？参考了连续边界来源、相邻离散单元与真实多壁的区分、移动前完整路径预测、失败后细分物理时间。没有复制 2D wall mobility、adhesion、RBC 或大型模块。
4. 3D 怎样组合 triangles？按共享流形边、局部定向、二面角和最近特征的局部 one-ring 分簇；区域标签只作诊断。三角形 ID 根据顶点坐标规范化。原 STL 完全不改。
5. 什么时候用一个 normal？候选属于一个光滑簇时只施加一个约束。面内用面法向；边/顶点保留加权 mesh pseudonormal 作定向，同时用球心到共同最近特征的距离梯度作为有限半径接触法向。两者明确区分，不是两个 wall constraints。
6. 什么时候用多个 normal？真正不同局部曲面保持独立法向，解最小欧氏改动的平移速度投影；枚举最多 rank 3 的 active set，不跨 sharp feature 平均。
7. 有没有改变 Rq=b？没有。仍是 R_bulk + R_BB_excess，twist pending；原系数、PCG、pair hard constraint、Frozen Flow 核心源码哈希不变。
8. 有没有加入 wall force？没有，也没有加入 wall resistance、wall lubrication 或 wall mobility。
9. LONG 有没有越过 0.494 s？有。MPI1/MPI4 都达到 **{long['time_s']:.16g} s**，不再因 triangle tie 停止。
10. 有没有达到 0.8 s？没有。原有 pair swept safety 无法接受后续试算，明确停止于 **PAIR_SAFETY_FAILED**，未修改 pair physics 或放宽容差。
11. 有没有自然 outlet？没有。端口与删除/迁移 fixture 通过，不把人工 fixture 当作自然流出事件。
12. MPI1/MPI4 是否一致？一致。source draw、IDs、候选/簇 ID、分类、V_used、轨迹、端口事件和物理终态完全一致；仅实际 owner rank、cross-rank 计数及 mpi_ranks 按并行分区不同。

**结论：surface-contact 模块 PASS；整体 workflow PARTIAL。** 全部是几何/运动学约束，WALL_HYDRODYNAMICS_STATUS = NOT_IMPLEMENTED。生产浓度仍 UNSPECIFIED，流场仍 ENGINEERING_TRANSIENT_FIELD_ONLY，不能作为真实科研生产模拟已就绪的声明。

**实际运行结果**

所有正式 case 从 t=0 开始，各运行 MPI1/MPI4。以下为独立回读 MPI1 的值；MPI4 的所有物理记录逐项相同。

| Case | 终止时间 s | 接受步数 | 终止原因 | draw/reject/admit/active | 最小壁面间隙 m | 最小 pair 扫掠间隙 m |
|---|---:|---:|---|---|---:|---:|
{chr(10).join(rows)}

MEDIUM 和 LONG 的最小 pair gap 略为负值，分别约 −0.9974 pm 和 −0.99996 pm，均在继承的 −1 pm 数值安全容差内；这里不把它们写成“严格非负”。壁面独立复核容差保持 2e−14 m，运行时仍使用原严格扫掠判定及 1e−10 m margin。后续失败试算不提交、不推进时间、不静默删粒子。

**旧故障的可复现证据**

在修改正式算法前，先用旧 stage 的只读 library，从 step 18134 的两粒子精确已接受状态重放。重现 `STOP_WALL_NORMAL_AMBIGUITY`：stage 1、particle 2、dt=1e−4 s、radius=9.208782268416533e−7 m。见 [OLD_LONG_FAILURE_REPLAY.json](validation/OLD_LONG_FAILURE_REPLAY.json) 和 [完整根因分类](validation/OLD_AMBIGUITY_ROOT_CAUSE.json)。

旧顺序 triangle 28422 / 28261 对应规范 ID 11572 / 11574。它们共享边顶点 5745 / 5908，均为 BIOLOGICAL_CORE；两个最近点约相差 1.8e−13 m，但最近距离只相差约 2.4e−20 m。九组预设数值参数全部归为一个 smooth patch。ID 只用于稳定标签，不用于任选最终法向。

**有限半径法向的必要区别**

第一版直接把加权面法向用于有限半径球的边接触，曾使 LOW 在约 0.2363 s 进入微小步长：mesh pseudonormal 与距离梯度相差 1.32755°，投影后真实 gap 仍以 2.29208e−6 m/s 减小。这是一次不合格开发试跑，已停止并完整保存在 `runs/DEV_V1_LOW_INTERRUPTED` 和 `provenance/development_v1`，不计入通过结果。

最终实现保留 area-weighted edge / angle-weighted vertex mesh pseudonormal，用它判断局部曲面定向；单一光滑特征的有限半径约束使用 `normalize(X − weighted_closest_feature_point)`。这个 offset-distance gradient 才与球心 gap 的一阶变化一致。它不是任意挑 triangle，也不是 wall force。真实多个曲面继续保持各自法向。该修正是对推荐 pseudonormal 配方的明确补充，见 [设计对照](design/2D_TO_3D_CONTACT_TRANSLATION.md)、[开发回归诊断](validation/LOW_DEVELOPMENT_REGRESSION_DIAGNOSIS.json) 和 [最终合同](contracts/SURFACE_CONTACT_CONTRACT.json)。

不声称未修改的 STL 是解析光滑曲面：辅助 mesh pseudonormal 本身按特征定义，可在面/边间有限变化。最终 offset contact normal 在测试的凸面、边、顶点 Voronoi 区域连续；平面跨特征法向和速度完全连续。曲面轨迹通过细化物理步长收敛，所有接受路径仍以原 triangle swept certificate 为真值。

**拓扑、参数和 2D 来源**

WALL_ONLY：67,071 triangles、33,633 exact-coordinate vertices、100,706 edges，其中 100,507 manifold、199 boundary、0 nonmanifold，1 个连通分量，0 定向冲突，0 退化三角形。199 条边界边均对应被排除的端口帽。区域接缝按几何连续性处理，标签重命名不改变结果。证据：[拓扑](validation/WALL_SURFACE_TOPOLOGY_AUDIT.json)、[端口边界](validation/WALL_BOUNDARY_PORT_AUDIT.json)、[区域接缝](validation/REGION_SEAM_AUDIT.json)、[流体侧法向](validation/FLUID_SIDE_NORMAL_AUDIT.json)。

默认 smooth_dihedral_threshold=15°；敏感性为 10°/15°/25°。tie distance tolerance=`factor × 64 × eps × max(norm(X),1e−6 m)`，factor=0.5/1/2；barycentric feature tolerance=1e−9。这些都是数值网格参数，不是 wall physics 参数，事先记录于设计。真实 LONG 九次完整 0→0.51 s 积分的候选数、簇数、模式、使用速度和轨迹完全相同，最大位置差和速度差均为零。没有靠极窄容差掩盖问题。

2D reference 固定于 commit `8c60b2620509dda748df9e958ff3d1251fdc2ab9`，分支 `sync/ulm-microbubble-traj-gen-2d-20260917`。独立 bare clone + git show 读取，不 checkout 用户 3D 仓库；读取文件和 SHA256 见 [REFERENCE_2D_SOURCE.json](provenance/REFERENCE_2D_SOURCE.json)。REFERENCE_ONLY=YES，CODE_COPIED_DIRECTLY=NO。2D 参考中的真实多壁仍会显式停错，其 mobility-space correction 没有搬入；本阶段重新实现 3D 欧氏 multi-normal projection。

**验证范围与实测误差**

合成资格检查包含 9,072 个 C++/独立 NumPy 查询及 QP 对照：随机平面、三种平面三角化、face/edge/vertex、三角形重排、三角化圆柱/球面、30°/60°/90°/120°/150° 楔角、三壁角点、故意非流形边。包含 216 个 true multi-surface 查询和 36 个预期 FAIL_NONMANIFOLD。另有 3,003 个连续特征路径查询。远端同源 C++ 实际运行 14 组轨迹，再由 NumPy/VTK 回读扫掠安全；平面路径跨原点顶点和多个三角形。

圆柱/球面采用 dt=1e−4、5e−5、2.5e−5、1.25e−5 s。相对最细档的圆柱端点误差为 7.0223e−14 / 1.6800e−14 / 3.3678e−15 m；球面为 6.1610e−14 / 3.4159e−14 / 1.7449e−15 m，均递减。平面 single-normal 极限与旧公式的实测相对误差为 0。形式化 QP 校验将几何浮点误差和优化误差分别测量；独立重建法向，再独立解 KKT，检查全部约束、目标值和 tight active set。

| Case | 独立壁面阶段记录 | 独立 dense raw solves | 最大 scaled raw error m/s | 接受位置修正 | 最大独立法向差 |
|---|---:|---:|---:|---:|---:|
{chr(10).join(checks)}

两次 RK2 的位置重构误差均为 0；Omega 在全部正式阶段的 raw/used 输出逐位一致。独立 validator 重建 source draw、FIFO、flux integral、inlet admissibility、LAMMPS 原生粒子数、全部接受轨迹/扫掠、pair non-overlap、端口事件和 near-wall accounting。原 Rq=b 系数/PCG/采样器源码及其 7 个稳定组件哈希不变；共检查 1,333 个独立 dense raw solves。

LOW/MEDIUM 不需要位置修正。LONG 有 104 次接受修正，占阶段记录 0.634%；总修正 2.43465e−10 m，占总粒子路径 2.84596e−6，最大单次 2.21051e−11 m < margin/4=2.5e−11 m。仅在真实 wall retry 之后、单簇、微小残差时允许一次修正，并重新认证整段路径；multi-surface 不作位置修正。修正不是主要运动算法。见 [修正使用审计](validation/BOUNDED_CORRECTION_USAGE.json)。

将 LONG 的 max_retry 从正式 24 临时降为 1 的独立诊断 fixture，会在 step 3712 停止；`WALL_CONTACT_FAILURE.json` 在终态文件之前已 flush，保存半径、位置、raw/used velocity、candidate 坐标/最近点/面法向/邻接/区域、簇和法向集、dt、RK2 stage、原始及投影路径最小 gap。其投影和扫掠值均已独立重算。正式配置能通过该步并继续。见 [失败诊断资格](validation/FATAL_DIAGNOSTIC_QUALIFICATION.json) 和 [实际失败文件](runs/DIAGNOSTIC_LONG_RETRY1_mpi1/WALL_CONTACT_FAILURE.json)。

主要验证入口：[合成资格](validation/SYNTHETIC_SURFACE_QUALIFICATION.json)、[特征连续性](validation/FEATURE_NORMAL_CONTINUITY_AND_SINGLE_NORMAL_LIMIT.json)、[远端合成轨迹独立验证](validation/REMOTE_NATIVE_SYNTHETIC_INDEPENDENT.json)、[生产几何 C++/Python 对照](validation/PRODUCTION_SURFACE_CPP_PYTHON.json)、[真实容差敏感性](validation/REAL_LONG_TOLERANCE_SENSITIVITY.json)、[正式回读](validation/INDEPENDENT_WORKFLOW_VALIDATION.json)、[MPI 对照](validation/MPI_DETERMINISM.json)、[原生 injection/ports/lifecycle](validation/NATIVE_INHERITED_CHECKS.json)。

**图**

![TRIANGLES ≠ PHYSICAL WALLS](visualization/continuous_surface_contact_3d_scheme.png)

![LONG old versus new](visualization/LONG_old_vs_new_trajectory.png)

其余 11 张图均在 `visualization/`：旧故障局部几何、候选与簇、smooth-edge pseudonormal、true multi-surface、二面角分布、簇数分布、LONG gap、约束模式、dt 历史、平面重三角化不变性、九组容差敏感性。合计 13 张，均标注 TEST_ONLY 浓度免责声明。

**保护、同步与限制**

本地复核 {pres['local_protected_files']} 个历史/原始文件，含前一 wallslide stage 全部 292 文件与 bcflux 全部 545 文件；原 STL、Frozen Flow、SonoVue 和 authoritative Q 来源不变。两个 3D 仓库的 branch、commit、完整 dirty/untracked 状态与接手前一致。无 git add/commit/push/reset/clean/merge/rebase。最终 remote/local 全文件 SHA256、远端 1,036 个历史输入及 9 个原 CPU build/library 文件的复核结果见 [同步收据](validation/REMOTE_LOCAL_SHA256_CHECK.json)；文件清单见 [SYNC_MANIFEST.json](provenance/SYNC_MANIFEST.json)。

当前新限制仍是 pair safety：wall 投影后继续执行原 pair swept check，无法接受共同试算时停止，没有新增 wall–pair 联合动力学求解器。无自然出口、无 0.8 s 完成声明。无 RMBW、CF2003 production correction、Falade–Brenner、Rallabandi、Higdon–Muldowney、wall lubrication、adhesion、RBC、GPU 优化或 large-N production。

本地 stage：`{S}`  
远端 results：`{P['remote_stage']}`  
远端 work：`{P['remote_work']}`  
2D 本地 reference：`/home/lzy/projects/compre_output/reference/ulm_microbubble_2d_20260917_003848`  
2D 远端 reference：`/workspace/reference/ulm_microbubble_2d_20260917_003848`

**最终状态**

| 字段 | 值 |
|---|---|
{statusrows}

NOT EXPERIMENTAL CONCENTRATION。REAL_SCIENTIFIC_PRODUCTION_READY = NO。
'''
name='MICROBUBBLE_VASCULAR_WORKFLOW_V0_SURFACE_CONTACT_REPORT.md';(S/name).write_text(report);(S/'report/INDEX.md').write_text(f'[完整报告](../{name})\n\n[机器可读状态](../WORKFLOW_STATUS.json)\n')
(S/'README.md').write_text(f'# CONTINUOUS_SURFACE_CONTACT_3D_V0\n\nSurface contact module PASS; overall workflow PARTIAL. LONG = {long["time_s"]} s, PAIR_SAFETY_FAILED.\n\n[报告]({name}) · [状态](WORKFLOW_STATUS.json) · [设计](design/2D_TO_3D_CONTACT_TRANSLATION.md)\n\nNOT EXPERIMENTAL CONCENTRATION. Wall hydrodynamics NOT_IMPLEMENTED.\n')
executed=[]
for p in sorted((S/'runs').glob('*/RUN_STATE.json')):
 d=p.parent;j=json.loads(p.read_text());scope='FORMAL' if d.name in [n+f'_mpi{k}' for n in cases for k in [1,4]] else 'SENSITIVITY' if d.name.startswith('SENS') else 'DIAGNOSTIC_RETRY_BUDGET_FIXTURE';executed.append(dict(name=d.name,scope=scope,time=j['time_s'],reason=j['reason'],config=str(d/'case.cfg')))
(S/'configs/EXECUTED_CASES.json').write_text(json.dumps(executed,indent=2)+'\n');(S/'contracts/CASE_INDEX.json').write_text(json.dumps(dict(formal=['LOW','MEDIUM','LONG_TRANSPORT'],formal_mpi_ranks=[1,4],sensitivity_runs=9,diagnostic_runs=2,development_interrupted='DEV_V1_LOW_INTERRUPTED; excluded from all acceptance claims',all_test_only=True),indent=2)+'\n')
p=S/'contracts/WORKFLOW_CONTRACT.json';j=json.loads(p.read_text());j['status']='PARTIAL';j['stop_reason']='PAIR_SAFETY_FAILED in MEDIUM and LONG; surface module qualified';p.write_text(json.dumps(j,indent=2)+'\n')
terminal=['MICROBUBBLE_VASCULAR_SURFACE_CONTACT_3D_V0_COMPLETE']+[f'{k}: {v}' for k,v in status.items() if k!='disclaimer']+[f'REPORT: {S/name}',f'LOCAL_STAGE: {S}',f'REMOTE_STAGE: {P["remote_stage"]}','REFERENCE_2D_PATH: /home/lzy/projects/compre_output/reference/ulm_microbubble_2d_20260917_003848']
(S/'report/FINAL_TERMINAL_STATUS.txt').write_text('\n'.join(terminal)+'\n');print('REPORT_READY',S/name)
