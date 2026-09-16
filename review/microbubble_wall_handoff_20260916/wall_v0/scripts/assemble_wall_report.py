from pathlib import Path
import json,hashlib,csv,itertools,sys,shutil
import numpy as np
R=Path(__file__).resolve().parents[1];P=json.loads((R/'provenance/TASK_PATHS.json').read_text());W=Path(P['local_work'])
def load(p):return json.loads((R/p).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=load('validation/OFFLINE_INDEPENDENT_FINALIZER.json');run=load('validation/RUNTIME_INDEPENDENT_FINALIZER.json');supp=load('validation/SUPPLEMENTAL_INDEPENDENT_FINALIZER.json');mpi=load('validation/MPI_WALL_MODEL.json');tab=load('contracts/RMBW_WALL_TABLE_CONTRACT.json');cwall=load('validation/C_WALL_CONVERGENCE.json');geom=load('validation/REAL_GEOMETRY_AUDIT.json');handoff=load('validation/FAR_FIELD_HANDOFF_AUDIT.json');contact=load('validation/CONTACT_INDEPENDENT_FINALIZER.json');integrity=load('validation/REMOTE_TO_WSL_INTEGRITY.json');vis=load('VISUALIZATION_PROVENANCE.json')
assert all(x['status']=='PASS' for x in [a,run,supp,mpi,integrity])
for name,digest in load('provenance/RUNTIME_SOURCE_SHA256.json').items():assert sha(R/name)==digest,name
selected=[]
for d in sorted((R/'cases').glob('*')):
 if not (d/'RUNTIME_SELECTION.json').exists():continue
 x=json.loads((d/'RUNTIME_SELECTION.json').read_text());assert x['selected_C_wall']==cwall['selected_C_wall'];assert (d/'case.selected.cfg').read_text()==(d/'case.cfg').read_text().replace('c_wall 0.4\n','c_wall 0.2\n');assert sha(d/'case.selected.cfg')==x['selected_config_sha256'];assert sha(d/'CASE_CONTRACT.json')==x['original_case_contract_sha256'];receipt=json.loads((d/'RUN_RECEIPT.json').read_text());assert 'case.selected.cfg' in receipt['command'];selected.append(d.name)
(R/'validation/RUNTIME_FROZEN_SELECTION_VERIFICATION.json').write_text(json.dumps({'status':'PASS','cases':selected,'original_contracts_unchanged':True,'only_derived_C_wall_selection_applied':.2},indent=2)+'\n')
old=load('provenance/PREVIOUS_REAL_8_CASE_CONTRACT.json');x=np.array(old['initial_positions_m']);rad=np.array(old['radii_m']);jpair=min(np.linalg.norm(x[i]-x[j])-rad[i]-rad[j] for i,j in itertools.combinations(range(8),2));jgap=run['checks']['J_real_8']['min_initial_gap_m'];jstate=load('cases/J_real_8/RUN_STATE.json');sem={'status':'BLOCKED_LOCAL_PLANE_VALIDITY','accepted_steps':0,'initial_state_only':True,'velocity_columns':'zero placeholders because no solve was performed; not physical results','wall_active_column':'candidate h/a range flag, not a successfully applied wall model','constraint_active_column':'no solve performed','raw_min_pair_gap_sentinel':1e100,'independent_initial_pair_gap_m':jpair,'trajectory_comparison_at_0p02s':'UNAVAILABLE'};(R/'cases/J_real_8/INITIAL_STATE_SEMANTICS.json').write_text(json.dumps(sem,indent=2)+'\n')
# Aggregate only actual timestep-model observations, keeping intentional controls labelled.
events=[]
for d in sorted((R/'cases').glob('*')):
 if not (d/'WALL_INTERACTION_EVENTS.csv').exists():continue
 state=json.loads((d/'RUN_STATE.json').read_text())
 if state.get('accepted_steps',0)==0:continue
 with (d/'WALL_INTERACTION_EVENTS.csv').open() as f:
  for row in csv.DictReader(f):events.append({'case':d.name,'known_penetrating_control':state['known_penetrating_control'],**row})
with (R/'WALL_INTERACTION_EVENTS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(events[0]));w.writeheader();w.writerows(events)
(R/'validation/WALL_EVENT_AGGREGATION.json').write_text(json.dumps({'rows':len(events),'Case_J_excluded':'blocked before first solve; raw initial diagnostics retained separately','controls_explicitly_labelled':True},indent=2)+'\n')
finalizer={'status':'PASS','scientific_stage_status':'BLOCKED_LOCAL_PLANE_VALIDITY','definition':'Independent math, runtime evidence, gate enforcement, MPI identity and input integrity validated; not a PASS of the blocked real trajectories','parts':['OFFLINE_INDEPENDENT_FINALIZER.json','CONTACT_INDEPENDENT_FINALIZER.json','RUNTIME_INDEPENDENT_FINALIZER.json','SUPPLEMENTAL_INDEPENDENT_FINALIZER.json'],'source_sha256':sha(R/'src/finalize_wall_hydrodynamics_v0.py'),'production_CPP_called':False};(R/'validation/INDEPENDENT_FINALIZER.json').write_text(json.dumps(finalizer,indent=2)+'\n')
normalgate=load('validation/FLAT_WALL_MODE_AUDIT.json');tablegate=load('validation/TABLE_INTERPOLATION_AUDIT.json');querymax=max(a['geometry']['REAL']['max_distance_error_m'],a['geometry']['CURVED']['max_distance_error_m']);states=[json.loads(d.read_text()) for d in (R/'cases').glob('*/RUN_STATE.json') if d.parent.name!='I_real_single'];safe=[s for s in states if not s['known_penetrating_control']];constraint_count=sum(s['wall_constraint_activations'] for s in safe)
S={
'MICROBUBBLE_WALL_HYDRODYNAMICS_V0':'BLOCKED_LOCAL_PLANE_VALIDITY',
'CORE_IMPLEMENTATION':'PASS_WITH_REFERENCE_LIMITATIONS',
'MODEL':'RMBW_LOOKUP_LOCAL_PLANE_RESISTANCE_V0','PRIMARY_REFERENCE':'RMBW_LUBRICATION','RMBW_COMMIT':tab['RMBW_commit'],'RMBW_RUNTIME_DEPENDENCY':'NO',
'RMBW_TABLE_FILE':str(R/tab['table_file']),'RMBW_TABLE_SHA256':tab['table_sha256'],'TABLE_GRID_POINTS':tab['grid_points'],'TABLE_EPSILON_MIN':.001,'TABLE_EPSILON_MAX':20,'NEAR_FIELD_HIGH_CONFIDENCE':'0.001–0.2','QUALIFIED_REFERENCE_RANGE':'0.001–5','LIMITED_REFERENCE_RANGE':'5–20',
'TABLE_MAX_INTERPOLATION_ACTION_ERROR':tablegate['max_qualified_action_error'],'TABLE_MAX_RECIPROCITY_ERROR':a['max_reciprocity_error'],'TABLE_MIN_SCALED_EIGENVALUE':a['minimum_mobility_scaled_eigenvalue'],'TABLE_MIN_SCALED_EIGENVALUE_MATRIX':'MOBILITY; full grid','BULK_DOUBLE_COUNTING':'NO',
'WALL_QUERY':'PASS','MAX_WALL_DISTANCE_ERROR_M':querymax,'LOCAL_PLANE_MODEL':'ON','LOCAL_PLANE_RMS_GATE':'0.10 a','LOCAL_PLANE_NORMAL_SPREAD_GATE_DEG':20,'LOCAL_PLANE_VALIDITY':'4/10000 valid; Case I no valid d50 near-wall start; Case J 8/8 invalid',
'WALL_RESISTANCE':'ON','NORMAL_WALL_RESISTANCE':'ON','TANGENTIAL_WALL_RESISTANCE':'ON','ROTATIONAL_WALL_RESISTANCE':'ON_WITH_REFERENCE_LIMITATIONS','TRANSLATION_ROTATION_WALL_COUPLING':'ON_WITH_REFERENCE_LIMITATIONS','AMBIENT_STRAIN_INDUCED_WALL_FORCE':'PENDING',
'HARD_WALL_NONOVERLAP':'ON','WALL_OVERLAP_TOL_M':1e-12,'C_WALL_TESTED':'0.4,0.2,0.1','SELECTED_C_WALL':.2,
'CASE_A_TABLE_LOOKUP':'PASS','CASE_B_NORMAL_WALL':'PASS','CASE_C_TANGENTIAL_WALL':'PASS','CASE_D_TR_COUPLING':'PASS','CASE_E_ROTATION':'PASS_WITH_REFERENCE_LIMITATION','CASE_F_HARD_WALL':'PASS','CASE_G_CURVED_QUERY':'PASS','CASE_H_REAL_STATIC_GEOMETRY':'PASS','CASE_I_REAL_SINGLE':'BLOCKED_NO_VALID_NEAR_WALL_START','CASE_J_REAL_8BUBBLE':'BLOCKED_LOCAL_PLANE_VALIDITY',
'CASE_J_PHYSICAL_TIME_S':0,'CASE_J_TARGET_TIME_S':.02,'CASE_J_MIN_WALL_GAP_UM':jgap*1e6,'CASE_J_MIN_BUBBLE_BUBBLE_GAP_UM':jpair*1e6,'CASE_J_GAP_STATISTICS_SCOPE':'INITIAL_STATE_ONLY; no accepted trajectory',
'WALL_CONSTRAINT_ACTIVATIONS':constraint_count,'WALL_CONSTRAINT_ACTIVATIONS_SCOPE':'all accepted safe synthetic runs, including three Cwall candidates; no real activations',
'NONFINITE_COUNT':sum(s['nonfinite_count'] for s in safe),'LOST_ATOMS':sum(s['lost_atoms'] for s in safe),'INVALID_QUERY_COUNT':sum(s['invalid_query_count'] for s in safe),'ACCEPTED_WALL_OVERLAP_COUNT':sum(s['accepted_wall_overlap_count'] for s in safe),'ACCEPTED_BUBBLE_OVERLAP_COUNT':sum(s['accepted_pair_overlap_count'] for s in safe),'SAFETY_COUNT_SCOPE':'hard-wall enabled runs; intentionally penetrating controls separately reported',
'MPI_WALL_MODEL':'PASS','KOKKOS_COMPATIBILITY_SMOKE':'PASS','GPU_PERFORMANCE_READY':'NO','INDEPENDENT_FINALIZER':'PASS',
}
for p in sorted((R/'visualization').glob('VIS_*.png')):S[p.stem]='GENERATED_BLOCKED_CASE_ANNOTATION' if p.stem in ['VIS_REAL_SINGLE_BUBBLE_WALL','VIS_REAL_SINGLE_BUBBLE_GAP','VIS_REAL_8BUBBLE_WALL_COMPARISON','VIS_REAL_WALL_ACTIVITY'] else 'GENERATED'
S.update(VISUALIZATION_GENERATION='PASS',HUMAN_VISUAL_REVIEW='PENDING',PALABOS_BASELINE_MODIFIED='NO',LAMMPS_CORE_MODIFIED='NO',PASSIVE_TRANSPORT_MODIFIED='NO',BUBBLE_BUBBLE_BASELINE_MODIFIED='NO',SONOVUE_SAMPLER_MODIFIED='NO',RMBW_SOURCE_MODIFIED='NO',FULL_EXACT_WALL_HYDRODYNAMICS='PENDING',MICROBUBBLE_ADHESION='PENDING',BUBBLE_BUBBLE_TWIST='PENDING',BUOYANCY='OFF',LIFT='OFF',RBC='OFF',REMOTE_TO_WSL_INTEGRITY='PASS',REMOTE_RESULT_DIR=P['remote_result'],LOCAL_RESULT_DIR=str(R),REPORT=str(R/'MICROBUBBLE_WALL_HYDRODYNAMICS_V0_REPORT.md'),NEXT_STAGE='USER_REVIEW_BEFORE_MICROBUBBLE_ADHESION_V0')
(R/'FINAL_TERMINAL_SUMMARY.txt').write_text(''.join(str(k)+' = '+str(v)+'\n' for k,v in S.items()));(R/'FINAL_STATUS.json').write_text(json.dumps(S,indent=2)+'\n')
case_table='| 核查 | 实测结论 |\n|---|---|\n'+''.join('| '+k+' | '+str(v)+' |\n' for k,v in S.items() if k.startswith('CASE_') and k not in ['CASE_J_GAP_STATISTICS_SCOPE'])
cwtable='| C_wall | 最大位置差 (m) | 峰值归一角速度差 | 积分转角差 (rad) | 判定 |\n|---|---:|---:|---:|---|\n'+''.join(f"| {r['C_wall']} | {r['position_error_m']:.8g} | {r['omega_relative_peak_error']:.8g} | {r['integrated_rotation_error_rad']:.8g} | {r['status']} |\n" for r in cwall['rows'])
report=f'''RMBW 先离线告诉我们：一个刚性球在不同离墙距离下有多难向墙靠近、沿墙滑动和旋转。

正式时域仿真不会每一步运行 RMBW，而是根据当前 h/a 快速查一张提前生成并验证好的阻力表。

**本轮总状态：BLOCKED_LOCAL_PLANE_VALIDITY。** 查表、C++ 阻力装配、统一硬约束、连续扫掠、独立核查及 MPI/Kokkos 技术测试已经实现并通过；真实血管 Case I 未找到符合冻结门限的 d50 近壁起点，Case J 原 8 个起点全部违反局部平面门限，所以本阶段不满足完整 PASS 条件。没有自动进入 adhesion 阶段。

**数据与科学定义。** 主参考为 [RigidMultiblobsWall 官方仓库](https://github.com/stochasticHydroTools/RigidMultiblobsWall)的 Lubrication 实现，固定 commit `{tab['RMBW_commit']}`。上游明确为球提供润滑修正路径；本轮只离线调用该既定路径，保留原始矩阵输出，不把上游代码并入生产 timestep。冻结参考核查、量纲/符号约定及来源位于 `provenance/reference_audit/`。RMBW 的 GPL-3.0 来源/版权及独立 Eigen 头文件许可记录保留；生产驱动没有运行时 RMBW、Python 或 subprocess 依赖。

表文件为 [{tab['table_file']}]({tab['table_file']})，SHA256 `{tab['table_sha256']}`。三种半径直接沿用 SonoVue d10/d50/d90；流体黏度为 0.001 Pa·s。真实场 SHA `{sha(R/'fields/FROZEN_FLOW_FIELD_V0.h5')}`，STL SHA `{sha(R/'provenance/closed_geometry_m.stl')}`。原 8 泡的尺寸、ID、起点和目标时长 0.02 s 均保持原值。真实流场仍为 **ENGINEERING_TRANSIENT_FIELD_ONLY**。

**查表与误差。** 表含 {tab['grid_points']} 个节点，覆盖 h/a=0.001–20，包含原始参考节点和源公式分支附近节点。先冻结网格及 holdout，再生成；每个尺寸有 5,500 个独立 holdout，每点用 128 个归一 generalized vectors 核查阻力 action 与 mobility action。生产 C++ 查表与独立 NumPy 查表一致，合格区间最大 action error 为 {tablegate['max_qualified_action_error']:.10g}，低于 1e-3。全表最大互易误差 {a['max_reciprocity_error']:.10g}；最小 scaled mobility eigenvalue {a['minimum_mobility_scaled_eigenvalue']:.10g}，总阻力最小 scaled eigenvalue {a['minimum_total_resistance_scaled_eigenvalue']:.10g}。矩阵未经平滑、拟合或对称化；eigenvalue 审核中的对称部分只用于诊断。

R_SI = D⁻¹ R_scaled D⁻¹，其中 D=diag(sqrt(M_T,bulk),sqrt(M_R,bulk))。HDF5 内完整保存 6×6 TT、RR、TR、RT，运行时缓存全表，用 log(h/a) 分段线性插值。总装配为 R_bulk + R_pair,excess + R_wall,excess；bulk 只计算一次。wall RHS 加 R_wall,excess q_inf，q_inf=[U_inf, 0.5 curl U_inf]，沿用已验证梯度路径。

135 组平壁 normal/tangent/torque 响应与 RMBW 直接矩阵一致，最大非零分量相对误差 {normalgate['maximum_component_error']:.8g}；10,000 随机法向的独立旋转误差 {a['rotation_error']:.8g}。C/D/E 的通过代表对冻结 RMBW 参考的实现正确，不消除参考自身 RR/TR 的精度限制。单个较小 TR 系数的相对误差可以大于 action 误差（最坏 holdout 约 0.547%），保留在完整 CSV 中，没有隐藏。

**参考与闭合限制。** h/a≤0.2 为近场高置信区，≤5 为 qualified 区，5–20 为 limited 区。这些标签不意味着 finite-gap RR/TR 已经有独立精确解背书。保留上轮的 RR_parallel 常数 0.3817 与文献 0.3709 差异，以及 h/a>9.018296 上游 TR=0、RR 回到 bulk 的分支限制。本模型只处理相对局部背景的阻力；**任意背景应变引起的额外 wall stresslet forcing 未实现**。因此单泡若 V=U_inf 且 Ω=Ω_inf，壁面 excess 项可为零，中心轨迹未必改变。

h/a<0.001 时只钳制查表阻力，不视为防穿墙机制。h/a>20 将 excess 置零，明确为 DEVELOPMENT BULK HANDOFF：总 scaled normal TT 跳变 {handoff['normal_TT_jump']:.8g}，parallel TT 跳变 {handoff['parallel_TT_jump']:.8g}，最大矩阵 action 跳变上界 {handoff['max_matrix_action_jump_upper_bound']:.8g}。实际外移短测在相邻采样点观察到 normal velocity 相对跳变 {handoff['normal_velocity_relative_jump']:.8g}，parallel 跳变 {handoff['parallel_velocity_relative_jump']:.8g}；这含很小的相邻 gap 变化，精确 h/a=20 矩阵跳变另列 [handoff 审核](validation/FAR_FIELD_HANDOFF_AUDIT.json)。没有把有限距离的壁面作用宣称为真实物理零。

**几何、局部平面及真实场阻断。** C++ 保留 BVH 最近三角形查询；n=(X−P)/|X−P|，三角形原始法向、朝球调整法向与 PCA 法向分别留存。切向基取与 n 最不平行的全局轴，经 Gram–Schmidt 构成右手系；相邻步按前一步切向符号保持确定性。

冻结 PCA patch 定义是：最近点周围半径 2a 球相交的完整三角形，顶点权重 area/3，按原始三角形 ID 排序。完整三角形的顶点可以延伸到 patch 球外，结果依赖此明确的离散几何解释。RMS/a≤0.10、面积加权法向角 P95≤20°，最近面法向方向歧义门为 20°。未用 RMS 推导未经校准的曲率。

真实几何 10,000 个安全内部点的同 STL C++/VTK 最近距离误差至多 {geom['max_distance_error_m']:.8g} m；PCA 分类无分歧。仅 4 个点通过，均为 d10、位于几何平坦端盖；其余 9,996 点未通过。端盖也在原闭合 STL 中，不能把这些点解释成已验证的生理血管近壁区域。该比例是本次固定抽样和三个尺寸下的结果，不是连续血管体积的精确有效比例。

Case I 先搜索全部 182,694 个 frozen field 节点的 h/a=0.1 投影，再对 67,262 个三角形中心、两个法向符号及 h/a=0.05/0.2 检查额外 269,048 个静态候选。没有找到同时满足 d50、局部平面、流体 8 角插值 stencil 的起点。**这不是对连续域不存在可行点的数学证明**，但没有已验证起点可用于本轮时域测试。Case J 原起点的 RMS/a 约 0.237–0.352，P95 约 40.18–66.70°；在第一步之前停止，t=0，目标仍为 0.02 s。没有更换尺寸、移动原 8 泡起点、减小时长或放宽门限。

{case_table}

Case J 原始 CSV 的零速度列是“未求解初态”的占位，不是物理解；`wall_active` 是 gap 区间候选标记，不表示已成功应用壁面模型。其 `min_pair_gap_m=1e100` 是未推进的 sentinel，不能用于科学摘要。独立重算的初始最小壁面间隙为 {jgap*1e6:.10g} µm，最小两泡间隙为 {jpair*1e6:.10g} µm。解释见 [初态语义](cases/J_real_8/INITIAL_STATE_SEMANTICS.json)。不存在可比较的新 0.02 s 真实轨迹。

**硬约束与时间步。** 在原 active-set 框架中按粒子 ID/三角形 ID 排 wall rows，再排 sphere-pair rows；约束为 h_base+dt n_base·V≥0。RK2 midpoint 和完整 proposal 均做连续 swept sphere/triangle 验证，拒绝后减半子步；BVH 查询结合距离的 1-Lipschitz 下界递归证明安全，达到递归深度上限时拒绝，不接受未证明区间。

200 个独立扫掠测试中，100 个安全、100 个拒绝，包含 50 个端点均安全而中间穿越的案例。两墙约束加一个两泡约束同时激活的静态 12×12 KKT 检查通过，速度误差 {supp['joint_velocity_error_m_s']:.8g} m/s，角速度误差 {supp['joint_omega_error_rad_s']:.8g} rad/s。时域两泡与壁面共同装配短测另保留 20 个 accepted steps，pair 最小间隙约 1.435 nm。

三个 C_wall 候选均运行到 0.003 s，对齐到三组时间节点的并集比较；C_adv=0.25、pair C_gap=0.4 不变。冻结位置/间隙门 1e-8 m、峰值归一角速度差 2%、积分转角差 1e-3 rad。

{cwtable}

因此选择 **C_wall=0.2**。0.4 的角速度差略超门，没有四舍五入为 PASS。初始 case contract 与候选配置保持原文件；后续 `.selected.cfg` 由冻结选择规则只替换 C_wall，另存 SHA 和选择记录，独立核对其余参数逐字不变。

硬壁组 accepted 最小间隙为 0，穿墙计数 0；无壁面对照最小间隙 −58.102 nm，仅阻力对照 −8.418 nm，都是刻意关闭硬约束的负对照。汇总安全计数不把它们混入硬约束通过组，完整负对照记录仍保留。根目录 WALL_INTERACTION_EVENTS.csv 带 case 和 control 标签，并排除未求解的 Case J 初态占位。

**MPI、Kokkos 与独立 finalizer。** 128 粒子真实 LAMMPS MPI1/MPI4 测试到 0.002 s；MPI4 发生 {mpi['MPI4_owner_changes']} 次实际 owner 变化，记录 {mpi['wall_active_observations']} 次 accepted wall-active 观测。位置/速度/角速度差分别为 {mpi['position_error_m']}/{mpi['velocity_error_m_s']}/{mpi['omega_error_rad_s']}。Kokkos 使用现成 22Jul2025 Update6 GPU build，仅做 8 粒子、2 步兼容性测试；没有做性能优化，GPU_PERFORMANCE_READY=NO。

独立 Python finalizer 不加载生产 C++，重算查表误差、正定/互易性、全部静态几何/PCA、坐标旋转、RK2 两阶段解、连续最小间隙、MPI 身份，以及真实案例阻断是否正确。大 MPI 案例先确认无近场 pair，再独立求解各 6×6 block；小耦合系统使用完整 NumPy KKT active-face enumeration。`INDEPENDENT_FINALIZER=PASS` 指证据与门执行核查通过，不能替代 Case I 的成功时域结果。

**图像与 ParaView。** 已生成 12 张 PNG 和 4 个规定 VTP，另加一个按查询覆盖标注的真实表面 VTP。Case I 的轨迹 VTP 为空且注明 BLOCKED；Case J 只含 8 个初始点，无伪造轨迹线。它们的速度/角速度用 NaN 表示不适用，不属于 solver nonfinite。真实图 9–12 明确展示阻断原因或原结果与新初态，不能冒充完整轨迹图。静态 surface 上灰色表示未查询，validity 随粒径变化，不外推未查询面。`HUMAN_VISUAL_REVIEW=PENDING`，生成和机器读取不替代人工审核。

**完整性。** 远端输入 SHA 验证、编译、运行 receipt、终态与逐文件结果清单均归档；WSL 已实际 `sha256sum -c REMOTE_SHA256SUMS` 通过 {integrity['files_verified']} 文件。远端原 LAMMPS 源码完整检查 15,195 文件无变化；本地 7 个既有基线清单和 RMBW 固定 commit/已跟踪源码前后核查通过。Palabos、LAMMPS 核心、passive transport、bubble-pair baseline、SonoVue、真实场和几何未修改。新增驱动、环境和输入副本都属于本轮独立目录。

保留了两项运行前工程错误：source archive 没有 .git，改用原有全源码 SHA 清单；合成零场 integer dtype 导致加 affine float 梯度失败，在创建 case contract 和任何运行之前修正为 float。编译和正式远端 case 均一次成功启动/退出，没有隐藏失败运行或重试。

**复核入口。** 阅读 [最终摘要](FINAL_TERMINAL_SUMMARY.txt)、[冻结 stage contract](contracts/MICROBUBBLE_WALL_HYDRODYNAMICS_V0_CONTRACT.json)、[表 contract](contracts/RMBW_WALL_TABLE_CONTRACT.json)、[独立 finalizer](src/finalize_wall_hydrodynamics_v0.py)、[图像来源](VISUALIZATION_PROVENANCE.json) 和 [复现说明](README.md)。SHA256SUMS 覆盖最终正式产物；远端下载清单单独保留。

**仍待用户审阅。** 本轮说明相同 RMBW 查表与约束实现可在技术场景中通过核查，但当前真实几何/尺寸与冻结局部平面模型之间存在严重适用性冲突。下一步需要人工审查几何有效性图和局部平面模型适用范围，再决定后续模型工作。不能将这个阻断转换成完整 wall hydrodynamics PASS，也不能自动开启 adhesion。AMBIENT_STRAIN_INDUCED_WALL_FORCE、FULL_EXACT_WALL_HYDRODYNAMICS、MICROBUBBLE_ADHESION、BUBBLE_BUBBLE_TWIST 均 PENDING；buoyancy/lift/RBC 均 OFF。
'''
(R/'MICROBUBBLE_WALL_HYDRODYNAMICS_V0_REPORT.md').write_text(report)
readme=f'''本目录保存 MICROBUBBLE_WALL_HYDRODYNAMICS_V0 的代码、冻结输入、实际结果和独立核查。

总状态为 BLOCKED_LOCAL_PLANE_VALIDITY；技术实现通过，Case I/J 未获得完整真实时域结果。先读 [中文报告](MICROBUBBLE_WALL_HYDRODYNAMICS_V0_REPORT.md)。

复核只读取结果，不执行求解器：

```bash
cd {R}
sha256sum -c SHA256SUMS
OPENBLAS_NUM_THREADS=1 {W}/env/bin/python -B src/finalize_wall_hydrodynamics_v0.py --phase offline
OPENBLAS_NUM_THREADS=1 {W}/env/bin/python -B src/finalize_wall_hydrodynamics_v0.py --phase runtime
OPENBLAS_NUM_THREADS=1 {W}/env/bin/python -B src/finalize_wall_hydrodynamics_v0.py --phase supplemental
```

独立 finalizer 会重写 validation 结果文件；若要保留原始归档，请先复制目录再复核。环境版本见 provenance/LOCAL_VALIDATION_REQUIREMENTS.txt。离线 RMBW 生成使用之前冻结的 Python3.8 参考环境；不是这个新 Python3.12 核查环境。生成脚本有 FIRST 输出保护，不覆盖原始 RMBW 输出。

生产驱动位于 src/rigid_lammps_main.cpp；仅链接既有 LAMMPS 静态库，不重新编译其核心。远端 cmake 配置/编译命令见 provenance/control_scripts/build_remote.py，源与二进制 SHA 见运行 receipt 和 provenance/RUNTIME_SOURCE_SHA256.json。全部现有 case 已到终态；run_wall_phase.py 拒绝隐式重跑。不要直接重启现有 case 来冒充新证据。

静态表与 contract 在 tables/、contracts/；几何/PCA10k×2 在 raw/、validation/；全部实际 timestep 和两阶段记录在 cases/；12 张图和 ParaView 数据在 visualization/。Case I VTP 为空；Case J VTP 只有初始点，NaN 表示未求解量。

远端原始数值结果：{P['remote_result']}
远端独立构建：{P['remote_work']}
本地独立项目：{W}

RMBW 源码/GPL 来源和 Eigen 头文件许可分别在 provenance/reference_audit/、provenance/EIGEN_COPYRIGHT.txt 中。生产未复制 RMBW 实现源码。
''';(R/'README.md').write_text(readme)
for sub in ['src','scripts','contracts','tables','validation','raw']:
 p=W/sub
 if not p.exists():p.symlink_to(R/sub,target_is_directory=True)
(W/'README.md').write_text('本轮项目代码、表、输入与核查归档位于：\n'+str(R)+'\n\n'+readme)
print('REPORT_AND_SUMMARY_WRITTEN',S['MICROBUBBLE_WALL_HYDRODYNAMICS_V0'])
