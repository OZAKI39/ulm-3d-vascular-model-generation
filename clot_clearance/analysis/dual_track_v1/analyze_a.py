from analysis.dual_track_v1.common import *
import csv

def main():
 topo=json.loads((A/'PERIODIC_TOPOLOGY_AUDIT.json').read_text());g=json.loads((A/'PVCE_ENERGY_RELEASE_RESULTS.json').read_text());low=[json.loads((A/f'A{i}/EQUILIBRIUM.json').read_text()) for i in range(3)];rec=json.loads((A/'STATIC_SOLVER_RECOVERY.json').read_text());high=json.loads((A/'lps/wide_A0_search_1/RESULT.json').read_text());fd=json.loads((A/'READ_ONLY_FD_SCALE_DIAGNOSIS.json').read_text());oldtwo=6.1859550723487025e-6;oldthree=6.170305164660184e-6
 summary=dict(geometry=topo,low=low,G=g,old_G=dict(two_J_m2=oldtwo,three_J_m2=oldthree),relative_G_change=dict(two=g['G_two_J_m2']/oldtwo-1,three=g['G_three_J_m2']/oldthree-1),recovery=rec,high=high,read_only_FD_diagnosis=fd,formal_equilibrium_attempts=5,measured_low_high_wall_s=sum(v['wall_seconds'] for v in low)+high['wall_seconds'],recovery_wall_upper_bound_s=rec['wall_upper_bound_s'],measured_low_high_CPU_s=sum(v['cpu_seconds'] for v in low)+high['cpu_seconds'],peak_observed_RSS_MiB=max(topo['peak_RSS_MiB'],high['peak_RSS_MiB'],*[v['peak_RSS_MiB'] for v in low]),recovery_CPU_RSS='Unavailable due reporting-only exception; do not replace by zero',delta_over_dimension=dict(old_W=.1505/.6,wide_W=.1505/5,H=.1505/1,old_B=.1505/.6,wide_periodic_B=.1505/.4),scientific_stop='Low3 accepted, old saved high recovery certified; one wide high candidate rejected on frozen coordinate plateau and1pm component FD. No further formal A solves. Read-only1nm FD diagnoses scale/roundoff without regrading acceptance.')
 save(A/'ANALYSIS_SUMMARY.json',summary)
 status=dict(WIDE_GEOMETRY_VALIDITY='PASS_16000_POINTS',PERIODIC_BOND_TOPOLOGY='PASS_SINGLE_RANK_STATIC_Y_PERIODIC',LPS_FORCE_ENERGY_CONSISTENCY='LOW_PASS_NATIVE_AND_FD_HIGH_COMPONENT_FD_GATE_NOT_PASSED',LOW_LOAD_PVCE_EQUILIBRIUM='PASS_3_OF_3',LOW_LOAD_G_DRIVE='PASS_SPECIMEN_DRIVING_FORCE',STATIC_SOLVER_RECOVERY='PASS_SAVED_OLD_STATE_CERTIFIED_METADATA_LIMIT_REPORTED',NEAR_ONSET_LOCATION='UNACCEPTED_HIGH_CANDIDATE_GRIP_DOMINANT_NO_ONSET_IDENTIFIED',NEAR_ONSET_G_DRIVE='NOT_RUN_ENTRANCE_GATE_FAILED',MATERIAL_GC='NOT_IDENTIFIED',NUMERICAL_SCOPE='FINITE_XZ_PERIODIC_Y_SMALL_PURE_SHEAR_LIKE_STATIC_DIAGNOSTIC_NO_BIOMATERIAL_VALIDATION')
 save(A/'FINAL_STATUS.json',dict(run_id=RUN,states=status,formal_equilibria=5,remaining_budget_unused=3,stopped=True))
 rows=[]
 for i,v in enumerate(low):rows.append(dict(model='wide_periodic_Y',case=f'A{i}',opening_m=2e-6,crack_length_m=.0015+i*.00005,U_J=v['U_J'],G_two_J_m2=g['G_two_J_m2'] if i==0 else '',G_three_J_m2=g['G_three_J_m2'] if i==0 else '',accepted=v['accepted'],maximum_free_force_N=v['max_native_free_force_N'],wall_s=v['wall_seconds']))
 rows.append(dict(model='old_finite_small_V1',case='A0_A1_A2_reused',opening_m=2e-6,crack_length_m=.0002,U_J='',G_two_J_m2=oldtwo,G_three_J_m2=oldthree,accepted=True,maximum_free_force_N='',wall_s='NOT_RERUN'))
 rows.append(dict(model='wide_periodic_Y',case='high_A0_candidate',opening_m=high['opening_m'],crack_length_m=.0015,U_J=high['U_J'],G_two_J_m2='',G_three_J_m2='',accepted=False,maximum_free_force_N=high['max_native_free_force_N'],wall_s=high['wall_seconds']))
 with (A/'CASE_COMPARISON.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 text=f'''# LPS_WIDE_PURE_SHEAR_PVCE_FEASIBILITY_V1 科研审查

运行：`{RUN}`。低载三态通过；旧小试样高载静态恢复通过；宽薄高载候选未通过全部冻结门槛，已停止。本报告不以总 PASS 代替分项科学状态。

## 用户要求的 12 项回答

1. **实际几何、节点与键。** 名义 X×Y×Z＝5×0.4×1 mm，dx＝50 µm，δ＝150.5 µm；100×8×20＝16,000 点，909,392 条完整参考键。体积 2 mm³、质量 2 mg。Y 周期厚度固定为 0.4 mm，宏观 Y 应变为零，X 自由，Z 上下端最外层各 800 点仅控制轴向位移；1,600 个受控自由度、46,400 个自由自由度。该几何是缩小 Pure-Shear 类周期材料层，不是 Liu 2024 的完整实验或真实有限厚度血栓。

2. **周期认证。** 共 131,196 条跨 Y 周期键。独立整数格点枚举与最小镜像 KD 邻居键逐一相等，无自键、重复或遗漏；Y 平移一层后 m 和三个切口集合不变。原子 ID、完整参考键 ID、周期影像和长度均保存于 `reference_topology.npz`。原生单 MPI rank ghost 通信已在完整态及 A0、参考态及周期非仿射位移上核对；未宣称多 MPI rank 或任意大变形周期拓扑已认证。

3. **m、θ、力与能量。** 所有状态使用完整未切参考 m；仅在任务目录把既有 XY 周期储能评估器改为 Y 周期。原生力与独立 −∇U 在非仿射预检查中最大相对差 {max(x['force_relative'] for x in topo['checks']):.6e}，θ 最大相对差 {max(x['theta_relative'] for x in topo['checks']):.6e}。四个方向差分均通过 1e−6 门槛。A0 参考态独立接触检测中的极微小数值接触能约 3.9e−37 J，原生接触为零，属于浮点距离舍入，未出现物理接触。没有用原生不完整 bond PE 代替独立储能。

4. **三态平衡。** A0/A1/A2 均通过相对自由力≤1e−6、最大自由力≤1e−13 N、末十次能量跨度≤1e−8、位置变化≤1e−10 m、反力不平衡≤1e−6、原生/独立力≤1e−9、接触及固定拓扑门槛。原生最大自由力依次为 {', '.join(f"{x['max_native_free_force_N']:.6e}" for x in low)} N。加载总开度均为 2 µm。

5. **PVCE 数值。** 从负 X 面切入 a＝1.50/1.55/1.60 mm，z＝0 参考交点按半开区间删除；删键数 16,416/16,976/17,536，严格嵌套，每增量 560 条。BΔa＝2.0e−8 m²。U0/U1/U2＝{', '.join(f"{x:.12e}" for x in g['U_J'])} J。G_two＝{g['G_two_J_m2']:.12e} J/m²，经验误差估计 ±{g['empirical_error_G_two_J_m2']:.6e}；G_three＝{g['G_three_J_m2']:.12e} J/m²，经验误差估计 ±{g['empirical_error_G_three_J_m2']:.6e}。估计合并残余力做功、能量平台和舍入，**不是严格离散误差界、材料误差或置信区间**。

6. **与旧小试样对照。** 旧低载 G_two/G_three＝{oldtwo:.12e}/{oldthree:.12e} J/m²。本轮分别提高 {summary['relative_G_change']['two']*100:.4f}% / {summary['relative_G_change']['three']*100:.4f}%。旧结果只读复用，没有重算。

7. **可否单独归因。** 不可以。δ/W 从 0.250833 降至 0.0301，δ/H 仍为 0.1505；宽度、厚度、裂纹长度和 Y 边界身份同时改变，当前不是分离这些因素的对照试验。两种 G 的差别属于试样相关驱动力差异，不应当作材料本构改变的证据。

8. **旧高载静态恢复。** 原 73.5605577 µm A0 保存态仅进行一次局部 Hessian-vector 残差修正，不改变材料、载荷、参考键或门槛。原最大自由力 1.129107e−12 N，恢复后 {rec['final']['max_force_N']:.6e} N，所有原门槛通过。原最终双精度状态已保存；随后 stdout JSON 序列化错误导致报告写出失败，现已从该保存态独立重建验证报告，没有再次求解。该次 CPU、峰值 RSS、HVP 计数未恢复，明确保留为未知；墙钟上界为 {rec['wall_upper_bound_s']:.6f} s。

9. **宽薄近临界检查及停止。** 满足入口条件后，只计算一个 A0 候选，开度 {high['opening_m']*1e6:.9f} µm。L-BFGS 浮点终止后使用预先记录的同类局部修正；最终自由力通过，但末十次位置变化 {high['last']['position_change_m']:.9e} m，大于 1e−10 m；1 pm 分量差分误差 {high['finite_difference_scaled_error']:.9e}，大于 1e−6。故该候选**未接受**，不再计算第二候选或高载 A1/A2。

10. **首达阈值位置。** 未建立真实首达阈值区间，也未运行自然断键。未接受候选的最大 s/s00＝{high['risk']['max_utilization']:.9f}，位于有限 X 正端与上下夹持面的交界附近；相应参考键与原子 ID 见 `FIRST_FAILURE_AUDIT.json`，只是空间诊断，不能称为可信起裂事件。Y 没有自由侧面。裂尖距负 X 面 1.5 mm，裂尖带与 X 面交集为空；另在 `SAVED_CANDIDATE_SPATIAL_REGIONS.json` 按 δ 宽度单独报告缺口口部与负 X 边界的交界区，及缺口前缘、夹具、其他 X 自由面和内部五区，不据此改变任何验收门槛。只读复核将差分步长设为固定 1 nm 后，相同保存态误差为 {fd['scaled_error']:.6e}，支持 1 pm 检查受舍入影响的解释；**未据此重评接受态，位置门槛失败仍保留**。

11. **能否称材料 Gc。** 不能。当前是已给定裂纹、固定参考 m、理想化周期试样的 `WIDE_PERIODIC_SPECIMEN_PVCE_G_DRIVE`；没有可信目标缺口临界分支、自然新裂纹面积或生物实验认证。

12. **是否支持修改 LPS 本构。** 没有。当前发现指向边界、有限尺寸和静态数值接受问题，不构成新增软化、Yeoh、Mullins 或修改 s00 的证据。

## 成本、审查图与复现

正式静态调用 5/8：低载 3、旧候选恢复 1、宽薄高载候选 1；没有新增动态算例。低载及宽薄高载记录墙钟合计 {summary['measured_low_high_wall_s']:.3f} s，CPU 合计 {summary['measured_low_high_CPU_s']:.3f} s；另加恢复墙钟上界 {rec['wall_upper_bound_s']:.3f} s。周期预检查 {topo['wall_s']:.3f} s，所记录的峰值 RSS 最大 {summary['peak_observed_RSS_MiB']:.2f} MiB。恢复元数据缺项不以零补齐。

`OPEN_RESULTS.html` 包含真实点阵、周期影像、三种切口、收敛/能量/G、恢复与停止证据、高载实际键位置以及三态真实平衡图。所有原始状态、原生输入、完整键表与日志留在本目录；材料点不是医学外壳。

正式运行前方案：`DESIGN_FREEZE.json`；高载数值步骤冻结：`HIGH_NUMERICAL_FREEZE.json`。复现入口和本轮接口修正见共同目录 `REPRODUCE.txt`。原生源码/库与历史结果是否保护，以共同 `PROTECTION_AUDIT.json` 的完整 SHA 复核为准。
'''
 (A/'REVIEW_LPS_WIDE_PURE_SHEAR_PVCE_FEASIBILITY_V1.md').write_text(text)
if __name__=='__main__':main()
