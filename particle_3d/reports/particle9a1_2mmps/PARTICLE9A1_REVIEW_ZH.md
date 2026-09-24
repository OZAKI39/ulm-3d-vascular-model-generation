# 一句话结论

自动状态 **P9A1_AUTOMATED_VALIDATION_PASS**；人工状态 **PENDING_USER_REVIEW**。P9-A.1 same-12 完成 11/12；30 条：EXECUTED；500 条：EXECUTED。

# 这轮修了哪三个问题

1. **入口附近错误接管**：保留真实 FEM 背景载荷，局部平面剪切只提供经过 canonical benchmark 标定的仿射修正；真实开放边界 rim 使用 P6.5 fallback。
2. **2 nm 约束反复开关**：在当前已解出的运动路径上寻找第一次到达 handoff 下限的时刻，先推进到事件，再重算剩余真实时间。
3. **重复接触条件**：用数值秩和非负锥冗余检查移除重复方向，再检查全部原始约束，保留真实不同的不可穿透信息。

# 为什么没有改 mobility 系数

原单平面 benchmark 没有显示 mobility 系数错误。M_eff、R_eff、delta_R、平移/旋转耦合常数、互易投影、法向 P6.5 阻力和 2 nm 下限均保留。源码逐项 AST/字节比较见 `data/unchanged_contracts.json`；原连续证书函数、P6.5 装配、admission、population source 与所有原积分守卫实现均未改写。
本地 84 项测试、服务器 24 项、共享核心 96 项测试通过。新旧 P9-A canonical 最大差异：Vt = 2.16840434e-19 m/s，Omega = 2.84217094e-14 s^-1。
相对原始 2D 四位小数非完全互易 reference：Vt 最大差异 3.1769044e-09 m/s，Omega 0.00996369103 s^-1；均在原有解析互易投影误差界内，没有扩大测试容差。

# 新的 bulk-preserving 方法是什么意思

“墙可以改变球怎么滑和转，但不能因为局部剪切梯度很小，就假装原本存在的血流消失了。”

设 $q_{lin}=[T^T(Hs),\;aT^T(n\times s)/2]^T$，$q_{ref}=(1-w)q_{lin}+wq_{wall}$。使用

$$\Delta b=R_{eff}q_{ref}-R_{free}q_{lin},\qquad R_{eff}q=R_{free}q_{bulk,actual}+\Delta b.$$

等价地，$q=q_{ref}+M_{eff}R_{free}(q_{bulk,actual}-q_{lin})$。canonical 流中真实 bulk 等于 q_lin，恢复旧解；远壁 w=0 恢复 free/P6.5。真实 bulk 的差额始终进入载荷。法向自由度没有加入此块。
ID 22 平均路径速度从 1.6173317e-08 m/s 变为 0.0018697796 m/s。这里平均值是实际路径长度/已积分物理时间。

# inlet rim 怎么处理

开放边界和血管壁的交线不是已验证的无限平面墙。按原始 GlobalNodeID，求 WALL perimeter edges 与各 cap perimeter edges 的严格交集。INLET/OUTLET_01/OUTLET_02/OUTLET_03 分别有 29/22/19/23 条 rim edges。全部 cap 边缘必须匹配，文件 SHA 随拓扑记录保存。
用原最近距离 kernel 的 EDGE/VERTEX barycentric witness 对照该拓扑；只有真实 rim witness 才令 P9 切向 delta_R 和仿射 delta_b 为零。普通 WALL shared edge 保持 P9 修正。未增加坐标距离带，P6.5 normal lubrication/handoff 保持工作。

| ID | FEM 向内速度 (m/s) | 旧 P9-A | 新 P9-A.1 | 新轨迹结果 |
|---|---:|---:|---:|---|
| 3 | 0.00240193063 | -1.39336009e-09 | 0.00240193062 | COMPLETED |
| 15 | 0.00240193062 | -1.02544293e-08 | 0.00240193062 | OTHER |

# handoff 怎么修

预测当前路径会到达 2 nm 边界时先处理首次到达事件，而不是越界后无限缩小时间步。有限三角形距离的凸性用于寻找步内最小距离，因此两端都清楚但中间穿过的情况也会被检查。首次到达用原有限几何求根，真实时间推进至事件后重新装配和求解剩余时间；没有位置投影、gap clamp、下游偏移或虚假时间。radius-dependent 原下限继续适用。
另一个数值问题是安全的边缘掠过路径：真实最小间隙高于下限，但原端点支撑平面证明过于保守。对此仅在最小距离处分割证明区间，每段重新调用完全未改动的 wall_handoff_certificate。只有全部原证书通过才接受同一条路径，运动路径、速度和位置都不改变。原 swept_clearance_certificate 仍检查完整区间。真实穿越测试必须拒绝。
根搜索上限 64 次来自 binary64 分辨率；事件及证明分段使用原 MAX_REFINEMENT_DEPTH=48 实现守卫。原 provider/progress/stationary/horizon 守卫均未放宽。新增失败会保留，不会为了出血管而调参数。
开发证据也保留：development_trial_01 是仅端点括根版本，development_trial_02 增加步内最小值检测。最终版本进一步组合原连续证书；两轮开发记录不是正式结果。

# ID 7 怎么修

“两条几乎重复的接触命令只保留独立信息。”在白化约束 W 上做 SVD，因 dual 使用 W W^T，数值可分辨阈值取 sqrt(max(W.shape)*eps)*norm(W,2)，不使用固定角度或距离。还要求被移除行能由保留行非负组合表示，避免错误删除相反方向的接触。解后用保留系统的严格预算检查所有原始行；若失败即拒绝，不能利用原大条件数扩大接受预算。
实际 ID 7 原矩阵 fixture：约束 2→1，旧 solver numerical rank 2→1；同一 Gram 分辨率下 rank before/after 均为 1，说明只删冗余信息。condition 1.1949323e+13→1。最小法向速度 −6.76861815e−7→-1.08420217e-19 m/s。保留/删除 ID、奇异值、残差、非负系数全部见 JSON。

# same-12 结果

严格复用原 500 admission ledger 中前 12 个 accepted ID 的 birth time、radius、position、quaternion、SonoVue draw、seed 和 admission。新阶段 ledger 只增加阶段 provenance，事件内容逐项相同。dt=0.00025 s，horizon=1.5 s，2.0 mm/s FEM SHA=129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d，Q=1.5513591604402322e−14 m³/s，原出口 classifier 不变。

| ID | P6.5 原版本 | P9-A 原版本 | P9-A.1 | P6.5 + 新 handoff |
|---|---|---|---|---|
| 3 | COMPLETED | INLET_ESCAPE | COMPLETED | COMPLETED |
| 5 | COMPLETED | COMPLETED | COMPLETED | COMPLETED |
| 7 | HANDOFF_OR_SAFETY_STOP | HANDOFF_OR_SAFETY_STOP | COMPLETED | COMPLETED |
| 10 | COMPLETED | COMPLETED | COMPLETED | COMPLETED |
| 11 | HANDOFF_OR_SAFETY_STOP | COMPLETED | COMPLETED | COMPLETED |
| 12 | HANDOFF_OR_SAFETY_STOP | HANDOFF_OR_SAFETY_STOP | COMPLETED | COMPLETED |
| 14 | COMPLETED | HANDOFF_OR_SAFETY_STOP | COMPLETED | COMPLETED |
| 15 | HANDOFF_OR_SAFETY_STOP | INLET_ESCAPE | OTHER | OTHER |
| 17 | HANDOFF_OR_SAFETY_STOP | HANDOFF_OR_SAFETY_STOP | COMPLETED | COMPLETED |
| 18 | HANDOFF_OR_SAFETY_STOP | COMPLETED | COMPLETED | COMPLETED |
| 19 | COMPLETED | COMPLETED | COMPLETED | COMPLETED |
| 22 | HANDOFF_OR_SAFETY_STOP | TIME_LIMIT | COMPLETED | COMPLETED |

完成数量：P6.5 5/12；旧 P9-A 5/12；新 P9-A.1 11/12；P6.5 + 新 handoff 11/12。

| ID | 模型 | 失败原因 | handoff events | rejected trials |
|---|---|---|---:|---:|
| 15 | P9A1 | EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS_NOT_PHYSIOLOGICAL_TRAPPING_PROOF | 27 | 30 |
| 15 | P65_NEW | EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS_NOT_PHYSIOLOGICAL_TRAPPING_PROOF | 22 | 77 |

ID 7 实际新轨迹中的最大 contact condition 为 14.774131038863084；ID 12/14/17 新轨迹均完成，原 HANDOFF_ENDPOINT_BELOW_LOWER 拒绝次数均为 0。仍有原几何 SWEEP_NOT_CERTIFIED 导致的真实时间细分，不能解释成所有细分均已消除。

# 30-track smoke

EXECUTED。outcomes: `{'COMPLETED': 29, 'OTHER': 1}`。仅 same-12 gate 通过才会运行；选择原 ledger 前 30 个事件，已经计算的 12 条只在 source/event identity 一致时复用。
30 条中的唯一未完成是 same-12 已出现的 ID 15：三个独立 FACE 接触，rank=3，所有乘子为正，contact condition≈7.16，平移速度在原 KKT 数值预算内为零；P6.5 同位置也是三个独立接触。新增的 18 条全部完成。因此按“无新的系统性数值错误”的门槛通过，不把这条静止路径改标为 completed，也不以 30/30 为门槛。证据见 data/stationary_stop_audit.json。

# 500-track formal

**EXECUTED**。启动时门槛记录独立保存在 `data/production_launch_gate.json`，之后的完整验证另计。
正式批次复用本轮同版本 smoke 的 30 条记录，另外积分 470 条；原 P9-A 旧轨迹没有作为本阶段轨迹缓存使用。正式记录 500 条，完成 486 条；outlet distribution：`{'OUTLET_01': 0, 'OUTLET_02': 477, 'OUTLET_03': 9, 'NO_EXIT': 14}`；failure counts：`{'EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS_NOT_PHYSIOLOGICAL_TRAPPING_PROOF': 14}`。所有未完成 ID：`[15, 87, 100, 202, 221, 330, 378, 389, 440, 475, 504, 523, 597, 705]`。

最终状态重装配支持三个独立接触下静止的未完成数：14；尚未解释的失败 ID：`[]`。这只审核保存的最后状态，没有重新积分、改变粒径或强迫通过。详见 `data/production_outcomes.csv`、`data/production_failed_final_states.json`。

# 当前限制

仍是 single nearest planar wall、no curvature、no multi-wall hydrodynamics、no shear lift、tetra-local gradient、independent MB。没有新增 RBC 闭合、MB 变形、CFD 或经验平滑。open-rim 单点拓扑 gate 没有额外 smoothing band；在网格特征切换处可能不连续。硬接触可能同时检查多个真实方向，这不代表加入 multi-wall hydrodynamics。
原守卫停止或静止不等于证明生理捕获；未完成的轨迹不能作为完整通过数据使用。30/500 gate 不通过时不得把 same-12 或开发轮次称作正式 500 数据。

自动门槛：

- unit_regressions: PASS
- canonical: PASS
- p65_normal: PASS
- id3_no_inlet_escape: PASS
- id15_no_inlet_escape: PASS
- id22_no_old_collapse: PASS
- id7_contact_rank: PASS
- old_handoff_chatter_removed: PASS
- no_wall_penetration: PASS
- no_new_solver_pathology: PASS
- smoke30_no_systematic_failure: PASS
- evidence_and_time_identity: PASS
- formal_all_500_recorded: PASS
- formal_no_penetration: PASS
- formal_no_unresolved_solver_failure: PASS

无穿透核对：0 条 WALL penetration；0 条 handoff 下限违规。时间/位置/输入/source/hash 独立校验：True。离线验证读取保存的逐步位置与 gap 并核对哈希；连续路径安全由在线原证书逐段保证，未声称离线重新计算全部几何查询。
服务器 `/workspace/particle9a1_2mmps_20260924T001605Z`；6 CPU workers，OMP/OPENBLAS/MKL/NUMEXPR 均为 1 线程；GPU 未用于积分。最终各批耗时（秒）：`{'same12': 116.76379561424255, 'smoke30': 47.7155487537384, 'production': 1155.6788318157196}`。开发轮次另有日志和逐轨迹 receipt，不能计作正式结果。
保护核对：341 个原报告/数据/冻结输入/用户已有修改文件的 SHA 已比较；mismatches=0；原 git index 未改：True。分支 dev/particle9a1-production-compatibility，未 commit/push。

修改/新增科学源码：contact_redundancy.py, nearfield_handoff.py, open_boundary_rim.py, particle65_motion.py, particle9a1_audit.py, particle9a_motion.py, physical_time_refinement.py, planar_wall_hydrodynamics.py, resistance_solver.py, wall_geometry.py。`logs/git_diff_stat.txt` 是完整工作区 diff，包含本轮开始前已有修改；本轮原始状态见 `logs/git_before.json`。

图与可直接重绘的来源：

- [01_same12_before_after.png](figures/01_same12_before_after.png) / [PDF](figures/01_same12_before_after.pdf)；来源：data/same12.csv, data/same12.json
- [02_inlet_3_15.png](figures/02_inlet_3_15.png) / [PDF](figures/02_inlet_3_15.pdf)；来源：data/inlet_3_15.csv, data/inlet_3_15.json
- [03_canonical_planar.png](figures/03_canonical_planar.png) / [PDF](figures/03_canonical_planar.pdf)；来源：data/canonical.csv, data/canonical_summary.json
- [04_handoff_before_after.png](figures/04_handoff_before_after.png) / [PDF](figures/04_handoff_before_after.pdf)；来源：data/handoff.csv, data/same12.csv
- [05_contact_redundancy_id7.png](figures/05_contact_redundancy_id7.png) / [PDF](figures/05_contact_redundancy_id7.pdf)；来源：data/contact_id7.csv, data/id7_fixed_fixture.json, reference/id7_duplicate_contact.json
- [06_real_same12_trajectories.png](figures/06_real_same12_trajectories.png) / [PDF](figures/06_real_same12_trajectories.pdf)；来源：data/real_paths.csv, data/trajectory_geometry.json, data/same12.json

入口：[OPEN_RESULTS.html](OPEN_RESULTS.html)；机器记录：[gates.json](data/gates.json)；逐轨迹安全/时间验证：[independent_validation.csv](data/independent_validation.csv)。

# 人工审核

- [ ] inlet 修复合理
- [ ] canonical planar physics 未破坏
- [ ] handoff 不再 chatter
- [ ] duplicate contact 处理合理
- [ ] same-12 轨迹合理
- [ ] 30-track smoke 合理
- [ ] 若已运行，500-track 数据合理
- [ ] 同意进入下一阶段
