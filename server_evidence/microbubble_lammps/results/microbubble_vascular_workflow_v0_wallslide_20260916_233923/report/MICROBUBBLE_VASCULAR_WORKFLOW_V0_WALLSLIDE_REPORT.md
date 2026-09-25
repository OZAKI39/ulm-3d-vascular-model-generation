# MICROBUBBLE_VASCULAR_WORKFLOW_V0 — GEOMETRIC_WALL_SLIDING_V0

1. **之前为什么贴壁卡死？** 原始速度仍指向墙内，旧版只会拒步和减小 dt；有限半径气泡最终无法前进。
2. **现在怎么处理？** 在 RK2 的 START 和 MIDPOINT 求解后，仅当 raw 预测线段不安全且速度指向墙内时，删除向壁分量。曲面微小越界在先拒步后允许严格有界的几何位置修正。
3. **有没有 wall force？** 没有。矩阵、RHS、PCG、BB 阻力均未改动，也没有 wall hydrodynamics。
4. **沿墙速度保留了吗？** 已验证记录中保留；独立检查的是切向速度向量，而不仅是速率。
5. **离开墙的运动保留了吗？** 保留；outward 和纯切向平墙测试满足相对变化 ≤ 10⁻¹²。Omega raw/used 完全相同。
6. **越过原停滞时刻了吗？** LOW：0.235279081861 → 0.250000000000 s，达到 0.25 s 终点。MEDIUM：0.070947410132 → 0.075815810579 s，随后由 PAIR_SAFETY_FAILED 停止。
7. **有没有新增穿墙或泡间重叠？** 已核验 accepted 状态及连续线段没有超过继承容差。MEDIUM 最小泡间间隙约 −9.99×10⁻¹³ m，属于原 10⁻¹² m 数值容差；不能写成严格几何间隙始终非负。
8. **有没有自然出口？** 未观察到。LONG 的最后落盘 accounting 时间为 0.494122765648 s，2 个 ACTIVE 气泡，随后触发 STOP_WALL_NORMAL_AMBIGUITY；设定的 0.8 s 没有完成。
9. **MPI1/MPI4 一致吗？** LOW/MEDIUM 完整运行一致，LONG 已落盘 CSV 逐字节相同，两者触发同一停止原因。LONG 没有终态文件，不能声称恢复了完整故障终态。

## 本轮结论：BLOCKED

已实现几何约束并解决两个指定的旧 normal-component stall；真实曲壁长期运行尚未合格。LONG 在 step 18134 后的试算触发用户规定的硬停止条件。随后停止所有模拟和模型修改，只进行已有输出的只读核验、报告与镜像归档。未执行排在 LONG 后面的 CAPACITY_STRESS_LIMIT 和 FLUX_WEIGHTED 新版 release；旧版对应证据只作继承来源。

法向判据检测到最近壁点候选在浮点距离容差内并列且点位置不一致。本阶段尚未区分真实几何非唯一与数值近似并列。fatal exception 未写出触发粒子、查询点和候选三角形；不对根因作超出证据的归因，也未改变判据后重跑。

## 求解与几何约束

第一步仍为 `R q_raw = b`，其中 `R = R_bulk + R_BB_excess`，twist pending。第二步使用 WALL_ONLY 最近点 `P`，`n=(X-P)/|X-P|`。raw stage segment 不安全且 `V_raw·n<0` 时，`V_used=V_raw−min(V_raw·n,0)n`；其余保持原值。第三步由 midpoint RK2 使用 V_used 更新位置，Omega_used=Omega_raw。START 的基点是 X_n、步长 dt/2；MIDPOINT 在 X_mid 查询法向，但最终 proposal 的基点仍是 X_n、步长 dt。

**KINEMATIC WALL CONSTRAINT / WORKFLOW_V0 APPROXIMATION，不是 PHYSICAL WALL HYDRODYNAMICS。** 被投影后的 q_used 一般不再满足原来的 Rq=b；独立 dense solve 验证的是保留的 q_raw。未引入 wall force、R_wall、滚动或旋转约束。

raw 求解、Frozen Flow、adaptive injection、SonoVue、lifecycle、ports 继承 BCFLUX 源码。安全距离仍为 gap≥10⁻¹⁰ m；独立几何回读容差 2×10⁻¹⁴ m，pair swept 容差 10⁻¹² m。没有新增 h/a 接触阈值。

旧递归 Lipschitz 检查无法证明恰好贴平墙的切向线段安全；新壁面路径判据计算线段—三角形最小距离，覆盖端点、边和面相交，并用 BVH 剪枝。原有距离 API 和递归检查保留，工作流两阶段采用精确连续距离证书，未降低安全距离。位置投影后重新检查整个 base→mid 和 base→final 线段，并检查两段的 pair swept safety。

## 有界位置投影及曲面验证

优先测试 velocity-only。内侧圆柱/球面接触的直线切向 midpoint 在任意有限 dt 下都有 O(dt²) 越界；证据保存在 `validation/VELOCITY_ONLY_CURVED_EVIDENCE.json`。据此启用 BOUNDED_GEOMETRIC_OFFSET_PROJECTION：先实际拒步减小 dt；只有该 stage 已删除 inward velocity 且 endpoint deficit 很小时，沿 endpoint 的最近壁法向修回原 margin 加浮点舍入保护量。每 stage 最多一次，修正总量上限固定为 2.5×10⁻¹¹ m（原 margin/4，含 roundoff cushion）。较深 proposal 仍拒步减半，不接受深穿透。

`GEOMETRIC_PROJECTION_EVENTS.csv` 保存每个已完成 attempt 的修正，包括 accepted=false 的回滚尝试。致命试算没有完整日志，不将其纳入投影资格声明。每粒子统计只累计 accepted 修正；wall_stall_count 仅统计可归属的完整记录，致命查询没有 particle ID，不能归属到某个粒子；total_time_under_constraint 是两个 stage 的 dt/2 加权采样时间，不代表真实接触驻留物理。

| 曲面 | dt_max (s) | 步数 | 轨迹误差 (m) | 最大修正 (m) | 累计修正 (m) |
|---|---:|---:|---:|---:|---:|
| cylinder | 0.0001 | 213 | 6.4828e-14 | 7.8127e-13 | 1.5573e-10 |
| cylinder | 5e-05 | 411 | 1.6242e-14 | 1.9532e-13 | 7.7997e-11 |
| cylinder | 2.5e-05 | 811 | 4.0648e-15 | 4.8829e-14 | 3.9031e-11 |
| sphere | 0.0001 | 213 | 6.4828e-14 | 7.8127e-13 | 1.5573e-10 |
| sphere | 5e-05 | 411 | 1.6242e-14 | 1.9532e-13 | 7.7997e-11 |
| sphere | 2.5e-05 | 811 | 4.0648e-15 | 4.8829e-14 | 3.9031e-11 |

两个曲面均持续运动 0.01 s、沿壁路径约 1 µm；dt 减半时误差约按四倍下降，最大及累计修正均下降。此合成资格不能替代 LONG 在真实 STL 上的硬停止。LONG 后段 dt 约为 3.90625×10⁻⁷ s，连续安全约束显著增加计算成本。

## 实际重放与独立核验

所有运行采用 TEST_ONLY direct number flux，**NOT EXPERIMENTAL CONCENTRATION**。LOW λ=18.046472132892898/s、MEDIUM λ=97.75172405316988/s 保持上一阶段配置、source draws 和 seeds（仅 stage 路径变更）。LONG λ=6/s，max_time=0.8 s，adaptive injection 始终开启。权威入口 Q=2.7369132390905703×10⁻¹⁵ m³/s，继承已审计来源；不以 position support 面积分替代。

| Case | 已验证步数 | 时间 (s) | admitted / active | hydrodynamic observations | 原因 |
|---|---:|---:|---:|---:|---|
| LOW | 2559 | 0.250000000000 | 3 / 3 | 0 | MAX_PHYSICAL_TIME |
| MEDIUM | 893 | 0.075815810579 | 7 / 7 | 350 | PAIR_SAFETY_FAILED |
| LONG persisted prefix | 18134 | 0.494122765648 | 2 / 2 | 0 | STOP_WALL_NORMAL_AMBIGUITY；无 RUN_STATE |

| 核验范围 | 速度投影次数 | accepted 位置修正次数 | dense solves | 最大 scaled raw q 误差 (m/s) | 最小 wall gap (m) |
|---|---:|---:|---:|---:|---:|
| LOW_release_mpi1 | 245 | 1 | 28 | 1.6263e-19 | 1.0997e-10 |
| MEDIUM_release_mpi1 | 129 | 2 | 290 | 3.1266e-17 | 1e-10 |
| abort_prefix_mpi1 | 48456 | 5 | 240 | 2.1684e-19 | 1e-10 |

MPI1 唯一轨迹计数（不把 MPI4 重复计入）：已记录位置修正尝试 26896 次，accepted 8 次；accepted 最大 2.31933770548e-11 m，累计 8.378849122e-11 m。所有已记录 attempt 最大 2.46957690354e-11 m。

独立 Python 使用 VTK 最近壁点和 NumPy 连续几何证书重算法向、切向分解、unsafe 激活条件、outward/Omega 不变性、wall/pair swept safety、RK2 raw/used/offset 位置关系，并对所有 BB stage 和稀疏非耦合 stage 做独立 resistance assembly + dense solve。注入前缀、FIFO、source/admitted 语义、权威流量积分、LAMMPS 实际粒子计数继续检查。

LONG 原始中止文件保持原样。`validation/abort_prefix_mpi1/` 是明确标注的只读派生核验视图，不是重跑，不是恢复的终态；common complete CSV prefix 与最后落盘 accounting 都到 step 18134。其 accepted 状态可以核验，致命 attempt 的诊断和完整终态不能恢复；生物 near-wall 结束事件不宣称通过。`MPI_CONSISTENCY.json` 证明完整 LOW/MEDIUM 与 LONG 已落盘 CSV 一致。

## 限制与状态

本模型不能预测真实 lubrication slowdown、wall-induced rotation、rolling velocity、wall hydrodynamic force、wall residence physics 或 adhesion probability。唯一用途是在不可穿透约束下推进 Workflow V0。入口位置支持仍是 WORKFLOW_V0_APPROXIMATION；Frozen Flow 仍为 ENGINEERING_TRANSIENT_FIELD_ONLY；production bubble concentration 为 UNSPECIFIED；真实科研生产就绪为 NO。

| 状态 | 结果 |
|---|---|
| ADAPTIVE_INJECTION_STATUS | PASS |
| RAW_RESISTANCE_SOLVER_STATUS | PASS |
| KINEMATIC_WALL_CONSTRAINT_STATUS | PARTIAL |
| TANGENTIAL_PRESERVATION_STATUS | PASS |
| OUTWARD_MOTION_STATUS | PASS |
| CURVED_WALL_SLIDING_STATUS | PARTIAL |
| HARD_WALL_NONPENETRATION_STATUS | PASS |
| PAIR_NONOVERLAP_STATUS | PASS |
| RK2_STATUS | PASS |
| PARTICLE_ACCOUNTING_STATUS | PASS |
| MPI_STATUS | PASS |
| OLD_LOW_STALL_RESOLVED | YES |
| OLD_MEDIUM_STALL_RESOLVED | YES |
| NATURAL_OUTLET_STATUS | NOT_OBSERVED |
| POSITION_PROJECTION_STATUS | USED_AND_QUALIFIED |
| WORKFLOW_TECHNICAL_STATUS | BLOCKED |
| WALL_HYDRODYNAMICS_STATUS | NOT_IMPLEMENTED |
| FLOW_PHYSICS_STATUS | ENGINEERING_TRANSIENT_FIELD_ONLY |
| PRODUCTION_BUBBLE_CONCENTRATION_STATUS | UNSPECIFIED |
| REAL_SCIENTIFIC_PRODUCTION_READY | NO |
| STOP_REASON | STOP_WALL_NORMAL_AMBIGUITY |
| RMBW_USED | NO |
| CF2003_USED | NO |
| FALADE_BRENNER_USED | NO |
| WSS_USED | NO |
| ADHESION_USED | NO |

PASS 的范围严格限于上述可核验记录与合成测试；KINEMATIC / CURVED 为 PARTIAL，整体 BLOCKED。POSITION_PROJECTION USED_AND_QUALIFIED 只表示合成 refinement 和已记录 accepted prefix 的数值资格，不批准失败试算或后续运行。

## 文件与来源

直接来源为只读 BCFLUX stage，其 545 个文件全量哈希冻结；原始 stable rigid 仅作为水动力来源。见 `provenance/BCFLUX_IMMUTABLE_FILES.json`、`SOURCE_DIFF.patch`、`PROJECT_BUILD_PROVENANCE.json`、`INPUT_MANIFEST.tsv`。原始求解/流场/注入/lifecycle/ports 关键源码逐字节不变，新增代码只在本 wallslide stage。

编译及执行身份包含 compiler、MPI、LAMMPS release/commit/static library、主机、日期、所有源码、binary/library SHA256。`validation/ORIGINAL_INPUT_PRESERVATION.json` 记录本地 794 个保护文件与 Git 完整状态；最终远端检查还覆盖 744 个保护输入和 2 个原 CPU build 文件。

已生成要求的 11 张图，并补充 `synthetic_curved_refinement.png`。所有图标记 TEST_ONLY / NOT EXPERIMENTAL CONCENTRATION。

同步证据：`provenance/SYNC_MANIFEST.json` 与 `validation/REMOTE_LOCAL_SHA256_CHECK.json`；报告中的同步完成声明以最终 SHA256 回执为准。

本地：`/home/lzy/projects/compre_output/microbubble_vascular_workflow_v0/wallslide_20260916_233923`
远端 results：`/workspace/microbubble_lammps/results/microbubble_vascular_workflow_v0_wallslide_20260916_233923`
远端 work：`/workspace/microbubble_lammps/work/microbubble_vascular_workflow_v0_wallslide_20260916_233923`

没有 Git add/commit/push/checkout/reset/clean。没有自动进入近壁理论、adhesion、RBC、GPU 或 large-N 阶段。
