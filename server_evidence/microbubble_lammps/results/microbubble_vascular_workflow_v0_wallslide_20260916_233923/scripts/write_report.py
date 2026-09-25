from pathlib import Path
import csv,json,hashlib,datetime,difflib
S=Path(__file__).resolve().parents[1];B=S.parent/'bcflux_20260916_224835';J=lambda p:json.loads(p.read_text());out=S/'report/MICROBUBBLE_VASCULAR_WORKFLOW_V0_WALLSLIDE_REPORT.md'
validation=J(S/'validation/INDEPENDENT_WORKFLOW_VALIDATION.json');assert validation['status']=='PARTIAL';checks={r['case']:r for r in validation['cases']};assert checks['abort_prefix_mpi1']['status']=='PASS_PREFIX_ONLY'
assert J(S/'validation/LOG_LINKAGE_VALIDATION.json')['status']=='PASS'
mpi=J(S/'validation/MPI_CONSISTENCY.json');assert mpi['status']=='PASS';syn=J(S/'validation/SYNTHETIC_WALL_VALIDATION.json');assert syn['status']=='PASS';ev=J(S/'validation/LONG_TRANSPORT_ABORT_EVIDENCE.json')['cases'][0]
low=J(S/'runs/LOW_release_mpi1/RUN_STATE.json');med=J(S/'runs/MEDIUM_release_mpi1/RUN_STATE.json');oldlow=J(B/'runs/LOW_release_mpi1/RUN_STATE.json');oldmed=J(B/'runs/MEDIUM_release_mpi1/RUN_STATE.json')
assert low['time_s']>oldlow['time_s'] and med['time_s']>oldmed['time_s']
status=dict(ADAPTIVE_INJECTION_STATUS='PASS',RAW_RESISTANCE_SOLVER_STATUS='PASS',KINEMATIC_WALL_CONSTRAINT_STATUS='PARTIAL',TANGENTIAL_PRESERVATION_STATUS='PASS',OUTWARD_MOTION_STATUS='PASS',CURVED_WALL_SLIDING_STATUS='PARTIAL',HARD_WALL_NONPENETRATION_STATUS='PASS',PAIR_NONOVERLAP_STATUS='PASS',RK2_STATUS='PASS',PARTICLE_ACCOUNTING_STATUS='PASS',MPI_STATUS='PASS',OLD_LOW_STALL_RESOLVED='YES',OLD_MEDIUM_STALL_RESOLVED='YES',NATURAL_OUTLET_STATUS='NOT_OBSERVED',POSITION_PROJECTION_STATUS='USED_AND_QUALIFIED',WORKFLOW_TECHNICAL_STATUS='BLOCKED',WALL_HYDRODYNAMICS_STATUS='NOT_IMPLEMENTED',FLOW_PHYSICS_STATUS='ENGINEERING_TRANSIENT_FIELD_ONLY',PRODUCTION_BUBBLE_CONCENTRATION_STATUS='UNSPECIFIED',REAL_SCIENTIFIC_PRODUCTION_READY='NO',STOP_REASON='STOP_WALL_NORMAL_AMBIGUITY',RMBW_USED='NO',CF2003_USED='NO',FALADE_BRENNER_USED='NO',WSS_USED='NO',ADHESION_USED='NO')
status.update(validation_scope='PASS flags refer to synthetic tests, complete LOW/MEDIUM, and the independently verified persisted LONG prefix through step 18134. Fatal trial has no complete terminal state or query identity. No global wall-sliding release PASS.',position_projection_scope='Native synthetic refinement, strict cap and every recorded completed attempt; accepted-prefix corrections validated. Fatal in-flight attempt is outside qualification.',MPI_scope='Complete LOW/MEDIUM and byte-identical persisted LONG CSVs plus identical stop; no recovered terminal state',timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),report=str(out),**J(S/'provenance/STAGE_PATHS.json'))
(S/'WORKFLOW_STATUS.json').write_text(json.dumps(status,indent=2,ensure_ascii=False)+'\n')
lines=['# MICROBUBBLE_VASCULAR_WORKFLOW_V0 — GEOMETRIC_WALL_SLIDING_V0','',
'1. **之前为什么贴壁卡死？** 原始速度仍指向墙内，旧版只会拒步和减小 dt；有限半径气泡最终无法前进。',
'2. **现在怎么处理？** 在 RK2 的 START 和 MIDPOINT 求解后，仅当 raw 预测线段不安全且速度指向墙内时，删除向壁分量。曲面微小越界在先拒步后允许严格有界的几何位置修正。',
'3. **有没有 wall force？** 没有。矩阵、RHS、PCG、BB 阻力均未改动，也没有 wall hydrodynamics。',
'4. **沿墙速度保留了吗？** 已验证记录中保留；独立检查的是切向速度向量，而不仅是速率。',
'5. **离开墙的运动保留了吗？** 保留；outward 和纯切向平墙测试满足相对变化 ≤ 10⁻¹²。Omega raw/used 完全相同。',
f"6. **越过原停滞时刻了吗？** LOW：{oldlow['time_s']:.12f} → {low['time_s']:.12f} s，达到 0.25 s 终点。MEDIUM：{oldmed['time_s']:.12f} → {med['time_s']:.12f} s，随后由 PAIR_SAFETY_FAILED 停止。",
'7. **有没有新增穿墙或泡间重叠？** 已核验 accepted 状态及连续线段没有超过继承容差。MEDIUM 最小泡间间隙约 −9.99×10⁻¹³ m，属于原 10⁻¹² m 数值容差；不能写成严格几何间隙始终非负。',
f"8. **有没有自然出口？** 未观察到。LONG 的最后落盘 accounting 时间为 {ev['last_flushed_accounting_time_s']:.12f} s，2 个 ACTIVE 气泡，随后触发 STOP_WALL_NORMAL_AMBIGUITY；设定的 0.8 s 没有完成。",
'9. **MPI1/MPI4 一致吗？** LOW/MEDIUM 完整运行一致，LONG 已落盘 CSV 逐字节相同，两者触发同一停止原因。LONG 没有终态文件，不能声称恢复了完整故障终态。','',
'## 本轮结论：BLOCKED','',
'已实现几何约束并解决两个指定的旧 normal-component stall；真实曲壁长期运行尚未合格。LONG 在 step 18134 后的试算触发用户规定的硬停止条件。随后停止所有模拟和模型修改，只进行已有输出的只读核验、报告与镜像归档。未执行排在 LONG 后面的 CAPACITY_STRESS_LIMIT 和 FLUX_WEIGHTED 新版 release；旧版对应证据只作继承来源。','',
'法向判据检测到最近壁点候选在浮点距离容差内并列且点位置不一致。本阶段尚未区分真实几何非唯一与数值近似并列。fatal exception 未写出触发粒子、查询点和候选三角形；不对根因作超出证据的归因，也未改变判据后重跑。','',
'## 求解与几何约束','',
'第一步仍为 `R q_raw = b`，其中 `R = R_bulk + R_BB_excess`，twist pending。第二步使用 WALL_ONLY 最近点 `P`，`n=(X-P)/|X-P|`。raw stage segment 不安全且 `V_raw·n<0` 时，`V_used=V_raw−min(V_raw·n,0)n`；其余保持原值。第三步由 midpoint RK2 使用 V_used 更新位置，Omega_used=Omega_raw。START 的基点是 X_n、步长 dt/2；MIDPOINT 在 X_mid 查询法向，但最终 proposal 的基点仍是 X_n、步长 dt。','',
'**KINEMATIC WALL CONSTRAINT / WORKFLOW_V0 APPROXIMATION，不是 PHYSICAL WALL HYDRODYNAMICS。** 被投影后的 q_used 一般不再满足原来的 Rq=b；独立 dense solve 验证的是保留的 q_raw。未引入 wall force、R_wall、滚动或旋转约束。','',
'raw 求解、Frozen Flow、adaptive injection、SonoVue、lifecycle、ports 继承 BCFLUX 源码。安全距离仍为 gap≥10⁻¹⁰ m；独立几何回读容差 2×10⁻¹⁴ m，pair swept 容差 10⁻¹² m。没有新增 h/a 接触阈值。','',
'旧递归 Lipschitz 检查无法证明恰好贴平墙的切向线段安全；新壁面路径判据计算线段—三角形最小距离，覆盖端点、边和面相交，并用 BVH 剪枝。原有距离 API 和递归检查保留，工作流两阶段采用精确连续距离证书，未降低安全距离。位置投影后重新检查整个 base→mid 和 base→final 线段，并检查两段的 pair swept safety。','',
'## 有界位置投影及曲面验证','',
'优先测试 velocity-only。内侧圆柱/球面接触的直线切向 midpoint 在任意有限 dt 下都有 O(dt²) 越界；证据保存在 `validation/VELOCITY_ONLY_CURVED_EVIDENCE.json`。据此启用 BOUNDED_GEOMETRIC_OFFSET_PROJECTION：先实际拒步减小 dt；只有该 stage 已删除 inward velocity 且 endpoint deficit 很小时，沿 endpoint 的最近壁法向修回原 margin 加浮点舍入保护量。每 stage 最多一次，修正总量上限固定为 2.5×10⁻¹¹ m（原 margin/4，含 roundoff cushion）。较深 proposal 仍拒步减半，不接受深穿透。','',
'`GEOMETRIC_PROJECTION_EVENTS.csv` 保存每个已完成 attempt 的修正，包括 accepted=false 的回滚尝试。致命试算没有完整日志，不将其纳入投影资格声明。每粒子统计只累计 accepted 修正；wall_stall_count 仅统计可归属的完整记录，致命查询没有 particle ID，不能归属到某个粒子；total_time_under_constraint 是两个 stage 的 dt/2 加权采样时间，不代表真实接触驻留物理。','',
'| 曲面 | dt_max (s) | 步数 | 轨迹误差 (m) | 最大修正 (m) | 累计修正 (m) |','|---|---:|---:|---:|---:|---:|']
for case in syn['curved']:
 for r in case['refinement']:lines.append(f"| {case['shape']} | {r['dt_max']:.3g} | {r['steps']} | {r['trajectory_error_m']:.5g} | {r['max_correction_m']:.5g} | {r['total_correction_m']:.5g} |")
lines+=['','两个曲面均持续运动 0.01 s、沿壁路径约 1 µm；dt 减半时误差约按四倍下降，最大及累计修正均下降。此合成资格不能替代 LONG 在真实 STL 上的硬停止。LONG 后段 dt 约为 3.90625×10⁻⁷ s，连续安全约束显著增加计算成本。','',
'## 实际重放与独立核验','',
'所有运行采用 TEST_ONLY direct number flux，**NOT EXPERIMENTAL CONCENTRATION**。LOW λ=18.046472132892898/s、MEDIUM λ=97.75172405316988/s 保持上一阶段配置、source draws 和 seeds（仅 stage 路径变更）。LONG λ=6/s，max_time=0.8 s，adaptive injection 始终开启。权威入口 Q=2.7369132390905703×10⁻¹⁵ m³/s，继承已审计来源；不以 position support 面积分替代。','',
'| Case | 已验证步数 | 时间 (s) | admitted / active | hydrodynamic observations | 原因 |','|---|---:|---:|---:|---:|---|',
f"| LOW | {low['accepted_steps']} | {low['time_s']:.12f} | {low['admitted']} / {low['active']} | {low['hydrodynamic_pair_observations']} | MAX_PHYSICAL_TIME |",
f"| MEDIUM | {med['accepted_steps']} | {med['time_s']:.12f} | {med['admitted']} / {med['active']} | {med['hydrodynamic_pair_observations']} | PAIR_SAFETY_FAILED |",
f"| LONG persisted prefix | {ev['common_diagnostic_prefix_step']} | {ev['common_diagnostic_prefix_time_s']:.12f} | 2 / 2 | {checks['abort_prefix_mpi1']['hydrodynamic_observations']} | STOP_WALL_NORMAL_AMBIGUITY；无 RUN_STATE |",'',
'| 核验范围 | 速度投影次数 | accepted 位置修正次数 | dense solves | 最大 scaled raw q 误差 (m/s) | 最小 wall gap (m) |','|---|---:|---:|---:|---:|---:|']
for name,r in checks.items():
 w=r['wall_constraints'];lines.append(f"| {name} | {w['active_velocity_projections']} | {w['position_corrections_accepted']} | {r['independent_dense_solves']} | {r['max_dense_reference_scaled_error_m_s']:.5g} | {r['min_wall_gap_m']:.5g} |")
 for stat in w['particles']:stat['case']=name
 with (S/'events'/f'{name}_WALL_PARTICLE_SUMMARY.csv').open('w') as f:
  wr=csv.DictWriter(f,fieldnames=list(w['particles'][0])+['disclaimer']);wr.writeheader();wr.writerows([{**r,'disclaimer':'NOT EXPERIMENTAL CONCENTRATION'} for r in w['particles']])
accepted_dist=[];all_dist=[]
for case in ['LOW','MEDIUM','LONG_TRANSPORT']:
 for r in csv.DictReader((S/'runs'/f'{case}_release_mpi1/GEOMETRIC_PROJECTION_EVENTS.csv').open()):
  if r.get('disclaimer')!='NOT EXPERIMENTAL CONCENTRATION':continue
  all_dist.append(float(r['correction_distance']))
  if r['accepted']=='1':accepted_dist.append(float(r['correction_distance']))
lines+=['',f"MPI1 唯一轨迹计数（不把 MPI4 重复计入）：已记录位置修正尝试 {len(all_dist)} 次，accepted {len(accepted_dist)} 次；accepted 最大 {max(accepted_dist):.12g} m，累计 {sum(accepted_dist):.12g} m。所有已记录 attempt 最大 {max(all_dist):.12g} m。",'',
'独立 Python 使用 VTK 最近壁点和 NumPy 连续几何证书重算法向、切向分解、unsafe 激活条件、outward/Omega 不变性、wall/pair swept safety、RK2 raw/used/offset 位置关系，并对所有 BB stage 和稀疏非耦合 stage 做独立 resistance assembly + dense solve。注入前缀、FIFO、source/admitted 语义、权威流量积分、LAMMPS 实际粒子计数继续检查。','',
'LONG 原始中止文件保持原样。`validation/abort_prefix_mpi1/` 是明确标注的只读派生核验视图，不是重跑，不是恢复的终态；common complete CSV prefix 与最后落盘 accounting 都到 step 18134。其 accepted 状态可以核验，致命 attempt 的诊断和完整终态不能恢复；生物 near-wall 结束事件不宣称通过。`MPI_CONSISTENCY.json` 证明完整 LOW/MEDIUM 与 LONG 已落盘 CSV 一致。','',
'## 限制与状态','',
'本模型不能预测真实 lubrication slowdown、wall-induced rotation、rolling velocity、wall hydrodynamic force、wall residence physics 或 adhesion probability。唯一用途是在不可穿透约束下推进 Workflow V0。入口位置支持仍是 WORKFLOW_V0_APPROXIMATION；Frozen Flow 仍为 ENGINEERING_TRANSIENT_FIELD_ONLY；production bubble concentration 为 UNSPECIFIED；真实科研生产就绪为 NO。','',
'| 状态 | 结果 |','|---|---|']
for k,v in status.items():
 if k.isupper():lines.append(f'| {k} | {v} |')
lines+=['','PASS 的范围严格限于上述可核验记录与合成测试；KINEMATIC / CURVED 为 PARTIAL，整体 BLOCKED。POSITION_PROJECTION USED_AND_QUALIFIED 只表示合成 refinement 和已记录 accepted prefix 的数值资格，不批准失败试算或后续运行。','',
'## 文件与来源','',
'直接来源为只读 BCFLUX stage，其 545 个文件全量哈希冻结；原始 stable rigid 仅作为水动力来源。见 `provenance/BCFLUX_IMMUTABLE_FILES.json`、`SOURCE_DIFF.patch`、`PROJECT_BUILD_PROVENANCE.json`、`INPUT_MANIFEST.tsv`。原始求解/流场/注入/lifecycle/ports 关键源码逐字节不变，新增代码只在本 wallslide stage。','',
'编译及执行身份包含 compiler、MPI、LAMMPS release/commit/static library、主机、日期、所有源码、binary/library SHA256。`validation/ORIGINAL_INPUT_PRESERVATION.json` 记录本地 794 个保护文件与 Git 完整状态；最终远端检查还覆盖 744 个保护输入和 2 个原 CPU build 文件。','',
'已生成要求的 11 张图，并补充 `synthetic_curved_refinement.png`。所有图标记 TEST_ONLY / NOT EXPERIMENTAL CONCENTRATION。', '',
'同步证据：`provenance/SYNC_MANIFEST.json` 与 `validation/REMOTE_LOCAL_SHA256_CHECK.json`；报告中的同步完成声明以最终 SHA256 回执为准。', '',
f"本地：`{S}`",f"远端 results：`{status['remote_stage']}`",f"远端 work：`{status['remote_work']}`",'',
'没有 Git add/commit/push/checkout/reset/clean。没有自动进入近壁理论、adhesion、RBC、GPU 或 large-N 阶段。']
out.write_text('\n'.join(lines)+'\n')
(S/'README.md').write_text('# GEOMETRIC_WALL_SLIDING_V0 — BLOCKED\n\nSTOP_WALL_NORMAL_AMBIGUITY in LONG_TRANSPORT (MPI1 and MPI4). Do not resume simulation or relax geometry gates as part of this completed stop-and-archive action. See report/MICROBUBBLE_VASCULAR_WORKFLOW_V0_WALLSLIDE_REPORT.md and WORKFLOW_STATUS.json.\n\nCopied preparation/audit scripts are retained for provenance; scripts inherited from BCFLUX are not a new instruction to run earlier stages. Authoritative executed configs are runs/*/case.cfg. Frozen source is src/, direct donor snapshots are provenance/bcflux_source/.\n')
diff=[]
for p in sorted((S/'src').iterdir()):
 old=B/'src'/p.name;a=old.read_text().splitlines(True) if old.exists() else [];diff.extend(difflib.unified_diff(a,p.read_text().splitlines(True),fromfile='bcflux/src/'+p.name,tofile='wallslide/src/'+p.name))
(S/'provenance/SOURCE_DIFF.patch').write_text(''.join(diff))
(S/'validation/OLD_STALL_COMPARISON.json').write_text(json.dumps(dict(status='PASS',LOW=dict(old_time_s=oldlow['time_s'],new_time_s=low['time_s'],resolved=True),MEDIUM=dict(old_time_s=oldmed['time_s'],new_time_s=med['time_s'],resolved=True),same_seed_source_stream_flux=True,only_executed_config_change='stage absolute path'),indent=2)+'\n')
print('REPORT_WRITTEN',out)
