# 当前科学代码逻辑审核：最终证据

源级审查在新实现前于 2026-09-24T23:18:39.390469+00:00 完成；原始20条发现和源码SHA保存在 data/code_logic_audit.json。指定科学路径未确认真实bug，existing科学源码修改0。以下区分源码审查、执行过的测试及尚未覆盖的风险，不把legacy行为归为bug。

baseline115通过；最终portable115 + 新35 + legacy/open-cap15 =165通过，0失败/0跳过。测试命令与开发执行问题见 REPRODUCE.md、BUG_FIX_LEDGER_ZH.md。

## 1. particle_3d/src/particle_3d/injection_population.py

- **function/class**：InjectionScheduler.pop / PopulationSource.next_mb
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：整数累计通量阈值，MB 使用完整 SonoVue；独立旧 RNG STREAMS 固定不变。
- **预期行为**：旧重放保持；P9-A.4 独立 Poisson、条件尺寸源。
- **最小复现/检查方法**：固定 Q,C 对 mb_clock.time_at(1:4) 作差；恒定 1/(CQ)。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：test_poisson_inverse_transform.py; test_poisson_fixed_seed_replay.py; source-level legacy clock audit

## 2. particle_3d/src/particle_3d/injection_admission.py

- **function/class**：FiniteSizeAdmission.check / attempt
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：check 单次；attempt 同一事件重抽位置；active 实现真实 pair exclusion。WALL 使用包围球证书或精确壁距及 lower handoff。
- **预期行为**：新源只调用 check(active={}) 一次，拒绝事件永久记账。
- **最小复现/检查方法**：用始终拒绝的 checker 和计数 sampler 证明 attempt 多次 draw；新单次测试建立 call count。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：test_rejected_event_never_retries_position.py; test_rejected_event_never_retries_size.py; rejection call-count tests

## 3. particle_3d/src/particle_3d/injection_method_c.py

- **function/class**：MethodCSource.event / TruncatedSonoVue
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：全局不可行 size 重抽，可行 size 固定后条件位置重抽；确定时钟表示进入事件率。条件 CDF 实现 D≤4，无截断堆积。
- **预期行为**：旧 event 保留；新源仅复用无状态 TruncatedSonoVue.sample，不调用 MethodCSource。
- **最小复现/检查方法**：tests/particle9a2_inlet 重放、尺寸冻结、全局不可行及无 4um 堆积测试。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：test_old_method_c_replay_unchanged.py; test_truncated_sonovue_no_4um_spike.py; legacy particle9a2_inlet

## 4. particle_3d/src/particle_3d/particle82a_admission.py

- **function/class**：common_event / method_b / sphere_check
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：B 固定 anchor 后重抽全分布尺寸，已明确不是原 SonoVue accepted PDF；旧 method_c 是另一历史内移诊断。
- **预期行为**：保留两种 legacy C 的命名上下文，不将其视作 P9-A.4。
- **最小复现/检查方法**：现有 Method B 参考出生记录重放及固定 anchor 检验。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：test_old_method_b_replay_unchanged.py; 500-request inlet-only diagnostic

## 5. particle_3d/src/particle_3d/inlet_flux.py

- **function/class**：positive_pieces / InletFluxSampler / frozen_boundary_flux
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：精确 q=0 裁剪、P1 正向通量、Dirichlet 混合；GlobalNodeID -1；owner tetra 朝向使 INLET inward、OUTLET outward。非有限数据与无正向通量拒绝。
- **预期行为**：原样复用整个正向 cap，不提前条件化。
- **最小复现/检查方法**：原 flux 测试、混合正负/零顶点解析积分、P1 采样均值。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：test_flux_position_sampler_unchanged.py; complete 100k worker replay; baseline Network flux checks

## 6. particle_3d/src/particle_3d/inlet_size_capacity.py

- **function/class**：InletClearanceTree
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：距闭 WALL 集的 1-Lipschitz 上下界；全局尺寸可行和条件 proposal 供 legacy C。
- **预期行为**：本轮 direct thinning 不依赖树；诊断用独立直接 Monte Carlo。
- **最小复现/检查方法**：旧容量与条件通量 tests；P9-A.4 路径不得引用此类。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：legacy particle9a2_inlet capacity tests; source inspection confirms absent from new production path

## 7. particle_3d/src/particle_3d/sonovue_adapter.py

- **function/class**：read_sonovue / diameter_um_to_radius_m
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：冻结源码/直方图 SHA、路径边界、单位 um 到 m，执行原 sampler；单验证 MB 读 full distribution。
- **预期行为**：新源显式只有条件分布合同。
- **最小复现/检查方法**：SHA 校验及 D/r 单位测试；端点/条件 CDF parity。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：conditional CDF parity / source SHA / exact D-radius tests

## 8. particle_3d/src/particle_3d/field.py

- **function/class**：FrozenFEMField.locate / sample
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：float64 几何；闭域边界；候选 cell ID 排序选最小 owner；域外量为 NaN 并标 invalid。
- **预期行为**：原样，cap 中心必须通过该查询。
- **最小复现/检查方法**：真实 cap/共享边/顶点重复与反序调用，确认稳定 owner；当前 paired 测试。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：test_inlet_cap_not_solid.py: actual cap centroids, shared vertices and edges

## 9. particle_3d/src/particle_3d/wall_geometry.py

- **function/class**：WallGeometry.from_frozen
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：只把 manifest WALL.vtp 建入 BVH；cap 仅用于 open-rim 拓扑，不加入 solid triangle。
- **预期行为**：保留 open inlet，球可 straddle cap。
- **最小复现/检查方法**：cap/WALL GlobalNodeID 三角集合不相交；开放/封闭 cap 对照和 rim overlap。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：test_inlet_cap_not_solid.py; test_open_inlet_straddling_sphere_accepts.py; test_real_wall_overlap_rejects.py

## 10. particle_3d/src/particle_3d/validation_boundary.py

- **function/class**：ValidationBoundaryClassifier.first_event
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：有序有限三角首交；rim WALL 优先。当前轨迹实例只包含正式 OUTLET，保持原终止定义。
- **预期行为**：不改变 classifier；新 smoke 额外事后检查 INLET 逃逸。
- **最小复现/检查方法**：原 outlet 复判测试；新 smoke 根据原 segment 事后独立分类。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：unchanged original analyze.one + analyze_smoke.py all 30 saved paths; outlet matches and 0 outward INLET crossings

## 11. particle_3d/src/particle_3d/particle7_checkpoint.py

- **function/class**：write_checkpoint / read_checkpoint
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：原二进制与 sidecar/hash/schema/model identity 双重验证；恢复 pending 和全部使用 ID。新增 module 会改变全源码 identity，因此旧 checkpoint 应由原版本重放。
- **预期行为**：不放宽科学 identity；新 population checkpoint 独立保存 next source ID、已接受 ID、time 和全 ledger。
- **最小复现/检查方法**：旧 checkpoint tests；新 checkpoint 丢拒绝行、ID 跳跃、合同不符应失败。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：source-level original checkpoint audit; new test_restart_preserves_next_source_event.py; old checkpoint migration intentionally not performed

## 12. particle_3d/src/particle_3d/particle82_checkpoint.py

- **function/class**：pending_intervals / restore_stepper_recording
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：保存接受状态，恢复原二分 interval 剩余路径；拒绝 trial 无物理状态。
- **预期行为**：不与 source rejection 混淆。
- **最小复现/检查方法**：tests/particle82/test_saved_checkpoint_restore.py；新 population rejection restore 单独验证。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：source-level original saved-interval restore audit; no claim of new cross-version trajectory checkpoint replay

## 13. particle_3d/src/particle_3d/particle82a_integration.py

- **function/class**：integrate_admitted
- **分类**：CORRECT_KEEP_UNCHANGED
- **当前行为**：保存 birth/hash/config；匹配结果才能复用；先原 check，再未改 stepper；npz 临时替换后 metadata 替换。
- **预期行为**：新 output 必须独占目录；部分输出不能被认作完整 smoke。
- **最小复现/检查方法**：source-level 审核；独占新目录和已存在输出拒绝 tests；未来跨进程同 ID 写入仍标风险，不宣称有锁。
- **严重性**：NONE
- **是否修改**：False
- **永久测试/执行证据**：30 isolated outputs complete +154 returned file hashes; source-level writer review; new population partial/overwrite tests

## 14. particle_3d/src/particle_3d/particle9a_motion.py

- **function/class**：Particle9AStepper / augment_planar_system
- **分类**：SCIENTIFIC_MODEL_LIMITATION
- **当前行为**：最近局部平壁近似，open rim 回退 P6.5；现有 contact/handoff/rotation/时间推进。
- **预期行为**：按 P9-A.1 保护，所有原式和参数不改。
- **最小复现/检查方法**：当前 12 项 P9-A.1 回归；新 smoke 使用同一源码 SHA。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：portable P9-A.1 tests + 30 smoke / 14592 continuous certificates; no equations changed

## 15. particle_3d/src/particle_3d/particle81_simulation.py

- **function/class**：environment / SavedTrajectoryStepper / prepare_population
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：默认 frozen_reference/require_current_flow 是历史场；绝对 physical time 与局部积分 elapsed time 分开。
- **预期行为**：P9-A.4 必须专用 NEW loader，拒绝 OLD SHA，避免调用默认环境。
- **最小复现/检查方法**：新 OLD SHA 负向测试；旧默认入口不重写。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：test_old_flow_rejected_for_p9a4.py; test_new_contract_uses_network_h0_flow.py

## 16. particle_3d/scripts/run_particle9a1.py

- **function/class**：main / job / worker merge
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：稳定 ID；旧 runner 的 flow/ledger/model 身份各自绑定；Network runner 显式 SHA。收集完成顺序可能不同，正式 results 按 ID 排序。
- **预期行为**：旧 runner 不改；新入口不覆盖、不复用旧输出；新 merge 先完整校验 source ID，再生成连续 accepted ID。
- **最小复现/检查方法**：运行现有 paired tests；新 worker 1/3/6 和乱序 merge/重复 ID 负向测试。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：test_worker_independent_proposal_ledger.py; full 100k 1/3/6 byte equality

## 17. particle_3d/scripts/run_particle9a2.py

- **function/class**：main / job / worker merge
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：稳定 ID；旧 runner 的 flow/ledger/model 身份各自绑定；Network runner 显式 SHA。收集完成顺序可能不同，正式 results 按 ID 排序。
- **预期行为**：旧 runner 不改；新入口不覆盖、不复用旧输出；新 merge 先完整校验 source ID，再生成连续 accepted ID。
- **最小复现/检查方法**：运行现有 paired tests；新 worker 1/3/6 和乱序 merge/重复 ID 负向测试。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：test_old_method_c_replay_unchanged.py; ordered merge and duplicate-ID tests

## 18. particle_3d/reports/network_derived_flow_mb_validation_v1/scripts/runner.py

- **function/class**：main / job / worker merge
- **分类**：LEGACY_BEHAVIOR
- **当前行为**：稳定 ID；旧 runner 的 flow/ledger/model 身份各自绑定；Network runner 显式 SHA。收集完成顺序可能不同，正式 results 按 ID 排序。
- **预期行为**：旧 runner 不改；新入口不覆盖、不复用旧输出；新 merge 先完整校验 source ID，再生成连续 accepted ID。
- **最小复现/检查方法**：运行现有 paired tests；新 worker 1/3/6 和乱序 merge/重复 ID 负向测试。
- **严重性**：MODEL_SCOPE
- **是否修改**：False
- **永久测试/执行证据**：portable Network paired regression + unchanged runner used for smoke; original source SHA match

## 19. particle_3d/contracts/PARTICLE7_INLET_POPULATION_V0.json + PARTICLE9A2_PRODUCTION_V1.json

- **function/class**：concentration / flow contract
- **分类**：SCIENTIFIC_MODEL_LIMITATION
- **当前行为**：8.5e12 为 post-bolus 全体名义浓度；P9-A.2 的 conditional source 文字未体现 F(4)；旧 flow SHA 明确为 OLD。
- **预期行为**：保留旧合同；新合同用 C_modeled=C_total*F_original(4)，steady continuous infusion 是建模假设，不是该 bolus 文献所测。
- **最小复现/检查方法**：2e7/(0.0326*72)*1e6；读取冻结 CDF(4)；验证 NEW SHA。
- **严重性**：HIGH_IF_REUSED_FOR_NEW_MODEL
- **是否修改**：False
- **永久测试/执行证据**：test_new_contract_uses_network_h0_flow.py; concentration_literature_verification.json; original contracts SHA unchanged

## 20. particle_3d/src/particle_3d/particle82a_integration.py

- **function/class**：output writer concurrency
- **分类**：RISK_NEEDS_TEST
- **当前行为**：没有多写者锁；当前调度每 ID 仅一个 worker，旧 runner 外部并发同目录不在合同内。
- **预期行为**：新入口用独占新目录或完整 checkpoint 恢复，禁止重复 ID。
- **最小复现/检查方法**：未执行破坏性多写者竞争；新执行器永久测试 duplicate ID / existing output。
- **严重性**：MEDIUM_OUTSIDE_CURRENT_CONTRACT
- **是否修改**：False
- **永久测试/执行证据**：test_restart_preserves_next_source_event.py / duplicate-ID tests; simultaneous external writers NOT TESTED, risk retained

## 共同单位、边界与限制

SI: m、m/s、Pa、m³/s、number/m³；直径与半径转换在source标记处明确完成。INLET inward、OUTLET outward，GlobalNodeID→index减1，float64稳定owner；q=0/backflow按原positive clipping。新流场metadata及实文件双重SHA，禁止OLD替代与默认旧路径fallback。
永久测试覆盖inverse-transform的NaN/Inf、负/零/极小rate、CDF端点4 µm、正负通量裁剪、ID排序和拒绝事件恢复；负/非有限浓度及非正/非有限Q的构造器拒绝路径做了源级审查，未另计为已运行参数化测试。真实rim/WALL与open cap对照测试通过。active={}刻意排除many-body pair conflict；旧pair算法只做源级审查，未宣称本轮验证simultaneous人口模型。
输出单ID单worker及独占新目录是当前合同；外部同时写同目录/同ID没有锁的风险保留。旧轨迹checkpoint应在旧版本身份下重放；本轮不放宽旧身份来迁移。数学采样测试覆盖固定seed和宽松预设容差，不以出口配额判断correctness。
浓度语义已解析为full anchor×F4，但恒定稳态绝对浓度仍是假设。P1局部守恒、局部平壁动力学近似、supported stationary的生理解释属于模型限制。既有科学代码0个CONFIRMED_BUG只描述本轮证据范围，不能推导全历史代码无缺陷。
