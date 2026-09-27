# 一句话结论

新 P9-A.4 direct Poisson thinning 已通过入口联合分布、全量 worker 重放和 30 泡 smoke 审核，达到 **P9A4_POPULATION_INLET_READY_FOR_LARGE_SAMPLE_REVIEW**。绝对时间使用明确的名义稳态浓度假设；本轮没有启动正式 500，也没有 push 或 merge。

实际基线 `abb3ab05fb8dcddcf120765702b682f23a63a71d` 与预期相同。新 worktree 为 `/home/lzy/projects/ulm_particle_population_inlet_p9a4`，分支 `dev/p9a4-poisson-finite-size-flux-inlet-20260925`。旧科学源码、流场、合同、报告和输出均未修改。

# 一秒钟里微泡是怎么来的

`lambda = C Q`：入口每秒通过的体积 Q 乘以每立方米中的泡数 C，给出平均 source proposal 数。本例 C 为 D≤4 µm 群体的 `8.41287899745e+12 m⁻³`，Q 为 `1.551359160440232e-14 m³/s`，所以平均每秒产生 `0.130513969` 个候选泡。它是期望值，不是每秒强行放入固定数量。

# 为什么时间不是等间隔

每次间隔独立取 `-log1p(-U)/lambda_source`，然后按 source ID 顺序累加。100,000 次实测平均间隔 7.648288 s，方差 58.533982 s²。固定窗口 source Fano=1.002337，entering Fano=1.009734，仅作描述性检查。旧 deterministic scheduler 完全保留。

# 一颗 source bubble 怎么决定大小

唯一正式 source 为 `SONOVUE_D_LE_4UM_CONDITIONAL`；使用冻结 SonoVue 直方图和原分箱内均匀分布，条件化而非把超过 4 µm 的泡裁到 4 µm。每个 source ID 只抽一次 D，范围 0.75–4 µm。没有 4 µm 堆积。新合同 SHA：`345d529f216df9aefc867ee17b557202ef8ee74c21e2b82636ce0ab5e702e48f`。

# 为什么入口位置按血流量抽

直接复用原 `InletFluxSampler`，在整个真实 FEM INLET cap 上按 `(u·n)+/Q_in` 抽一次。原 P1 正向通量裁剪与三角采样不改，不预先限制到某一尺寸的可行区域，不用未来出口标签。

# 为什么只检查一次

现实中的泡不会因为这个位置放不下就瞬间传送到另一个位置。每个候选泡的时间、D、位置生成后，调用原 `FiniteSizeAdmission.check(active={})` 一次。失败就保留一次拒绝；不重抽位置或尺寸，也不延后重新插入。这里只代表独立单泡轨迹群体，不包含 simultaneous many-body pair exclusion。

# rejected source proposal 是什么

它是到达入口、但不满足真实 WALL 或原 lower-handoff 条件的物理候选事件，**不是 solver failure**。全部 100,000 行记录保存在 `data/inlet100k/proposal_ledger.jsonl`：19,221 accepted，80,779 rejected；拒绝率 80.779%。其中 WALL_REJECTED=80,672，WALL_NEARFIELD_REJECTED=107。每行都有 source ID、时间、D、位置和原因；只有 accepted 才有连续 particle ID。

# source distribution 和 entering distribution 为什么不同

大泡需要更大的空间，所以实际进入群体向较小直径偏移。必须分别解释两种 PDF。

| 分布；单位 µm | mean | p10 | p50 | p90 |
| --- | --- | --- | --- | --- |
| source | 2.046020 | 1.306021 | 1.927612 | 2.992582 |
| entering | 1.591172 | 1.158126 | 1.543728 | 2.101201 |

独立控制样本预测的 entering CDF 与实测最大差为 0.004869；16 个尺寸×位置联合单元的最大差为 0.005994。这些统计验证支持目标联合分布，不等于对任意分辨率分布的一致性证明。

# finite-size accessible flux 是什么

Q_acc(D) 是尺寸 D 的球在真实壁距与 handoff 规则下能进入的那部分流量。它影响进入概率，不参与预先改变 source 位置分布。独立 20,000 个全通量位置估计如下，800 个位置/尺寸组合与原 checker 的判断完全一致。

| D (µm) | Q_acc/Q_in | 点态 95% Wilson CI |
| --- | --- | --- |
| 0.5 | 0.84060 | [0.83546, 0.84561] |
| 0.8 | 0.66055 | [0.65396, 0.66708] |
| 1 | 0.55235 | [0.54545, 0.55923] |
| 1.2 | 0.45815 | [0.45125, 0.46506] |
| 1.6 | 0.29030 | [0.28405, 0.29663] |
| 2 | 0.16215 | [0.15711, 0.16732] |
| 2.6 | 0.03465 | [0.03220, 0.03728] |
| 3 | 0.00175 | [0.00126, 0.00243] |
| 3.5 | 0.00000 | [0.00000, 0.00019] |
| 4 | 0.00000 | [0.00000, 0.00019] |

零次接受的 MC 点不单独证明数学上严格为零；保留其上置信限。表中 CI 为逐尺寸点态区间，使用共同位置，因此不同尺寸估计相关。

# 和 Method B 有什么区别

B 固定 anchor 后反复抽原完整 SonoVue 尺寸，accepted PDF 是位置条件化分布。旧 replay 保持。本次 500 个入口诊断请求得到 306 个 accepted，共 108,361 次尺寸抽样；这是 legacy 尝试预算下的统计，不是 physical Poisson acceptance。

# 和 Method C 有什么区别

C 对全入口不可能的尺寸重抽，随后为可行尺寸重抽位置；其时钟定义为 deterministic entering `C_MB Q`。本次 500 个入口诊断请求得到 500 个 accepted，542 次 source 尺寸抽样和 859 次位置抽样。P9-A.4 保留 source rejection，进入率自然降低。B/C 的 500 只是入口统计，未运行 500 条轨迹。

# open inlet 为什么不是墙

真实 INLET 的 179 个三角不在 WALL BVH 中，cap/WALL 面重叠数为 0。球心位于 cap，球可以向上游跨出；不要求整个球进入体积。开放/封闭 cap 对照、rim overlap、真实 cap 场查询及边/顶点 owner tests 均通过。没有 inlet-plane 距离≥radius 的人工限制。[几何证据](OPEN_INLET_AUDIT_ZH.md)。

# absolute concentration 是不是 ground truth

不是。历史 8.5e12 m⁻³ 来自全泡 bolus dose/blood-volume 名义锚点，不是 D≤4 µm 浓度，也不是连续灌注实测值。新 C_modeled = 8.5e12 × F_original(4 µm)，F=0.989750470289。已澄清 total 与 conditional 语义；把这一锚点维持为恒定稳态属于 `MODEL_ASSUMPTION`。[追溯与文献](CONCENTRATION_PROVENANCE_ZH.md)。

# 当前发现了哪些真实代码 bug

在本轮指定科学路径和永久测试覆盖范围内，CONFIRMED_BUG=0，修复=0，existing scientific files 修改=0。deterministic timing、retry、旧 flow contract 都按 legacy/model scope 处理，没有把新模型偏好包装成旧 bug。新开发测试/分析脚本的执行问题及一个旧回归工具的相对路径限制另列于 [ledger](BUG_FIX_LEDGER_ZH.md)，没有隐去失败日志。未对所有 legacy 任意输入作无 bug 保证。

# NEW Network-H0 inlet audit

实际 NEW SHA `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`；OLD SHA 被新入口拒绝。100k source 的物理时钟跨度 764828.768732 s（约 8.8522 天），这是恒定模型的大样本统计窗口，不是建议一次真实 bolus 实验运行这么久。

实测进入率 0.025131115 s⁻¹，独立预测 0.025536308 s⁻¹，预测 MC 95% CI [0.025049982, 0.026022633]。两者差为合并标准误 1.514 倍。accepted 平均间隔 39.791310 s。acceptance 95% CI [0.189780, 0.194664]。

accepted clearance mean/p10/p50/p90 为 258.390, 39.048, 215.544, 543.684 nm。按直径和入口三角的拒绝统计见 `data/rejections_by_diameter.csv`、`data/rejections_by_inlet_triangle.csv`。1/3/6 worker 的完整 100k ledger 和 19,221 个 births 逐字节一致；拒绝事件 checkpoint/restart、重排合并和重复 ID 拒绝测试通过。

# 30 smoke

首 30 个 accepted births 在独立 Vast.ai 目录、6 CPU worker、每库单线程下运行，P9-A.1 动力学源码冻结。结果 **28 completed / 2 stationary / 0 solver fail**。stationary ID 4、8 有连续 32 步 rank-3 接触与 KKT 支持，不解释成生理滞留。穿透、handoff 违规、入口逃逸、NaN/Inf 均为 0（按原 roundoff 容差）。

MB 出口 O1/O2/O3=0/6/22；同一已接受样本的 point basin 为 1/12/14，noexit=3（30 s point horizon）。只用于描述，既不证明 population split，也不用于 quota 或接受门控。[完整 smoke 审核](SMOKE30_REVIEW_ZH.md)。

# P1 local conservation limitation

NEW P1/P1+VMS 合法 root section max≈3.454212%、RMS≈1.789924% 的局部守恒偏差保留。本轮只使用真实 FEM INLET，不修改 CFD、Taylor-Hood、grad-div 或 P9-A.3 internal section。

# 是否可以进入正式 500

**建议人工审核 source/entering PDF、80.779% rejection 的物理含义、名义浓度与进入率、两条 supported stationary 及 bug audit 后，再决定启动正式 500。** 本阶段仅 `READY_FOR_P9A4_LARGE_SAMPLE_REVIEW`，没有自动授权或执行 500。最终 165 tests pass、0 fail、0 skip；46,670 个旧发布文件 SHA 不变。

[图件总览](OPEN_RESULTS.html) · [数学合同](P9A4_MATHEMATICAL_CONTRACT_ZH.md) · [代码审核](CURRENT_CODE_LOGIC_AUDIT_ZH.md) · [复现](REPRODUCE.md) · [机器摘要](data/final_summary.json)
