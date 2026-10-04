# 1. 这一轮做了什么

本轮在原有工程上新增了独立的 STREAMING_EROSION 模式，未修改旧 cyclic 路径。
损伤由 mean flow 下 bond 的末端拉伸应变推动，按独立暴露时间积分。
演示参数事先固定，没有为得到断裂而调整载荷、损伤系数或阈值。
三个算例均完成 0.1 s 暴露，未出现断键或脱离。
首个 0.002 s 机械窗口没有通过预设的准稳态检查，所以这些损伤轨迹保留了尚未充分稳定的末端响应。
本轮同时保留了小算例和三档演示中损伤导致储能略增的证据，未将它隐藏为正耗散。
这是一轮软件与数值开发，不是血栓侵蚀的物理验证。

# 2. 模型公式

```text
stretch_ij = |x_j - x_i| / |X_j - X_i|
epsilon_plus_ij = max(stretch_ij - 1, 0)
rate_ij = C_t_s_inv × max(epsilon_plus_ij / epsilon_ref - threshold_ratio, 0)^m_t
D_new = min(1, D_old + erosion_exposure_step_s × rate_ij)
g = 1 - D; active iff D < 1
```

`X` 是参考位置，`x` 是实际当前位置；`stretch`、`epsilon_plus` 和损伤 `D` 都无量纲。
`C_t_s_inv = 0.01 s⁻¹` 是开发用速率，`epsilon_ref = 10⁻⁴` 是参考拉伸应变，`threshold_ratio = 0` 是无量纲阈值，`m_t = 1` 是指数。
`erosion_exposure_step_s` 是慢时间步长，基线为 0.01 s；`D_break = 1`，只有 D 达到 1 才断键。
所有侵蚀参数都是 **DEVELOPMENT / UNCALIBRATED PLACEHOLDERS**，Gc coupon calibration 不等于 streaming erosion calibration。

`mechanical_time_s` 只累计真正积分的机械窗口；`erosion_exposure_time_s` 只累计上式的慢暴露步长。
每个宏步机械窗口固定 0.002 s，机械 dt = 10⁻⁶ s，最后 20% 记录动能、最大和平均拉伸应变；以末端构型的应变更新 D，**不以时间窗 max-min 作为损伤驱动**。
漂移取“各 active bond 在检查窗内最大变化 / 窗内全局峰值应变”以及标量最大、平均应变的相对变化三者中的最大值，分母下限为 10⁻¹²，要求不超过 5%。
这只是预设的驱动稳定性检查；动能全部保存，检查通过也不等于真实物理准静态得到证明。
若未通过，则记录 `MECHANICAL_STEADY_STATE_NOT_ESTABLISHED`，仍以明确标为未验证准稳态的末端应变完成诊断积分；不延长窗口，不重调参数。

只借鉴了形变影响时间依赖退化的思路；没有移植论文的细胞模型、参数、空间腐蚀或生化过程。[Bajd & Serša (2013) 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3870808/)

# 3. 代码改了哪些文件

以下均为新增文件；旧文件修改数为 0，最终核对见 `provenance/PRESERVATION_CHECK.json`。

- `configs/streaming_erosion_demo_coarse.json`：固定几何、材料、mean-only 流场和预先声明的侵蚀参数。
- `pd_clot/regularization/erosion.py`：实现单独的拉伸应变率损伤，并将末段漂移检查与损伤驱动分开。
- `pd_clot/regularization/erosion_runner.py`：复用原力学和流场，在独立机械时钟与暴露时钟下运行并保存完整状态。
- `pd_clot/regularization/erosion_energy.py`：保留损伤前后储能变化的符号，明确报告储能增加而不冒称正耗散。
- `tests/test_streaming_erosion_damage.py`：验证解析解、不可逆性、阈值、恒定驱动步长一致性和带符号能量记录。
- `tests/test_streaming_erosion_runner.py`：验证 mean-only 与静态导入流场、双时钟、旧周期参数独立性和粒子安全约束。
- `scripts/verify_streaming_erosion.py`：执行三档固定参数计算，并独立核对保存的状态和共同暴露时刻差异。
- `scripts/review_streaming_erosion.py`：从保存的状态生成五组 PNG/PDF、浏览页和本报告。
- `results/streaming_erosion_demo_v1/.gitignore`：保留本轮状态与 VTP 文件，不修改旧忽略规则。
- `verification/streaming_erosion_v1/.gitignore`：保留本轮验证输入和失败诊断状态，不修改旧忽略规则。

新计算写在 `results/streaming_erosion_demo_v1/`；`coarse` 即 0.01 s 基线，`dt_erosion_0p01` 是指向它的相对链接，不是重复计算。
新验证、日志、比较、PNG/PDF、输入 SHA256 和本报告写在 `verification/streaming_erosion_v1/`；浏览页写在 `visualization/streaming_erosion_v1/`。
每档计算保存 CONFIG、IDENTITY、SUMMARY、history JSON/CSV、states NPZ、逐步 VTP、慢时间 PVD、source_snapshot 和 SHA256；图中每个点都来自这些保存数据。

# 4. 测试结果

- Old tests：48 / 48 passed；开始开发时单独回归通过，最终全量测试仍全部通过。
- New tests：30 / 30 passed。
- Total：**78 passed**，最终 failures = 0；命令为 `python -m pytest tests -q -p no:cacheprovider --basetemp build_streaming_erosion_v1/full_tests`。
- 覆盖零驱动、常驱动解析解、压缩不损伤、不可逆上限、D_break 强制为 1、三档慢步长解析一致性、DeltaN 与频率元数据独立性、mean-only 小算例、静态 resolved field 的末端应变驱动和双时钟。
- 初次新测试有 23 passed、5 fixture errors：同一小算例在固定位置损伤后储能增加 1.21846804×10⁻¹⁸ J，旧正耗散账本拒绝继续；原日志和失败状态完整保存在 `tests/NEW_TESTS_FIRST.log` 与 `tests/first_attempt_energy_rejection/`。
- 修复仅增加侵蚀专用的**带符号储能账本**及两项能量测试；未改材料、侵蚀参数、流场、旧账本和能量残差容限；能量增加仍会发出 `ENERGY_DISSIPATION_NOT_ESTABLISHED`，没有被截断或消除。
- 最终 50 条 warnings 中 47 条来自 VTK/NumPy 弃用提示，3 条是明确保留的准稳态或能量诊断；三档完整状态独立审计均 PASS。

测试日志：`tests/BASELINE_OLD_TESTS.log`、`tests/FULL_TESTS.log`；数值审计：`comparisons/STATE_AUDIT.json`。
恒定驱动解析图的最大绝对误差为 2.168e-19；这项结论独立于机械积分。

# 5. 数值结果

mean/max damage 指**粒子相邻原始 bond 损伤的算术平均**后的粒子均值/最大值，不等同于最大 bond D。
同一等体积 360 粒子、12592 bond、60 固定粒子的粗网格用于全部三档；基线最大 bond D = 3.27487562e-03。

| 暴露步长 (s) | 总暴露 (s) | 平均粒子 D | 最大粒子 D | 断键数 | 脱离体积 (m³) | 最大位移 (m) | 最大拉伸应变 | Settling |
|---|---|---|---|---|---|---|---|---|
| 0.02 | 0.100 | 2.91265914e-04 | 4.99650530e-04 | 0 | 0.000e+00 | 2.5538506e-07 | 3.2813194e-04 | 未通过（首步） |
| 0.01 | 0.100 | 2.91191188e-04 | 4.98908990e-04 | 0 | 0.000e+00 | 2.5541622e-07 | 3.2819039e-04 | 未通过（首步） |
| 0.005 | 0.100 | 2.91153813e-04 | 4.98538125e-04 | 0 | 0.000e+00 | 2.5543181e-07 | 3.2821963e-04 | 未通过（首步） |

**NO FAILURE / NO DETACHMENT OBSERVED**。粒子数量与总体积保持，固定粒子位置误差精确为 0；没有人为切键或粒子删除。
基线最低 J = 0.99985491，最大能量账本相对残差 = 1.378263e-04。
基线未通过的宏步为 [1]，最差漂移 8.481045% > 5%；后续窗口通过不能撤销首步失稳对已累计损伤的影响。

| 暴露步长 (s) | 实际机械累计时间 (s) | 损伤处带符号累计储能下降 (J) | 最终储能 (J) | 半径 0.525 mm 内损伤占比 | 损伤加权离源距离 (mm) |
|---|---|---|---|---|---|
| 0.02 | 0.010 | -2.09188827e-17 | 9.97502339e-14 | 0.206706 | 0.823644 |
| 0.01 | 0.020 | -2.09610704e-17 | 9.97620464e-14 | 0.206566 | 0.823780 |
| 0.005 | 0.040 | -2.09821581e-17 | 9.97679570e-14 | 0.206496 | 0.823848 |

储能下降定义为固定位置更新 D 前的储能减去更新后的储能；负值就是储能增加，不能称耗散。
该量只记录实际机械路径，没有乘以慢时间/机械时间比例，也没有采用 coupon Gc；逐 bond 分摊仅按 ΔD 比例作诊断，不是断裂能标定。
基线累计储能增加为 2.09610704e-17 J，超过原 roundoff 阈值的增加次数为 10，状态为 `ENERGY_DISSIPATION_NOT_ESTABLISHED`。
本轮没有碎片，因此没有治疗时间中的碎片输运结论；原 transport 只随实际机械窗口执行。

# 6. Erosion time-step comparison

**INCONCLUSIVE**。
三档总慢暴露都是 0.1 s，机械窗口都为 0.002 s；宏步数不同使机械累计时间分别为 0.010、0.020、0.040 s，它们不能解释成治疗时间。
在共同慢时刻 0.02、0.04、0.06、0.08、0.10 s，以 0.005 s 档作为数值参考，计算 `abs(case-reference) / max(abs(reference), floor)`。
损伤、应变、比例分母下限为 10⁻¹²，位移/距离为 10⁻¹⁵ m，能量为 10⁻³⁰ J；初始零状态不参与相对差比较。

| 指标 | 共同慢时刻最大相对差 | 最终时刻最大相对差 |
|---|---|---|
| `mean_damage` | 0.184494% | 0.038502% |
| `maximum_damage` | 1.134017% | 0.223133% |
| `damage_stored_energy_drop_J` | 1.612081% | 0.301568% |
| `stored_energy_J` | 0.575474% | 0.017764% |
| `maximum_displacement_m` | 0.716724% | 0.018302% |
| `maximum_epsilon_plus` | 1.359772% | 0.026717% |
| `fraction_of_damage_inside_influence_region` | 0.522575% | 0.101723% |
| `damage_weighted_distance_from_source_m` | 0.130216% | 0.024788% |

所有指标、共同慢时刻的最大相对差为 **1.612081%**，出现在 `damage_stored_energy_drop_J`、步长 0.02 s、暴露 0.02 s。
即使曲线接近，三档首个机械窗口均未通过预声明检查，所以按预声明规则不能据此断言慢时间积分已收敛，也不宣称网格客观性。
详表见 `comparisons/time_step_differences.csv` 与 `comparisons/TIME_STEP_STUDY.json`。

# 7. 图怎么看

1. [Erosion rate law](figures/erosion_rate_law.png)：线是恒定拉伸驱动的解析解，空心点是数值更新；正常情况是零驱动不损伤，其余点落在解析线上。当前最大误差 2.168e-19，支持率方程的软件实现。
2. [Erosion time step consistency](figures/erosion_time_step_consistency.png)：左图比较三档平均损伤，右图放大相对 0.005 s 档的差异，横轴都是慢暴露时间。当前共同慢时刻平均损伤最大相对差为 0.184494%，但首步未稳定，因此结论仍为 INCONCLUSIVE。
3. [Mechanical settling](figures/mechanical_settling.png)：看每个 2 ms 窗口最后 20% 的动能与最大应变，以及末段漂移是否低于 5% 虚线；这里横轴明确是窗口内机械时间。当前第 1 步漂移 8.4810%，后续窗口通过，所以初始瞬态影响没有被误报为已排除。
4. [Spatial erosion driver](figures/spatial_erosion_driver.png)：颜色是原始影响权重加权的相邻 bond 末端拉伸应变，白色菱形标出 clot 左侧的预设流场中心，方块是固定区；合理的图应与实际加载及锚定约束一起判断，不能仅凭热点判物理正确。当前范围 1.76403e-06–5.42733e-05，按真实位置和 1× 位移绘制，没有箭头或蓝色微泡。
5. [Final damage map](figures/final_damage_map.png)：看真实粒子损伤及色条数值，亮色只代表本图的最大值而非 D=1。当前真实 min/max = 4.83989e-06/4.98909e-04，最大 bond D 也远低于 1，未发生断裂。

五图均同时提供同名 PDF；保持黑底、Arial 常规字重、英文标注，浏览页为 `visualization/streaming_erosion_v1/index.html`。
图中的“源”是预设局部声流场中心，管流牵引仍是分布载荷；它不是求解得到的微泡表面，也不是画在 clot 上方的合力中心。

# 8. 当前科学状态

```text
STREAMING_EROSION v1: SOFTWARE / NUMERICAL DEVELOPMENT
PHYSICAL CALIBRATION: NOT PERFORMED
MESH OBJECTIVITY: NOT ASSESSED IN THIS STEP
WALL-ATTACHMENT MESH CONSISTENCY: STILL OPEN
CYCLIC_FATIGUE: UNCHANGED
COMBINED EROSION + FATIGUE: NOT IMPLEMENTED
REAL MICROBUBBLE CFD: NOT USED IN THIS VERIFICATION
PARAMETERS: DEVELOPMENT / UNCALIBRATED PLACEHOLDERS
EROSION TIME-STEP STUDY: INCONCLUSIVE
MECHANICAL_STEADY_STATE_NOT_ESTABLISHED
NO FAILURE / NO DETACHMENT OBSERVED
FRAGMENT RESOLUTION: NOT_ASSESSED_NO_DETACHMENT
EROSION ENERGY: ENERGY_DISSIPATION_NOT_ESTABLISHED
```

旧科学状态保持：SOFTWARE VERIFIED；fracture-energy coupon 为 ENERGY REGULARIZED FOR SPECIFIED OPENING PROTOCOL ONLY；旧 main streaming case 为 NOT CONVERGED；旧 N=25,000 算例为 NO FAILURE / NO DETACHMENT；旧 forced-failure case 为 FRAGMENT_RESOLUTION_INADEQUATE。
本轮未运行真实微泡 CFD、材料标定、空间腐蚀反馈、组合疲劳或生化过程；能量账本闭合、平滑图像和测试通过都不能证明局部物理精度。

# 9. 下一步建议

只做一个独立的**机械松弛窗口敏感性检查**：冻结当前 mean flow、材料和侵蚀参数，比较预先指定的较长窗口与本轮 0.002 s，检验首步末端应变能否稳定；本轮结果保持不变。
