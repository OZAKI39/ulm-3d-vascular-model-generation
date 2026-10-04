# 1. 这一轮为什么做

上一轮从未变形 clot 突然开流，第一步结束时还在晃动，应变漂移约 8.48%，不能据此认为损伤驱动已经平衡。
本轮先问“剩余的力是否足够小、运动是否停下、拉伸应变是否稳定”，并保持材料、mean flow、机械 dt、阻尼和侵蚀参数不变。
**目前只完成了无需 ADR 论文的验证部分；整轮任务尚未完成，状态为 NEED_PAPER_INPUT。**

# 2. 我们采用什么已有方法

- **Kilic & Madenci ADR：未实现。** 指定论文 DOI 为 [10.1016/j.tafmec.2010.08.001](https://doi.org/10.1016/j.tafmec.2010.08.001)。出版社接口只返回元数据，作者机构页面只有摘要，OpenAlex 与 Semantic Scholar 未提供开放全文；本地还检查了 24 个非生成、去重 PDF 的前两页，没有找到该文。检索证据见 `resources/LOCAL_PDF_SEARCH.json`、`resources/openalex`、`resources/semantic_scholar`、`resources/elsevier`；没有根据二手文章或记忆填写虚拟质量、阻尼或更新公式。
- **Peridigm：设计参考，没有移植求解器。** 核查 BSD-3-Clause 许可证，阅读官方 `src/core/Peridigm.cpp` 的 QuasiStatic/NOX、自由度处理和力残差，以及 ElasticCorrespondenceQSTension 算例；固定源码提交为 `f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903`。其 quasi-static 路径求解静态力平衡，显式路径使用 velocity-Verlet；damped Newton 不能冒称 Kilic–Madenci ADR。[官方源码](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/core/Peridigm.cpp)、[2023 论文 §2.6](https://doi.org/10.1007/s42102-023-00100-0)。
- **SciPy：直接调用已有 least_squares。** 使用 trust-region reflective、三点数值 Jacobian 和 linear loss，只优化 36 粒子小算例中的 72 个自由坐标；固定粒子不进入未知量。位移除以网格间距、残差除以初始外力尺度只是数值缩放，不更改力平衡根。[官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)

SciPy 的 xtol、ftol、gtol 均为 10⁻¹²，max_nfev=500，配置运行前已冻结。
独立性是“不同求解算法求同一力平衡方程”，不是第二个材料模型，也不是物理验证。
未复制整个 Peridigm 或 Trilinos；保存的许可证和示例用于审计，不作为本项目求解代码。

# 3. 什么叫“平衡”

自由粒子的 `R = F_internal + F_surface` 应接近零；两项在本项目中已是 N，不能再乘一次粒子体积。
记录 L2、RMS 和最大粒子残余力，RMS 定义为 L2 / sqrt(自由粒子数)。
相对残差为 `||R_free|| / max(||F_external_free||, 1e-12 N)`，预设门槛为 **10⁻⁴**。

对 fixed damping 同时要求物理质量计算的动能 / max(|储能|,10⁻²⁴ J) ≤ **10⁻⁶**，末段驱动漂移 ≤ **1%**。
每 50 个机械步检查一次，在这 50 步的应变范围中同时检查 bond 向量和最大/平均标量漂移；分母下限为 10⁻¹²。
需要连续 **5 次**通过全部条件，不只检查应变变化。
SciPy 是代数求根，不产生物理速度、机械积分时间或 ADR pseudo time；保存为 null 而不是假造动能史，并单独重算最终力残差。

配置见 `configs/mechanical_equilibrium_verification_v1.json`，预声明见 `provenance/PREDECLARED_PLAN.json`；没有跑完后放宽容限。
ADR 如后续实现，只能标记为 NUMERICAL EQUILIBRATION DEVICE，不能解释成血栓黏度，也不能把虚拟质量写入物理密度。

# 4. 代码改了什么

以下均为新增；旧源码、旧 78 项测试和历史结果未修改。

- `pd_clot/regularization/equilibrium.py`：提供统一平衡接口、冻结拓扑力回调、原固定阻尼法和小规模 SciPy 基准；ADR 调用明确抛出 NEED_PAPER_INPUT。
- `configs/mechanical_equilibrium_verification_v1.json`：冻结联合门槛、数值缩放和 2/4/8/16 ms 窗口设置。
- `tests/test_equilibrium_reference.py`：增加 15 项前置验证，核对参考根、固定点、静态输入、时钟及旧路径逐元素一致性。
- `scripts/verify_mechanical_equilibrium.py`：按照硬件预算并行运行 10 个无损伤算例，保存逐步历史、独立日志、状态、源码快照和哈希。
- `scripts/review_mechanical_equilibrium.py`：生成实际已完成部分的四组 PNG/PDF 和本报告。
- `verification/mechanical_equilibrium_v1/.gitignore`：保留本轮状态和静态输入文件。

未修改旧 `erosion_runner.py`；由于 ADR 依据和验证未完成，还没有创建或接入正式 equilibrated erosion runner。
硬件记录：13th Gen Intel(R) Core(TM) i7-13700HX，24 个逻辑 CPU、12 个可用物理核心，RAM 7.62 GiB，初始可用 4.67 GiB；GPU 仅记录信息，现有 kernel 仍用 CPU。
两个 sweep 都使用 3 个进程，每个 BLAS/OMP/MKL 线程数为 1；按 1 GiB/worker 和至少 20% 内存余量预估，各 job 有独立目录和日志，实际峰值 RSS 保存于 SUMMARY。

# 5. 测试

- Old tests：**78 passed**；在新增测试前单独回归，最终全量回归仍通过。
- New tests：**15 passed**。
- Total：**93 passed**；failures = 0；未删除、改弱、skip 或 xfail 旧测试。
- 新测试实际覆盖：零载荷残余力、两方法固定坐标、SciPy 独立根判据、legacy 与上一轮保存 x/v 逐元素相同、两方法 D=0/暴露时间=0、未通过条件的显式状态、缺论文时拒绝 ADR、禁止 SciPy 自动生产使用、位置相关表面载荷、DeltaN/频率元数据独立性、static resolved field，以及物理密度/输入不变。
- 尚未编写或执行 ADR 收敛、ADR 对照、ADR 迭代数独立性，以及正式 pre-relaxation/erosion runner 的损伤推进门控测试；这些不是已通过项目，也没有用 skipped tests 充数。

日志：`tests/BASELINE_78_TESTS.log`、`tests/NEW_TESTS_FIRST.log`、`tests/FULL_TESTS.log`。
窗口未通过平衡条件是被保留的数值观察，不是删除或掩盖的 pytest 失败。

# 6. 小算例独立验证

36 粒子、12 固定粒子、72 自由坐标、550 bonds；几何尺寸沿用上一轮 small fixture，材料和流场与原设置一致。
SciPy manufactured reference：**SCIPY_REFERENCE_FORCE_BALANCE_PASSED**；相对力残差 **3.71329103e-12**，L2 **1.04584382e-18 N**。
Termination reason：``xtol` termination condition is satisfied.`；SciPy nfev=6，计入数值 Jacobian 的真实残差评估次数=726，最终 Jacobian 条件数=337.996。
Static resolved-field reference 也通过，相对力残差 **5.47882345e-10**；它只用于静态力平衡，不解释成 cyclic fatigue。

| fixed window (ms) | RMS 位移差 (m) | 最大位移差 (m) | RMS差/参考最大位移 | 最大应变相对差 | 联合比较 |
|---|---|---|---|---|---|
| 2 | 1.06913e-09 | 1.98291e-09 | 1.077381e+00% | 1.202482e+00% | NOT PASSED |
| 4 | 2.00102e-11 | 4.32026e-11 | 2.016458e-02% | 5.064184e-02% | NOT PASSED |
| 8 | 2.73748e-14 | 4.78247e-14 | 2.758600e-05% | 1.156603e-05% | PASS |
| 16 | 9.66982e-19 | 2.21201e-18 | 9.744425e-10% | 0.000000e+00% | PASS |

RMS 和 L2 比较只统计自由粒子；门槛为相对残差≤10⁻⁴、最大应变相对差≤1%、RMS差/参考最大位移≤1%。
8/16 ms 在这个小问题上满足预声明对照条件，但 **ADR vs SciPy 尚未进行**；小算例不能直接替代 360 粒子的独立基准。

| 方法 | 最大位移 (m) | 平均拉伸应变 | 储能 (J) | 平均牵引 (Pa) | 最大牵引 (Pa) |
|---|---|---|---|---|---|
| SciPy | 9.9234354e-08 | 2.6655626e-05 | 7.8322095e-15 | 3.6695293 | 4.2414443 |
| 2 ms | 1.0033218e-07 | 2.6973940e-05 | 7.9036813e-15 | 3.6695288 | 4.2414442 |
| 4 ms | 9.9276071e-08 | 2.6661209e-05 | 7.8393780e-15 | 3.6695293 | 4.2414443 |
| 8 ms | 9.9234371e-08 | 2.6655639e-05 | 7.8322099e-15 | 3.6695293 | 4.2414443 |
| 16 ms | 9.9234354e-08 | 2.6655626e-05 | 7.8322095e-15 | 3.6695293 | 4.2414443 |

完整比较含 relative L2 displacement difference，见 `COMPARISONS.json`。

# 7. Fixed-window sensitivity

每个算例都从同一未变形 X、v=0、D=0 开始，dt=1 µs、damping=3500 s⁻¹；只改变窗口长度，不接续前一档最终构型。
全部窗口记录每个机械步的残余力、动能、最大/平均应变和最大位移，见各 case 的 `history.npz` 与 `history.csv.gz`。

| 粒子数 | 窗口 (ms) | 最大拉伸应变 | 相对残余力 | 最终动能 (J) | 联合门槛 |
|---|---|---|---|---|---|
| 36 | 2 | 4.3295056125e-04 | 9.0829920e-03 | 5.61796e-18 | NOT PASSED |
| 36 | 4 | 4.2802291538e-04 | 7.9473917e-04 | 3.25554e-21 | NOT PASSED |
| 36 | 8 | 4.2780631590e-04 | 4.8242827e-07 | 2.99437e-27 | PASS |
| 36 | 16 | 4.2780626642e-04 | 4.2897583e-12 | 2.81391e-35 | PASS |
| 360 | 2 | 3.2330422929e-04 | 1.7245191e-02 | 5.90680e-17 | NOT PASSED |
| 360 | 4 | 3.2759680293e-04 | 5.2803789e-04 | 8.42406e-20 | NOT PASSED |
| 360 | 8 | 3.2767573725e-04 | 5.3202194e-07 | 5.51699e-26 | PASS |
| 360 | 16 | 3.2767585692e-04 | 6.3561998e-12 | 3.48439e-34 | PASS |

360 粒子 16 ms 最终力残差 L2=4.27971230e-18 N，RMS=2.47089305e-19 N，max=5.23909993e-19 N。
本轮没有改变窗口直到“好看”；2/4 ms 的未通过结果全部保留，8/16 ms 是事先声明的两档。

图的阅读方式：

- [Fixed window convergence](figures/fixed_window_convergence.png)：看延长等待后应变稳定到哪里；左图是小算例并有 SciPy 参考线，右图是原 360 粒子基线，不能拿左图参考值去判断右图。
- [Force residual history](figures/force_residual_history.png)：看残余力是否降到 10⁻⁴ 虚线以下；横轴是原阻尼积分的机械时间，**不是 ADR iteration**，图明确标为 legacy only。
- [Reference vs fixed windows](figures/reference_vs_fixed_windows.png)：只比较同一个小算例的位移和应变，各用一个纵轴；没有画未计算的 ADR 柱。
- [No-damage startup diagnostics](figures/no_damage_startup_diagnostics.png)：看基线 clot 的动能、残余力和应变如何稳定；竖线是 legacy 条件首次连续通过，不是侵蚀时钟启动线，整段暴露时间都为 0。

四图均有同名 PDF，黑底、Arial 常规字重、英文标注；没有修改科学数据来平滑曲线。

# 8. Pre-relaxation 结果

**正式 pre_erosion_equilibration：NOT RUN。**
无损伤基线 first_converged_iteration=5400，对应 5.40 ms；这是原 fixed damping 法的诊断，不是 ADR 的虚拟时间或正式 pre-relaxation 验收。
本轮所有已跑状态 D=0、erosion_exposure_time_s=0；没有启动侵蚀时钟。
启动条件仍为：完成 ADR/小算例参考验证，并且 pre-relaxation 的力残差、运动量和驱动稳定性全部通过；随后每个损伤更新前均需通过同样门槛。

# 9. Equilibrated erosion rerun

**NOT RUN。** 因 ADR 方法来源与验证仍缺失，未写入任何声称“已平衡”的侵蚀结果。
旧 `results/streaming_erosion_demo_v1/` 完整保留；没有将旧 mean/max D 改名为本轮结果，也没有继续调用旧“未稳定仍更新损伤”的路径。
正式新 runner 的停机与保存状态分支尚待实现，本轮不宣称该分支已验证。

# 10. 时间步比较

**本轮 equilibrated erosion slow-step study：NOT RUN / NOT ASSESSED。**
未比较新的 0.02/0.01/0.005 s 侵蚀曲线，因此最大相对差异为 N/A；上一轮的 1.612081% 不属于本轮结果。
ADR residual 图、EROSION CLOCK STARTS HERE 图、equilibrated slow-step 图和最终 driver/damage maps 尚未生成；不会用占位数值或 D=0 验证图冒充正式结果。

# 11. Energy warning 是否仍存在

**ENERGY_DISSIPATION_NOT_ESTABLISHED 保持。**
本轮只做 D=0 力学验证，没有损伤跳变，故没有新的 U_before_damage、U_after_damage 或 delta_U_damage 结果。
上一轮基线累计损伤处储能增加为 2.09610704e-17 J；本轮没有修改损伤方程、NOSB 材料、integrity 更新或带符号能量审计。
机械稳定和力残差变小不能消除这一旧警告；区分启动瞬态与损伤跳变的能量贡献，需要后续真正完成 equilibrated erosion 后才可判断。

# 12. 当前状态标签

```text
TASK: PARTIALLY COMPLETED; NEED_PAPER_INPUT
ADR: NOT IMPLEMENTED — PRIMARY METHODS UNAVAILABLE
SCIPY REFERENCE: VERIFIED FOR SMALL NO-DAMAGE CASE ONLY
LEGACY WINDOW COMPARISON: COMPLETED; 2/4 ms NOT PASSED, 8/16 ms PASSED
MECHANICAL_EQUILIBRIUM: NOT VERIFIED
PRE_EROSION_EQUILIBRATION: NOT RUN
EQUILIBRATED EROSION / SLOW-STEP STUDY: NOT RUN
PHYSICAL CALIBRATION: NOT PERFORMED
MESH OBJECTIVITY: NOT ASSESSED / STILL OPEN
WALL ATTACHMENT: STILL OPEN
REAL MICROBUBBLE CFD: NOT USED
EROSION ENERGY: ENERGY_DISSIPATION_NOT_ESTABLISHED
CYCLIC_FATIGUE / DAMAGE LAW: UNCHANGED
```

需用户提供 Kilic & Madenci (2010), DOI 10.1016/j.tafmec.2010.08.001 的 PDF 本地路径或上传论文。
这是用户提示词第 1、19 节的明确前置条件；取得原文前暂停 ADR 及其依赖步骤，不凭记忆补公式，不进入新的 clot physics 开发。
