# 1. 这一轮为什么做

之前开流后的第一段 2 ms 还在晃动，应变漂移约 8.48%。这时直接开始记损伤，会把启动过程混入侵蚀。
现在先让力学达到平衡，之后每次损伤更新前再次检查平衡。计算机迭代次数不再计入治疗曝光秒数。
**结论：本开发算例的机械平衡与慢步比较通过；损伤能量耗散问题仍未解决。**

本报告是 `adr_completion/` 新增的完成版。上级同名报告保留为当时缺论文的历史记录，不能作为最新结论。
原 BraVa CFD、微泡批次、原 cyclic fatigue、原 erosion_v1 及全部旧结果均未改动。

# 2. 我们采用什么已有方法

直接依据用户提供的 Kilic & Madenci (2010) 全文 §2.1，印刷页 196–197，公式 (11)–(14)、(18)–(19)。
已提取 11 页正文并渲染核对公式；PDF 路径和 SHA 在 `resources/PAPER_IDENTITY.json`。
[论文 DOI](https://doi.org/10.1016/j.tafmec.2010.08.001)。

原文 (14) 规定虚拟质量应覆盖刚度绝对行和。当前 NOSB 无现成切线矩阵，故对现有 free 总力作中央差分，
逐列累加行和，包括随位置变化的流体牵引。没有使用不适用于当前 NOSB 的论文键模型公式 (15)。
虚拟质量取行和界的 1.1 倍，差分扰动 h×10⁻⁵，结束时以半扰动重新检查行和界；没有修改力学 kernel。
这里采用 (14) 的坐标系相关形式，不声称具有论文基于 (15) 的框架不变性质。

首次严格按正文写出的 probe 出现负局部 Rayleigh 商，在迭代 16 后停止；完整失败数据保留在 `initial_small_probe/`。
随后核实 **Madenci & Oterkus (2014) 作者提供的 Fortran 附录**：非正商使阻尼取零；dtau=1 时 c>2 则使用1.9。
据此补上原作者保护分支，而非擅自取绝对值。来源、文件 SHA、行号和实现前声明见
`DAMPING_GUARD_ADDENDUM.md`、`resources/SPRINGER_CODE_IDENTITY.json`。
[Springer 官方附录](https://extras.springer.com/?query=978-1-4614-8464-6)。
正式标准小算例中负商保护 64 次、上限保护 32 次；并不掩盖这些触发。

Peridigm 仅作成熟软件设计参照：核对了 BSD-3-Clause 许可、固定自由度及 QuasiStatic/NOX 残量设计、
correspondence 准静态示例；既未 vendor Peridigm/Trilinos，也未复制其求解器。
[源码固定版本](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/core/Peridigm.cpp)，
[2023 论文](https://doi.org/10.1007/s42102-023-00100-0)。
SciPy 直接调用现有 [least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)，
只作为小问题独立求根算法；不是第二个材料模型，更不是实验验证。

# 3. 什么叫“平衡”

同时满足三项：自由粒子的未平衡总力足够小、数值运动量足够小、拉伸应变不再明显变化。
每 50 次迭代检查，连续 5 次通过：相对 free 总力 ≤10⁻⁴、末 50 步 driver drift ≤1%、
伪动能诊断/|储能| ≤10⁻⁶。分母下限分别为10⁻¹² N、10⁻¹² strain、10⁻²⁴ J，最大50000步。
所有收敛容差先声明再运行，首个失败后也没有放宽。

ADR 使用虚拟时间步1及虚拟质量，按论文中央差分更新位移，并以局部刚度估计阻尼。
它是 **NUMERICAL EQUILIBRATION DEVICE**，不能解释为 clot 黏度。
物理密度始终1056 kg/m³；pseudo half velocity 不当作物理速度输出，physical kinetic energy 为 null。
ADR 全过程 D、integrity、拓扑、曝光和周期数均冻结。具体公式、单位和局限见 `METHOD_AUDIT.md` 及保护补充。

# 4. 代码改了什么

**此前已有文件修改数：0。** 新增文件：

- `pd_clot/regularization/equilibrium_adr.py`
- `pd_clot/regularization/erosion_equilibrium_runner.py`
- `configs/streaming_erosion_equilibrium_v1.json`
- `configs/streaming_erosion_equilibrium_guarded_v1.json`
- `tests/test_equilibrium_adr.py`
- `tests/test_erosion_pre_relaxation.py`
- `scripts/verify_adr_completion.py`
- `scripts/audit_adr_completion.py`
- `scripts/review_adr_completion.py`

新接口支持 ADR，并转发旧 legacy 和小规模 SciPy 路径；旧 `equilibrium.py` 保留为当时缺论文的历史快照。
新 runner 每个阶段保存 x、D、力、逐迭代历史、虚拟质量及独立时钟。最后一次损伤更新后还做一次平衡，
仅用于最终状态和图，曝光仍停在0.1 s。任一未通过的阶段保存状态并停止损伤推进。

硬件：13th Gen Intel(R) Core(TM) i7-13700HX，24逻辑CPU，12可用物理核心，
RAM 7.62 GiB。当前内存预算使两个新 sweep 各用
2 / 2 个进程；每进程BLAS/OMP/MKL线程数1，
按1 GiB/worker并保留至少20%可用内存估算。逐job wall/CPU/peak RSS 已保存；无GPU重写。
旧窗口上一阶段已按冻结条件并行重跑，本次复用其不可变结果，没有重复耗费计算。

# 5. 测试

旧78项全部保留通过，前置独立参考阶段新增15项保留通过，本次新增19项通过。
**全量 112 passed，0失败、0跳过。**
覆盖了小问题解析解/独立参考、残余力下降、固定点、输入不变、质量界、作者保护分支、
额外迭代、静态导入场、旧路径逐元素复现、初始D和曝光为零、周期元数据独立、
初始和中途失败阻止后续损伤更新、同坐标 signed energy 重算。
日志见 `NEW_TESTS_FIRST.log`、`FULL_TESTS.log`。
另对所有 41 个保存平衡状态重算了残余力，
对所有35次正式损伤跳变重算了能量及损伤增量，见 `STATE_AUDIT.json`。

# 6. 小算例独立验证

36粒子，12固定，72自由坐标，550 bonds；材料和mean flow均未调整。
SciPy：`xtol` termination condition is satisfied.，相对力残差 3.7132910e-12。
ADR：350次迭代，相对力残差 5.4549430e-11。
最大应变相对差 1.0380615e-12；RMS位移差/参考最大位移
5.3174004e-10，均小于预设1%。
RMS位移差 5.2766879e-17 m，最大位移差 7.3197178e-17 m，
relative L2位移差 1.0423680e-09。

| 方法 | 最大位移(m) | 最大拉伸应变 | 平均拉伸应变 | 储能(J) | 牵引均值/最大值(Pa) |
|---|---|---|---|---|---|
| ADR | 9.92343539e-08 | 4.27806266e-04 | 2.66556260e-05 | 7.83220951e-15 | 3.669529 / 4.241444 |
| SciPy | 9.92343539e-08 | 4.27806266e-04 | 2.66556260e-05 | 7.83220951e-15 | 3.669529 / 4.241444 |
| 2 ms | 1.00332176e-07 | 4.32950561e-04 | 2.69739403e-05 | 7.90368131e-15 | 3.669529 / 4.241444 |
| 4 ms | 9.92760708e-08 | 4.28022915e-04 | 2.66612092e-05 | 7.83937800e-15 | 3.669529 / 4.241444 |
| 8 ms | 9.92343712e-08 | 4.27806316e-04 | 2.66556395e-05 | 7.83220989e-15 | 3.669529 / 4.241444 |
| 16 ms | 9.92343539e-08 | 4.27806266e-04 | 2.66556260e-05 | 7.83220951e-15 | 3.669529 / 4.241444 |

额外500次迭代和static resolved-field变体也通过；完整对照在 `SMALL_REFERENCE_GATE.json`。
小问题的参考值只用于同几何比较，不拿来替代360粒子问题的参考解。

# 7. Fixed-window sensitivity

全部从同一未变形状态、v=0、D=0出发，dt=1 µs、原阻尼3500 s⁻¹，不接续前一档构型。

| 粒子数 | 窗口(ms) | 最大拉伸应变 | 相对力残差 | 联合条件 |
|---|---|---|---|---|
| 36 | 2 | 4.3295056125e-04 | 9.0830e-03 | NOT PASSED |
| 36 | 4 | 4.2802291538e-04 | 7.9474e-04 | NOT PASSED |
| 36 | 8 | 4.2780631590e-04 | 4.8243e-07 | PASS |
| 36 | 16 | 4.2780626642e-04 | 4.2898e-12 | PASS |
| 360 | 2 | 3.2330422929e-04 | 1.7245e-02 | NOT PASSED |
| 360 | 4 | 3.2759680293e-04 | 5.2804e-04 | NOT PASSED |
| 360 | 8 | 3.2767573725e-04 | 5.3202e-07 | PASS |
| 360 | 16 | 3.2767585692e-04 | 6.3562e-12 | PASS |

360粒子 ADR 无损伤平衡 driver=3.2767585700e-04。
2/4 ms 未通过的历史数据完整保留；8/16 ms趋向ADR平衡，图中小算例另有SciPy线。

# 8. Pre-relaxation 结果

360粒子基线从未变形状态开始，**450次ADR迭代**后通过全部门槛；这是虚拟迭代，不是450秒治疗。
最终相对残差 1.2029018e-09，L2=8.0992945e-16 N。
整段D=0、erosion exposure=0。只有连续5次联合检查通过之后，才允许第一次非零曝光损伤更新。
图 `pre_relaxation_driver` 明确标出 EROSION CLOCK STARTS HERE。

# 9. Equilibrated erosion rerun

新结果在 `results/streaming_erosion_equilibrium_v1/`，原v1结果不覆盖。
C_t=0.01 s⁻¹、m=1、epsilon_ref=10⁻⁴、threshold=0、D_break=1，mean flow及材料不变。
基线慢步0.01 s，10次更新，曝光0.1 s；更新后最终平衡不增加曝光。
最终粒子平均D=2.9109696751e-04，最大粒子D=4.9812387726e-04，
最大bond D=3.2793352033e-03；broken bonds=0，
detached fraction=0。
**没有观察到断裂或脱离**，这不是“模型失败”，也不能称为已验证断裂功能。

最终相对总力残差 1.1133051e-09；
L2=7.4961626e-16 N，RMS=4.3279115e-17 N，
max=1.0611346e-16 N。
现有ADR范围是锚定且拓扑不变的开发算例；遇到断裂导致范围改变会保存并停止，不冒充脱落后的物理动力学。

# 10. 时间步比较

三档曝光都严格结束于0.1 s，ADR迭代没有混入秒数。
比较所有共同曝光时刻的六项预声明指标：平均/最大粒子D、最大bond D、最大driver、最大位移、储能。
以0.005 s为参考，最大相对差 **0.0131107%**，低于预声明5%；
仅说明这三档、本算例的开发一致性，不构成网格客观性或一般时间收敛阶证明。

| 慢步(s，p表示小数点) | 平均粒子D | 最大粒子D | 最大bond D | ADR总迭代 | 六指标最大相对差 |
|---|---|---|---|---|---|
| 0p02 | 2.9107751636e-04 | 4.9808057792e-04 | 3.2790485822e-03 | 1700 | 0.0131107% |
| 0p01 | 2.9109696751e-04 | 4.9812387726e-04 | 3.2793352033e-03 | 2950 | 0.0043709% |
| 0p005 | 2.9110669758e-04 | 4.9814553379e-04 | 3.2794785449e-03 | 5450 | 0.0000000% |

六组图均为黑底、Arial常规字重、英文、无箭头，PNG和PDF各一份：

- [Fixed window convergence](figures/fixed_window_convergence.png)：等更久后稳定到哪里。
- [Force residual history](figures/force_residual_history.png)：残余力如何降至门槛以下。
- [ADR versus reference](figures/adr_vs_reference.png)：同几何不同算法的位移与应变。
- [Pre-relaxation driver](figures/pre_relaxation_driver.png)：何时才允许启动曝光。
- [Erosion time step consistency](figures/erosion_time_step_consistency_equilibrated.png)：平衡后的慢步曲线。
- [Final driver and damage maps](figures/final_driver_and_damage_maps.png)：真实1×位移，各色条使用各自实际min/max。

driver图为邻接bond拉伸应变的influence权重平均；粒子D是邻接bond损伤的算术平均，
都沿用已有聚合定义，不把两者的色值混为一谈。颜色亮不代表D=1。

# 11. Energy warning 是否仍存在

**仍存在：ENERGY_DISSIPATION_NOT_ESTABLISHED。** 每次保存同一x下的U_before、U_after及
delta_U=U_before−U_after，保留负号；没有将虚拟阻尼耗散记成物理耗散，没有修改damage law。

| 曝光步长 | 储能增加次数/更新数 | 累计signed delta_U(J) |
|---|---|---|
| 0p02 | 5 / 5 | -2.101751147e-17 |
| 0p01 | 10 / 10 | -2.101047404e-17 |
| 0p005 | 20 / 20 | -2.100689527e-17 |

在更新前已经通过机械平衡时，损伤跳变仍导致储能增加，因此这一现象不能仅归因于启动振荡。
它表明当前NOSB完整性更新的热力学耗散尚未建立；不因此否认已验证的力平衡，也不借平衡通过掩盖能量问题。
本轮不修复它，下一轮也没有自动启动。

# 12. 当前状态标签

```text
MECHANICAL_EQUILIBRIUM: VERIFIED_FOR_CURRENT DEVELOPMENT_CASE
ADR: VERIFIED FOR CURRENT DEVELOPMENT CASE, WITH DOCUMENTED AUTHOR GUARDS
SCIPY REFERENCE: PASSED FOR SMALL NO-DAMAGE CASE
EROSION UPDATES: ALL 35 GATED BY EQUILIBRIUM
SLOW-STEP COMPARISON: PASSED PREDECLARED DEVELOPMENT ACCEPTANCE
EROSION ENERGY: ENERGY_DISSIPATION_NOT_ESTABLISHED
FAILURE / DETACHMENT: NOT OBSERVED
PHYSICAL CALIBRATION: NOT PERFORMED
MESH OBJECTIVITY: NOT ASSESSED / STILL OPEN
WALL ATTACHMENT: STILL OPEN
REAL MICROBUBBLE CFD: NOT USED
```

本轮交付限于这一机械平衡问题。没有增加材料、壁面、疲劳、声流或断裂物理。
源码、配置、历史数据和图的数据来源均可追溯；保留检查见 `provenance/PRESERVATION_CHECK.json`。
