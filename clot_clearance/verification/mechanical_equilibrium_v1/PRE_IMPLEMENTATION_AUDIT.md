# 开发前源码审计

本文件在新增平衡求解代码之前生成。当前分支为 `clot-clearance/energy-regularization-20260930`，上一轮 STREAMING_EROSION v1 尚是工作树中的新增文件；本轮把它们全部作为需要保护的既有文件。没有发现适用的 AGENTS.md。全部既有文件的身份见 `provenance/PREEXISTING_SHA256.json`。

## 现有流程

1. **力评估接口**：`fast_mechanics.evaluate_supported(cloud, x, integrity, material, safety, shape)` 返回粒子力 N、F、J、弹性能密度及稳定化能密度；`prepare_supported_shape` 根据当前 integrity 计算 shape。`ConsistentFragmentFluid.surface_force(...)` 返回表面力 N、法向和牵引 Pa。自由粒子残余力应为这两项力之和，不能重复乘体积。
2. **固定粒子**：`erosion_runner.run` 使用参考 z 坐标严格小于固定层上界（减去很小的网格容差）的粒子；每个机械步强制其 x=X、v=0，并将加速度置零。基线有 60 个固定粒子。固定粒子上的反力不应进入自由粒子的收敛范数，SciPy 也不应把它们列为未知量。
3. **表面力依赖当前位置**：暴露掩码在固定拓扑窗口内不变，但法向由当前位置重新计算；provider 在当前位置取速度梯度与压力。因此必须在每次 residual evaluation 时重新求流体表面力，不能冻结初始牵引向量。
4. **mean-only 时间性**：manufactured provider 仍读取代表频率，但当 oscillatory_amplitude=0、mean_fraction=1 时，其时变项及 du/dt 都为零；管流也是静态的。ResolvedStreamingFieldProvider 为静态输入。当前字段是软件验证场，不是真实微泡 CFD。
5. **阻尼位置**：旧 `erosion_runner.py` 每个 velocity-Verlet 步前后各乘 exp(-damping × dt/2)，附着粒子 damping=3500 s⁻¹；它影响机械迭代而不是侵蚀率。脱离粒子的速度松弛输运另有 half step，本轮 D=0、全体附着，故这项不触发。ADR 的虚拟质量和阻尼以后必须分开保存，不能覆盖物理密度或称为血栓黏度。
6. **2 ms 定义**：relaxation_time_s=0.002，dt_s=10⁻⁶，每个宏步固定 2000 次机械积分；最后 20% 保存动能及拉伸应变，用应变范围检查漂移。首步约 8.48% 超过原 5% 门槛。旧模式会保留诊断标记后继续更新损伤；新流程必须在不平衡时停止推进损伤。
7. **损伤冻结**：旧机械窗口内 D、integrity、shape 和连通拓扑不变；只有窗口结束时才用末端 tensile strain 更新 D，再更新 shape、拓扑与暴露面。`erosion.py` 不使用时间窗应变幅值、DeltaN 或 coupon Gc；本轮保留此损伤方程及 `erosion_energy.py` 的带符号能量审计。
8. **两类时间保存**：旧 mechanical_time_s 由机械 tick × dt 累计；erosion_exposure_time_s 由宏步 × exposure_step 单独累计，保存在 history 与 states。新 equilibrium_iteration / equilibrium_pseudo_time 必须另列；无损伤窗口和 SciPy 基准的侵蚀时间始终为零。

已阅读：`erosion_runner.py`、`erosion.py`、`flow.py`、`fast_mechanics.py`、`fragment_mechanics.py`、`mechanics.py` 和上一轮 `REVIEW_STREAMING_EROSION.md`。

## 文献与实现边界

- Peridigm 的 BSD-3-Clause 许可证已核查。检查了官方 `src/core/Peridigm.cpp` 的 QuasiStatic/NOX 求解及残余力范数，以及 `test/verification/ElasticCorrespondenceQSTension/ElasticCorrespondenceQSTension.xml`。借鉴力残差、边界自由度处理和独立容限的设计；不会复制其求解器，也不会 vendor Trilinos。
- 该实现中的 quasi-static Newton/NOX 与显式积分入口是不同路径；不能把 Peridigm 的 damped Newton 称作 Kilic–Madenci ADR。已查看代码中没有找到可直接用于此项目的该论文 ADR 实现。
- [Peridigm 官方源码](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/core/Peridigm.cpp)、[2023 论文](https://doi.org/10.1007/s42102-023-00100-0)、[许可证](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/LICENSE.md)。
- SciPy `optimize.least_squares` 仅作为小问题裁判，使用其已有 trust-region reflective 实现。未知量只含自由粒子位置，残差是总力向量；采用位移和残差数值缩放，不改变材料或载荷；保存 termination reason，并独立检验最终残余力。[官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)
- Kilic & Madenci (2010), DOI 10.1016/j.tafmec.2010.08.001：本地项目、references、上传附件和用户目录按文件名未找到该文；出版社返回书目元数据而无方法正文，作者机构页面仅有摘要；OpenAlex 与 Semantic Scholar 均未提供开放 PDF。检索证据保存在 `resources/`，已向用户请求 PDF 本地路径。

**当前 ADR 门槛：NEED_PAPER_INPUT。未取得完整方法前不实现 ADR 公式，不借用二手文章冒充原文。** 按提示词第 19 节，先保留并完成无损伤 fixed-window sensitivity 与 SciPy reference；ADR、正式 pre-relaxation、equilibrated erosion 和其 slow-step study 依赖文献及基准验收，不能提前推进。

## 预定最小实现

增加独立的 `equilibrium.py`，先支持 legacy_fixed_damping 与仅验证用途的 scipy_reference，并让 adaptive_dynamic_relaxation 明确拒绝执行直至文献可核对。复用现有力学 kernel，不复制或改写它。旧 erosion_runner 完整保留；取得 ADR 依据后再用新 runner 副本连接统一接口，不覆盖上一轮结果。

无损伤小算例沿用上一轮 4×3×3 fixture，其材料、网格间距、源位置、mean flow、机械 dt 和 damping 不变。小算例验证完成后，在 360 粒子基线作 2/4/8/16 ms 固定窗口比较；每档从同一未变形状态出发，独立目录、进程和日志。收敛与 SciPy 配置会在运行前冻结，不根据结果放宽。
