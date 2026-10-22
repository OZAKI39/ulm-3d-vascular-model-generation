# 原生 P2/P1 表示验证与曲面装配停止点

2026-10-06。**P2 的二次场表示优势已证实，但实际 P2 CFD 对微泡轨迹的价值仍为 NOT_EVALUABLE。** A、B、C、D1 通过；D2 直接复现原生曲面混合装配不一致，触发 Phase E 9.1 的运行前停止条件。没有启动 GPU pipe/BraVa CFD 或粒子轨迹。

| 问题 | 结果 |
|---|---|
| native svMultiPhysics P2/P1 fluid path 是否真的存在并可运行？ | 原生基函数和3D单元装配确实可执行；完整 P2 GPU CFD 未运行，不能称 production 可用 |
| P2 是否解除 388769-like P1 单自由度/rank限制？ | 是。解析 fixture 中 P1 的 G rank=1，P2 的 G rank=3，精确复现二次场 |
| 真实 138584 topology 在P2下是否仍可能全零？ | 不再被拓扑强制全零：有1个自由边中点、3个标量速度DOF；但单节点控制仍使 G rank≤1。实际解是否非零未知 |
| quadratic Poiseuille exact-P2 interpolation 是否达到预期？ | 是。960原探针+4096随机点，u/G/s/curl 均达到预声明舍入误差门槛 |
| actual P2/P1 pipe CFD 是否比actual P1更可靠？ | NOT_EVALUABLE；曲面混合几何语义未通过，未运行 |
| 29.8/59.7 µm discrete-solution差异是否改善？ | 未知；只证明解析插值精确，未生成实际 P2 离散解 |
| 是否运行BraVa P2 pilot？为什么？ | 否。前置曲面装配门槛阻断，pipe CFD 也未运行 |
| 若运行，ID1189 local continuity / G/shear 是否改善？ | NOT_RUN，不声称改善 |
| ID822固定空间急增是否变得更稳定？ | 未重新评估。原 M0/medium/fine 增速36.228/4.507/16.858 mm/s的非稳定性结论不变 |
| 是否值得进入P2 particle sampler/trajectory阶段？ | **NOT_EVALUABLE**；不能凭表示优势越过曲面实现阻断点 |

## 本轮身份与旧 P2 的区别

完整 [原生路径审计](NATIVE_P2P1_PATH_AUDIT_ZH.md)、[结果身份](RESULT_IDENTITY.json)、[输入 SHA](INPUT_SHA256.json) 已保存。读取并核对实际 production retained source，commit `c3f0bb892b765b718f61069ecd9726dbc6d177fd`；生产 solver binary SHA 为 `43cd7d7b9ae9bde91275982225c6bf084859141cdaf2ec1030e1b83c3b94691d`，PETSc 身份及 build flags 保留。仅新建独立 validation 目录、测试 harness 和 CPU sampler，生产 solver/particle 源码均未修改。

退休 local P2 使用 DOLFINx、局部 Domain A、中心线网络边界；中/细档近壁 u 差约20–22%、G差约19–22%、WSS差约15.95%，bulk u 差1.30%，故停止于 `LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED`。旧方案已经使用 P2 导数，并非把 P2 导出成 P1。因此本轮的区别是现有 svMultiPhysics 原生 mixed path、原拓扑、先解析 gate 和同一生产 GPU 路线的条件性计划，不是重新包装旧采样能力。本轮没有恢复旧实现。

原18 mL/min历史校准 BraVa 场和旧1565轨迹身份不变：458/527/437出流、141观察截止、2数值失败；未把它们称为新固定灌注协议结果。未来新边界设计仍遵从中心线阻力方法；本轮历史对照不改出口压力，也未创建新的 BraVa BC。time-consistent VMS 已降级和 RT0 仅诊断的上轮结论继续有效，本轮未实施相关改动。

## A–C：原生单元、真实拓扑与采样器

原生 fixture 复用 production object files，实际执行 `nn::get_gnn`、原生 mesh/function-space 初始化、`fluid_3d_m/c` 和 `fs::thood_val_rc`。速度10节点、压力4角点；单元有34个独立未知量和40个存储槽位。局部线性场连续性 residual 误差2.776e-17。完整 GPU solve、VTU/restart round-trip 未执行；源码连接与单元验证不冒充完整算例认证。

B1 在一般仿射 tetra 上采用各分量均含混合二次项的解析速度。4106点（10节点+4096随机内部点）上，P2 的 u/G/curl/strain 最大缩放误差约5.7e-16/5.0e-16/7.5e-16/6.0e-16，低于预声明1e-11；P1 对照误差非零。B2 采用 `u=(z(1+x),z(2+y),z²)`，z=0壁面六节点精确零速，三个内部边中点保持解析值，P1 rank=1、P2 rank=3。这是表示 fixture，不是无散 CFD 解析解。

真实138584的4个角点均为 WALL。根据原 WALL 三角形的 GlobalNodeID，逐边结果如下，节点编号均为**1-based GlobalNodeID**：

| 边端点 | WALL三角形数 | 中点速度 |
|---|---:|---|
| 18691–18692 | 2 | 固定0 |
| 16114–18691 | 2 | 固定0 |
| 16114–18692 | 2 | 固定0 |
| 16115–18692 | 2 | 固定0 |
| **16115–18691** | **0** | **自由** |
| 16114–16115 | 2 | 固定0 |

对应原输入局部边[1,3]，中点 local node8。该边连接两个wall角点，但不属于任何wall三角形；未虚构内点或解除真实 BC。令这一自由节点具有任意非零速度后，原生 P2 basis 在4096内部点给出非零表示，说明它有解除全零单元的能力。与此同时，`u=N8 U8` 仍有 `G=U8⊗∇N8`，rank≤1；不能称拓扑缺陷已全部消除。真实388769也已逐边保存。详见 [图1](01_P1_VS_P2_REPRESENTATION.png) 及 [JSON/NPZ](P2_REPRESENTATION_FIXTURES.json)。

新增 [CPU/NumPy sampler](p2_sampler.py) 保留TET10全部速度DOF、TET4角点压力，返回原始u/G/curl/strain、pressure、inside_lumen和element ID。采用原生basis与解析Jacobian；多个包含单元取最小输入cell ID。封闭边界包含只允许已记录的浮点残差容差，不移动点、不截断重心坐标、不外推、不平滑、不恢复梯度。域外保留 NaN/false/-1。曲面候选边界采用 Bernstein control-point bounds；曲面反解仅通过小fixture，不等于认证整个solver的混合几何。25项 [sampler检查](P2_FIELD_SAMPLER_TESTS.json) 通过，含共享面/边/顶点、线性压力、刚体运动、二次场、有限值和保留梯度跳变。本轮没有 GPU particle kernel。 P2速度仍是C0空间，G仍可跨面跳变；采样器保留这些跳变。高阶表示不自动保证局部守恒或轨迹稳定。

## D1：冻结 V2 上的精确 Poiseuille 表示

保留36211个原角点、174785个tetra；按排序后的223119条唯一边建立算术中点，得到259330个TET10速度节点。只按native reader规则调整局部定向，cell顺序和拓扑不变。边界24248个三角面继承为6节点面，标签1/2/3不变。没有 mesher、refinement 或完整曲面网格生成。

参数保持 R=0.0005967792299824433 m、L=0.00716135075978932 m、Umean=0.1 m/s、rho=1056 kg/m³、mu=0.00345312 Pa·s。全部节点填解析Poiseuille，全部960个原探针和4096个随机内部点严格包含。最大物理尺度归一化误差：u=4.163e-16，G/curl=2.900e-14，s=1.984e-14，全部低于1e-11。u尺度为2Umean，导数尺度为4Umean/R；表中的相对L2另以各组真实解析值为分母。

| 壁距 µm | exact-P1 u相对L2 % | exact-P2 u相对L2 % |
|---:|---:|---:|
| 0.377 | 9.758956 | 1.139e-11 |
| 1 | 4.980433 | 3.381e-12 |
| 2 | 3.479448 | 1.915e-12 |
| 5 | 2.352807 | 8.226e-13 |
| 10 | 1.611492 | 4.516e-13 |
| 29.838961 | 0.380179 | 1.137e-13 |
| 59.677923 | 0.272039 | 7.707e-14 |
| 149.194807 | 9.455478 | 2.444e-14 |
| 298.389615 | 8.181824 | 1.963e-14 |

[图2](02_P1_VS_P2_POISEUILLE_REPRESENTATION.png) 显示u/G/s/curl全部10个原距离组，0.1µm仍仅诊断。CSV/NPZ保留图源。算术中点位于圆柱弦面内部，其解析速度一般非零，**该纯表示场不是施加wall no-slip后的CFD解**；需要曲面几何 gate 才能进入下一步。因此上述表不能替代原actual P1结果，也不能证明29.8/59.7µm的离散解差异已经解决。

## D2：停止原因来自实际原生曲面装配

仅取冻结V2中两个实际tetra：中央壁面cell143654、入口rim cell35265。分别比较算术中点和按解析圆柱投影wall中点、保持入口cap平面的单元；角点不动，共享边标签继承。两者全reference域的正Jacobian有充分证书：仿射基准det>0，`max_corner ||J0⁻¹(J−J0)||₂` 分别0.012740和0.012606，均<1；由于J扰动对reference坐标是仿射的，凸性给出全域界。不是只用几个采样点声称无翻转。

直接链接并调用生产对象里的 `fluid::construct_fluid`，对 `u=x−x_corner3`（解析G=I、div=3）截取真实局部连续性残差，没有时间推进、线性求解或替代装配。另调用原生basis/gnn解释它的控制流，并以完全相同的连续性quadrature、仅统一几何Jacobian作对照；不把quadrature误差混入该对照。

`fluid.cpp:591–607` 构造10节点速度坐标xwl与4角点压力坐标xql。动量循环先算压力再算速度，积分权重采用高阶速度Jacobian。连续性循环顺序相反（705–725行）：g=0最后执行压力gnn覆盖Jacobian；后续g因压力lShpF=true跳过该调用，权重又使用高阶速度Jacobian。压力形函数梯度始终来自四角点映射。因此原生速度几何确实是高阶，但mixed变量没有一致的曲面映射/积分权重。

| 单元fixture | 仿射对照残差相对差 | 曲面残差相对差 | 压力形函数梯度相对差 | 两映射最大物理位置差 |
|---|---:|---:|---:|---:|
| 中央wall | 1.568e-15 | **0.0302553%** | **0.270357%** | 35.3825 nm |
| 入口rim | 9.017e-15 | **0.0567942%** | **0.298948%** | 35.3825 nm |

这里“相对差”比较原生装配与一致曲面几何的同quadrature诊断；**不是实际CFD的速度/压力误差百分比**。原生残差与源码控制流预测的相对最大差≤3.53e-15，明确支持该实现定位。连续性第一个点的used/velocity Jacobian分别0.9990426066和0.9982015910，后三点均1；仿射对照全部约1。见 [图3](03_CURVED_GEOMETRY_NATIVE_ASSEMBLY.png)、[原生几何审计](P2_CURVED_GEOMETRY_AUDIT.json) 和 [逐点数据](P2_CURVED_GEOMETRY_SAMPLES.npz)。

虽然这些单元级差值不大，仍不能在用户要求“曲面语义可信”的前提下把它们自动视作合格；也不能由它们推算未运行CFD的误差，更不能用新P2分支的缺陷解释旧P1场。按Phase E 9.1，“native solver不能在这种curved TET10 geometry上可靠运行”即停止；状态为 **P2_CURVED_GEOMETRY_NOT_SUPPORTED**，范围限定为本轮所需的mixed曲面语义。未patch solver，未换成线性壁面绕过gate，未扫参数。

## 验证、资源与交付

[最终58项工程检查](TEST_RESULTS.json) 全部通过，其中“检测曲面不一致”通过表示成功定位阻断，不表示科学gate通过。原P1场、历史输入/协议、2431个retained源码文件、生产binary及PETSc身份保持。

启动前CPU配额7.68核、RTX4090、约15.1GiB空闲显存；cgroup上限约42.1GiB，接近上限的占用主要包含可回收文件缓存，不能只用memory.current判断可用内存。其它GPU任务保持原状。监督记录覆盖全部分析执行，包含两次已修正的harness失败，精确累计耗时及1秒采样进程树内存见 [receipt](EXECUTION_RECEIPT.json)。监督计算耗时不包含源码阅读、开发与报告整理，也不是耗时预测。native fixture仅1个MPI进程，BLAS/OMP=1；没有CFD，未为提高GPU利用率安排无关任务。receipt保留限额、初始内存、退出码、失败记录及最终图件SHA/300dpi核查。所有本项目监督任务已退出、TASK移除，无自动后续。

新增 GPU pipe CFD=0、BraVa CFD=0、population/46-track/sentinel/dt sweep=0、完整曲面网格=0。conditional field/run identity JSON 均明确 NOT_RUN；没有checkpoint或steady认证。三张图均300dpi，原始NPZ/CSV、输入/输出SHA、源码、native harness和复现入口均保留。[统一可视化入口](OPEN_RESULTS.html)。

## 唯一主建议

**P2P1_NOT_EVALUABLE_DUE_TO_IMPLEMENTATION_BLOCKER**。保留已验证的表示与采样证据，停在原生曲面混合装配语义阻断点；本轮不进入particle迁移或轨迹验证。
