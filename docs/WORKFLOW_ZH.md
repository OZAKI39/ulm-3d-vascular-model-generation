# 科学流程与版本边界

2026-10-22 更新：当前仅保留 ROI-only 正式结果。下述 H0/best-feasible 比较是方法来源说明及初次整理时的历史核验，相关整套旧结果已按用户指令删除；当前安装与测试状态见 SOLVER_INSTALL_AND_CLEANUP_ZH.md。

## 几何和网络

实际源样本为 `fMOST_0_5_6_0_0_6_0001_02_01`，不是 Schmid/NW1 模型。所用处理记录是 `all_run_20260825_133152`、`20260825_133201_radius_plus_structure_k5` 与 `ultraliser_anchor003274_20260825_133350`。原始 SWC、图像/掩模、采样结果、ROI 的源节点/半径、表面和元数据均在 `inputs/geometry/`。

当前几何函数从实际引用的 `utils.rodent_vasculature`、`sampling`、`cfd_lumen`、`cfd_preprocess`、`cfd_surface_prepare` 递归提取；打印支撑、模具和 BraVa 处理不构成新工程依赖。几何入口 `retained_roi()` 实际读取保留 ROI，已验证 anchor 与样本身份。原 ROI 裁切端口记录和 CFD 延长/封口后四个流体端口不能仅按数量混同；完整转换由保存的映射、面标签和延长段报告追踪。

1D/0D 使用原始半径、SI 长度和 `linear_radius_resistance`。`roi_only_hydraulics.load_roi_cache` 从完整 A-H0 图中抽取已标识的内部边，保留源索引和流向系数；固定入口体积流量、给定三个出口压力后求解节点守恒。`regress_h0` 重构保存的 H0 端口状态；`roi_boundary_design.design` 保留直接压力基底解、迭代交叉校验与 gauge 不变性检查，没有使用 CFD 反向校准。

H0 是全网络理想化边界模型；best-feasible 是其受限设计结果；ROI-only 是独立的 ROI 出口压力设计。`parameterized_fem_handoff` 与 `fem_h0_case` 保留真实切口到 CFD cap 的延长段压力转移与 XML 边界序号。三个方案的 3D 出口比例不能按 0D 目标覆盖。

## FEM、冻结场和 WSS

网格生成保持原 SimVascular TetGen API 与 frozen mesh policy。已有 CFD 为 P1/P1 + VMS；Python 负责准备、启动、监控和验收，不重写原生求解器。原生产过程以零初值开始，按原 policy 检查场变化、流量和连续稳态区间，通过 STOP_SIM 停止协议正常结束，并校验 native checkpoint。`solver_support/native_run.py` 保留这些步骤，而不是仅等待可执行文件退出。

`flow_field.load_case` 验证相应场的 SHA、规范网格点序/单元顺序、GlobalNodeID 和边界面，再复用原 `FrozenFEMField`。Stage-Q 曾是较低入口速度历史参考，其整套结果现已删除；默认材料/几何相同不代表它与 2 mm/s 场是同一个场。

WSS 采用原 `solver_support.wss` 与 `wss_case`。计算使用 `mesh_arrays.npz` 中的原始 float64 坐标，保留法向、父单元、面积、切向牵引和原节点面积加权显示值；本轮在 ROI-only 上重新只读计算，45,221 个壁面面的保存量与结果逐项比较。原图中使用的显示平均不能替代原始 facet WSS，也不改变物理量。

## 源、准入与有限尺寸运动

SonoVue 原 histogram 与 sampler contract 都有原 SHA。新的 `sonovue_adapter` 用显式 import 加载工程内唯一 sampler；原 sampler 的身份字节作为数据保存，验证时只读取/哈希，绝不 exec。D≤4 µm 截断仍是原分布条件化，浓度由原名义浓度乘保留质量；它是模型假设，非本轮测得的灌注浓度。

`ContinuousInfusionSource` 沿用按 source event ID 分角色的 PCG64、Poisson 到达间隔、一次粒径抽样、一次完整入口流量加权位置抽样、一次有限尺寸准入。失败的 source event 保留；仅接受事件分配 particle ID。`source_for_case` 复用这些数值规则和保存种子，并赋予新的代码身份，不冒充原 source identity。本轮对 ROI-only/best-feasible 的最小源事件与原 ledger 实际对比。

`formal_dynamics_p9a5`、`particle9a_motion`、`particle65_motion` 等旧名称仍是真实有效的模块。动力学保留有限尺寸、近壁阻力、接触约束、连续 handoff 证书、细分和事件顺序；名义步长分别为 0.5 ms 或 1 ms，观察扩展为 3/6/12 s。出口穿越由原面分类器确认；接触支持静止需要原支撑/速度条件，观察截止与失败不能改作静止。

ROI-only 和 best-feasible 的正式生产使用私有函数绑定中的 C++ 最近点内核；`native_adapter` 的有限私有 import 代理只替换当前包内 handoff 查询函数，保留原绑定方式和编译浮点选项。它不动态加载旧目录代码，不修改全局 `sys.modules`。本轮在真实冻结 WALL 三角形上逐字节对比 Python 与 C++ 最近点/重心坐标结果。

RBC 相关基础类在真实闭包内保留；当前正式单泡模型没有 MB–MB 或 RBC 耦合。LAMMPS bridge 的存储/邻居/重启用途不等于 LAMMPS 代替正式粒子动力学。

## 身份、保存与展示

原轨迹的 COMPLETE、event hash、identity、原始 NPZ/gzip、metrics、support、receipt 全部保持原字节。新程序的 SHA 列在 `configs/runtime_source_identity.json`；它与旧 114-source 契约并列记录，绝不把修改后的文件哈希写回旧生产契约。

`validation.py` 从数据重建统计并核验完成标记；它不调用新积分。`visualization.py` 生成当前汇总图；`presentation/` 保留原流场渲染和点诊断方法。点流线中的按出口分层选线是展示规则，不能作为粒子来源选择或 CFD 分流比例。原视频按轨迹年龄对齐，不能当成同时注入的观测。

已知局部 P1 表示限制、历史时间步/网格研究的适用范围、未完成 H0 交付和 `SKIPPED_BY_USER` 项都保留；本次迁移和一致性检查没有给它们新增物理精度结论。
