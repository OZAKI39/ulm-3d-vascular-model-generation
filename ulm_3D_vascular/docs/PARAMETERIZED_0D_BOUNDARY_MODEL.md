# 可参数化几何水力网络：纯 0D 辨识与冻结

本实现是新的 `PARAMETERIZED_GEOMETRY_HYDRAULIC_0D_V1`，保留 H0 源码、报告和结果。参数辨识只使用既有几何、解析阻力、完整 0D 网络及合成/实验目标。它不读取 CFD 场或微泡/RBC 轨迹，不调用 3D 求解器，不使用 3D surrogate/admittance，不在未来 CFD 后反向调参。

## 数学模型与状态变量

局部阻力唯一权威实现为 `network_1d0d/hydraulic_resistance.py::linear_radius_resistance()`：

\[
R_e(\theta)=\frac{8\mu L_e}{\pi}
\frac{r_{0,e}(\theta)^2+r_{0,e}(\theta)r_{1,e}(\theta)+r_{1,e}(\theta)^2}
{3r_{0,e}(\theta)^3r_{1,e}(\theta)^3}.
\]

它是半径沿段长线性变化时对 $r^{-4}$ 的解析积分，等半径时连续退化为 Poiseuille 公式。旧 `one_d_flow.edge_resistance()` 的标量接口现在委托给它；没有新增第三套阻力公式。原旧函数在近等半径分支存在减法消差及等半径近似误差，统一实现消除了两套实现的漂移。

用每段 $g_e=1/R_e$ 组装整个网络的拉普拉斯导纳矩阵，执行

\[
G_{ff}(\theta)P_f=-G_{fb}(\theta)P_b,\qquad
Q_e=(P_u-P_v)/R_e(\theta).
\]

P、Q 是联立求解得到的状态，不是优化变量。本例是 7422 节点、7421 边的完整 A 网络，包含 ROI 之外保存的支路；不把三个出口各自当成独立的 `P=RQ` 黑箱。

源点仍是 H0 明确指定的结构节点 2410（模型假设，不是已证实的动脉标注）；其他结构叶节点为共同的表压参考终端。相对共同参考压力先令源压为 1 Pa、终端为 0 Pa，求出正确方向的 ROI 入口流量 $Q_{in}^{unit}>0$，再以

\[
\lambda=Q_{target}/Q_{in}^{unit},\quad
P=P_d+\lambda P^{unit},\quad Q=\lambda Q^{unit}
\]

匹配规定的 ROI 总流量。不会用 `abs(Qin)` 隐藏错误方向；不会把 $P_d$ 随 λ 放大。当前配置只支持**共同**远端参考压力，异质远端压力会被明确拒绝。

## 参数、单位与来源

| 参数/量 | 单位 | 定义与来源 |
|---|---|---|
| `s_O1`、`s_O2` | 无量纲 | 对现有 O1/O2 downstream component 节点半径的有效尺度；默认优化 `log(s)` |
| `terminal_resistance_O3_pa_s_m3` | Pa·s/m³ | O3 显式终端闭合参数；可选优化 `log(R/R_ref)` |
| `R_ref_pa_s_m3` | Pa·s/m³ | 对数坐标的单位参考，默认 1e17；不是实验观测或辨识约束 |
| `mu_pa_s` | Pa·s | 默认 0.00345312；可配置，但不是本版优化自由变量 |
| 密度 | kg/m³ | 默认 1056，记录在配置/冻结文件中，不进入当前稳态电阻网络方程 |
| `roi_target_flow_m3_s` | m³/s | 默认 1.551359160440232e-14，沿用用户指定 H0 工作点；不是新的 CFD 调参输入 |
| `distal_reference_pa` | Pa | 表压参考；不是测得的真实组织压力 |
| 长度、半径 | m | 保存的 SI 管网几何，求解内部不使用 μm/mm |
| 节点压力、边流量、阻力 | Pa、m³/s、Pa·s/m³ | 网络状态或由几何积分导出的量 |

`GeometryHydraulicCache` 预存拓扑、段长、原始端点半径、默认黏度阻力、ROI 边掩码、O1/O2 下游节点及边掩码、源/终端/端口索引和端口符号。数组复制到不可写快照；每次 evaluation 构造新的有效半径数组，再对所有几何边向量化调用权威阻力函数。

下游掩码通过移除 ROI 内部边后的无向连通分量确定，和 H0 `external_radius_variant()` 一致，不依赖 SWC 边的存储朝向。**real-cut 节点不缩放**，因此第一条外部边的一端固定、另一端缩放，自然形成 taper transition。不能对第一条边统一乘 `s^-4`。端口流量符号复用 `port_signs()`：INLET 向 ROI 内为正，O1/O2/O3 向外为正；输出顺序固定为 INLET/O1/O2/O3。

## O3 的两个模式

`LEGACY_REFERENCE`：现有 O3 节点直接固定在 $P_d$，严格复现 H0，不声称存在有限几何下游阻力。

`TERMINAL_RESISTANCE`：添加独立的虚拟 distal 节点及一个阻力边，移除 O3 自身的 Dirichlet 条件，只固定新 distal 节点。R 来源标为 `MODEL_CLOSURE_PARAMETER`，不标为 geometry-derived。虚拟项只进入 augmented-network metadata，原始 SWC ID 数组不伪造节点/边编号。求解后核对 $P_{O3}-P_d=R_{O3}Q_{O3}$；此式是检查，不替代全网求解。固定闭合阻力不会自动随黏度缩放。

## 目标、残差与可辨识性

配置字段显式区分 `value_dimensionless/sigma_dimensionless`、`value_m3_s/sigma_m3_s` 和 `value_pa/sigma_pa`。支持出口分流、平均出口 Q、real-cut 压力、命名端口压差及对数参数 prior。未知配置字段会报错，避免拼写或单位错误被忽略。

残差向量包含 `(model-target)/sigma` 和 `(x-x_prior)/sigma_prior`，交给 `scipy.optimize.least_squares(jac='3-point')`，不只返回标量目标。sigma 是标准差/残差尺度，权重是 `1/sigma²`。当前采用对角加权，不实现测量协方差；分流观测的相关性不能据此忽略。暂不输出参数置信区间。

三出口分流只有两个独立自由度；固定 Qin 时，出口 Q 和 fraction 也不能重复增加约束数。结构计数不足立即抛出 `UnidentifiableParameterizationError`，给出自由参数数、独立流量约束、压力约束和 prior 约束。固定 anchor 是**从自由参数中移除参数**，不是凭空增加数据。legacy O3 已固定的压力不作为辨识信息。

计数足够后，还实际用纯 0D 前向差分检查初始和最终残差 Jacobian 的秩。数值阈值为 `max(1e-8, largest_singular_value*1e-7)`；这只是以当前 sigma 缩放后的局部敏感性门槛，不是全局唯一性证明。prior 可以补足秩，但报告明确写 `PRIOR_REGULARIZED_NOT_DATA_IDENTIFIED`。数值秩不足会拒绝返回伪唯一结果。参数边界不当作独立观测。

`SYNTHETIC_TARGET` 由同一个纯 0D 模型生成或显式声明，不是实验真值。`EXPERIMENTAL_TARGET` 接受未来独立实验数据与来源说明；本轮没有真实实验校准。程序拒绝 CFD/particle target 类型，但不能自动鉴别人为误标成实验数据的外部数字；科研审阅仍需核实来源。

## 记录、冻结和未来 3D

每次 optimizer residual evaluation、初始/最终秩检查都写入 `parameter_trace.csv`；包含物理参数、P/Q、signed fractions、质量残差、所有标准化 residual 和范数。`nfev` 是 SciPy 的计数，另记 `total_forward_evaluations`，后者包含有限差分调用，不能混淆。失败评估会保留 FAILED 行并抛错，不用任意罚值或裁剪掩盖。

`fit_summary.json` 记录初末参数与残差、目标、prior、Jacobian 秩、优化状态、代码/配置/输入哈希。`network_state_si.npz` 保存所有节点压力、边流量和半径/阻力；虚拟节点如存在，通过单独元数据标明。没有反流时另给 positive fractions；一旦存在出口负流，保留 signed state 和反流审计，positive-only summary 为 unavailable，而不是把负流裁为零。

冻结命令重核配置、几何、模型代码、固定参数和收敛/秩门槛，只在**固定最终参数**上再执行一次 0D forward 校验，绝不再优化。输出 `FROZEN` YAML，独占创建、不覆盖，文件为只读权限并有 canonical SHA256 内容封印；读取时验证封印。具有文件系统所有权的人仍可修改权限，所以这里的 immutable 指“不可覆盖 + 可检测篡改”的工作流契约，不是密码学身份认证或不可篡改存储。

工作流固定为：`几何参数 → 解析 R → 全网 G P=b → 0D residual → 纯 0D optimizer → frozen θ → 未来独立 3D forward prediction`。**禁止从后续 3D 差异反向 retune，也禁止从微泡轨迹反推水力参数。**本轮 CFD 调用次数为 0。

optional FEM exporter 本版未实现。已阅读 `extension_transfer.py`，其 `measure_extensions()` 隐式使用默认黏度、`pressure_transfer()` 写死 H0 reference 描述；直接接入会错误表达新模型来源或可配置黏度。将来如单独实现，只能读取冻结参数及表面几何，以 `cap_pressures()` 和统一 `shift_pressure_gauge()` 分别输出 physical real-cut、raw cap 和 gauge-shifted cap pressure，不能读取 CFD 结果或重启拟合。本版不创建 FEM case。

## 文件复用与隔离

| 既有文件 | 本轮用法 |
|---|---|
| `hydraulic_resistance.py` | 唯一权威 local law；保持文件未改 |
| `utils/cfd_preprocess/one_d_flow.py` | 检查旧 `edge_resistance()`/`solve_global_flow()`；仅将阻力标量接口委托到权威实现，旧全网求解不改 |
| `network_solver.py` | 原样复用 `conductance_matrix` → `solve_network` 和正向 `operating_point_scale` |
| `idealized_h0.py` | 文件不改；复用 `H0Domain`、`h0_boundary`、`port_signs`；旧 `solve_operating_point`/`external_radius_variant` 仅作独立 regression oracle |
| `boundary_conditions.py` | 复用通过 `h0_boundary` 调用的显式 source gate；沿用表压语义；本轮不做 cap transfer |
| `extension_transfer.py` | 已审查，未在 optimizer 中导入/调用，未修改 |
| H0 两份中文报告 | 只读理解历史工作流；其中 3D 数字不进入拟合目标 |
| `roi_operating_point_H0.json` | 仅作 legacy 0D 回归金标准，不作为实验数据 |
| `downstream_resistance_audit.json` | 只读核对 O1/O2 既有下游和 O3 N/A 的语义，不作可独立优化的三个黑箱 R |

便携几何加载使用已保存的 `analysis_A_H0_graph_si.npz` 和 v1 `roi_ports_in_a.json`；只保留后者的名称、节点 ID、法向几何字段，不使用 FEM/pressure 元数据。原 `load_h0()` 仍可独立运行，但会检查含历史 FEM 文件的绝对路径清单，因此新优化器不调用它。

回归优先读取当前 H0 0D JSON；仅当该金标准文件缺失时，使用 `configs/parameterized_hydraulics/legacy_h0_snapshot.json` 的明确快照，报告标为 fallback。快照附有原文件 SHA256、来源提交和“非实验真值”说明；缺失几何不会触发该回退。

## 使用与复现

从仓库根目录执行，需 Python + NumPy + SciPy + PyYAML + pytest。项目 Python 环境可作为 `python` 使用。

```bash
python ulm_3D_vascular/scripts/fit_a_parameterized_0d.py \
  --config ulm_3D_vascular/configs/parameterized_hydraulics/legacy_h0_regression.yaml \
  --output /tmp/parameterized_0d_legacy_new

python ulm_3D_vascular/scripts/fit_a_parameterized_0d.py \
  --config ulm_3D_vascular/configs/parameterized_hydraulics/synthetic_radius_fit_example.yaml \
  --output /tmp/parameterized_0d_synthetic_new

python ulm_3D_vascular/scripts/freeze_a_parameterized_0d.py \
  --config ulm_3D_vascular/configs/parameterized_hydraulics/synthetic_radius_fit_example.yaml \
  --fit-summary /tmp/parameterized_0d_synthetic_new/fit_summary.json \
  --output /tmp/parameterized_0d_synthetic_new/frozen_parameterized_boundary.yaml

PYTHONPATH=ulm_3D_vascular python -m pytest \
  ulm_3D_vascular/tests/network_parameterized \
  ulm_3D_vascular/tests/network_1d0d \
  ulm_3D_vascular/tests/network_h0 -q
```

输出目录必须是新目录。已有示例数值、日志和冻结文件见 `reports/parameterized_0d_v1/`；不要覆盖这些审计记录。旧 H0 测试保留其原有绝对路径依赖，在其他机器需具备其历史输入；新增参数化测试仅用仓库内已有几何/报告，可独立运行。

## 科学限制

当前仅为稳态、不可压缩 Newtonian 圆管等效电阻网络；无 compliance、inertance、RCR、脉动传播或 Womersley impedance。SWC 几何、根/叶边界和圆管近似并不等于真实脑循环。s 表示有效半径修正，不自动解释为测得的管径；O3 R 是未解析下游的闭合假设。几何半径与黏度/闭合阻力等可能相关，尤其全局黏度无法仅由 fixed-Q 的 flow split 独立辨识。合成参数恢复证明软件自洽，不证明真实生理参数可唯一校准或 3D 预测准确。
