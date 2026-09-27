# SOURCE_CFD_STEADY_REFERENCE：当前源码及已存审计复核

本文件只读取现有输出、核验代数运算；没有运行 Musubi 或读取 rho/u 作为 restart。
当前 accepted iteration 为 598755，physical time 为 0.001220703363914373 s。
所有 source 文件副本和 SHA256 见 source_cfd_steady_reference.json。

## 准确数学定义

`utils/cfd_flow/steady_state.py::trapezoidal_mean` 对严格递增时间（固定 dt 下用 iteration）作梯形积分，除以窗口总时长。`audit_steady_window` 只使用 a=long-start、b=short-start、c=end 三个等间隔 checkpoint；short 用 b,c，long 用 a,b,c。

| 指标 | 当前 source 精确定义 | 迁移属性 |
|---|---|---|
| R_mass(W) | abs(mean_W(Qin) - sum mean_W(Qout_i)) / max(abs(mean_W(Qin)), smallest positive normal float64) | SOLVER_INDEPENDENT；是体积流量闭合，不是离散 LBM mass identity |
| physical_volume_closure | abs(Qin(c) - sum Qout_i(c))/abs(Qin(c))；`physical_port_flux.py::evaluate_physical_port_fluxes` | SOLVER_INDEPENDENT 代数定义；source 使用末时刻值，本轮主 gate 改为规定的窗口平均，另报告此 source 对照指标 |
| R_velocity | abs(mean_speed(c)-mean_speed(b))/max(abs(mean_speed(c)),1e-12 m/s)；mean_speed 是各 fluid node 的速度模平均 | SOLVER_INDEPENDENT；只是标量均速变化，不能代替全场 L2 残差 |
| R_pressure | max over inlet gauge pressure and three inlet-to-outlet drops of abs(P(c)-P(b))/max(abs(P(c)),1 Pa) | 代数 SOLVER_INDEPENDENT；inlet pressure 在 source 来自 adaptive controller，取得方式 SOLVER_SPECIFIC |
| R_inlet | abs(short_mean_Qin-Qtarget)/Qtarget | SOLVER_INDEPENDENT |
| flow_fraction_drift | 每个出口在 a,b,c 三个 checkpoint 的 instantaneous flow_fraction 最大值减最小值，再取出口最大 | SOLVER_INDEPENDENT；本轮另要求 short/long 及 confirmation 的窗口均值比例变化，不混为同一个量 |
| significant_backflow | 对 short、long 任一窗口，任一 mean_Qout<0 且 abs(mean_Qout)>0.05*abs(mean_Qin) | SOLVER_INDEPENDENT；本轮要求所有 outlet mean>0，比 source 的 5% 豁免更严格 |

## 本轮采用方式和差异

本轮保留 source 的可迁移标量均速、闭合及显著回流诊断作为伴随列，但所有主流量 gate 使用密集采样的物理时间窗口平均。窗口起点采用线性插值后梯形积分，不把旧 iteration 数移植。

按任务要求增加真正的中央控制体 dM/dt + outward mass flux；相对总质量漂移仅作为单独诊断。速度主 gate 使用逐节点全场 L2，压力主 gate 使用固定 gauge 压力场 L2，标为 HEMOCELL_SPECIFIC_CONVERGENCE_METRIC；不把弱得多的 source 标量指标冒称为全场收敛。

source inlet controller target、controlled link flux、full-timestep Musubi population identity、streaming slot 及 pressure_eq/wall_libb 实现均为 SOLVER_SPECIFIC，不移植为 HemoCell BC 或守恒实现。HemoCell 继续使用已冻结 Guo regularized off-lattice no-slip、VelocityPlugProfile3D 与 DensityNeumann。source 的 wall_libb 与 Guo 是不同数值实现。

阈值来自当前 `validated_contract.py` 和 `configs/cfd_flow.yaml`：mass 0.01、velocity 0.01、pressure 0.005、inlet 0.01、fraction drift 0.01。新增 CV 0.01、多面 spread 0.02、dx 对 dx/2 0.01 来自本轮用户要求，仅称 STEP3B_CROSS_SOLVER_COMPARISON_GATE。

## 复核结果

已从当前 checkpoint history 独立重算 short/long 四端口均值、R_mass、R_velocity、R_pressure、R_inlet、fraction drift 与末 checkpoint physical closure，全部与当前 final audit 一致（允许浮点末位舍入）。精确数值及逐项结果见 JSON。HemoCell 尚未自动收敛时，这些数字只作 reference，不进行 cross-solver agreement 判定。
