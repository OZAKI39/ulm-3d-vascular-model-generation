# 有限设计范围内的三出口分流优化

本工作流为 `DESIGN_OPTIMIZATION`，参数为 `DESIGN_VARIABLES`，没有实验数据辨识主张（`data_identification_claim: NONE`）。它与原有 calibration 工作流并列，未修改其验收标准。

## 为什么改变研究问题

完整 A 管网、固定 ROI、被动正阻力、O3 为 `legacy_reference` 且远端共同参考压力为零时，已有不可达性证书给出必要条件 `Q_O1/Q_O3 <= 0.35320878510001563`。等分要求此比值为 1，因此当前模型不能实现三出口精确等分。原 `BALANCED_0D_TARGET_NOT_REACHED`、拟合拒绝记录和证书保留；本轮新成功状态为 `BEST_FEASIBLE_BALANCE_DESIGN_PASS`，二者不矛盾。

参数辨识询问数据是否唯一约束参数，设计优化则在声明的空间中选择目标函数最小的参数。`parameter_fit.py` 的 `UnidentifiableParameterizationError`、Jacobian 秩检验及其他辨识保护保持原样。没有实验观测时，不得称本轮参数被数据唯一识别；设计冻结也不要求数据 Jacobian 满秩。

## 正式设计问题

- 仅优化 `x=(log(s_O1),log(s_O2))`，两项原生 bounds 都是 `[log(0.5),log(2)]`。
- 来源为 `SYNTHETIC_DESIGN_ENVELOPE`：人为声明的数值设计范围，不是生理范围、实测半径范围或置信区间。有限上界从算法上禁止半径向无穷大逃逸。
- 每次评价调用原 `solve_parameterized_operating_point()`：完整 A 几何 → 解析边阻力 → 全图 `G P=b` → ROI 端口。没有 ROI 代理模型或独立端口阻力调整。
- 固定 `mu=0.00345312 Pa·s`、`Pd=0 Pa`、`O3=legacy_reference`、O3 终端附加阻力 `None`、网络 ROI 流量 `1.551359160440232e-14 m³/s`。ROI 及 J2–O1/O3 几何不变。
- 主目标 `J=sum((f_i-1/3)^2)`，在守恒条件下等于三倍分流比例的总体方差。另报极差、总体标准差、最大等分偏差。`1/3` 是目标中心，精确达到不作为验收条件。

采用确定性的 `scipy.optimize.least_squares`，在 log 坐标使用原生 bounds、三点数值 Jacobian 和 `ftol=xtol=gtol=1e-12`；不通过 clipping 伪造约束。从 `(1,1)` 出发，辅以 25×25 对数网格，网格最佳点再作为第二初值；重复首个优化检查完整评价 trace 一致。网格和两个初值一致仅为交叉检查，不是全局最优的数学证明。

最优点落在上界标为 `BOUND_ACTIVE_OPTIMUM`，不是失败。它表示结果依赖声明的设计范围，不代表生理真实最优。上界 1.25、1.5、2、3 的纯 0D 敏感性仅用于表格展示；进入唯一 CFD 的始终是预先声明的上界 2。更大必要可行域的等方差投影 `[0.1616970483999667,0.38050842546960284,0.4577945261304305]` 不是有限范围内已实现的设计。

## 冻结、单向传递与唯一 3D 预测

冻结前检查原生优化成功、有限正半径/阻力、有限且在 bounds 内的参数、全管网质量守恒、无出口倒流、目标和极差同时改善、网格/局部优化一致、确定性重复一致。冻结文件为 `reports/parameterized_0d_v1/best_feasible_balance/frozen_balance_design.yaml`，使用内容封印与输入/代码 SHA256；禁止覆盖。

`parameterized_fem_handoff.py` 只读冻结设计和当前 FEM 延伸段几何，不读 CFD 或轨迹。延伸段按现有解析渐变半径阻力积分法，显式使用上述黏度：

```
R_ext = (8*mu/pi) * integral(ds/r(s)^4)
P_cap_raw[i] = P_realcut[i] - R_ext[i]*Q[i]
P_cap_applied[i] = P_cap_raw[i] - min(P_cap_raw)
```

所有出口共同平移参考零点保持两两压差，原始负表压无需裁剪。出口牵引压力不是每个边界节点的流体压力 Dirichlet 条件。`s_O1=s_O2=1` 的解析 handoff 回归必须先复现旧 H0。

`prepare_best_balance_fem.py` 只注册一个指定算例。原 H0 XML 仅允许三个出口压力 Value 不同，其余材料、网格、壁面、入口、P1/P1+VMS、PETSc、时间策略与停止门槛保持原值；入口仍是 production `1.551359160885543e-14 m³/s`。`launch_best_balance_once.py` 校验 preflight 与所有输入，通过排他创建 dispatch 文件阻止再次启动，调用字节不变的 H0 运行器一次；禁止自动重启。

CFD 只作后续 forward prediction，结果不反馈给 0D。独立重读原生日志、每个保存快照、最终场和 checkpoint，按原门槛验收后才导出 frozen flow。分流按带符号的 `Qout/Qin` 计算，不裁剪或强制归一化。与实际旧 H0 VTU 重积分所得指标比较；即便 3D 未改善也不 retune，不进行第二次 CFD。

## 入口、结果与解释范围

纯 0D 入口：`scripts/optimize_balanced_0d_design.py`；准备：`scripts/prepare_best_balance_fem.py`；唯一启动：`scripts/launch_best_balance_once.py`；独立验收：`scripts/validate_best_balance_fem.py`。完整执行命令、数值及文件哈希见 `reports/best_feasible_balance_forward_v1/BEST_FEASIBLE_BALANCE_FORWARD_REPORT_ZH.md`。

本轮不修改当前 Particle/RBC 生产流场，不运行轨迹；冻结流场标为 `particle_production_promoted: false`。0D 的圆管近似、刚性牛顿模型、外部有效半径参数化、单向延伸段传递、固定单张 P1/P1 网格都会限制解释。求解收敛、质量守恒或 WSS 图平滑均不证明网格无关、生理真实性或实验参数可辨识性。
