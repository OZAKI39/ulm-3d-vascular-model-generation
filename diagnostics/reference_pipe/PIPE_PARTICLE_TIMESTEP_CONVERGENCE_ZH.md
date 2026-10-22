# Particle timestep convergence

2026-10-04。状态：**NOT_RUN_FIELD_GATE_FAILED**。GPU manufactured 稳态通过，但 field 科学门槛失败；按用户规则停止，时间步测试执行数为 **0**。

`0.5 / 0.25 / 0.125 ms` 的解析场与 CFD 场 canonical 对照均未运行。因此没有 `0.25 vs 0.125 ms` 收敛证据，也没有 `0.5 vs 0.25 ms` 支持或否定证据。不能写 DT_0P5MS_SUPPORTED、DT_0P5MS_NOT_SUPPORTED 或 TIMESTEP_NOT_CONVERGED；机器字段采用 null 表示尚无结论。

正式 BraVa 配置仍为 **dt=0.5 ms**。没有自动改为 0.25 ms，没有运行 0.0625 ms；本轮结果不构成修改正式步长的依据。固定灌注 1.240815660 s、随后观察 10 s、浓度 4×10⁹ m⁻³保持，未运行新正式 population。

阻断证据：[field 报告](MANUFACTURED_TRAJECTORY_FIELD_VALIDATION_ZH.md)。机器状态：[TIMESTEP_CONVERGENCE.json](TIMESTEP_CONVERGENCE.json)。
