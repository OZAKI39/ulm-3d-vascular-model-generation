# Frozen SimVascular FEM Baseline — Stage SV1.3Q

状态：**ILU_REBUILD_POLICY_WINNER_FOUND**。FEM DEVELOPMENT: FROZEN。下表由本分支 frozen artifact fresh-read 和 accepted Stage Q JSON 自动生成，生成器是 `scripts/fem_freeze_sync/write_handoff.py`。

| 项目 | 实际值 |
|---|---:|
| rho (kg/m³) | 1056.0 |
| mu (Pa s) | 0.00345312 |
| Q target (m³/s) | 2.7369132390905703e-15 |
| dt (s) | 1.5040600916882215e-07 |
| 首次正式稳态 step | 70 |
| safe-stop step | 71 |
| 最终 TimeValue (s) | 1.0678826650986389e-05 |
| wall time (s) | 1571.0034701200202 |
| 最大速度 (m/s) | 0.0012522159476823856 |
| 最低压力 (Pa) | -0.22374142984762152 |
| 最高压力 (Pa) | 434.33471073169557 |
| Qin (m³/s) | 2.736913239090574e-15 |
| OUTLET_01 (m³/s) | 1.1653417879825574e-16 |
| OUTLET_02 (m³/s) | 2.3319920670259264e-15 |
| OUTLET_03 (m³/s) | 2.883869932663918e-16 |
| mass error |Qout−Qin|/Qtarget | 0.0 |
| inlet error | 1.2970356614863377e-15 |
| 最大壁面速度 (m/s) | 0.0 |
| KSP attempts | 161 |
| KSP iterations | 78874 |
| ILU numerical rebuilds | 15 |
| ILU reuse | 146 |
| linear unrecovered failures | 0 |
| nonlinear failures | 0 |

正式 SimVascular 网格有 70363 节点、371402 四面体和 45704 外表面三角形，单位 m。不是旧 DOLFINx 网格。网格 SHA256：`1a7a69dcba52f0475b2280775921e4dedcdfefa67b7965096d86cb9dba38bdf9`。

BC：INLET 为 steady Dirichlet flat profile、impose_flux、zero_out_perimeter；原 XML 入口值 `-2.736913239876193e-15`；三个出口为零 Neumann，WALL 为零 Dirichlet。ρ、μ、Q、dt 和 generalized-alpha spectral radius 0.5 保持实际运行值。

稳态门槛：每 10 步保存，连续 5 区间 E_u≤1e-05、E_Q≤1e-06，mass gate≤1e-06。step 70 触发停止请求；step 71 是允许完成的 in-flight timestep，final VTU 与 native checkpoint 时间相同。这里的 mass error=0 是当前浮点积分输出，不是物理误差恰好为零的断言。

背景流：[frozen_reference/flow/steady_flow_stage_sv1_3q.vtu](frozen_reference/flow/steady_flow_stage_sv1_3q.vtu)，8997629 bytes，SHA256 `373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1`。`Velocity`/`Pressure` 是 Float64 节点数据，均 finite，wall no-slip PASS；速度梯度、涡量和应变率没有存储。

Mesh convergence: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

Time-step sensitivity: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

CPU/GPU field-level equivalence: **DEFERRED BY USER DECISION**。

这些是 USER-ACCEPTED PROJECT DECISIONS，不是通过的验证。FEM performance tuning: frozen after Stage Q。Stage Q is the frozen particle-development background-flow baseline；不能解释为网格收敛、时间步收敛、CPU/GPU 科学等价或新的 FEM production approval。

[接口入口](PARTICLE_HANDOFF.md) · [Stage Q 原报告](reports/sv1_3q/REPORT.md) · [fresh-read 结构化摘要](frozen_reference/baseline_summary.json) · [原始验收](reports/sv1_3q/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json)
