# 有限设计范围内最佳可实现三出口分流：0D 设计与唯一 3D 前向验证

## 结论

固定完整 A 模型与 `[0.5,2.0]²` 合成设计范围内，得到 `s_O1=2`、`s_O2=0.898280224172206`。0D 主目标 J 比 `(1,1)` 改善 51.78994731%，O1 上界激活。唯一新 CFD 已实际运行且通过独立日志、原稳态门槛和最终场检查；状态为 **BALANCE_IMPROVED_IN_3D**，3D J 相对实际旧 H0 改善 40.00537202%。没有 CFD 回调优化，没有第二次 CFD，没有 Particle/RBC 运行或生产晋升。

精确等分：**NO / NOT REQUIRED**。旧不可达性证书及 `BALANCED_0D_TARGET_NOT_REACHED` 原样保留。本轮结果是声明范围内经网格搜索与双初值检查一致的设计，不是生理最优、实验参数辨识或网格无关性证明。

## 1. 版本、方法与科学约束

- 实际 starting commit：`c5635dc0b0a2404f8eba9b347784d0429071f795`；历史基线 `f43c09b1e55cbd9702a65c600a464154fcf346f3` 未回退。
- 工作树：`/home/lzy/projects/temp_storage/github_sync_20260927`；分支 `sync/current-h0-dt1ms-wss-audit-20260927`。
- 完整图：7422 节点、7421 边，原始几何与端口映射哈希已封存。
- 所有优化 forward 调用原 `parameterized_hydraulics.solve_parameterized_operating_point()`；完整图解析阻力 → `G P=b` → ROI 提取。无 ROI-only surrogate、无独立端口 R/P 调参。
- `workflow_kind=DESIGN_OPTIMIZATION`，`parameter_interpretation=DESIGN_VARIABLES`，`data_identification_claim=NONE`。旧辨识保护未修改。
- `mu=0.00345312 Pa·s`，`Pd=0 Pa`，`O3=legacy_reference`，附加 O3 terminal resistance 为 `None`，ROI/J2→O1/J2→O3 几何固定。
- 0D ROI Qin `1.551359160440232e-14 m³/s`；production FEM Qin `1.551359160885543e-14 m³/s` 原样保留，不把网络流量的末位浮点差异写回入口。
- bounds 来源 `SYNTHETIC_DESIGN_ENVELOPE`：两项 `[0.5,2.0]`，仅为当前数值设计范围，**不是生理/测量范围或置信区间**。
- 旧证书状态 `EQUAL_SPLIT_IMPOSSIBLE_UNDER_PRESCRIBED_MODEL`，哈希核验及对应 pytest 通过。必要约束 `Q_O1/Q_O3 <= 0.35320878510001563` 排除了 1/3 等分。更大必要区域投影 `[0.1616970483999667, 0.38050842546960284, 0.4577945261304305]` 不代表有限参数可实现解。

## 2. 纯 0D 优化与基线

目标 `J=sum((f_i-1/3)^2)`；守恒时等于 `3*population_variance(f)`。优化器 `scipy.optimize.least_squares`，log 坐标、原生有限 bounds、三点 Jacobian，`ftol=xtol=gtol=1e-12`，没有参数 clipping。NumPy 2.5.3；SciPy 1.18.1。

| 工况 | O1 % | O2 % | O3 % | J | 极差 | 总体 std | 最大等分偏差 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0D baseline | 6.30232927 | 49.59000950 | 44.10766123 | 0.111104084196 | 0.432876802261 | 0.192444004146 | 0.270310040589 |
| 0D primary design | 14.96125190 | 38.68966709 | 46.34908101 | 0.0535633375268 | 0.313878291133 | 0.133620529269 | 0.183720814373 |
| 旧 H0 3D（实际 VTU 积分） | 6.54082456 | 44.82012154 | 48.63905390 | 0.108404991171 | 0.420982293383 | 0.190092075208 | 0.267925087731 |
| 新 3D（唯一 forward） | 15.51494709 | 32.91168032 | 51.57337259 | 0.0650371711635 | 0.360584254993 | 0.147238096478 | 0.18240039255 |

主初值 `(1,1)`；优化器 `ftol termination condition is satisfied.`，nfev=9、njev=7、optimality=2.76994e-12。实际主优化 forward 评价 38 次，primary 全流程含网格/第二初值/重复共 726 次。`nfev` 不含全部有限差分调用，因此不同于 forward 总次数。

`BOUND_ACTIVE_OPTIMUM`：`s_O1: UPPER`，`s_O2: NONE`。在有限范围内落在上界是允许的结果，结果依赖设计范围。J 改善 0.0575407466690238，极差降低 0.118998511127909，std 降低 0.0588234748769973。

全网质量审计：最大相对残差 1.7460901e-13，ROI 相对残差 2.4611237e-14；无端口倒流，所有有效半径、边阻力有限且为正。

25×25 对数网格共 625 点，网格最好 J=0.0536588288800794；从其启动第二次 bounded 局部优化。两个局部结果参数绝对差 `[0.0, 2.271137500287068e-09]`，J 差 1.179612e-15，交叉检查 **CONSISTENT**。重复主优化的 38 次评价载荷逐项完全一致。此检查不等同于全局最优数学证明。

### 纯 0D 上界敏感性（仅表格）

下界始终 0.5。1.25、1.5、3 的结果仅为 0D 诊断；唯一 CFD 始终使用预先声明的上界 2，未根据结果更换 envelope。

| s_max | 最佳 s_O1 | 最佳 s_O2 | O1/O2/O3 % | J | 极差 | active bounds | 用途 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.25 | 1.25 | 0.9206458195 | 10.581691/41.326710/48.091599 | 0.0799337697936 | 0.375099076074 | UPPER/NONE | 0D_DIAGNOSTIC_ONLY |
| 1.5 | 1.5 | 0.9082140165 | 12.898728/39.868145/47.233127 | 0.0653481114553 | 0.343343987442 | UPPER/NONE | 0D_DIAGNOSTIC_ONLY |
| 2.0 | 2 | 0.8982802242 | 14.961252/38.689667/46.349081 | 0.0535633375268 | 0.313878291133 | UPPER/NONE | PRIMARY_DESIGN |
| 3.0 | 3 | 0.8940897079 | 15.903373/38.188369/45.908257 | 0.0485503600358 | 0.300048842153 | UPPER/NONE | 0D_DIAGNOSTIC_ONLY |

冻结文件：`/home/lzy/projects/temp_storage/github_sync_20260927/ulm_3D_vascular/reports/parameterized_0d_v1/best_feasible_balance/frozen_balance_design.yaml`。状态 `FROZEN_DESIGN`，SHA256 `e2ebf0127b11f9015f2261b973a6161b20e647bb588b1f2c339a2d671633f4a7`。设计在 FEM handoff 与唯一 CFD 启动前已封印；后续没有改写。

## 3. 单向 FEM handoff

从冻结设计取得 real-cut P/Q；只读取当前延伸段几何，用显式 `mu=0.00345312 Pa·s` 和现有渐变半径解析阻力积分。沿每段使用原 H0 的 40 个截面中点取样及端部半格常半径封闭规则。

```text
R_ext = (8*mu/pi) integral(ds/r(s)^4)
P_cap_raw[i] = P_realcut[i] - R_ext[i]*Q[i]
P_cap_applied[i] = P_cap_raw[i] - min(P_cap_raw)
```

| 端口 | 0D Q (m³/s) | real-cut P (Pa) | R_ext (Pa·s/m³) | raw cap P (Pa) | applied cap P (Pa) |
| --- | --- | --- | --- | --- | --- |
| O1 | 2.32102751805388e-15 | 52.019814677518 | 5.38266722060601e+16 | -72.9133727180135 | 237.679049278879 |
| O2 | 6.00215694616186e-15 | 3671.11749297271 | 5.57850462591372e+16 | 3336.28689007647 | 3646.87931207336 |
| O3 | 7.1904071401862e-15 | 0 | 4.31953874017835e+16 | -310.592421996893 | 0 |

共同减去 `-310.592421996893 Pa`（即加上 `310.592421996893 Pa`）；平移前后两两压差最大误差 `0 Pa`。raw cap 负值是参考零点下的表压，不裁剪。对每个出口使用相同常数平移，保持驱动压差。

纯解析 H0 handoff 回归已运行：R_ext、raw cap P、applied cap P 的最大绝对差均为 **0**；没有为此运行 CFD。详见 `handoff_H0_regression.json`。

## 4. 唯一 CFD：配置、GPU 与运行证据

本地算例：`/home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-best-feasible-balance-v1`。

服务器：`root@50.115.148.16:4159`；算例 `/workspace/flow_best_feasible_balance_20260928/mean-2p0-mmps-A-best-feasible-balance-v1`；private supervisor `best_balance_gpu`，autorestart=false。`reports/one_scientific_case_dispatch.json` 排他登记启动，完成记录为 `reports/dispatch_completion.json`。

**NEW SCIENTIFIC CFD RUN COUNT = 1**。其他 envelope 的 CFD、校准重试 CFD、粒子与 RBC 调用：均 0；**CFD-based retuning = NO**。

预检通过后才启动：H0 XML 仅许可 3 个出口 traction `Value` 字段变化；实际 O1 从 585.2864332598318 到 237.679049278879 Pa，O2 从 2932.0152710782013 到 3646.87931207336 Pa，O3 两者均 0。其余 XML 字段、体网格/外表面/壁面、入口、policy、PETSc options 均逐文件/逐字段核验不变。

求解器为原 `svmultiphysics`，P1/P1+VMS，刚性不可压缩牛顿流体，`rho=1056 kg/m³`，无滑移壁面；网格 70363 节点、371402 四面体。RTX 4090，MPI=1、OMP=1；沿用原 PETSc CUDA 矩阵和原预条件配置。没有为本任务改变并行数、后端或性能配置。

实际 dt 为 `1.5040600916882215e-07 s`，取自固定原算例 policy（分支名中的 dt1ms **不代表本算例步长**）。每 10 步采样，连续 5 个间隔满足速度变化 ≤1e-05、流量变化 ≤1e-06、全局质量残差 ≤1e-06，线性/非线性原门槛不变。从零场开始，无重启或缩放旧场。

实际 exit code `0`，wall time `1644.023508 s`（27.400 min），停止器合格 step `70`，最终原生输出 step `71`，最终物理时间 `1.0678826651e-05 s`。两者可能因写 STOP_SIM 的安全退出多一个时间步而不同。

GPU 使用证据：solver.log 的矩阵类型 `seqaijcusparse`；5 秒采样最大 GPU 利用率 98%，显存峰值 3840 MiB。实际 PETSc profile 摘录：

```text
Using PETSc Release Version 3.25.5, Aug 30, 2026
   GPU Mflop/s: 1e-6 * (sum of flop on GPU over all processes)/(max GPU time over all processes)
   GPU %F: percent flops on GPU in this event
MatMult           168420 1.0 4.2433e+01 1.0 5.34e+12 1.0 0.0e+00 0.0e+00 0.0e+00  3 14  0  0  0   3 14  0  0  0 125843   148506    159 2.04e+04    0 0.00e+00 100
MatLUFactorNum        15 1.0 8.7051e+01 1.0 2.01e+11 1.0 0.0e+00 0.0e+00 0.0e+00  5  1  0  0  0   7  1  0  0  0  2313       0      0 0.00e+00    0 0.00e+00  0
KSPSolve             159 1.0 1.1697e+03 1.0 3.92e+13 1.0 0.0e+00 0.0e+00 0.0e+00 71 99  0  0  0  92 99  0  0  0 33488   34804   167983 2.12e+04 82998 3.91e+02 100
PCSetUpOnBlocks      159 1.0 9.3832e+01 1.0 2.01e+11 1.0 0.0e+00 0.0e+00 0.0e+00  6  1  0  0  0   7  1  0  0  0  2146       0      0 0.00e+00    0 0.00e+00  0
PCApply           167506 1.0 1.0556e+03 1.0 2.67e+13 1.0 0.0e+00 0.0e+00 0.0e+00 64 68  0  0  0  83 68  0  0  0 25314   25519      0 0.00e+00    0 0.00e+00 100
Using PETSc directory: /workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3n/external/petsc325/install_gpu13
Using PETSc arch:
```

这是 CPU/GPU 混合实现：FEM 装配及 ILU 因子构建仍有 CPU 工作；profile 中 MatMult、PCApply 的 GPU 浮点运算占比为 100%，MatLUFactorNum 和 PCSetUpOnBlocks 为 0%。后端保持原 GMRES(100)、右预条件 ASM overlap=2、子域 ILU(2)、CUDA matrix/vector 配置。不能把显存占用或峰值利用率解读为整个 CFD 都在 GPU 上执行。完整执行 command、host、软件/库哈希及采样保存在算例 `reports/` 和 `run/petsc_profile.txt`。

## 5. 独立验收与实际 3D 分流

本地重新解析原生 solver.log，与远端历史逐项一致，线性/非线性门槛通过。重新读取所有保存的稳态快照，独立重算间隔，合格 step `[70]`。最终网格坐标/拓扑与输入一致，checkpoint step/time 验证通过。

- 实际 Qin `1.551359160440232e-14 m³/s`；入口相对目标误差 `2.8704574e-10`。
- 全局质量残差 `2.0339865e-16`，补偿求和残差 `1.0169933e-16`；各出口积分流量均为正，无净倒流，未据此声称每个面片都不存在局部回流。
- 最大速度 `0.00706521457226 m/s`，压力范围 `[-1.9701796326175063, 5579.380779850303]` Pa；全部有限。
- 壁面最大速度 `0 m/s`，无滑移通过。
- 最终场相对停止合格场的速度 L2 变化 `4.4623213e-15`、归一化端口流量变化 `2.0339865e-16`，均通过原门槛。

| 端口 | 实际 Qout (m³/s) | 实际 Qout/Qin % | 3D−0D 百分点 | cap 面积平均流体 P (Pa) |
| --- | --- | --- | --- | --- |
| O1 | 2.4069255290344207e-15 | 15.51494709 | +0.55369519 | 241.22806695982177 |
| O2 | 5.105783675399739e-15 | 32.91168032 | -5.77798677 | 3655.1780850602668 |
| O3 | 8.000882399968158e-15 | 51.57337259 | +5.22429158 | 9.478022864899993 |

分流采用带符号 Qout/实际 Qin，没有 clipping 或强行归一化。表中 cap 面积平均流体压力由解场积分，区别于施加的法向牵引压力；两者不必逐点相同。3D–0D 差异如实保留，既不是调参信号，也不作为必须为零的验收门槛。

实际旧 H0 3D 分流由其正式 frozen VTU 重新积分，未硬编码四舍五入数字。比较状态 **BALANCE_IMPROVED_IN_3D**：仅在新 J 和极差均低于旧值时才标记改善。J 相对降低 40.00537202%，极差减少 0.0603980383906，std 减少 0.04285397873（正值为改善，负值为恶化）。旧/新完整指标见第 2 节和 `all_balance_metrics.csv`。

新 frozen flow：`/home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-best-feasible-balance-v1/frozen_flow/steady_flow_mean_2p0_mmps_A_best_feasible_balance.vtu`。其 `manifest.json` 包含设计 envelope、0D 参数/指标、cap pressures、实际 3D 指标、0D–3D 差异、旧 H0 比较、输入/输出哈希及最终 step；`particle_production_promoted=false`。

## 6. 已有 WSS 链路复用

复用并归档已验证的 `flow_solver_support/wss.py` 与 `wss_case.py` 原字节源码，哈希见 `evidence/wss_core/provenance.json`，旧解析张量及回归证据同时保留。没有重跑圆管 CFD。本轮对旧/新实际 FEM 场重新运行该链路：从壁面相邻四面体的 P1 完整速度梯度张量，计算单位外法向的切向黏性牵引模长：

```text
tau = mu * (grad(u) + grad(u).T) @ n
WSS = norm(tau - dot(tau,n)*n)
```

只处理壁面 tag 1，排除入口/出口截面；原坐标 m、速度 m/s、黏度 Pa·s 得到 Pa。统计使用原始面片量，面积加权均值与经验累积面积分位数。J1/J2 区域按面片中心落入半径 5 μm 的球选择，球心分别 `(92,49,111)` 和 `(130.04,82.04,87.18)` μm。全壁面统计含人工延伸段，不能称其为仅真实血管统计。

| 工况 | 区域 | 面片数 | 面积 μm² | 面积均值 Pa | 面积 P5/P50/P95 Pa | 最小/最大 Pa |
| --- | --- | --- | --- | --- | --- | --- |
| old_H0_3D | whole_wall_including_extensions | 45221 | 2024.19 | 12.531 | 2.03148/14.7158/26.1013 | 0.414767/46.5202 |
| old_H0_3D | J1_ball_5um | 2211 | 98.181 | 23.3055 | 15.3199/20.4385/35.7881 | 2.76347/46.5202 |
| old_H0_3D | J2_ball_5um | 2136 | 95.8368 | 10.6897 | 2.08486/6.30428/22.3984 | 0.414767/30.2411 |
| best_feasible_3D | whole_wall_including_extensions | 45221 | 2024.19 | 14.1133 | 2.21022/14.6146/30.8396 | 0.811451/46.2614 |
| best_feasible_3D | J1_ball_5um | 2211 | 98.181 | 23.4838 | 11.6244/24.7151/35.9028 | 1.32113/46.2614 |
| best_feasible_3D | J2_ball_5um | 2136 | 95.8368 | 13.5769 | 4.79507/8.10362/26.6284 | 0.811451/36.8608 |

`figures/wss_raw_comparison.png/pdf` 使用相同正交视角、相同网格、原始面片值与线性色标 0–55 Pa，无平滑。VTP 另外保留原链路的面积加权节点显示值，但本图和 CSV 不使用它。0D 网格目标与真实 3D 分流比较见 `figures/balance_design_and_forward.png/pdf`；bound sensitivity 仅表格。

## 7. 测试、执行与复核

启动前实际执行命令（从仓库根目录）：

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=ulm_3D_vascular \
/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B -m pytest \
ulm_3D_vascular/tests/network_parameterized \
ulm_3D_vascular/tests/network_1d0d \
ulm_3D_vascular/tests/network_h0 -q -p no:cacheprovider
```

结果 **138 passed in 10.23s**。覆盖旧完整测试与新增的不可达性保留、辨识 guard 不变、原生 bounds、无限半径拒绝、设计改善、方差关系、trace 确定性、双初值一致、bound-active、冻结 roundtrip/篡改/覆盖拒绝、零 CFD 依赖、H0 handoff 回归、XML 受控改变及唯一启动保护。早期一次新增测试发现 NumPy scalar YAML 序列化问题，已在正式冻结前修复；失败日志保留，不冒充首轮全部通过。

完整输入和下述只读后处理入口随提交保存；独立验收的 frozen 输出禁止覆盖，复核时应在独立副本中执行。**不要重启 scientific CFD**；唯一 dispatch 文件存在时启动器会拒绝重复执行。

```text
scripts/optimize_balanced_0d_design.py --config configs/parameterized_hydraulics/best_feasible_balance_design.yaml --output <new-audit-dir> --freeze
scripts/prepare_best_balance_fem.py --help
scripts/validate_best_balance_fem.py --case <case-copy> --source-case <old-H0> --flow-root <FEM-root> --report <report-copy>
scripts/postprocess_best_balance.py --case <case> --source-case <old-H0> --flow-root <FEM-root> --design <design-dir> --report <report-dir>
scripts/report_best_balance_forward.py --case <case> --source-case <old-H0> --design <design-dir> --report <report-dir>
```

数值证据：`design_trace.csv`、`coarse_grid.csv`、`bound_sensitivity.csv`、`all_balance_metrics.csv`、`final_3D_ports.csv`、`wss_comparison.csv`；独立间隔检查见 `independent_validation.json`。原生中间快照留在本地和服务器 case 的 `run/1-procs/`；Git 保存正式最终场、checkpoint、完整求解日志和验收摘要，不重复提交所有中间大文件。

## 8. 关键哈希

| 文件/对象 | SHA256 |
| --- | --- |
| frozen_balance_design.yaml | `e2ebf0127b11f9015f2261b973a6161b20e647bb588b1f2c339a2d671633f4a7` |
| new frozen VTU | `fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2` |
| old H0 frozen VTU | `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4` |
| volume mesh | `1a7a69dcba52f0475b2280775921e4dedcdfefa67b7965096d86cb9dba38bdf9` |
| solver binary | `0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7` |
| PETSc library | `b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef` |
| case solver XML | `afe21f47bc4b53990c8f32a3d03044b36dc9a78038e79429eb4966ec2a302388` |
| solver.log | `313c9632a3d2b1d570208871d67b8e0c8fa414d66e55e2571a4e0a929cd7a0be` |
| parameter_fit.py unchanged | `4bad7b5760c6976acfa5bb9550c477101dc4a9c59e0d38e5e55011cbe5b56610` |
| feasibility certificate unchanged | `2aece51d379d392533ddf2b43e1ec393590c6c57d3347eee1478904d016e209c` |

## 9. 尚未解决的科学限制

1. 这是给定完整 A、O3 参考闭合和合成 bounds 下的设计，不是生理真值或实验 calibration；O1 上界激活表明结果依赖 envelope。
2. 25×25 网格和双初值交叉检查增加数值可信度，不构成严格全局最优证明；已有不可达性证书的必要区域投影也不是已实现的有限设计。
3. 0D 圆截面/充分发展阻力近似、延伸段阻力补偿和 3D 分叉几何不同，因此允许分流偏差；本轮没有耦合迭代或 CFD 反馈校正。
4. 保持原单张网格、P1/P1+VMS 和相同时间策略，只完成受控边界改变的 forward 比较；没有新的网格独立性、内部截面质量误差研究、时间步或模型不确定性量化。全局质量守恒不能证明局部梯度/WSS 准确。
5. WSS 面片差异受 P1 梯度、曲面离散和原网格分辨率影响；本轮不据图像平滑程度诊断物理正确性，也不为颜色均匀调压或平滑。
6. 本轮未与实验流量/压力/WSS 比较，未运行 CORE500、microbubble 或 RBC，未把新流场晋升为轨迹生产输入。
