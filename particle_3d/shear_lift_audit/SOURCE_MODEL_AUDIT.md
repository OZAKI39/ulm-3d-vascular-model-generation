# 正式动力学源码核对

核对工作树 HEAD：`775adc019536585a4ee2f4dbf2c083c959e0a339`。具体工作树内容由 `data/readonly_baseline.json` 的逐文件 SHA256 固定，不能仅用提交号代表有未提交文件的工作树。

正式保存轨迹入口是 `particle81_simulation.SavedTrajectoryStepper`，P82A 和 PPT 数据也保存同样列结构和 `UNCHANGED_PARTICLE65STEPPER_NO_LAMMPS_NO_MB_MB_OR_RBC_COUPLING` 元数据。调用链：

`FrozenFEMField.sample` → `u, 0.5 curl(u)` → `Particle65Stepper.step_to` → `v1_trial` → `assemble_v1` → `solve_resistance` → `refine_interval` 接受子步 → `SavedTrajectoryStepper._accepted`。

`hydrodynamic_resistance.sphere_self_diagonal` 的平动、转动块分别为 6πμa 和 8πμa³。`particle65_motion.assemble_v1` 只把各球 self diagonal 乘以背景自由广义速度作为右端项 `diag*u`；矩阵为 self block 加 `zeta*row.T@row`。法向行只作用于平动分量。

`nearfield_regularization.NearFieldRegularizationV1.evaluate` 复用 P5 球壁 leading 6πμa²/h 系数，乘既有 C1 权重，分母使用既有 handoff 下限。`nearfield_handoff.handoff_constraints` 保留原壁面几何，构建法向不穿透速度约束。`resistance_solver.solve_resistance` 求准静态速度并施加约束；没有质量乘加速度的独立惯性力积分。

核对结论：上述单 MB 调用路径中没有 Saffman、wall-induced inertial lift、deformation lift、Magnus lift、Brownian、重力/浮力或声辐射力的右端项。背景有剪切和球体随涡量旋转，并不表示已有剪切升力或 Magnus 力。P6.5 的既有报告 `particle_3d/reports/particle6_5/PARTICLE6_5_REVIEW.md` 也明确没有加入惯性、摩擦、表面力或壳层力学。

## 保存时序

`SavedTrajectoryStepper._accepted` 先接受新位置，再把本步保持的速度和 `evaluated_at_time_s=old.time_s` 一起保存。第 k 行位置是 x_k，速度是前一区间 V_[k−1,k]。因此：

- 每行完整保留 endpoint 位置、速度、endpoint 流速以及 `endpoint_slip_speed_m_s`，满足原始保存状态诊断。
- 主要力比较使用 x_(k−1) 的流速、梯度、法向与间隙，配合同一行 V_[k−1,k]。验证时间、dt 及 x_k=x_(k−1)+dt V 三个恒等式。
- 第一行是初始自由流诊断，明确标记 `initial_diagnostic`，dt=0；不是已完成一次阻力平衡的证据。
- 曲壁条件下用最接近的原有限三角形计算球壁法向；对齐间隙与正式保存值逐行交叉检查。没有用光滑替代壁面。

本阶段没有编辑以上源码。新增代码在本目录 `code/`，正式模块不导入该目录。永久测试验证球体自由流右端项、对齐关系、输入哈希及正式源码没有导入后处理审计。远程程序只读已保存数组，不调用 `step_to`、`integrate_one`、`v1_trial` 或轨迹姿态推进。
