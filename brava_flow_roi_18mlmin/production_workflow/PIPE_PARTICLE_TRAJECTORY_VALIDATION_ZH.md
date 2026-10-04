# Production particle trajectory 验证

2026-10-04。状态：**NOT_RUN_FIELD_GATE_FAILED**。本轮 GPU manufactured CFD 已达稳态，但 `REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC` 阻断后续阶段。

0.377 µm 带 u 相对 L2=8.519477%（要求≤5%），三个 bulk u 组也超过 5%。五个 particle bands 的 G/s/curl 均通过，仍不满足完整 field gate。详见 [field 报告](MANUFACTURED_TRAJECTORY_FIELD_VALIDATION_ZH.md) 和 [FIELD_VALIDATION.json](FIELD_VALIDATION.json)。

| 验证项目 | 实际执行 | 结论 |
|---|---:|---|
| Analytic production-particle oracle | 0 cases | 未运行，不能判 PASS/FAIL |
| 四个 sentinel / 完整十个独立场景 | 0 cases | 未运行 |
| CFD-vs-oracle 相同物理时间比较 | 0 cases | 未运行，无一致性结论 |
| 三个代表场景 CPU reference parity | 0 cases | 未运行 |

停止依据是用户本轮第 21 条“科学 field gate FAIL 立即停止”。因此没有用失败场计算轨迹，也没有用旧 CPU step10 代替本轮场。production particle physics、wall handoff/contact、自适应积分、accept/reject、FP64/GPU backend 均未修改；冻结合约覆盖的 40 个生产源码 hash 已在本地与服务器核对一致。

这次停止不能写成 PARTICLE_ORACLE_FAIL，也不能据此判断原生产粒子求解器正确或错误。恢复任何后续计算需要用户新的任务指令；本轮没有启动 BraVa 新 CFD、population 或旧队列。

机器状态：[PARTICLE_ORACLE_VALIDATION.json](PARTICLE_ORACLE_VALIDATION.json)。总入口：[OPEN_RESULTS.html](OPEN_RESULTS.html)。
