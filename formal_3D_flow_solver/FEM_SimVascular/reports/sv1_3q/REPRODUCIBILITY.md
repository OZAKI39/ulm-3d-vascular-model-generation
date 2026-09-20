# Stage Q 证据复查入口

复查使用现有证据，不能因 pytest 或日志解析而重新启动 CFD。

- `USER_REQUEST.txt` 保存本阶段要求，`configs/sv1_3q/policy.json` 是运行前固定的规则。
- `reference_freeze.json` 和 `baseline_r1.json` 保存 Stage P 交付哈希、R1 晚期窗口、完整计算与早期历史区间。R1 没有重跑；早期历史区间的时间分辨率与启动开销差异单独说明。
- `early_checkpoint.json` 保存 Stage P 原生 step10 检查点的独立结构检查与 SHA256；晚期仍为原 Stage N step60 检查点。
- `patches/sv1_3q/ilu_rebuild_policy.patch` 是相对冻结 Stage P 的完整补丁；只增加重建策略及原生停止确认时间戳。
- `source_patch.json` 保存修改前后源码逐文件哈希。新二进制、构建缓存与远端库链接保存在 `native_artifact_mirror.json` 和 `remote/svmp_reuse_build.json`。
- `configs/sv1_3q/runplans/` 是逐次不可覆盖的运行计划；`remote/*_execution.json` 保存实际命令、全部求解原因、进程耗时、原生输出哈希、GPU 采样与退出状态。
- `logs/sv1_3q/remote/` 保存原始日志，包括失败尝试。`SV13Q_BEGIN/END/RECOVERY` 逐次记录年龄、参考迭代、重建原因、恢复及计数。
- `remote/adaptive_recovery_probe.json` 是使用同一生产头文件的独立 GPU 合成故障测试。它证明一次恢复和失败后停止分支，不是血管稳定性或性能证据。初版测试程序的调用/清理错误日志保留。
- 每份验收证据核对原始 PETSc `MatLUFactorNum` 与 `KSPSolve` 次数，独立进程重新读取最终 VTU，并检查原生 checkpoint。PETSc 嵌套事件不相加。
- 只对晚期前两名运行早期窗口；`winner.json` 根据两窗口总实测耗时选优。只有至少 10% 收益且两窗口健康才运行一次完整血管计算。
- 完整运行沿用原 `SteadyStopMonitor` 和阈值；只合并压缩文件传输，并在判断通过后立即发 STOP_SIM。`stop_latency.json` 保存检测、请求、原生确认三个时间点。

只读测试命令：

```bash
.venv/bin/python -B -m pytest -q -p no:cacheprovider tests/test_sv13q_*.py
```

完整测试的历史失败单列，Stage P 和 CPU production 均保持不变。各策略只有一次开发观测，不是正式性能统计；CPU/GPU 科学等价验证继续推迟。

最终源码复现应以完整补丁 `ilu_rebuild_policy.patch` 和最终构建清单为准。中间增量构建与测试程序修正脚本仅保留工程过程，不能在完整补丁上重复叠加应用。
