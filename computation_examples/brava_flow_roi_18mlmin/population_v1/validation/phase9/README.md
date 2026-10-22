# Phase 9：实际性能、生产与交付

`benchmark.py` 使用真实 BraVa FEM 中选出的内部非交互测试群体，N=1/10/100/1000/10000，CPU/GPU 同一状态、同一完整 0.5ms 推进。夹具是性能输入，不冒充实际灌注；没有活跃球对时 pair physics 成本为零。阶段 CUDA events 预热取中位数，完整步单次计时，不作稳定吞吐置信区间。

`benchmark_sparse_batch.py` 16 个独立 realization 各 1 泡，逐组/批处理端点 parity。`launch_production.py` 记录正式 run 的命令、PID 和资源遥测；实际结果用自己的统计，不套用夹具速度。

两次开发运行揭示共享时钟进展保护的责任范围问题，保留失败记录。最终修复增加两项回归检查，实际近共线接触测试补充后26测试通过。正式 run `realization_0000_contactfix` 仍必须区分每条轨迹的正常/失败/静止/删失，不因全局覆盖达标把科学完整性写成 PASS。

复现命令见根目录 `REPRODUCE.md`；全部数值与结论在主报告 `POPULATION_VALIDATION_REPORT_ZH.md`。动画使用真实 global physical time、英文界面、既有样式和精确1.30倍正交主图。
