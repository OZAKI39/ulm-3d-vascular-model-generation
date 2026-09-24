# Taylor-Hood 验证停止记录

状态：`TAYLOR_HOOD_VALIDATION_STOPPED_RESOURCE_COST`。

本轮由用户决定停止，原因是完整规模验证的求解器与资源工程成本已经超出当前项目阶段的合理投入。这不是 Taylor-Hood 的科学失败。

已经验证：原生 TET10 / P2 velocity / P1 pressure、nFs=2、VMS 关闭；TET4→TET10 几何与全局共享边拓扑；CUDA 稀疏矩阵；官方制造解短验证。实际全模型为 534,979 节点、371,402 四面体、约 235,823,344 nnz（P1 的 14.74 倍）。原 ILU2 构建超过 25 分钟且主存达到约 36 GiB，已按特定 PID 安全停止；RCM + ILU2 随后发生 CUDA allocation OOM。

用户停止指令到达前，最后一次已启动的 RCM + shift + ILU1 短 smoke 自行正常退出：1 个时间步、1 次 Newton 迭代，17 次线性迭代，true relative residual 约 4.40e-11。该结果仅证明此次单步代数求解成功，不能称为完整稳态 P2 流场，也不能证明局部守恒改善。

尚未验证：完整 P2 稳态求解、完整物理验收、P2 局部守恒改善与资源成本收益。最后一次 smoke 的 result_001 和 checkpoint 只作原始存档，不用于 P2/P1 科学比较。不再启动任何 Taylor-Hood attempt、预条件器调试或 grad-div。

2026-09-24 21:00 UTC 核对本机与 Vast.ai 的 PID、命令、cwd、PPID、启动 ticks：无属于本轮 Taylor-Hood workspace 的活动进程，因此此时无需发送新信号。原始失败记录和全部脚本、报告、网格及结果均保留，不删除、不覆盖。

本机报告：本目录。服务器工作区：`/workspace/flow_mean_2p0_mmps_A_H0_TaylorHood_20260924T194503Z`。最新原始 smoke 位于 `case_int32/smoke_shift_rcm_ilu1`。其完整只读副本归档到新 MB 报告的 `archive/taylor_hood_last_completed_smoke.tar.gz`。新任务的 `logs/taylor_hood_termination.json` 与保护检查记录保存停止依据。

未来仅在另行授权恢复后，从现有 source capability / native function-space / mesh elevation audit 及真实资源记录开始评估新的硬件或求解器工程预算；不得把单步 smoke 当完整科学结果。已存在的 `solve_p2_remote.py` 是未执行的完整求解入口，本次不运行。
