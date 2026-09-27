# P9-A.5 复现入口

本轮名义积分步长为 **1.0 ms = 0.001 s**，依据用户后续明确指令，覆盖最初请求中的 0.25 ms。原 P9-A.1/P6.5 动力学文件保持原字节；正式调用显式传入 `dt=DT`。近壁接触及连续安全证书仍可触发更小的内部子步。

原 0.25 ms 批次已停止，原服务器目录、已完成轨迹、未完成调试证据及当时源码均保留。其位置见 `data/dt_revision.json`。新旧步长使用相同 CORE500，但不同运行身份，完成记录不能跨步长复用。

## 当前运行与输入

- 当前本地、服务器和归档目录：`data/run_context.json`。
- 源版本：`37c40d7c08ca6527c8379c0ffa81436ae708f59f`。
- NEW 流场实际 SHA：`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。
- 永久 CORE500 SHA：`851c3a4893031334526942fed8b8d430c95245da324eafc7ac7b8600033cfdb9`。
- 准入数据来自已审核的 P9-A.4 accepted births；本轮不生成新 population。
- 原始入口 ledger 逐事件关联核验：`data/population_ledger_linkage_audit.json`。
- 生产合约：`../../contracts/P9A5_FORMAL_PRODUCTION_V1.json`。

## 环境

本地 Python：`/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`。
服务器 Python：`/root/particle8_2_runs/env/bin/python`。
所有积分与测试均设置 `PYTHONDONTWRITEBYTECODE=1`，并将 `OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS`、`MKL_NUM_THREADS`、`NUMEXPR_NUM_THREADS` 设为 `1`。服务器实际硬件、cgroup 配额、依赖版本保存在 `logs/server/`。

## 执行顺序

在工作树根目录执行本目录 `scripts/` 下的脚本。另行复现实验必须先创建全新的 worktree、归档和 UTC 服务器目录，并建立新的 `data/run_context.json`；不要在已发布的原始证据目录覆盖文件。

1. 原有165项回归：`sync_metadata/p9a4_flow_particle_rbc_20260925/run_checks.py --output <新的日志目录>`。原路径迁移只发生在该工具创建的临时测试副本。
2. `prepare.py` 核验输入和历史1.5秒上限，冻结 CORE500、固定24条benchmark队列及本轮正式合约。
3. `pytest particle_3d/tests/particle9a5_formal_trajectories` 执行新增永久测试，保存 JUnit XML 和文本日志。
4. `inventory_server.py` 创建新服务器目录并保存资源清单。
5. `deploy.py` 一次部署并启动固定24条轨迹的1/2/4/6/8/12/16 workers测速。原积分器以相同1.0 ms步长另算参考轨迹，检查调度适配器逐点一致性。完成前不得选定正式worker数量。
6. `start_production.py` 在回归与benchmark通过后，依据实测资源预先声明最大N；复用新步长、同一身份、逐文件SHA匹配的24条benchmark轨迹，再完成CORE500。缺少自然出口则依次加入连续accepted births，到2000前每批250，之后每批500，止于全部出口出现或预定资源上限。
7. 同一个活动积分器仅在需要时延长3→6→12秒，并保存不变的旧前缀。完成轨迹按SHA复用；被中断的未完成轨迹保存到 `interrupted_or_mismatched/`，随后从原始出生状态重算该条。
8. 启动器在最终MB cohort冻结后自动执行 `run_points.py`，对所有相同出生中心进行事后point诊断。point不影响正式成员。
9. `collect.py` 增量收集并验证全部服务器输出SHA；计算中可使用 `--progress` 收集日志及已完成批次摘要。
10. `compare_dt.py` 生成固定24条新旧步长的配对差异记录；它只审查变化，不构成时间收敛证明。
11. `analyze_formal.py`、`render_formal.py`、`render_animations.py` 生成实际统计、PNG/PDF和视频。
12. `final_regression.py` 重跑全部旧测试与新增测试；`verify_protection.py` 检查所有原文件SHA和Git差异。刷新分析并运行 `build_review.py` 生成中文报告及HTML。最终摘要中的保护状态只可在核验通过后写入。

原始NPZ、全量gzip诊断、逐轨迹metrics/support/receipt、批次manifest、完成标志及失败调试证据均属于可审核资产。CORE500与FULL须分别报告；达到时间上限不能自动标记为stationary。绘图抽稀仅作用于显示路径顶点，所有正式成员及原始数据保留。
