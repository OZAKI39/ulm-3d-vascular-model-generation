# 当前结果入口

## 最新微泡：ROI-only-balanced-pressure-v1，dt = 0.5 ms

2026-09-28 已完成 1500 条：O1=787、O2=61、O3=467。终态：COMPLETED=1315、SUPPORTED_STATIONARY=185；执行失败 0。
保持上一批原有粒径分布及固定来源顺序，未按出口补样。数值检查、CUDA 复核、本地独立完整性核验均通过。

- [英文动画（主体放大 30%）、图件与数据入口](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/OPEN_RESULTS.html)
- [新增三个视角 MP4：侧向、反向、斜俯视](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/additional_views/OPEN_RESULTS.html)（同一批轨迹，英文界面，保留 1.30 倍显示）
- [中文结果报告](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/MICROBUBBLE_RESULTS_ZH.md)
- [0.5 ms 位置数据](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/data/trajectories_dt0p5ms.npz) · [逐轨迹 CSV](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/data/trajectory_catalog.csv) · [出口统计](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/data/outlet_summary.csv)
- [1500 条原始轨迹及审计](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/tracks/)
- [本地独立核验](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/LOCAL_VERIFICATION.json) · [复现说明](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/REPRODUCE.md)

8 CPU 进程执行有限尺寸积分，RTX 4090 用于 CUDA FP64 场/轨迹复核和 EGL 动画渲染；视频编码为 CPU libx264。
动画按轨迹年龄对齐，不代表同时注入；1.30 倍为相对原始显示比例，未重复叠加。所有既有结果继续保留。

## 上一批微泡：best-feasible-balance-v1，dt = 0.5 ms

2026-09-28 已完成 1500 条：O1=0、O2=133、O3=1182，接触支持静止 185 条；执行失败 0。保持原有粒径分布，按来源顺序取首 1500 个有效样本，没有按出口补样。

- [英文动画（主体放大 30%）、图件与数据入口](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/OPEN_RESULTS.html)
- [中文结果报告](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/MICROBUBBLE_RESULTS_ZH.md)
- [0.5 ms 位置数据](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/data/trajectories_dt0p5ms.npz) · [逐轨迹 CSV](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/data/trajectory_catalog.csv) · [出口统计](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/data/outlet_summary.csv)
- [1500 条原始轨迹及审计](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/tracks/)
- [本地独立核验](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/LOCAL_VERIFICATION.json) · [复现说明](particle_3d/reports/microbubble_best_balance_dt0p5ms_n1500/REPRODUCE.md)

8 CPU 进程执行有限尺寸积分；RTX 4090 执行 FP64 场/轨迹复核和动画渲染。视频按轨迹年龄对齐，不代表同时注入。9000 个结果文件及 114 个原有源码文件通过哈希核验；该检查不构成时间步收敛证明。

## 保留基线：Network-H0，dt = 1.0 ms

- [500 条正式轨迹](particle_3d/reports/particle9a5_formal_trajectories/outputs/formal/core500/tracks/)
- [逐轨迹目录](particle_3d/reports/particle9a5_formal_trajectories/data/trajectory_catalog.csv)
- [当前结果摘要](particle_3d/reports/particle9a5_formal_trajectories/data/final_summary.json)
- [当前图件](particle_3d/reports/particle9a5_formal_trajectories/figures/)
- [固定 CORE500](particle_3d/reports/particle9a5_formal_trajectories/data/CORE500_COHORT.json)

该基线结果为 500 条：O1=0、O2=139、O3=291，另有 70 条接触支持静止；没有计算失败。清理保留了所有轨迹和现有图件。原 OPEN_RESULTS.html 仍是较早的阶段页；正式动画和最终交付核验此前没有完成，本次清理不把它们标记为已完成。

## RBC 代码、轨迹与展示

| 内容 | 路径 |
|---|---|
| RBC 几何分布与轨迹 | [particle2](particle_3d/reports/particle2/) |
| 受限血管内形状/轨迹模型 | [particle3](particle_3d/reports/particle3/) |
| RBC 与微泡局部相遇 | [rbc_mb_interaction](particle_3d/reports/rbc_mb_interaction/) |
| 正常 RBC 相遇展示 | [rbc_mb_normal](particle_3d/reports/rbc_mb_normal/) |
| 共流展示 | [rbc_mb_coflow](particle_3d/reports/rbc_mb_coflow/) |
| 运动与旋转展示 | [rbc_mb_flow_rotation](particle_3d/reports/rbc_mb_flow_rotation/) |
| HemoCell 代码和结果副本 | [hemocell_restore](server_evidence/hemocell_restore/) |

这些 RBC 结果使用各自保存的参考场/理想化模型；不能解释为已完成新 H0 流场中的 RBC 耦合计算。代码位于 `particle_3d/src/particle_3d/`，生成和渲染入口位于 `particle_3d/scripts/`。
