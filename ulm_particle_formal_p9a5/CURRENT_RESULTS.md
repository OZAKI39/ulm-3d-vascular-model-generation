# 当前结果入口

## 微泡：Network-H0，dt = 1.0 ms

- [500 条正式轨迹](particle_3d/reports/particle9a5_formal_trajectories/outputs/formal/core500/tracks/)
- [逐轨迹目录](particle_3d/reports/particle9a5_formal_trajectories/data/trajectory_catalog.csv)
- [当前结果摘要](particle_3d/reports/particle9a5_formal_trajectories/data/final_summary.json)
- [当前图件](particle_3d/reports/particle9a5_formal_trajectories/figures/)
- [固定 CORE500](particle_3d/reports/particle9a5_formal_trajectories/data/CORE500_COHORT.json)

现存结果为 500 条：O1=0、O2=139、O3=291，另有 70 条接触支持静止；没有计算失败。清理保留了所有轨迹和现有图件。原 OPEN_RESULTS.html 仍是较早的阶段页；正式动画和最终交付核验此前没有完成，本次清理不把它们标记为已完成。

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
