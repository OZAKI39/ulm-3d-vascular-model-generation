# 小鼠 A：Network-H0 基线

迁移后的目录：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-H0-pressure-v1`。

原状态：`RETAINED_QUALIFIED_CFD_AND_MATCHED_TRAJECTORIES`。本次只迁移整理，未重算流场或轨迹。

500 条正式微泡，dt=1 ms；保留其既有交付状态。

- 冻结流场：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-H0-pressure-v1/flow/frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu`。
- 轨迹与日志：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-H0-pressure-v1/particle_3d/reports/particle9a5_formal_trajectories`。
- 可视化：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-H0-pressure-v1/visualization/OPEN_RESULTS.html`。
- 服务器原位置：`/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1`。

流场 SHA256：`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。

`code/` 指向已迁入 `_shared/` 或对应 BraVa 算例的实际代码，避免为每个算例复制同一套库。旧路径兼容链接用于保留历史配置和哈希记录的可读性。

[全部算例与迁移说明](../README.md)
