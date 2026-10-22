# 小鼠 A：ROI-only 平衡压力

迁移后的目录：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1`。

原状态：`RETAINED_QUALIFIED_CFD_AND_MATCHED_TRAJECTORIES`。本次只迁移整理，未重算流场或轨迹。

1500 条微泡，dt=0.5 ms；保留支持静止与实际出口覆盖。

- 冻结流场：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/flow/frozen_flow/steady_flow_mean_2p0_mmps_A_ROI_only_balanced_pressure.vtu`。
- 轨迹与日志：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500`。
- 可视化：`/home/lzy/projects/computation_examples/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/visualization/OPEN_RESULTS.html`。
- 服务器原位置：`/workspace/flow_roi_only_balance_20260928/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1`。

流场 SHA256：`fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12`。

`code/` 指向已迁入 `_shared/` 或对应 BraVa 算例的实际代码，避免为每个算例复制同一套库。旧路径兼容链接用于保留历史配置和哈希记录的可读性。

[全部算例与迁移说明](../README.md)
