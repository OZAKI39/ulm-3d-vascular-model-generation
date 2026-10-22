# BraVa 等流量出口校准

迁移后的目录：`/home/lzy/projects/computation_examples/brava_calibration_equal_split`。

原状态：`SOLVED_CALIBRATION_NOT_FINAL_PRESSURE_BC`。本次只迁移整理，未重算流场或轨迹。

强制等流量出口用于校准压力；没有配套生产微泡，不能把它当作自然分流结果。

- 冻结流场：`/home/lzy/projects/computation_examples/brava_calibration_equal_split/flow/frozen_flow/steady_flow.vtu`。
- 服务器原位置：`/workspace/brava_flow_roi_18mlmin_20260928/cases/calibration_gpu8_production`。

流场 SHA256：`f5d7a2ccf7cdcc08cbfebd725d0a17e23e5780df23e4ef6488a26ffa5983d4d2`。

`code/` 指向已迁入 `_shared/` 或对应 BraVa 算例的实际代码，避免为每个算例复制同一套库。旧路径兼容链接用于保留历史配置和哈希记录的可读性。

[全部算例与迁移说明](../README.md)
