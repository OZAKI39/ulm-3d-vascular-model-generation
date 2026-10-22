# WSS 验证：pipe_nr16_cpu_mpi8_halfdt

迁移后的目录：`/home/lzy/projects/computation_examples/validation_cases/pipe_nr16_cpu_mpi8_halfdt`。

原状态：`PASS`。本次只迁移整理，未重算流场或轨迹。

独立数值验证算例；没有生产微泡轨迹。细档及 O2 压力扰动仍为 SKIPPED_BY_USER。

- 冻结流场：`/home/lzy/projects/computation_examples/validation_cases/pipe_nr16_cpu_mpi8_halfdt/flow/frozen_flow/steady_flow.vtu`。
- 服务器原位置：`/workspace/wss_validation_v2_20260927T1230Z/stage2/pipe_nr16_cpu_mpi8_halfdt`。

流场 SHA256：`f3444feddb2a210a0a26a043fb4ec6579e3ba71264d414cf50428fb6bff9aeda`。

`code/` 指向已迁入 `_shared/` 或对应 BraVa 算例的实际代码，避免为每个算例复制同一套库。旧路径兼容链接用于保留历史配置和哈希记录的可读性。

[全部算例与迁移说明](../../README.md)
