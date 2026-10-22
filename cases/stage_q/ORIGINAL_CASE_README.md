# 早期 FEM Stage-Q 冻结参考

迁移后的目录：`/home/lzy/projects/computation_examples/fem_reference_stage_q`。

原状态：`FROZEN_REFERENCE_NOT_REVALIDATED_AS_CURRENT`。本次只迁移整理，未重算流场或轨迹。

保留历史网格、配置与冻结结果；不是当前 ROI-only 或 BraVa 流场。

- 冻结流场：`/home/lzy/projects/computation_examples/fem_reference_stage_q/flow/flow/steady_flow_stage_sv1_3q.vtu`。

流场 SHA256：`373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1`。

`code/` 指向已迁入 `_shared/` 或对应 BraVa 算例的实际代码，避免为每个算例复制同一套库。旧路径兼容链接用于保留历史配置和哈希记录的可读性。

[全部算例与迁移说明](../README.md)
