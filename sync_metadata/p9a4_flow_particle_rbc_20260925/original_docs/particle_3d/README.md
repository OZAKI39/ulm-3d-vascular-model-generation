# 当前 Particle 科学代码

当前任务是 Network H0 新旧流场的 30 泡配对验证，状态 `NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY`。OLD/NEW 各 29 条完成、1 条 stationary、0 条数值失败。

- 运行入口：[runner.py](reports/network_derived_flow_mb_validation_v1/scripts/runner.py)。
- 正式报告：[NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md](reports/network_derived_flow_mb_validation_v1/NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md)。
- 科学实现：`src/particle_3d/`；原科学快照与固定输入：`reports/network_derived_flow_mb_validation_v1/server_bundle/`。
- 核心调用：`Particle9AStepper → particle82a_integration.integrate_admitted → particle6_stepper / particle65_motion → 阻力与接触求解`。轨迹在 CPU 上按泡并行。

`FrozenFEMField` 进行四面体场查询，壁面/阻力/接触共享模块及历史数学回归均保留。点示踪对照使用 `NativePointTracer`。旧阶段名称仍包含实际依赖，未按编号删除科学核心。

从仓库根目录运行 [run_checks.py](../sync_metadata/network_h0_particle_20260925/run_checks.py)，可完成包含当前 Particle 24 项的完整 115 项回归；它在临时目录迁移历史路径，不修改本仓库数据。

新计算必须指定新的输出名。CPU 配对复现可使用：

```bash
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1 \
python particle_3d/reports/network_derived_flow_mb_validation_v1/scripts/runner.py \
  --root particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle \
  --label NEW --workers 6 --ids all --output replay_new_001
```

该命令会实际重新积分，入口会验证科学快照哈希并拒绝覆盖已有输出。本次同步仅运行回归与完整性检查。本次验证的确切依赖版本见 [requirements-verified.txt](../sync_metadata/network_h0_particle_20260925/requirements-verified.txt)。

各阶段历史数据、RBC 示例与独立升力审核保留其原科学适用范围。跨仓库入口见[主 README](../README.md)。
