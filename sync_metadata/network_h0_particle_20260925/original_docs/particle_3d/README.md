# Particle：当前科学代码与运行入口

当前状态：**NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY**。Network H0 配对验证使用相同的 30 个初始状态，OLD/NEW 各 29 条完成、1 条 stationary、0 条数值失败。尚未启动 500 条新流场生产计算。

- 当前入口：[reports/network_derived_flow_mb_validation_v1/scripts/runner.py](reports/network_derived_flow_mb_validation_v1/scripts/runner.py)。
- 正式报告：[NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md](reports/network_derived_flow_mb_validation_v1/NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md)。
- 回归测试：`tests/network_flow_mb_validation_v1/`、`tests/particle9a1_2mmps/`。
- 共享实现：`src/particle_3d/`；固定科学快照与配对输入：`reports/network_derived_flow_mb_validation_v1/server_bundle/`。

## 核心调用链

`runner.mb_job → Particle9AStepper + particle9a1_audit.instrument → particle82a_integration.integrate_admitted → particle6_stepper / particle65_motion → 阻力与接触求解`。

`field.FrozenFEMField` 查询四面体并插值背景流；`wall_geometry`、`wall_gap`、`wall_contact` 处理壁面几何与间隙；`hydrodynamic_resistance`、`resistance_assembly`、`resistance_solver` 和 `planar_wall_hydrodynamics` 提供阻力/近壁模型。点示踪对照入口是 `particle82_point_native.NativePointTracer`。

轨迹在 CPU 上按泡分配给最多 6 个进程；各泡求解自身运动与接触。旧阶段模块仍被复用，文件名前缀不代表废弃。

## 复现与检查

默认 `formal_3D_flow_solver/FEM_SimVascular/frozen_reference` 仍是历史 2 mm/s 场；H0 场由当前 adapter 显式加载，不修改默认绑定。

在本工作树根目录运行当前回归：

```bash
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1 .venv/bin/python -m pytest -q -p no:cacheprovider particle_3d/tests/network_flow_mb_validation_v1 particle_3d/tests/particle9a1_2mmps
```

服务器复现目录为 `/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z`。使用 `/root/particle8_2_runs/env/bin/python scripts/runner.py --root <该目录> --label NEW --workers 6 --ids all --output <新的输出名>`，并设置上面的单线程环境变量。入口检查科学快照 SHA，且拒绝覆盖已有输出；输出名必须是该目录下尚不存在的子目录名。

历史 P9A1 的复现入口 `scripts/run_particle9a1.py`、原回放/渲染脚本和测试保留。当前分析与绘图代码在配对报告的 `scripts/` 下。旧部署、迁移和一次性收尾脚本已归档，按需从[清理归档](/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z)恢复。

跨仓库与服务器依赖详见[当前工作流索引](/home/lzy/projects/ACTIVE_VASCULAR_WORKFLOW.md)。
