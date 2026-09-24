# SimVascular / svMultiPhysics 3D FEM

当前主线是 Network A/H0 压力边界驱动的冻结 3D 流场。完整路径见[当前工作流索引](/home/lzy/projects/ACTIVE_VASCULAR_WORKFLOW.md)。

| 内容 | 入口 |
| --- | --- |
| H0 算例生成 | [/home/lzy/projects/ulm_3D_vascular/network_1d0d/fem_h0_case.py](/home/lzy/projects/ulm_3D_vascular/network_1d0d/fem_h0_case.py) |
| 当前远程求解 | [solve_a_h0_fem_remote.py](/home/lzy/projects/ulm_3D_vascular/scripts/solve_a_h0_fem_remote.py) |
| 当前物理验收 | [validate_a_h0_fem.py](/home/lzy/projects/ulm_3D_vascular/scripts/validate_a_h0_fem.py) |
| 测量、收敛、checkpoint | `src/sv_validation/` |
| 当前日志解析 | `scripts/sv13q/flow_parser.py` |
| 2 mm/s 算例、WSS 与流线派生量 | [flow_2mmps](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/scripts/flow_2mmps) |
| 当前流场与算例 | [/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases) |
| 旋转图/动画及原渲染接口 | [rotate_visualization](/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html) |

当前 H0 场在 `mean-2p0-mmps-A-H0-pressure-v1/frozen_flow/`，物理验收是该算例的 `reports/physics_validation_H0.json`。
旧 2 mm/s 场与各工作树的冻结输入有不同角色，使用 manifest 和 SHA 区分。

必要工具链重建脚本位于 `scripts/sv13n/` 与 `scripts/sv13l/`；Stage Q 运行脚本、source patch、C++ helper、配置与第三方求解器保留。旧阶段一次性试装/重试/交付脚本已退役，原始文件可从[清理归档](/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z)恢复。
`scripts/fem_freeze_sync/validate_frozen.py` 只对应早期 Stage Q 冻结交付，不能作为 H0 新场的验收入口。

当前派生量/物理门禁回归在流场工作树运行：`PYTHONPATH=src /home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python -m pytest -q -p no:cacheprovider tests/flow_2mmps`。
历史全量测试保留其历史语义；部分测试依赖原服务器路径或原始未镜像产物，不能替代当前入口的验收。
