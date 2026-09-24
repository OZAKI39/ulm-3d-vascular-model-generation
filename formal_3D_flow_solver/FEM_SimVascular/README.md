# 3D FEM / SimVascular

当前是 Network A/H0 出口压力驱动的冻结 3D 血管流场。

- H0 算例生成：[fem_h0_case.py](../../vascular_network/network_1d0d/fem_h0_case.py)。
- 远程求解：[solve_a_h0_fem_remote.py](../../vascular_network/scripts/solve_a_h0_fem_remote.py)。
- 物理校验：[validate_a_h0_fem.py](../../vascular_network/scripts/validate_a_h0_fem.py)。
- 共享实现：`src/sv_validation/`；当前日志解析：`scripts/sv13q/flow_parser.py`。
- OLD 与 H0 算例：[flow_cases](flow_cases/)；当前派生量入口：`scripts/flow_2mmps/`。
- 原旋转可视化接口：[OPEN_RESULTS.html](rotate_visualization/OPEN_RESULTS.html)。
- 固定求解器源码：[vendor/svMultiPhysics_stage_q](vendor/svMultiPhysics_stage_q/)，其三个修改文件与本次服务器实际源文件逐字节相同；运行记录见 [server_evidence](../../server_evidence/README.md)。

`frozen_reference` 供当前 Particle 使用；`upstream_stage_q_reference` 保留 Flow 工作树的旧冻结输入。两者角色不同。独立回归脚本在临时目录选择对应输入，避免改写任何已有冻结目录。
Stage N 的 PETSc/CUDA 与 Stage L 的 MPI/Fortran 构建脚本、启动 wrapper 和 Stage Q reuse patch 仍有实际用途，已保留。

完整性和当前 42 项流场测试包含在[主 README](../../README.md)的统一检查命令中。Taylor–Hood 验证已因资源成本停止，代码/日志作为历史证据保留，不代表已有完整 P2/P1 稳态场。
