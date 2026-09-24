# Network A/H0 → 3D FEM → Microbubble

这是截至 2026-09-25（柏林时间）的工作流同步分支，包含代码瘦身后的本地实现、必要输入、冻结流场、轨迹、可视化、正式报告，以及服务器源码和执行记录。

当前结果为 **NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY**：同一组 30 个初始状态，OLD/NEW 各 29 条完成、1 条 stationary、0 条数值失败。500 条新流场生产计算尚未启动。

| 内容 | 仓库入口 |
| --- | --- |
| 血管 A、ROI、端口、Network H0 | [vascular_network](vascular_network/README.md) |
| H0 3D FEM 算例、求解、校验 | [FEM_SimVascular](formal_3D_flow_solver/FEM_SimVascular/README.md) |
| 当前 OLD/NEW 配对轨迹与科学代码 | [particle_3d](particle_3d/README.md) |
| 配对验证正式报告 | [NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md](particle_3d/reports/network_derived_flow_mb_validation_v1/NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md) |
| 压力、WSS、流线、局部矢量原可视化接口 | [rotate_visualization](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html) |
| 服务器保留源码、配置、日志、执行证据 | [server_evidence](server_evidence/README.md) |
| SonoVue 固定分布 | [sonovue_size_distribution_v0](sonovue_size_distribution_v0/README.md) |
| 历史 FEM / RBC / 2D 来源 | `formal_3D_flow_solver/FEM/`、`particle_3d/reports/`、`source_dependencies/` |
| 同步范围、清理记录、逐文件来源与排除原因 | [同步报告](sync_metadata/network_h0_particle_20260925/SYNC_REPORT_ZH.md) |

## 下载与检查

```bash
git clone --single-branch --branch sync/network-h0-particle-cleanup-20260925 \
  https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git
cd ulm-3d-vascular-model-generation
python3 sync_metadata/network_h0_particle_20260925/verify_snapshot.py
```

完整性检查只需 Python 标准库；本分支使用普通 Git 对象。发布验证使用 Python 3.13.11，依赖版本见 [requirements-verified.txt](sync_metadata/network_h0_particle_20260925/requirements-verified.txt)。在具备这些依赖的环境中运行：

```bash
python sync_metadata/network_h0_particle_20260925/verify_snapshot.py --scientific-inputs
python sync_metadata/network_h0_particle_20260925/run_checks.py --output /tmp/ulm-workflow-checks
```

发布前在临时迁移目录中完成 **115 项回归，全部通过**（Particle 24、Network 49、流场 42）。测试使用仓库内文件；只在临时副本中迁移 5 份 JSON 的历史绝对路径，源码、数组、公式、参数及测试断言保持原字节。原科学清单中的 92 个解释器缓存条目在临时验证清单中排除，全部 114 个实际科学源文件仍校验原 SHA。原报告和原清单保留。

## 流场角色

| 路径 | 角色 |
| --- | --- |
| `formal_3D_flow_solver/FEM_SimVascular/frozen_reference/` | Particle 默认使用的历史 2 mm/s 场 |
| `formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/frozen_flow/` | OLD，SHA `129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d` |
| `formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow/` | NEW/H0，SHA `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4` |
| `formal_3D_flow_solver/FEM_SimVascular/upstream_stage_q_reference/` | Flow 原工作树的早期 Stage Q 冻结交付，用于历史流场回归 |
| `formal_3D_flow_solver/FEM_SimVascular/legacy_inputs/` | Particle 迁移前输入的历史角色 |

流场求解使用固定的 svMultiPhysics/PETSc GPU 工具链；微泡轨迹由 CPU 进程按泡并行。必要构建脚本、源补丁、版本与服务器记录已保留；已安装环境与二进制可按记录重建。

原始全脑数据、虚拟环境/构建树、多余渲染帧、中间 CFD 状态、重复部署包及非必要大文件未纳入。依赖当前回归的约 29.5 MiB 流线候选池作为重要文件保留。历史全流程脚本仍有原机器路径，重新计算前需按[迁移说明](sync_metadata/network_h0_particle_20260925/PORTABILITY_ZH.md)配置；本次没有重算 CFD 或轨迹。
