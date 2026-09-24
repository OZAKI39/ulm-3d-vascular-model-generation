# 3D vascular flow, microbubble and RBC workflow

本分支整理截至 2026-09-23 的血管几何、FEM/SimVascular 流场、微泡轨迹、RBC 运动、可视化和科学审核结果，保留对应代码、必要输入、正式数据、服务器运行记录及测试。它是当前工作成果的精简快照，不是全部原始数据和临时计算目录的镜像。

**完整同步范围、验证记录与排除说明：** [中文同步报告](sync_metadata/flow_mb_rbc_20260923/SYNC_REPORT_ZH.md)。

## 从这里开始

| 内容 | 仓库位置 |
|---|---|
| 血管几何生成与有效 ROI 输入 | `model_generate.py`、`utils/`、`test_data/`；[原几何工程说明](sync_metadata/flow_mb_rbc_20260923/source_README.md) |
| FEM/SimVascular 求解代码、配置、构建补丁、服务器记录 | [formal_3D_flow_solver/FEM_SimVascular](formal_3D_flow_solver/FEM_SimVascular/) |
| Particle 使用的旧冻结流场与交接契约 | [frozen_reference](formal_3D_flow_solver/FEM_SimVascular/frozen_reference/)；[PARTICLE_HANDOFF](formal_3D_flow_solver/FEM_SimVascular/PARTICLE_HANDOFF.md) |
| 新的 2.0 mm/s 完整求解、最终场、残差/Excel、压力、WSS | [flow_cases/mean-2p0-mmps](formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/) |
| 最新全血管流线、局部速度矢量、压力与 WSS 旋转显示 | [rotate_visualization/README_ZH.md](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/README_ZH.md) |
| 微泡及 RBC 的运动算法、永久测试、阶段脚本 | [particle_3d](particle_3d/)：`src/`、`tests/`、`scripts/` |
| 正式微泡轨迹与约 500 条轨迹的 PPT 数据 | [particle_3d/outputs](particle_3d/outputs/) 中 `particle8_1`、`particle8_2a`、`particle8_2a_ppt` 等 |
| RBC/微泡共同流动及旋转 | [rbc_mb_flow_rotation](particle_3d/reports/rbc_mb_flow_rotation/)，另保留 `rbc_mb_coflow`、`rbc_mb_normal`、`rbc_mb_interaction` |
| 各阶段 CSV/JSON、图表、动画、中文审核报告 | [particle_3d/reports](particle_3d/reports/) |
| 独立剪切升力量级审核 | [shear_lift_audit](particle_3d/shear_lift_audit/) |
| 冻结 SonoVue 尺寸分布 | [sonovue_size_distribution_v0](sonovue_size_distribution_v0/) |
| 早期 FEM、来源迁移与研究总结 | [formal_3D_flow_solver/FEM](formal_3D_flow_solver/FEM/)、[source_dependencies](source_dependencies/)、[THESIS_SUMMARY](THESIS_SUMMARY/) |

## 流场与模型版本

- `frozen_reference` 是现有正式微泡轨迹使用的旧冻结流场。新的 `flow_cases/mean-2p0-mmps` 是入口平均速度 2.0 mm/s 的独立求解结果；不能把旧场轨迹当作新场轨迹。
- 正常形态 RBC 与微泡共同移动/旋转的最新演示使用较宽理想化血管中的 Poiseuille 场；它不等同于真实血管内解析 RBC 膜变形的双向流固耦合。
- 近壁和几何变形代理仍沿用原阶段模型定义。独立剪切升力审核没有向正式轨迹方程加入升力。
- 本次同步未修改运动方程，也未重新进行完整 CFD 求解或大批轨迹积分。

## 下载与验证

```bash
git clone --single-branch --branch sync/flow-microbubble-rbc-workflow-20260923 \
  https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git
cd ulm-3d-vascular-model-generation
python3 sync_metadata/flow_mb_rbc_20260923/verify_snapshot.py
```

上面的 SHA256 检查只需 Python 标准库。使用原有科学 Python 环境后，还可验证冻结契约和运行本次检查：

```bash
python sync_metadata/flow_mb_rbc_20260923/verify_snapshot.py --scientific-inputs
python sync_metadata/flow_mb_rbc_20260923/run_checks.py
```

Particle 包要求 Python >= 3.11；请按各模块的 `pyproject.toml`、依赖锁文件和 README 准备环境。旋转可视化的精确依赖版本保存在 `rotate_visualization/requirements.txt`。这里保留的是经过验证的代码和依赖记录，不包含本地 `.venv` 或服务器已编译的运行环境。

## 数据边界

正式轨迹、冻结网格/流场、必要检查点、最终 MP4/图片与结果摘要已保留。大型准入缓存、重复汇总 CSV、升力审核逐状态派生数组、中间 CFD 状态、多余渲染帧、部署压缩包和原始全脑数据集未同步。详见[逐文件清单](sync_metadata/flow_mb_rbc_20260923/file_inventory.csv)和[恢复说明](sync_metadata/flow_mb_rbc_20260923/REBUILD_OMITTED_ZH.md)。

历史报告中的绝对路径、原始 SHA256 清单和阶段验收结论按原样保留，因此其中可能出现本次排除的文件。当前快照以本次 `workflow_files.sha256` 为准；迁移到新机器后，重跑历史脚本前需更新配置中的工作目录和服务器路径。
