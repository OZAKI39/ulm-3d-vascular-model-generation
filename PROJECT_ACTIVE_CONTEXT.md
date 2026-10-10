# 当前项目入口 — Y 血管血栓断裂与早期背景流场耦合

更新：2026-10-10。当前分支是血栓项目的新服务器迁移及 GitHub 精选同步，不是早期 PD2 开发终点。

- [当前完整说明、结果与复核方法](clot_clearance/github_sync_20261010/README.md)
- [当前科学状态](clot_clearance/verification/y_clot_point_cloud_and_early_flow_coupling_v1/FINAL_STATUS.json)
- [最新 4K 视频](clot_clearance/results/y_clot_point_cloud_and_early_flow_coupling_v1/whole_vessel_4k_v1/y_clot_whole_vessel_force_directions_4k.mp4)
- [原生 LAMMPS 修改](clot_clearance/github_sync_20261010/native_patches/lammps/SOURCE.json)

服务器独立目录 `/workspace/y_clot`，环境入口 `source /workspace/y_clot/env.sh`。完整逐节点力审计与原生依赖在服务器保留；GitHub 保留关键完整精度轨迹、网格流场、源代码、配置、紧凑证据与可视化。

方法边界：定性工程断裂；原生无血栓 Y 背景流场与碎片单向耦合；附着核心在 2.33 ms 后保持参考快照。未新增生理准确性或双向 FSI 结论。旧 PD2 边界工作及历史报告保留于分支祖先；上版根上下文存于 `clot_clearance/github_sync_20261010/PREVIOUS_ROOT_CONTEXT.md`。
