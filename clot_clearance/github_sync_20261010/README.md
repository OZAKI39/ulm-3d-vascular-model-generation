# 当前 Y 血管血栓项目同步（2026-10-10）

这是当前 LAMMPS 点云断裂与早期背景流场耦合项目的代码、关键数据及证据快照。完整工作环境位于新服务器 `/workspace/y_clot`。本目录的迁移验证记录是已经完成的结果；GitHub 不包含软件环境和原生依赖二进制。

## 当前结果

- [最新 4K 视频](../results/y_clot_point_cloud_and_early_flow_coupling_v1/whole_vessel_4k_v1/y_clot_whole_vessel_force_directions_4k.mp4)
- [加载与阻力方向图](../results/y_clot_point_cloud_and_early_flow_coupling_v1/whole_vessel_4k_v1/y_clot_whole_vessel_force_directions_4k_load_and_drag.png)
- [三视角离散化图](../results/clot_discretization_combined_v2/Fig_clot_discretization_three_views.png)
- [当前科学评审](../verification/y_clot_point_cloud_and_early_flow_coupling_v1/REVIEW_Y_CLOT_POINT_CLOUD_AND_EARLY_FLOW_COUPLING_V1.md)
- [当前最终状态](../verification/y_clot_point_cloud_and_early_flow_coupling_v1/FINAL_STATUS.json)

保留 1,927 个血栓材料点、1,561 个固定锚点、157,637 条初始键。当前已接受实验终态脱离 217 个节点，形成 208 个脱离分量。501 帧的全部节点坐标和拓扑保持原数据精度，不做降采样。原生 Y 网格含 26,712 个节点、146,230 个四面体。

这是定性工程断裂/输运模型；流场为无血栓背景 CFD。仍不代表生理血栓参数标定、血栓阻塞后流场或双向流固耦合。传输阶段附着主体沿用 2.33 ms 的固定参考快照。

## 同步范围

当前完整精度 `accepted_trajectory.npz`（21.6 MB）和 `FRAGMENT_TRANSPORT_TRAJECTORIES.npz`（3.1 MB）是必要数据，故保留。几何、断键、准确事件帧、交接信息、Y 原生网格与速度/压力场、输入、配置、源代码、分析脚本、紧凑报告/日志和最新 4K 展示一并保存。

不上传近 1 GB 的逐节点 `saved_frame_force_audit.bin`、其他大型力审计、重复 LAMMPS 文本轨迹、大型 CSV、重复 ZIP、旧版视频、编译产物、运行环境、第三方完整源码树和商业字体。它们的当前完整版本保留在新服务器。较早实验的完整原始数据不都在此 GitHub 快照中，原有报告可能引用未上传的文件。

[源文件选择清单](SOURCE_FILES.json) 记录每个选定文件的源路径、字节数及 SHA-256；[未上传文件清单](EXCLUDED_FILES.json) 记录项目文件大小和理由。此分支继承旧血栓分支的历史，本轮不重写历史大文件。

原生 LAMMPS 修改以 [上游提交与文件哈希](native_patches/lammps/SOURCE.json)、补丁及三个修改后文件保存，避免上传整个第三方仓库。CFD 已有修补与相关工作流位于 `../github_dependencies/`。依赖目录的路径适配是本轮唯一生产 Python 变更，[记录在此](PATH_ADAPTATION.json)；科学算法和已接受数据保持不变。

## 小型受力显示缓存

[DISPLAY_FORCE_CACHE.npz](DISPLAY_FORCE_CACHE.npz) 保存已有渲染器在所有 501 个已接受帧使用的合外力、作用位置、代表节点 482 的阻力，以及逐帧坐标哈希。[提取来源](DISPLAY_FORCE_CACHE_PROVENANCE.json) 和 [提取代码](export_force_cache.py) 均保留。

缓存来自现有原始力记录和原有输运力表达式，不包含全部逐节点受力，不能替代完整力学审计。它允许在不下载近 1 GB 审计文件的情况下复核并重绘已有力箭头，没有新增求解、插值或数值舍入。

## 环境与复核

Python 3.13.11 与确切依赖列表见 [锁定文件](migration/requirements_new_server.lock.txt)。新电脑应单独创建虚拟环境，并自行准备需要的原生 LAMMPS/MPI/PETSc/SimVascular 依赖；GitHub 快照不是一键完整环境镜像。

在已配置的新服务器上，可以使用现有独立环境验证任意新 checkout：

```bash
source /workspace/y_clot/env.sh
cd /path/to/checkout/clot_clearance
export PYTHONPATH="$PWD"
export Y_CLOT_FONT_DIR=/workspace/y_clot/deps/fonts
python github_sync_20261010/replay_compact.py --output /tmp/clot_replay_new
```

输出目录必须不存在。该命令核对全部 501 帧和 15 个原生资产，然后按已有风格绘制第 52 帧（4K）。`--verify-only` 跳过绘图，`--frame N` 可选择其他已接受帧。完整逐节点审计没有上传，因此原完整 `render.py` 入口需要从服务器取回所列大文件；使用这里的 compact replay 可直接重绘选定帧。

`Y_CLOT_DEPENDENCY_ROOT` 可覆盖依赖路径。默认使用本 checkout 的 `github_dependencies`，正常部署仍支持原 `/workspace/y_clot/deps` 布局。`Y_CLOT_FONT_DIR` 可指定合法 Arial 字体目录；缺少 Arial 时若系统安装 Liberation Sans，会使用对应替代字体。原始上传图片/视频保持原字体和字节不变。

[迁移验证](migration/) 包括 MPI、LAMMPS run 0、原生网格 API、CPU/GPU 各两步流场对照和 4K 可视化。两步验证只确认环境/后端一致性，未重新认证稳态收敛或物理精度。构建脚本记录新服务器的实际路径与构建参数，外部源码/运行库需要另行准备。

## GitHub 推送与本地备份清理

科学快照提交 `eadd653f2c686ae8d5e699d542238d86611aa6b0` 已推送，并通过 GitHub 远端提交/树核验（[记录](REMOTE_SNAPSHOT_VERIFICATION.json)）。随后按用户要求删除本地完整迁移备份和迁移暂存副本，合计约 40.7 GB（[删除记录](LOCAL_BACKUP_DELETION.json)）。原本地开发目录、新分支 checkout、新服务器完整项目及轻量审计记录保留。此删除记录作为后续小型提交附在同一分支。
