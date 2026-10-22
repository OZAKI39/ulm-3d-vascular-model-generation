# 实验控制与实时视觉平台 — 2026-10-22 同步快照

本分支包含实物实验相关软件：Chronos 相机传输、Tektronix 信号发生器、PZT 通道映射与状态、实时 GUI/录像、人工标注、SAM/Cellpose 辅助及微调、OpenCV 校准、PIV 与群体运动研究。数值仿真位于 [sync/microbubble-simulation-f10ag-20261022](https://github.com/OZAKI39/ulm-3d-vascular-model-generation/tree/sync/microbubble-simulation-f10ag-20261022)。

这是独立快照分支，保留实际工作区的未提交修复。视觉根目录的旧 HEAD 不能单独代表当前代码：逐文件源哈希见 [COPIED_FILES.json](sync/20261022/COPIED_FILES.json)。相机和 AFG 采用实际运行的 `F:\CHRONOS_control`、`F:\Tektronix_control` 最新文件；WSL 不同版本的控制日志另存为证据。

## 使用入口

| 部分 | 入口 |
|---|---|
| Chronos 接收、调焦与黑校准 | [CHRONOS_control](external_connect/CHRONOS_control/) |
| AFG 控制和状态查询 | [Tektronix_control](external_connect/Tektronix_control/) |
| 实时观测、PZT 布局和自动录像 | [realtime_collective_observation_v1](external_connect/vision_tracking_research/realtime_collective_observation_v1/) |
| 人工标注、候选生成和训练 | [smart_online_tracking_v2](external_connect/vision_tracking_research/smart_online_tracking_v2/) |
| 新显微镜标签、校准与 PIV 结果 | [microscope_20261019_v1](external_connect/vision_tracking_research/smart_online_tracking_v2/microscope_sessions/microscope_20261019_v1/) |
| 群体视觉研究 F10U–F10X | [vision_tracking_research](external_connect/vision_tracking_research/) |
| 一段真实输入视频 | [experiment_overlay_v1.mp4](external_connect/sample_videos/v2/experiment_overlay_v1.mp4) |

## 保留内容与状态

人工实例标签、审核区域、数据划分、候选接受/拒绝记录、已采用的 SAM decoder、最新显微镜候选权重、训练与独立验证记录均有保留。历史未通过候选的日志也保留，不能将其误认为已启用模型。OpenCV 离线候选的准入状态保持原记录，未改变实时默认检测器。

保留最新 OpenCV/PIV 数值数组、逐帧记录、时间戳、显示规则与脚本、最终预览，以及相机黑校准记录、PZT 配置、读回和控制日志。历史 `observation_status` 与查询请求仅归档在 `sync/20261022/controller_state_evidence`，没有放回活动命令邮箱，本次没有查询或更改仪器输出。

本地绝对目录链接已实体化；显微镜会话的共享模型目录改为等价相对链接。基础 SAM、Cellpose、RAFT 等可下载权重、环境、VLC/驱动安装包、合成训练大数组和完整长录像没有整批上传。模型身份与下载信息见原模块说明及排除清单。长实验录像仍是外部原始资料，本分支不是完整录像灾备。

## 校验与继续开发

```bash
python sync/20261022/restore_payload.py
python sync/20261022/restore_payload.py --restore
```

第二条命令按 SHA-256 还原无损压缩的文本/数据原文件名。环境需按各模块依赖文件重新创建；Windows 控制器与 WSL 视觉程序分别配置。原脚本中的 `F:`、`/home/lzy/projects`、视频和模型路径需对应实际安装位置。README 不意味着全新电脑已经即插即用；本次进行了离线打包/回归检查，没有在新设备上执行硬件联调。

详见 [EXPORT_POLICY.json](sync/20261022/EXPORT_POLICY.json)、[排除清单](sync/20261022/EXCLUDED_FILES.json)、[来源身份](sync/20261022/SOURCE_IDENTITIES.json)、[VALIDATION.json](sync/20261022/VALIDATION.json)。旧报告中的绝对路径和历史清单保持原样，以 `sync/20261022` 定义当前快照。
