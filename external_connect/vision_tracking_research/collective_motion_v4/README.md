# F10X-PIV — 全通道微泡运动观测

已完成现成 CPU/GPU 后端部署、937 帧全通道离线分析、实际 PTS 因果回放、9 组科研图和 3 段视频。结论为**工程原型完成、科学准入部分通过**：空间范围和同条件后端速度提高，但尚未可靠达到 30 fps，仍有低对比度漏检及非仿射错误接受。

EXPLORATORY · PHYSICS_UNCALIBRATED · training_ready=false。输出是二维可见微泡图案的表观运动；不是水流真值、物理融合、纯声学响应或控制准入。

- [中文开发报告：33 项回答](F10X_DEVELOPMENT_REPORT.md)
- [最终推荐](FINAL_RECOMMENDATION.md)；[机器可读 Gate](FINAL_GATE.json)
- [科研图](figures/figure_manifest.json)；[完整视频](videos/full_video_manifest.json)；[后端比较视频](videos/backend_video_manifest.json)
- [活动代码索引](ACTIVE_CODE_INDEX.md)；[源码复用与许可证](SOURCE_REUSE_AUDIT.md)
- [恢复清单](RECOVERY_MANIFEST.json)；[输入输出哈希](INPUT_OUTPUT_SHA256.json)

正式离线结果：`results/full_video_verified_v5/results/offline_regular/`；正式定速结果：`results/replay_final_v6/results/`。`full_video_fast_v2` 和 `full_video_verified_v3` 是未通过全量数值等价检查的加速试验，不得作为正式结果。科研冻结配置、独立程序种子及最终测试在 `results/research_v1/`。

本地：`/home/lzy/projects/computation_examples/_archives/mevo_vision_f10x_20261017/external_connect/vision_tracking_research/collective_motion_v4`。服务器：`/workspace/f10x_piv_20261017`。本轮分支 `dev/f10x-full-channel-piv-20261017`，不推送 GitHub。

本地检查：`.venv/bin/python -m pytest -q tests`。依次重建报告：`scripts/analyze_delivery.py`、`scripts/publication_figures.py`、`scripts/render_videos.py`、`scripts/render_backend_video.py`、`scripts/write_reports.py`。绘图/编码只读保存的正式计算结果，不连接硬件。

服务器环境 `/workspace/f10x_piv_20261017/.venv`。用 `scripts/remote_job.py` 建立当前 Git 的不可变快照和唯一 Supervisor 标签；新计算必须使用新的输出目录。数据环境变量参见 `src/xcommon.py`，不能覆盖既有结果。
