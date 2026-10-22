# 代码与历史审计

实际从F10V完整交付80440a8fc5f15716a2a5835b8eb2d509e20a30c1创建隔离worktree和dev/f10w-pretrained-full-video-20261017，未用较旧F10U同步覆盖最新工作。原F10V工作树启动时干净；服务器物理仓库66d4ceb未改。适用目录未发现额外AGENTS；服务器运行前已读/etc/vast-agents-guide.md并使用独立Supervisor标签。

F10U CandidateDetector/distribution/OpenPIV/pycpd/affine/state直接只读复用；F10V只复用图块splat、M4 v2、资源和视频编解码辅助。新代码只提供官方模型适配、数据协议、运行编排、结果审核与可视化；没有重写RAFT、CoTracker、Sinkhorn或视频编解码器。third_party被Git排除。冻结执行快照及代码SHA见各COMPLETE/EXECUTION文件。

原片、人工标注、F10U/F10V、相机/AFG和物理历史共5214个受保护文件逐一复核，结果见evidence/final_integrity.json。科研总进度文档仅在全部实测后按用户要求备份并添加F10W段落，未删历史结论。没有GitHub push。
