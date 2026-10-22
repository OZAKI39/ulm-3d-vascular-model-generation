# F10W-Vision

EXPLORATORY · PHYSICS_UNCALIBRATED · training_ready=false

已完成全937帧冻结PIV、三种官方预训练模型+POT/CPD的实际对照、四路全文件定速回放、11组PNG/PDF、3段MP4与审查。结论：暂不替换PIV。

- [26项最终中文报告](FINAL_RECOMMENDATION.md)
- [完整录像报告](FULL_VIDEO_PIV_REPORT.md)
- [模型对照](ZERO_SHOT_MODEL_COMPARISON.md)
- [代码与复现入口](ACTIVE_CODE_INDEX.md)
- [机器Gate](FINAL_GATE.json)
- 本地视频位于`videos/`；两段full_duration均937帧，比较视频明确为64帧selected clip。

本地目录：`/home/lzy/projects/computation_examples/_archives/mevo_vision_f10w_20261017/external_connect/vision_tracking_research/collective_motion_v3`；分支`dev/f10w-pretrained-full-video-20261017`。服务器独立目录`/workspace/f10w_vision_20261017`。没有GitHub push，第三方完整仓库、权重、大缓存和原片不纳入Git；恢复索引见REMOTE_ARTIFACTS.json和models/INDEX.json。

隔离环境复现：本地`.venv/bin/python -m unittest discover -s tests -v`；已下载结果下运行`figures_full.py`、`analyze_models.py`、`full_videos.py`、`comparison_video.py`重绘。所有脚本位于scripts。使用新的F10W_OUTPUT_ROOT运行新实验，不能覆盖原结果。模型推理需服务器CUDA环境；最终锁定依赖在evidence/server_bundle/setup/requirements_final.lock。NumPy1.26使用PyArrow19，不能直接升到26。
