# 活跃代码与复现入口

| 范围 | 入口 | 说明 |
|---|---|---|
| 路径/身份/只读复用 | src/wcommon.py | 环境变量F10W_OUTPUT_ROOT/VIDEO/THIRD_PARTY/MODELS |
| 原片审计 | scripts/audit_source.py | PyAV完整解码、PTS、逐帧SHA |
| 完整离线PIV | scripts/full_video.py | 冻结历史参数，5配置，每50帧保存可验证chunk |
| 读取大记录 | src/full_records.py | 分块SHA校验后流式读取 |
| 官方模型 | src/motion_models.py | LRAFT/UnLite/CoTracker thin adapters；POT官方函数 |
| 新素材协议 | configs/model_protocol.json | 46/92/46/138，选型与最终测试隔离 |
| 数据 | scripts/generate_model_data.py | 解析运动+真实图块，禁止称物理真值 |
| 统一评价 | src/benchmark.py / scripts/model_benchmark.py | 先开发、验证、校准冻结，再final；已有freeze拒绝重调 |
| 实际定速回放 | scripts/causal_replay.py | 全原片、latest slot1、断帧重置、M4 v2→v1 |
| 历史/真片/官方例 | scripts/extra_diagnostics.py | real/historical/ot/official；历史仅诊断 |
| 图表 | scripts/figures_full.py / analyze_models.py | 保存数据重绘，无预测修饰 |
| 视频 | scripts/full_videos.py / comparison_video.py | 复用PyAV/既有Writer，源PTS、完整解码 |
| 合同测试 | tests/test_contracts.py | 数学、越界、OT失败、Schema、TTL、ZMQ |
| 环境和资产 | scripts/remote_asset_audit.py / audit_server_environment.py | 官方SHA、许可证、锁定依赖、有限I/O探针 |
| 发布快照 | scripts/deploy_snapshot.py | Git archive不可变目录、独立Supervisor标签，无GitHub push |
| 最终核查 | scripts/verify_delivery.py | 全帧、manifest、分割、历史保护、视频/图表 |
| 报告/进度 | scripts/write_reports.py / update_progress.py | 原始结果出表，保留中文历史及diff |

恢复大结果：`python scripts/restore_assets.py fullvideo_delivery.tar.gz`及`python scripts/restore_assets.py models_complete.tar.gz`；脚本先校验SHA，只补缺失文件，不覆盖已有不同内容。

大数据可读取入口：pyarrow.parquet.read_table('results/full_video_piv_per_frame.parquet')；局部向量用full_video_vectors.parquet（dx/dy为区间像素位移，速度需除实际PTS差）。完整候选/图像质量/各方法对象在full_video_chunks，远程映射与SHA见REMOTE_ARTIFACTS.json。原始录像直接使用指定本地文件，模型权重用models/INDEX.json恢复。F10V历史依赖仍在相邻collective_motion_v2，不需要拷贝或改写。
