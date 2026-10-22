# F10T-Vision-2 活跃入口

| 入口 | 用途 |
|---|---|
| start_annotation.sh / scripts/annotate.py | 已修复的持续窗口、napari原生绘图、类别与版本化保存 |
| scripts/sam_worker.py | 可取消原SAM后台辅助；不含SAM2 |
| src/annotation_store.py | 五类掩码、确认关口、SHA、Cellpose导出 |
| scripts/prepare_candidates.py / acquire_live_candidates.py | 预登记录像候选、独立直播采集；拒绝覆盖 |
| scripts/review_labels.py | 核查人工标签；完成后冻结具体修订和导出 |
| src/annotation_snapshot.py | 撤回误标整帧审核而保留像素；带哈希的图像/标签快照、仅从快照导出、数据就绪盘点 |
| scripts/correct_partial_reviews_20261016.py | 本次用户明确的第 1–5 帧审核范围修正；带窗口锁和备份，已执行后拒绝重复覆盖 |
| scripts/windows_frame_server.py / src/live_bridge.py | 原Windows接收器的有界原始图像IPC |
| src/vision_engine.py | 两条候选检测路线、Norfair、图像矩、稀疏LK |
| src/runtime.py / scripts/run_live.py | 最新帧、过期拒绝、mock、资源与真实直播日志 |
| scripts/compare_routes.py | 固定同帧同ROI成对测试、严格因果离线短回放 |
| scripts/render_delivery.py | 四组真实PNG/PDF和两段MP4；--figures-only只更新图 |
| tests/test_contracts.py | 标注I/O、冲突、导出、形态、时间/session、mock契约 |
| scripts/write_delivery_reports.py | 从真实日志生成交付文档与分层关口 |
| scripts/verify_delivery.py | 最终文件身份、视频、标签保护及交付检查 |

本地分支dev/f10t-vision2-20261016；基线9c22d207be5318f9d45b014caaf52bd04d345328。原相机 / AFG工程不修改，历史未清理。旧版GUI仅保留evidence/annotation_repair_20261016作故障证据，不再是启动入口。

2026-10-16 图层 / 局部 SAM 修复：`src/sam_prompt_geometry.py` 负责原图坐标和候选范围检查；`scripts/check_annotation_v3.py` 是当前独立 Qt 回归入口；`tests/test_sam_prompt_geometry.py` 检查坐标和候选保护。有限训练提示对照入口为 `scripts/benchmark_sam_annotation_v3.py` 与 `scripts/check_sam_point16_v3.py`，现有输出目录不可覆盖。修复证据及当前源码哈希见 `evidence/annotation_repair_v3_20261016/`；此前交付清单保留原快照含义。
