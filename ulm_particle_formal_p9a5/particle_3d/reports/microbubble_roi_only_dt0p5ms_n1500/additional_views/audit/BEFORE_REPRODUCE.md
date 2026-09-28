# 本轮复现与数据位置

服务器：`vast4090`，独立目录：
`/workspace/microbubble_roi_only_dt0p5ms_n1500_20260928`。

本地目录：
`/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500`。

本轮使用的新流场：`fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12`。
当前请求授权该流场进入新微泡批次；此前 CFD manifest 中“未推广到粒子生产”是上一阶段的历史状态，不修改该历史记录。

## 原始依赖

- 服务器原有源码与几何依赖：`/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms`，本轮只读。
- 新流场：`/workspace/flow_roi_only_balance_20260928/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/frozen_flow/steady_flow_mean_2p0_mmps_A_ROI_only_balanced_pressure.vtu`。
- 114 个受保护源码哈希：`data/preparation.json`；本輪不覆盖其源码。
- 配置：`config.json` 中保存服务器绝对路径。在另一台机器复现时，须映射这些依赖到实际副本路径，且流场、网格及源码验证必须继续通过。

## 命令（在新建目录中复现，不能覆盖已有完成记录）

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/root/particle8_2_runs/env/bin/python scripts/campaign.py prepare
/root/particle8_2_runs/env/bin/python scripts/campaign.py pilot
/venv/main/bin/python scripts/gpu_validate.py mesh
/root/particle8_2_runs/env/bin/python scripts/test_native.py
/root/particle8_2_runs/env/bin/python scripts/verify_native_bytes.py
/root/particle8_2_runs/env/bin/python scripts/production_entry.py
/venv/main/bin/python scripts/gpu_validate.py tracks
/root/particle8_2_runs/env/bin/python scripts/summarize.py
/root/particle8_2_runs/env/bin/python scripts/finish_delivery.py
```

生产前必须通过本机实际几何查询与完整轨迹逐位对照；编译核针对本机 NumPy 的浮点累加顺序，不承诺跨平台逐位相同。未通过时不能直接绕过门槛。

当前运行由独立 supervisor 管理，`autorestart=false`；`production` 仅生成当前 1500 条有限尺寸微泡轨迹，不启动 CFD。`delivery` 只等待本批结果并渲染；`resources` 只读采样资源；本地 `collect_live.py` 增量回传，不启动计算。

原始轨迹采用安全判据要求的自适应子步。`data/trajectories_dt0p5ms.npz` 另存按 0.5 ms 采样的位置；原始接受步、事件和审核日志全部保留。

`finish_delivery.py` 在最终汇总存在后调用渲染并写入交付哈希清单。完整回传本地后，以含 NumPy、Pillow 和 imageio-ffmpeg 的 Python 执行 `scripts/verify_collected.py`；它逐条核验原始记录、来源身份、交付哈希及视频解码，并生成 `LOCAL_VERIFICATION.json`。

## 动画显示更新

2026-09-28 按用户要求将动画、预览页及配套图件文字改为英文。动画主体用正交相机放大至 1.30 倍，微泡标记由 7 像素变为 9.1 像素；保留配色、背景、视角方向、回放时间及全部轨迹数据。相机平移仅用于避开标题和页脚。实际投影尺寸比和画面边界见 `data/render_manifest.json`，对比核验见 `logs/presentation_english_zoom130_verification.json`。旧展示文件及校验清单保存于 `logs/presentation_before_english_zoom130/`。本次只重新启动独立 supervisor 中的 `delivery` 渲染任务。
