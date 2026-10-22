# Candidate 0：第 5 视角微泡动画

本次仅对已有轨迹重新绘图，没有启动 CFD 或轨迹求解，没有改变原始数据。

- 固定 `view_05_az240_el25.0_t6p0s.png` 的观察方向（相对方位角 +240°、仰角 25°），动画中不旋转相机。
- 去除白色包围框与框体底色。直接复用 `visualization/style.py` 的 `CoordinateGrid` 和 `AnnotationPlan` 绘制坐标参考平面、坐标轴、毫米刻度；固定相机投影和文字可用区域按本界面适配。
- 血管屏幕尺寸相对该参考静态图放大 **41.79%**。完整血管保持在画面内，保留入口及三个出口的真实端口标注。
- 黑色背景、蓝色微泡、透明灰色血管、英文四行标题保留；无底部说明文字。
- 显示物理时间 **0–6 s**，24 fps，共 **1306 帧**，视频时长 **54.417 s**。沿用此前放慢 80% 的播放速率；因整帧取整，物理时间采样间隔有极小差别。标题中的 `dt = 0.5 ms` 指原始模拟时间步，不是视频帧间隔。

## 显示筛选

此前的 2 条求解失败轨迹和 58 条数值风险轨迹全程隐藏。对每个视频时刻 t，读取已接受轨迹段的真实三维速度，在 `[t−0.05, t+0.05] s` 内按时间加权平均速度模长，隐藏平均速度低于 2 mm/s 的微泡。出生或终止附近仅使用实际存在的时间段，不补零、不外推。该规则由用户确认逐帧执行，因此微泡跨过阈值时会显示或隐藏，这不代表新的注入、消失或改变运动状态。

6 s 末帧显示 **124** 个微泡，与此前静态图的 ID 完全一致。全过程共有 1505 条轨迹至少显示过一帧。筛选只用于展示，不能据此认定保留轨迹均已通过物理正确性验证。

## 文件与复现

- `animations/microbubble_view05_physical_time.mp4`：动画。
- `figures/view05_t6p0s.png`：6 s 原始渲染预览；`video_start/middle/end.png` 为编码后抽帧。
- `data/frame_summary.csv`：每帧物理时间、存活数、原有排除数、额外低速排除数及显示数。
- `data/animation_samples.npz`：真实插值位置、平均速度、有效窗口时长、逐帧显示掩码及粒子 ID。
- `data/layout.json`、`coordinate_annotations.json`：相机、坐标轴、端口坐标及标签落点。
- `data/render_manifest.json`、`QC.json`：来源、参数和检查结果。
- `data/track_hashes.json`、`protected_files_sha256.json`：源轨迹及既有成果的校验值。

服务器使用 RTX 4090 / VTK EGL 绘图，CPU 编码；脚本位于 `population_v1/scripts/`：

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 env/bin/python -u -B scripts/render_view05_flow_axes.py --end-time 6
OPENBLAS_NUM_THREADS=2 env/bin/python -u -B scripts/check_view05_animation.py
```

本地新成果路径：`/home/lzy/projects/brava_flow_roi_18mlmin/population_v1/candidate_0/view05_flow_axes/`。
服务器镜像：`/workspace/brava_flow_roi_18mlmin_20260928/population_v1/candidate_0/view05_flow_axes/`。
既有动画、静态图、审计记录及 candidate 15 均保留。
