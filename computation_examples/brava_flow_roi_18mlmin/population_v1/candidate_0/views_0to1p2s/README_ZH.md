# Candidate 0：0–1.2 s、四个观察角度

已生成 4 段固定视角动画，包括当前视角。黑色背景、蓝色微泡、透明灰色血管、无包围框的流场坐标轴、英文标题和端口标注沿用上一版。所有角度均自动放大到当前预留显示范围内；当前 240° 版本的相机与上一版一致。

| 相对方位角 | 绝对方位角 | 仰角 | 说明 |
|---|---|---|---|
| 240° | 16.218875° | 25° | 保留的当前视角 |
| 330° | 106.218875° | 25° | 新增视角 |
| 60° | 196.218875° | 25° | 新增视角 |
| 150° | 286.218875° | 25° | 新增视角 |

右上角 `Azimuth offset` 表示绕候选姿态坐标系 Z 轴、相对原始基准视角的方位角增量。基准绝对方位角为 136.218875235°；绝对方位角从 +X 朝 +Y 测量。`Elevation` 为相对于 XY 平面的仰角。镜头在每段动画中固定。

## 时间和显示筛选

- 物理时间严格为 **0–1.2 s**，每段 **262 帧、24 fps、10.916667 s**。沿用此前放慢 80% 的播放速度，帧数按整数取整；`dt = 0.5 ms` 是原始模拟时间步。
- 全程隐藏此前 2 条求解失败轨迹和 58 条数值风险轨迹。
- 对每帧时刻 t，按真实三维速度模长在 `[t−0.05, t+0.05] s` 内的时间加权平均值筛选；小于 2 mm/s 时隐藏。在实际轨迹出生或终止处仅使用存在的时间段。
- 视频结束于 1.2 s，但其末帧速度窗口仍为 1.15–1.25 s：计算使用已存在的后续轨迹，未在视频边界截断窗口，也未新增求解。
- 4 个视角的每帧位置、时间、平均速度、可见性掩码均完全一致；末帧各显示 **1027** 个微泡，整个片段共有 **1442** 条轨迹至少显示过一次。
- 微泡跨过速度阈值时会出现或隐藏，这是显示筛选行为；不改变运动轨迹，也不代表注入、消失或求解状态改变。

## 文件

每个 `azXXX_el25p0/` 子目录包含：

- `animations/microbubble_azXXX_el25p0_0to1p2s.mp4`：动画。
- `figures/snapshot_end.png`、`video_end.png`：1.2 s 原始渲染及编码后抽帧。
- `data/render_manifest.json`、`layout.json`：来源、渲染参数、相机及端口标签。
- `data/frame_summary.csv`：每帧活跃、隐藏及显示数量。
- `data/animation_samples.npz`：真实插值位置、时间、速度、粒子 ID 和显示掩码。

总目录提供 `OPEN_RESULTS.html`、`CONTACT_SHEET.png`、`data/view_summary.csv` 和 `data/QC.json`。

## 已执行检查

4 段视频均完成全片解码；分辨率、帧率、帧数、时长正确。相机角度与右上角标注一致；4 组逐帧数据完全相同。利用原始已接受轨迹段的位移/时长，对 5 个时刻的 3630 组有效粒子窗口独立核算，平均速度最大差为 7.08×10⁻¹⁵ m/s，位置差为 0。当前视角相机与上一版的最大数值差为 8.53×10⁻¹⁴。源轨迹、流场、既有动画、静态图和 candidate 15 经哈希校验保持不变。

显示筛选仅用于展示，不等于已经证明保留轨迹的物理准确性。本轮没有启动 CFD 或轨迹重算。

## 复现

服务器使用 RTX 4090 / VTK EGL 渲染、CPU 编码：

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 env/bin/python -u -B scripts/render_population_four_views.py
OPENBLAS_NUM_THREADS=2 env/bin/python -u -B scripts/check_population_four_views.py
```

本地根目录：`/home/lzy/projects/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s/`。
服务器镜像：`/workspace/brava_flow_roi_18mlmin_20260928/population_v1/candidate_0/views_0to1p2s/`。
