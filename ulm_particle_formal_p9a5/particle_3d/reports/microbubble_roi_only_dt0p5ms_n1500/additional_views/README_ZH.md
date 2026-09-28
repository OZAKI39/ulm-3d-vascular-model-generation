# 当前微泡动画：三个额外视角

复用 ROI-only-balanced-pressure-v1 的现有 1500 条有限尺寸微泡轨迹，dt = 0.5 ms；没有重新积分轨迹或运行 CFD。原动画、原绘图代码、原报告及轨迹保持不变。

[打开三个视角的预览页面](OPEN_RESULTS.html) · [原动画与数据](../OPEN_RESULTS.html)

## 视角及样式

| 视角 | 方位角 | 仰角 | MP4 |
|---|---:|---:|---|
| Side view | -133.781125° | 6.294494° | `animations/microbubble_age_aligned_side.mp4` |
| Reverse view | -43.781125° | -6.294494° | `animations/microbubble_age_aligned_reverse.mp4` |
| Elevated oblique view | -178.781125° | 45.000000° | `animations/microbubble_age_aligned_oblique.mp4` |

方位角从 +X 朝 +Y 计算，仰角从 XY 平面计算；相机方向向量指向从目标到相机。三种视角均固定 Z 轴朝上、采用正交投影。相对原视角，侧向绕 Z 转 90°，反向取原方向的相反方向，斜俯视绕 Z 转 45°并将仰角设为 45°。

保持原来的英文标题、图例、说明文字、字体、背景、血管透明度和颜色。投影比例统一为原默认尺寸的 1.30 倍，没有在上一轮 1.30 倍基础上再次放大；微泡显示点直径为 9.1 像素。不同视角的几何投影外观尺寸会因朝向而改变，物理显示比例相同。正交相机 parallel_scale 均为 47.042127755972054 μm。

每段 MP4 均为 1920×1080、24 fps、288 帧、12 秒，轨迹年龄范围为 0–146.029479 ms。三个视角的每一帧对应相同年龄和相同粒子位置。按各轨迹起始年龄对齐回放，不表示同时注入；已完成出口穿越的微泡在其计算结束后消失，接触支持静止微泡保持最后位置。

图例表示整批轨迹终态计数：O1=787、O2=61、O3=467、接触支持静止=185；不是每一帧中可见的微泡数。

## 渲染及验证

- GPU：服务器 NVIDIA GeForce RTX 4090，NVIDIA EGL / PyVista 渲染；编码为 CPU libx264。
- `render_views.py` 直接导入原 `../scripts/render_results.py` 的样式常量和放大函数，仅追加相机定义；坐标插值和出口后显示规则与原代码相同。
- 渲染前后核验 1500 条轨迹的 3000 个直接使用的轨迹/统计文件；45 个原交付文件（包括原 MP4、原代码和原页面）保持哈希不变。
- `verify_views.py` 实际完整解码三段视频，检查每段 288 帧、尺寸、时长、统一比例、不同相机和文件哈希。结果见 `VALIDATION.json`。
- 已人工查看每个视角的第 24、96、160、240 帧拼图：英文标题、图例和脚注可读，血管位于内容区内。拼图见 `figures/*_frames.png`。
- 第 96 帧预览另存 PNG（300 dpi 元数据）和 PDF。它们是实际渲染帧，不是另一套科学计算。
- `RENDER_MANIFEST.json` 保存相机坐标、投影边界、时间轴哈希、后端信息和视频哈希；`logs/` 保存本次执行记录。

## 路径与复现

本地目录：
`/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/additional_views`

服务器目录：
`/workspace/microbubble_roi_only_dt0p5ms_n1500_20260928/additional_views`

实际执行入口由本目录 `supervisord.conf` 指定；仅包含 `render_views`，autostart=false、autorestart=false。日志记录正常退出码 0；本次专用 supervisor 在渲染完成后已正常关闭。

本地验证命令：
```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -B /home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/additional_views/verify_views.py
```

如需重现，应在同一批次根目录下建立一个新的兄弟输出目录，复制 `render_views.py` 及 `audit/SOURCE_LOCK.json`，建立空的 `animations/`、`figures/` 子目录。该脚本按自身父目录定位批次数据，并拒绝覆盖已有 MP4。在服务器使用 `/root/particle8_2_runs/env/bin/python -B <新目录>/render_views.py`；长任务仍通过专用 supervisor 执行。输入流场路径取父目录 `config.json`。无需启动粒子计算或 CFD。
