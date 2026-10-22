# Stage 3 三维微泡运动动画

本次只读取已接受的 Stage 3 圆管轨迹并生成显示帧，不重新运行声场或粒子求解器。

## 打开结果

本地输出目录：`results/stage3/visualization/reference_style_v1/`。

- [播放页面](../results/stage3/visualization/reference_style_v1/OPEN_RESULTS.html)
- [C3 主动画](../results/stage3/visualization/reference_style_v1/stage3_C3_motion_3d_60fps.mp4)
- [C1–C4 同步对照](../results/stage3/visualization/reference_style_v1/stage3_C1_C4_comparison_3d_60fps.mp4)

两段均为 1920 × 1080、60 fps、1321 帧，时长约 22.017 秒。
正式物理时间 0–2 ms 映射为 20 秒运动过程，开头和结尾另有停留。

## 样式与数据含义

参考视频为
`/home/lzy/projects/clot_clearance/visualization/streaming_fragmentation_demo/streaming_traction_v4/particle_damage_streaming_traction_60fps.mp4`。

沿用黑底、Arial 常规字重、大号英文标注、深蓝色三维坐标网格、无外围面板边框和固定正交相机
（仰角 24°、方位角 −68°）。微泡使用蓝色，青绿色曲线表示最近 0.25 ms 的轨迹，不是力矢量。
没有把参考视频中的损伤量、流牵引力或固定基底移植到微泡结果中。

这是泡群局部特写，单位 µm；圆管全长 20 mm、半径 1.5 mm，管壁不在特写范围内。
所有场景使用相同相机、空间界限和物理比例。所有 64 个粒子及其身份保留。
粒子中心使用原始世界坐标，位移缩放为 1。
圆形投影按真实平衡半径 1.78 µm 换算屏幕尺寸，边缘细线用于辨识；不显示未保存的径向脉动。

| 场景 | 物理开关 | PZT 相位 |
| --- | --- | --- |
| C1 | 仅 Poiseuille 流动 | 全部关闭 |
| C2 | 流动 + 主辐射力 | 0°, 0°, 0°, 0° |
| C3 | 流动 + 主辐射力 + 次辐射力 | 0°, 0°, 0°, 0° |
| C4 | 同 C3 | 0°, 90°, 0°, 0° |

## 来源及显示插值

正式运行父目录为 `results/stage3/cylinder/20261002T121713774073Z/`。
具体目录从 `results/stage3/cylinder/STATE.json` 读取，子目录分别以 C1、C2、C3、C4 结尾。
直接读取 `particles/trajectories.h5` 的 `time_s`、`positions_m`、`radii_m`、`particle_id`。

逐帧时钟采用统一物理时间，在线性相邻保存状态间插值；没有时间外推、空间平滑、接触投影修正或位移放大。
代码验证在原始保存时刻能逐位重建原始坐标，并验证起终点一致。尾迹包含窗口内的每一个原始保存顶点。

线性显示插值可能使两个接触粒子的间隙比保存状态略小：本次 C3 最小约 −10.7 nm，C4 约 −4.34 nm。
它是显示插值造成的弦线误差，不能作为新增计算结果使用；未通过移点修补外观。
确切数值及来源哈希保存在 `RENDER_MANIFEST.json`。
Stage 3 原有近接触、准稳态平移及历史力敏感性限制继续适用，动画不增加物理验证证据。

## 复现与检查

渲染使用现有共享 Python 环境，不需要启动服务器计算任务：

```bash
cd /home/lzy/projects/computation_examples/mevo_acoustic_microbubbles
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /home/lzy/projects/computation_examples/_shared/python_environment/bin/python \
  -m scripts.stage3.render_motion \
  --output results/stage3/visualization/reference_style_replay
```

先检查关键帧可追加 `--preview-only`。已有 MP4 的目录不会被覆盖。
Arial 默认来自 `/mnt/c/Windows/Fonts/arial.ttf`，可通过 `--font` 指定相同字体文件。
依赖已有 NumPy、h5py、Matplotlib、Pillow、imageio 与 imageio-ffmpeg；不导入生产求解器。
新环境可通过项目的 `visualization` 可选依赖安装视频编码组件；本次直接使用已有环境。

输出包括：

- 两段 MP4、10 张关键帧、预览图和三时刻联图；
- `display_samples.npz`：显示坐标、粒子身份、半径、视频及物理时钟；
- `frame_map.csv`：每组每帧对应的左右原始状态和插值系数；
- `provenance/`：渲染器副本及输入文件 SHA-256；
- `RENDER_MANIFEST.json`：相机、颜色、屏幕边界、视频规格和输入未变更检查；
- `DELIVERY_SHA256.json`：输出文件哈希。

运行中检查轨迹完成标记、时间单调性、坐标有限性、初始云一致性、原始时刻重建、画面内粒子范围，
并在编码后核对实际帧数、分辨率、帧率及首中末帧解码。渲染前后检查所有直接输入的哈希不变。
