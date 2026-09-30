# 粒子损伤与局部声流式牵引：v4

本版将用户指定的**局部声流式牵引**直接画成表面箭头，移除上一版位于血栓上方的空间包络中心菱形。该箭头表达表面单位面积上的牵引，单位为 Pa；方向、大小均来自现有模型的局部分量。

结果入口：[播放页面](OPEN_RESULTS.html)、[60 fps MP4](particle_damage_streaming_traction_60fps.mp4)、[20 fps GIF](particle_damage_streaming_traction_20fps.gif)、[初态预览](preview.png)。

## 图中元素

| 元素 | 含义 |
| --- | --- |
| 圆形粒子及颜色 | 原计算的粒子位置与损伤；损伤色条范围 0–1 |
| 白色方形轮廓 | 固定基底粒子，内部颜色仍表示损伤 |
| 浅绿色箭头 | 当前仍受载表面上的局部声流式牵引方向和强度 |
| X、Y、Z 网格 | 原模型空间坐标，单位 mm，位移比例为 1 |

保持纯黑背景、Arial 常规字重、英文标题、放大的单主图、固定视角和损伤色条。无中心标记、微泡球体、N 标记或右侧 Force 面板。原版本、科学计算配置及结果全部保留；本次没有重跑动力学模拟。

## 牵引方向与数据分量

数据读取 `results/streaming_fragmentation_demo`，实际解析到 `runs/fragmentation_pilot_003`，包含 26 个保存状态和 360 个粒子。

原 `pd_clot/streaming.py` 的 `AnalyticTestStreaming.traction` 给出本算例的局部分量：

```text
t_local = A(x) · g(t) · 100 Pa · [e_x − (e_x · n)n]
A(x)    = exp[−|x − c|² / (2L²)]
g(t)    = 1 + 0.9 sin(2π · 500 Hz · t)
L       = 0.35 mm
```

`n` 取每个保存状态 VTP 中的实际表面法向，`x` 取同一状态的粒子中心。配置的法向附加牵引为 0 Pa，因此已保存状态的局部牵引是切向的。它以 **+X 管轴方向在局部表面上的投影**为方向：平坦上表面的箭头主要沿 +X，倾斜及变形表面会改变箭头方向。上方的包络中心仅影响空间强弱分布，不表示从上往下施力。本版不再绘制该中心。

读取原配置后仅在内存副本中关闭 `include_pipe_traction`，用于提取局部分量；磁盘配置没有变化。逐个原状态核对：

```text
局部声流式牵引 + 背景管流压力/黏性牵引 = 原 VTP 保存的总表面牵引
```

26 个状态的最大重构误差为 **0 Pa**；局部牵引与法向的点积绝对值最大约 **2.49×10⁻¹⁴ Pa**；保存状态的局部牵引峰值约 **81.794 Pa**。完整数值记录于 `RENDER_MANIFEST.json` 和 `QC.json`。

箭头不包含背景管流压力/剪切、内部弹性力或脱落碎片的速度松弛项。粒子的已有运动仍来自原完整模型，不能归因于画出的单一分量。本模型仍是未标定的解析声流式牵引验证模型，不是已求解、实验验证的超声声流场。

## 箭头采样与比例

- 受载资格沿用原求解器规则：已暴露、仍连接基底、有有效表面法向且非固定粒子。
- 在曾满足资格的粒子中，按初始空间覆盖选取 32 个固定材料 ID，避免每帧重新随机采样。某 ID 脱离受载资格时，对应箭头消失；不会给已脱落碎片另加表面载荷。
- 所有帧共用线性比例 **0.005 mm/Pa**，即 50 Pa 对应 0.25 mm 的三维箭头长度。不逐帧归一化、不截断幅值；屏幕长度还受三维视角投影影响。
- 为让箭头在表面附近清晰可见，箭尾沿显示法向外移半个粒子间距，即 **0.0625 mm**。这一偏移仅用于显示；牵引仍在原粒子中心计算，没有改变力学作用点。
- 图例用于识别箭头类别，其印刷长度不作为空间比例尺。

## 60 fps 显示与时间边界

MP4 为 **1920×1080、60 fps、751 帧、12.5167 s**；GIF 为 **1280×720、20 fps** 便捷预览。相邻原状态映射为 0.5 s 播放时间，播放时间不是统一比例的物理时间或治疗时间。

位置在相邻原状态之间线性插值，仅用于平滑显示；每第 30 帧与原保存位置精确一致，固定粒子始终保持原坐标。损伤不插值，保持前一个原状态的值直到下一保存状态。

牵引矢量只对**两个相邻原状态都受载**的同一材料粒子做显示插值。否则保持前一个保存值；受载掩码也保持前一状态，直到原状态端点才更新。显示用法向可插值并归一化，仅用于放置箭尾。每第 30 帧的牵引矢量与提取的原状态分量精确一致。

插值矢量不是新的求解器采样；法向发生变化时，中间显示矢量也不保证严格满足瞬时切向条件。26 个原端点满足切向条件。没有用插值重建未保存的 500 Hz 周期内振动，也没有人为添加脉动。颜色和受载资格仍会在原宏步更新，这是离散记录的变化，不是视频缺帧。

## 文件及核验

- `traction_samples.npz`：原状态局部牵引、背景管流牵引、原总牵引、法向、受载资格、材料 ID、代理时间和坐标。
- `display_samples.npz`、`frame_map.csv`：逐帧位置、损伤、牵引、箭头及原状态映射。
- `keyframes/`：全部 26 个原状态对应的图像。
- `QC.json`：独立公式核验、分量重构、切向性、原端点一致性、插值规则、字体、文件保护以及 MP4/GIF 解码检查。
- `provenance/`：新渲染器、校验脚本、上一版渲染器副本、输入哈希及渲染日志。
- `DELIVERY_SHA256.json`：交付目录文件哈希，不包含清单自身。

复现时从项目根目录运行，输出目录必须是新目录：

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD" MPLCONFIGDIR="$PWD/build_fragmentation/matplotlib" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B scripts/render_streaming_traction_v4.py --run results/streaming_fragmentation_demo --output visualization/streaming_fragmentation_demo/streaming_traction_v4_replay
PYTHONDONTWRITEBYTECODE=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B scripts/check_streaming_traction_v4.py --output visualization/streaming_fragmentation_demo/streaming_traction_v4_replay
```

字体读取本机 `/mnt/c/Windows/Fonts/arial.ttf`，可用 `--font` 指向另一台机器已安装的 Arial。未分发系统字体文件。
