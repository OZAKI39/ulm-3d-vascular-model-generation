# 粒子损伤与力源：v3

本版依据用户最后的更正，顶部使用 **Arial 常规字重，不加粗**。从 v2 渲染器创建副本继续调整，原视频、代码和计算结果均保留。

主要结果：[60 fps MP4](particle_damage_load_source_60fps.mp4)、[20 fps GIF](particle_damage_load_source_20fps.gif)、[结果页面](OPEN_RESULTS.html)、[末态图片](preview.png)。

## 画面变化

- 保持纯黑底、英文、损伤色条、固定视角、等空间单位及真实位移比例。
- 删除右侧 Force 曲线、所有力数值和 N 标记。
- 删除全部粒子受力箭头及力箭头比例图例。
- 主图使用释放出的右侧空间。白色空心菱形表示原配置的局部载荷中心；白色方形轮廓表示固定基底。
- 标题为 `Particle damage & load source`，实际使用本机 `/mnt/c/Windows/Fonts/arial.ttf`，不使用替代字体。字体路径、字重和 SHA256 记录在 `RENDER_MANIFEST.json`。

## 力源标记的定义

空心菱形的位置来自原配置 `streaming.bubble_center_m = [0.0001, 0, -0.0002] m`，即 `[0.1, 0, -0.2] mm`。它表示原合成模型中局部载荷包络的中心，既不是有实际大小的微泡球面，也不是新增的超声换能器或单个粒子上的力。

没有增加新的力学源，没有将菱形大小或亮度解释为力幅值，也没有人为添加源的脉动。原模型中的背景管流及表面牵引继续存在于已有计算结果中；分布式背景载荷不被伪造为另一个点源。

## 流畅度处理与数据边界

原视频虽然编码为 12 fps，但只有 26 个不同状态，每个状态静止显示 0.5 s，因此实际位置更新仅约 2 次/秒。这是主要的阶跃观感来源。

新 MP4 为 **1920×1080、60 fps、751 帧、12.5167 s**。每相邻保存状态之间生成 30 个均匀显示间隔，逐帧渲染粒子位置，不重复插入静止帧。固定基底始终保持精确原坐标，镜头固定，数据不外推。

**平滑只用于显示位置，不能视为新增的求解器轨迹。** 粒子位置在相邻真实状态之间做线性插值；每第 30 帧与原保存位置逐位一致。线性插值无过冲，但不恢复未保存的周期内振动或真实曲线路径。

**损伤不插值。** 每个宏步区间保持前一个已保存损伤值，到下一个原状态时才更新。因此局部颜色仍可能按宏步变化，这不代表视频丢帧，也不改变原模型的失效时刻。原始位置、损伤、断键和碎片结果文件均未修改。

每两个原状态仍映射为 0.5 s 播放时间；原机械代理时间的首间隔为 4 ms、其余为 2 ms。播放时间不是统一比例的物理时间或治疗时间。逐帧映射见 [frame_map.csv](frame_map.csv)，显示坐标及所用原始状态索引见 [display_samples.npz](display_samples.npz)。GIF 为 1280×720、20 fps 的便捷预览，流畅观看以 60 fps MP4 为准。

## 核验与复现

[QC.json](QC.json) 检查原始关键帧一致性、所有粒子保留、固定基底、线性插值公式、损伤不插值、字体身份、原文件哈希、MP4/GIF 完整解码和视频时间戳间隔。只做渲染，不重新求解或修改任何科学模型。

在项目根目录运行，输出目录必须是新目录：

```bash
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$PWD/build_fragmentation/matplotlib" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B scripts/render_damage_source_v3.py --run results/streaming_fragmentation_demo --output visualization/streaming_fragmentation_demo/damage_source_v3_replay
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B scripts/check_damage_source_v3.py --output visualization/streaming_fragmentation_demo/damage_source_v3_replay
```

另一个系统上可用 `--font /path/to/arial.ttf` 指定已安装的 Arial 字体。渲染器、上一版脚本副本和输入哈希保存在 `provenance/`，未捆绑分发系统字体。
