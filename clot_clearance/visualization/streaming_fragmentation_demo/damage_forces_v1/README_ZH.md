# 粒子损伤与外力：新可视化版本

本目录使用既有 `runs/fragmentation_pilot_003` 的 26 个保存状态重新绘图，未重新积分轨迹，未修改原计算结果或原两部视频。

打开 [OPEN_RESULTS.html](OPEN_RESULTS.html) 或 [MP4](particle_damage_external_forces.mp4)。另有 [GIF](particle_damage_external_forces.gif)、[末态预览](preview.png)、[外力历史 CSV](force_history.csv)、[逐粒子力数据 NPZ](force_samples.npz) 和 [核验记录](QC.json)。

## 画面含义

- 单一三维视图，所有粒子均按损伤 0～1 着色；方形白色轮廓表示固定基底。保留真实坐标与位移，固定镜头。
- 不绘制微泡标记、碎片编号/分类面板、活动键/断键连线或清除平面。
- 青色箭头表示抽样材料粒子所受的**外部流体力**，箭头尾端位于粒子的真实坐标。方向来自力矢量，长度正比于力的大小。
- 全动画使用同一线性比例：**1 μN 对应坐标中的 0.35 mm 箭头长度**。这是力矢量的显示比例，不是粒子位移。没有逐帧归一化、对数缩放或截断大箭头。
- 箭头固定选取 54 个非基底材料粒子，按初始空间位置覆盖选点；严格零力的箭头不画。粒子本身不抽样，360 个全部显示。右侧统计与曲线使用全部粒子，不受箭头抽样影响。
- 右侧 Fx、Fy、Fz 是全部粒子外部流体力的矢量和，单位 μN；虚线和圆点标记当前状态。正负号对应坐标轴方向。`Peak particle force` 是当前所有粒子中的最大外力模长。
- `Force scale: 1 μN` 是固定参考箭头，不是额外施加的力，也不代表微泡。

## 外力的数据定义

附着、暴露且可移动的材料使用保存的表面牵引 `surface_traction_Pa`，乘以计算中相同的受载面积 `h² × area_factor`。脱离的可移动材料使用保存的流体速度与粒子速度：

```
F_surface = traction * h² * area_factor
F_relaxation = rho_s * V * (u_fluid - v_particle) / tau_h
```

附着材料的松弛力为零，脱离材料的表面牵引为零；固定基底的上述外力均为零。没有将内部 PD 键力、数值阻尼或固定约束反力混入“外部流体力”。因而图中的合力不是用于计算整块材料加速度的全部力之和。

每帧对应宏步损伤/连通性更新后的保存状态，显示该状态下的瞬时模型外力。脱离粒子的松弛力由原模型公式在保存状态上计算，不使用宏步冲量代替瞬时力。保存的速度场及表面牵引均用身份哈希匹配的原提供器只读重算核对，误差为零。

这些状态均处于规定的 500 Hz 慢载荷同一相位。动画可反映几何、暴露面和附着状态变化带来的外力变化，**不能恢复未保存的周期内振荡**。曲线连接保存点仅用于阅读，没有补造中间轨迹、方向或损伤值。

## 播放与核验

视频 1920×1080、12 fps、156 帧、13 秒；每个保存状态重复 6 帧。播放时间不是实际物理时间。已完整解码 MP4 与 GIF，并检查初始、首次脱离及末态画面；原始位置、损伤、力求和、统一箭头比例及原文件哈希核验见 `QC.json`。本版本只改展示，原模型的未标定、低维支撑等科学限制继续适用。

从项目根目录重新生成时使用新的输出目录：

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD" MPLCONFIGDIR="$PWD/build_fragmentation/matplotlib" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B scripts/render_damage_forces_v1.py --run results/streaming_fragmentation_demo --output visualization/streaming_fragmentation_demo/damage_forces_replay
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B scripts/check_damage_forces_v1.py --output visualization/streaming_fragmentation_demo/damage_forces_replay
```

渲染器基于原 `scripts/render_fragmentation.py` 另存开发；原脚本副本、本版脚本及输入哈希保存在 `provenance/`。
