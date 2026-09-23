# 全血管压力场与 WSS：放大旋转可视化

打开 `OPEN_RESULTS.html` 可播放两段动画并下载 4K 图片。

- 与最新版全血管流线采用相同取景，固定平行投影尺度为 45.31911463464213 µm。相对旧压力/WSS 动画，单位物理长度的显示比例放大 33.436%。
- Z 轴始终竖直，绕平行于原始 Z 方向、通过固定视图中心的轴旋转。坐标轴、物理坐标刻度、网格和 µm 单位保留。
- Inlet / Outlet 01–03 在血管旁侧显示，以箭头指向对应端口。文字没有黑色底框；换位及显隐采用约 0.875 秒渐变，包含动画末尾与开头的衔接。
- 压力顶部为 `Full vessel | Pressure field`；底部仅为 `Steady FEM pressure on the vessel surface`。
- WSS 顶部为 `Full vessel | Wall shear stress`；底部仅为 `Derived from FEM velocity gradient`。
- 颜色条只保留标题和数字。压力为 `Pressure (Pa)`，范围 −5～2500 Pa；WSS 为 `WSS (Pa)`，范围 0～55 Pa。色表、范围与旧压力/WSS 图相同。
- 两段动画均为 1920 × 1080、24 fps、432 帧、18 秒；18 秒为稳态场的展示旋转时长。4K 图片为 3840 × 2160，选用单一清晰标签位置，便于放入 PPT。

压力沿用完整血管表面的 `Pressure_Pa`；WSS 沿用物理血管壁上的 `WSS_display_Pa`。后者由 FEM 速度梯度计算壁面切向黏性应力，再按三角形面积加权得到节点显示值。未重新求解、修改数值、裁剪数值或平滑几何；曲面法向光照仅影响显示。

独立验证确认：压力与原 FEM 节点值完全一致；WSS 面积加权显示值重构误差为 0；入口出口文字与血管投影无重叠；所有帧完整取景、Z 轴竖直、透明度平滑变化。颜色条附加文字清除通过逐帧像素检查。数值、相机与文字记录见 `MEDIA_VALIDATION.json`、`INDEPENDENT_VALIDATION.json` 和两个视图的 `*_camera.json`、`*_annotations.json`。

完整代码在 `../../render_surface_fields.py`，公共相机与文字代码在 `../../render_visualization.py`。在工程目录运行 `./run.sh --surface-fields` 可重现；运行 `./run.sh --validate-surface-fields` 可独立验证。
