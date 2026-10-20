# V3 红色血栓与载荷场附加可视化

只读取 `clot_fragment_surface_coupling_v3/20261018T170400Z` 已接受数据，未启动动力学。原有视频、正式源码冻结和原始数据保留。

结果目录：`results/clot_fragment_surface_coupling_v3/red_z_vertical_20261019/`。

- `VIDEO_V3_RED_Z_VERTICAL.mp4`：原对比视频右侧的新 V3 数据单独显示，保留相同实际帧时序；1920×1080、24 fps。红色材料点和真实存活键，固定全景相机、Z 轴竖直、位移 ×1。仅保留时间、坐标、简短图例和 2 µN 力标尺。
- `FIG_V3_PEAK_LOAD_FIELD.png/.pdf`：2.4 ms 名义驱动峰值平台末端，实际表面牵引场。
- `FIG_V3_L2_FRAGMENT_LOAD_FIELD.png/.pdf`：2.650125 ms 的真实 46 节点游离分量，101 条内部存活键、4 个有效 L2 面。

静态图为 3200×1800 PNG；PDF 内嵌相同高分辨率栅格图，非矢量几何。彩色面显示保存的实际切向牵引模值（Pa）；蓝色箭头显示实际单面合力（N），世界长度比例为 2 µN → 200 µm。全景仅抽取空间分散的真实箭头以保持可读性，未改变载荷值；全部有效面数据在对应 `_DATA.csv`。未绘制的面或箭头不代表零物理阻力。

保存的 `surface_traction` 已包含开放权重，故实际面力＝该牵引向量×几何面积。CSV 的 `effective_area_m2` 是几何面积×开放权重，不能再与该已加权牵引相乘而重复施加权重。色值不是未加权的 β×相对滑移，也不是名义 1000 Pa 的均匀填色。

视频不再用多色区分碎块，仅改变颜色表达；碎块、掉点和破裂仍来自同一真实原始材料图。全景绘制存活短键，静态特写绘制该分量全部存活内部键；不绘制失效键轨迹，不构造平滑外壳。表面依然是原方法的部分近似求积，不是新增流场求解。

在原隔离环境内运行：

```bash
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python analysis/clot_fragment_surface_coupling_v3_red_visuals/render.py --preview
# 预览审阅记录保存为 verification/clot_fragment_surface_coupling_v3_red_visuals/20261019/PREVIEW_REVIEW.json 后：
python analysis/clot_fragment_surface_coupling_v3_red_visuals/render.py
```

`VALIDATION.json` 记录真实帧来源及 SHA256、Z 轴屏幕投影误差、视频完整解码与尺寸帧率校验。附加文件有独立清单，不重写上一轮归档。
