# Candidate 0：黑色背景、端口标签、坐标网格及显示筛选

仅更新 candidate 0。保持蓝色微泡、透明血管、原视线方向、英文单列顶部信息、24 fps 和 120 秒慢放时长。未运行 CFD 或重新计算轨迹。

## 画面与指定的右图风格

- 背景为纯黑色；微泡保持统一蓝色 `#3399ff`，血管的材质、透明度、原视线方向和顶部英文信息继续沿用。
- Inlet、Outlet 1、Outlet 2、Outlet 3 在整段动画中显示。箭头由官方 `mesh/SV_MESH/mesh-surfaces/INLET.vtp` 和 `OUTLET_0*.vtp` 截面中心经过 candidate 0 刚体变换后投影得到，不按视觉位置猜测编号。
- 坐标风格参考用户指定命令 `cd /home/lzy/projects/vascular_printing && ./run_python.sh s1-2_swc_roi_generate_human.py configs/swc_roi_generate_human_spatial.yaml` 的**右侧子图**。仅阅读该流程的显示代码与已有截图，没有执行它的几何、ROI 或采样流程。
- 实际代码链：`spatial_display_adapter.py` → `compact_annotations.py` → `utils/rodent_vasculature/interactive.py`。抽取其中坐标显示函数为 `scripts/spatial_axes_reference.py`，不修改参考项目。
- 使用同样的红色三维包围框 `#FF5C5C`、4 px 框线、填充 opacity = 0.025；六个面的灰色网格为 `#BFC7D5`、1 px、opacity = 0.28。
- X/Y/Z 每轴分别 4/4/5 个刻度，每个轴仅显示一组标签。刻度贴在实际三维框边并朝外，毫米数值保留两位小数，Arial 14 px 刻度和 16 px 标题。没有保留中间预览中的二维 u/v 网格。
- 此处数据已经是 candidate 0 的 **mm 坐标**，没有再套用参考项目的 μm→mm 除以 1000。数值是姿态变换后的真实 X/Y/Z 坐标。
- 此红框表示当前 CFD 血管的显示包围范围，不是参考项目所选择的 ROI，也不是新生成的物理边界。
- 保留当前相机的 Z-up 视向；参考子图是 Y-up，复制的是网格、字体及贴边刻度样式，不重排血管姿态。Z 轴标题额外外移 18 px，避免与中部刻度数值相碰。
- 左下角采用参考图的彩色 X/Y/Z 方向轴及白色单位标签。为避免当前视角下 X、Y 标签相碰，方向箭头显示长度不同；方向小图不表示距离。
- 主图按**完整三维包围框、刻度和标题共同占据的范围**最大适配到顶部文字以下，所有元素不裁切。真实 XYZ 包围框比血管投影更宽，因此其适配缩放与先前没有坐标框的版本不同；不能以裁掉框边或刻度换取血管放大。
- 相机、安全显示范围、血管投影范围、所有实际坐标标签位置、官方端口位置与源代码哈希保存在 `data/render_manifest.json` 和 `data/axes_style_reference.json`。

## 用户选择的显示筛选

按本轮明确选择，整个隐藏以下轨迹（动画及静态轨迹图使用同一个名单）：

1. `trajectory_catalog.csv` 中 `SOLVER_FAILURE`：2 条，ID 440、1505。
2. 已有 `cfd_wall_audit/acceleration_overlap.csv` 中 `window_any_eta_gt_0p5=True`：58 条。这表示此前筛出的显著加速事件前后各约 0.06908 s 的窗口中访问过 η>0.5 单元，其中 η=|div(u)|/||grad(u)||F；使用已有数据，不重新设定阈值。

两类不重叠，共隐藏 60 条，显示 1505/1565 条。没有因为正常加减速、低速、近壁运动或到达观察终点而额外删除轨迹。剩余 143 条显著加速候选也未仅凭加速现象被隐藏。

**58 条属于数值风险筛查，不是逐条证实的物理错误；保留显示的轨迹也不能因此被称为已独立证明正常。**这是展示筛选，不是数据修复或新的流场精度结论。原始 1565 条轨迹、终态分类及统计完整保留。

显示轨迹的原始终态为：Outlet 1 = 452，Outlet 2 = 483，Outlet 3 = 429，观察终点截尾 = 141。此为筛选后的显示计数，不能替代原始模拟流出统计。

## 文件与复现

- `animations/microbubble_physical_time.mp4`、`OPEN_RESULTS.html`：当前动画及预览。
- `data/display_selection.csv`：全部 1565 个 ID 是否显示及隐藏原因。
- `data/render_manifest.json`：端口、相机、坐标、筛选来源哈希和渲染信息。
- `data/animation_samples.npz`：完整未删减的慢放采样数组；渲染时才应用显示名单。
- `data/render_qc.json`：本次完整视频、筛选名单、坐标及数据一致性检查；检查代码为 `scripts/check_render_candidate0.py`。
- `revisions/before_black_labels_maxfit_20260929/`：本轮调整前的完整可视化备份。

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
/workspace/brava_flow_roi_18mlmin_20260928/population_v1/env/bin/python -u -B \
/workspace/brava_flow_roi_18mlmin_20260928/population_v1/scripts/render_population.py \
--run-name realization_0000_contactfix --pose 0 --playback-speed 0.2 --hide-audited-anomalies
```

只调整显示，积分 dt 仍为 0.5 ms。视频时间仍映射到原始物理时间 0–13.240815660385584 s。入口平均速度仍按入口网格面积 Q/A 得到 0.2695615893 m/s，画面显示 0.270 m/s（18 mL/min）。服务器 RTX 4090 执行 VTK EGL 渲染，CPU 执行 libx264 编码。
