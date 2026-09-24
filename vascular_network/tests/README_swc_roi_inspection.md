# 原始 SWC 的三个 ROI 交互检查

此界面用于回答“局部半径变化是否已经存在于 SWC 文件中”。默认读取最近一次成功的人脑采样结果，按 `selection_rank` 选取前三个代表 ROI；本次对应 `BG001.CNG.swc` 的锚点 **8、155、647**。

在项目目录运行：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tests/inspect_swc_rois.py
```

三个视窗可独立旋转、缩放。把鼠标放在血管所在视窗中操作：

| 操作 | 功能 |
| --- | --- |
| 左键拖动、滚轮 | 旋转、缩放 |
| 右键单击血管，或将鼠标放在血管上按 `P` | 选取最近的原始节点，显示编号、SWC 行号、半径、直径及父子节点半径；橙色点表示被选节点中心 |
| `1` | 管面：仅在原始节点设置原半径截面，节点之间直线连接 |
| `2` | 原始中心线和节点；节点标记大小固定，仅用于定位，不代表半径 |
| `3` | 每个原始节点放置一个球，球半径严格取自 SWC 第六列；球的重叠不是额外血管 |
| `C` | 在统一颜色与半径着色之间切换；三个 ROI 使用同一色域 |
| `R` | 恢复三个初始视角 |
| `S` | 保存当前截图 |

界面采用英文标注，避免本机缺少中文字体导致文字无法显示。底部 `r` 是半径，`d` 是直径，`ratio` 是所选节点半径与父节点半径之比。父子节点信息来自完整源文件，可能包含当前 ROI 外的节点。

也可以在启动时更换三个 ROI，例如：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tests/inspect_swc_rois.py --anchors 8 155 647
```

用 `--run-dir outputs/human_brava/sampling/20260923_185451_radius_plus_structure_k5` 固定本次采样批次。锚点编号应来自该批次 `manifests/candidate_rois.csv`，且三个 ROI 必须属于同一份源 SWC。

## 数据含义与观察边界

采样结果只提供 ROI 的原始节点编号集合。坐标、半径、父子连接均重新读取原始 SWC，单位依据该批次保存的源配置；本次为 **mm**。界面不读取处理后的 ROI 几何，不做重采样、中心线平滑、半径平滑、半径修补或隐式表面重建，也不绘制方向箭头、ROI 方框、坐标网格和其他背景血管。

为保持“只显示原始 SWC”，不显示裁剪时新增的边界节点，也不显示只有一个端点在该 ROI 内的边。因此与主程序相比，末端可能略短，节点数分别是 **139、42、54**，少于包含人工裁剪节点的 ROI。保留采样时选定的连通成分，不额外加入同一包围盒内其他未选中的血管。

SWC 给出离散位置、半径和连接，并没有提供真实血管壁表面。管面视图只是这些样本的分段连接，弯折或分叉处可能出现管面相交；不能把这种相交理解为真实管腔。判断节点半径时，应结合节点读数、球视图及导出的原始数值。球面和圆截面通过有限多边形显示，但半径参数不经过缩放或修正。

原始记录中确有明显变化。例如节点 **2768 → 2769 → 2770** 的半径为 **0.31 → 0.62 → 0.31 mm**，节点 **16 → 17** 为 **0.9796 → 0.31 mm**。严格经过这些样本的平滑显示仍会产生粗细变化；单靠渲染不能消除这些变化，同时又保持全部原值。这些记录说明源数据存在变化，并不能证明真实血管具有相同程度的局部狭窄，更不能据此判定生理正常或异常。

## 可复核输出

每次运行在 `tests/swc_roi_inspection_output/时间戳/` 保存源文件路径与 SHA-256、三份原始节点表以及三份按变化倍数排序的相邻半径表。节点表记录源文件行号，便于逐行核对。截图也保存在同一目录。`source_audit.json` 说明单位和节点选择规则。

无窗口生成三个视图的检查图：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tests/inspect_swc_rois.py --off-screen
```

本次预览和数值记录位于 [swc_roi_inspection_output/verified_preview](swc_roi_inspection_output/verified_preview)。自动验证覆盖原始 ID 和坐标保留、真实节点半径、分叉管面截面半径、球半径、节点拾取对应关系、三种视图切换及源文件不变性：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python -m pytest -q tests/test_swc_roi_inspector.py
```
