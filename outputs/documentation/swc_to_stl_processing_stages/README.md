# 4.4.2 整体处理流程：三维可视化

本目录继续使用 anchor 3274 的真实保存 ROI 和当前 `swc_stl_model_generate.py` 配置，解释 `polylines-with-spheres → solid fill → DMC → mesh processing`。

## 六个面板分别表示什么

1. **A：带球端折线。** 折线负责连接实际中心线采样点；每个采样点的球由其 feed radius 决定大小。球和相邻线段共同标记血管应占据的空间，而不等于最终 STL 三角面。
2. **B：壳与实体。** 表面体素化先得到封闭壳；`--solid --voxelization-axis xyz` 再把壳内部标为实体。图中橙色是一层教学性壳单元，蓝绿色是剖切后的内部单元。该体素展开使用真实输入和 `6.0` voxels/µm，但不是 Ultraliser 导出的内部 volume。
3. **C：DMC。** 官方 `--isosurface-technique dmc` 从实体域的内外分界产生初始网格，共 `145172` 个三角形。
4. **D：Laplacian。** 当前程序先执行 `10` 次 Laplacian 处理。它主要根据邻域移动顶点，缓和局部高频起伏；此处三角形数仍为 `145172`。
5. **E：自适应优化与平滑。** `--adaptive-optimization --optimization-iterations 5 --smooth-iterations 5` 调用的优化过程包含法向处理、表面平滑、细化以及稠密/平坦区域简化。网格减少到 `47230` 个三角形，而不是简单地对所有三角形做同一种移动。
6. **F：watertight。** 最终阶段检查并修复为统一水密表面。捕获的 watertight STL 与当前验收 STL 的 SHA-256 完全一致，因此 C–F 确实来自生成当前验收结果的同一参数组合。

## 文件真实性边界

- `official_stages/*.vtp` 和 `official_local_views/*.vtp`：由官方 `ultraVessMorpho2Mesh` 的 DMC、Laplacian、optimized 和 watertight STL 转换得到，属于真实阶段几何。
- `explanatory_surface_shell.vtu`、`explanatory_solid_interior.vtu`：根据同一实际中心线、feed radius 与体素间距构造的教学性体素拆分；用于说明“壳内填充”，不冒充 Ultraliser 内部 volume。
- `official_invocation.txt`：本次阶段捕获使用的完整命令。
- `processing_stage_metadata.json`：阶段计数、参数、SHA-256 核对和文件路径。
