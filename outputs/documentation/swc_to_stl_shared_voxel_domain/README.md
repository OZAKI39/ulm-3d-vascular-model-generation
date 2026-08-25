# SWC 到统一 STL：共享体素域可视化

本目录由 `doc_visualize_v2.py` 基于当前建模代码实际使用的数据生成。

## 如何阅读总览图

1. **A：实际输入。** 保存的 ROI 包含中心线节点、连接边和每个节点的源半径；正式送入 Ultraliser 的半径为源半径乘以 `0.91`。
2. **B：反例式对照。** 三条分支各自生成管面时，虽然都位于正确的三维坐标，它们仍是三个相交或重叠的表面对象；分叉处可能含接缝、内表面或重叠关系。
3. **C：共享占据域。** 每条中心线段及其两端半径决定一组占据单元，所有分支对同一个布尔体素数组执行逻辑并集，因此分叉处只保留“内部/外部”的整体分类。这里的单元间距是 `0.166667 µm`（`6.0` voxels/µm）。
4. **D：统一边界。** 正式 Ultraliser 输出从整体形态提取一个三角形外边界；完整表面有 `47230` 个三角形、`1` 个连通分量，并且 watertight=`True`。

因此，原文所说的“统一内腔边界”不是指把若干分支 STL 文件摆到一起，而是先在同一个空间域中决定整个血管树哪里属于管腔内部，再只提取这一个整体的外边界。

## 文件

- `swc_to_stl_shared_domain.png`：四阶段总览图。
- `panels/*.png`：四个可单独放大的三维面板。
- `source_centerline_and_radius.vtp`：实际 ROI 中心线，包含 `source_radius_um`、`feed_radius_um`、`degree` 和 `local_node_id` 点数据。
- `explanatory_shared_occupancy.vtu`：局部解释性占据单元，可在 ParaView 中旋转和裁切。
- `explanatory_occupancy_boundary.vtp`：上述解释性占据的外层方格边界。
- `accepted_ultraliser_surface_local_view.vtp`：正式表面在图示立方体内的裁切视图；仅用于显示，不用于 QC。
- `visualization_metadata.json`：数据来源、参数、计数和解释边界。

## 重要边界

`explanatory_shared_occupancy.vtu` 是根据实际中心线、实际 feed 半径和正式体素间距编写的教学性胶囊并集。它说明“全部分支进入一个共享占据域”这一数据结构和几何含义，但不是 Ultraliser C++ 程序导出的内部体素文件，不能用它替代正式 Ultraliser 结果或正式 QC。
