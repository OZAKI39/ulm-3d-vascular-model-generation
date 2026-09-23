# 人工平面端口契约 v2 的修订依据

Stage 1.5 的正式结论永久保留为 **FAIL**。其三个候选的旧 3D scalar-area 相对误差约 6.3e-10～3.0e-9，超过冻结的 1e-12；longdouble 独立检查确认这不是求和精度问题。依据：[Stage 1.5 area analysis](../stage01_5/area_invariant_analysis.json)。

WALL 是真实血管表面，坐标、连接和标签完全冻结。CAP 是人工 CFD 边界，不是 anatomical surface。原始 float32 VTP 留下约 6e-12～1.5e-11 m 的 rim 非共面性。同一条未移动的 3D rim，在不同内部三角剖分下可以产生不同的 3D 标量面积；它不适合作为人工平面端口的身份定义。

v2 端口定义为固定 3D rim 顶点及边、原始 x0、原始 outward normal 和固定正交基上的投影多边形。正式面积 A_port 是投影多边形面积，正式中心是该平面多边形面积中心。由有序 3D rim 计算的向量面积是独立方向及面积不变量。Stage 0 source_contract.json 保持原样，v2 引用其 SHA256。

原 rim 坐标不投影替换，新增内部点才放在原平面。原 3D scalar triangle area 与 scalar-area centroid 保留为 legacy provenance；每个候选仍公开其标量面积变化。修订的是人工端口面积定义，不是放宽旧误差门槛，更不是改变 vascular anatomy。投影面积和向量面积相对误差门槛仍为 1e-12。

全部质量、密度、四面体数和 P2 拓扑成本预算在候选运行前冻结。只检查网格；不创建 FEM 空间，不求解流动。Stage 1.5 FAIL、Stage 2 solver core 和所有历史证据不变。
