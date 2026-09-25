# HemoCell 几何兼容性只读结论

NOT_READY_FOR_HEMOCELL_VOXELIZATION_TEST

官方 examples/pipeflow 的 STL 入口为 getFlagMatrixFromSTL，读取 geometry/refDir/refDirN，交给 TriangleSet(DBL) → DEFscaledMesh → TriangleBoundary3D → VoxelizedDomain3D(inside)。证据：[pipeflow.cpp:55](</home/lzy/projects/hemocell_starter/examples/pipeflow/pipeflow.cpp:55>)、[voxelizeDomain.cpp:93](</home/lzy/projects/hemocell_starter/helper/voxelizeDomain.cpp:93>)。

该 loader 不读取 STL 的物理单位标签。Palabos 取 referenceDirection 的 bbox跨度 deltaX，用 refDirN/deltaX 归一化到 LU，先移到最小坐标原点，再增加 margin+extraLayer；返回原始坐标单位下的 dx=deltaX/refDirN。物理配置 dx 则独立定义 m/LU。因此未来必须使 refDirN×物理dx 与选轴的实际米制长度一致，并保存坐标变换；本轮未选择 refDirN 或 dx，也没有把旧 CFD dx 搬过来。证据：[triangularSurfaceMesh.hh:1568](</home/lzy/projects/hemocell_starter/palabos/src/offLattice/triangularSurfaceMesh.hh:1568>)。

该 helper 调用 inflate，默认量为1e-3 LU；随后无条件在 X 两端复制邻近切片以打开端部。证据：[triangularSurfaceMesh.h:332](</home/lzy/projects/hemocell_starter/palabos/src/offLattice/triangularSurfaceMesh.h:332>)、[voxelizeDomain.cpp:140](</home/lzy/projects/hemocell_starter/helper/voxelizeDomain.cpp:140>)。当前几何入口外法向约+Z、另有三个不同方向出口，不能直接视为这个“两端X管道”接口已经兼容。STL也不携带 VTP 的 port/region labels；本轮没有实现边界映射。

该 closed-domain inside 路径要求有可识别内外的封闭表面；源码注释和流程支持此预期，但未找到这里自动检查watertight后才执行的显式完整验收。源几何已有封闭/单连通证据，仍不能宣称实际voxelization通过。

本轮结论依据：几何候选已固定且轻量拓扑检查通过；现有helper方向假设、实际物理尺度/refDirN契约、人工评审与部分QC阈值差异仍需评审。没有执行HemoCell、改变core/Palabos/examples、生成case或启动Step2。
