# 单位与坐标合同

FINAL_STL_UNIT = m

最终文件为 `/home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/geometry/cfd_surface_axis_aligned_inlet_m.stl`。米制结论来自 s4 显式导出、其 meter_scale_qc、下游刚体变换合同和 Base 网格的实际 geometry/hash 引用；没有根据坐标范围猜单位。旋转本身 scale=1，保留长度单位。历史轴对齐导出脚本在当前工作区未找到，精确导出代码位置为 UNVERIFIED；变换矩阵、结果和消费关系有实际文件证明。

| contract | value |
|---|---|
| SOURCE_GEOMETRY_UNIT | Raw fMOST xyz is treated as voxel indices; raw radius physical calibration UNIT_UNVERIFIED. Saved normalized SWC/ROI/lumen uses code-contract um. |
| S3_INPUT_UNIT | um for normalized SWC, ROI and reference lumen; source calibration is a separate unresolved issue |
| S3_OUTPUT_UNIT | geometry metadata um; normals dimensionless; pressure Pa; volumetric flow m3/s; s3 does not create a new STL |
| S4_INPUT_UNIT | um |
| S4_INTERNAL_UNIT | um |
| S4_OUTPUT_UNIT | um STL/VTP plus explicit m STL copy |
| FINAL_STL_UNIT | m |
| FINAL_COORDINATE_SYSTEM | axis_aligned_CFD |
| HEMOCELL_TARGET_UNIT | Lattice units inside voxelizer; configured physical dx is m/LU. refDirN/reference-axis extent controls mesh scaling; no run-specific mapping selected. |


| operation | object | before | after | factor | evidence |
|---|---|---|---|---|---|
| raw xyz multiplied by configured spacing | SWCData.points_um | raw voxel index | program-contract um | [1, 1, 2] | [swc_io.py:185](</home/lzy/projects/ulm_3D_vascular/utils/rodent_vasculature/swc_io.py:185>) |
| radius compensation and H5 diameter adapter | feed_radii and points fourth column | source_radius_um | diameter_um | 2 * 0.91; positions not rescaled here | [ultraliser_backend.py:250](</home/lzy/projects/ulm_3D_vascular/utils/cfd_lumen/ultraliser_backend.py:250>) |
| s2 creates meter mirror | mesh_m.vertices | um | m | 1e-06 | [ultraliser_qc.py:51](</home/lzy/projects/ulm_3D_vascular/utils/cfd_lumen/ultraliser_qc.py:51>) |
| s3 one-dimensional resistance calculation only | root_radius_m / edge length and radii | um | m | 1e-06 | [one_d_flow.py:114](</home/lzy/projects/ulm_3D_vascular/utils/cfd_preprocess/one_d_flow.py:114>) |
| s4 meter export after cap | meter_mesh.vertices, copied from final_mesh | um | m | 1e-06 | [vmtk_qc.py:748](</home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/vmtk_qc.py:748>) |
| historical rotation about inlet area-centroid pivot | whole final surface and all port patches | s4 frame (m or um form) | axis_aligned_CFD (same length unit) | 1.0 | [anatomical_to_cfd_transform.json:109](</home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/transform/anatomical_to_cfd_transform.json:109>) |
| future HemoCell loader scaling; read only, not applied | TriangularSurfaceMesh in memory | STL coordinates | LU | refDirN / deltaX_reference | [triangularSurfaceMesh.hh:1568](</home/lzy/projects/hemocell_starter/palabos/src/offLattice/triangularSurfaceMesh.hh:1568>) |


`s3` 对长度/半径换算只用于阻力与流量计算，不缩放或生成 STL。s4 输出 um、m 两份，最终包已是 m，未来不得重复应用 1e-6。H5 的 0.91 是半径重建补偿，不是全局单位转换，也不修改原始 SWC 半径。

保存变换角度为 42.264604756361436°，以 inlet 面积重心为 pivot，入口外法向映射到 +Z，内法向为 −Z。绕非零 pivot 的 4×4 矩阵有非零平移项；这不是额外局部变形。m 与 um 两套齐次矩阵不可混用。`metadata/anatomical_to_cfd_transform.json` 包含原始数值及逆变换。

原始 fMOST spacing=[1,1,2] um 来自配置/保存 manifest；独立采集标定和原始半径校准仍为 UNIT_UNVERIFIED。`anatomical_fMOST` 是历史文件的坐标系名称，不构成解剖轴、手性、全脑配准或世界原点的实测证明。
