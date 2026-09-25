# 24项完成标准逐项回答

| question number | answer |
|---|---|
| 1 | 后续候选为实际冻结CFD使用的axis-aligned完整m STL；当前HemoCell helper适配尚未就绪。 |
| 2 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/geometry/cfd_surface_axis_aligned_inlet_m.stl |
| 3 | 历史axis_aligned_inlet geometry run的刚体变换阶段；精确导出器当前源码位置UNVERIFIED。 |
| 4 | fMOST原始标注SWC→单连通analysis→ROI003274→Ultraliser lumen。 |
| 5 | 坐标按spacing归一；ROI裁切；radius feed0.91/H5直径；Ultraliser重建；局部切口/TPS延长/active remesh；cap；整体刚体旋转。 |
| 6 | m，来自代码/保存单位合同及下游消费证据。 |
| 7 | 存在um→m=1e-6；旋转scale=1；本轮未变换。 |
| 8 | s4 utils/cfd_surface_prepare/vmtk_qc.py:meter_mesh.vertices；另有s2 meter mirror和s3阻力内部SI换算。历史轴对齐exporter位置UNVERIFIED。 |
| 9 | 使用s4含延长/cap表面的刚体旋转派生版；非原始lumen，也非未旋转s4表面直接作为当前坐标。 |
| 10 | 1 |
| 11 | 3 |
| 12 | cut_000→inlet；cut_001→outlet_01；cut_002→outlet_02；terminal_000→outlet_03，全部ASSUMED生理角色。 |
| 13 | 四个人工延长段的远端cap，完整最终坐标见boundary_contract.json。 |
| 14 | 见boundary_contract.json center：最终cap面积加权重心；保留s3/s4其他中心定义。 |
| 15 | 见boundary_contract.json normal：从最终cap三角形计算的单位外法向。 |
| 16 | 见boundary_contract.json area(m²)，不是s3 πr²名义面积，也不是canoND矩形面积。 |
| 17 | s3源图方向和cut/terminal映射；s4 Hungarian匹配cap到计划端点，写VTK标签；旋转保留ID。 |
| 18 | s3 summary/geometry_reference + port_classification.csv + extension_plan.csv + BC JSON，s4按port_id交叉核对。 |
| 19 | YES，4个人工延长段。 |
| 20 | YES，4个cap共191面属于最终STL。 |
| 21 | 来源SWC重建的CORE45722面；其中FAR_CORE45512面保留、core-side collar210面已重网格；非原始实验表面逐点真值。 |
| 22 | EXTENSION21349面，入口7414、出口4502/4371/5062面；归属详情见region_contract.json。 |
| 23 | 已有拓扑/刚体/半径QC，本轮另有轻量只读检查；逐指标限制与未完成的人工评审保留。 |
| 24 | 最终SHA256=840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb；provenance.json/source_paths.json/SHA256SUMS保存关键来源、代码、配置、metadata和副本校验。 |
