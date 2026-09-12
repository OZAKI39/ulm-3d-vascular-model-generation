# s4 调用链与实际产物链

| function/stage | input | output | unit | boundary role | evidence |
|---|---|---|---|---|---|
| s4.main → load_surface_prepare_config | 当前 configs/cfd_surface_prepare.yaml | 固定 s3_run；严格 VMTK schema | um | 4个边界 | [s4_cfd_surface_prepare.py:42](</home/lzy/projects/ulm_3D_vascular/s4_cfd_surface_prepare.py:42>) |
| load_surface_inputs | s3 summary/geometry_reference/BC/classification/extensions | SurfaceInputs；跨文件核对 port_id | um + SI | 只接受1入3出和1 TRUE_TERMINAL | [io.py:79](</home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/io.py:79>) |
| load_original_surface | geometry_reference 的 lumen_surface_um.vtp | 只读几何的内存工作副本 | um | 原始封闭 lumen | [vmtk_pipeline.py:369](</home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/vmtk_pipeline.py:369>) |
| local_plane_cut | 每个边界的局部邻域与切面 | 四开口 surface | um | 真实重建管腔与人工延长的近端界面 | [vmtk_pipeline.py:372](</home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/vmtk_pipeline.py:372>) |
| run_official_vmtk | open surface + boundarynormal/TPS config | 人工延长段；官方 VMTK 1.5.0历史记录 | um | 5D语义；ExtensionRatio10×meanRadius | [vmtk_pipeline.py:409](</home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/vmtk_pipeline.py:409>) |
| assign_cross_seam_active_entities / entity_remesh_official_vmtk | CORE/EXTENSION、核心侧2层 collar | FAR_CORE保留；ACTIVE跨界重网格 | um | interface会改变，不能说全部core未动 | [guarded_remesh.py:123](</home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/guarded_remesh.py:123>) |
| 实际历史 recovery cap_only | 固定201129 open VTP/STL SHA | 221611封帽结果；没有再次extension/remesh | um | 4个CAP | [cap_request.json:2](</home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/vmtk/cap_request.json:2>) |
| tag_and_export_final_surface | capped VTP + s3 boundary predicted endpoints | CellEntityIds/port_id/region/manifest、STL | um与m | 匈牙利一对一匹配；保留ASSUMED角色 | [vmtk_qc.py:218](</home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/vmtk_qc.py:218>) |
| 下游历史刚体变换 | s4冻结带标签表面 | axis_aligned_inlet_m.stl + tagged VTP + transform | m/um分开保存 | 标签不变；法向/中心随整体旋转 | [anatomical_to_cfd_transform.json:3](</home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/transform/anatomical_to_cfd_transform.json:3>) |
| 当前 Base mesh消费 | 轴对齐wall/outlet、理想数值入口平面 | 冻结mesh；production replay继续用它 | m | 与s4原坐标不同 | [seeder.lua:9](</home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/seeder.lua:9>) |


当前 s4 源码定义完整 local cut → VMTK TPS extension → active collar remesh → cap 路径。当前实际被下游溯源引用的 `vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611` 是恢复封帽结果，其 input/frozen_open_geometry_reference.json 指向 `/home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_remesh_anchor003274_20260826_201129`，记录不重建、不重新延长、不重新重网格。这两条事实需要同时保留，不能伪造“当前 s4 在本 WSL 从零生成该历史结果”。

当前 `load_surface_inputs` 直接 Path.resolve 保存的 Windows geometry_reference 路径，没有本地前缀映射；当前 VMTK 配置也仍指向 Windows 环境。历史快照 YAML 含已删除的 legacy 字段，当前严格 loader 不一定接受；均未修正。

| region | actual evidence | meaning |
|---|---|---|
| ORIGINAL_VASCULAR_CORE | SurfaceRegionId=0: 45722 faces | 源SWC重建管腔；含210个已重网格的核心侧collar面 |
| FAR_CORE | RemeshEntityId=1: 45512 faces | s4 remesh/cap前后局部形状与连接保留，随后所有几何整体旋转 |
| ARTIFICIAL_INLET_EXTENSION | port0: 7414 faces | 人工入口延长壁；无实测血管身份 |
| ARTIFICIAL_OUTLET_EXTENSION | port1/2/3: 4502/4371/5062 faces | 三个人工出口延长壁 |
| CAP | SurfaceRegionId=2: 191 faces | 4个远端封帽，属于最终完整STL |
| INTERFACE_REGION | CORE∩CROSS_SEAM_ACTIVE: 210 faces | 近端核心侧collar；语义由距原始区域最近关系恢复，不是精确实验分界 |


区域计数从实际 VTP arrays 读取；每个 extension 的 port 归属由 open VTP CrossSeamPortIndex 与未旋转三角形精确匹配，再核对最终旋转 VTP 的面连接/标签相同得到。没有写回或变换几何。详情见 metadata/region_contract.json。
