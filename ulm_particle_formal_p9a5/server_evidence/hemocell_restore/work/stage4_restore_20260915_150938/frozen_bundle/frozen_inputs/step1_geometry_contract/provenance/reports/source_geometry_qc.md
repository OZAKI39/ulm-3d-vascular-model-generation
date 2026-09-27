# 已有 QC 与本轮独立读取检查

EXISTING_PROJECT_QC：下表来自已有源记录及明确标注的阈值比对；没有重跑 repair/remesh/self-intersection。

| item | status | finding | source |
|---|---|---|---|
| watertight | PASS | s4 topology true；axis rigid QC same；本轮m STL exact-edge incidence亦闭合 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json |
| single component | PASS | 既有component_count=1；独立exact-coordinate adjacency=1 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json |
| boundary edges | PASS | 既有0；本轮0 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json |
| nonmanifold edges | PASS | 既有0；本轮0 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json |
| self intersection | PASS | 既有0/5402候选；本轮未复算；刚体变换保持几何交关系 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json |
| degenerate triangles | PASS | 既有0；本轮精确零面积0，未重跑历史近零面积阈值 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json |
| normal consistency | PASS | 既有winding_consistent；本轮共享边方向一致；入口m STL外法向≈+Z | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json |
| triangle angle ≥20° for every active face | FAIL | 活动区域 minima=[15.99819806657744, 17.527033106729867, 18.963475614794085, 21.574683072286202]；历史 quality hard_gate=false，只标 finite_metrics，并非逐面20°通过 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/crossseam_mesh_quality_qc.json |
| aspect ratio ≤5 for every active face | FAIL | outlet_01 max=5.196656573395647；历史bad fraction=0.00022016732716864817；并不等同全表面硬阈值失败 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/crossseam_mesh_quality_qc.json |
| cap planarity at nominal 1e-5 um | UNVERIFIED | 既有distal开口误差=[6.885416564239577e-06, 1.2079953231120683e-05, 5.4833143296253994e-06, 1.3014544371850109e-05] um；其中两项超过配置1e-5，但open_profile_qc未把planarity加入checks；不能继承整体PASS作为cap硬阈值验收 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/extension_geometry_qc.json |
| surface distance: accepted original-side collar | PASS | P95=0.007965553277467063 um ≤0.05；FAR_CORE max motion=0；只对声明区域和方法有效 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/active_collar_original_side_distance_qc.json |
| surface distance: legacy bidirectional collar | FAIL | 原始诊断P95=0.09527295926669675 um、max=0.20809758402773543 um；已被历史流程改为DIAGNOSTIC_ONLY，不能隐藏失败 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/legacy_reclassified_collar_bidirectional_distance.json |
| radius fidelity | PASS | 既有radius_fidelity.status=PASS；采样截面证据，非所有位置最小口径证明 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/radius_fidelity.json |
| manual final capped review | UNVERIFIED | 保存状态MANUAL_REVIEW_REQUIRED；没有在本轮范围内确认最终人工评审完成 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/run_summary.json |


INDEPENDENT_READ_ONLY_CHECK：`independent_geometry_check.json` 单独记录二进制 STL 的原始坐标测量。三角面=67262，按精确坐标等价计数的顶点=33633，连通体=1，边界边=0，非流形边=0，重复面=0，共享边方向不一致=0。该计数与 VTP 的 73416 存储点数量不同，原因是 VTP 存在重复坐标点、STL逐面保存顶点；没有焊接或修改文件。

最终m STL cap平面误差=[1.060932738459606e-11, 1.691960281621276e-11, 8.200240986385859e-12, 1.6642996502559376e-11] m。STL float32量化及不同中心/法向定义会影响末位。只提供测量；未宣称全套HemoCell体素化、边界或几何质量验收通过。
