# 自动一致性检查

| id | status | description | evidence |
|---|---|---|---|
| A | PASS | 冻结 s3 input_manifest/geometry_reference 与报告的 August lumen 来源、SHA一致；不是断言当前null选择器重跑同一批次。 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/input/input_manifest.json; /home/lzy/projects/ulm_3D_vascular/outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/input/geometry_reference.json |
| B | PASS | 当前 s4 配置、历史open run引用同一s3包和原始lumen。 | /home/lzy/projects/ulm_3D_vascular/configs/cfd_surface_prepare.yaml; /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_remesh_anchor003274_20260826_201129/input/original_surface_reference.json |
| C | PASS | s4→刚体变换→实际冻结网格链一致；直接s4 STL与当前final STL字节不同，差异为已记录旋转。 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/input/source_provenance.json; /home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/qc/full_fluid_center_containment.json; /home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/qc/input_manifest.json |
| D | PASS | 单位由s4缩放代码、meter QC和m/um齐次变换共同证明；历史旋转导出器具体源码未找到。 | /home/lzy/projects/ulm_3D_vascular/utils/cfd_surface_prepare/vmtk_qc.py; /home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/transform/anatomical_to_cfd_transform.json |
| E | PASS | s3/s4/最终VTP、cap STL和合同均为1入口3出口；port_id、global edge和terminal身份一致。 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/roi/port_classification.csv; /home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/boundaries/boundary_manifest.csv; boundary_contract.json |
| F | PASS | 源STL与复制STL SHA256完全一致；未重写或变换。 | /home/lzy/projects/ulm_3D_vascular/outputs/cfd_flow/axis_aligned_inlet_geometry_anchor003274_20260829_111451/geometry/cfd_surface_axis_aligned_inlet_m.stl; /home/lzy/projects/compre_output/step1/20260912_215759/geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl |


附加明确差异：s4未旋转STL与当前最终STL直接字节相等 = FAIL（预期不同）；历史Windows路径直接用于当前Linux = 不可用；当前FRESH_STEADY入口未实现自动新稳态执行。这些差异已解释且未修正原工程。

只读前后完整目录元数据与证据文件SHA核对见 read_only_verification.json；没有执行源工程测试或solver。
