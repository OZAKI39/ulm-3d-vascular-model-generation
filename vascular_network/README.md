# vascular A → ROI → Network H0

本目录保留血管预处理、Network A/H0 的源码、配置、测试和必要源数据。

| 阶段 | 入口 |
| --- | --- |
| 小鼠 ROI | `s1-1_swc_roi_generate_mouse.py`、`utils/schmid_pkl/`、`utils/sampling/` |
| Ultraliser 表面 | `s2_swc_stl_model_generate.py`、`utils/cfd_lumen/` |
| 端口与 CFD 表面 | `s3_cfd_1D_data_preprocess.py`、`s4_cfd_surface_prepare.py` |
| Network A 基础图与精确切口 | `scripts/run_a_network_1d0d_boundary_v1.py`、`network_1d0d/` |
| H0 假设与求解 | `scripts/run_a_network_h0_v2.py`、`network_1d0d/idealized_h0.py` |
| 3D H0 算例与验收 | `network_1d0d/fem_h0_case.py`、`scripts/solve_a_h0_fem_remote.py`、`scripts/validate_a_h0_fem.py` |

H0 仍依赖 v1 图/端口，因此保留 v1 的完整必要输入哈希闭包。只选择本案例的源 SWC、ROI、几何和配置，没有复制完整原始脑数据集。
当前报告见 [A_NETWORK_H0_REVIEW_ZH.md](reports/a_network_1d0d_boundary_v2_idealized/A_NETWORK_H0_REVIEW_ZH.md) 和 [ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md](reports/a_network_1d0d_boundary_v2_idealized/ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md)。

原 JSON 中的绝对路径按原样保留以保存来源；[统一回归工具](../sync_metadata/network_h0_particle_20260925/run_checks.py)会在临时目录迁移这些路径并运行全部 49 项 Network 回归。

表面重建仍使用 Ultraliser，`feed_radius_um = source_radius_um * 0.91`，H5 第四列为直径 `2 * feed_radius_um`，CFD 坐标从微米按 `1e-6` 缩放为米。配置与说明在 `configs/`、`docs/`；第三方构建环境按记录重建。
其余 LBM、human/MeVO/TopBrain 的已有源码保留各自范围；并行进行的 TopBrain/BRAVA 大型输出不属于本次血流/Particle 数据快照。
