# 小鼠脑微血管：vascular A → ROI → Network H0

当前跨仓库流程与服务器入口见[当前工作流索引](/home/lzy/projects/ACTIVE_VASCULAR_WORKFLOW.md)。

| 阶段 | 当前代码 |
| --- | --- |
| 小鼠源数据与 ROI | `s1-1_swc_roi_generate_mouse.py`、`utils/schmid_pkl/`、`utils/sampling/` |
| ROI 表面重建 | `s2_swc_stl_model_generate.py`、`utils/cfd_lumen/` |
| 1D 数据/端口预处理 | `s3_cfd_1D_data_preprocess.py`、`utils/cfd_preprocess/` |
| CFD 表面与延伸段 | `s4_cfd_surface_prepare.py`、`utils/cfd_surface_prepare/` |
| Network A 基础图与精确切口 | `scripts/run_a_network_1d0d_boundary_v1.py`、`network_1d0d/` |
| Network A H0 | `scripts/run_a_network_h0_v2.py`、`network_1d0d/idealized_h0.py` |
| H0 FEM 算例与物理验收 | `network_1d0d/fem_h0_case.py`、`scripts/solve_a_h0_fem_remote.py`、`scripts/validate_a_h0_fem.py` |

H0 仍读取 v1 的图与端口数据，v1 的生成入口和共享模块保留。
H0 结论见 [A_NETWORK_H0_REVIEW_ZH.md](reports/a_network_1d0d_boundary_v2_idealized/A_NETWORK_H0_REVIEW_ZH.md) 与 [ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md](reports/a_network_1d0d_boundary_v2_idealized/ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md)。

## 表面重建约定

Ultraliser `ultraVessMorpho2Mesh` 是当前表面重建后端。输入为保存的 ROIRecord，源坐标、源半径和 CUT_PORT 元数据保持原始定义。
`feed_radius_um = source_radius_um * 0.91`；H5 `points` 第四列是 `2 * feed_radius_um` 的直径。
`lumen_surface_um.stl` 用于微米坐标检查，`lumen_surface_m.stl` 坐标按 `1e-6` 缩放用于 CFD。
从本目录运行 `python s2_swc_stl_model_generate.py configs/swc_stl_model_generate.yaml`，所用 Python 环境需具备该阶段依赖。
细节见 [ULTRALISER_PIPELINE.md](docs/ULTRALISER_PIPELINE.md)、[CFD_PREPROCESS.md](docs/CFD_PREPROCESS.md)、[CFD_SURFACE_PREPARE.md](docs/CFD_SURFACE_PREPARE.md)。

## 回归与其他分支

```bash
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python -m pytest -q -p no:cacheprovider tests/network_h0 tests/network_1d0d
```

`s5_cfd_flow_solve.py` / `utils/cfd_flow/` 为现存 LBM 分支；其他 human、MeVO、TopBrain 预处理和可视化也保留。它们有独立入口或测试，本次没有足够证据判定为废弃代码。
原 README 和历史一次性脚本的归档入口为 [vascular_workflow_cleanup_20260924T214356Z](/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z)。
