# 当前方法和数据入口

下面路径均相对于本工程。各算例共用 `mouse_vessel_flow/`，没有复制完整求解器。ROI-only 仍为默认；其他已接受流场有各自的身份、压力条件、队列和时间步。

| 功能或算例 | 当前源码与函数 | 独立配置/数据 |
|---|---|---|
| 几何 | `geometry.py:retained_roi`；`geometry_tools/rodent_vasculature/` 与 `geometry_tools/cfd_lumen/` | `configs/geometry/retained_roi.json`，`inputs/geometry/` |
| ROI-only 压力设计 | `network.py:roi_cache, validate_design`；`network_models/roi_boundary_design.py:regress_h0, design`；`boundary.py:roi_only_design` | `cases/roi_only/case.json`、`results/roi_only/network/` |
| Best-feasible 压力方案 | `network_models/balance_design.py`、`parameter_fit.py`、`parameterized_hydraulics.py` 与 `parameterized_fem_handoff.py` | `cases/best_feasible/case.json`、`results/best_feasible/network/`；本轮核对保存结果，没有重新拟合 |
| Network-H0 | `network_models/idealized_h0.py`、`network_solver.py`、`fem_h0_case.py:serialize_fixed_pressure_xml, create_case` | `cases/network_h0/case.json`、`results/network_h0/network/`；保留 500 条正式记录 |
| 原始 2 mm/s、Stage-Q | 与其他算例共用 `flow_field.py:load_case`、`fem.py:inspect_case` | `cases/mean_2mmps/case.json`、`cases/stage_q/case.json`；只作为历史参考，不冒称新的当前验收 |
| 网格与 CFD | `meshing.py:launch_embedded`、`mesh_tools/tetgen_native.py:generate`；`fem.py:prepare_case, run_case`；`solver_support/native_run.py:run_native` | 各算例原 `solver.xml`；原生安装根目录 `new_vessel_flow_solvers`，软件身份在 `configs/environment.json` |
| 微泡 | `particles/injection.py:saved_events, source_for_case`；`particles/dynamics.py:bind_job`；`particles/trajectories.py:run_saved_cohort` | ROI-only/Best-feasible 为 0.5 ms；H0 为 1 ms。各自保留原流场哈希、粒径、入口队列、契约和完整轨迹 |
| WSS | `wss.py` 导出 `solver_support/wss.py` 的 `p1_gradients, tangential_traction` 等；`validation.py:wss_check` | 原双精度网格、ROI-only 保存壁面量、`validation/wss_audit/` 和 `wss_validation_v2/` |
| 只读核验 | `validation.py:validate`、`particles/analysis.py:verify_tracks`、`scripts/check_independence.py` | 数据身份与原始文件哈希；不进行新正式积分 |

完整函数索引见 `CURRENT_SOURCE_FUNCTIONS.csv`。受原正式结果保护的旧阶段名称（P6/P9 等）仍是共同算法依赖，不能按名称“旧”就删除。`inputs/reference/source_identity/` 中原始源码 ZIP 是字节身份证据，不被导入或执行。原生第三方源码和 WSS 已执行验证的源码证据不冒充第二套活跃 Python 工程。

已移除的新工程旧入口：`scripts/migrate_local.py`、`migrate_sources.py`、`finalize_audit.py`、`inventory_data_review.py`。以前复制源码的同步工作区已注销；正式结果目录没有恢复 `.py/.sh/.ps1` 旧执行脚本。原始历史工程仍在待确认清单中，未宣称它们已全部删除。

本轮未重新执行最佳压力拟合、完整几何重建、网格生成、CFD 或正式微泡批次。原生安装核验读取先前实际小规模试算的身份凭据。旧 Ultraliser 极小模型问题继续按用户要求暂缓处理；细档血管与 O2 敏感性仍为 `SKIPPED_BY_USER`。
