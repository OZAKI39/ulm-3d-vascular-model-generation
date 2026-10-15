# clot 科研快照（2026-10-15）

本分支 `clot-clearance/research-sync-20261015` 从 2026-10-11 科研分支继续，补齐之后的原生修复、静态 CFD 接口、能量/软化可行性、事件积分及真实 V3 集成的代码、必要输入、数值证据、执行日志和正式图表。同步过程没有运行新的物理工况。

## 最新结果

变分渐进软化已接入真实 V3 血栓，保留 1927 材料点、1561 固定锚点和 157637 原始键。唯一历史 C1 Gaussian 工况完成 0–2.33 ms，状态 `FULL_PROTOCOL_NUMERICAL_PASS`，科学性质 `EXPLORATORY_UNCALIBRATED`。

- 本体/界面完全失效键分别为 9567/420；180 材料点脱离，占初始材料质量 9.34%，形成 174 个游离分量。
- 最大最终游离分量 3 点，离散体积等效直径约 89.5 µm。没有 ≥200 µm 游离块；1–2 节点碎屑占脱离质量 98.33%，不能称为临床微栓子。
- 10502 接受步，1182 拒绝试步，2307 接受细分步；最大归一化能量余额 0.0814%。主段墙钟 3130.263 s，峰值 RSS 0.996 GiB，单 rank、单线程、无 GPU。
- 最大粒子接触重叠 0.436dx；峰值来自 0.376 ms 两个仍附着材料点的非原始邻接自接触，接触几何精度未认证。本次未触发 Y 壁罚接触。
- 本体宏观强度与实际 Gc 未认证；界面软化尾部是工程假设。旧 G0 理论/离散值仍是旧瞬断模型标称参考。背景 COM 阻力为单向近似，人工 Pa 不是超声功率；没有接入真实微泡 streaming、双向 FSI 或正式 AI 训练。

[当前开发入口](../ACTIVE_DEVELOPMENT_ENTRYPOINTS.md) · [28 项中文报告](../verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/REVIEW_LPS_VARIATIONAL_SOFTENING_V3_INTEGRATION_AND_SINGLE_PILOT_V1.md) · [最终状态](../verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/FINAL_STATUS.json) · [AI 探索性记录](../verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/V3_SOFTENING_PILOT_AI_RECORD.json)

[科研说明（相对链接）](clot_rupture.md) · [WSL 原文快照](clot_rupture.original.md) · [文档来源与 SHA](DOCUMENT_SNAPSHOT.json)

## 六组科研图

| 内容 | 文件 |
|---|---|
| V3 材料、锚点与混合键 | [PNG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig1_V3_material_and_bonds.png) · [SVG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig1_V3_material_and_bonds.svg) |
| 首次软化、失效与脱离 | [PNG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig2_first_damage_events.png) · [SVG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig2_first_damage_events.svg) |
| 真实点云演化 | [PNG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig3_actual_V3_evolution.png) · [SVG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig3_actual_V3_evolution.svg) |
| 新旧模型碎片统计 | [PNG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig4_historical_fragment_comparison.png) · [SVG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig4_historical_fragment_comparison.svg) |
| 受力、损伤与独立能量账本 | [PNG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig5_load_motion_energy.png) · [SVG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig5_load_motion_energy.svg) |
| 计算成本与 AI 使用边界 | [PNG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig6_cost_and_AI_admission.png) · [SVG](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/Fig6_cost_and_AI_admission.svg) |

[HTML 浏览入口](../results/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/OPEN_RESULTS.html)。GitHub 可直接查看 Markdown、PNG、SVG；HTML 请下载后打开，或在 `clot_clearance` 下执行 `python -m http.server 8000 --bind 127.0.0.1`。

## 本次补齐的科研报告

历史失败、预算中止与最终通过的运行均保留各自结论，不把中间失败改写为通过。

| 任务 / Run ID | 报告 |
|---|---|
| `lps_native_state_identity_and_fracture_bugfix_v1` / `20261011T203015Z` | [中文报告](../verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z/REVIEW_LPS_NATIVE_STATE_IDENTITY_AND_FRACTURE_BUGFIX_V1.md) |
| `lps_final_native_bugfix_full_v3_regression_v1` / `20261011T211309Z` | [中文报告](../verification/lps_final_native_bugfix_full_v3_regression_v1/20261011T211309Z/REVIEW_LPS_FINAL_NATIVE_BUGFIX_FULL_V3_REGRESSION_V1.md) |
| `lps_fracture_timelayer_and_timestep_sensitivity_v1` / `20261011T221502Z` | [中文报告](../verification/lps_fracture_timelayer_and_timestep_sensitivity_v1/20261011T221502Z/REVIEW_LPS_FRACTURE_TIMELAYER_AND_TIMESTEP_SENSITIVITY_V1.md) |
| `y_clot_native_static_cfd_load_interface_v1` / `20261012T185858Z` | [中文报告](../verification/y_clot_native_static_cfd_load_interface_v1/20261012T185858Z/REVIEW_Y_CLOT_NATIVE_STATIC_CFD_LOAD_INTERFACE_V1.md) |
| `lps_3d_fracture_surface_energy_identification_v1` / `20261012T215617Z` | [中文报告](../verification/lps_3d_fracture_surface_energy_identification_v1/20261012T215617Z/REVIEW_LPS_3D_FRACTURE_SURFACE_ENERGY_IDENTIFICATION_V1.md) |
| `lps_3d_fracture_surface_energy_identification_v1` / `20261012T215617Z` | [中文报告](../verification/lps_3d_fracture_surface_energy_identification_v1/20261012T215617Z/source_inputs/REVIEW_LPS_WIDE_PURE_SHEAR_PVCE_FEASIBILITY_V1.md) |
| `lps_controlled_mode_i_separation_energy_v1` / `20261012T230100Z` | [中文报告](../verification/lps_controlled_mode_i_separation_energy_v1/20261012T230100Z/REVIEW_LPS_CONTROLLED_MODE_I_SEPARATION_ENERGY_V1.md) |
| `lps_mode_i_low_load_admission_and_solver_recovery_v1` / `20261013T122000Z` | [中文报告](../verification/lps_mode_i_low_load_admission_and_solver_recovery_v1/20261013T122000Z/REVIEW_LPS_MODE_I_LOW_LOAD_ADMISSION_AND_SOLVER_RECOVERY_V1.md) |
| `lps_strength_energy_decoupled_softening_feasibility_v1` / `20261013T181824Z` | [中文报告](../verification/lps_strength_energy_decoupled_softening_feasibility_v1/20261013T181824Z/REVIEW_LPS_STRENGTH_ENERGY_DECOUPLED_SOFTENING_V1.md) |
| `lps_variational_softening_native_isolated_v1` / `20261013T201800Z` | [中文报告](../verification/lps_variational_softening_native_isolated_v1/20261013T201800Z/REVIEW_LPS_VARIATIONAL_SOFTENING_NATIVE_ISOLATED_V1.md) |
| `lps_native_softening_dynamic_commit_and_small_crack_v1` / `20261013T214800Z` | [中文报告](../verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/REVIEW_LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1.md) |
| `lps_native_softening_full_protocol_timestep_admission_v1` / `20261013T224838Z` | [中文报告](../verification/lps_native_softening_full_protocol_timestep_admission_v1/20261013T224838Z/REVIEW_LPS_NATIVE_SOFTENING_FULL_PROTOCOL_TIMESTEP_ADMISSION_V1.md) |
| `lps_native_softening_event_resolved_integration_v1` / `20261013T235945Z` | [中文报告](../verification/lps_native_softening_event_resolved_integration_v1/20261013T235945Z/REVIEW_LPS_NATIVE_SOFTENING_EVENT_RESOLVED_INTEGRATION_V1.md) |
| `lps_native_softening_domain_retention_and_failure_safe_output_v1` / `20261014T150825Z` | [中文报告](../verification/lps_native_softening_domain_retention_and_failure_safe_output_v1/20261014T150825Z/REVIEW_LPS_NATIVE_SOFTENING_DOMAIN_RETENTION_AND_FAILURE_SAFE_OUTPUT_V1.md) |
| `lps_native_softening_small_nonaffine_crack_validation_v1` / `20261014T210608Z` | [中文报告](../verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T210608Z/REVIEW_LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1.md) |
| `lps_native_softening_small_nonaffine_crack_validation_v1` / `20261014T214309Z` | [中文报告](../verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T214309Z/REVIEW_LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1.md) |
| `lps_native_softening_small_nonaffine_crack_validation_v1` / `20261014T223001Z` | [中文报告](../verification/lps_native_softening_small_nonaffine_crack_validation_v1/20261014T223001Z/REVIEW_LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1.md) |
| `lps_variational_softening_v3_integration_and_single_pilot_v1` / `20261014T232612Z` | [中文报告](../verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/REVIEW_LPS_VARIATIONAL_SOFTENING_V3_INTEGRATION_AND_SINGLE_PILOT_V1.md) |

早期资料入口保留在 [2026-10-11 快照](../github_sync_20261011/README.md)。其中的“当前/下一步”属于当时历史阶段，以本页和最新开发入口为准。

## 文件范围与校验

[来源与导出哈希](SOURCE_FILES.json) · [本次快照清单](SNAPSHOT_MANIFEST.json) · [排除汇总](EXCLUSION_SUMMARY.json) · [逐文件排除清单（gzip）](EXCLUDED_FILES.json.gz) · [导出转换记录](EXPORT_TRANSFORMS.json) · [校验结果](VALIDATION.json)

保留自定义 C++/Python、配置与材料卡、原始 ID/参考拓扑、关键 FP64 状态、完整标量曲线、损伤记录、必要日志、失败诊断、原生版本/构建记录及正式 PNG/SVG。排除约 4.17 GB 的重复上游/构建树、库与对象、密集流式数组/重启、优化器逐次检查点、草稿和论文全文。**原始服务器与本地科研工作区没有删除数据；这是精选科研快照，不是完整服务器镜像。** 没有改写历史 Git 提交。

大于 1 MiB 的 CSV、日志、JSON/JSONL 等文本以无损 gzip 保存，记录原始字节 SHA、导出 SHA 和还原 SHA；不做数值舍入或下采样。全部新文件小于 32 MiB。历史完整运行 manifest 保持原字节，仍描述包含未上传原生文件的完整环境；请用本次 `SNAPSHOT_MANIFEST.json` 核验 GitHub 快照。

在 `clot_clearance` 目录执行：

```bash
python github_sync_20261015/verify_snapshot.py
python github_sync_20261015/restore_compressed.py --verify-only
# 分析脚本需要未压缩的原路径时：
python github_sync_20261015/restore_compressed.py
```

还原器拒绝覆盖不同内容。历史脚本保留当时的绝对路径与原始库 SHA；在其他机器计算必须重新建立独立环境、按锁定上游与补丁重建库并适配新目录。部分完整数值重放仍需取回排除清单中的原生数据。保留旧 `portable_paths.py` 的 GitHub 依赖/字体定位适配，没有用服务器专用版本覆盖它。

最新材料卡与结果已冻结。下一项建议仅为只读追踪 0.376 ms 附着区域自接触/局部塌陷及其附近的失效顺序，不由本次同步启动任何新实验或自动调大载荷。
