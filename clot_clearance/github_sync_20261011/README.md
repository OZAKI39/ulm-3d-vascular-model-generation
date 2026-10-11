# 当前 clot 科研快照（2026-10-11）

本分支保存自上次迁移同步以来的代码、输入配置、正式图表、关键数值数组、完整曲线与日志，以及通过和未通过的科学审查结果。当前成果入口如下；旧目录仍保存当时的阶段结论。

## 当前结论

- **含栓 CFD → 材料点载荷**：已完成离线保力、保力矩映射。1927 个材料点中 580 个受载，1561 个固定锚点排除。归一化代数误差约 10⁻¹⁶；压力、剪切、总项 L2 修正分别约 0.618%、9.015%、0.772%。通过本轮局部工程筛查，尚未接入生产碎裂，也没有证明逐点网格收敛。
- **LPS 表面修正与失效**：M5C 的变分诊断通过力—能量一致性检查；V3 出现 15 个负权重点、59 条非正平均权重键，不能启用该修正。逐相互作用路径功可计算，但不能据此移植独立逐键弹性势能或已校准能量阈值。生产瞬断规则不变，Gc 仍未识别。
- **历史轨迹**：原点云断裂及碎片输运仍采用无栓背景 CFD 和工程加载；本次新含栓解没有追溯替换旧轨迹。当前没有双向 FSI、微泡 streaming 或新增 AI 训练结果。

[科研说明（GitHub 相对链接版）](clot_rupture.md) · [WSL 原文快照](clot_rupture.original.md) · [最新共同审查](../verification/y_clot_traction_mapping_and_lps_failure_feasibility_v1/20261011T184653Z/COMBINED_REVIEW.md) · [最新分项状态](../verification/y_clot_traction_mapping_and_lps_failure_feasibility_v1/20261011T184653Z/FINAL_STATUS.json)

## 科研过程与证据

| 阶段 | 报告/静态索引 |
|---|---|
| 宏观材料与时间步审查 | [索引](../results/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z/OPEN_RESULTS.html) |
| 剪切边界与缺口定位 V2 | [索引](../verification/lammps_shear_boundary_and_notch_localization_v2/20261011T001536Z/OPEN_RESULTS.html) |
| 夹持传力与缺口断裂 V3 | [索引](../verification/lammps_notch_fracture_grip_transfer_v3/20261011T115626Z/OPEN_RESULTS.html) |
| PVCE 能量释放 V1 | [索引](../verification/lps_pvce_energy_release_feasibility_v1/20261011T132103Z/OPEN_RESULTS.html) |
| 连续体与近起裂 V2 | [索引](../verification/lps_pvce_continuum_and_near_onset_audit_v2/20261011T142003Z/OPEN_RESULTS.html) |
| 宽薄 PVCE 与真正含栓 CFD | [两线索引](../verification/y_clot_cfd_and_wide_pvce_dual_track_v1/20261011T154854Z/OPEN_RESULTS.html) |
| 保守映射 B 线 | [中文报告](../verification/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/REVIEW_Y_CLOT_CFD_TO_LAMMPS_CONSERVATIVE_TRACTION_MAP_V1.md) · [索引](../verification/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/OPEN_RESULTS.html) |
| 表面修正与能量失效 A 线 | [中文报告](../verification/lps_surface_correction_and_energy_failure_feasibility_v1/20261011T184653Z/REVIEW_LPS_SURFACE_CORRECTION_AND_ENERGY_FAILURE_FEASIBILITY_V1.md) · [索引](../verification/lps_surface_correction_and_energy_failure_feasibility_v1/20261011T184653Z/OPEN_RESULTS.html) |

GitHub 可直接查看 Markdown、PNG 和 SVG。HTML 索引请下载后用浏览器打开；或在 `clot_clearance` 目录运行 `python -m http.server 8000 --bind 127.0.0.1`。

## 最新七组图与数据

| 组 | PNG |
|---|---|
| B1 实际几何与局部支持 | [查看](../results/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/B1_geometry_and_support.png) |
| B2 原 CFD 压力、剪切与法向 | [查看](../results/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/B2_native_surface_tractions.png) |
| B3 守恒与局部修正 | [查看](../results/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/B3_conservation_and_local_quality.png) |
| B4 材料点分项载荷方向 | [查看](../results/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/B4_material_point_loading.png) |
| A1 真实邻域 | [查看](../results/lps_surface_correction_and_energy_failure_feasibility_v1/20261011T184653Z/A1_actual_neighborhoods.png) |
| A2 权重与变分诊断 | [查看](../results/lps_surface_correction_and_energy_failure_feasibility_v1/20261011T184653Z/A2_weights_and_variational_diagnostic.png) |
| A3 路径功与失效限制 | [查看](../results/lps_surface_correction_and_energy_failure_feasibility_v1/20261011T184653Z/A3_path_work_and_failure_limits.png) |

[节点力 N](../verification/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/CLOT_NODAL_TRACTION_FORCES.npz) · [ID/单位/来源](../verification/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/CLOT_NODAL_TRACTION_METADATA.json) · [稀疏映射与约束修正](../verification/y_clot_cfd_to_lammps_conservative_traction_map_v1/20261011T184653Z/CFD_TO_LAMMPS_MAPPING_MATRIX.npz)

## 同步范围与文件核对

[SOURCE_FILES.json](SOURCE_FILES.json) 记录逐文件源/导出 SHA、大小及无损压缩方式；[EXCLUDED_FILES.json](EXCLUDED_FILES.json) 记录未上传项。[EXPORT_TRANSFORMS.json](EXPORT_TRANSFORMS.json) 记录导出副本的凭据清理，原始服务器与 WSL 科研数据不改。[SNAPSHOT_MANIFEST.json](SNAPSHOT_MANIFEST.json) 核对本次新增和更新的 Git 文件，历史科研清单继续代表原来的完整运行环境。

保留原精度的 NPZ 状态、CFD 表面牵引/紧凑网格、节点映射、CSV 曲线、求解日志、失败诊断和所有正式 PNG/SVG。超过 1 MiB 的 CSV/日志以 `.csv.gz` / `.log.gz` 保存，未下采样或舍入。不会把压缩文件误当作原始文件的字节哈希。

完整原生 VTU/VTP、重复的 body.data、逐节点大审计/重启、环境、共享库、重复传输包和论文全文留在服务器或现有本地论文目录。原生库须按已有构建记录重建；部分完整复现需要另取未上传的原生 CFD 文件。论文全文不随公开仓库分发；DOI、数学审查和论文文件 SHA 保留。旧提交中的既有大文件不重写历史。

在 `clot_clearance` 目录执行：

```bash
python github_sync_20261011/verify_snapshot.py
python github_sync_20261011/restore_compressed.py --verify-only
# 需要原路径的完整 CSV/日志时：
python github_sync_20261011/restore_compressed.py
```

还原器核对哈希，拒绝覆盖不同内容。重新计算应使用新输出目录和独立环境，参考各轮 `REPRODUCE.txt`、环境快照及 [原迁移构建记录](../github_sync_20261010/README.md)。原科学脚本保留运行时的服务器路径；本次同步没有重跑仿真或将其改成跨机器完整环境镜像。

下一轮唯一优先建议：现有 V3 完整态、禁止断键的原生静态受载接口验证，核对单位、合力、约束反力/力矩、小变形响应和重复计力。本次仅同步，不自动启动。
