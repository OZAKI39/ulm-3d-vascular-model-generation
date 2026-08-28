# CFD 生产入口瘦身审计

## 正式可达路径

- `cfd_preprocess.py`：YAML 配置 → 已保存 ROI/全局模型 → 全局一维流动 → ROI 边界传递 → 几何引用与 readiness 校验 → 边界条件、端口几何和运行摘要。
- `cfd_surface_prepare.py`：已保存 PASS 输入 → 局部切口 → 官方 VMTK TPS `BOUNDARY_NORMAL` 延长 → cross-seam 实体分配与局部重网格 → OPEN 全量 QC → 封口 → final QC、边界映射、半径、压力修正和米制副本。

## 下游兼容性

Surface 入口仍读取 preprocess 的 `qc/run_summary.json`、`roi/boundary_conditions.json`、`input/geometry_reference.json`、`roi/port_classification.csv`、`roi/port_extension_plan.csv` 和 `roi/port_planes.vtp`。这些正式产物及其所需摘要字段均保留。

## 已删除的开发期内容

- CLI 选项：`--resume-crossseam-open`；两个正式入口现在均只接受一个可选 YAML 路径。
- 状态码：从旧 surface CLI 白名单删除 36 个不可达名称，主要属于 cap-only/direct-cap、旧 entity-only、guarded/proximal-guard、恢复入口及旧 raw/global 实验。
- preprocess 重复文件日志、长篇运行报告、重复流量合计和固定 `false` 审计字段。
- surface 的恢复 CLI、全局重网格、cap-only、旧 entity-only、guarded/proximal-guard、centerline、历史 run 对比和旧自定义 extension/refinement/junction 路径。
- `seam_quality_comparison.csv`、历史对比图和固定历史 run ID 依赖。
- 仅维护旧自定义几何实现的测试和配置键。

## 保留的生产证据

- preprocess：全局一维表格/VTP 与 QC、端口传递与 readiness QC、边界条件、extension plan、port planes、几何引用、输入清单、运行摘要和当前结果图。
- surface：原始/OPEN/remeshed/final STL/VTP、边界 STL 与映射、官方 VMTK 子进程请求/环境/日志、当前 cross-seam 与 final QC、压力修正、米制副本、`crossseam_interface_closeups.png` 和 `final_surface_review.png`。

## 行为冻结结论

- 科学参数修改：否。
- 几何算法修改：否；仅将已经验证的 original-side collar P95、三站截面保真和 final cap/QC 合并到唯一主路径。
- 本次审计运行 VMTK：否。
- 本次审计运行真实 ROI：否。
- 删除历史输出目录：否。
