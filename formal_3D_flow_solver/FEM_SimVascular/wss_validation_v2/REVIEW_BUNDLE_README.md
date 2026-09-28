# 第二轮 WSS 精简审阅包

请先阅读 `WSS_VALIDATION_REPORT_V2.md`。`data/` 是带单位、区域、统计权重和比较基准的机器可读数据；`figures/` 保存5张诊断图的PNG/PDF；`evidence/production_fix.patch` 是本轮最小生产代码差异。`CASE_STATUS.csv` 区分实际合格解、原解复用、主动停止、失败和未执行输入，避免把性能试验当作科学结果。

`evidence/solver_log_excerpts/` 保留带原行号的求解日志摘要和原始日志SHA256；每例原配置、来源锁、独立检查和收敛记录按stage目录排列。完整原始VTU、NPZ、重启文件和原始日志在 `ARTIFACT_INDEX.csv` 所列本地/服务器路径，未塞入此精简包。路径为空的服务器栏不表示服务器副本已核实存在。

可以不访问整个项目，直接用NumPy重做包中J1/J2实际速度子集的生产WSS恢复：

```bash
python -B scripts/reproduce_local_wall_evidence.py
```

脚本逐一核对核心代码和数据SHA256，并与已保存原始面片向量/模长比较，结果写到 `data/local_subset_reproduction.json`。子集保留实际CFD节点值、原坐标和父单元，**不是一个闭合的新CFD域**，也不能代替完整求解或证明物理准确性。完整CFD复现需要原网格/几何、对应求解器和运行环境，见 `COMMANDS.md`；不要覆盖既有算例重跑。

`REVIEW_BUNDLE_MANIFEST.json` 给出包内主要文件的大小和SHA256；外置 `.sha256` 文件给出整个ZIP的SHA256。正式旧图、旧网格和上一轮审计保持原样，最终完整性核对记录在 `evidence/protected_inputs_after.json`。

本轮最终范围：三档圆管完整CFD保留；真实血管仅原网格与中档的两档敏感性分析。细档真实血管CFD及整个出口压力敏感性阶段均因用户限制时间成本取消，状态`SKIPPED_BY_USER`，不是PASS。没有压力响应图或扰动结果CSV；不得根据缺失结果推断不敏感。细档网格几何指标可以保留，但不能作为细档CFD证据。
