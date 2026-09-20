# FEM_SimVascular — Stage SV1

在 WSL 原生运行官方 SimVascular：Stage 1.7 tagged surface → PolyData model →
surface + tetra volume mesh → native fluid solver → VTU/VTP → Python QC 与可视化。

唯一工作目录：`/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular`。
旧 `/home/lzy/projects/formal_3D_flow_solver/FEM` 始终只读。

现有官方 Linux distribution 保留在 `external/SimVascularDistribution/`，
独立 Python 后处理环境为 `.venv/`。本工程只使用一条原生求解路线。
输入与物理量均采用 SI；总流量来自冻结的数值参考条件，不是实验泵流量。

清理记录见 `reports/sv1/cleanup_summary.json`；上一阶段报告保留于
`reports/sv0/REPORT.md`，其中原始附件已按本阶段要求精简。
新阶段结果、日志和报告分别写入 `outputs/sv1/`、`logs/sv1/`、`reports/sv1/`。

当前数值验收状态以 [SV1 报告](reports/sv1/REPORT.md) 为准。保存的失败瞬态场
仅用于诊断，不能作为稳态参考结果。测试会如实报告未通过的数值验收门限。

已有产物的复核命令（不启动求解或重新生成网格）：

```bash
.venv/bin/python -B scripts/postprocess_flow.py
.venv/bin/python -B scripts/reload_solution.py
.venv/bin/python -B scripts/flow_figures.py
.venv/bin/python -B scripts/audit_cleanup.py
.venv/bin/python -B scripts/check_old_fem.py
.venv/bin/python -B -m pytest -q --junitxml=reports/sv1/pytest_results.xml
.venv/bin/python -B scripts/write_report.py
```

运行配置固定于 `configs/sv_reference.yaml`、`configs/mesh_policy.json`、
`configs/time_policy.json` 和 `configs/sv_flow.xml`，其运行前哈希记录在报告目录。
原生构建的源版本、依赖及完整命令见 `native_solver.json` 与 `native_dependencies.json`。
只读旧工程审计包含 `.git` 与原有未提交改动；压缩清单可用标准 gzip 解压读取。
真实 VTU、检查点和完整日志保留在本地 `outputs/sv1/` 与 `logs/sv1/`，不纳入 Git。
