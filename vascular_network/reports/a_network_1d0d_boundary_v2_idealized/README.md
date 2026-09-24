# A-Network Idealized Hydraulic Baseline Model v2

本轮采用理想化模型：2410 为主动指定的水力源点，其余 122 个结构叶节点为 p_ref=0 Pa 的参考压力终端；0 Pa 是表压基准。这些定义均不是生理标注，结果不代表真实小鼠脑循环 ground truth。

状态：**A_NETWORK_H0_3D_VALIDATED**。

- [主审阅报告](A_NETWORK_H0_REVIEW_ZH.md)
- [明确的 H0 假设](H0_MODEL_ASSUMPTIONS_ZH.md)
- [ROI 工作点与守恒](ROI_H0_OPERATING_POINT_ZH.md)
- [真实切点到 FEM cap](NETWORK_TO_FEM_PRESSURE_TRANSFER_ZH.md)
- [终端和半径敏感性](TERMINAL_RADIUS_SENSITIVITY_ZH.md)
- [新 3D 独立验证](ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md)
- [机器摘要](data/final_summary.json)
- [测试日志](logs/pytest.log)：49 项，0 个失败/错误（v1 30 + v2 19）。

所有图有 PNG 与 PDF；Figure 06 在真实新 FEM 验证通过后生成。网络计算本地 1 worker，实测 5.845 s，峰值 RSS 143.37 MiB。新 FEM 使用独立服务器目录，不覆盖原 production。

# 复现入口

从 `ulm_3D_vascular` 根目录使用项目已有科学 Python 环境：

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/run_a_network_h0_v2.py
python -m pytest tests/network_1d0d tests/network_h0 -q
python scripts/render_a_h0_v2.py
python scripts/write_a_h0_reports.py
```

网络脚本只写本 v2 目录。若完整重跑，网络脚本会将摘要重置为 3D pending；必须再用已保留的实际 FEM 输出执行 `scripts/validate_a_h0_fem.py --case <新case> --flow-root <FEM_SimVascular根目录>`，再生成图表和报告。

创建新 FEM 输入使用 `network_1d0d.fem_h0_case.create_case(original, fresh_target, bc_json)`，目标存在就拒绝覆盖。服务器启动命令与环境参数保存在 `data/3d_remote_run.json`，runner 为 `scripts/solve_a_h0_fem_remote.py --case <remote_case> --reference-root /workspace/flow_mean_2p0_mmps_20260922`。实际 executable/PETSc SHA、命令、MPI/OMP、完整求解历史记录在新 case 的 reports/execution.json。长期保存本地结果，不依赖服务器永久保留目录。

当前分支 `dev/a-network-1d0d-idealized-baseline-v2`，本轮仅新增代码、测试、报告和独立 flow case；不提交或清理已有大量历史修改。Git 快照与最终 SHA 复验见 logs/。

最终 tracked diff-stat：`1135 files changed, 1701876 insertions(+), 1700623 deletions(-)`。并发工作区变化记录：{'file': '.gitignore', 'added_lines': ['data/external/TopBrain/', 'data/TopBrain'], 'action': 'Preserved without modification; not written by this task'}。保留已有修改，未回退其它工作的 `.gitignore` 规则。文件保护结果见 [最终输入审核](logs/final_input_protection_audit.json)，工作区差异见 [Git 审核](logs/git_workspace_final_audit.json)。
