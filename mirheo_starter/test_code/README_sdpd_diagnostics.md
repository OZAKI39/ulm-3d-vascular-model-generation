# SDPD 原因排查

本地根目录 `/home/lzy/projects/mirheo_starter`。所有入口使用 `.venv/bin/python -B`。输入由 `py_scripts/sdpd_diagnostics.yaml` 中的明确路径及 SHA-256 固定，不扫描“最新”结果目录。

```bash
cd /home/lzy/projects/mirheo_starter

# 已有数据、源码及实际探针的 CPU 诊断；不启动 GPU
.venv/bin/python -B -m py_scripts.diagnose_sdpd \
  --config py_scripts/sdpd_diagnostics.yaml --analyze-only

# 离线页面；不启动 Mirheo/GPU 实验
.venv/bin/python -B -m test_code.review_sdpd_diagnostics \
  --config py_scripts/sdpd_diagnostics.yaml --open

# 实际 Windows Chrome / CDP 检查，使用软件绘制；已通过的同哈希页面复用记录
.venv/bin/python -B -m test_code.review_sdpd_diagnostics \
  --config py_scripts/sdpd_diagnostics.yaml --browser-test

# 本次经过筛选的探针已完成；同版本再次调用只校验并复用，绝不重跑
.venv/bin/python -B -m py_scripts.diagnose_sdpd \
  --config py_scripts/sdpd_diagnostics.yaml --execute-approved-probes

# 全部 CPU 回归（合成用例不作为物性证据）
.venv/bin/python -B -m unittest discover -s test_code -p 'test_*.py' -v
```

诊断 CLI 无参数返回用法错误，`--help` 正常退出；两者都不启动 GPU。查看页面默认只读。依赖现有 Python/Plotly、Windows Node 和 Chrome；没有安装包、启动服务或重编译计算库。

## 结果与来源

- 正式结果输出到 `data/sdpd_diagnostics/result_<identity>/`；命令打印确切路径，重复查看复用该包。
- 页面输出到 `test_code/outputs/sdpd_diagnostics/<result>/<viewer_identity>/sdpd_diagnostics_review.html`。
- 每个包有 `package_sha256.json`；页面的 `html_record.json` 绑定包和 HTML。浏览器实测结果位于同目录 `browser_final/browser_checks.json`，绑定 HTML 哈希。人工验收始终 PENDING。
- 原报告修复前精确复现、原文件保护哈希、先失败后通过的测试和 CPU 日志在 `test_code/outputs/sdpd_diagnostics/setup_20260908T221052_971138Z/`。
- 新 GPU 运行位于 `runs/sdpd_diagnostics/thermal_cause_20260909/`。首次观测通道命名错误的失败和修正观测版本分目录保留。两次费用均计入原 `runs/fluid_calibration/shared_budget_pool.json`，不新建授权额度。
- `CPU_validation_attempt01.json`、`CPU_validation_corrected_observer.json`、`CPU_validation_final.json` 保留各次代码哈希和测试记录。`CPU_validation.json` 是当前验证记录；历史版本仍在，不覆盖旧实验报告。
- 最终共 **153 项 CPU 测试通过**，包括缓存配置/库版本失配和八块预算的采样网格取整。当前记录 `CPU_validation_delivery.json`，日志 `CPU_delivery_tests.log`。页面和检查器由单独的实际浏览器记录验证。
- 之前的预算草稿、失败测试、观测初始化失败、首次浏览器默认可见性断言错误及全部验证版本均保留。最终交付路径与验证哈希由 setup 目录的 `delivery_record.json` 明确绑定。

## 解读

`diagnosis_summary.json` 包含实际新探针、通道布局、CPU 前后相位复核、初态与后续真实快照的结构审查；`evidence_matrix.json` 明确区分已证实、支持但未证实、检查范围内排除、证据不足、未测。

`sampling_audit.csv` 记录计划/分配/实际步数、物理时长、截段、ACF、完整块及尾部；`temperature_audit.json` 保留完整轨迹和全部窗口。`profile_audit.json` 保留所有空间点、目标曲线、拟合、残差和加权平均差异。`matched_physical_time_audit.json` 比较相同物理阶段，未把短轨迹拼成稳态数据。

`pressure_audit.json` 将输入 EOS、机械压力、共同参考压差及协方差假设分开。`native_source_contract.json` 定位本地源码链；`time_scale_audit.json` 和 `snapshot_audit.json` 是独立公式推导，不能视为另一套流体求解器的验证。

新探针执行使用 `actual_parameters.json` 顶层 `steps=160000`、`dt_star=5e-7`，以及 `diagnostic_plan` 的 0.08* 固定物理时长。嵌套 `task.desired_time_star=0.25` 是继承的旧任务模板元数据，未被本次 worker 用于推进；以顶层执行值和 `worker_completion.json` 的实际 160000 步为准。该记录保持原样，没有事后改写。

本次对已有正式 Python 代码的修复只在 `py_scripts/fluid_physics/analysis.py`：均值/剖面使用与 CI 相同的完整块，另存全窗均值和剖面。不改原始 CSV、不删除尾部、不改阈值或单位。新增探针的 `__forces` 修正仅影响新观测代码。

没有修复数值模型，也没有选出合格液体；`selection=null`。后续只优先无驱动热平台实验，具体步数、固定观测窗、预计成本和所缺授权见 `proposed_experiments.json → next_priority_experiment`。当前没有追加预算，也不会自动继续运行。
