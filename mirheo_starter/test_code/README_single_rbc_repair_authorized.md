# 已授权修复对照：实际交付入口

Campaign：`rbc_repair_20260910T131105Z`。

隔离编译和 A0–A6 七项 GPU 对照已经执行。半步长完成零剪切 t*=30、Γ=0；面积与准备残余标准仍失败，没有新的 Γ=4 比较或合格速度比。费用与科学状态以 `data/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_results.json` 为准，浏览器状态以绑定 HTML 哈希的 `delivery_receipt.json` 为准。

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.single_rbc_repair.authorized_result_report --help
.venv/bin/python -B -m py_scripts.single_rbc_repair.authorized_result_report --verify
.venv/bin/python -B -m py_scripts.single_rbc_repair.authorized_result_report --open
```

这些查看命令不运行求解器。`--open` 打开现有 HTML；`--render` 仅从已有分析重新生成页面，重新渲染后浏览器验收必须重新绑定新哈希。

本次公共查看入口是 `authorized_result_report`。原 `repair_single_rbc_benchmark --review` 及 `test_code.review_single_rbc_benchmark` 属于先前冻结的页面生成代码，保留用于旧状态复核，不用于重新生成本次交付页面。原 `--execute` 队列已包含失败，不能用它重试。各次实际求解命令保存在执行记录中，A4/A5/A6 的新计划分别记录连续调用、延长准备、半步积分的改变；没有自动重试授权。

原始结果位于 `runs/single_rbc_repair/rbc_repair_20260910T131105Z/gpu/`，包含 HDF5、CSV、NPZ、原生失败快照、库哈希、执行日志与输出摘要。CPU 转换只写入 `data/.../runtime_reviews/`，没有改写这些原始文件。

重点证据：`run_boundary_audit.json`、`A5_failure_step_geometry.json`、`A6_full_population_membership.json`、`A6_preparation_quality.json`。安装库未被覆盖；前后 20,978 个旧文件哈希通过。旧 GitHub 分支未被本次工作改写，本轮没有提交或推送。
