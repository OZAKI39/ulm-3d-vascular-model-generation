# 本轮复现与数据说明

本目录是科学审核，不是新的 P9-A.1 production 数据。全部原始 formal 文件、admission、FEM、physics 保持只读；保护清单见 `data/protected_before.json`。默认不需要重新跑任何积分，可直接用已有证据重绘。

## 只读分析、重绘与测试

在 `/home/lzy/projects/ulm_particle_3d_particle0` 下：

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=particle_3d/src
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python particle_3d/scripts/audit_point_stagnation_field.py
.venv/bin/python particle_3d/scripts/verify_local_radius_passage.py
.venv/bin/python particle_3d/scripts/analyze_routing_stationary_audit.py
.venv/bin/python particle_3d/scripts/render_routing_stationary_audit.py
.venv/bin/python -m pytest -q particle_3d/tests/particle9a1_audit --junitxml=particle_3d/reports/particle9a1_routing_stationary_audit/logs/audit_tests.xml
.venv/bin/python particle_3d/scripts/finalize_routing_stationary_audit.py
.venv/bin/python particle_3d/scripts/write_routing_stationary_review.py
```

8 张 PNG 同时提供 PDF；每张图的数据文件在 `data/figure_manifest.json` 中逐项列出。`OPEN_RESULTS.html` 首先展示 Figure 08 总结图。CSV 中向量/嵌套列是 JSON 字符串；无出口永远保留在分母。

## 计算 provenance 与独立目录

实际服务器根目录见 `data/remote.json`，本轮为 `/workspace/particle9a1_routing_stationary_audit_20260924T064607Z`。服务器 Python 为 `/root/particle8_2_runs/env/bin/python`，实际 CPU quota 为 7.68。P6.5 使用 6 workers；另一个单 worker 依次处理 point、局部 radius、point 分辨率复核，计算并发峰值 7，BLAS 单线程。无需 GPU。

`reference/birth_ledger.json` 是正式 500 的原事件；`reference/all_scheduled_admission_records.json` 是产生它们的全部 834 候选；`reference/production_source_manifest.json` 固定生产源码。运行器启动时校验所有生产 Python 文件及正式 flow SHA。构建新的独立审计目录时，应复制同版生产源码和 Frozen reference，补充原 `particle82_point_native.cpp`，再放入本轮新增审核工具；不要将输出写入旧 formal 目录。

只有确实需要复算对照时，才在独立的 `/workspace/particle9a1_routing_stationary_audit_<timestamp>` 根目录执行：

```bash
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=particle_3d/src
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/root/particle8_2_runs/env/bin/python particle_3d/scripts/run_routing_stationary_audit.py p65 --workers 6
/root/particle8_2_runs/env/bin/python particle_3d/scripts/run_routing_stationary_audit.py point --workers 1
/root/particle8_2_runs/env/bin/python particle_3d/scripts/run_stationary_radius_audit.py --workers 1
/root/particle8_2_runs/env/bin/python particle_3d/scripts/check_point_routing_resolution.py
```

不要重新运行 P9-A.1 正式 500。`run_stationary_radius_audit.py` 明确读取旧 formal 只读参考路径，从已保存 pre-jam 状态进行 `DIAGNOSTIC_ONLY` 局部续算，不调用 admission、不修改正式 event，也不产生新正式轨迹。半径比只按 1、0.99、0.95、0.90 搜索并在局部通过后停止。

## 局部通过判据的审核修订

本轮初版局部诊断采用原卡点 FEM 方向的固定投影。少量弯曲路径已经远离原接触，但固定投影仍不足 3 个半径。保留全部原始计算，使用 `verify_local_radius_passage.py` 对保存路径按局部 FEM 方向累计前进、离开 3r 球域、连续移动及原接触间隙复核，没有重新积分。最终 `stationary_radius_audit.py` 已直接使用该判据以便复现时及时停止。

旧判据结果：`data/stationary_audit_fixed_direction_reference.json`；新判据：`data/radius_passage_policy.json`；超过最终最早通过比值、已在初版计算过的有限诊断记录保留在 supplementary 字段及 NPZ 中，不混进正式路径。

## 时间与哈希

`*_config.json`、`*_completed.json`、逐轨迹 receipt 和日志保存服务器、workers、源码/流场 SHA 及实际耗时。Point 分辨率复核复用原 P1 RK23 tracer，未替换采样方法。原轨迹 NO_EXIT=13（其中 accepted=2）通过解析壁面单元证据解释；这不是将其事后分配至某出口。

Git 工作区在本轮开始前已有 staged 与 unstaged 修改。`logs/git_before.json` 保存原状态；`logs/git_diff_stat.txt` 是整个工作区 diff，新增未跟踪审核文件另列于 `logs/new_audit_files.txt`。本轮未 commit/push。
