# 复现本轮截面门槛

当前状态是 `NO_VALID_INTERIOR_INJECTION_SECTION`。只有前置截面搜索/审核已实现；以下命令不会运行 2000 事件或任何微泡轨迹。

在仓库根目录（WSL）执行：

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=particle_3d/src
.venv/bin/python -m pytest particle_3d/tests/particle9a3_inlet -q
.venv/bin/python particle_3d/scripts/run_particle9a3_section_gate.py --workers 6
.venv/bin/python particle_3d/scripts/render_particle9a3_section_gate.py
```

服务器使用 `/root/particle8_2_runs/env/bin/python` 替换 `.venv/bin/python`。只在新的 P9-A.3 目录执行；搜索输出含实测 host/runtime，重新运行会更新该阶段的输出，不会覆盖旧 B/C 结果。审核报告生成器还需要保存的 `test_results.json` 和 `protected_verification.json`，它不会自行声称测试通过。

源几何追溯在 reference/source_contract.json 和 roi_core.swc。2D 参考 blob 与 SHA 在 data/two_d_source_identity.json。输入校验调用原 read_frozen，候选相对流量容差直接读取冻结 policy.mass_limit；不能用修改阈值的方式重现 PASS。
