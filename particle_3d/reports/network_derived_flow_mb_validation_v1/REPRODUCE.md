# 可复现入口

本地原科学包：`/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/src/particle_3d`（只读）。新增 adapter/分析/绘图位于本目录 `scripts/`；永久测试位于 `particle_3d/tests/network_flow_mb_validation_v1/`。

服务器工作区：`/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z`。输入包 SHA：`f2669f198c4a071d04119c1dbed26fd3ff17fe2cbe50431dd4f1be86276e0a94`。该工作区每次运行前逐一验证 `data/code_snapshot_manifest.json`、flow SHA、共享 cohort。所有 30 泡事件固定，不执行 sampler。

已执行的命令形态：
```
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1 TMPDIR=/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/tmp /root/particle8_2_runs/env/bin/python /workspace/particle_network_flow_mb_validation_v1_20260924T210047Z/scripts/runner.py --root /workspace/particle_network_flow_mb_validation_v1_20260924T210047Z --label NEW --workers 6 --ids all --output NEW
```
既有 output 名存在时 runner 主动拒绝覆写。复跑必须使用新的、独立的 workspace/output；不要重跑历史路径。determinism 用 `--ids 3,15,22 --workers 1/3`；point 用 `--mode point --workers 1`。本轮 OLD 通过源码/flow/event/样本 SHA 核验复用，记录在 cohort_provenance.json。

本地只读永久检查：
```
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python -m pytest -q -p no:cacheprovider /home/lzy/projects/ulm_particle_3d_particle0/particle_3d/tests/network_flow_mb_validation_v1
```
本目录 `archive/` 保存原始传输包、Taylor-Hood 最后已完成单步 smoke，以及修正图标签之前的原图。`server_bundle/` 保留科学源码快照与 frozen_reference。所有正式输出均由未经修改的 Particle 科学模块生成。
