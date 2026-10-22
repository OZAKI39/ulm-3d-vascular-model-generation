# 持续灌注工作流复现

这些命令只运行微泡模块/验证/后处理，不启动任何 CFD。复算会创建新输出目录，禁止覆盖已有算例。

## 权威路径

- 本地源码：`/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/src/particle_3d/population_v1/`
- 本地测试：`/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/tests/population_v1/`
- 本地全部新交付：`/home/lzy/projects/brava_flow_roi_18mlmin/population_v1/`
- 服务器 SSH 别名：`vast4090`
- 服务器全部新交付：`/workspace/brava_flow_roi_18mlmin_20260928/population_v1/`
- 服务器源码：上述目录 `source/particle_3d/src/`
- 服务器隔离环境：上述目录 `env/bin/python`，沿用服务器预装 Torch 的 `.pth`，其余科学包只装在该隔离环境。
- 最终生产：`runs/realization_0000_contactfix/`；三个开发中止批次 `realization_0000/`、`realization_0000_guardfix/`、`realization_0000_final/` 均不能作为有效完成数据。
- 最终生产积分源码 commit：`98d705f13babdfa0cba6c12331623292a588cb77`，本地 Git 工作树 `/home/lzy/projects/temp_storage/github_sync_20260928_latest`，分支 `dev/brava-population-gpu-20260929`。后续分析/文档提交不改变该运行身份。

## 在服务器运行测试

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
export PYTHONPATH="$PWD/source/particle_3d/src"
export POPULATION_VALIDATION_DIR="$PWD/validation"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env/bin/python -B -m pytest source/particle_3d/tests/population_v1 -q
```

真实 FEM/壁面/原积分器 parity 在 `scripts/validate_gpu_field.py`、`validate_gpu_wall.py`、`validate_real_reference.py`、`validate_real_gpu_dynamics.py`；真实源统计在 `validate_real_source.py`。这些都是独立微泡验证，不需要 CFD。

## 复算同一实现（仅在需要时执行）

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
export PYTHONPATH="$PWD/source/particle_3d/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env/bin/python -u -B -m particle_3d.population_v1.run \
  --environment data/environment.json \
  --output runs/reproduction_new_directory \
  --backend cuda --concentration 4e9 --seed 2026092909 \
  --minimum-exits 20 --post-infusion 12 --dt 0.0005 \
  --git-commit 98d705f13babdfa0cba6c12331623292a588cb77
```

与该 commit 完全相同的科学源码和相同软件/设备用于位级复现；后续源码版本须记录自己的 commit。输出身份还保存逐文件 SHA，不靠 Git 标签猜测执行版本。`launch_production.py` 只包装上述命令并记录 PID/资源监测，不恢复旧已停止任务。

## 仅分析/重新绘图

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env/bin/python -u scripts/analyze_population.py --source-root source --run-name realization_0000_contactfix
env/bin/python -u scripts/render_population.py --run-name realization_0000_contactfix --pose 0
env/bin/python -u scripts/render_population.py --run-name realization_0000_contactfix --pose 15
env/bin/python scripts/plot_summary.py --run-name realization_0000_contactfix
env/bin/python scripts/final_qc.py --run-name realization_0000_contactfix
env/bin/python scripts/build_report.py --run-name realization_0000_contactfix
```

动画要求有效 NVIDIA EGL；默认 `VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow`。编码固定 libx264 CPU 两线程。相机和所有文本布局由 `style_reference.py` 和 `render_population.py` 保存。

依赖清单为 `validation/server/environment_packages_final.txt`，动画需 `imageio`、`imageio-ffmpeg` 及 `/usr/bin/ffmpeg`。所有新增包仅安装在本任务隔离环境中。下载与本地清单命令：

```bash
cd /home/lzy/projects/brava_flow_roi_18mlmin/population_v1
bash scripts/collect_results.sh
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python scripts/build_report.py --run-name realization_0000_contactfix
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python scripts/inventory_delivery.py --run-name realization_0000_contactfix
```

`data/environment.json` 保存实际服务器路径；`data/environment.local.json` 是只改路径并保留哈希的 WSL 只读适配。WSL 环境无 CUDA Torch，不能把本地 skip 写成 GPU 通过。

## 数据语义与任意物理时刻重建

NPZ 原始状态字段：particle_id/source_event_id/realization_id、birth_time、age、x/v/omega/q、radius/diameter、wall_gap、tetra_id、status。单位 m、s、m/s、rad/s；status=0 表示快照中的活跃泡，终态分类在 `terminals.jsonl.gz` 中以稳定 ID 关联。pair IDs/count/contact 在逐步 `interaction_audit.jsonl.gz` 按物理时间和 ID 关联，无事件的步为零。

`state_time_s[j]` 为全局物理时刻，`state_offsets[j:j+2]` 为该时刻所有活跃粒子的行切片。出生本身也写快照；同一 ID 同时刻重复快照可取最后一个。保存每次接受端点，`v` 为结束在该端点的前一段 held velocity（出生行是初始自由速度）。用真实出生≤t<终态时间筛选，再在相邻端点之间线性插值；不外推、不按 age 把不同出生对齐。精确终态时刻可显示端点，但之后必须消失。

出生行 wall_gap 为单次 admission 的保守下界，tetra_id=-1 是初始未定位 hint。`analysis/.../tracks/` 的派生导出补算出生的精确 gap/tetra；原始状态不改写。接受推进后的 gap 和 tetra 是实际查询结果。

`SOURCE`、`ACCEPTED`、`ACTIVE`、`EXITED`、`SUPPORTED_STATIONARY`、`SOLVER_FAILURE`、`OBSERVATION_CENSORED` 是不同类别，不能统称已完成微泡。粒子中止不代表生理栓塞或溶解。源停止取决于实际出口覆盖，不能将这组计数当作固定时窗的无条件出口概率。

## CPU reference backend

同一 CLI 支持 `--backend cpu`；正式参数全相同时数值 CPU reference 成本很高。轻量执行例使用明确的短窗口验证覆盖，不把它当正式灌注结果：

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
export PYTHONPATH="$PWD/source/particle_3d/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env/bin/python -B -m particle_3d.population_v1.run \
  --environment data/environment.json --output validation_runs/cpu_example_new \
  --backend cpu --concentration 4e9 --seed 2026092909 \
  --minimum-exits 20 --fixed-infusion 0.002 --post-infusion 0.003 --dt 0.0005 \
  --git-commit 98d705f13babdfa0cba6c12331623292a588cb77
```

此例可能没有出口完成；覆盖结果为 false 是正确短窗口语义。`--fixed-infusion` 仅用于明确标记的测试；正式运行不提供它。

## 多独立 realization

`population_v1.gpu_batch.perfuse_batch_until(engines, sources, end_times)` 接收多个 `GPUPopulationEngine` / `PerfusionSource`，每个使用独立 realization_id。共有静态 `GPUField` / `GPUWall` 可复用；严禁不同 realization 互相碰撞。可执行参考为 `tests/population_v1/test_phase8_batch_gpu.py` 和 `scripts/benchmark_sparse_batch.py`。不会额外执行用户没有要求的正式群体。
