# BraVa 四端口算例的执行与复现

本工作区独立于原打印几何及小鼠结果。`cases/` 中可能包含失败诊断；只有 `reports/ACTIVE_FLOW.json` 指定且独立检查通过的结果可以进入正式展示或微泡批次。已存在求解日志的算例不允许原地覆盖。以下命令不表示应重复启动已完成的 CFD。

## 环境

服务器别名 `vast4090`；工作区 `/workspace/brava_flow_roi_18mlmin_20260928`。

- 科学 Python：`/root/particle8_2_runs/env/bin/python`。
- CUDA FP64 检查：`/venv/main/bin/python`，需要 PyTorch CUDA。
- CFD：`gpu_solver_fix/build.json` 记录本次独立可执行文件、编译/链接命令和源文件哈希。`scripts/build_gpu_sync_fix.py` 保留具体构建过程，不覆盖原求解器。
- 原始微泡模块：`/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms`；保护清单有114个文件。粒径分布使用其中的原 SonoVue 数据。
- 本地科学 Python：`/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`；源微泡项目 `/home/lzy/projects/ulm_particle_formal_p9a5`。
- CPU cgroup 配额约7.68核，采用8个 MPI 进程或8个轨迹进程；OMP/OpenBLAS/MKL 单线程。EGL 真实 NVIDIA 渲染，libx264 CPU 编码。

## 本次已保存的步骤

1. `scripts/prepare_geometry.py`：复制已验收四端口芯，读取制造端口身份，保存 mm→m 和打印姿态变换。
2. `scripts/generate_mesh.py`：官方 SimVascular/TetGen 0.2 mm 网格。原始生成日志在 `mesh/generation/`。
3. `scripts/audit_mesh.py`：闭合性、连接、正体积、单元质量、端口和源几何距离；生成 `mesh/SV_MESH/`。
4. `scripts/design_boundary.py`：中心线阻力及管长近似初步设计，属于辅助估计，不能代替3D分流。
5. `scripts/prepare_case.py`：基础配置。`cases/calibration_gpu8_production/` 保存真正执行的迭代预算和 CUDA8 选项；早期模板或失败算例不能替代该配置。
6. `scripts/validate_gpu_backend.py`：已有 CPU8 与独立 GPU8 两步结果的逐场/WSS比较，不启动新的 CFD。
7. `scripts/run_flow_sequence.py`：校准→独立检查→提取出口压力→全压力出口正式求解→独立检查。此脚本只用于初次顺序执行；现有完整算例不允许覆盖重跑。若需要新实验，应创建独立算例目录并冻结其输入哈希。

## 已有正式场的后处理与媒体

在服务器工作区运行；需要 `ACTIVE_FLOW.json` 和对应冻结数据：

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow
/root/particle8_2_runs/env/bin/python -B scripts/export_print_pose_fields.py
/root/particle8_2_runs/env/bin/python -B scripts/prepare_flow_visualization.py
/root/particle8_2_runs/env/bin/python -B scripts/render_flow_parallel.py --stills-only
# 检查两姿态预览、出口标签和色标后：
/root/particle8_2_runs/env/bin/python -B scripts/render_flow_parallel.py
/root/particle8_2_runs/env/bin/python -B scripts/check_internal_sections.py
/root/particle8_2_runs/env/bin/python -B scripts/check_internal_divergence.py
```

流线读取正式 P1 FEM 场，原始粒子库的 native RK23 仅用于流线绘制；不把它当成有限尺寸微泡积分。9条路径进行空间步长减半核对。每姿态5段动画，真实字段：流线速度、局部速度矢量、压力、节点显示 WSS、原始壁面三角形 WSS。媒体哈希、帧数、色标和 GPU renderer 记录在各自 `MEDIA_VALIDATION.json`。

## 新批次微泡

`scripts/run_microbubble_sequence.py` 首先验证两姿态的全部流场动画已完成并匹配正式流场哈希，然后依次：入口采样→CUDA 梯度检查→原 Python 完整单条轨迹→原/编译几何核及完整轨迹一致性→1500条积分→CUDA 轨迹检查→CSV/报告导出。

```bash
/root/particle8_2_runs/env/bin/python -B scripts/run_microbubble_sequence.py
BRAVA_POSE=0 /root/particle8_2_runs/env/bin/python -B microbubble/scripts/render_results.py --preview
BRAVA_POSE=15 /root/particle8_2_runs/env/bin/python -B microbubble/scripts/render_results.py --preview
# 检查图像后生成完整动画：
BRAVA_POSE=0 /root/particle8_2_runs/env/bin/python -B microbubble/scripts/render_results.py
BRAVA_POSE=15 /root/particle8_2_runs/env/bin/python -B microbubble/scripts/render_results.py
/root/particle8_2_runs/env/bin/python -B scripts/seal_microbubble_delivery.py
/root/particle8_2_runs/env/bin/python -B microbubble/scripts/verify_collected.py
/root/particle8_2_runs/env/bin/python -B scripts/write_results_report.py
```

`config.json` 内的源路径是服务器路径；本地完整性核查使用已记录的本地保护源目录回退，不改变冻结身份。正式轨迹的动态事件及近壁审计在 `microbubble/tracks/mb_*/`；每条 `COMPLETE.json` 验证本条文件及出生事件哈希。

## 字段与单位

- `steady_flow.vtu`：原生 double 节点 `Velocity` (m/s)、`Pressure` (Pa)，坐标 m。
- `flow_arrays_si.npz`：原始规范四面体连接、边界身份及相同节点场。
- `wall_wss_si.vtp`：原始面片 `WSS_raw_Pa`、切向应力矢量、单位外法向、父四面体编号；另有面积加权节点显示 `WSS_display_Pa`。
- `trajectories_dt0p5ms.npz`：轨迹年龄、出生后物理时间、SI位置；实际终点保留其真实出流时间。内部接受子步另存于每条原始轨迹。
- 打印姿态只影响显示副本；不修改原始 SI 物理数据。没有重力的当前模型可以复用两姿态，实际含重力的灌注实验不在此假设内。
