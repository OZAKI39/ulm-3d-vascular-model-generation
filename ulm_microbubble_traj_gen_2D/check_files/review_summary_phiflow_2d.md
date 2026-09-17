# Review Summary - Face-Flux 2D Lumen Flow Update

## 1. 本次实现目标

本次目标是把 2D lumen flow 从 cell-centered 伪投影改为 face-flux 有限体积投影，并让入口/出口、散度、压力矩阵、最终通量诊断使用同一套离散 face flux。

固定要求的实现状态：

1. 非物理合格流场仍通过 `FlowConvergenceError` 中止，不会进入 Taichi 粒子推进。
2. solver 身份改为 `phiflow_viscous_fv_projection_2d`；metadata 明确 PhiFlow 只做 explicit diffusion，压力投影是 `scipy_face_flux_finite_volume_cg`。
3. 入口/出口不再是 lumen 内部厚速度带，改为 mask 边界上的 open edge。
4. solid wall face 通量在投影系统中固定为 0；不再执行投影后的 wall cleanup。
5. divergence、pressure matrix、final diagnostics 均来自最终 face flux。
6. 每个 outlet 单独通过 open edge 积分检查 `u dot n * edge_length`。
7. 最终保存的 velocity 是 `face_flux_after_projection` 重建的 cell-centered velocity。

## 2. 修改了哪些文件，每个文件为什么改

| 文件 | 原因 |
|---|---|
| `ulm_microbubble_traj_gen/utils/vessel_rasterizer.py` | 使用 Shapely `LineString.buffer` / junction disk / `unary_union` 生成 lumen polygon，再用 Rasterio `rasterize` 烧录 lumen mask；保留后续 nearest vessel segment 属性填充。 |
| `ulm_microbubble_traj_gen/utils/flow_boundaries.py` | 重写入口/出口构造为 mask boundary open edge；新增 edge normal、edge length、face index、signed open flux、per-label target。 |
| `ulm_microbubble_traj_gen/utils/face_flux_projection.py` | 新增 face-flux 有限体积投影：fluid-fluid face 可压力修正，wall face 固定 0，open face 固定目标通量。 |
| `ulm_microbubble_traj_gen/utils/phiflow_solver.py` | 主求解链改为 PhiFlow diffusion -> face flux -> open edge flux -> FV projection -> physical gate；移除投影后壁面清理和入口/出口速度硬覆盖。 |
| `ulm_microbubble_traj_gen/utils/types.py` | `FlowField` 增加 face flux、open boundary flux、edge length、per-label actual flux 字段。 |
| `ulm_microbubble_traj_gen/utils/field_io.py` | 将新增 face-flux / open-edge 诊断数组写入 NPZ。 |
| `ulm_microbubble_traj_gen/utils/flow_diagnostics.py` | outlet CSV 改用最终 face flux 积分后的 per-label actual flux；stage divergence 改为真实执行阶段。 |
| `ulm_microbubble_traj_gen/utils/config.py` | 默认 solver mode 改为 `phiflow_viscous_fv_projection_2d`。 |
| `ulm_microbubble_traj_gen/configs/physics_flow_config.yaml` | 同步默认 solver mode。 |
| `ulm_microbubble_traj_gen/test_files/test_phiflow_2d_validation.py` | 测试从旧 cell velocity flux correction 改为 open-edge + face-flux projection 独立验证。 |
| `ulm_microbubble_traj_gen/check_files/review_summary_phiflow_2d.md` | 更新本次实现记录。 |

没有修改 `ulm_vascular_model_generator/vessel_generation.py`：该文件只是 DCCO 生成器 CLI/API re-export。lumen mask 和流场求解发生在 microbubble pipeline 中，DCCO graph 数据通过导出的 `Vessel` 列表被 downstream 模块读取。

## 3. 核心逻辑变化，按调用链说明

```text
generate_microbubble_trajectories.py
-> runner.py::run_generation
-> rasterize_vessels(...)
   -> Shapely buffered vessel corridors + junction disks
   -> unary_union / validity repair
   -> Rasterio rasterize -> lumen_mask
-> solve_velocity_field(...)
   -> build_flux_boundaries(...)
      -> root/terminal Vessel -> mask boundary open edges
      -> edge normal / length / face index / target signed flux
   -> build_initial_velocity_field(...)
   -> PhiFlow explicit diffusion
   -> cell_velocity_to_face_flux(...)
   -> apply_boundaries_to_face_flux(...)
   -> project_face_flux(...)
      -> same face flux builds divergence RHS
      -> same fluid-fluid faces build pressure matrix
      -> same face flux is corrected and diagnosed
   -> physical convergence gate
-> write_flow_diagnostics(...)
-> advect_particles(...)
```

## 4. 行为变化：改动前 vs 改动后

| 项目 | 改动前 | 改动后 |
|---|---|---|
| lumen mask | 基于逐段距离和形态学平滑 | Shapely corridor polygon + junction disk + unary_union + Rasterio rasterize |
| 边界 | lumen 内部薄/厚 boundary cells 上的速度值 | mask 边界 open edge 上的固定 face flux |
| 投影变量 | cell-centered velocity divergence | face flux finite-volume divergence |
| 压力矩阵 | cell 邻接近似投影 | fluid-fluid internal face graph Laplacian |
| wall no-penetration | 投影后尝试清理壁面法向速度 | wall face 从一开始固定 0，最终无后处理 |
| flux diagnostic | cell velocity dot normal 估算 | 最终 face flux 上逐 edge 积分 |
| outlet 检查 | 总量为主 | 每个 outlet label 单独检查 |
| 最终保存场 | 可能经过投影后再修改 | `face_flux_after_projection`，metadata 标记 `post_projection_velocity_modified=false` |

## 5. 新增/修改了哪些测试

修改 `ulm_microbubble_traj_gen/test_files/test_phiflow_2d_validation.py`，当前覆盖 4 个独立测试：

1. open boundaries 必须位于 mask edge，并且 face flux 积分匹配 target。
2. face-flux projection 使用同一套 flux 做 divergence 和 diagnostics，投影后 divergence 低于阈值，wall penetration 为 0。
3. in-plane wall shear proxy 仍匹配 2D Poiseuille 量级。
4. `physical_converged=False` 的 flow 会被 Taichi 粒子推进拒绝。

## 6. 已运行的检查命令及结果

```powershell
d:\anaconda3\envs\pmp\python.exe -m compileall ulm_microbubble_traj_gen\utils ulm_microbubble_traj_gen\test_files
```

结果：通过。

```powershell
d:\anaconda3\envs\pmp\python.exe -m unittest discover -s ulm_microbubble_traj_gen\test_files -v
```

结果：4 tests，OK。

```powershell
d:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\generate_microbubble_trajectories.py --quick-test
```

结果：通过并生成轨迹，输出目录：

```text
ulm_microbubble_traj_gen/results/20260709_144137
```

关键 quick-test 指标：

```text
solver_mode=phiflow_viscous_fv_projection_2d
saved_velocity_stage=face_flux_after_projection
post_projection_velocity_modified=false
post_projection_wall_cleanup_performed=false
final_normalized_divergence_error=7.114391206247104e-07
interior_normalized_divergence_error=5.617263177266306e-07
final_flux_relative_error=1.031721411169771e-15
final_outlet_flux_max_relative_error=1.8812800491501592e-16
wall_penetration_max_um_s=0.0
pressure_cg_info=0
```

```powershell
rg -n "correct_boundary_flux|initialize_boundary_velocity|_remove_wall_normal_velocity|div_after_wall_cleanup|div_after_final_boundary_apply" ulm_microbubble_traj_gen\utils ulm_microbubble_traj_gen\configs ulm_microbubble_traj_gen\test_files
```

结果：无匹配；旧的后处理/硬校正路径没有残留引用。

## 7. 需要人工重点看的 5 个位置

1. `ulm_microbubble_traj_gen/utils/flow_boundaries.py`
   - open edge 筛选逻辑，尤其是斜向 vessel endpoint 的 edge normal / lateral profile 是否符合预期。
2. `ulm_microbubble_traj_gen/utils/face_flux_projection.py`
   - face flux 正负号约定、pressure correction 符号、solid wall face 排除。
3. `ulm_microbubble_traj_gen/utils/phiflow_solver.py`
   - PhiFlow diffusion 和 FV projection 的连接处，以及 physical gate 使用的指标。
4. `ulm_microbubble_traj_gen/utils/vessel_rasterizer.py`
   - Shapely/Rasterio rasterization 的 grid transform、junction disk 半径和 `all_touched=True` 对细血管宽度的影响。
5. `ulm_microbubble_traj_gen/utils/flow_diagnostics.py`
   - outlet CSV 是否满足你需要的 per-outlet 检查格式；当前 actual 来自最终 face flux，而不是 cell velocity。

## 8. 是否有任何超出需求范围的改动

没有改动 SWC 导出、DCCO growth/hemodynamics、Taichi 粒子积分核心规则、MAT/声场相关模块。

本次没有修改 `ulm_vascular_model_generator/vessel_generation.py`，原因是它不是 lumen mask 或流场计算的执行位置。实际功能改动集中在读取 `Vessel` 列表后的 microbubble trajectory generation pipeline。

运行 quick-test 会新增结果目录 `ulm_microbubble_traj_gen/results/20260709_144137`；这是验证输出，不属于算法代码改动。
