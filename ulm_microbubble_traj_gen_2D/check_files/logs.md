# 当前流场计算工作流

本文档记录 `ulm_microbubble_traj_gen` 当前版本的流场计算流程、关键配置、收敛判据和最近一次验证结果。

## 1. 当前目标

当前代码的目标不是用经验公式直接给微泡速度，而是先根据导出的血管模型生成一个网格流场，然后让微泡在这个流场中插值运动。

核心思想如下：

```text
血管模型
-> 网格化 lumen 区域
-> 1D 血管阻力网络求压力
-> Stokes / creeping-flow 网格求解
-> 速度场和壁面剪切力场
-> Taichi 微泡轨迹推进
```

当前默认使用：

```yaml
field:
  solver_mode: stokes
```

## 2. 关键调用链

完整调用链如下：

```text
ulm_microbubble_traj_gen/generate_microbubble_trajectories.py
-> utils/runner.py::run_generation
-> utils/vascular_io.py::load_physics_input
-> utils/grid_domain.py::build_domain_from_vessels
-> utils/vessel_rasterizer.py::rasterize_vessels
-> utils/connectivity.py::validate_fluid_connectivity
-> utils/phiflow_solver.py::solve_velocity_field
-> utils/stokes_solver.py::solve_stokes_flow
-> utils/network_coupling.py::solve_network_pressures
-> utils/taichi_particles.py::advect_particles
```

其中最重要的流场计算入口是：

```text
utils/phiflow_solver.py::solve_velocity_field
```

当 `field.solver_mode` 为 `stokes` 时，会进入：

```text
utils/stokes_solver.py::solve_stokes_flow
```

## 3. 输入数据

当前流场计算使用 `ulm_vascular_model_generator` 已经导出的血管模型。

默认配置为：

```yaml
input:
  model_name: ulm_xz_planar_dcco_tree_seed_105
  swc_dir: ulm_vascular_model_generator/vessel_swc_models
```

实际读取的主要文件包括：

| 文件             | 用途                                                     |
| ---------------- | -------------------------------------------------------- |
| `.swc`         | 保存血管树的基础几何结构                                 |
| `.vessels.npz` | 保存更完整的 vessel 信息，例如半径、流量、黏度、父子关系 |

当前流场计算主要依赖 `.vessels.npz` 中的完整 vessel 信息，而不是只依赖 SWC。

## 4. 网格化流程

`rasterize_vessels(...)` 会把每一段血管转换到 X-Z 平面网格上。

生成的关键网格字段包括：

| 字段                    | 含义                       |
| ----------------------- | -------------------------- |
| `lumen_mask`          | 哪些网格单元在血管腔内     |
| `wall_mask`           | 哪些网格单元接近血管壁     |
| `vessel_id`           | 每个网格单元对应哪一段血管 |
| `radius_um`           | 网格单元所在血管的半径     |
| `flow_rate_um3_s`     | 网格单元所在血管的流量     |
| `viscosity_mpas`      | 网格单元所在血管的黏度     |
| `direction_xz`        | 血管在 X-Z 平面上的方向    |
| `distance_to_wall_um` | 网格单元到血管壁的距离     |

当前版本不再用经验速度公式直接生成微泡速度。

## 5. 连通性检查

在求解流场之前，代码会先检查 `lumen_mask` 是否为一个连续的流体区域。

调用位置：

```text
utils/connectivity.py::validate_fluid_connectivity
```

如果网格太粗，血管会被离散成多个断开的区域。此时程序会终止，不会继续求解流场。

典型失败原因：

```text
domain.grid_spacing_um 太大
domain.min_lumen_radius_cells 太小
```

默认完整配置使用：

```yaml
domain:
  grid_spacing_um: 2.0
```

该配置在最近一次验证中通过了连通性检查。

注意：`--quick-test` 使用 50 um 粗网格，可能会导致 lumen 断裂。这是 quick-test 网格太粗导致的，不代表 Stokes 求解器本身失败。

## 6. 1D 血管阻力网络

当前版本会先在血管树上建立一个 1D 阻力网络。

调用位置：

```text
utils/network_coupling.py::solve_network_pressures
```

每一段血管的阻力近似来自 Poiseuille 关系：

```text
R = 8 * mu * L / (pi * r^4)
```

其中：

| 符号   | 含义       |
| ------ | ---------- |
| `R`  | 血管段阻力 |
| `mu` | 动态黏度   |
| `L`  | 血管段长度 |
| `r`  | 血管半径   |

1D 网络使用：

```yaml
network_terminal_pressure_pa: 0.0
network_root_flow_scale: 1.0
```

也就是说：

```text
root 端使用导出血管模型中的 root flow
terminal 端使用固定 terminal pressure
```

最近一次验证中：

```text
network_root_flow_um3_s: 1692540.0679861226
network_pressure_min_pa: 0.0
network_pressure_max_pa: 104.15216240889335
network_terminal_nodes: 35
```

## 7. 1D 压力到网格边界的映射

1D 网络算出 root 和 terminal 的压力后，会映射到 lumen 网格边界。

当前模式：

```yaml
boundary_condition_mode: one_d_pressure
```

对应的边界风格为：

```text
one_d_pressure_with_source_sink_flow_prior
```

实际含义是：

```text
1D root pressure
-> root inlet 网格边界

1D terminal pressure
-> terminal outlet 网格边界

root / terminal flow
-> source/sink 通量软约束
```

最近一次验证中：

```text
root_inlet_boundary_cells: 320
terminal_outlet_boundary_cells: 1601
network_pressure_boundary_cells: 1921
```

## 8. Stokes 网格求解

当前流场求解器是 lumen-only 的 Stokes / creeping-flow 近似求解器。

调用位置：

```text
utils/stokes_solver.py::solve_stokes_flow
```

它只在 `lumen_mask` 内部建立稀疏方程，不在血管外部背景区域建方程。

压力投影求解器配置为：

```yaml
pressure_solver: pcg
pressure_preconditioner: pyamg
```

最近一次验证中：

```text
pressure_projection_backend: scipy_lumen_only_pcg
pressure_preconditioner: pyamg_smoothed_aggregation
pressure_cg_info: 0
```

`pressure_cg_info: 0` 表示线性压力投影子问题正常收敛。

## 9. 当前不再限制速度

当前版本不再对速度做硬上限裁剪。

默认配置为：

```yaml
max_speed_um_s: 0.0
```

含义是：

```text
不限制速度
只记录最大速度作为诊断指标
```

最近一次验证中：

```text
max_speed_um_s: 1799.5701904296875
configured_max_speed_um_s: 0.0
speed_limit_enabled: false
```

也就是说，当前速度不是被 `max_speed_um_s` 裁剪出来的。

## 10. 边界通量软校正

仅靠 1D pressure 映射到网格边界时，实际入口/出口通量可能明显偏离血管模型中的目标流量。

因此当前版本加入了边界通量软校正：

```yaml
boundary_flux_relaxation: 0.5
```

它的作用是：

```text
根据当前入口/出口通量误差
沿血管方向轻微调整 root inlet 和 terminal outlet 的速度
让实际边界通量接近导出血管模型中的 flow_rate
```

这不是速度上限，也不是经验速度场。

它只作用于入口/出口边界通量误差。

最近一次验证中：

```text
relative_flux_error: 0.01184109530480767
boundary_flux_relative_error: 0.01184109530480767
max_allowed_relative_flux_error: 0.05
```

说明当前入口/出口通量误差低于 5% 阈值。

## 11. 内部散度检查

开边界附近本身存在流体进入和流体流出，因此不能把入口/出口附近的数值层当作普通内部区域来检查 `divergence = 0`。

当前版本的做法是：

```yaml
divergence_boundary_exclusion_cells: 8
```

含义是：

```text
计算内部 divergence_rms 时
排除 root inlet 和 terminal outlet 附近 8 个网格单元的开边界数值层
```

在默认 2 um 网格下：

```text
8 cells ~= 16 um
```

最近一次验证中：

```text
divergence_rms_s_inv: 0.9423559363555104
raw_divergence_rms_s_inv: 6.915585119896035
max_allowed_divergence_rms_s_inv: 1.0
```

解释：

| 指标                         | 含义                                          |
| ---------------------------- | --------------------------------------------- |
| `divergence_rms_s_inv`     | 排除开边界数值层后的内部散度，用于收敛判据    |
| `raw_divergence_rms_s_inv` | 不排除开边界数值层的全 lumen 散度，只作为诊断 |

因此，当前真正用于物理收敛的是 `divergence_rms_s_inv`，不是 `raw_divergence_rms_s_inv`。

## 12. 收敛判据

当前严格物理收敛需要同时满足：

```text
relative_change <= tolerance
divergence_rms_s_inv <= max_allowed_divergence_rms_s_inv
relative_flux_error <= max_allowed_relative_flux_error
wall_penetration_max_um_s <= max_allowed_wall_penetration_um_s
pressure_cg_info == 0
```

默认阈值：

```yaml
tolerance: 1.0e-3
max_allowed_divergence_rms_s_inv: 1.0
max_allowed_relative_flux_error: 0.05
max_allowed_wall_penetration_um_s: 1.0
```

最近一次验证中：

```text
physical_converged: true
relative_change: 0.0009962971110523552
divergence_rms_s_inv: 0.9423559363555104
relative_flux_error: 0.01184109530480767
wall_penetration_max_um_s: 4.76837158203125e-07
pressure_cg_info: 0
```

## 13. 最近一次成功运行

最近一次完整默认运行成功完成。

输出目录：

```text
ulm_microbubble_traj_gen/results/20260709_021305
```

关键结果：

| 指标                             |                      数值 |
| -------------------------------- | ------------------------: |
| `physical_converged`           |                  `true` |
| `iterations`                   |                   `217` |
| `relative_change`              | `0.0009962971110523552` |
| `divergence_rms_s_inv`         |    `0.9423559363555104` |
| `raw_divergence_rms_s_inv`     |     `6.915585119896035` |
| `relative_flux_error`          |   `0.01184109530480767` |
| `boundary_flux_relative_error` |   `0.01184109530480767` |
| `wall_penetration_max_um_s`    |  `4.76837158203125e-07` |
| `max_speed_um_s`               |    `1799.5701904296875` |
| `speed_limit_enabled`          |                 `false` |
| `pressure_cg_info`             |                     `0` |

输出文件：

| 文件                                   | 用途                                    |
| -------------------------------------- | --------------------------------------- |
| `velocity_and_wall_shear_field.npz`  | 保存速度场、壁面剪切力场、lumen mask 等 |
| `microbubble_field_trajectories.npz` | 保存微泡轨迹                            |
| `domain_metadata.yaml`               | 保存本次网格、求解器和收敛诊断信息      |
| `run_config.yaml`                    | 保存本次运行使用的配置                  |

## 14. 当前版本仍需注意的地方

1. 当前仍是 X-Z 平面上的网格流场，不是完整 3D 体网格流场。
2. `one_d_pressure_with_source_sink_flow_prior` 是 1D 网络和网格流场的混合边界方法，不是完整双向 1D-3D 强耦合。
3. `flow_prior_relative_error` 仍然较高，最近一次为 `0.9743624737813034`。当前真正用于通量收敛的是 `boundary_flux_relative_error`。
4. `raw_divergence_rms_s_inv` 较高主要来自开边界数值层，因此当前使用排除边界层后的内部散度作为收敛指标。
5. 如果更换血管模型、网格间距或流量尺度，需要重新检查 `divergence_boundary_exclusion_cells` 和 `boundary_flux_relaxation` 是否仍然合适。

## 15. 已运行检查

已运行编译检查：

```powershell
d:\anaconda3\envs\pmp\python.exe -m py_compile `
  ulm_microbubble_traj_gen\utils\stokes_solver.py `
  ulm_microbubble_traj_gen\utils\config.py `
  ulm_microbubble_traj_gen\utils\phiflow_solver.py
```

结果：

```text
通过
```

已运行默认完整生成：

```powershell
d:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\generate_microbubble_trajectories.py
```

结果：

```text
通过
完成 Stokes 流场求解
完成 Taichi 微泡轨迹推进
结果保存至 ulm_microbubble_traj_gen/results/20260709_021305
```
