# generate_microbubble_trajectories.py 当前输出检查说明

当前版本只生成正常微泡流动轨迹，不再生成病灶黏附、状态分类、事件码或兼容字段。

## 数据来源

| 数据 | 来源 |
|---|---|
| 血管几何 | `ulm_vascular_model_generator/vessel_swc_models/<model>.swc` |
| 血管流量、压力、黏度等字段 | `ulm_vascular_model_generator/vessel_swc_models/<model>.vessels.npz` |
| 微泡仿真参数 | `ulm_microbubble_traj_gen/configs/bound_plane_default.yaml` |

## bubbles_*.mat 核心字段

| 字段 | 含义 |
|---|---|
| `bubbles` | 每个微泡的轨迹结构数组。 |
| `stats` | 轨迹库半径和直径分布统计。 |
| `transport_network` | 每段血管的半径、长度、流量、压力梯度、速度和壁面剪切力。 |
| `vascular_source` | 本次读取的 SWC 和 `.vessels.npz` 路径。 |

## 单个 bubble 核心字段

| 字段 | 含义 |
|---|---|
| `XYZ_laminar` | 微泡在每个采样点的三维坐标，单位 um。 |
| `vf_array` | 每个采样点的轴向速度，单位 um/s。 |
| `ecg_array` | 每个采样点使用的脉动流倍率。 |
| `closest_nodes` | 每个采样点最近的血管图节点。 |
| `radii` | 考虑微泡半径后，微泡中心可用的局部半径。 |
| `transport_vessel_id` | 每个采样点所在血管段 ID。 |
| `transport_flow_rate_um3_s` | 每个采样点所在血管段的流量。 |
| `transport_wall_shear_stress_pa` | 每个采样点所在血管段的壁面剪切力。 |
| `microbubble` | 该微泡的半径、注入时间、初始截面位置等个体属性。 |

## frames_*.mat 核心字段

| 字段 | 含义 |
|---|---|
| `frames` | 每帧微泡矩阵：Bubble ID、轨迹采样点、X、Y、Z。 |
| `frames_velocities` | 每帧每个微泡的速度。 |
| `frames_radii` | 每帧每个微泡的局部可用半径。 |
| `frames_poiseuille` | 每帧每个微泡的 Poiseuille 速度因子。 |
| `frames_ecg` | 每帧每个微泡的脉动流倍率。 |
| `frames_vessel_id` | 每帧每个微泡所在血管段 ID。 |
| `frames_flow_rate_um3_s` | 每帧每个微泡所在血管段流量。 |
| `frames_wall_shear_stress_pa` | 每帧每个微泡所在血管段壁面剪切力。 |

## quick test 检查项

| 检查项 | 目的 |
|---|---|
| `frames` 尺寸 | 确认帧矩阵行数和列数符合配置。 |
| `frames_velocities` 尺寸 | 确认速度矩阵与帧数、每帧微泡数一致。 |
| 速度有限且非负 | 避免 NaN、Inf 或负速度进入下游。 |
| 单条 bubble 序列长度一致 | 确认坐标、速度、半径和节点序列可一一对应。 |
