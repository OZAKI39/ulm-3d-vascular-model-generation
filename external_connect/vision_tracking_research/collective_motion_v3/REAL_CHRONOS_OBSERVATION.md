# 真实Chronos录像观测

EXPLORATORY · PHYSICS_UNCALIBRATED · training_ready=false

原始937帧都已完成PIV，详细见FULL_VIDEO_PIV_REPORT。整段因果回放均以源PTS定速，实际处理/丢弃可追溯。未连接真实相机，LIVE_CAMERA_PASS=false。

固定预登记selected clip为360–423帧，CoTracker只看当前及之前16帧，不使用未来；PIV/OT/LRAFT使用过去第2帧至当前。64个展示帧中各方法质量/风险门槛接受数：{'PIV': 64, 'LRAFT_ST': 19, 'CoTracker3': 64, 'POT_unbalanced': 4}。CoTracker局部轨迹与可见泡群运动相符是观察性描述，不是已核实单泡材料身份。LRAFT或OT拒判时视频显示null，不把其低可信原始预测绘成可信速度。

| Method | Processed /937 | Dropped | Full ready p50 / p95 ms | Fresh valid states | Fresh six Hz |
|---|---:|---:|---:|---:|---:|
| PIV | 860 | 77 | 36.26 / 45.88 | 667 | 21.400 |
| LRAFT_ST | 771 | 166 | 46.35 / 55.63 | 173 | 5.550 |
| CoTracker3 | 881 | 56 | 79.66 / 85.56 | 50 | 1.604 |
| POT_unbalanced | 922 | 15 | 7.34 / 12.14 | 66 | 2.118 |

REAL_SIX_PARAMETER_ACCURACY=NOT_VERIFIED，跨实验泛化未验证。候选中心、密度和2D散度均不代表3D物理团收缩率。
