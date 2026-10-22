# 当前权威方法与边界

自有运行实现统一在 `src/`。不同工况只保留各自配置和数据；`cases/*/reports/executed_*.py`、诊断目录内原脚本和服务器记录是历史执行证据，不是另一套日常源码。

| 功能 | 当前源码或身份位置 | 使用边界 |
|---|---|---|
| BraVa 几何与端口 | `src/human_vessel_flow/geometry.py`、`transforms.py` | 已封闭四端口几何，mm→m；0/15 为姿态变换 |
| 网格 | `meshing.py`、`mesh_audit.py`、`native/tetgen_mesh.py` | SimVascular/TetGen；保留原权威网格，本轮未生成新网格 |
| 新边界设计 | `boundary.py`、`fem.py` | 中心线阻力估算出口压力，未作新三维分流验证；不改旧校准压力 |
| FEM 支持 | `native_runner.py`、`src/flow_solver_support/` | 原日志、稳态和 checkpoint 门槛；生产 P1/P1 与 VMS |
| 原生生产求解器 | `native/svmultiphysics/`、`native/production_patch/`、`scripts/build_native.sh` | GPU 同步生产补丁；独立安装在 `new_vessel_flow_solvers/brava_svmultiphysics`，P2 实验补丁未混入 |
| 冻结场读取 | `flow_field.py`、`src/particle_3d/field.py` | 两个场各自 SHA，原节点/单元身份；SI 原始场不按打印姿态复制为新物理解 |
| WSS / 流量 | `wss.py`、`postprocess.py`、`src/flow_solver_support/wss.py` | 原网格梯度与原面片切向牵引；节点显示值不能替代原值 |
| 微泡动力学 | `particles.py`、`src/particle_3d/` 的真实依赖闭包 | 有限尺寸、壁面/泡间接触、阻力、CPU 回退与 FP64 GPU；不替换为点流线 |
| 固定时间窗 | `src/particle_3d/population_v1/run_fixed_window_population.py`、`configs/population/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json` | 灌注 1.240815660 s，观察 10 s，dt=0.5 ms，4e9 m⁻³；没有新正式群体 |
| V2 / P1 诊断 | `src/human_vessel_flow/reference_pipe/`、`diagnostics/reference_pipe/` | 保留场误差分解、exact/P1 敏感性及冻结 V2 生成来源；科学 FAIL 不变 |
| P2/P1 研究 | `diagnostics/p2p1/`、`p2_geometry.py`、`p2_gpu_budget.py` | 原曲面、边界、输出及资源证据；STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT，不是正式生产路线 |
| 可视化 | `visualization.py`、`results/visualization/` | 保留原四视角黑底蓝泡动画；逐帧低速筛选只影响显示，不改账本 |

表内省略包前缀的文件位于 `src/human_vessel_flow/`。实际 Python 模块和符号来源仍见 `MODULE_MAPPING.csv`、`FUNCTION_CLASS_MAPPING.csv`。

历史正式场 SHA 为 `73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`，历史强制等流量校准场 SHA 为 `f5d7a2ccf7cdcc08cbfebd725d0a17e23e5780df23e4ef6488a26ffa5983d4d2`。两者不能混称。

旧群体仍是 1565 接受、458/527/437 出流、141 观察截止、2 数值失败。REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC、MIXED_FIELD_ERROR 和 P1_FIELD_TRAJECTORY_IMPACT_LARGE 保留，V3 仅提案。全局守恒、CPU/GPU 一致或展示流畅都不表示这些科学问题已经解决。
