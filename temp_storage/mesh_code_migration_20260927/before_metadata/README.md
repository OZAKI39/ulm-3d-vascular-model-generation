# 当前网格工具、共享环境和流场可视化

- [旋转可视化](rotate_visualization/OPEN_RESULTS.html)：Network-H0 新边界条件流线、速度矢量、压力和 WSS，保留原有渲染风格。
- `rotate_visualization/input_data/`：当前可视化输入；`results/`：视频、4K 图片和矢量数据。
- `scripts/`、`src/`、`inputs/`、`configs/`、`outputs/sv1/`：仍需使用的网格准备工具与输入。
- `external/SimVascularDistribution/`：官方内嵌 Python / TetGen 网格工具。
- `.venv/`：被当前共享 Python 环境引用的依赖，不能作为旧阶段缓存删除。

当前 FEM 算例位于 [ulm_flow_mean_2p0_mmps](/home/lzy/projects/ulm_flow_mean_2p0_mmps/README.md)。旧求解阶段、构建环境、性能试验、日志和旧渲染副本已删除。
