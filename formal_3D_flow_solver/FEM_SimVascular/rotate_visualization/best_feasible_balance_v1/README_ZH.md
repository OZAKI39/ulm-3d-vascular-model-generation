# 最佳可实现分流：新三维流场旋转可视化

源算例：`mean-2p0-mmps-A-best-feasible-balance-v1`，唯一 CFD 的最终第 71 步。流场 SHA256：`fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2`。

打开 [OPEN_RESULTS.html](OPEN_RESULTS.html) 查看四段动画；18 秒是稳态场的展示旋转时长，不是非稳态血流模拟时间，也不是微泡运动时间。

| 内容 | 文件 |
| --- | --- |
| 全血管流线 | [streamlines_full_vessel_enlarged.mp4](results/animations/streamlines_full_vessel_enlarged.mp4) |
| 局部速度矢量 | [velocity_vectors_branch_detail_enlarged.mp4](results/animations/velocity_vectors_branch_detail_enlarged.mp4) |
| 全血管压力 | [pressure_full_vessel_enlarged.mp4](results/surface_fields/animations/pressure_full_vessel_enlarged.mp4) |
| 全血管 WSS | [wss_full_vessel_enlarged.mp4](results/surface_fields/animations/wss_full_vessel_enlarged.mp4) |

四段均为 1920×1080、24 fps、432 帧、18 秒；静态图片为 3840×2160，分别位于 `results/figures/` 和 `results/surface_fields/figures/`。

沿用已有放大视图、固定 Z 轴旋转、6° 仰角、旁置箭头标签、透明文字和渐变换位。共享流线/矢量渲染器与原文件字节一致；表面渲染器全部函数和类保持不变，只将压力色标从 −5–5000 Pa 扩展到 **−5–6000 Pa**。速度为 0–7.5 mm/s，WSS 为 0–55 Pa。样式核查见 `audit/STYLE_PRESERVATION.json`。

实际压力范围为 -1.970180–5579.380780 Pa；原始壁面面片 WSS 为 0.811451–46.261386 Pa，面积均值 14.113265 Pa。没有裁剪或归一化物理数据。WSS 动画继续使用原界面的面积加权节点显示值（并由渲染器插值显示）；原始面片 WSS、牵引向量和面积均保存在 `input_data/field_diagnostics/data/wall_wss_si.vtp`。

流线由新 P1 速度场通过原 double-precision RK23 点流线算法重新积分。原 768 个可视化候选中，自然终止事件为 `{'OUTLET_02': 245, 'OUTLET_01': 105, 'OUTLET_03': 387, 'POINT_TIME_HORIZON_30S': 31}`；仅选取自然抵达出口且通过几何步长加密检查的 96 条，最大路径加密差 0.00308900 μm。保留原显示配额 O1/O2/O3=16/56/24，**显示流线数量不代表出口流量比例**。真实 3D 分流为 O1=15.51494709%、O2=32.91168032%、O3=51.57337259%。时间上限候选未被伪装为完成路径。未运行 CFD、有限尺寸微泡或 RBC。

兼容接口文件名 `steady_flow_mean_2p0_mmps.vtu` 和 `physics_validation_H0.json` 分别存放新流场和新独立验收报告的完整副本，不能按旧文件名理解为 H0 数据。当前来源与别名关系见 `BUILD_CONTEXT.json`、`INPUT_MANIFEST.json`。`render_manifest.json` 只继承局部观察位置，旧 `MEDIA_VALIDATION.json` 只提供显示放大比例参考。

## 复现与检查

```bash
cd /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/best_feasible_balance_v1
./run.sh
./run.sh --surface-fields
./run.sh --validate
./run.sh --validate-surface-fields
sha256sum --check --quiet SHA256SUMS.txt
```

`run.sh` 自动使用现有共享渲染环境 `/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`。未安装新软件。FEM 环境缺少 imageio_ffmpeg 的首次导入失败已保留在 `logs/render_environment_initial_failure.log`，随后改用已有完整渲染环境。

重新生成流线需在独立新输出副本中运行本目录 `prepare_streamlines.py`，该脚本拒绝覆盖已完成流线。表面数据生成入口是父目录 `prepare_surface_data.py --case <上述源算例> --output <独立输出目录>`，使用现有 `solver_support` WSS 核心。

原 H0 可视化数据和动画保留在父目录原有 `input_data/`、`results/` 中；本轮所有新数据和图像位于本目录。完整逐帧校验、标签检查和输入哈希分别保存在新 `results/` 及 `results/surface_fields/` 下，最终汇总为 `VISUALIZATION_VALIDATION.json`。
