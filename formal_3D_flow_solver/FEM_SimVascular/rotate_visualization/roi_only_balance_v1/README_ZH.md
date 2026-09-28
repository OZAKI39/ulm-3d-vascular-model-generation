# ROI-only 新边界流场旋转可视化

打开 [总预览 OPEN_RESULTS.html](OPEN_RESULTS.html)。数据来自 `mean-2p0-mmps-A-ROI-only-balanced-pressure-v1` 的最终第 71 步；冻结场 SHA256 `fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12`。

生成内容为全血管流线、局部速度矢量、全血管压力和 WSS：1920×1080、24 fps、每段 432 帧/18 秒，配套 3840×2160 静态图。18 秒是稳态场的展示旋转时长，不是瞬态 CFD 或微泡运动时间。

| 内容 | 路径 |
| --- | --- |
| 全血管流线 | [MP4](results/animations/streamlines_full_vessel_enlarged.mp4) |
| 局部速度矢量 | [MP4](results/animations/velocity_vectors_branch_detail_enlarged.mp4) |
| 全血管压力 | [MP4](results/surface_fields/animations/pressure_full_vessel_enlarged.mp4) |
| 全血管 WSS | [MP4](results/surface_fields/animations/wss_full_vessel_enlarged.mp4) |
| 流线/矢量 4K 图片 | [figures/](results/figures/) |
| 压力/WSS 4K 图片 | [surface_fields/figures/](results/surface_fields/figures/) |

## 数据和显示方式

沿用原放大视图、固定 Z 轴旋转、6° 仰角、黑色背景、英文图内标签、旁置箭头及透明文字渐变。主渲染器与父目录代码字节一致；表面渲染器的全部函数/类相同，仅压力色标数值更新到 **−5–7000 Pa**，覆盖实际最大值 6445.547523 Pa。速度色标仍为 0–7.5 mm/s，WSS 为 0–55 Pa。`audit/STYLE_PRESERVATION.json` 保存核查。

当前压力范围 -1.017965–6445.547523 Pa。原始壁面面片 WSS 0.920194–46.139020 Pa，面积均值 15.429464 Pa。WSS 根据最终 FEM 完整 P1 速度梯度计算壁面切向黏性牵引，mu=0.00345312 Pa·s，排除入口和出口截面。动画沿用面积加权节点显示值及渲染插值；原始面片值、面积、法向和牵引向量均保留在 `input_data/field_diagnostics/data/wall_wss_si.vtp`。没有裁剪物理数据或新增平滑，也没有由本次绘图证明网格无关。

从新速度场重新积分 768 个可视化候选，实际终止事件 `{'OUTLET_02': 192, 'OUTLET_01': 257, 'OUTLET_03': 290, 'POINT_TIME_HORIZON_30S': 29}`。选出的 96 条均自然抵达真实出口并通过步长减半检验；最大路径差 0.001719647 μm。保留原显示配额 O1/O2/O3=16/56/24，**显示条数不代表流量分配**。实际 CFD 分流为 34.50368463% / 27.27289025% / 38.22342512%。时间上限候选未作为已完成路径。

原始场及网格位置：`/home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1`。兼容旧接口的 `steady_flow_mean_2p0_mmps.vtu` 与 `physics_validation_H0.json` 实际是本次 ROI-only 新场及验收报告副本，不能凭文件名判断为 H0。`render_manifest.json` 仅继承观察中心；输入目录旧 `MEDIA_VALIDATION.json` 仅提供历史显示比例参考。本次验证在 `results/` 和 `results/surface_fields/`，最终汇总为 `VISUALIZATION_VALIDATION.json`。

## 复现

服务器渲染目录 `/workspace/roi_only_flow_visualization_20260928`，使用现有 NVIDIA EGL/RTX 4090 绘制，CPU libx264 编码。实际 OpenGL 设备写入各 MEDIA_VALIDATION.json。没有启动 CFD、有限尺寸微泡或 RBC 计算；只运行新场派生数据处理及绘图。

```bash
cd /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1
./run.sh
./run.sh --surface-fields
./run.sh --validate
./run.sh --validate-surface-fields
sha256sum --check --quiet SHA256SUMS.txt
```

上述命令使用已准备的数据。重新生成流线时应使用独立输出副本；`prepare_streamlines.py` 会拒绝覆盖已完成流线。表面数据可由父目录原 `prepare_surface_data.py --case <source case> --output <new output>` 重建；源脚本及 WSS 核心哈希保存在计算记录。原 CFD 输入/结果逐文件保护清单位于 `audit/SOURCE_CASE_LOCK.json`。

旧 H0 与 full-A best-feasible 可视化继续保留在父目录和 `best_feasible_balance_v1/`，不会被新数据覆盖。
