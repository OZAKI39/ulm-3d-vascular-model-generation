# 当前最佳可实现分流：三维流场可视化

[总预览 OPEN_RESULTS.html](OPEN_RESULTS.html) 已更新为 `mean-2p0-mmps-A-best-feasible-balance-v1` 最终第 71 步流场。

本轮新文件位于 [best_feasible_balance_v1/](best_feasible_balance_v1/)，数据说明、复现命令与验证结果见其 [README_ZH.md](best_feasible_balance_v1/README_ZH.md)。四段动画分别展示流线、局部速度矢量、压力和 WSS，均为 1080p、24 fps、18 秒；配套 4K 图片已生成。

沿用原放大显示、固定 Z 轴旋转、旁置标签和渐变文字。新压力最高约 5579 Pa，因此新结果压力色标扩展到 −5–6000 Pa，速度与 WSS 色标保持原范围。18 秒为稳态流场的展示旋转时长；显示流线条数不是流量占比。

原 H0 数据/结果仍位于本目录原有 `input_data/`、`results/`，可打开 [OPEN_H0_RESULTS.html](OPEN_H0_RESULTS.html)。本目录原始渲染代码保持不变，新结果使用子目录中的可复现副本。

## 保留的 H0 对照说明

打开 [OPEN_RESULTS.html](OPEN_RESULTS.html) 查看全血管流线、局部速度矢量、压力和 WSS。H0 对照输入来自 `mean-2p0-mmps-A-H0-pressure-v1` 的最终第 71 步冻结流场。审核说明见 [REVIEW_NEW_FLOW_ZH.md](REVIEW_NEW_FLOW_ZH.md)。

沿用原版放大显示、固定 Z 轴旋转、旁置箭头标签、透明文字及渐变切换。压力色标为用户已确认的 −5–5000 Pa，速度色标为 0–7.5 mm/s，WSS 为 0–55 Pa。

2026-09-27 网格输入目录更名为 `input_data/solver_mesh/`，读取脚本和当前校验清单已同步。视频、4K 图片、相机参数、标签布局和数据数值保持不变。`PATH_RENAME_PROVENANCE.json` 记录仅涉及路径的渲染代码变更，历史审核快照保留原文。

| 内容 | 路径 |
| --- | --- |
| 全血管流线动画 | [streamlines_full_vessel_enlarged.mp4](results/animations/streamlines_full_vessel_enlarged.mp4) |
| 局部速度矢量动画 | [velocity_vectors_branch_detail_enlarged.mp4](results/animations/velocity_vectors_branch_detail_enlarged.mp4) |
| 全血管压力动画 | [pressure_full_vessel_enlarged.mp4](results/surface_fields/animations/pressure_full_vessel_enlarged.mp4) |
| 全血管 WSS 动画 | [wss_full_vessel_enlarged.mp4](results/surface_fields/animations/wss_full_vessel_enlarged.mp4) |
| 流线、矢量 4K 图片 | [results/figures](results/figures/) |
| 压力、WSS 4K 图片 | [results/surface_fields/figures](results/surface_fields/figures/) |
| 新输入数据副本 | [input_data](input_data/) |
| 局部矢量数值及几何 | [results/data](results/data/) |

四段动画均为 1920×1080、24 fps、432 帧、18 秒；4K 概览图为 3840×2160。18 秒是相机旋转时长，流场保持稳态。

`input_data/frozen_flow/steady_flow_mean_2p0_mmps.vtu` 是兼容原接口的文件名，内容是新 `steady_flow_mean_2p0_mmps_A_H0.vtu` 的完整副本，SHA-256 为 `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。输入映射和逐文件校验值见 [INPUT_MANIFEST.json](INPUT_MANIFEST.json)。原始冻结清单保存在 `SOURCE_FROZEN_MANIFEST.json`，其中名称保持源文件名。

流线已按原算法从新速度场重新积分；压力表面和 WSS 已按原方法重新导出；68 个局部速度矢量取自新速度场。渲染阶段只读取这些已准备的数据，因此渲染报告中的 `new_streamline_integrations=0` 和 `source_fields_recomputed=false` 只描述渲染阶段，不代表沿用旧流场的派生数据。

`input_data/reports/render_manifest.json` 保留原局部观察中心，`input_data/field_diagnostics/MEDIA_VALIDATION.json` 只提供历史显示比例参考，不是本次结果的验证报告。本次报告位于 `results/`、`results/surface_fields/` 和 `audit/`。

在本目录复现渲染与校验：

```bash
./run.sh
./run.sh --surface-fields
./run.sh --validate
./run.sh --validate-surface-fields
sha256sum --check --quiet SHA256SUMS.txt
```

`run.sh` 自动寻找包含 `requirements.txt` 所列依赖的 Python，也可设置 `FEM_ROTATE_PYTHON`。重渲染使用本目录输入即可，无需重跑 FEM 或流线积分。重建派生输入的脚本为 `prepare_streamlines.py`、`prepare_surface_data.py`，其依赖的本机原算法路径记录于 `BUILD_CONTEXT.json` 和派生数据计算报告。

替换前的完整可视化包保存在 `../rotate_visualization_history/before_A_H0_data_20260924T190826Z/`。本轮没有修改 FEM 求解器、入口采样算法或 Particle 动力学，没有运行有限尺寸微泡生产模拟。
