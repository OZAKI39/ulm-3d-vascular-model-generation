# 血栓模拟开发归档 / Clot-clearance development snapshot

本目录同步自现有 `/home/lzy/projects/clot_clearance`，包含直管初始原型、旧强制碎裂演示，以及最新能量正则化开发。原始源码、配置、日志、科学数据与历史报告按字节保留；本地算例未修改。仓库中其余血管建模文件未改变。

## 最新结果入口

- [中文报告（GitHub 相对链接版）](github_sync/reports/REGULARIZATION_REPORT_ZH.md)
- [English report (relative links)](github_sync/reports/REGULARIZATION_REPORT_EN.md)
- [最新粒子损伤动画 MP4](visualization/streaming_regularized_demo/regularized_damage_transport.mp4)
- [GIF](visualization/streaming_regularized_demo/regularized_damage_transport.gif)、[预览图](visualization/streaming_regularized_demo/preview.png)
- [13 组 PNG/PDF 科学图](visualization/streaming_regularized_demo/figures/)
- [本地 HTML 查看器](visualization/streaming_regularized_demo/index.html)：下载后用浏览器打开，或从本目录运行 `python -m http.server 8000 --bind 127.0.0.1`。
- [18 组算例汇总 CSV](verification/regularization/comparisons/study_summary.csv)
- [当前科学状态](verification/regularization/FINAL_STATE.json)、[原测试记录](verification/regularization/FINAL_TESTS_SCOPED.log)
- [本次复制后测试记录](github_sync/COPY_TESTS.log)、[副本校验](github_sync/SNAPSHOT_CHECK.json)

当前结论：48 项原测试通过，18 组计算完成；Gc_demo=0.01 J/m² 的指定 coupon 标定通过。**主算例没有断键或脱落，模型 NOT CONVERGED，碎片分辨率 NOT_ASSESSED_NO_DETACHMENT。** 不将制造流场称为真实微泡 CFD，不将旧强制碎裂结果当作新模型结果。

动画保留黑底、Arial 常规字重、无箭头及真实位移；源位于 clot 左面上游 0.2 mm。实际损伤色标为 0–0.0005，黄色不表示 D=1。原科学报告和哈希不因 GitHub 展示而改写；相对链接报告仅是附加导航副本。

## 同步范围与排除

完整文件清单、源 SHA-256、分片 SHA-256 和逐文件排除原因见 [SNAPSHOT_MANIFEST.json](github_sync/SNAPSHOT_MANIFEST.json)。

保留全部求解器/分析/可视化脚本、测试、配置、输入场、源码快照、执行日志、失败诊断、标定曲线、能量/拓扑/清除账本、全部算例的科学 NPZ 状态，以及正式图表和 MP4。旧 results 兼容链接以相对符号链接保留。

排除约 672 MB：构建/缓存/临时安装、可由保留的算例重导出的逐帧 VTK 及其 PVD 索引、布局草稿、可重新插值得到的 display_samples、两个大型旧版重复 GIF、第三方论文全文及其全文提取副本。输入 VTU/VTP 与导入器验证夹具保留。论文 DOI 与 [实现边界说明](references/PAPER_READING.md) 保留；论文为 [Karmakar et al. (2026)](https://doi.org/10.1016/j.apples.2026.100347)。

保留的原始文件载荷约 477 MB（十进制，不含 Git 压缩）。没有科学数据降采样或数值改写。历史保护清单仍记录完整原工作区；其中被排除的衍生文件不属于 GitHub 包完整性失败，应使用本次 snapshot manifest 校验。

## 还原大文件与校验

细网格 `results/streaming_regularized_demo/fine/states.npz` 原始大小 109,514,198 字节，按最多 40 MiB 的原始字节块拆为三片，避免 GitHub 单文件限制。无需 Git LFS、外部存储或额外 Python 库。

从仓库根目录：

```bash
cd clot_clearance
python github_sync/verify_snapshot.py
python github_sync/restore_large_files.py
```

还原器逐片校验、检查拼接后原始 SHA-256，拒绝覆盖不同数据。还原出来的 `fine/states.npz` 已被局部 `.gitignore` 排除，避免再次提交超限文件。`--verify-only` 可只核对分片。已还原且校验一致的文件可重复使用。保留的 18 组正式算例数据可在还原后读取；批量 VTK 导出需在新输出目录重新运行对应算例。

## 运行环境与测试

本子项目使用与仓库原血管几何工具独立的环境。记录的环境是 Linux/WSL、Python 3.13.11，主要库见 [provenance/ENVIRONMENT.json](provenance/ENVIRONMENT.json)，本次环境版本见 [ENVIRONMENT.json](github_sync/ENVIRONMENT.json)。建议为该子项目单独建立环境；依赖列表是 [requirements-recorded.txt](github_sync/requirements-recorded.txt)。这些是实际使用的版本记录，不表示已在所有平台验证。

在具备所列依赖的 Python 环境中，从 `clot_clearance/` 执行：

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PWD"
export NUMBA_CACHE_DIR="$PWD/build_github/numba_cache"
export MPLCONFIGDIR="$PWD/build_github/matplotlib"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -B -m pytest tests -q -p no:cacheprovider --basetemp build_github/pytest_tmp
```

必须指定当前 `tests` 目录，避免收集 `provenance` / `source_snapshot` 内冻结的同名测试副本。

最小新算例与能量标定：

```bash
python -B -m pd_clot.runner --config configs/minimal.json --output runs/minimal_new
python -B -m pd_clot.regularization.coupon \
  --config configs/streaming_regularized_demo_coarse.json \
  --spacing 0.000125 --Gc-demo 0.01 \
  --output verification/regularization/calibration_replay
python -B -m pd_clot.regularization.runner \
  --config configs/streaming_regularized_demo_coarse.json \
  --calibration verification/regularization/calibration_replay/CALIBRATION.json \
  --output results/streaming_regularized_demo_replay/coarse
python -B -m pd_clot.regularization.campaign --study mesh --output-root results/streaming_regularized_demo_replay
python -B -m pd_clot.regularization.campaign --study cycle --output-root results/streaming_regularized_demo_replay
```

使用新的输出目录，避免覆盖历史结果。初始直管和强制碎裂流程见 [README.md](README.md)、[README_FRAGMENTATION.md](README_FRAGMENTATION.md)。

## 可移植性边界

科学计算和测试使用项目相对路径；历史日志、审计清单、报告生成器及旧保护检查保留原工作区绝对路径，用于追溯，不能直接当作新电脑上的全工作区检查。原可视化脚本使用本地 Arial 文件 `/mnt/c/Windows/Fonts/arial.ttf`；未上传商业字体。重新渲染时，支持 `--font` 的渲染器可指定自己的合法 Arial 文件；其余绘图入口需要在副本中适配字体路径。现有 MP4/PNG/PDF 可直接查看。

下一个科学步骤是导入实际微泡解析 VTU（u、p、明确单位、压力基准、坐标配准和足够的表面覆盖），执行导入诊断与一宏步对照。该快照没有新增或恢复任何 BraVa 批次或服务器计算。
