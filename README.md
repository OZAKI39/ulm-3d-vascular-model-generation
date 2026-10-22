# 微泡聚集与声场仿真 — F10AG 同步快照

本分支只包含数值仿真项目：PZT 声场、微泡声响应、逐泡运动、接触求解、控制探索、数值验证与科研结果。实验相机、信号发生器和实时 GUI 位于独立分支 [sync/experiment-control-vision-20261022](https://github.com/OZAKI39/ulm-3d-vascular-model-generation/tree/sync/experiment-control-vision-20261022)。

源工作区：`/home/lzy/projects/computation_examples/_archives/mevo_capture_f10ag_20261019`；原始 HEAD：`4e80e6b18e6de014acb9bbb611ed6c4c3bcb61a2`。这是独立快照分支，不导入庞大的开发 Git 历史；原历史分支不变。继承的源码、报告和原执行身份保持原样，本次没有重新计算或改写科研结论。

## 阅读入口

- [项目代码、配置与依赖](computation_examples/mevo_acoustic_microbubbles/)
- [F10AG 完整报告](computation_examples/mevo_acoustic_microbubbles/results/straight_tube_exp1/f10ag/F10AG_DEVELOPMENT_REPORT.md)
- [F10AG 图表](computation_examples/mevo_acoustic_microbubbles/results/straight_tube_exp1/f10ag/figures/)
- [F10AG 固定 ID 数值数据](computation_examples/mevo_acoustic_microbubbles/results/straight_tube_exp1/f10ag/dataset/)
- [共享项目流程说明（历史原文）](docs/shared_project_overview_20261022.md)
- [同步范围](sync/20261022/EXPORT_POLICY.json)、[来源身份](sync/20261022/SOURCE_IDENTITIES.json)、[实际文件清单](sync/20261022/COPIED_FILES.json)、[排除清单](sync/20261022/EXCLUDED_FILES.json)

## 数据与恢复

保留历史代码、关键几何/响应输入、F10AG/F10AF 数值数据、验收与失败证据、科研图、运行日志。补入服务器上的 APECSS 响应、初态及 2.4/2.5/2.6 MHz 粗网格缓存。

较大的文本和 HDF5 文件采用 **gzip 无损封装**；未抽样、未修改数值。原路径和原 SHA-256 保存在清单中。运行前按需还原：

```bash
python sync/20261022/restore_payload.py             # 校验压缩文件及解压后原始哈希
python sync/20261022/restore_payload.py --restore   # 在当前克隆中还原原文件名
```

环境没有打包，按项目 `pyproject.toml` 和各阶段依赖锁定文件重建；历史脚本可能包含原 WSL/服务器路径，需要按本机布局调整。不要用原历史输出目录重跑。同步检查记录在 [VALIDATION.json](sync/20261022/VALIDATION.json)。

完整逐步 HF 原始数据、可再生成的细网格声场缓存、大型 LF 诊断轨迹、重复绘图数组及视频未整批进入 Git。服务器恢复位置和哈希见 [REMOTE_RAW_INDEX.json](computation_examples/mevo_acoustic_microbubbles/results/straight_tube_exp1/f10ag/REMOTE_RAW_INDEX.json)；本分支不是这些外部大文件的完整灾备。历史清单可能引用未打包文件，以 `sync/20261022` 为本次同步依据。

科学适用范围保持：`EXPLORATORY`、`HOMOGENEOUS_WATER_CONTROL_PROXY`、`PHYSICS_UNCALIBRATED`、`training_ready=false`。局部捕获、条带与持续旋转仍须分别验收，不能因代码同步而扩大原结论。
