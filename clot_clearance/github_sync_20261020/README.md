# Clot 项目同步快照 — 2026-10-20

目标仓库：`OZAKI39/ulm-3d-vascular-model-generation`  
新分支：`clot-clearance/project-sync-20261020`  
开发基线：`5036268`。本次汇总本地与服务器的项目代码、计算日志、研究数据、报告和可视化，未启动新仿真、修改材料或覆盖接受态。

## 当前方法与结果入口

- [项目开发入口](../ACTIVE_DEVELOPMENT_ENTRYPOINTS.md)
- [原 Y 血管血块上的 V4 方法说明](../analysis/y_clot_irregular_current_method_v4/REPRODUCE.md)
- [当前 Y 案例科学报告](../analysis/y_clot_irregular_current_method_v4/REVIEW_Y_CLOT_IRREGULAR_CURRENT_METHOD_V4.md)
- [最新固定视角碎裂视频与静态图](../results/y_clot_fracture_reference_style_v2/20261020T094240Z/OPEN_RESULTS.html)
- [最新可视化审查报告](../analysis/y_clot_fracture_reference_style_v2/REVIEW.md)
- [本次离线导航页](OPEN_RESULTS.html)
- [开发笔记原文快照](clot_rupture.md)

当前接受轨迹为 `y_clot_irregular_current_method_v4/20261020T041100Z`，计算至 3.6 ms。1927 个材料点、1561 个固定锚点、153986 条参考键；412 个内部材料点扰动，边界参考坐标保留。积分体积 0.240875 mm³。采用非均匀体积变分 LPS、双相不可逆渐进软化、当前表面载荷近似与壁接触。1000 Pa 为机械剪切代理的名义输入；没有新增声场、微泡声流或双向 FSI。材料与微结构仍未标定，不能作为已验证的临床血栓预测模型。

最新视频固定相机、Z 轴竖直，显示实际碎裂过程；149 个接受态、501 视频帧，物理时间 0–3.6 ms。视频未放大位移或插值节点。四个最新视频均包含在分支中，包括两种视野及 PPT 版本。

## 同步范围

保留既有 Git 历史及历次已审查快照，补齐 10 月 15 日以来逆向载荷、碎块机制、双相材料、表面耦合、V4 不规则点云、原 Y 血块计算及最新可视化的源代码、配置、C++ 扩展/补丁、数值验证、失败诊断、完整标量曲线、关键状态、报告和结果图。最新本地开发笔记按原字节复制，原路径链接作为历史记录保留。

服务器与本地 Python/C++/配置源码的共同文件已做 SHA-256 比较；唯一差异为 `pd_clot/portable_paths.py`，继续保留支持 GitHub 依赖和字体定位的本地版本。

另核对了 4837 个本地与服务器共同的待导出文件，其中 4834 个字节完全相同。三个差异是开发入口文档及旧 Y v1 的两张预览图；保留本地接受档案，并在 `server_variants/` 单独保存服务器版本及其 SHA，不覆盖历史结果。

遵循此前“不上传无关紧要的大文件”的要求：

- 大于 1 MiB 的文本数据可采用确定性 gzip，解压后的 SHA-256 与源文件完全相同；不做数值舍入或下采样。
- 密集状态采用 `compact_states/` 独立路径导出。保留坐标、损伤、原始 ID、时刻、速度、能量及已有外力等指定数组，每个 `.npy` 成员逐字节校验。未保留的数组逐文件列出。**这些是研究/回放数据，不是完整原生重启检查点。**
- 最新 Y 碎裂视频所用的全部 149 个接受态均有对应数据；原视频记录中的 153 个源文件身份均核对。原参考拓扑、非均匀积分体积和材料相区继续保留。
- 排除虚拟环境、共享库/对象、编译缓存、冗余重启和密集诊断流、过大的重复视频。单个新增文件小于 32 MiB。源数据仍保存在原工作区；未删除本地或服务器科学数据。

这是可核验的科研快照，**不是完整服务器镜像，也不是无需构建的独立求解器发行包**。原生重建仍需要对应上游依赖、锁定版本和构建记录。旧脚本中绝对路径与原库 SHA 是历史证据，没有批量改写。历史完整运行 manifest 描述原始全量档案；不要用它判定 Git 精选快照缺失。

## 校验与还原

从 `clot_clearance` 目录执行：

```bash
python github_sync_20261020/verify_export.py
# 需要 NumPy；只读取数据，不进行动力学计算：
python github_sync_20261020/verify_current_replay.py
# 仅还原被完整 gzip 压缩的文本；不同内容的已有文件拒绝覆盖：
python github_sync_20261020/verify_export.py --restore-text
```

`SNAPSHOT_MANIFEST.json` 覆盖本快照的项目文件，排除自身及校验结果以避免自引用。`EXPORT_RECORDS.local.json.gz` 与 `EXPORT_RECORDS.server.json.gz` 记录逐文件源 SHA、导出 SHA、压缩/字段选择和排除原因。`SOURCE_COVERAGE.json`、`SYNC_SUMMARY.json` 和 `VALIDATION.json` 提供最终统计。

`compact_states/<原相对路径>` 的保留数组可直接用 `numpy.load(..., allow_pickle=False)` 读取。字段裁剪后的 NPZ 不冒充原完整 NPZ，也不能直接覆盖被冻结的原路径。最新回放验证器给出如何按记录查找对应文件的示例。

一个继承的历史快捷链接 `analysis/clot_fragment_surface_coupling_v3/states` 仍指向未上传的完整原生状态目录，见 `INHERITED_LINK_LIMITATIONS.json`。该案例的紧凑状态位于本快照的 `compact_states/verification/clot_fragment_surface_coupling_v3/...`；未把紧凑状态冒充完整状态来替换旧链接。最新回放导航页的链接已单独核验。

服务器源根为 `/workspace/y_clot/project`；本地源根为 `/home/lzy/projects/temp_storage/clot_github_sync_20261015/clot_clearance`。完整源目录索引以 gzip 留存。此前服务器上已排除的历史上游/构建树和原生大数据沿用 [10 月 15 日快照说明](../github_sync_20261015/README.md) 及其排除清单，不重复导入。
