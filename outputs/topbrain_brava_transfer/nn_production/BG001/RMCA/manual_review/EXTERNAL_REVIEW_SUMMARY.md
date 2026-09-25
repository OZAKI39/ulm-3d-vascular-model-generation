# BG001 RMCA 候选边界审查摘要

## A. 来源与范围

审查对象为已冻结的 BG001 RMCA、TYPE 4。映射状态 DERIVED_SOURCE_SUPPORTED，仅适用于当前指定文件。raw SHA256：`5b69c6c406e4724ffc15a2a4c6b5958cbb169b822a05098d6ebedeecbbb3afbe`；ColorCoded SHA256：`733b255aeeb7795de45799f64e73b19a93ada45ca63518e9aa479c203658ef0e`。分支 CSV SHA256：`9f67c1accba48469610f93d27a028ad140c76e0f1c2f60134eefe8683c373758`。完整 ROI manifest 源绑定及下载版本说明见包内 evidence。

本次导出只读取生产结果。registration、ranking、NN、UNKNOWN 校准、分支聚合、ROI 提取、VascularMD 建模及表面生成均未执行；真实导出时相应入口被禁止调用断言保护。所有审查状态仍为 UNREVIEWED。

## B. RMCA 总体统计

共 797 个原始节点、55 条原生拓扑分支：M1 8 条，MeVO 44 条，UNKNOWN 3 条。严格 MeVO 为 5 个组件，合计 687 个原始节点、44 条分支。完整分支保留共享端点，因此与点级 MeVO 的 678 个节点数量不同。供体为 MRA024、013、015、022、011；既有支持阈值为 agreement≥0.60、有效 donor≥3、归一化距离≤12.7949772843。

## C. 候选事件计数

共 **44 个事件**，对应 **31 个原始 SWC 位置**。同一组件分叉根可有多条首层分支，也可同时产生不同类型事件。

| 类型 | 数量 |
|---|---:|
| PROXIMAL_M1_TO_MEVO | 6 |
| DISTAL_MEVO_TO_UNKNOWN | 0 |
| DISTAL_MEVO_TO_TERMINAL | 26 |
| UNKNOWN_TO_MEVO | 2 |
| MEVO_TO_M1_REVERSE | 0 |
| MEVO_COMPONENT_ROOT | 8 |
| DETACHED_MEVO_COMPONENT | 2 |

## D. 五个组件与分离原因

| 组件 | 节点/分支 | 首层 MeVO 分支 | 直接上游标签 | 与其他组件的分离原因集合 |
|---|---:|---|---|---|
| roi_part01 | 142 / 11 | 95_96 | M1 | DIFFERENT_PROXIMAL_PARENT;UNKNOWN_GAP |
| roi_part02 | 44 / 1 | 107_111 | M1 | DIFFERENT_PROXIMAL_PARENT;UNKNOWN_GAP |
| roi_part03 | 330 / 24 | 112_113;112_128 | UNKNOWN | UNKNOWN_GAP |
| roi_part04 | 130 / 6 | 137_138;137_143 | M1 | DIFFERENT_PROXIMAL_PARENT;UNKNOWN_GAP |
| roi_part05 | 41 / 2 | 144_145;144_146 | M1 | DIFFERENT_PROXIMAL_PARENT;UNKNOWN_GAP |

分离原因根据既有 branch-parent 路径描述，未修改任何连接。完整两两关系、连接路径和各路径标签见 `tables/component_separation.json`。UNKNOWN_GAP 指相关组件之间的原生路径含 UNKNOWN；DIFFERENT_PROXIMAL_PARENT 指它们具有不同的直接近端父分支；其余原因的定义见该 JSON。

## E. 全部边界简表

短 ID 链接到独立四面板图；完整 ID、原始节点号和全部数值在边界 CSV 中。

| 短 ID | 类型 | part | 上游 → 下游 |
|---|---|---:|---|
| [B0276d7d6a2](boundaries/BG001_RMCA_B0276d7d6a2/BG001_RMCA_B0276d7d6a2_review.png) | PROXIMAL_M1_TO_MEVO | 01 | 94_95 → 95_96 |
| [Be9d570d1ed](boundaries/BG001_RMCA_Be9d570d1ed/BG001_RMCA_Be9d570d1ed_review.png) | MEVO_COMPONENT_ROOT | 01 | 94_95 → 95_96 |
| [B766d3abd5b](boundaries/BG001_RMCA_B766d3abd5b/BG001_RMCA_B766d3abd5b_review.png) | PROXIMAL_M1_TO_MEVO | 02 | 95_107 → 107_111 |
| [Bb8f7da57b7](boundaries/BG001_RMCA_Bb8f7da57b7/BG001_RMCA_Bb8f7da57b7_review.png) | MEVO_COMPONENT_ROOT | 02 | 95_107 → 107_111 |
| [B59f02eeabe](boundaries/BG001_RMCA_B59f02eeabe/BG001_RMCA_B59f02eeabe_review.png) | DETACHED_MEVO_COMPONENT | 03 | 94_112 → 112_113 |
| [B6bd8df27e7](boundaries/BG001_RMCA_B6bd8df27e7/BG001_RMCA_B6bd8df27e7_review.png) | DETACHED_MEVO_COMPONENT | 03 | 94_112 → 112_128 |
| [B84196289a3](boundaries/BG001_RMCA_B84196289a3/BG001_RMCA_B84196289a3_review.png) | MEVO_COMPONENT_ROOT | 03 | 94_112 → 112_113 |
| [Bc8983008b2](boundaries/BG001_RMCA_Bc8983008b2/BG001_RMCA_Bc8983008b2_review.png) | UNKNOWN_TO_MEVO | 03 | 94_112 → 112_113 |
| [Bf1145d6339](boundaries/BG001_RMCA_Bf1145d6339/BG001_RMCA_Bf1145d6339_review.png) | MEVO_COMPONENT_ROOT | 03 | 94_112 → 112_128 |
| [Bf349521b9e](boundaries/BG001_RMCA_Bf349521b9e/BG001_RMCA_Bf349521b9e_review.png) | UNKNOWN_TO_MEVO | 03 | 94_112 → 112_128 |
| [B03ae50e0e3](boundaries/BG001_RMCA_B03ae50e0e3/BG001_RMCA_B03ae50e0e3_review.png) | PROXIMAL_M1_TO_MEVO | 04 | 93_137 → 137_138 |
| [B328beeea6a](boundaries/BG001_RMCA_B328beeea6a/BG001_RMCA_B328beeea6a_review.png) | PROXIMAL_M1_TO_MEVO | 04 | 93_137 → 137_143 |
| [B7e9b954694](boundaries/BG001_RMCA_B7e9b954694/BG001_RMCA_B7e9b954694_review.png) | MEVO_COMPONENT_ROOT | 04 | 93_137 → 137_143 |
| [Bef0c71cf9b](boundaries/BG001_RMCA_Bef0c71cf9b/BG001_RMCA_Bef0c71cf9b_review.png) | MEVO_COMPONENT_ROOT | 04 | 93_137 → 137_138 |
| [B14fd1a2758](boundaries/BG001_RMCA_B14fd1a2758/BG001_RMCA_B14fd1a2758_review.png) | MEVO_COMPONENT_ROOT | 05 | 92_144 → 144_145 |
| [B359e8f830b](boundaries/BG001_RMCA_B359e8f830b/BG001_RMCA_B359e8f830b_review.png) | PROXIMAL_M1_TO_MEVO | 05 | 92_144 → 144_146 |
| [Bb88d50380b](boundaries/BG001_RMCA_Bb88d50380b/BG001_RMCA_Bb88d50380b_review.png) | PROXIMAL_M1_TO_MEVO | 05 | 92_144 → 144_145 |
| [Bf58fac9b08](boundaries/BG001_RMCA_Bf58fac9b08/BG001_RMCA_Bf58fac9b08_review.png) | MEVO_COMPONENT_ROOT | 05 | 92_144 → 144_146 |
| [B086ac3307a](boundaries/BG001_RMCA_B086ac3307a/BG001_RMCA_B086ac3307a_review.png) | DISTAL_MEVO_TO_TERMINAL | 01 | 101_105 → TERMINAL |
| [B203ba3aa68](boundaries/BG001_RMCA_B203ba3aa68/BG001_RMCA_B203ba3aa68_review.png) | DISTAL_MEVO_TO_TERMINAL | 01 | 102_103 → TERMINAL |
| [B647b716f2c](boundaries/BG001_RMCA_B647b716f2c/BG001_RMCA_B647b716f2c_review.png) | DISTAL_MEVO_TO_TERMINAL | 01 | 96_106 → TERMINAL |
| [Be6a082a350](boundaries/BG001_RMCA_Be6a082a350/BG001_RMCA_Be6a082a350_review.png) | DISTAL_MEVO_TO_TERMINAL | 01 | 98_100 → TERMINAL |
| [Beb6b35bc52](boundaries/BG001_RMCA_Beb6b35bc52/BG001_RMCA_Beb6b35bc52_review.png) | DISTAL_MEVO_TO_TERMINAL | 01 | 98_99 → TERMINAL |
| [Bfb28a793d7](boundaries/BG001_RMCA_Bfb28a793d7/BG001_RMCA_Bfb28a793d7_review.png) | DISTAL_MEVO_TO_TERMINAL | 01 | 102_104 → TERMINAL |
| [Bbeaa20b6c9](boundaries/BG001_RMCA_Bbeaa20b6c9/BG001_RMCA_Bbeaa20b6c9_review.png) | DISTAL_MEVO_TO_TERMINAL | 02 | 107_111 → TERMINAL |
| [B24e8a4201e](boundaries/BG001_RMCA_B24e8a4201e/BG001_RMCA_B24e8a4201e_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 125_127 → TERMINAL |
| [B36f87bb7b8](boundaries/BG001_RMCA_B36f87bb7b8/BG001_RMCA_B36f87bb7b8_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 122_124 → TERMINAL |
| [B3df620f7bc](boundaries/BG001_RMCA_B3df620f7bc/BG001_RMCA_B3df620f7bc_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 117_118 → TERMINAL |
| [B3f21f48138](boundaries/BG001_RMCA_B3f21f48138/BG001_RMCA_B3f21f48138_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 125_126 → TERMINAL |
| [B67d7bacc8a](boundaries/BG001_RMCA_B67d7bacc8a/BG001_RMCA_B67d7bacc8a_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 115_116 → TERMINAL |
| [B6af69813b7](boundaries/BG001_RMCA_B6af69813b7/BG001_RMCA_B6af69813b7_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 120_121 → TERMINAL |
| [B76bb5667af](boundaries/BG001_RMCA_B76bb5667af/BG001_RMCA_B76bb5667af_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 122_123 → TERMINAL |
| [B7a4b65ef1d](boundaries/BG001_RMCA_B7a4b65ef1d/BG001_RMCA_B7a4b65ef1d_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 128_129 → TERMINAL |
| [B94f95566d7](boundaries/BG001_RMCA_B94f95566d7/BG001_RMCA_B94f95566d7_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 132_133 → TERMINAL |
| [Baa5c06b3cc](boundaries/BG001_RMCA_Baa5c06b3cc/BG001_RMCA_Baa5c06b3cc_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 131_135 → TERMINAL |
| [Bb8aa111741](boundaries/BG001_RMCA_Bb8aa111741/BG001_RMCA_Bb8aa111741_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 117_119 → TERMINAL |
| [Bc4dad6a3e6](boundaries/BG001_RMCA_Bc4dad6a3e6/BG001_RMCA_Bc4dad6a3e6_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 130_136 → TERMINAL |
| [Be8a71d5760](boundaries/BG001_RMCA_Be8a71d5760/BG001_RMCA_Be8a71d5760_review.png) | DISTAL_MEVO_TO_TERMINAL | 03 | 132_134 → TERMINAL |
| [B14a4147b08](boundaries/BG001_RMCA_B14a4147b08/BG001_RMCA_B14a4147b08_review.png) | DISTAL_MEVO_TO_TERMINAL | 04 | 137_143 → TERMINAL |
| [B36db3ce01f](boundaries/BG001_RMCA_B36db3ce01f/BG001_RMCA_B36db3ce01f_review.png) | DISTAL_MEVO_TO_TERMINAL | 04 | 140_141 → TERMINAL |
| [B8e0a3e3f84](boundaries/BG001_RMCA_B8e0a3e3f84/BG001_RMCA_B8e0a3e3f84_review.png) | DISTAL_MEVO_TO_TERMINAL | 04 | 138_139 → TERMINAL |
| [Ba6382479c7](boundaries/BG001_RMCA_Ba6382479c7/BG001_RMCA_Ba6382479c7_review.png) | DISTAL_MEVO_TO_TERMINAL | 04 | 140_142 → TERMINAL |
| [B27750a65a0](boundaries/BG001_RMCA_B27750a65a0/BG001_RMCA_B27750a65a0_review.png) | DISTAL_MEVO_TO_TERMINAL | 05 | 144_145 → TERMINAL |
| [B7578ebefeb](boundaries/BG001_RMCA_B7578ebefeb/BG001_RMCA_B7578ebefeb_review.png) | DISTAL_MEVO_TO_TERMINAL | 05 | 144_146 → TERMINAL |

## F. 高优先级审查

共有 **18 项 HIGH_PRIORITY_REVIEW**。包括全部近端 M1→MeVO、UNKNOWN→MeVO、每个组件首层分支，以及存在时的逆序。支持较低的 MeVO/UNKNOWN 转换按缓存 known fraction、confidence 排在同类前部；不因排序改变其标签或审查状态。对应图片索引见 `HIGH_PRIORITY_REVIEW.md`。

特别应检查 part03 的两个 UNKNOWN→MeVO 入口及其共同分叉位置。此处只提出审查问题：供体支持、相邻节点投票与组件入口的关系是否足够清楚？没有自动认定该候选正确或错误。

## G. VascularMD 状态

part01、02、04、05 已在此前完成原生建模及 smooth SWC/VTK/STL。part03 的原生入口检查失败：其根为分叉且没有符合既有规则的直接上游 M1 context。**VascularMD entry failure does NOT imply semantic boundary failure.** 本次不重跑模型，不改变 context，不产生新血管表面。QC 和原运行日志已独立复制入包。

## H. 已知限制

这些是 TopBrain-informed candidate labels，不是 BraVa 人工分段真值；Pilot 指标不能当作 BG001 准确率。既有 OOD 支持阈值未被本次重新优化。类型4的左右命名是 BG001 特定来源推断。PNG 是中心线证据图，不能用于精确量取管腔外表面；半径以原始数值为准。所有方向箭头只代表 native parent→child。终端也不意味着已确认的 M3 末端。

模型与制造适用性均不代替语义审查。所有真实事件均为 UNREVIEWED。缺数据或渲染错误的完整列表位于 `bundle_generation_summary.json`；本次 warnings：无。Errors：无。
