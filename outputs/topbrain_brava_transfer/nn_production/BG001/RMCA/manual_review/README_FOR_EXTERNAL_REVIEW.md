# BG001 RMCA MeVO 候选边界：外部审查入口

本包包含一个锁定版本的 BraVa BG001 右侧大脑中动脉及其 TopBrain 引导的语义候选。MeVO = M2 OR M3；M1 为近端类别，MeVO 为合并后的 M2/M3 候选，UNKNOWN 表示当前支持不足或分支判别未满足固定阈值。标签不是人工 ground truth。本次只整理已经冻结的结果，没有重新配准、投票、聚合分支、提取 ROI 或建模。

TYPE 4 = RMCA 采用 **BG001-specific DERIVED_SOURCE_SUPPORTED** 约定：由官方 BG001 具名长度表推断，并由用户明确选定。仅绑定本包原始和 ColorCoded SWC 的哈希，不推广至其他 BraVa。生产记录中的历史左右命名冲突仍完整保留在 evidence；本次审查状态命名不表示获得了新的官方直接 TYPE 字典。

## 建议阅读顺序

1. 阅读 `EXTERNAL_REVIEW_SUMMARY.md`，再查看 `global/` 中的三类总览。边界总览另提供相差 90° 的第二视图。
2. 优先阅读 `HIGH_PRIORITY_REVIEW.md`。每项链接到固定四面板图片和机器可读证据。
3. 图 A：完整 RMCA 为淡灰色，当前组件突出，红色圆形 marker 为当前边界。图 B：最多两条上游分支及两层下游分支；图 C：同一区域第二视角。图 D：直接来自缓存的概率、支持距离、半径、长度与供体证据。
4. PROXIMAL 表示靠近原生 RMCA 根的上游侧，DISTAL 表示沿原生 parent→child 的下游侧。箭头遵从 SWC 父子方向，不能等同于已经测量的生理血流方向。
5. 对照 `tables/` 的 CSV 和缩进拓扑。必要时使用单边界的 metrics/local_topology JSON；VTP 仅是可选的第二级证据。

每张图有固定颜色图例：M1 橙色，MeVO 青色，UNKNOWN 紫色，当前边界红色，背景 context 灰色。组件总览另以五种颜色区分 part01–part05。中心线按原样绘制，线宽是显示用宽度；真实半径以 CSV/JSON 数值为准，不以 PNG 线宽测量。

## 边界类型

- **PROXIMAL_M1_TO_MEVO**：上游 M1 分支与下游 MeVO 分支之间的候选入口。
- **DISTAL_MEVO_TO_UNKNOWN**：MeVO 结束并进入 UNKNOWN 的候选远端。
- **DISTAL_MEVO_TO_TERMINAL**：MeVO 分支在原始 SWC 终端结束；不等同于已确认的 M3 末端。
- **UNKNOWN_TO_MEVO**：UNKNOWN 下游进入 MeVO；需审查支持空缺及入口依据。
- **MEVO_TO_M1_REVERSE**：MeVO 下游重新出现 M1 的逆序候选；最高优先级。
- **MEVO_COMPONENT_ROOT**：每个严格连通组件的首层 MeVO 分支；分叉根可对应多条首层分支。
- **DETACHED_MEVO_COMPONENT**：组件入口没有直接上游 M1；不表示原始血管断裂，也不判定语义错误。

一个物理分叉点可同时具有入口、组件根等多个事件；分叉根的每条首层 MeVO 分支也单独保留事件。ID 由 subject、side、上游分支、下游分支和类型确定，与 CSV 行号无关。因此事件数不等于空间位置数。总览 marker 后的 `+` 表示同一位置有多个事件，旁侧索引列出同位置的完整短 ID；每个事件都有独立图片。

## 数值和状态的边界

概率、known fraction、confidence、分支支持距离和有效 donor 数全部来自生产缓存。组件平均概率仅是对这些既有分支指标按已存长度作描述性汇总，不改变任何分支标签。radius ratio = 下游分支原始样本半径中位数 / 上游对应中位数。根距离及最近分叉距离均沿原生中心线弧长测量；最近分叉允许包含边界点自身。

供体 before/after 证据读取共享节点前、后最近的**不同原始 SWC 样本**，同时给出节点编号；不是重新发起 NN 查询。它可能不同于全分支的平均支持。终端不存在下游分支，相关字段为 null，CSV 对应空单元；这不是零概率。任何缺字段会明确列为 FIELD_NOT_AVAILABLE，不猜值。

所有事件的 review_status 保持 **UNREVIEWED**。高优先级仅表示应优先审查，不表示 ACCEPT、REJECT 或已确定的 M2 起点。

**VascularMD entry failure does NOT imply semantic boundary failure.** part03 是有效保留的语义候选区域，其建模入口检查失败与语义边界是否合理完全分开。本包不会把建模失败自动转为语义否定。所有几何仍为 ANATOMICAL_TRANSFER_CANDIDATE，制造状态为 MANUFACTURING_NOT_VALIDATED。

## 文件结构与离线使用

```text
README_FOR_EXTERNAL_REVIEW.md / EXTERNAL_REVIEW_SUMMARY.md / HIGH_PRIORITY_REVIEW.md
review_manifest.json / bundle_generation_summary.json / bundle_generation.log
tables/       边界、组件、分支 CSV；55 分支拓扑 TXT/JSON；根与组件分离原因
global/       标签、组件、边界总览及边界第二视角
components/   五个组件总览
boundaries/<boundary_id>/  review PNG、metrics JSON、local_topology JSON、VTP+显示说明
evidence/     两个 BG001 源 SWC、原始预测和映射、小型支持 NPZ、严格/modelable SWC、原生 QC 与日志
```

生产表中的绝对路径用于溯源；查看本包无需这些路径存在。严格/modelable SWC 的 review 表使用包内相对路径。大型 VascularMD surface/STL 仅记录原保存路径，未复制；其状态由包内 QC 和日志说明。未包含 TopBrain NIfTI、供体点云、其他 BraVa 个体或完整 VascularMD 缓存。PNG、CSV、JSON 和源中心线足够支持第一轮边界审查。

`review_manifest.json` 对每个普通文件记录 SHA256、用途、优先级和关联事件/组件。manifest 本身使用“将自己的 SHA256 字段设为 null 后的 canonical JSON”的摘要规则，避免自引用哈希不可能固定的问题；完整 ZIP 摘要位于包旁的 `.sha256` 文件。ZIP 内路径均为相对路径。
