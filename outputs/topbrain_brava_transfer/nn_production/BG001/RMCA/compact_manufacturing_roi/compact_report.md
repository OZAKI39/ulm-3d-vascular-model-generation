# 单一冻结 MeVO 组件的紧凑制造 ROI 研究

## 摘要

本研究将既有大范围制造模型保留为 FULL_CONTEXT_REFERENCE，并从单一冻结语义组件生成 MINI、BALANCED 与 RICH 三个紧凑派生物。实际来源为 part04，默认候选为 BALANCED。未重新执行配准、供体排序、NN、UNKNOWN 校准、分支语义聚合或近端精修。所有保留坐标、TYPE 与有向边均来自真实 BraVa；半径补偿和打印旋转只属于制造层。当前状态为 PRINT_READY_CANDIDATE，未进行实体打印。

## 方法

首先检查配置优先的 part04，并同时读取已有 semantic_refined 与 diameter_075 单根文件。以每轮三个候选均具有最小分叉结构为来源可行性条件，优先使用配置中的组件顺序，并以冻结语义支持减去需补偿中心线比例作为工程评分。程序未对原语义进行重新投票或聚合。

从精修根沿原始拓扑向远端展开，分支按冻结语义支持、原始中位直径及原长度的 0.45/0.35/0.20 权重排序。最多保留两个女儿，最多跨两次解剖拓扑分叉；分叉深度仅描述制造简化，不对应 M2/M3 分段。路径在预算以内的最后一个原始样本截断，不插值、不连接其他组件。删除短于 12 mm 且原始中位直径低于 1 mm 的末端 twig，其余必要连接路径保留。每条未完全保留的原语义分支均在 branch_selection.csv 中解释。

近端只选择距目标 10 mm 最近、位于 8–15 mm 内的真实上游前缀，不增加兄弟分支。根据用户后续明确授权，45/55/65 mm 路径长度及 60/80 mm 全局长度仅保留为参考，不作为硬截断。算法对原始采样点距离与全部拓扑边界进行确定性搜索，在单组件、分支/出口/深度、200 mm 总中心线、最大空间范围等硬限制内，权衡空间目标、代表性分叉、总长度及补偿比例，自动确定各候选截点。RICH 的空间范围仍服从全局最长边 110 mm、对角线 140 mm 上限。

制造半径下限为 0.60 mm，使用已有约 5 mm 局部 smoothstep 过渡；每点原始半径、补偿半径和增量保存在映射表。补偿比例按边长及两端补偿指示量的梯形权重统计，另列采样点比例。半径倍数超过 1.75 时标注 LARGE_MANUFACTURING_COMPENSATION，不自动抹平原始半径差异或判为解剖改变。

## 来源可行性与尺寸冲突

真实 part04_01 的第一处分叉距离精修根约 58.62 mm，第二处分叉约 123.05 mm；固定 55 mm 上限会在首处分叉之前截止。用户因此明确允许自适应路径长度，并要求优先保持合理复杂度与空间范围。本研究以该授权为准，未偷偷延长默认路径或更改语义边界。

自适应搜索后实际选择 part04_01 的 semantic_refined 层，使用 part04=True，回退 part03=False。既有 d075 也被检查，但其短细截断末端进一步筛除后不能提供三个模式都需要的最低分叉结构；因此采用同一组件的冻结语义版，并只在制造派生物中局部截断。

本次尺寸目标作为偏好，allow_undersized_candidate=True。结构、空间上限与切片可行性分别检查；不通过坐标缩放、跨组件连接或增加长 M1 前缀扩大模型。

| 来源 | 层 | BALANCED 分支/分叉/出口 | 结构评价 |
|---|---|---|---|
| part01 | semantic_refined | 3/1/2 | MINI_MINIMUM_STRUCTURE_UNAVAILABLE |
| part01 | diameter_075 | 1/0/1 | BIFURCATION_BUDGET, OUTLET_BUDGET, MINI_MINIMUM_STRUCTURE_UNAVAILABLE, RICH_MINIMUM_STRUCTURE_UNAVAILABLE |
| part02 | semantic_refined | 1/0/1 | BIFURCATION_BUDGET, OUTLET_BUDGET, MINI_MINIMUM_STRUCTURE_UNAVAILABLE, RICH_MINIMUM_STRUCTURE_UNAVAILABLE, LOW_FROZEN_SEMANTIC_SUPPORT |
| part03 | semantic_refined | 5/2/3 | 可选 |
| part03 | diameter_075 | 3/1/2 | MINI_MINIMUM_STRUCTURE_UNAVAILABLE |
| part04_01 | semantic_refined | 5/2/3 | 可选 |
| part04_01 | diameter_075 | 1/0/1 | BIFURCATION_BUDGET, OUTLET_BUDGET, MINI_MINIMUM_STRUCTURE_UNAVAILABLE, RICH_MINIMUM_STRUCTURE_UNAVAILABLE |
| part04_02 | semantic_refined | 1/0/1 | BIFURCATION_BUDGET, OUTLET_BUDGET, MINI_MINIMUM_STRUCTURE_UNAVAILABLE, RICH_MINIMUM_STRUCTURE_UNAVAILABLE |
| part05 | semantic_refined | 3/1/2 | 可选 |

## 三种紧凑候选的实际结果

| 候选 | 分支 | 分叉 | 出口 | 长度 mm | AABB mm | 最长边 mm | context mm | 补偿长度比例 |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| MINI | 3 | 1 | 2 | 146.239 | 29.14 x 27.90 x 68.20 | 68.200 | 10.343 | 49.40% |
| BALANCED | 5 | 2 | 3 | 173.694 | 29.14 x 27.90 x 70.68 | 70.680 | 10.343 | 54.62% |
| RICH | 5 | 2 | 3 | 173.694 | 29.14 x 27.90 x 70.68 | 70.680 | 10.343 | 54.62% |

MINI 原生建模：success；拟合后最小直径 1.1703 mm；最大补偿倍数 2.299，超过 1.75 的采样点数 37；自适应根至末端截断距离 121.382 mm，加入 context 后最长路径 131.725 mm。


BALANCED 原生建模：success；拟合后最小直径 1.1703 mm；最大补偿倍数 2.299，超过 1.75 的采样点数 46；自适应根至末端截断距离 136.539 mm，加入 context 后最长路径 146.882 mm。


RICH 原生建模：success；拟合后最小直径 1.1703 mm；最大补偿倍数 2.299，超过 1.75 的采样点数 46；自适应根至末端截断距离 136.539 mm，加入 context 后最长路径 146.882 mm。


旧模型有 39 条拓扑分支、19 个分叉、20 个出口，包围盒为 62.93 x 80.29 x 91.14 mm。BALANCED 为 5 条分支、2 个分叉、3 个出口，总长度减少 87.78%；出口减少 17，分叉减少 17。

全部候选为单连通、单入口。BALANCED/RICH 总长度超过 160 mm 优选区间但仍低于 200 mm 硬上限；这是保留第二处分叉及其真实末端支路的结果。若 RICH 与 BALANCED 的补偿 SWC 哈希相同，会在清单记录 same_geometry_as，不虚构第三处分叉。

## 打印方向与实际切片

实际打印配置为 Bambu Lab P1S，喷嘴 0.4 mm，材料为任务明确指定的 ABS，打印体积为 256.00 x 256.00 x 250.00 mm，软件版本 02.07.01.57。原有 GUI 预设保持不变。

仅对 BALANCED 表面搜索 73 个确定性方向，保存前三名；MINI、RICH 仅输出表面，不自动切片。方向评分增加入口高度、近端主干与床面 30–60° 的偏好，以及远端长水平段惩罚。最优方向 XYZ 欧拉角为 [86.76752772701319, 7.499175533260423, -138.9437244233769] 度，打印包围盒为 54.29 x 57.85 x 34.39 mm，入口高度 0.433 mm，主干角度 26.55°。

向下需支撑面积启发式为 100.250 mm²；该指标不是切片器支撑体积。实际只调用 BALANCED Top1 切片，状态为 BAMBU_SLICED。

切片记录：状态 BAMBU_SLICED；时间 897.0 秒；总耗材 1.8 g。支撑信息为 {'support_used': 'true', 'support_type': 'tree(auto)', 'filament_usage_g': None, 'quantity_source': 'Separate support quantity not exposed by this archive', 'feature_tags': ['Support', 'Support interface']}；未提供的独立支撑质量不补造。实际原始元数据以 BambuStudio/slice_results.json 为准。

## 讨论与可复核性

三种模型用于比较同一局部血管网络的制造复杂度，不能解释为不同解剖标签。所选 part04_01 本身只有两个分叉，因此 BALANCED 与 RICH 可能收敛为相同拓扑和几何；不会为制造一个更大的 RICH 而连接 part04_02、part03 或其它组件。三次二叉分叉加一个入口主干通常需要七条拓扑分支，与最多六条的约束也不相容，因此三分叉是上限，不是强制目标。自适应路径可以超过旧固定参考长度，但总长度与空间范围仍受硬上限限制。阈值为本项目工程启发式，NOT manufacturer guaranteed limits. Must be calibrated experimentally.

旧 manufacturing_roi、semantic_refined、diameter_075、既有 VascularMD 输出和上游库均保持原文件哈希；旧大模型的 FULL_CONTEXT_REFERENCE 身份记录在新清单内，不改写旧清单。s1-2 可视化代码未修改。没有向打印机提交任务，不能称为 MANUFACTURING_VALIDATED。

本次明确记录的工程提示：MINI: LARGE_MANUFACTURING_COMPENSATION; BALANCED: LARGE_MANUFACTURING_COMPENSATION; RICH: LARGE_MANUFACTURING_COMPENSATION; MINI: VMD_PREFERRED_FLOOR_UNDERSHOOT; hard minimum checked separately; BALANCED: VMD_PREFERRED_FLOOR_UNDERSHOOT; hard minimum checked separately; BALANCED: ABOVE_PREFERRED_TOTAL_LENGTH; below hard 200 mm cap; RICH: VMD_PREFERRED_FLOOR_UNDERSHOOT; hard minimum checked separately; RICH: ABOVE_PREFERRED_TOTAL_LENGTH; below hard 200 mm cap; STEM_ANGLE_OUTSIDE_PREFERRED_RANGE: soft preference traded against inlet height and support。

本轮核验 1077 个受保护文件，哈希全部保持一致。测试及逐点几何核验的完整记录见 QC/test_results.json。

实际回归共 202 项通过，其中既有 190 项、新增 12 项，最终失败 0 项。既有 scikit-image/NumPy 弃用提示 61 条，新增测试无警告。六个导出 SWC 均通过原坐标、TYPE、有向边和原始半径映射核验；三套原生模型均成功且保持拓扑，三套打印表面均为单连通、封闭且无非流形边。

s1-2_swc_roi_generate_human.py 的 SHA256 为 `3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729`，与保护快照一致。旧大模型的输出及其九个实现文件哈希也保持不变。
