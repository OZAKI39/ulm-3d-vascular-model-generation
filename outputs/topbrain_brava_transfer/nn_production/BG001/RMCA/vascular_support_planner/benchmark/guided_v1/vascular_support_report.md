# BG001 RMCA 血管专用支撑规划报告

当前状态：**READY_FOR_HUMAN_REVIEW**。本轮输入为最新 2 mm 五面盒与既有血管、端口；上游模型 SHA256 复核结果：True。

## 为什么原支撑需要调整

细血管的长悬空段可能在自动支撑中被遗漏；增加自动支撑密度又可能让树枝沿血管连续包裹。本轮将接触位置与支撑结构分开处理：先确定危险位置，再让成熟切片器生成支撑。真实 ABS 打印和机械拆除仍是后续实验。

## 检测与接触位置

检测到 **13 段危险悬空区域**，最长 **19.98 mm**。804 个候选位置被筛选为 45 个计划接触点。5 mm 是本次工程阈值，并非通用物理定律；4、5、6 mm 的敏感性结果另存为表。

Bambu 原生切片约有 **51 个弧长接触区间**；Guided Orca 实际约有 **59 个**。此处是从挤出刀路投影到血管中心线后估计的断续接触区间，不等于切片器内部树叶节点数；连续包裹和网格接触片数量另外记录，避免混淆。

按这一共同估计口径，本次**未达成接触数量下降的目标**。筛选后的计划点数不能替代实际切片的接触区间数；本次改进主要体现在减少连续接触面积，不能宣称所有拆除指标均已改善。

## 必要位置是否仍有遗漏

Guided Orca 是否仍有危险段没有覆盖：**NO**。实际覆盖率 100.00%，最大未覆盖危险弧长 0.00 mm。计划覆盖与刀路覆盖分开计算；仅计划通过不能作为生产通过。

## 接触与拆除负担

Guided Orca 的接触面积相对 Bambu 变化：减少 57.48%。该面积由下表面三角形与实际支撑挤出轨迹的邻近关系估计；所有方案使用同一估计方法。

分叉 2 mm 避让区附近接触区间数：5。这表示靠近分叉，不自动等于支撑与分叉实体相交。预计被困支撑数：0。探针可达性是离散几何估计，能抓到树干不等于一定能完整抽出。

分叉完全避让尚未实现。计划中进入避让区的点逐项标记为覆盖例外：移除后会漏掉危险要求，保留理由见 support_contacts.csv。真实切片接触范围可能大于计划点，因此另行核验实际接触面；QC/08 从下方着色显示实际接触，不能把“碰撞检查通过”解读成“分叉没有接触”。

## 同条件比较

|方案|接触区间数（估计）|接触面积 mm²|危险覆盖率 %|被困风险数|切片时间估计 s|硬门槛失败|
|---|---:|---:|---:|---:|---:|---|
|BAMBU_NATIVE_BASELINE|51|237.14|99.83|0|35628.00|CRITICAL_COVERAGE_INCOMPLETE, SUPPORT_COLLISION_NOT_CLEARED|
|ORCA_NATIVE_BASELINE|38|235.71|99.83|1|34831.00|CRITICAL_COVERAGE_INCOMPLETE, TRAPPED_SUPPORT_RISK, SUPPORT_COLLISION_NOT_CLEARED|
|GUIDED_ORCA|59|100.84|100.00|0|34675.00|无|
|CUSTOM_CLEVER_TREE|尚未取得有效结果|尚未取得有效结果|尚未取得有效结果|1|尚未取得有效结果|CUSTOM_RESEARCH_NOT_PRODUCTION, SELF_SUPPORT_NOT_CONFIRMED, CUSTOM_TIP_COLLISION|

表中的 SUPPORT_COLLISION_NOT_CLEARED 表示保守挤出包络的间距下界未通过，不等同于已经证明发生实体碰撞；失败候选均未被推荐。接触面积与支撑质量均为刀路估计。自定义树的名义接触圆盘面积不可与这些实际刀路面积直接比较，所以其实际接触面积留空。

Orca 使用其 X1C 0.2 mm 官方配置结构承接本机保存的 Vessel 工艺值。层高、首层、喷嘴、线宽、顶部和 XY 间隙均逐项保持；不同引擎的启动宏使用各自官方实现，转换明细见 profiles/orca_translation_audit.json。

## 下一次打印选择与剩余风险

通过全部软件硬门槛后的排序推荐：**GUIDED_ORCA**。Guided Orca 项目文件可人工审核，STL 中的 enforcer/blocker 是控制区域，不是直接打印实体。自定义 Clever Tree 和短桥均为研究候选，其失败结果也保留。

最大剩余风险是软件中可支撑、可触达的细血管，在 ABS 冷却收缩和真实拆支撑受力时仍可能断裂。没有执行真实打印、上传或启动打印，也不声明制造验证完成。

本次推荐参数：{"support_type": "tree(manual)", "tree_support_tip_diameter": "0.5", "tree_support_top_rate": "20%", "tree_support_branch_diameter_organic": "2.0", "tree_support_angle_slow": "15"}。具体候选目录：`/home/lzy/projects/vascular_printing/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner/slicer/guided_00`。支撑估计用料 0.68 g；整件打印时间估计 34675.00 s，并非本次计算耗时。

## 方法与复现

四种检测证据分别来自表面法向、按实际层高切片的上下层承托、FreeFloating 图覆盖与局部低点、SWC/已冻结端口路径的弧长分析。Greedy 覆盖选择沿同一分支弧长计算，弯曲度和局部半径会缩短影响范围；不是固定欧氏球。各指标分别输出，候选先过硬门槛，再按接触面积等字典序排序。

`./run_python.sh s1-8_vascular_support_planner.py --stage all`

输入来源见 source_discovery.json；检测、接触、覆盖见 tables/ 和 JSON；12 个以内 Orca 搜索候选见 slicer/；阶段图见 QC/；输入哈希见 protected_before.json 和 source_protection.json；永久测试位于 tests/test_s1_8_vascular_support_planner.py。
