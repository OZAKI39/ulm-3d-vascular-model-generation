# 学位论文式总结的来源与取舍记录

本文件属于内部证据追溯，不是论文正文。整理截止于本轮已保存材料。正文以真实物理问题组织，不按 Particle 编号编年；本表保留阶段名称以便查证。本文的“已建立”指在注明范围内已有实现和验证证据，不将自动检查通过改写为所有生理目标或用户人工验收完成。

## 证据优先级与本轮边界

优先采用当前正式轨迹队列的根级 VALIDATION 与其数据目录，其次采用各阶段科学合同、正式机器记录及配套 REVIEW；开发试算、参考重放和历史较小队列不混入正式统计。若最新记录仍写有未收敛、未建立或用户审核待定，正文保留这些限制。旧结果只在解释模型修正时使用，并说明比较对象。

本轮只整理既有研究和制作说明性图，不重新求解血流，不运行新的粒子轨迹，不恢复已停止的远程任务。Particle-8.2 尚未整体完成；其中间计算、性能数字和延长计算结果均未用于正文“已完成”结果或图 16。第 17 节仅陈述待完成的并行与多出口研究。源目录现存三处未提交代码改动及未跟踪的 Particle-8.2 报告均保留。

原始文献保留项目采用时的版本。文中具体研究数值来自本地证据，不从此前对话记忆抄录。正文对结果作合理舍入；高精度值与来源字段已另存 [extracted_values.json](audit/extracted_values.json)。

## 章节与数值证据

下列路径均相对于仓库根目录。配套 [证据散列清单](audit/evidence_inventory.json)包含 32 份阶段审核、验证、文献和合同文件的实际散列。字段名称仅在本文件出现。

| 章节 / 数值 | 主要来源 | 取舍及解释 |
|---|---|---|
| 第 2 节：70,363 节点、371,402 四面体 | `particle_3d/reports/particle0/PARTICLE0_VALIDATION.json`：`node_count`、`tetra_count` | 采用正式冻结体网格，不用可视化抽稀后的网格 |
| 第 2 节：45,221 壁面三角形 | `particle_3d/reports/particle3/PARTICLE3_VALIDATION.json`：`wall_triangle_count`；P8.1 `full_geometry_wall_triangles` | 两份记录一致；开放端盖另外保留，不计为固体墙 |
| 第 2、16 节：入口流量与出口分配 | `particle_3d/reports/particle7/PARTICLE7_VALIDATION.json`：`inlet_Q_m3_s`、`outlet_Q_m3_s`、`mass_balance`；`particle_3d/outputs/particle8_1/PARTICLE8_1_VALIDATION.json`：`clock`、`outlet_frozen_volume_flux_m3_s`、`outlet_frozen_volume_flux_fraction` | P7 与 P8.1 仅末位舍入差；论文流量比例及图 16 以 P8.1 正式记录为准，整体残差采用 P7 原积分记录，不重新调整出口流量 |
| 第 2、9 节：黏度 0.00345312 Pa·s | `particle_3d/reports/particle5/PARTICLE5_VALIDATION.json`：`dynamic_viscosity_pa_s`；冻结 `baseline_summary.json` 与 `run/solver.xml` | 沿用原流场物性，无新参数选择 |
| 第 3 节：解析误差及真实节点回采样误差 | P0 VALIDATION：`max_affine_velocity_error`、`max_affine_pressure_error`、`max_affine_gradient_error`、`max_node_velocity_error` | 人工解析真值与真实节点恢复分开，不能当作 FEM 生理真值误差 |
| 第 4 节：1.176 µm；91.75/91.83/91.90 ms；69.10 µm | `particle_3d/reports/particle1/PARTICLE1_VALIDATION.json`：`single_mb_diameter_um`、`real_trajectory_cases` | 表述为入口邻近内部起点；路程与代表值采用三个验证步长中最细一档，未称其为生产步长 |
| 第 5 节：直径 6.79±0.93 µm；体积均值 47.9 fL、CV 0.148 | `particle_3d/contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json`；`particle_3d/reports/particle2/literature/sources.json` | 直径和 MCV 来自项目采用的 Moss 预印本；CV 是项目选择。文献 2.6 fL 为鼠间 MCV 离散度，不是单细胞 SD；未将其代入体积分布 |
| 第 5 节：100,000 样本、接受后 6.80 µm / 47.87 fL | `particle_3d/reports/particle2/PARTICLE2_VALIDATION.json`：`validation_population_N`、`accepted_D_statistics`、`accepted_V_statistics` | 与潜在均值不同是条件化后的模型样本差异；不“纠正”回输入参数 |
| 第 6 节：姿态关系与形状比 | 上述 RBC 合同；`particle_3d/src/particle_3d/rbc_orientation.py`；P2 VALIDATION 的静止/旋转/剪切记录 | 正文把原代码方向记号 p 改为数学符号 s_R，避免与压力 p 冲突；短轴形状周期与有向向量周期分开 |
| 第 7、8 节：壁面及六类形状接触 | P3、P4 REVIEW 与 VALIDATION；`wall_gap.py`、`pair_geometry.py`、`wall_contact.py`、`kinematic_contact.py` | 硬约束与实际表面几何；未把椭球变成等效球，也未称接触乘子为真实测量力 |
| 第 9 节：3.233 nm 局部样本、速度比例 0.00255 | P5 VALIDATION：`real_near_wall_mb_gap_m`、`real_near_wall_attenuation_ratio` | 采用靠近壁面样本；更小的原始间隙可能对应离墙，不能混作靠近验证 |
| 第 10 节：0.01/0.05、2 nm、10⁻³ 相对尺度 | `particle_3d/contracts/NEAR_FIELD_REGULARIZATION_V1.json`；P6.5 文献审计 | 最新正则化合同优先于旧 P5 的未截断 1/h；参考半径与阻力约化半径分开 |
| 第 10 节：0.0235 nm 历史 / 0.0815 nm 同初态 / 约 2 nm 新模型 | `particle_3d/reports/particle6_5/PARTICLE6_5_VALIDATION.json`：`real_old_min_gap_m`、`real_matched_p5_min_gap_m`、`real_v1_min_gap_m` | 三者的案例身份明确区分；不把历史双球和重放单球写成同一初态对照 |
| 第 10 节：基准 6 次、细步长 0 次、2.038/3.233 nm | P6.5 VALIDATION 的 `real_replay`、`sensitivity_results` 与 REVIEW 第 8、9 节 | 安全性成立但事件时间未收敛，事件次数是判定带诊断 |
| 第 11 节：H_D=0.45、C_MB=8.5×10¹² m⁻³ | `particle_3d/contracts/PARTICLE7_INLET_POPULATION_V0.json`；P7 `literature/sources.json` | 入口体积输送比例及名义充分混合浓度；不是管内瞬时压积或实测 PK |
| 第 11 节：0.0233 s⁻¹、42.99 s、1.23×10⁻¹⁵ m³/s | P7 VALIDATION：`mb_number_rate_s_inv`、`mean_mb_event_interval_s`、`rbc_volume_rate_m3_s` | 由已存通量关系导出；计划事件与接纳分开 |
| 第 12 节：状态承载及重启一致 | P6 REVIEW、VALIDATION；P7 `restart`；P8 VALIDATION 的 `restart` | P6 负责已有运动的一致性与保存恢复；P7/P8 再验证动态人口，不扩写成新增悬浮液物理 |
| 第 13 节：1110/1/1109 RBC，1/0/1 MB | `particle_3d/reports/particle8/PARTICLE8_VALIDATION.json`：`real_geometry.real_mixed_inlet_smoke_counts` | 采用明确标为静态接纳障碍的真实入口场景，不用人工开放区的百万级人口证明真实 passage |
| 第 14 节：2200/1969/1249、231、720 | P8.1 根 VALIDATION：`scheduled`、`admitted`、`completed`、`end_reasons` | 根级正式队列；排除 `reference_*` 与步长敏感性重算，保持两个群体守恒等式 |
| 第 14 节：停留时间及路程统计 | 同上：`completed_residence_time_s`、`completed_path_length_m` | 仅完成样本；从 SI 换算为 ms/µm 后舍入，不把停止时刻当退出时刻 |
| 第 14 节：26.27 h、0.25 ms | 同上：`acquisition_birth_window_s`、`integration_config.dt_s` | 采集供给时钟不是壁钟耗时或单次注射有效期；步长属于本数据集验证设置 |
| 第 15 节：150×150 格、前一区间时间加权 | `particle_3d/src/particle_3d/particle81_figures.py` 第 99–104 行；P8.1 REVIEW | 实际为 XY/XZ/YZ 二维直方图，值为 s/bin；原请求建议的三维密度仅作连续理想解释，不能误写为已经计算了三维体素场 |
| 第 16 节：0/1249/0 与非零出口流量 | P8.1 根 VALIDATION | 最新完成队列仍只有 1/3 出口覆盖，不引用未完成 P8.2 来改善结论 |

## 历史阶段如何进入论文

| 阶段 | 论文去向 | 本文认可的已完成范围 | 保留的边界 |
|---|---|---|---|
| Particle-0 | 2–3 节 | 冻结输入与任意点插值验证 | 无新 FEM 生理验收 |
| Particle-1 | 4 节 | 单微泡基础运动及完整验证路径 | 验证步长、入口邻近内部初始化 |
| Particle-2 | 5–6 节 | RBC 统计几何与自由姿态 | 无真实膜变形及有限尺寸通行证明 |
| Particle-3 | 7、13 节 | 壁面间隙、连续安全性、简化变形可行性 | RBC 真实 passage 未建立 |
| Particle-4 | 8 节 | 有限尺寸配对与同时接触验证 | 几何运动学，不是完整悬浮液 |
| Particle-5 | 9–10 节 | 球形法向近场与阻力响应验证 | 旧亚纳米结果由新版解释边界替代；非球形未建立 |
| Particle-6 | 12 节 | 同物理运动在状态承载、邻近查询及重启中的一致性 | 无新增粒子物理，无多核扩展结果宣称 |
| Particle-6.5 | 10 节 | 平滑激活及连续介质下限安全性 | 交接时间、次数和步长稳定性未建立 |
| Particle-7 | 11–13 节 | 通量驱动人口与生命周期控制 | 人工长程开放区不冒充真实血管通行 |
| Particle-8 | 12–15 节 | 既有状态的混合场景与生命周期回放 | 真实 RBC 仍是受限接纳 |
| Particle-8 可视化增强 / full3d | 15 节 | 完整几何、三维旋转与保存状态回放 | 相机动画不增加轨迹或物理证据 |
| Particle-8.1 | 14–16 节 | 1249 条正式完成轨迹及累积可视化 | 用户审核仍待定；生理完成未建立；分支不完整 |
| Particle-8.2 | 17 节 | 不作为已完成成果 | 已暂停的未完成专项，只列后续研究 |

P8.1 正式状态明确为自动检查通过、显示有局限、完整网络分支覆盖不完整、生理完成未建立，且 `user_acceptance=PENDING_USER_REVIEW`。本轮未将该状态改成用户已验收。

## 文献核查与适用性

- RBC 参数采用本地已保存的文献表格索引及对应合同，不因新写总结而重新选择参数；预印本身份按项目采用版本保留。
- 本轮补核 ULM 背景的 Errico 原始 Nature 页面、Ness 论文的爱丁堡大学原文 PDF，以及 Patel、Dauba 的原始研究题录。没有用搜索摘要替代项目数值结果。
- Jeffery 出版社页面本轮返回访问错误；题录 DOI 及项目已实现公式相互核对，未宣称本轮读取全文。部分 PMC 页面出现浏览器验证，因此沿用项目既有表格位置及提取记录，不声称本轮重新获取全部实验表。
- 2 nm 为项目选择。Ness 只支持相对尺度先例；受限水、脂质水合和 SonoVue 壳层文献只支持适用尺度讨论，不支持直接测得血管 2 nm 截断。
- 参考文献 [4]、[5] 为 RBC 体积尺度交叉核对，不是新的混合品系统计总体；正文未用它们修改已经采用的分布。

## 图像与审计范围

主图 16 张中，14 张逐字节复制历史结果图；图 1 为本次方法路线示意，图 16 根据 P8.1 既有流量与完成计数重绘。没有生成新数值轨迹。历史图中的阶段标识和符号保留，中文图注阐明科学内容，符号换记见配套符号表。

本轮对文稿、公式、图像路径、来源散列和两种正文版本一致性进行检查，没有重跑粒子科学测试。历史测试通过数只能归属于原运行；详见 [验证状态索引](audit/stage_validation_index.json)。LaTeX 无可用编译器，静态检查不等于 PDF 编译通过。
