# F10R：仿真四PZT控制序列与XY导航交付

已完成三种仿真控制模式、原X精确重放、真实成对/分段实验、有限相位/初态对照、新Y稳定弱控制及四组图件。X短目标可达；新Y数值稳定但未达目标；没有真实硬件执行。

本报告来自实际保存的2 ms轨迹；单分辨率观察、完整数值准入、实验标定分别统计。`EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED`；`training_ready=false`。


## 原始身份与历史保护

最新公开同步HEAD为467b17411deeb13a77370565b7382f071b41a774，本地完整科学源码HEAD为41ae8ed2b88f2cbbed834fffccaa0b5de6e91589，950个源码/配置/测试/依赖文件逐项SHA相同。完整本地F10Q工作区干净，公开包与科研包的不同来自归档范围；最新清单位于sync/20261014，根SYNC_MANIFEST为Stage5旧清单。独立分支dev/fast-pzt-f10r-simulation-control-navigation-20261014，无安全废弃对象可删，清理审计单独提交3e55371324701bb32e9eec229ebd226cff46eb25。所有旧实现、失败及轨迹保留，独立只读核验无旧跟踪文件修改或删除。

| 控制 | 终点COM (µm) | 对应目标误差 (µm) | COM数值 | 形状数值 | 固定ID | 目标准入 |
|---|---|---|---|---|---|---|
| OFF | 1.58818678e-16, 0, 0 | 0.2 | PASS | PASS | PASS | 否 |
| 原X全精度重放 | 0.200896324, -1.09787383e-06, 1.25731453e-16 | 0.000896325125 | PASS | PASS | PASS | 是 |
| 同轴X对 | 0.131843603, -9.26442286e-16, -6.6174449e-18 | 0.0681563975 | PASS | PASS | PASS | 是 |
| 单片依次X−/X+ | 5.93980663e-09, 0, 1.32348898e-17 | 0.199999994 | NOT_RUN | NOT_RUN | NOT_RUN | 仅单分辨率/不准入 |
| 两个1 ms段 | 0.200674583, 5.86855223e-07, 1.41613321e-15 | 0.000674582901 | PASS | PASS | PASS | 是 |
| 新低压Y | -1.8507483e-09, 0.0129440325, 0 | 0.187055968 | PASS | PASS | PASS | 否 |

## 用户要求的22项回答

1. 是，从最新公开F10Q同步分支的准确HEAD建立独立worktree开发。
2. 已核实公开与本地完整科研版本的950个科学源码/配置等文件字节相同，提交不同不代表科学差异。
3. 无可证明安全删除对象，未强行删代码；审计和活跃入口索引已独立提交。
4. F1–F10Q历史证据完整保留；本轮170个远端文件回传并逐项SHA核验，原始失败不覆盖。
5. 已实现独立段号、起止、频率、启停、幅相、场身份及实际事件日志，非网格事件缩短步长准确对齐。
6. ONE_HOT_PZT、AXIS_PAIR_CONTROL、FOUR_CHANNEL_INDEPENDENT均已通过接口测试并实际用于轨迹。
7. 仍为2.5 MHz均匀水圆片代理与原合法APECSS，缓存复用，无频率外推或静默压力裁剪。
8. 严格区分Pa-equivalent与Vpp；原1500源上限不代表功放安全电压，未建立虚构换算。
9. 没有开发/运行真实仪器通信，也没有USB或实际电压输出。
10. 原X控制全精度重放，全部2001帧、32ID坐标逐位一致，最大差0 m。
11. OFF静止；原X、X对、两段X都达到原目标且通过完整分层检查；单片反向序列基本相消，未达目标，数据见上表。
12. ±5°团初始刚性旋转后观察到X约0.200894 µm，均在容差；仅两种确定性扰动、未做各自dt/cache重复，不泛化鲁棒性。
13. 找到低压力四路Y候选，COM/形状/拓扑/固定ID都通过实际半步与细缓存对照。
14. Y仅移动0.012944032 µm，目标误差0.187055968 µm，未达到0.2 µm目标。
15. 存在：历史F10Q Y的COM通过，但形状/ID失败。新弱Y通过并不抹去历史失败。
16. 本例两段控制比原固定控制目标误差少约0.221742 nm，形态代价相近；不声称分段普遍占优。
17. IDEAL_SCHEDULED_SEGMENTS；相位公共时钟和理想同步，USB延迟/同步偏差未测，不虚设仪器参数。
18. 最长实际完成2 ms；没有5/10/30 ms、长时间陷阱或任意Z导航。
19. 全部科学数值、分析、测试及独立核验255.274639 s；其中远端四进程229.789284 s、实际积分224.291758 s（已包含在远端时间中，不重复相加）。CPU QP 155.818568 s、阻力装配31.842679 s均为积分内分项；峰值RSS 646.902 MiB，纯GPU核时未独立剖析。
20. 四组PNG/PDF在本工作区results/straight_tube_exp1/f10r/figures/；完整绝对路径见文末。
21. 已完成备份及最小插入更新，原文全部字节及标题结构保留；证据progress/UPDATE_AUDIT.json，备份名microbubble_aggregation.md.bak_F10R_20261014。
22. 优先验证相位移动X陷阱的三维力稳定性并有限扩大X距离；Y需避免把代理近零点结论外推实验，先标定/确认场再优化低风险Y。任意团形状控制证据不足，不优先扩大复杂形状任务。

## 补充要求的实际结论

单片场先于新Y轨迹完成：Y/X压力效率约1/202.177，Y零点距4.504167 mm、实际距4.5 mm，确认代理特有的Y低效。每片的复梯度及一次力完整保存在SINGLE_PZT.json和Y报告，全部Y标签HOMOGENEOUS_WATER_PROXY_SPECIFIC。

只改X+相位±1/±3/±6/±10°，局部一次声力敏感性3.73%/度，与6%/度预测不一致。−6/+6/+10°真正运行完整2 ms，均观察到容差内；+10°额外半步/细缓存通过，距边界仅约1.583 nm。−10°未运行轨迹，不能把场推算称实际越界。

等幅同轴对连续压力峰的位置曲线X约-0.842628 µm/度，Y强正侧支约-0.774094 µm/度；都是场计算，未模拟长时陷阱导航，也不是已证明的三维稳定陷阱。

## 失败、未运行与工程修订

25条本轮轨迹均完成各自登记窗口并通过基础质量，无积分拒绝；但这不代表25条导航成功。单片相消序列和新Y未达目标；弱单Y/Y对未进入完整窗口；原Y形态/ID失败保留；旧强控制0.884 ms活动集循环是历史失败，本轮不重启。

原Y压力峰最近中心选取混淆分支，原错误结果保留并明确作废，改为连续强分支审查。另发现输出循环变量覆盖RESULT.source文本（在积分完成后），输入、控制及轨迹均正确；派生报告改读不可变INPUT.source，并恢复原OFF零响应准入规则。原RESULT和远端SHA不改，SOURCE_METADATA_REVIEW.json保留修订前后。没有修改科学轨迹或阈值掩盖失败。

未运行项包括旋转/±6°的完整数值准入、−10/±1/±3°完整轨迹、长时移动陷阱、任意Z或复杂团形状、实验标定与实际硬件。不存在把NOT_RUN当PASS的情况。

## 冻结执行与图件

四个独立Supervisor标签mevo_f10r_pilot_v1、screens_v1、navigation_v1、robustness_v1，均已退出；各自固定源码commit与tar SHA在execution/和logs/*_EXECUTION.json，运行没有重用可变工作目录。四场缓存/APECSS输入及预登记SHA一并保存。

工作区：`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10r_20261014/computation_examples/mevo_acoustic_microbubbles`。图件目录：`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10r_20261014/computation_examples/mevo_acoustic_microbubbles/results/straight_tube_exp1/f10r/figures`。Fig_F10R_1_segmented_control、Fig_F10R_2_x_navigation_robustness、Fig_F10R_3_y_morphology_sensitivity、Fig_F10R_4_control_admission各有PNG/PDF。输入SHA、脚本身份和图件SHA在FIGURE_PROVENANCE.json，实际视觉审查在VISUAL_REVIEW.json。原始轨迹、所有逐步接触/真实互补性/事件、力采样与分析数组均可追溯。

图形渲染、网络传输与工程开发时间不计入科学600 s预算；没有把这些耗时称为GPU求解。所有物理仍未实验校准，数值可重复不等于实验准确或AI训练准入。
