"""Write evidence-based reports and append-only Chinese progress update."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib, json, difflib, re, shutil
import numpy as np
from .common import *

def main():
 assert read(BASE/'analysis/INDEPENDENT_REVIEW.json')['status']=='PASS'
 cfg=read(CONFIG);runs=read(BASE/'analysis/RUN_ANALYSIS.json');comparisons=read(BASE/'analysis/COMPARISONS.json');gates=read(BASE/'analysis/FINAL_GATE.json')
 costs=read(BASE/'analysis/COSTS.json');partial=read(BASE/'history/partial_analysis_before_transfer_completed/COSTS.json')
 items=dict(A_campaign_s=read(BASE/'a_navigation_v1/CAMPAIGN.json')['wall_s'],B_campaign_s=read(BASE/'b_independent_v1/CAMPAIGN.json')['wall_s'],
      diagnostic_C_D_s=read(BASE/'diagnostics/EXECUTION.json')['wall_s'],PDMS_preflight_s=read(BASE/'diagnostics/PDMS_PREFLIGHT.json')['wall_s'],
      partial_analysis_retained_s=partial['local_analysis_s'],complete_analysis_s=costs['local_analysis_s'],independent_review_s=read(BASE/'analysis/INDEPENDENT_REVIEW.json')['wall_s'])
 total=sum(items.values());write(BASE/'analysis/TOTAL_COSTS.json',dict(components=items,total_measured_scientific_s=total,budget_s=600,budget_pass=total<600,
     timing_scope='Measured scientific function/campaign wall clocks; campaign includes context initialization, checks, integration and output. Interpreter imports, deployment, SSH transfer, development and figure rendering are outside these timers. Not pure GPU time.',
     remote_campaign_s=costs['remote_process_s'],integration_s=costs['integration_s'],CPU_solver_s=costs['CPU_solver_s'],CPU_resistance_s=costs['CPU_resistance_s'],
     plot_render_last_pass_s=read(BASE/'figures/FIGURE_PROVENANCE.json')['render_wall_s'],max_particle_physical_time_s=.003))
 prefix='EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED；training_ready=false。\n\n'
 A=runs['X_PAIR_3MS_DT1'];minus=runs['X_MINUS_ONLY_DT1'];plus=runs['X_PLUS_ONLY_DT1']
 docA='''# F10S：X同轴对3 ms验证

从 F10R `navigation_v1/runs/X_PAIR_DT1/INPUT.json` 读取完整精度值：X−和X+均为1036.2016581041792 Pa-equivalent，相位分别0和3.6651914291880923 rad（210°），Y±关闭，2.5 MHz。原32个固定ID、R=1.78 µm，原3/1.5 µm复场缓存、APECSS、F2、润滑阻力、五次Sequential QP与R3壁面均未改。3 ms从原初态新鲜计算，没有拼接。

目标固定绝对坐标(0.2,0,0) µm；原容差0.075 µm，新增严格终点评价0.01 µm，均已预登记。

| 分辨率 | X终点/µm | 目标误差/µm | 基础质量 |
|---|---:|---:|---|
'''
 for key in ['X_PAIR_3MS_DT1','X_PAIR_3MS_DT05','X_PAIR_3MS_FINE']:
  v=runs[key];docA+=f"| {key} | {v['endpoint_m'][0]*1e6:.12f} | {v['target_error_initial_final_m'][-1]*1e6:.12f} | {v['status']} |\n"
 docA+='''
主结果0.197999446649 µm、终点误差0.002000553351 µm；比旧2 ms误差0.068156397464 µm下降约97.06%。原容差与严格评价均PASS；最接近目标在3 ms，未越过目标，保存的每个时间步目标误差没有反向增加。Y/Z仅机器精度舍入。

全窗口COM步长/缓存最大差0.000117420/0.005690938 nm；固定ID差0.001179081/0.005695114 nm。原COM、去质心形状（Rg、Q、长短比、匹配分布、可识别方向）、网络及ID分别PASS。32成员网络完整保持3 ms；原2 ms持续性定义不变，3 ms扩展仅修正覆盖判定。初始主轴在XY平面退化，不虚构方向控制。

Rg减少0.295179629 nm，Q相对范数变化约0.00999925%，长短比仍约2.23606798。变化超过步长/缓存差，但只是接触初态轻微压紧，不是指定形态控制或物理融合。泡处压力789.080–1207.482 Pa；一次力每泡平均模2.21418 pN、最大2.63964 pN；二次合力平均模170.847 pN、最大238.553 pN。力统计来自实际0–2.9 ms采样；精确终点APECSS/F2另存 `X_PAIR_3MS_DT1_FROZEN_FORCE.json`。

最大接触修正0.0381507 nm，墙面修正严格0，最大实际互补等效间隙约3.83e−10 nm；原17.8 nm重叠、178 nm修正、1 nm硬球几何、0.01 nm真实互补及残差阈值不变。质量通过不等于实验物理准确。

数据：`results/straight_tube_exp1/f10s/a_navigation_v1/`；全窗口指标、比较与分层判定：`analysis/`。
'''
 docB='''# F10S：独立X单源验证

各自重新加载同一原始32泡状态，400 Pa-equivalent、零相位、2.5 MHz，只启用一枚X端PZT，其余三枚关闭。不是旧先后切换序列的截取。每枚先0.5 ms独立筛选，实测方向和基础质量通过后，重新从原初态运行2 ms及半步/细缓存对照。

| 单源 | 2 ms ΔX/µm | 平均速度/(µm/s) | 最大COM数值差/nm |
|---|---:|---:|---:|
| X− | +0.015835804138 | +7.917902069 | 0.000489614 |
| X+ | −0.015835804138 | −7.917902069 | 0.000489614 |

运动信号约为最大实测数值差的32344倍；只读原OFF对照逐泡严格不动。全模型一次力、二次力、F2、润滑、真实接触反力和壁面参与了两条独立方向轨迹，没有人为翻转力。两组基础、COM、形状、网络、固定ID都PASS；本任务没有登记单源位置目标，TARGET_REACHED=NOT_RUN。

实际泡处压力742.614–742.998 Pa，一次Fx每泡时空平均分别±0.2656631 pN，最大一次力模0.2657900 pN。二次合力平均模92.3964 pN、最大99.4998 pN；接触反力平均模87.7772 pN。精确末态F2增量平均45.1251 pN、最大47.1663 pN。力模和带符号X分量分开记录，二次内力总和为零不能解释为团级外推力。

两组Rg均减少约0.295179607 nm（约0.005%），Q相对变化约0.010%，长短比基本不变；固定ID分布及原32泡网络保持。这个微小收缩超过所测数值差异，但不支持显著或任意形态控制。最大修正约0.0170414 nm、墙面修正0，原门槛不变。

数据：`b_independent_v1/`；短筛选和完整2 ms分别保留，所有时间步、反力、压力、末态F2均可审查。该定性方向证据仅适用于当前均匀水代理和初态，不是外部实验复现。
'''
 docC='''# F10S：相位压力峰与力稳定性

预登记X等幅对400/1036.2016581041792 Pa-equivalent，相位差−6/0/+6°，原2.5 MHz APECSS，原32泡几何；6个候选均在表有效范围内。直接评估原48阶圆盘积分及同一积分的解析梯度，独立96阶核查和压力有限差分通过；没有新增k-Wave或场外缓存外推。

压力最大值不严格等于一次力零点：0°由对称性重合，±6°约相差±18.94 nm。全部六个候选的一次力Jacobian三个实特征值均为负，因此在当前模型中存在局部三维恢复性单泡平衡；步长50/25 nm的Jacobian核查通过。1036.202 Pa-equivalent、0°的特征值约[−9.13777e−8,−1.37855e−9,−1.37855e−9] N/m。该局部稳定性不代表任意Z向独立导航；本轮平面布局对称初态的Z位移为零。

固定32泡内部几何的整体平移诊断使用 `Tᵀ ζ A T U=Tᵀ F`，A来自原成对润滑阻力。对于刚性平移，成对相对阻力和内部二次力的总和消去；这些力仍保留用于压紧风险评估。分别寻找团心平衡并用实际APECSS一次力求速度Jacobian，标记 **FROZEN_SHAPE_COM_DIAGNOSTIC**。

400 Pa-equivalent的X恢复时间约2.47253 s，1036.202约0.368446 s。以F10R实测局部峰斜率−0.842628 µm/度做1°阶跃，局部线性模型预测2 ms跟随位移仅−0.0006813/−0.0045616 µm；它是诊断预测，不是逐泡轨迹，也不是已经验证的相位斜坡跟随。时间尺度分别比2 ms大约1236/184倍，因此按预登记停止条件不跑长时间陷阱。

高压峰增加压紧风险：0°、400与1036.202源幅下，原几何处的二次合力平均模约369.239/2477.857 pN、最大397.208/2665.549 pN；F2增量平均180.240/1209.538 pN，来自实际响应而非P²缩放。自由速度诊断估计最强接近在1 µs内约0.0446/0.2994 nm；压紧virial为负。不能把这一局部自由估计当成接触预算通过，也不能宣称已观察到严重数值失败：未运行完整高峰接触轨迹，完整陷阱BASIC/形状/ID准入均NOT_RUN。

所有原坐标、复压、梯度、一次/二次/F2力、边网络、平衡位置、Jacobian、有限差分与简化预测保存在 `diagnostics/TRAPS.json`。未把声压峰的移动当作已实现团导航。下一步相位跟随需要低风险压力、内部形态松弛与闭环滞后研究；目前更适合短程直接推移后OFF停止。
'''
 pd=read(BASE/'diagnostics/PDMS_PREFLIGHT.json')
 docD='''# F10S：Y近场零点与PDMS计算适用性

所有Y结论为 **HOMOGENEOUS_WATER_PROXY_SPECIFIC**，几何偏移仅为探索性分析，不能当作用户装配误差。

固定1500 Pa-equivalent、2.5 MHz，Y−源距中心4.5 mm，理论轴上第二零点在4.504166667 mm。源距偏移−25/−5/0/+5/+25 µm，团心压力分别96.808/30.335/13.778/2.753/68.614 Pa；基线Y/X幅值比0.00494616（约1/202.18）。同样的±5 µm合法团心位移也跨越明显压力变化；一次Fy在零点两侧反号。单源公共相位0/90/180°不改变压力模和周期平均力，已检查。

2.475/2.5/2.525 MHz仅场级变化得到174.756/13.778/202.258 Pa。非2.5 MHz没有调用旧APECSS力或运动，保存NOT_RUN_FREQUENCY_OUTSIDE_APECSS。因此Y低效高度依赖波长和此几何的均匀水近场零点，不能外推真实PDMS。

参考材料：Guangyao Xu et al., *Acoustic Characterization of Polydimethylsiloxane for Microscale Acoustofluidics*, Physical Review Applied 13, 054069 (2020), DOI [10.1103/PhysRevApplied.13.054069](https://doi.org/10.1103/PhysRevApplied.13.054069)。Tables II/III的Sylgard184 10:1、85°C固化、25°C测量参考组，测量频带3–7 MHz；原F10M保存论文来源、摘要参数、PDF SHA和表格目视审查。F10S另核查作者/DTU公开记录。论文含Kramers–Kronig拟合，允许作为新假设研究外推，但本轮没有把2.5 MHz外推材料静默准入。

预检保留32.5×9×9 mm块、1.6 mm通道、两个1.5 mm顶部水孔、四7 mm圆片、X封闭几何。域含条件性1 mm外部水包层，BLI额外36网格支持余量，域外20网格PML；参考低频声速仅作条件性资源规划下界，不作为已验证2.5 MHz材料输入。水波长0.6 mm；CFL≤0.2；PDMS等效流体不包含完整剪切波，真实端封材料/厚度与外部耦合仍未校准。

| PPW | 网格间距/µm | 含PML格数 | 估计峰值/GiB | 单源估计/s |
|---|---:|---|---:|---:|
'''
 for r in pd['candidates']:docD+=f"| {r['PPW']} | {r['spacing_m']*1e6:.3f} | {'×'.join(map(str,r['expanded_shape']))} | {r['estimated_peak_GiB']:.3f} | {r['estimated_single_source_wall_s']:.1f} |\n"
 docD+='''
估算按F10L实测三源373.706 s、15100 MiB、378×192×192、3504步，按体素/时间步放缩，显存×1.2、耗时×1.5异质余量。是估计而非本轮实测；此源支持与均匀网格方案并非所有可能方法的最小成本。两候选都超20 GiB安全预算、24 GiB硬件及本轮120 s诊断预算。资源和材料准入未满足，**PDMS_Y_FIELD_NOT_ADMITTED / PDMS_FIELD_NOT_RUN**；没有求解多材料压力场，不伪造对照图。
'''
 for name,body in [('F10S_X_PAIR_3MS.md',docA),('F10S_INDEPENDENT_X_SOURCES.md',docB),('F10S_PHASE_TRAP_FEASIBILITY.md',docC),('F10S_Y_PDMS_SCOPE.md',docD)]:
  (ROOT/'docs'/name).write_text(body+'\n'+prefix)
 summary=f'''# F10S 开发交付（2026-10-15）

{prefix}完整基线为F10R交付提交 `9f1db411836e9a67dff1f77a344e5457674e6a81`，不是公开F10Q重建。分支 `dev/fast-pzt-f10s-x-pair-onehot-trap-scope-20261015`；工作区 `{ROOT}`。清理审查独立提交f22c48de3，没有证据充分的废弃代码，因此不删除；F1–F10R科学源码、历史、失败、轨迹、APECSS/场表及复现入口不变。

1. 源码基准：上述完整F10R HEAD；289个历史归档哈希在开发前核对。
2. 安全清理：仅审查与活跃索引；没有删除文件，历史全部保留。
3. X对确实从原始32泡状态新鲜积分到3 ms；1/0.5 µs与1.5 µm细缓存共三条。
4. 主结果ΔX=0.197999446649 µm，终点目标误差0.002000553351 µm。
5. 比旧2 ms误差0.068156397464 µm更近，误差降低约97.06%。
6. 原0.075 µm目标容差PASS。
7. 新0.01 µm严格终点评价PASS；目标未移动。
8. 完整3 ms的步长与缓存比较PASS，最大COM差0.005690938 nm，固定ID差0.005695114 nm。
9. 独立X−单片沿+X：2 ms为+0.015835804138 µm。
10. 独立X+单片沿−X：2 ms为−0.015835804138 µm。
11. 两者最大COM数值差约0.000489614 nm，运动信号约为其32344倍。OFF历史同模型逐泡严格静止。
12. 单源Rg仅减少0.295179607 nm，长短比不显著改变；微小变化可分辨但不属于指定形态控制。
13. 压力峰与力平衡在0°对称重合，±6°约差18.94 nm。
14. 六个候选存在当前代理模型中的局部三维恢复性单泡平衡；没有完整团的长期陷阱验证。
15. 固定内部几何诊断X跟随时间约2.47253/0.368446 s。忽略内部形态演变与相位移动过程中非线性/接触重排，不是真实刚体模型或训练轨迹。
16. 高峰二次力增强，最大每泡合力约2.666 nN且有压紧趋势；完整高压接触失败与否NOT_RUN，不虚称安全或已发生失败。
17. Y零点对几何和波长敏感；所有Y结论HOMOGENEOUS_WATER_PROXY_SPECIFIC。
18. 完成2.5 MHz PDMS单源条件性材料/域/网格/PML/显存/耗时预检。
19. 未运行PDMS源。参考材料有效带、真实密封/耦合未确定，估计60.0–109.5 GiB、32.6–74.7 min；PDMS_Y_FIELD_NOT_ADMITTED。
20. 本轮最长真实逐泡时间3 ms；11条新轨迹完成各自窗口，无长时间陷阱或硬件任务。
21. 已计时必要科学计算、分析与独立核查合计 **{total:.6f} s**；远端campaign {costs['remote_process_s']:.6f} s（A46.399/B67.258），包含保留的部分归档分析8.452 s。详见TOTAL_COSTS，非纯GPU时间，开发/传输/绘图/解释器导入不混入。
22. 四组英文PNG/PDF：`results/straight_tube_exp1/f10s/figures/`；输入SHA、绘图代码身份、轨迹commit与视觉审查齐全。
23. 科研进度文档先备份 `.bak_F10S_20261015`，仅插入F10S段落；before/after/diff/SHA及标题结构审查在progress_update/。
24. 下一轮优先有限扩大X指定目标并研究短推移后OFF停止/反馈，保持原数值预算；相位跟随时间太长且高峰压紧较强，暂不跑长陷阱。Y先做独立声场/实验参数校核；复杂形态控制和仪器通信继续后置。

| 独立状态 | X对3 ms | X−单源2 ms | X+单源2 ms |
|---|---|---|---|
| BASIC_NUMERICS_PASS | PASS | PASS | PASS |
| GROUP_POSITION_NUMERICALLY_ADMITTED | PASS | PASS | PASS |
| GROUP_SHAPE_NUMERICALLY_ADMITTED | PASS | PASS | PASS |
| GROUP_TOPOLOGY_NUMERICALLY_ADMITTED | PASS | PASS | PASS |
| PARTICLE_ID_NUMERICALLY_ADMITTED | PASS | PASS | PASS |
| TARGET_REACHED | PASS | NOT_RUN（未登记目标） | NOT_RUN（未登记目标） |

TRAP_FORCE_STABILITY_CHECKED=PASS（局部力检查）；完整移动陷阱轨迹NOT_RUN。PDMS_FIELD_ADMITTED=FAIL（预检拒绝），PDMS_FIELD=NOT_RUN。EXPERIMENTALLY_CALIBRATED=NOT_VERIFIED，TRAINING_READY=false。网络保持不等于新聚合或融合；初始32泡是人工合法接触演示初态。

工程核查：原17.8/178 nm预算、1 nm几何、0.01 nm实际互补及全部残差不变，墙修正严格0。2 ms旧判定回归与3 ms部分窗口负例通过；GPU压力/力与CPU原参考、APECSS非法范围拒绝和壁面拒绝由A的CHECKS.json实测。独立复核覆盖全部11轨迹、缓存范围、原初态、固定ID、真实球面距离、梯度与积分阶数。没有新的物理求解器。

保留问题：第一次单元测试使用了错误的历史比较字典键，修正测试夹具后通过，没有改科学结果；首次分析在回传尚未完成时正确把缺少X+重复标NOT_RUN，保留history目录，全部75个远端文件SHA验证后重做完整分析。B的部署辅助metadata中route=A来自未使用的旧默认选项；实际F10S_MODE=B，campaign.mode=B，执行代码按mode读取180 s B预算，不按该辅助字段决定计算。所有原始执行元数据保留。

复现入口（先读取预登记，远端标签和结果目录必须全新）：`python -m scripts.straight_tube_exp1.f10s.prepare`（只在新空结果目录）；`remote --mode A/B --label ...`；`archive`；`diagnostics`、`pdms_preflight`、`analyze`、`review`、`figures`。原输入为继承F10N/Q/R，不能无提示覆盖历史或重新prepare当前结果。主要科学任务A执行14731fd56、B/诊断a5f9e3701，最终分析3d21896da、核查d7bd7079f，图0d859d971；完整SHA与归档源码见execution/。局部诊断只用真实响应，未重算k-Wave；场缓存被逐泡试验复用。

详见F10S_X_PAIR_3MS、F10S_INDEPENDENT_X_SOURCES、F10S_PHASE_TRAP_FEASIBILITY、F10S_Y_PDMS_SCOPE及F10S_ACTIVE_CODE_MAP。
'''
 (ROOT/'docs/FAST_PZT_F10S_DEVELOPMENT.md').write_text(summary)
 # Preserve original bytes and heading structure; insert exactly one new passage.
 progress=Path('/home/lzy/projects/microbubble_aggregation.md');backup=Path(str(progress)+'.bak_F10S_20261015')
 addition=f'''\n**F10S 用真实3 ms验证了X同轴对的短距离导航。**在完整F10R提交9f1db4118上独立开发，有界清理无安全删除对象，所有旧物理和历史保留。严格读取原完整精度1036.2016581041792 Pa-equivalent、X+相对210°控制，从原32固定ID接触演示初态重新积分，不拼接2 ms轨迹。原绝对目标(0.2,0,0) µm及0.075 µm容差不变，运行前另登记0.01 µm严格终点评价。主结果ΔX=0.197999447 µm，终点误差0.002000553 µm，比旧2 ms误差0.068156397 µm更小，两项目标评价均通过；最接近时刻3 ms，未越过目标、误差无反向增加。半步及细缓存均通过，COM最大差0.005691 nm、固定ID0.005695 nm；基础、团位置、形状、32成员网络和ID分别通过完整窗口。Rg仅减少0.295180 nm、长短比不显著变，不能称指定形态控制或物理融合。原17.8 nm重叠、178 nm修正、真实互补及残差预算不变，墙面修正严格0。\n\n**两枚单源确实分别向相反方向推动团。**400 Pa-equivalent、2.5 MHz，每枚都先从同一原初态独立0.5 ms筛选，再新鲜2 ms及半步/细场：X−向+X移动0.015835804 µm，X+向−X移动同量，平均速度±7.917902 µm/s；不是先后切换片段。运动比实测COM数值差约大32344倍，原OFF逐泡严格静止。实际泡处压力742.614–742.998 Pa，一次Fx平均±0.265663 pN，原二次/F2、润滑和接触反力均保留。三组完整数值层次通过，Rg同样只微缩约0.29518 nm。\n\n**压力峰附近有恢复力，但相位陷阱不是毫秒级导航捷径。**少量预登记相位及实际APECSS一次力显示局部三维恢复性，0°压力峰与零力位置对称重合，±6°约相差18.94 nm。固定原32泡内部几何、整体平移并使用原阻力的FROZEN_SHAPE_COM_DIAGNOSTIC得到400/1036.202源幅下X跟随时间约2.473/0.368 s；1°峰阶跃的2 ms预测只跟随约0.000681/0.004562 µm，远慢于峰的0.843 µm变化，不是真实逐泡轨迹。高峰二次合力平均升至约2.478 nN、最大2.666 nN，有压紧风险；没有运行完整高峰接触，不能称已失败或已安全。按预登记停止长时间陷阱分支。\n\n**Y场低效仍是代理特例，复杂PDMS求解未获准。**源距−25至+25 µm的少量变化使团心压力从原13.778 Pa变到2.753–96.808 Pa，一次Fy在近场零点两侧反号；±1%频率只计算场，不外推2.5 MHz响应。Y全部标HOMOGENEOUS_WATER_PROXY_SPECIFIC。2.5 MHz真实尺寸域的条件性预检保留水腔、双顶孔和封闭X端，6/8 PPW、20网格PML、源支持余量估计约60.0/109.5 GiB、单源32.6/74.7 min；参考Sylgard材料测量带为3–7 MHz，密封和外耦合仍未知。因此PDMS_Y_FIELD_NOT_ADMITTED、PDMS_FIELD_NOT_RUN，没有伪造多介质图。\n\n本轮11条新轨迹完成各自窗口，最长3 ms；已计时科学计算、两次分析及独立复核共约{total:.3f} s，其中远端campaign约113.657 s，绘图/传输/开发另计。75个远端文件逐项SHA回传，四组英文PNG/PDF、原始轨迹、实际压力力场、输入/源码身份及视觉检查保存在独立F10S工作区results/straight_tube_exp1/f10s/，详见docs/FAST_PZT_F10S_DEVELOPMENT.md。下一步优先扩大有限X目标并研究短推移后OFF/反馈；相位跟随和Y材料声场先做有界校核，复杂形态与真实硬件继续后置。持续 **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED，training_ready=false**。\n\n'''
 before=progress.read_bytes();assert b'**F10S ' not in before,'Already updated; do not overwrite backup'
 assert not backup.exists(),'Existing backup must be reviewed before repeat update'
 shutil.copy2(progress,backup);marker='---\n\n### 尾声：故事里还没讲完的部分'.encode();assert before.count(marker)==1
 after=before.replace(marker,addition.encode()+marker,1);assert after.replace(addition.encode(),b'',1)==before
 out=BASE/'progress_update';out.mkdir(exist_ok=True)
 (out/'before.md').write_bytes(before);(out/'after.md').write_bytes(after)
 (out/'change.diff').write_text(''.join(difflib.unified_diff(before.decode().splitlines(True),after.decode().splitlines(True),fromfile=str(backup),tofile=str(progress))))
 headings=lambda b:re.findall(rb'^#{1,6} .+$',b,re.M)
 assert headings(before)==headings(after)
 progress.write_bytes(after)
 write(out/'AUDIT.json',dict(status='PASS',before_sha256=hashlib.sha256(before).hexdigest(),after_sha256=sha(progress),backup_sha256=sha(backup),backup=str(backup),heading_structure_identical=True,heading_count=len(headings(before)),historical_bytes_preserved_exactly=True,only_F10S_inserted=True))
 print('Reports and byte-preserving progress update complete',total)
if __name__=='__main__':main()
