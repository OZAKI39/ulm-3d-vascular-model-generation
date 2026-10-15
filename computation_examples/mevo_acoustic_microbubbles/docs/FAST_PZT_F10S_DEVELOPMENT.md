# F10S 开发交付（2026-10-15）

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED；training_ready=false。

完整基线为F10R交付提交 `9f1db411836e9a67dff1f77a344e5457674e6a81`，不是公开F10Q重建。分支 `dev/fast-pzt-f10s-x-pair-onehot-trap-scope-20261015`；工作区 `/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10s_20261015/computation_examples/mevo_acoustic_microbubbles`。清理审查独立提交f22c48de3，没有证据充分的废弃代码，因此不删除；F1–F10R科学源码、历史、失败、轨迹、APECSS/场表及复现入口不变。

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
21. 已计时必要科学计算、分析与独立核查合计 **132.966068 s**；远端campaign 113.657085 s（A46.399/B67.258），包含保留的部分归档分析8.452 s。详见TOTAL_COSTS，非纯GPU时间，开发/传输/绘图/解释器导入不混入。
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
