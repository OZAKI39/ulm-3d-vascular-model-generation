# 直管损伤—断裂—脱离—输运验证报告

日期：2026-09-30。工程：`/home/lzy/projects/clot_clearance`，在既有原型上增加独立扩展文件，原 NOSB-PD 实现和原算例保留。新入口为 [OPEN_FRAGMENTATION.html](OPEN_FRAGMENTATION.html)，复现命令见 [README_FRAGMENTATION.md](README_FRAGMENTATION.md)。

**结论：已跑通由载荷与经验损伤驱动的实际断键、连接丧失、脱离及流体驱动输运流程。结果只属于 qualitative forced-failure verification；所有促成失效的参数未经实验标定。** 新动画没有施加开裂位移、手动切键、手动分离粒子或人工设置碎片编号，也没有放大位移。

## 1. 结果身份与精确指标

选定配置：[configs/streaming_fragmentation_demo.json](configs/streaming_fragmentation_demo.json)。实际数据：`runs/fragmentation_pilot_003`，逻辑入口：`results/streaming_fragmentation_demo`。完成 25 个宏步、N = 25,000 个代表循环，机械代理时间 0.052 s；计算墙钟时间 191.854210388 s。原始配置、代码哈希、粒子和键数见 [IDENTITY.json](results/streaming_fragmentation_demo/IDENTITY.json)。

| 指标 | 结果 |
|---|---:|
| 首次真实断键 | N = 2,000，机械时间 0.006 s |
| 首次脱离 / 首次脱离分量 | N = 11,000，机械时间 0.024 s |
| 首次清除平面穿越 | N = 15,000，机械时间 0.032 s |
| 初始 / 最终保留粒子 | 360 / 360 |
| 初始键数 / 最终活动键 | 12,592 / 1,538 |
| 最终断键数 / 比例 | 11,054 / 0.877858958068615（87.7858958068615%） |
| 最终附着体积占比 | 0.39444444444444426（39.4444444444444%，142 个粒子） |
| 最终脱离体积占比 | 0.6055555555555556（60.5555555555556%，218 个粒子） |
| 最终分量数 | 178 = 1 个附着分量 + 177 个脱离分量 |
| 最大碎片质心净位移 | 0.0030876717753563687 m = 3.087671775356369 mm |
| 最终清除体积占比 | 0.4333333333333333（43.3333333333333%，156 个不重复粒子） |
| 最终平均 / 最大粒子损伤 | 0.8651173715376823 / 1.0 |
| 全部粒子总体积 | 7.031250000000002 × 10⁻¹⁰ m³ |

177 个脱离分量包含 **172 个单颗粒、4 个双颗粒分量和 1 个 38 颗粒分量**。最大的脱离多颗粒分量 ID = 103，最终净移动 1.504275399482664 mm，体积 7.421875000000003 × 10⁻¹¹ m³，其内部活动键仍参与受力。最大 3.08767 mm 指标包括单颗粒分量，不能解释为全部碎片均移动了这一距离，也不是累计路程。

## 2. 计算与统计如何形成

每个宏步在实际动力学轨迹上测量键伸长的半峰峰幅值 Q，再更新 `DeltaD = DeltaN * C_damage * max(Q/Q_ref - 1, 0)^m`。活动键由 `1-D` 削弱；达到 `D_break` 后完整性设为零且永不恢复。所有 11,054 条断键均有含 Q、损伤前后值、增量和阈值的 [failure_ledger.jsonl](results/streaming_fragmentation_demo/failure_ledger.jsonl)。

在活动键图中，与至少一个固定基底粒子连通的材料为附着材料；其余为脱离材料。稳定碎片 ID 由连通分量的分裂谱系产生。每次宏更新均重新判断流体暴露面，输出质心、质心速度、包围盒和体积。独立并查集核验了连通性、附着可达性、全部关键统计、损伤公式和每条断键的阈值证据。

附着暴露粒子在当前坐标上接受原解析载荷；脱离粒子继续保留质量和活动键，同时接受 `a_h = (u_f-v)/tau_h`，固定基底不移动。清除平面位于 x = 1 mm；分量质心向下游穿越时仅累计尚未记账的粒子体积。清除材料不删除，后续分裂不重复记账。因此“清除比例”表示诊断平面穿越量，不是化学溶解量或临床血管再通率。

额外力学归因对照从 N = 11,000 的实际脱离状态出发，保持两支的初始坐标、速度和现有键完全相同，冻结损伤仅用于这项 2 ms 对照。开启流体力时，首批 4 个脱离粒子的质心下游位移为 +0.189024143334417 mm；关闭时为 −0.0872842856983594 mm，两支质心位移向量之差为 0.276427440843857 mm。没有重置速度或施加分离命令。证据见 [TRANSPORT_FORCE_ATTRIBUTION.json](results/streaming_fragmentation_demo/TRANSPORT_FORCE_ATTRIBUTION.json)。正式主算例的损伤更新没有冻结。

## 3. 所有为失效演示选择的参数及其他变化

以下列出相对原温和配置的全部变化；机器可读差异为 [CONFIG_DIFF.json](provenance/fragmentation_extension/CONFIG_DIFF.json)。配置变化均在新配置文件内，试跑之间没有调整力学源代码常数。

| 配置项 | 原轻度损伤 | 新演示 | 作用与解释 |
|---|---|---|---|
| `damage.Q_ref` | 0.014 | **0.001** | 显著降低循环损伤启动尺度，促成失效 |
| `damage.D_break` | 无独立提前阈值 | **0.25** | 提前断键；新代码默认仍为 1 |
| `damage.C_damage` | 2×10⁻⁵ | **5×10⁻⁷** | 与降低后的 Q_ref 联合选择；C 本身降低，不能称为单独增大疲劳速率 |
| `simulation.number_of_macro_steps` | 8 | **25** | 总代表循环从 8,000 增至 25,000 |
| `clot.cells` | [12,8,4] | **[10,6,6]** | 增厚并缩短、收窄局部斑块，形成可失效的自由材料 |
| `clot.origin_m` | [−0.00075,−0.0005,−0.0011] | **[−0.000625,−0.000375,−0.0011]** | 配合新几何保持居中及基底位置 |
| `streaming.bubble_center_m` | [0.00015,0,−0.00045] | **[0.0001,0,−0.0002]** | 使局部载荷对应增厚后的上表面附近 |
| `streaming.traction_scale_Pa` | 140 | **100** | 重新选择载荷尺度；未为视觉效果增大应力 |
| `simulation.detached_damping_s_inv` | 未区分 | **0** | 自由材料不再受附着固体的数值阻尼；仍有流体松弛 |
| `transport.enabled` | 无 | **true** | 启用脱离后的替代输运耦合 |
| `transport.tau_h_s` | 无 | **0.0002 s** | 0.2 ms 强松弛使自由材料跟随规定流速；未经标定 |
| `transport.background_velocity_multiplier` | 无 | **1** | 背景流速不额外放大 |
| `transport.x_clearance_m` | 无 | **0.001 m** | 清除诊断平面 |
| `surface.weighted_bulk_fraction` | 原固定表面 | **0.8** | 邻域体积亏损触发动态表面暴露 |
| `surface.area_factor` | 无 | **1** | 每暴露粒子使用 h² 的近似受载面积 |
| `surface.normal_moment_tolerance` | 无 | **10⁻¹²** | 邻域方向矩低于阈值时不定义载荷法向 |
| `safety.maximum_displacement_m` | 0.0004 | **0.004** | 允许真实下游移动；只是中止阈值，不施加或缩放位移 |
| `safety.rank_deficient_policy` | 原严格中止 | **intrinsic_correspondence** | 显式低维支撑延续，限制见下一节 |
| `safety.support_rank_relative_tolerance` | 无 | **10⁻¹²** | 识别实际剩余支撑维数 |
| `safety.maximum_radius_ratio` | 无 | **1** | 超出管半径则中止；不等于壁面接触求解 |
| `name` / `verification_label` | 原内部算例名 / 无 | `streaming_fragmentation_demo` / `qualitative forced-failure verification; UNCALIBRATED` | 区分结果身份 |

保留而仍必须披露的经验数值为：`m_damage=2`、`DeltaN=1000`、流涡速度尺度 0.01 m/s、长度尺度 0.35 mm、慢载荷频率 500 Hz、振幅 0.9、均值系数 1、法向附加载荷 0、附着阻尼 3500/s。载荷与速度场分别指定，不是联立求出的声流应力与速度。`s1=1.1`、`s2=1.3` 保留在配置中，但此循环损伤分支不使用它们。

时间步 2 μs、粒子间距 0.125 mm、horizon/spacing = 3.01、G = 1000 Pa、K = 50000 Pa、稳定化参数 3、密度 1056 kg/m³ 均未改变。直管半径 1.25 mm、参考长度 6 mm、Q = 18 mL/min、黏度 0.00345312 Pa·s 均保留。1 MHz 是载波元数据，2 μm 是规定气泡半径，未求解气泡形变。最终最大粒子位移 3.269093151749531 mm，低于 4 mm 中止阈值。

## 4. 试跑、原算例回归与检验

| 数据目录 | 结局 |
|---|---|
| `runs/fragmentation_pilot_001` | 完成 35,000 循环；断键比例 0.0054002541296060995，但没有脱离，未满足新演示验收 |
| `runs/fragmentation_pilot_002` | C = 10⁻⁶，已有断键/脱离；第 25 宏步中触发 4 mm 位移保护，在 N = 24,000、机械时间 0.050602 s 停止。保留 FAILED.json 与部分结果，不列为成功完成 |
| `runs/fragmentation_pilot_003` | C 调为 5×10⁻⁷，完整完成；即本报告选定结果 |
| `verification/fragmentation_extension/mild_reproduction` | 原 runner 与原配置副本重新运行；全部已保存数值数组与旧结果逐位相等 |

原轻度损伤结果仍称为 `straight_pipe_damage_demo`，保留 384 粒子、12,804 条键、N = 8,000、零断键、零脱离。新结果不能覆盖它。

CTest 集成运行通过 26 项 pytest 检验，包括原 18 项及新增 8 项。新增项检查全维旧内核一致性、低维能量梯度与客观性、阈值不可逆性、稳定分量编号与清除去重、暴露检测、固定基底与流体松弛、速度场散度，以及存活的低维键受力。日志：[tests_final.log](logs/fragmentation/tests_final.log)。

全时序独立审核：[INDEPENDENT_AUDIT.json](results/streaming_fragmentation_demo/INDEPENDENT_AUDIT.json)。原文件哈希、VTK 与 NPZ 一致性、源代码身份、视频完整解码和本地链接检查集中在 [FINAL_QC.json](provenance/fragmentation_extension/FINAL_QC.json)。浏览器交互本身未做自动化端到端测试。

原工程的 352 个既有文件全部保持原哈希。另一个 2,640 项的保护清单中，2,639 项没有变化；共享 `BRAVA_CONTEXT_HANDOFF.md` 在本次工作期间由其他工作追加了网格研究进度。其原有 16,513 字节前缀仍与原 SHA256 完全一致，新增 2,692 字节原样保留，本扩展没有写入或回滚该文档。最终核验状态明确记为 `PASS_WITH_SHARED_CONTEXT_APPEND`，并保存观察副本。

## 5. 可视化交付

- [损伤与碎片双面板 MP4](visualization/streaming_fragmentation_demo/fragmentation_damage.mp4)、[对应 GIF](visualization/streaming_fragmentation_demo/fragmentation_damage.gif)。
- [碎片与损伤双面板 MP4](visualization/streaming_fragmentation_demo/fragmentation_fragments.mp4)、[对应 GIF](visualization/streaming_fragmentation_demo/fragmentation_fragments.gif)。
- [四阶段 PNG](visualization/streaming_fragmentation_demo/figures/four_panel_summary.png)、[PDF](visualization/streaming_fragmentation_demo/figures/four_panel_summary.pdf)：N = 0、2,000、11,000、25,000。
- [七项时间历程 PNG](visualization/streaming_fragmentation_demo/figures/time_histories.png)、[PDF](visualization/streaming_fragmentation_demo/figures/time_histories.pdf)。

两部 MP4 均为 1600×900、12 fps、156 帧、13 s；26 个真实保存状态各重复 6 帧，仅调节播放速度。镜头、坐标范围与单位比例固定；所有位置使用真实位移因子 1。显示每个分量的活动键生成树及每宏步最多 350 条新断键；完整键数据始终保留在 VTK 中。碎片按稳定 ID 着色，显式文字标注最大的脱离分量。灰方块为固定基底，蓝点为放大显示的规定气泡位置，绿色虚线为清除平面。已查看初始、首次断键、首次脱离和末态，末态可见材料向下游分离移动。

## 6. 适用范围与显著限制

最终粒子的支撑秩 0/1/2/3 分别为 172/18/20/150。全维仍使用既有 NOSB-PD；低维使用显式内禀支撑近似，零秩粒子无内部键力但保留流体力与运动积分。不能把这一延续方案当作通过实验验证的三维碎片本构。详细公式、表面规则、稳定编号及清除记账见 [FRAGMENTATION_MODEL.md](FRAGMENTATION_MODEL.md)。

本演示只验证“声流式载荷 → PD 损伤 → 真实断键 → 脱离 → 碎片输运”的数值流程。它不证明正确的超声频率响应、微泡形状振荡、临床清除率、实验标定疲劳寿命或定量正确的碎片尺寸。没有对该新失效算例声称网格/时间步收敛；没有双向流固反馈、碎片碰撞、化学溶解或真实壁面接触。轴向背景场作解析延拓，清除面不是出域删除面。

N/1 MHz = 0.025 s 只作周期账本；实际积分机械代理时间为 0.052 s，二者都不是已校准的治疗时间。供应论文不是本次超声损伤参数的实验依据。下一科学步骤仍是替换为解析充分或实验验证的微泡流场及真实流体—固体载荷耦合，再验证材料、损伤、碎片尺度和清除率。
