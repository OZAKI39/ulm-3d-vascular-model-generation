# F10P实际控制结果与23项回答

本轮完成独立原型、有限实验、分层验收和四组图件；没有完成强控制2 ms或XY导航。原强控制合法时间从B2的0.174 ms增加到硬球dt1的0.463 ms，但新真实间隙互补门槛失败，不能称完整强控制通过。固定2.5 MHz、32泡、R1.78 µm、四枚7 mm PZT、均匀水复声压缓存、F2/APECSS、原R3封闭X真实壁面和严格零流不变。

## 逐项回答

1. 开发前核验完整F10O的98cf17b8e和全部归档；本轮无新增安全删除对象。546f286f1为独立入口整理提交，M0、B2、所有历史源、数据、失败、场及表均保留。
2. 在当前法向下求受约束速度：min ½sᵀAs−bᵀs，Bs+g/R≥0，s=dt v/R。A为原成对阻力；反力满足非负、受力平衡和互补。只保证SPD局部QP的唯一速度，不宣称全多球最终位置唯一。
3. 是，硬球模块没有接触弹簧刚度，也没有事后PBD。B2刚度仅用于原参考对照，未修改。
4. 是，原F10G resistance原样复用，包括原润滑激活与间隙下限；并逐步做对称性和Cholesky检查。
5. 两泡/四泡六个独立参考全部通过：分离、挤压、释放、无接触、同时接触、网络平移；非法矩阵及迭代上限拒绝检查通过。
6. 是，150步及0.25 ms通过，越过原0.175 ms；这不是完整2 ms通过。
7. OFF和弱控制各完成dt1、dt0.5、细缓存的2 ms。强控制dt1在0.464 ms提案失败，合法0.463 ms；dt0.5诊断仅登记并完成0.5 ms，细缓存合法0.463 ms。强控制完整2 ms减半步长与细缓存没有运行。
8. 直接活动集各正式及诊断提案最大实际重叠约8.47e−13 nm，是浮点舍入量，不是实验精度。首轮OSQP拒绝提案最大约1.85e−9 nm；原17.8 nm及新1 nm要求都未放宽。
9. 强dt1拒绝提案局部平面互补等效间隙1.05e−15 nm、原始乘积约1.51e−33 N·m；真实球面等效间隙0.010056962 nm，超过预登记0.010000000 nm。物理受力残差约3.49e−23 N。不能用很小的局部QP残差掩盖真实几何反力不一致。
10. 原重叠、178 nm接触修正及墙预算均满足；全部墙修正0，未穿PDMS或封闭X端。但强控制完整轨迹未通过更严格、运行前登记的真实间隙互补检查，因此FULL_STEP_ACCEPTANCE失败。
11. 强控制共同0–0.463 ms同ID步长差21.64711846 nm，缓存差0.05214383 nm；前缀低于50 nm，不能提升为完整2 ms准入。弱控制2 ms同ID步长差1.0316757e-05 nm、缓存差0.00017910581 nm，通过。
12. 强前缀COM步长/缓存差4.2411653e-07/0.0017640246 nm；弱控制完整2 ms为7.9181899e-07/4.67898e-05 nm。OFF严格不动；弱控制质心位移模1.61910689 nm，大于所测差异，但不是微米级导航。
13. 弱控制Rg减少0.14685379 nm，X/Y从1到1.000049250923，全局长短比从2.2360679775到2.237041981373。Rg、Q、比值、匹配分布和可识别主轴分别通过登记阈值；初始XY简并，不能定义初始主轴转角。强前缀Rg步长差2.586263 nm、归一化Q差0.00099785525，仅为前缀证据。
14. OFF和弱控制初始32成员网络维持完整2 ms，无新持续合并/分裂。强控制只有短窗连通；dt0.5到0.5 ms发生两条瞬时边变化但仍一团，观察不足2 ms，不能验证持续事件。瞬时10 µs与持续2 ms验收独立。
15. 本轮硬球没有出现团级通过而单泡未通过的实际案例。代码另用明确标记的元数据测试验证两个门槛独立，绝不把该测试当真实科研样本。
16. 不完全一致。相同弱控制的B2 Rg减少0.184224 nm，硬球减少0.146854 nm，两者差0.037370 nm，大于双方测得步长/缓存误差；M0则增加11.1069 nm。标记CONTACT_MODEL_DEPENDENT，不能判断哪种模型更接近真实泡壳。
17. 原X=(0.2,0,0)、Y=(0,0.2,0) µm及0.075 µm容差完整保留。由于强控制基础质量失败，新XY候选优化没有运行；可选实现保持未准入，不制造优化轨迹。
18. 没有达到目标。强dt1合法终点位移(0.05196697,0.00944982,约0) µm，X/Y误差0.14833434/0.19750934 µm；弱2 ms误差0.19839314/0.19980847 µm。强前缀也没有完整位置准入，不能宣布TARGET_REACHED。
19. 强dt1实际平均整步2.389454 ms，原B2约4.989761 ms；两者有效时间不同，不称严格同窗加速比。硬球150步约3.469294 ms；弱dt1约2.467678 ms。成本包括场、力、约束、墙和传输。
20. 科学数值计算共44.855733 s，含首轮OSQP失败、当地小参考/冻结输入复核与四个服务器进程；远端44.801865 s。积分及逐步证据42.984579 s，CPU约束20.748126 s、阻力11.567359 s，峰值RSS536.996 MiB。未测独立GPU纯核时间，也不把部署、传输、图件时间冒充科研积分时间。
21. 四组英文PNG/PDF位于本工作区results/straight_tube_exp1/f10p/figures/，绘图源、真实数据、执行SHA、输入SHA、FIGURE_PROVENANCE及实际PNG/PDF预览检查一并保存。
22. 科研进度文档最小更新和备份由progress_document.py执行；before/after、diff、SHA及标题结构检查存于progress_document/。备份名microbubble_aggregation.md.bak_F10P_20261014。
23. 下一轮优先解决有限步接触反力与真实球面几何的一致性、核查强团内部迅速重排及步长差；再通过完整2 ms团级位置/形状准入后研究更远导航。不能优先用更大力、改容差、弹簧或更长轨迹掩盖该阻塞。

## 实际运行清单

| Run | Status | Legal time (ms) | Mean step (ms) | Integration (s) | Reason |
|---|---|---:|---:|---:|---|
| STRONG_PREFIX_DT05 | PASS | 0.5 | 3.1869 | 3.2206 | — |
| STRONG_PREFIX_FINE | FAIL | 0.463 | 3.2336 | 1.5195 | CONTACT_COMPLEMENTARITY |
| pilot_v1__STRONG_150 | FAIL | 0.003 | 13.2649 | 0.0596 | SOLVER_RESIDUAL |
| pilot_v2__STRONG_150 | PASS | 0.15 | 3.4693 | 0.5332 | — |
| pilot_v2__STRONG_SHORT | PASS | 0.25 | 3.4093 | 0.8660 | — |
| OFF_DT05 | PASS | 2 | 1.7115 | 6.9172 | — |
| OFF_DT1 | PASS | 2 | 1.6185 | 3.2727 | — |
| OFF_FINE | PASS | 2 | 1.8872 | 3.8186 | — |
| STRONG_DT1 | FAIL | 0.463 | 2.3895 | 1.1216 | CONTACT_COMPLEMENTARITY |
| WEAK_DT05 | PASS | 2 | 2.7983 | 11.2932 | — |
| WEAK_DT1 | PASS | 2 | 2.4677 | 4.9792 | — |
| WEAK_FINE | PASS | 2 | 2.6668 | 5.3833 | — |

三次失败是首轮OSQP上限和两次真实间隙互补拒绝，均保存LAST_STEP、INPUT、逐步HDF5及拒绝提案。77个远端文件8,216,449 bytes逐项SHA一致。未运行5/10/30 ms、AI训练、新CFD或新声场。

## 分层状态

ENGINEERING_PASS=PASS；OFF/WEAK的BASIC_NUMERICS_PASS、GROUP_POSITION_NUMERICALLY_ADMITTED、GROUP_SHAPE_NUMERICALLY_ADMITTED、GROUP_TOPOLOGY_NUMERICALLY_ADMITTED（仅初始网络生存范围）、PARTICLE_ID_NUMERICALLY_ADMITTED分别PASS。

STRONG完整2 ms BASIC=FAIL，后续完整窗口对照NOT_RUN且不准入；共同合法前缀单独报告，不能由COM小差异覆盖基础失败。TARGET_REACHED=false；SHAPE_TARGET_REACHED=false；EXPERIMENTALLY_CALIBRATED=NOT_VERIFIED；TRAINING_READY=false。

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / training_ready=false。
