# Particle-8.1：完整血管微泡轨迹可视化层——中文审核报告

Frozen flow full-vessel ULM-style microbubble trajectory simulation

已实际生成 **1249 条真实 INLET→OUTLET 完成轨迹**，保存全量接受子步、原始属性和失败记录，输出十张静态图、四组 MP4/GIF、逐帧账本与 CSV/JSON。达到 frozen-flow full-vessel ULM visualization 的可展示标准；不是 RBC-coupled physiological completion。
**覆盖范围：1/3 个出口分支有完成轨迹。完整 Frozen 网格均已显示，未观测分支不算轨迹重建完成。**

## 九项审核问题

1. **是否真正入口到出口？** 是。每条完成轨迹均用原 P6.5 积分器在官方 Frozen FEM 中推进，并独立复核末接受线段与真实 OUTLET 三角面相交。出口语义是球心首次穿过开放端盖，不宣称整个有限球体完全离开。无人工直线补齐、短轨迹延长或 synthetic 通道。
2. **完成数量？** scheduled=2200，admitted=1969，completed=1249。最低 300 条要求已达到；优选 1000 条状态：PASS。所有未完成 ID 保留原始尺寸和明确终止原因。
3. **采集时间？** 从 t=0 到末次计划出生为 94567.677818595 s（26.268799 h）；最终记账到 94567.771039076 s。相邻出生间隔 42.985308099 s。长时保持名义浓度是此数据集的假设，不是一次注射后的有效 PK 持续时间预测。
4. **出口数量？** OUTLET_01：0、OUTLET_02：1249、OUTLET_03：0。未分类的成功退出=0。终止/入场条件对样本有筛选，不能把此出口比例当作未截尾人群的生理分流概率。
5. **独立叠加？** 是。ULM replay is assembled from independently integrated full-vessel microbubble trajectories in a frozen flow field. 每个 MB 使用自身稳定 ID/出生时间/原始尺寸/朝向/入口随机流独立积分；不含 MB–MB、RBC–MB 同时水动力耦合。
6. **旋转与累积？** 已输出四组完整血管旋转、ULM 式累积、三条代表性完整 journey、出口着色动画；均由同一 scene 导出。完整网格包含全部 45,221 个 WALL 三角形及四个开放端盖。累积只纳入当时已经完成的轨迹；未发生的未来轨迹不提前出现。
7. **补上了什么？** 此前真实 P8 只有入口附近 1 ms 单 MB smoke replay。本阶段新增长采集出生账本、全血管独立轨迹积分、每 ID 终止审计、数百/上千条完整 passage 累积与完整单 MB journey。没有把此前 synthetic lifecycle 作为真实通行证据。
8. **还缺什么？** 真实 RBC passage、RBC-resolved suspension、完整 PK、新 CFD、Particle-7.5、新 RBC deformation 均未建立。P6.5 仍限 sphere-normal regularization，非球形 lubrication、多体 mobility、production dt/neighbors 与 handoff event-time 收敛未冻结。没有声学成像或真实定位噪声模型；未覆盖出口分支仍缺轨迹。
9. **是否可展示？** 是，作为真实 Frozen 三维几何中可追溯的独立微泡完整路径可视化层；应连同覆盖缺口、计算截尾与上述模型边界展示。完整生理物理完成、整棵血管所有分支重建不能由此自动推出。

## 数据与数值证据

- 共 819762 个原始保存状态；完成轨迹长度 min/median/max=68.406/69.716/71.306 µm；residence=88.917/95.286/149.586 ms。
- Q_IN=2.73691323909e-15 m³/s，C_MB=8.5e12 m⁻³，H_D=0.45。P7 deterministic cumulative flux clock 原样复用；未改成 Poisson。H_D 是 feed volume fraction，本任务不模拟 RBC，也不强设 tube Hct。
- 原 SonoVue inverse CDF 对每 ID 只抽尺寸一次；P7 isotropic orientation V0 明确是模型假设。位置采用原 P7 局部法向流量权重 proposal；FiniteSizeAdmission 最多 512 次重试位置，不重抽尺寸，不增加入口偏移。通过者分布已被有限尺寸条件化；Figure 02 分开显示初始 proposals 与 admitted 位置。
- 每 ID NPZ 保留所有接受子步，含 elapsed t、位置、上一接受区间使用的速度/角速度、q、h_geom、g_nf、h_lower 和 nearfield state。首行速度是入场自由流诊断；不能当作已执行位移。CSV.gz 另含物理绝对时间、速度范数和 stable ID。绝对时钟用于 acquisition，elapsed t 保留长采集下的子步精度。
- 原 P6.5 接受/拒绝、连续 wall 证书与二分时间细化保持原样；只用精确输入查询缓存减少重复计算。批量开发中扩展了相同 BVH 请求缓存，较早记录使用仅 wall-gap 缓存；两版均不改变物理算术，完整接受数组、终止原因与 provider 次数的对照见 cache_v2_parity.json。永久测试比较 native 与缓存适配器，原 full native ID 1/4 的完整原始数组也作为逐位一致对照；安全停止 ID 13 的全部保留状态与原生 1.5 s 长计算的前缀逐位一致，证明该保护只截断而未改变已接受轨迹。可选全局运动界优化已测但正式队列未开启。
- dt=0.25 ms 沿用 P7 验证配置，只是本数据集数值步长，不是 production freeze；每 MB elapsed horizon=1.5 s。确定性计算保护包括 16,000 provider calls、512 次无接受推进、32 个名义步完全不动、64 个接受子步仅有 roundoff 量级位移或总时间推进低于名义 dt 的 1%。保护只停止计算并保存此前状态，不改物理力、位置、尺寸或接受规则。
- 近壁计算停止不能解释为真实捕获或排除通行：原规则有 2 nm handoff 附近固定状态和纳秒级细分。诊断归档保留原完整 1.5 s 停滞例与开发参考；这些归档不混入正式队列。终止比例依赖计算保护，尚未建立 handoff 事件时间/轨迹完成率收敛。

| 终止类型 | 数量 |
|---|---:|
| OUTLET_02 | 1249 |
| INLET_ADMISSION_GUARD_EXHAUSTED_UNRESOLVED | 231 |
| INTEGRATION_SAFETY_STOP | 720 |

Frozen 的三个出口均有非零体积流量；某出口没有完成 MB 轨迹，不表示该分支没有流动。有限尺寸入场条件化和计算截尾使该队列不能直接替代无偏体积通量示踪。

| 出口 | Frozen 流量 (m³/s) | Q_OUT / Q_IN | 完成 MB |
|---|---:|---:|---:|
| OUTLET_01 | 1.16534179e-16 | 4.2579% | 0 |
| OUTLET_02 | 2.33199207e-15 | 85.2052% | 1249 |
| OUTLET_03 | 2.88386993e-16 | 10.5369% | 0 |

另对预先选定 ID 1、4（基准完成）及 7、13（基准计算停止）分别重算 dt/2 与 dt/4，结果见 [步长敏感性](data/timestep_sensitivity.json)。这些重算不计入正式队列。比较只在共同 elapsed 时间内线性插值计算误差，不外推。四例有限对照不能建立所有近壁事件时间或完成率收敛。
ID 1/4 在 dt/4 的驻留时间分别为 89.165189/95.704231 ms；相对基准变化 -0.355750/-0.747065 ms；共同时间轨迹最大偏差为 0.405972/0.822959 µm。两例仍从 OUTLET_02 退出，但不可把约亚微米偏差表述成零误差。ID 7/13 仍计算停止；其中 32 个名义步的静止保护随 dt 改变等待时长，因此停止时刻不是物理捕获时刻。

## 显示变换与回放

- 所有可视化只读同一 saved scene：线性位置插值、P8 SLERP 朝向插值只作用于显示，不外推、不回写物理状态。轨迹使用保存折线，无样条、无路径抽稀、无轴向拉伸。仅平移显示原点并把米换算为 µm，当前球保留原半径。MB ID 采用屏幕空间引线标注以避开端盖文字；引线是注释，不计入物理轨迹。
- 相机方位 30°→120°、俯仰 28°，正交等比例显示。1/4 动画慢放代表性真实出生到退出窗口并明确 CUT 空闲时间，最后切至最终记账；2 将完整采集窗口压缩为 240 帧；3 慢放并顺序显示三条独立的完整 birth→exit，并标记跨 MB 空闲时段切换；每次到达出口后保持已计算端点状态 12 个显示帧（0.8 视频秒），用于观察退出，不推进物理时间。单 MB 动画侧栏明确只统计当前选定 ID，其余三类侧栏统计完整 acquisition cohort。真实出生稀疏，未通过复制活跃 MB 制造粒子云。
- MP4 1600×1000、15 fps；GIF 960×600，每五帧及首/中/末关键帧取样，333 ms 预览节奏。物理时钟以逐帧 JSON 为准。视频 active 只统计处于有效保存计算区间的轨迹。停止后不再外推；这些微泡后续状态未知，不等同于物理死亡、排出或真实悬浮液中消失。历史线透明度仅是显示淡化，不改变历史样本；accumulation 只画已完成历史，不把未完成尾迹当作成功路径。
- XY/XZ/YZ 累积图按原接受 dt 加权，避免近壁大量细分造成虚假的计数密度。它是计算轨迹占据图，不是超声测量重建。

## 验证与审查范围

本轮实际运行 296 项永久测试，全部通过；包括 P8.1 与 P1/P6.5/P7/P8/此前 full3d 相关回归，并未声称重新执行全部历史阶段测试。锁定的 1768 个上游文件 SHA256 保持不变。
原始尺寸、isotropic q、每 ID 出生时刻、首次入口 proposal 由保存种子重新生成逐项精确比较；每条完成轨迹末段重新分类到官方出口。另对每 admitted MB 最多九个保存的非出口点复查原 FEM 内部性，这项是抽样，非全步额外点定位。
每个 MP4/GIF 已全帧解码；保存帧账本与 scene 状态逐帧比较，验证 IDs、位置、未来历史排除、完成/失败记账和时钟。重复建立两个独立 renderer 导出同一帧，RGB 与物理状态精确一致。这项像素一致性覆盖当前 software_environment.json 记录的运行环境。
代理实际查看十张图及四组已导出视频的首/中/末帧/关键帧，审查说明与文件散列见 visual_inspection.json；这不是用户人工验收，也未声称人眼看完每一视频帧。用户验收状态为 PENDING_USER_REVIEW。

代码提交：`7decc31301fbb27dcb5bc2a8497901c8861389cd`；所有代码与测试及产物绑定。开发中停止的试算保留在 data/reference_*，不计入正式数据量，最终以根 trajectories、trajectory_catalog 与本报告为准。

## 可直接打开的产物

- [本地可视化总览](index.html) · [验证 JSON](PARTICLE8_1_VALIDATION.json)
- [轨迹目录 JSON](data/particle8_1_trajectory_catalog.json) · [逐 MB CSV](data/trajectory_catalog.csv) · [全部样本 CSV.gz](data/trajectory_samples.csv.gz)
- [事件 JSON](data/particle8_1_events.json) · [事件 CSV](data/events.csv) · [物理审计](data/physics_audit.json)
- [测试记录](data/test_run.json) · [媒体解码](data/media_audit.json) · [视觉检查](data/visual_inspection.json)
- [完整文件列表](OUTPUT_FILES.json) · [SHA256 清单](ARTIFACT_SHA256.json)

- 01_full_vessel_replay：[MP4](animations/particle8_1_anim_01_full_vessel_replay.mp4) · [GIF](animations/particle8_1_anim_01_full_vessel_replay.gif) · [帧记录](frames/particle8_1_anim_01_full_vessel_replay_frames.json)
- 02_accumulation：[MP4](animations/particle8_1_anim_02_accumulation.mp4) · [GIF](animations/particle8_1_anim_02_accumulation.gif) · [帧记录](frames/particle8_1_anim_02_accumulation_frames.json)
- 03_single_journeys：[MP4](animations/particle8_1_anim_03_single_journeys.mp4) · [GIF](animations/particle8_1_anim_03_single_journeys.gif) · [帧记录](frames/particle8_1_anim_03_single_journeys_frames.json)
- 04_outlet_ensemble：[MP4](animations/particle8_1_anim_04_outlet_ensemble.mp4) · [GIF](animations/particle8_1_anim_04_outlet_ensemble.gif) · [帧记录](frames/particle8_1_anim_04_outlet_ensemble_frames.json)