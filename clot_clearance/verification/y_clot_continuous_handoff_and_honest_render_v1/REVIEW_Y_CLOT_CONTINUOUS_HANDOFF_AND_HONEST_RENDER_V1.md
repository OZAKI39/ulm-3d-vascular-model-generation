# Y_CLOT_CONTINUOUS_HANDOFF_AND_HONEST_RENDER_V1 审查

本轮已完成原生交接审查、原生断裂机制诊断、离线重制MP4/GIF及raw/medical对照。没有新增LAMMPS、输运或CFD求解；历史244帧视频、上轮数据与代码保持只读。结论仍是 **LOCAL_DEBRIS_SHEDDING_WITH_ATTACHED_CORE**，未实现“大块血栓破裂”。

## 直接查看

- [新版MP4](../../results/y_clot_continuous_handoff_and_honest_render_v1/animations/y_clot_fracture_to_flow_seamless_v2.mp4) / [GIF](../../results/y_clot_continuous_handoff_and_honest_render_v1/animations/y_clot_fracture_to_flow_seamless_v2.gif)
- [初态raw/medical](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/raw_vs_honest_medical_initial.png) / [交接态raw/medical](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/raw_vs_honest_medical_handoff.png)，均另有PDF。
- [速度与播放时钟诊断](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/handoff_speed_and_playback_audit.png) / [破坏机制诊断](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/failure_spatial_and_connectivity_diagnosis.png)
- [显示契约](DISPLAY_VS_PHYSICS_CONTRACT.md)、[最终状态](FINAL_STATUS.json)、[逐组件交接表](handoff_velocity_audit.csv)、[完整帧清单](RENDER_MANIFEST.json)。

## Gate A：没有真实速度重置

直接重读native LAMMPS最后5个保存帧（2.290–2.330 ms），逐组件核对handoff和transport t=0，并审查2471个早期真实accepted steps（每组件至少5步，含x、v、P1 CFD速度、阻力、墙约束标记）。208组件的位置最大跳变 **0 m**，速度最大跳变 **0 m/s**，包括6节点分量和207个1–2节点分量。其质量总计2.8125e−8 kg，成员225节点。

| 类别 | 数量 | handoff速度 min / median / max (m/s) | 所处CFD速度 min / median / max (m/s) |
|---|---:|---|---|
| 1–2节点质量包 | 207 | 0.095294 / 0.858503 / 1.568498 | 0.000320 / 0.004585 / 0.009877 |
| 6节点分量 | 1 | 0.599820 / 0.599820 / 0.599820 | 0.004657 / 0.004657 / 0.004657 |

原阻力在当前实际pre-step状态重建：F_start=m(u−v)/tau，平均力来自该冻结局部流速下指数松弛增量；它不是新增求解。未发生墙约束的早期步速度复现误差上界2.255e−16 m/s，3个带墙约束的步另标记，不能把约束冲量误算为drag。

在“各组件首个不早于500 µs的真实接受步”上，1–2节点速度/初速度比中位数0.006822，方向变化中位数18.727°；6节点比值0.112257、转向0.151°。**这些并非同一500 µs时刻**：采样时刻范围0.502–5.496 ms，中位数0.934 ms；逐行精确时刻在CSV。它们支持快速松弛，不能包装成更高时间分辨率的实际轨迹。

原因分开记录：TRUE_VELOCITY_RESET_BUG=NOT_OBSERVED；背景drag在handoff开启的力模型切换、100 µs工程代理松弛、旧视频2000倍播放率切换均存在。19组件处局部CFD低于1 mm/s，全部208组件单按当地CFD估算的100 µs位移低于旧局部视图1像素；这只是背景平流的像素界限，并不否认handoff初始惯性位移。

另外，1702附着节点的末态COM速度仍为0.054716 m/s，而Phase2未继续求解其弹性运动，只保留末态。主体“停住”有这个明确的模型截断原因，不能称为输运速度清零bug，也不能声称其已物理静止。

## Gate B/C：视频与医学显示的修正

旧视频只查到一处重复端点：零基帧117/118多出50 ms显示，总端点显示100 ms，没有长时间人为hold证据。新版删除重复端点，234个原生保存帧直接接250个输运显示时刻，25 fps共19.36 s；帧233=2.330 ms、帧234=2.340 ms，首个跨阶段步10 µs，与前段4000倍慢放一致。之后仅改变选帧时刻，画面持续给出绝对物理时间与每帧慢放/压缩倍数，不重置为0。不能把自然减速全部“消除”。

输运画面只对相邻真实accepted steps之间同一成员组件的COM线性插值；52000条组件/帧出处可追溯，全部225游离成员显示中心均经native域/端口合法性检查。没有LAMMPS位置或断键拓扑插值。所有镜头固定，流场箭头仅在Phase2出现，字幕明确未含血栓的单向背景流、附着主体末态快照及出口终止标记。

医学显示复用原各组件独立外壳，仅改深红材质、表面法向与柔和光照；关键网格逐顶点/逐面/体积与原重建完全相同。细小组件不拼接、不膨胀。橙色实心球改为真实节点上的2 px半透明质量标记；末态220标记包含219游离未解析节点和1附着单节点。加上1701节点主体和6节点壳恰好覆盖1927节点，质量未删。6节点参考体积等效直径约112.725 µm；该变形组件显示包围盒为266.872×141.507×192.018 µm，与原外壳逐点相同。最长方向大于200 µm不等于体积等效直径达到200 µm，也不是渲染新增膨胀，不能冒充几百微米体积的大块或经生物验证的组织碎片。

首次开发渲染的画面缓存未刷新，已拒收保留；修复后检查不含时间字幕的实际画面变化及编码后解码帧，避免“只有时间在走”。材料视觉QA只确认比上轮柔和且忠于数据；50 µm离散支撑仍限制轮廓，并非真实组织光学验证。

## Gate D：未触发额外求解

没有速度/单位/初始化bug，故不重算transport。也没有经过验证的全程detachment-triggered drag实现，未运行可选LAMMPS+transport。未来最小契约见FUTURE_BACKGROUND_DRAG_COUPLING_CONTRACT.md，包含原生脱离图、数值阻尼与流阻关系、质量/力/单位、墙约束和附着主体范围。本轮没有声称FSI。

## Gate E：为什么还没有成块破裂

11991/82270条bulk键、819/22040条interface键断裂。按参考键中点到真实初始露出面的距离，**99.2828% bulk断键位于100 µm内**；0–50 µm层断键率53.016%，50–100 µm层19.552%，100–150 µm层0.449%。因此浅表剥落不只是渲染错觉。该深度是初始表面距离代理，不是连续裂面厚度。

960节点初有interface连接，16节点失去全部直接interface键，未建立宏观脱粘前缘。总附着1702节点需精确区分为 **1701节点最大核心+1附着单节点**；上轮汇总的88.323819%是总附着质量。最大核心仍有70261条bulk键、21220条interface键并连接1282个anchor。三轴中央80%范围的45个截面至少仍有184条核心bulk键穿过，说明保留通路，不能据此声称任意方向不存在局部裂纹。

重建全部2030个原生断键时刻，得到291个精确成员集合生命周期；峰值游离块13节点、参考体积等效直径145.865 µm。末态208游离分量仍只有一个6节点，其余207个1–2节点。没有≥200 µm（至少34节点）或≥300 µm（至少114节点）的游离大块。7个极短生命周期没有匹配保存帧COM，明确置空；native触发断键位置、参考COM、成员、质量和生灭时刻保留，不伪造首次出现的COM。

目前Gaussian参考200 Pa、Liu-LPS critical stretch、曲面锚定组合没有形成隔离宏观体积的切割/剥离路径。连通图说明残存路径，不证明唯一因果；提高力幅也可能加重节点级粉化，不能保证大块。下一阶段只列一个独立候选：同一真实Y、固定材料上的“宽尺度梯度切向剥离对照”，先定义载荷身份及连贯前缘指标，本轮不运行。

## 冻结范围、开销和交付

K=33333.333333333336 Pa、G=940 Pa、s00_bulk=0.15、s00_interface=0.140624859375、contact_distance_factor=1.0dx、原三层anchor及native WALL均未改。V3连续体积0.271574388 mm³，离散1927节点/0.240875 mm³保持。B1/B2 TERMINATED，PD2 fracture STOPPED，FSI与bubble streaming未接入；含血栓CFD NOT_RUN。

原有近壁不确定、保留主体附近路径不确定、无流体反馈、无片间碰撞和无转动/变形等限制仍在。WALL_NO_PENETRATION=PRIOR_SCOPE_ONLY；本轮显示插值合法性检查不扩张上一轮碰撞验证范围。生物碎片尺寸仍NOT_ESTABLISHED，ONE_WAY_TRANSPORT=EXPLORATORY_ONLY。

记录的只读交接审查约169.0 s，机制审查49.2 s；通过版本渲染/静态图106.5 s，峰值renderer RSS 1.25 GiB（不含encoder子进程）。另保留被拒首次渲染92.0 s，不将开发成本隐去。均不是性能扫描；详见PERFORMANCE_RENDER_ONLY.json。

568个历史文件stat保护及关键数据SHA核对结果见FINAL_INTEGRITY_CHECK.json；源码、CSV、PNG/PDF、ParaView、报告与视频分别保留。evidence.zip不复制大型原轨迹、全部中间shell缓存或被拒视频；原路径与hash完整记录，MP4/GIF作为独立交付并列hash。本轮结束，等待人工审查。
