# CODEX TASK — Y_CLOT_CONTINUOUS_HANDOFF_AND_HONEST_RENDER_V1

请在已有项目工作树中实际执行，完成审查、必要修正、离线重新渲染、少量条件性求解与交付；不是只给建议或重新生成旧动画。

## 0. 新阶段定位与必须尊重的事实

本轮针对上一轮 `Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1` 的用户反馈：

1. 游离碎屑的医学视觉效果不自然，像散落的橙色塑料颗粒；
2. 视频在 LAMMPS 断裂结束后看起来“先停住，再被背景流带走”；
3. 用户希望看到具有体积的血栓块破裂、脱落与运动，但当前仍主要是细小碎屑剥落。

必须正视已保存的原始科学事实：

- V3 连续血栓体积约 0.271574388 mm³；LAMMPS 1927 个节点，dx=50 µm，体积约 0.240875 mm³；
- 上轮仅 C0 零载和 C1 Gaussian 200 Pa 两例 LAMMPS 动力学；
- C1 末态附着质量 88.323819%，游离质量 11.676181%；
- 208 个游离连通分量，共 225 节点；其中仅一个 6 节点分量，其余 207 个为 1–2 节点未解析碎屑；
- 全部原生断键时刻重建图得到的最大短暂游离分量仅 13 节点；从未出现 >=200 µm / >=300 µm 游离大块；
- LAMMPS 断裂只计算至 2.33 ms；随后将末态交给最长 5 s 的单向背景 CFD 输运；
- LAMMPS 断裂阶段未连续加载背景流体阻力；输运阶段才加入冻结的、未包含血栓的 CFD 场；
- 1–2 节点质量包的输运含工程数值速度松弛 100 µs，不能冒充物理颗粒尺度或真实流体动力学；
- 完整 Y 网格管壁碰撞本轮此前已经通过限定范围内的零越壁检查；旧 CFD 尚不是含血栓 CFD。

因此**本轮禁止通过放大粒子、聚合不连通组件、手动画裂口、造假大块或假造运动轨迹来让效果“更真实”**。必须把“显示质量修正”和“实际碎裂机制尚未产生大块”分开。

## 1. 从真实工作树开始，保护全部历史

检查上轮 REVIEW、FINAL_STATUS、`fragment_handoff.py`、`fragment_transport.py`、`medical_fragment_renderer.py`、`native_y_wall_contact.py`、视频逐帧 manifest、LAMMPS 原生轨迹、surviving bond/fragment lineage、输运轨迹与 raw/medical 对照。

记录 `STARTING_GIT_STATE.md`。禁止 `git reset --hard`、`git clean`、`git stash`、覆盖未提交工作、删除旧 outputs 或修改平行项目的 CFD 文件。旧 244 帧 MP4 及其全部数据只读保留。

继续冻结：

- Liu LPS K=33333.333333333336 Pa、G=940 Pa；s00_bulk=0.15、s00_interface=0.140624859375；
- `contact_distance_factor=1.0dx`；原三层 anchor 与真实 Y WALL ；
- 原背景 CFD 身份、坐标、SI 单位与端口；
- B1/B2 TERMINATED，PD2 fracture STOPPED，FSI 和 bubble streaming 未接入。

## 2. Gate A — 必做：定位“停住再运动”的真实原因（完全不新增 LAMMPS 求解）

从**原始求解数据，而非 MP4 像素**逐条核查：

A1. 对每个原生游离 component 以及已独立的 6 节点分量，记录：

- 最后 3–5 个 LAMMPS accepted frames 的 COM `x(t),v(t)`；
- 末态 2.330 ms handoff 的 COM、速度、质量、成员 atom ID；
- 输运 `t=0` 的 COM 和速度；
- 后续至少 5 个真实输运步的 `x(t),v(t),u_CFD(x),drag_force`。

A2. 计算并报告：

- `max_handoff_position_jump_m`；
- `max_handoff_velocity_jump_m_s`（分别按 >=3 节点和 1–2 节点质量包报告）；
- `min/max/median` 游离组件速度；
- LAMMPS 末态至输运早期的运动方向/速度变化；
- 100 µs 数值松弛是否造成速度迅速衰减/重定向；
- 近壁 CFD 很低的速度是否令最初位置变化低于像素分辨率；
- 视频在 phase handoff 附近是否插入重复静态帧、停顿时长和是否重置时间标签。

A3. 必须区分：

- `TRUE_VELOCITY_RESET_BUG`：读取/单位/初始化造成真实速度被清零或突跳；
- `PHYSICAL_MODEL_FORCE_DISCONTINUITY`：背景 drag 只在 handoff 后施加；
- `PRESENTATION_TIMEWARP`：前 2.33ms 放慢展示，后 5s 用另一时间压缩，造成感知停顿；
- `NEAR_WALL_SLOW_FLOW`：近壁本来速度小；
- `UNRESOLVED_DEBRIS_PROXY_RELAXATION`：100µs 工程模型；
  多个原因可同时存在，必须给证据。

若证据表明速度未被重置，不能为了迎合用户而声称存在 velocity reset bug。

## 3. Gate B — 必做：不用重算 CFD/LAMMPS 的无缝时间叙事与视频重制

用已保存的真实状态重制 MP4/GIF。原则：**画面播放可以变速，物理位置、速度、断裂拓扑绝不能伪造**。

B1. 保留两阶段真实时钟，明确显示：

- Phase 1 `LAMMPS fracture: physical t = 0–2.330 ms`；
- Phase 2 `ONE-WAY BACKGROUND FLOW: physical t = 2.330 ms + transport_elapsed_s`；
- `Playback time is not physical time`，分别标注演示阶段的时间压缩/慢放。

B2. 不在 handoff 处重复多帧末态/初态，不设置人为 hold。视觉镜头与实际碎片位置连续；如果两阶段同一位置确实连续，应直接承接。流线或箭头的出现应明确属于第二阶段的“载荷模型启用”，不能利用叠加淡入来暗示原 LAMMPS 已计算真实流动。

B3. 不能为了更流畅而任意插值断键或凭空移动碎块。如果采用仅用于渲染的相邻帧位置插值，必须只在同一拓扑、同一成分成员集合和合法的时间区间内进行；原始的断键发生时间不许被平滑/改变。相机运动也不许造成似乎碎块突然位移的错觉。

B4. 输出可直接验证的 `handoff_velocity_audit.csv`、`handoff_timeline_audit.json`、`no_artificial_pause_report.md`。若确有真实物理速度停滞，保留并如实报告，绝不能靠动画伪造运动。

## 4. Gate C — 必做：重做不虚假的医学展示（先不新增力学求解）

现有“每节点体素支持 → 合并外壳 + 小粒子橙色球”的画法对极小碎屑不自然。

采用两套同步、可切换的表示：

### C1. Attached main clot 与真正有几何支持的 resolved component

- 每个 surviving-bond connected component 独立重建外壳，绝不能跨组件自动 bridge；
- 附着主体用更柔和的表面法向/材质、真实的非均匀深红色、适度亚表面视觉效果；不产生显著的人工颗粒方格；
- 形态纹理可以是**显示层 shader/roughness**，不可改求解坐标、体积、质量、裂缝位置；
- component 若只有 3–6 个节点，不得因为圆滑渲染变成几百微米的“巨型软组织块”；
- Y 血管保持真实几何与当前 native wall 碰撞，适度半透明/剖切；
- 统一 1x 长度比例、同步局部/全 Y/分叉镜头，避免眼花缭乱的快速相机缩放。

### C2. 数值未解析碎屑（1–2 节点）

- 用**真实节点位置/质量/成员身份**绘制非常细小、半透明的 mass markers，或由真实离散质量位置构成的透明浓度云（仅医学展示层）；
- 取消像彩色糖粒/橙色弹珠的实心粗粒外观和超长拖尾；
- 浓度云必须标注 `DISPLAY-ONLY DENSITY OF NUMERICALLY UNRESOLVED MASS MARKERS`，不代表经过验证的微米级颗粒尺寸、扩散或实际分子浓度；
- 不许删除这些组件，不许把它们拼成大碎块；质量/位置/速度统计仍从原始 solver/transport 得到；
- 同帧 raw solver view 必须可审查全部原始点、surviving bonds 和分量。

### C3. 面向用户的清晰画面

优先把视频重点放在：

- 原始附壁不规则形态；
- 出现局部裂口的真实位置；
- 小规模 detached fragments / resolved vs unresolved 的区别；
- 随背景流向分叉运动（仅 ONE-WAY DEMO）。

如果仍没有 >=200µm 游离块，画面必须如实呈现，不准起标题“large-clot fragmentation achieved”。

可自由使用已有 PyVista/VTK/ffmpeg；只在必要时安装隔离依赖。渲染缓存场景、表面和同一分量 mesh，避免每帧完整重建无关对象；不做无意义 benchmark。

## 5. Gate D — 条件执行：真实流体交接不连续时，只做最小物理修正

只有 Gate A 发现可证实的 velocity reset、单位错位、插值初始化错误时，修复对应 bug，重新从原 handoff 检查点进行**输运重算**，不重跑 LAMMPS。

如果位置速度连续，但只是背景流拖曳在 handoff 才启用：

- 不把它冒充 bug；
- 明确记录 `BACKGROUND_FLOW_NOT_APPLIED_DURING_FRACTURE`；
- 本轮优先通过 B 阶段诚实的时间叙事改善感知；
- 可评估未来把背景单向拖曳提前到“component 真正脱离之时”所需的最小耦合契约，但**未经验证不能静默把两阶段模型说成全程流固耦合**。

只有当已有条件已经支持一个局部、明确且可验证的实现时，才允许**最多一次**与当前物理载荷身份明确区分的对照 LAMMPS+transport；否则停止并列为下一轮物理耦合需求。不能为让视频好看而增加重力、放大流速、修改拖曳、任意加速颗粒。

## 6. Gate E — 必做：对用户想要的“成块破裂”给出真实机制诊断，而不是假视觉

本轮不扩大力幅，不新增断裂材料，不重调 s00，不加入 VES/softening，不增加高保真 FEM/PD2。

直接从 C1 原生破坏轨迹定位：

- 断键是否集中于暴露表面的浅薄区域；
- bulk vs interface 断键空间分布；
- 可断 interface 是否形成宏观脱粘前缘；
- 是否出现横贯血栓/包围一块较大体积的连通裂面；
- 所有 detached-component 的首次出现位置、最大 nodes、质量、生命周期；
- 后期最大附着主体仍存在的原因（以连接拓扑证据说明，不得断言唯一因果）。

输出：
`CURRENT_FAILURE_MODE = LOCAL_DEBRIS_SHEDDING_WITH_ATTACHED_CORE`（若证据仍支持）。

解释为什么这不是“碎裂成几个大块”：目前 Gaussian 200 Pa 载荷、Liu-LPS critical stretch 和曲面锚定组合没有形成宏观切割/剥离路径。任何“更强力一定大块”的说法都不成立。

只给下一阶段一个**物理上明确且独立标记**的候选：例如在同样的真实 Y 几何下，研究具有宽尺度空间梯度的切向拉扯/界面剥离机制，或以后从真实微泡场确认存在的法向 lift/impulse；若需要人为预裂缝，必须单独明确 `PRESCRIBED_FLAW_DEMONSTRATION`，不得把它装成真实血栓自然破裂。**本轮不为产生大块而启动额外盲目求解。**

## 7. 验收与交付

新增目录：
`verification/y_clot_continuous_handoff_and_honest_render_v1/`
`results/y_clot_continuous_handoff_and_honest_render_v1/`

至少交付：

- `TASK_REQUEST.md`
- `STARTING_GIT_STATE.md`
- `SOURCE_AND_PRIOR_RESULT_REGISTER.json`
- `HANDOFF_CONTINUITY_AUDIT.json`
- `handoff_velocity_audit.csv`
- `VIDEO_TIMEWARP_AUDIT.json`
- `FRAGMENTATION_MECHANISM_DIAGNOSIS.md`
- `DISPLAY_VS_PHYSICS_CONTRACT.md`
- `RENDER_MANIFEST.json`
- `PERFORMANCE_RENDER_ONLY.json`
- `FINAL_STATUS.json`
- `REVIEW_Y_CLOT_CONTINUOUS_HANDOFF_AND_HONEST_RENDER_V1.md`
- 原始/医学对照图至少 2 张 PNG/PDF；
- 重制视频 `y_clot_fracture_to_flow_seamless_v2.mp4` 和预览 GIF；
- 轻量 evidence.zip（原始轨迹如大则保留原路径与 hash）。

最终状态分别给出：

- `HANDOFF_POSITION_CONTINUITY = PASS/FAIL`
- `HANDOFF_VELOCITY_CONTINUITY = PASS/FAIL/NOT_ESTABLISHED`
- `ARTIFICIAL_VIDEO_HOLD = FOUND_AND_FIXED/NOT_FOUND`
- `BACKGROUND_DRAG_DURING_FRACTURE = NOT_INCLUDED/INCLUDED_WITH_EVIDENCE`
- `UNRESOLVED_DEBRIS_RENDER_HONEST = PASS/FAIL`
- `CLOT_SURFACE_VISUAL_QA = PASS/FAIL`
- `WALL_NO_PENETRATION = PRIOR_SCOPE_ONLY/RECHECKED_IF_PHYSICS_CHANGED`
- `MACRO_FRAGMENT_200UM = NOT_OBSERVED/OBSERVED`
- `MACRO_FRAGMENT_300UM = NOT_OBSERVED/OBSERVED`
- `LARGE_BLOCK_FRACTURE_MECHANISM = NOT_ESTABLISHED/ESTABLISHED`
- `CLOT_PRESENT_CFD = NOT_RUN`
- `ONE_WAY_TRANSPORT = EXPLORATORY_ONLY`
- `BIOLOGICAL_FRAGMENT_SIZE_VALIDATION = NOT_ESTABLISHED`

## 8. 停止规则

本轮目标是**真实解释并消除能够修复的“停一下再走”问题、显著改善医学展示品质、诚实揭示细小碎屑与大块破裂之间的差距**。

不得为了得到想象中的大块，修改材料、掩藏粒子、把短期小组件渲染成很大的 clot chunks、造假裂面、制造人为轨迹、偷偷重新跑 CFD、做广泛参数扫描或恢复高保真 fracture。

完成后停止等待人工审查。