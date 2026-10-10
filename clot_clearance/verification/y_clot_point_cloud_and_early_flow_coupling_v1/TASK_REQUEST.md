CODEX TASK — Y_CLOT_POINT_CLOUD_AND_EARLY_FLOW_COUPLING_V1
你已经知道如何访问并使用服务器。请实际执行本轮所要求的软件修改、必要计算、图件和动画；不要只返回计划。
A. 本轮背景与最高目标
项目仍以 LAMMPS BPM/PERI LPS 低成本血栓碎裂 为主，不恢复 B1/B2 phase-field、不开发完整 PD2 fracture。
上一轮 Y_CLOT_CONTINUOUS_HANDOFF_AND_HONEST_RENDER_V1 已证明：
- HANDOFF_POSITION_CONTINUITY=PASS；HANDOFF_VELOCITY_CONTINUITY=PASS；TRUE_VELOCITY_RESET_BUG=NOT_OBSERVED。
- 旧重复播放端点已修复；新视频仍不符合用户视觉偏好。
- 背景 drag 直到 2.33 ms 全局交接后才启用：BACKGROUND_FLOW_NOT_APPLIED_DURING_FRACTURE=true。
- 207 个 1–2 节点未解析碎屑 + 1 个 6 节点碎块在末态游离；真实破坏仍是 LOCAL_DEBRIS_SHEDDING_WITH_ATTACHED_CORE，不是宏观块体破裂。
- 用户明确要求彻底放弃 medical shell 渲染，采用其上传的 Fig_clot_discretization.png 所示 点阵模式。也希望碎屑从实际脱离开始就受到流场影响，而不是到全局 handoff 才启动。
本轮仅完成两项：
1. 把点阵显示设为新的默认科学展示方式，实际渲染比上一轮更统一自然；
2. 在不改变原 LPS 断裂材料/几何/接触的前提下，尝试最小实现“按真实脱离时刻触发背景流体力”，并消除模型切换引起的突然开启阻力和视频时序错觉。
绝对不能把“立即受到流体力”错误解释为“速度瞬时被设成局部 CFD 速度”。
B. 保留科学身份和原资产
保持：
- 原生 Y 血管三维 WALL、入口与出口几何。
- V3 曲壁血栓、LAMMPS 1927 材料节点、三层锚点及全部原始连接。
- 同一 Liu 2024 小应变 LPS：K=33333.333333333336 Pa，G=940 Pa，rho=1000 kg/m³（工程密度）。
- bulk s00=0.15、interface s00=0.140624859375、alpha=0。
- 1.0dx 接触候选；完整 Y 壁面碰撞与原始流场查询；dx=50 µm；delta=3.01dx；dt=0.5 µs。
- 原 Gaussian 200 Pa 工程载荷、1.1605 ms 脉冲、2.33 ms 本构动力学窗口。
- 旧背景 CFD 仍是 UNOBSTRUCTED / NO CLOT PRESENT，不是包含血栓的 FSI 或真实血栓面 traction。
必须记录现有 git status 与全部关键历史目录。禁止 git reset --hard、git clean、git stash、覆盖当前工作区未提交结果。不得删除旧 medical render；仅将其标记为历史可选，不再用于新视频的默认输出。
C. 先做点阵渲染，不需要重新求解
优先只读利用上一轮已保存的 LAMMPS 轨迹、断键 events、component IDs、输运 trajectories，从中重新渲染一部点阵版 MP4，并生成 3 个代表性静态帧。
直接对齐用户提供的 Fig_clot_discretization.png 风格：
- 血栓原始材料节点：暗红/胭脂红、清晰的离散点阵；
- 固定血管壁锚点：浅蓝/青灰小点；
- 真正游离的 1–2 节点或更大分量：中灰或稍深灰点，不再显示橙色实心球；
- surviving bonds：极细、低透明度线条；允许分别开关 bulk/interface/wall 的显示，但所有被隐藏的内容必须在 legend 和 manifest 中说明，求解与统计不能受影响；
- 血管：原生 Y 网格淡灰、低透明度的轮廓参考，不能用不透明管壁遮住内部碎屑；
- 所有节点使用原始三维位置，绝不放大位移、物理碎片大小、合并真实未连接粒子；
- 点尺寸采用固定屏幕像素/合理的深度着色，是图形符号，不将其误称为真实微粒直径；
- 不再使用平滑实体 shell、红色大块蒙皮、琥珀色碎片球、视觉拟合裂缝、假散射体积。
动画同时展示：局部近景、整个 Y 管定位；坐标轴用 mm、真实流向可用少量细线/箭头显示，但必须标注 UNOBSTRUCTED BACKGROUND CFD。
减少冗余字幕。时间标签写明确物理时间，播放率变化必须显示。固定相机稳定，避免自动缩放掩盖细小真实位移。允许为讲清真实轨迹显示半透明历史轨迹线，轨迹每个点必须来自真实数值输出。
输出：y_clot_point_cloud_only_preview.mp4、PNG、显示模式说明、帧清单；用户要求的 point cloud 为默认而非备用。
D. 明确解释为什么之前看起来“先停后走”
先从实际旧数据复述并输出短报告：
- 交接处 x,v 是否连续；
- 1–2节点碎屑的 LAMMPS handoff 速度与 CFD 背景速度尺度；
- 原先 100 µs 数值松弛导致的减速时间尺度；
- 2.33 ms 之前 F_background_drag=0、之后才开启的载荷切换；
- 背景近壁无滑移造成的低流速；
- 旧视频中的真实慢流/快惯性以及播放倍率切换。
不能把“速度快速下降”误报成速度重置，也不能通过把碎片直接附着到流线上伪造物理运动。
E. 最小物理升级：流体力在脱离时立即生效
目标不是做 FSI、也不是改变血流。目标是将已经在 Phase-2 用的单向拖曳近似提早、连续地用于 LAMMPS 中真正脱离的节点/组件。
先检查现有 LAMMPS fix external 壁面回调、破键事件和 surviving-bond/anchor topology。开发时必须以当前计算步的真实已断连接为基础识别 DETACHED_FROM_ALL_WALL_ANCHORS，不能用预知未来的后处理标签判定当前断裂。
建议实现事件触发的脱离状态更新：只在连接拓扑变化后更新附着/游离组件身份，不必每个时间步从零重建全图。如果现有 API 不能可靠即时读取，可采取最小、隔离的原生 hook / callback，不能重新编写 PD 材料求解器。
对真正游离且与背景 CFD 查询域有效的组件：
- 最早从首次允许的物理时间步开始施加背景阻力；
- 1–2节点 NUMERICALLY_UNRESOLVED_DEBRIS 只能用明确标记的工程速度松弛代理，基线 tau=100 µs；不能宣称已知真实碎片水动力直径；
- =3节点组如复用已有等效球阻力模型，采用现有组件体积/质量与速度，正确处理组COM vs 节点力分配，并如实报告其限制；
- 同时保留真实壁碰撞，按当前坐标查询 CFD P1，域外不外推；
- 原 CFD 在原血栓内部和新近壁层的有效性有限，不得把在未求解的实际 clot solid 内的流速当可靠实测；需要标记这些查询为 CLOT_PRESENT_FLOW_UNCERTAIN / NEAR_WALL_UNCERTAIN；
- 已附着主体不应被误当为自由粒子套用 Stokes/Schiller–Naumann 阻力；其原动量与力学约束照常保留；
- 旧数值阻尼 -500m v 和新增 drag 要明确避免对游离节点的重复计入。必要时把数值阻尼限制在声明的部分节点/状态，记录这对物理模型的改变；绝不通过直接赋予流速进行速度重置。
仅当接入方式能可靠保存：detachment event、流体力首次生效步、作用节点、x/v、局部 CFD u、外力/壁力和对应数值阻尼时，才启动新力学 case。
若在本轮约束下无法可靠实现 event-triggered force，必须停在 EARLY_FLOW_COUPLING=NOT_ESTABLISHED，交付真实阻碍和最小接口需求；绝不通过视频插值、原始帧重排、伪造流动和提前剪接输运来冒充成功。
F. 新增力学求解严格限于一个真实 case
只允许重算：
V3 / Y vessel / GAUSSIAN / reference mean 200 Pa / 1.1605 ms pulse / 2.33 ms observation / 1.0dx contact
除新增即时流体作用所必需的项外，几何、材料、断裂、载荷与断后接触保持原配置。旧 case 作为历史对照读入，不重跑旧 baseline。
绝不追加 Gaussian/ Opposing 幅值扫描、材料校准、网格收敛、GPU 压测、VES、PD2、B1/B2 或含栓 CFD。
必须报告的科学差异
- 新旧首次 bulk/interface break 时刻；
- detached node count、native broken bonds、attached core mass、未解析碎屑质量；
- 每个代表碎块的 t_detach 与 t_drag_on（特别报告延迟 Δt）；
- 实际 x,v 是否连续，slip |v-u_f| 随脱离时间衰减；
- LAMMPS 终时交接到长时输运时是否仍存在额外阻力律跳变；
- 是否改变近壁/管壁接触、是否出现不允许穿壁；
- 实际 solver wall time 与渲染时间。
若新增流体耦合改变断键数量，按真实结果报告；不因此回调 s00 或 shear 幅值。
G. 连贯输运和视频要求
不能等所有碎屑到 2.33 ms 才第一次感觉到流体力。仅对真实脱离组件从脱离后开始施力，短时 LAMMPS 动力学结束后再进入轻量输运。
输运初始化必须来自新的 LAMMPS accepted terminal state，且使用同一阻力/速度松弛约定，位置与速度接续；时间用连续物理绝对时间。
后续仍采用最多 5 s 单向背景流输运演示，无 CFD 重算，不人为增加水流速度或任意缩短松弛常数以迎合视觉效果。
附着主体在后续阶段停止力学推进时，必须明确表示“FROZEN REFERENCE SNAPSHOT AFTER 2.33 MS”；不得称其已经物理静止。
新视频：
y_clot_point_cloud_early_flow_coupling.mp4
附一份短 GIF、真实首断时刻与脱离时刻关键帧，以及 velocity/flow timeline 图。
持续记录横跨 phase boundary 的 x,v,Fdrag,time，确保从源计算上消除突变，而不只是播放剪辑技巧。
H. 必需的最少软件核对
仅限：
1. 点阵渲染中 1927 个原 clot 节点都可追溯，anchor/clot/free 身份与当前 frame 的真实 topology 一致；
2. 只有 detach 后才允许触发流体阻力；触发延迟不超过已明确的受限算法步进粒度；
3. 不存在额外速度初始化重置/两次 drag/非法域外插值；
4. 真实 Y 管壁无非法穿透；
5. 视频帧与真实数据 manifest 一致、影像主体没有 shell、没有凭空扩大碎片。
不要重复额外 benchmark、测试矩阵或其他工程科学不相关工作。
I. 文档与软件来源
LAMMPS official:
https://github.com/lammps/lammps
fix external:
https://docs.lammps.org/fix_external.html
LAMMPS Python:
https://docs.lammps.org/Python_module.html
BPM/PERI bonds:
https://docs.lammps.org/bond_bpm_peri.html
BPM/PERI contact:
https://docs.lammps.org/pair_bpm_peri.html
PyVista:
https://docs.pyvista.org/
OVITO Python:
https://www.ovito.org/manual/python/
允许在隔离环境中安装/拉取任何必要包。不要污染平行项目或覆盖当前 LAMMPS 原安装；若增加 hook，需固定源码 SHA、独立 build 和可回退 patch。
J. 交付与状态
新目录建议：
verification/y_clot_point_cloud_and_early_flow_coupling_v1/
results/y_clot_point_cloud_and_early_flow_coupling_v1/
至少：
- TASK_REQUEST.md
- STARTING_GIT_STATE.md
- BASELINE_REASON_FOR_PERCEIVED_STOP.md
- POINT_CLOUD_RENDERING_CONTRACT.md
- EARLY_DETACHMENT_DRAG_DESIGN.md
- DETACHMENT_DRAG_EVENT_AUDIT.csv
- VELOCITY_AND_FLOW_COMPARISON.csv
- LAMMPS_CASE_CHANGE_MANIFEST.md
- FINAL_STATUS.json
- REVIEW_Y_CLOT_POINT_CLOUD_AND_EARLY_FLOW_COUPLING_V1.md
- 点阵 MP4、GIF、重要 PNG/PDF、对应帧来源 manifest
- 轻量 evidence.zip；大体积数据保留原位并保存 hash。
必须单列状态：
POINT_CLOUD_DEFAULT_RENDER=PASS/FAIL
SMOOTH_SHELL_DEFAULT=DISABLED
TRUE_VELOCITY_RESET_BUG=NOT_OBSERVED/OBSERVED
BACKGROUND_DRAG_DURING_FRACTURE=IMPLEMENTED/NOT_IMPLEMENTED
DETACHMENT_TRIGGER_CORRECT=PASS/FAIL/NOT_ESTABLISHED
HANDOFF_POSITION_CONTINUITY=PASS/FAIL
HANDOFF_VELOCITY_CONTINUITY=PASS/FAIL
BACKGROUND_FLOW=UNOBSTRUCTED_ONE_WAY_EXPLORATORY
CLOT_PRESENT_CFD=NOT_RUN
MACRO_FRAGMENT_200UM=OBSERVED/NOT_OBSERVED
LARGE_BLOCK_FRACTURE_MECHANISM=NOT_ESTABLISHED（除非实际新数据证明）
本轮核心停止条件：点阵视图成为正式默认；要么在一个真实 Y 血栓 case 中验证“脱离时开始承受背景流体力”，要么给出无法可靠接入的明确原因；不允许额外扩大任务。
最终用中文回答：用户为什么原来看到先停再动、现在是否从真实脱离时开始受流体影响、碎片速度为什么不可能瞬间等于水流速度、点阵显示是否已经和用户给的示例一致、还有什么物理限制。