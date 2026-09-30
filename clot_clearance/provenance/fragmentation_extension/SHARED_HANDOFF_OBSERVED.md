你现在接手我正在持续开发的脑血管 CFD/FEM → 微泡群体运动 → 科学可视化项目。请在现有工程上继续，不要从零重建。以下为截至 2026-09-29 的交接上下文；下一条消息会给出具体新任务。本条仅用于接手，不代表授权自动启动新的 CFD、微泡生产批次或恢复旧任务。

一、先读取的索引与路径

首先阅读：
/home/lzy/projects/brain_report.md

该文件已核对 276 个本地路径、17 个服务器入口，包含代码功能、数据、报告、动画及历史版本区别。按当前任务读取相关章节和证据，不要无差别扫描所有历史文件。

主要路径如下，后文用这些缩写：
B = /home/lzy/projects/computation_examples/brava_flow_roi_18mlmin
P = B/population_v1
C = /home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d/population_v1
L = /home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/src/particle_3d
F = B/cases/balanced_pressure_final
R = P/runs/realization_0000_contactfix

上轮已将本地保留的流场算例、配套轨迹、代码、数据、日志和展示迁入 /home/lzy/projects/computation_examples，按算例命名。共用代码在 _shared，既有开发/发布工作树在 _archives。原路径保留兼容符号链接，实体已迁移；不要删除链接、再次迁移，或批量改写冻结记录中的旧路径。

迁移报告：/home/lzy/projects/computation_examples/MIGRATION_REPORT.md
路径对照：/home/lzy/projects/computation_examples/PATH_MIGRATION.csv
原始几何/打印工程仍在 /home/lzy/projects/vascular_printing，指南为该目录的 开发与使用指南.md。

二、当前 BraVa 几何与正式 CFD

对象为 BraVa BG001 RMCA、BALANCED 已验收四端口血管芯，包含人工接管，排除倒模盒体。candidate 0、15 是同一物理几何的两种刚体打印姿态；当前无重力、刚壁、牛顿流体模型复用同一个流场，不是两个独立 CFD。

输入身份/单位/端口/变换：B/inputs/geometry.json
权威网格：B/mesh/SV_MESH/mesh-complete.mesh.vtu
边界：B/mesh/SV_MESH/mesh-surfaces/{INLET,OUTLET_01,OUTLET_02,OUTLET_03,WALL}.vtp
网格为 135852 节点、692313 四面体；源几何 mm，计算坐标 m。

当前有效场由 B/reports/ACTIVE_FLOW.json 指定：
- 冻结速度/压力：F/frozen_flow/steady_flow.vtu
- 规范数组：F/frozen_flow/flow_arrays_si.npz
- 原始 WSS：F/postprocess/wall_wss_si.vtp
- 独立检查：F/postprocess/INDEPENDENT_CHECK.json
- 实际设置：F/run/solver.xml、F/run/PETSC_OPTIONS.txt
- 执行和收敛：F/reports/execution.json、F/run/solver.log
- 出口统计：F/postprocess/boundary_summary.csv

流场 SHA256：73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc
网格 SHA256：e755a6ce04eca518fc02df5cd6139c98887e94617bdf9410d72f55b99f443c94

材料：rho=1056 kg/m³，mu=0.00345312 Pa·s，刚性无滑移壁面。
入口规定 Q=18000 μL/min=18 mL/min=3e-7 m³/s，入口压力不规定。
正式出口压力/traction BC：
O1=556.9706555488448 Pa；O2=0 Pa；O3=17.13742715356202 Pa。
求解后的入口流量占比：
O1=33.33342362562679%；O2=33.322943568142236%；O3=33.343632806230964%。
不要把截面平均压力误当规定 traction 数值。

流程：准备几何/网格 → 等流量校准 → 测量出口压力并统一平移参考值 → 三个自然压力出口正式复算。校准实体在 /home/lzy/projects/computation_examples/brava_calibration_equal_split/flow；强制校准分流不能冒充最终自由分流验证。

代码入口：B/scripts/run_flow_sequence.py、solve_remote.py、postprocess_flow.py。实际使用 svMultiPhysics + PETSc，8 MPI ranks、CUDA 向量/稀疏矩阵、FGMRES 与 CPU ASM/LU；CFD dt=0.01 s。GPU 同步修复、独立构建与 CPU/GPU 对照在 B/gpu_solver_fix。不要把构建时的 BUILT_NOT_YET_VALIDATED 状态当成后续没有完成对照验证。

三、必须区分“已有微泡结果”和“当前计算协议”

已有正式群体结果是 R，历史规则为每个出口至少 20 个实际球心穿越后停灌，再观察 12 s。
- C_MB=4e9 m⁻³，名义 dt=0.0005 s。
- 源事件/接受数均为 1565，最大同时活跃 1438。
- 实际停灌：1.2408156603855842 s；结束：13.240815660385584 s。
- O1/O2/O3 出流：458/527/437。
- 141 条观察截止、2 条数值失败、0 条支持静止。
- conservation_pass=true，但 scientific_integrity_pass=false；不能称所有轨迹成功或全部通过。

读取 R/SUMMARY.json、IDENTITY.json、POPULATION_CONTRACT.json、MANIFEST.json。
原始接受子步群体状态在 R/states/chunk_*.npz；source.jsonl.gz、terminals.jsonl.gz、interaction_audit.jsonl.gz 分别保存来源、终态和交互账本。
逐条轨迹与统计在 P/analysis/realization_0000_contactfix/，重点为 tracks/、trajectory_catalog.csv、population_timeline.csv、INDEPENDENT_AUDIT.json、failure_diagnostics.json。

当前代码已经改成固定时间协议，不能沿用旧出口覆盖停止规则：
- 正式入口：python -m particle_3d.population_v1.run_fixed_window_population
- C/run.py 仅兼容转发到新入口。
- 合同：P/fixed_window_protocol/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json
- T_infusion=1.240815660 s；停灌后观察 10.0 s；T_end=11.240815660 s。
- C=4e9 m⁻³，Q=3e-7 m³/s，源率1200/s，dt=0.5 ms，默认 seed=2026092909。
- 出生区间 [0,T_infusion)，停灌/结束不依赖出口结果；终点活跃泡保留末态并记右删失。
- 旧 --minimum-exits、--post-infusion、--fixed-infusion 参数已经禁用。

固定协议已完成小规模 CPU/GPU、调度、源重放和真实场零源 CLI 验证，但尚未运行新的正式大批量 BraVa 轨迹。不能将已有1565条旧协议结果改称新协议结果。

当前协议报告：P/fixed_window_protocol/FIXED_WINDOW_PROTOCOL_REPORT_ZH.md
旧正式群体验证报告：P/POPULATION_VALIDATION_REPORT_ZH.md
P/REPRODUCE.md 含旧协议命令，仅作历史参考；当前入口以固定协议报告及实际源码为准。

四、微泡模型及计算分工

当前为冻结 CFD 单向驱动、过阻尼有限尺寸刚性球、真实有限三角形壁面接触、持续 Poisson 灌注、共享物理时钟及微泡间 V1 法向近场/接触模型。保持原 SonoVue D≤4 μm 条件分布、真实入口正通量采样和一次准入，不按出口结果挑泡或篡改粒径。

主要代码：C/environment.py、source.py、fixed_window.py、engine.py、gpu_engine.py、gpu_field.py、gpu_wall.py、gpu_neighbors.py、gpu_resistance.py、components.py、output.py。
CPU 权威物理依赖位于 L 下，如 field.py、particle9a_motion.py、planar_wall_hydrodynamics.py、wall_geometry.py。依既有 identity 清单保护原科学源码，不因阶段名旧就删除。

GPU 实际执行 FEM/壁面查询、邻居搜索、批量阻力和主要推进；CPU 执行源/准入/严格出口、困难球对与有限边缘小组件回退及输出。旧正式批次约99.83%的接受动力学更新由GPU求解；不代表所有代码全在GPU。

测试：/home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d/tests/population_v1
没有 RBC、声学、壳破裂、聚并、变形泡或双向 CFD 耦合；其它 RBC 参考模型不是本 BraVa 批次。

五、已知近壁数值风险与加速原因

近壁全面审计：P/cfd_wall_audit/WALL_CFD_AUDIT_REPORT_ZH.md
加速诊断：P/analysis/acceleration_diagnosis/ACCELERATION_EXPLANATION_ZH.md
内部通量：B/reports/internal_sections/

全局边界守恒通过，但内部截面偏差最高约占入口Q的2.49%，局部P1近壁单元存在高散度/明显向壁速度分量。典型 tetra 388769 的 |div(u)|/||G||F 约0.995。当前只有单档 BraVa 网格，不能宣称网格无关。

原WSS按相邻四面体完整P1梯度计算 mu*(G+G.T)*n，再作切向投影；逐面片值 WSS_raw_Pa 与面积加权节点显示 WSS_display_Pa 分开。源码在 B/vendor/flow_solver_support/wss.py、wss_case.py。独立重算一致仅证明实现对应，不证明局部物理精度。

已有“先慢后快”主要对应近壁慢滑 → 面片边缘/离壁 → 间隙增大、近壁阻力减弱 → 进入较快FEM流区。不是动画后半段改变时间倍率，也不是惯性蓄能弹开。该机制有记录及CPU状态复核，但受到上述CFD/几何离散风险限制，不能直接解释成已验证的真实生理现象。

六、当前展示版本与必须保持的样式

流场预览：B/visualization/candidate_0/OPEN_RESULTS.html 和 candidate_15/OPEN_RESULTS.html。
每姿态包含流线、局部速度矢量、压力、节点显示WSS、原始面片WSS五段动画。candidate 0 流线已改为36 s/周、整周最大安全取景；其它流场视频保持各自原样。

最新微泡四视角：
P/candidate_0/views_0to1p2s_large_axes/OPEN_RESULTS.html

此版本：
- 只输出 candidate 0；物理时间0–1.2 s。
- 四个相对基准方位角240°、330°、60°、150°，仰角均25°；不是世界坐标绝对方位角。
- 黑色背景、蓝色微泡、灰色半透明血管、英文文字、端口标签、流场风格坐标轴/网格，无红白外框。
- 主图保持当前最大安全取景，修改展示时不得机械再乘1.3导致裁切。
- 右上角标角度；轴名称40 px、刻度30 px、刻度文字间距38 px。
- 1920×1080、24 fps、262帧、约10.916667秒视频；维持既有慢放映射。
- 隐藏2条失败及58条在加速窗口访问高散度区域的数值风险轨迹；58条并非逐条已证实物理错误。
- 每帧继续隐藏当前时刻前后各0.05 s内三维时间加权平均速度<2 mm/s的泡。
- 筛选只改变显示，不删除原轨迹、不修改源/终态账本或出口统计。

渲染：P/scripts/render_population_four_views.py → render_view05_flow_axes.py。
最新字号需参数 --axis-title-size 40 --axis-tick-size 30 --axis-tick-offset 38；默认参数仍是小字号。
核验：P/scripts/check_population_four_views.py。
风格/坐标公共代码：B/visualization/style.py。
各角度 data/render_manifest.json、animation_samples.npz、frame_summary.csv、layout.json 记录来源、掩码和相机；总 data/QC.json 为已完成验证。
不要把旧 candidate_0/animations 或 candidate_15 动画误认为最新四视角版本。

七、服务器、环境及版本

SSH别名：vast4090。使用前读取实例 /etc/vast-agents-guide.md；不要修改其他项目服务。
服务器根：/workspace/brava_flow_roi_18mlmin_20260928
服务器 population：上述根目录/population_v1
服务器核心：population_v1/source/particle_3d/src/particle_3d/population_v1
服务器环境：population_v1/env/bin/python
服务器 CFD 可执行：/workspace/brava_flow_roi_18mlmin_20260928/gpu_solver_fix/svmultiphysics
旧 CFD 后处理环境：/root/particle8_2_runs/env/bin/python
本地科学环境：/home/lzy/projects/computation_examples/_shared/python_environment/bin/python

既有服务器记录为 RTX 4090；实际运行前核对当前资源和进程。GPU 用于真实计算/渲染，视频编码为 CPU libx264；不把渲染或nvidia-smi利用率当作动力学GPU证据。本地环境不默认具备CUDA Torch。

相关仓库：OZAKI39/ulm-3d-vascular-model-generation。
保留开发工作树：/home/lzy/projects/computation_examples/_archives/github_sync_20260928_latest
当前核对分支：dev/brava-population-gpu-20260929
当前核对HEAD：0b808fb6c90e27674d79078c49358ced2e9258ff
旧正式群体生产commit：98d705f13babdfa0cba6c12331623292a588cb77
固定协议CODE_VERSION记录：9836c264863eca7f34806c22e020b4c174ef3a7c

当前HEAD、旧运行身份和协议开发版本不同；运行以逐源码哈希为准。目录迁移和brain_report整理没有提交或推送，不要声称GitHub已同步本地所有最新改动。

八、历史范围与继续工作的约束

- B/microbubble/ 是用户停止的旧独立批次，不得恢复其自动启动链。
- P/runs/realization_0000、realization_0000_guardfix、realization_0000_final 是中止开发记录，不能因目录名带final就选作正式数据。
- 旧小鼠WSS项目的细档真实血管和O2压力扰动仍为SKIPPED_BY_USER；不要恢复。旧三档圆管验证和两档真实血管敏感性结果可查阅，不能写成真实血管网格无关。
- 原小鼠ROI-only、best-feasible、H0、WSS及RBC参考案例已按名称收纳在 computation_examples；详见其README和CASE_INDEX，不要误当BraVa输入。
- 保留原始结果、失败、删失、身份哈希和已接受模型约束；显示筛选不能变成科学数据删改。
- 根据新任务在独立输出目录执行必要工作，优先复用已有结果；不要自动重复完整求解、性能实验或启动额外realization。
- 不为“颜色更均匀”“动画更顺”调压力、平滑物理场或挑选科学结果。诊断区分已运行验证、代码检查、推测和未验证。
- 用户偏好直接完成已授权工作、保持目录整洁；普通可逆实现选择自行判断，不反复索取确认。

接手后先简短说明你已识别的有效流场、历史群体结果、当前固定协议和最新显示入口。若下一条消息已给具体任务，直接按其范围继续；没有具体任务时不要主动启动计算。


## 2026-09-29 用户要求暂停开发（优先于此前继续执行计划）

用户为休眠电脑要求在安全点暂停。`mesh_particle_validation_v2` 已暂停：M2_medium_frame180 CFD 第 12 步（0.12 s）正常退出并保存通过审计的原生重启检查点；本研究服务器任务和自动队列均停止，本地已有检查点及结果副本。等待用户后续明确恢复指令。恢复前先完整读取 `/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/mesh_particle_validation_v2/PAUSE_RESUME_ZH.md` 和 `PAUSED_STATE.json`，不要直接重启旧队列。原任务尚未完成，不将中断 CFD 计为合格稳态结果。


## 2026-09-29 用户已明确恢复开发

用户后续指令“恢复开发”已解除上述暂停。先读 `mesh_particle_validation_v2/RESUME_REPORT_ZH.md` 与 `RESUME_AUTHORIZATION.json`。M2 中档已完成原生 checkpoint step 10→12 连续性验证（速度 relative L2 4.41e-15），从 step 12 在独立 segments/resume_0001 续算。服务器 supervisor `resume_campaign`、`replay_followup`、`report_watch` 已运行，恢复前须先查实际状态，不能重复启动。既有暂停段与原始 FAIL 保留；未宣称稳态或全部研究完成。


## 2026-09-29 持续监控与固定边界补充验证

用户要求开发停止前持续监控并显示进度；不能因服务器任务仍在后台运行而结束当前开发。状态脚本为 study/scripts/monitor_progress.py，记录于 monitoring/LATEST.json 和 status_history.jsonl。20:06 UTC 的 M2 中档已至 step 71，最近 step 70 达 4/5 稳态区间，尚未宣称完成。

新增 M1_fixedplc_medium / fine（762096 / 895647 tetra）已通过完整 QA，逐面保留 M0 的全部边界三角形；固定 1059 个采样点、4 个原异常位置及 4096 个原零速单元体积探针均在域内。实际尺度报告 FIXED_PLC_REFINEMENT_RESULTS_ZH.md；仍有全域内部重划分，不能等同于整体近壁加密。

新增 supervisor 程序 fixedplc_validation_followup 已排队，先等待 resume_campaign 的原 M2/M3 CFD 队列结束，再用相同 CFD 与 46 个初态逐档验证。与其他 CFD/replay 共用 HEAVY_COMPUTE.lock。其加入 ADDITIONAL_COMPUTE_PRODUCERS.json，report_watch 会等它完成。实际工作状态须从 supervisor 和结果记录核对；若再次要求暂停，应把此等待队列一并安全停止，不能遗漏而让其随后自动开算。没有更改原物理模型或发布新 production flow。


## 2026-09-29 20:23 UTC：M2 中档续算完成

M2_medium_frame180 已在 step 81 完成并通过独立后处理与 field audit，flow SHA256 `45575ad4fe989e11c500c190817360a8fc01734a6ffe415254c1729fea0d3ef3`；父目录 RESUME_INDEX 为 QUALIFIED_STUDY_FLOW，原暂停 FAIL 仍保留。全域严格零速单元为 0，全部冻结采样记录仍在域内；局部散度并非全面改善，详见 M2_MEDIUM_FIELD_REVIEW_ZH.md 及同物理截面图。M2_fine_frame180 已启动，M2 中档 replay 正等待共享重计算锁；当前科学任务仍未完成，用户持续监控指令仍生效。


## M2 细档与中档重放完成；M3 正在计算（更新 UTC：2026-09-29T23:02:58.367107+00:00）

M2_fine_frame180 在 step 81 完成 CFD，8676.749338 s，执行/独立后处理/全体场审计均完成。flow SHA256 `8269b9fc48935c8ac3c681459a32a49df07d8f9c62945563ebeda291cbdcbb4f`，mesh SHA256 `8dbb49224c21cf181adc9c93a709f6bfba6fb9c4bf1f9dc0e6d81a8e9140802d`。中细档均无严格零速四面体，1059 条冻结记录全部在域内。M2 中细档的场分量检查 FAIL：全体固定点 u/G/WSS 相对 L2 差为 0.476%/26.514%/4.347%；原高散度区域为 1.176%/55.708%/8.370%。最大 10 条梯度差记录占差平方和 86.012%，壁距约 11.75–34.25 µm，未删点。

入口 WSS 峰值位于同一壁面三角形 78835，距原 388769 中心约 62.34 mm。壁邻接单元三个壁节点为零，唯一非壁节点属于 INLET。高度 25.818→17.072 µm，切向速度 0.281413→0.277222 m/s；μ|u_t|/h 精确重现 WSS 37.6385→56.0739 Pa（相对差 <1e-12）。这是 Flat + Impose_flux + Zero_out_perimeter 下入口离散过渡的网格敏感性，不把该峰值误当原下游热点或真实剪切增强。原 BC 与速度场未改动。

M2_medium_frame180_dt0.00050000 的 46 条固定单泡全部 EXITED：O1/O2/O3=13/18/15，无失败/删失/有支撑停留，耗时 343.586009 s。仍有 16 次 slow-fast 事件，涉及 9 条轨迹；ID636 仍有一次，ID393/1189 未再筛出但出口转为 O1。相对同协议 M0，37 条两次均出流者中 24 条改变出口，另 9 条终止类别改变。故全部出流不等于收敛，也不能外推总体失败率；ID1505 在 M0 本就出流。M2 细档 replay 尚未执行完，不能把场分量 FAIL 混写成已完成全部轨迹 gate。

研究目录 mesh_particle_validation_v2 下新增报告 `M2_MEDIUM_FINE_FIELD_REVIEW_ZH.md`，机器记录 `04_field_audit/M2_MEDIUM_FINE_FIELD_REVIEW.json`，图 `figures/M2_fixed_point_convergence.png/pdf`，已集成 OPEN_RESULTS.html 与总报告。相同物理截面已包含 M0/M2中/M2细。最新本地收集为 4 个候选 CFD + 4 个 replay；HTML 60 个本地链接与 34 个图文件哈希已校验。45 个 study 脚本当前版本已归档并在服务器逐一验证；只改审计/报告工具，无 production 核心变更。

服务器 resume_campaign 当前运行 M3_medium_frame180，M2 细档 replay 等待 HEAVY_COMPUTE.lock，固定边界补充队列仍等待原 CFD 队列。持续监控指令仍生效，不应结束开发或另起重复任务。后续以实时 supervisor/monitor 文件为准；没有放行 production、dt study 或 P2。若再暂停，须连同 fixedplc_validation_followup 等自动队列安全停止。
