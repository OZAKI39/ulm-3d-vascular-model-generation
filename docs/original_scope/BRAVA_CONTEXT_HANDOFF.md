<!-- BRAVA_SERVER_MIGRATION_20261009_START -->
## 2026-10-09：BraVa CFD 已迁至新服务器，独立环境验证通过

当前 BraVa CFD 服务器为 `root@95.3.33.46:15283`；使用本机已有 `~/.ssh/vast_step3b_ed25519` 密钥。独立根目录 `/workspace/brava_cfd`，项目 `/workspace/brava_cfd/project`，Python `/workspace/brava_cfd/envs/cfd`。先 `source /workspace/brava_cfd/activate.sh`，当前求解器/MPI 入口为 `bin/brava-solver` / `bin/brava-mpi`。旧 `vast4090` 及下文旧服务器路径是历史记录，不是当前运行入口。

本轮仅迁移 CFD 工作集约 2.42 GB；旧 45.38 GB 完整服务器备份仍位于本地 `server_migration/brava_vast4090_20261009T094459Z`。没有迁移声学项目、旧微泡批次/动画、废弃加密批次或重复工具链。当前本地 Git `a4b238b27437e1760b4f74fc219ede4f0180ee8f` 的相关代码优先，原服务器执行证据按身份补充；远端工作集不含完整 Git 历史。

新机 RTX4090、驱动550.120、CPU配额15.36核、内存42.07 GiB；原 CUDA13 binary 保留身份，新入口按相同 production 源码与 PETSc3.25.5 重建为 CUDA12.6/SM89。文件校验、Python独立性、MPI/PETSc CPU/CUDA、自有P1 sampler和已有小管道单步CPU1/GPU1/GPU2检查通过；速度相对L2差约3.5e-15。仅 PETSc 自动生成路径配置及16条内部链接按新目录适配，原科学数据/源码不变。

科学状态以项目 `ACTIVE_DEVELOPMENT_SCOPE.md` 的最新2026-10-06记录为准，仍为 `STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT` 等原门槛；下方更早的阶段状态仅作历史。环境验证不是局部场精度认证，没有新 BraVa/P2/轨迹/群体计算。安装与验证作业及专用 supervisor 均已退出，后续等用户具体开发指令；新边界仍使用中心线阻力估算出口压力。

[迁移与验证报告](server_migration/brava_new_server_20261009/MIGRATION_AND_VALIDATION_ZH.md) · [当前服务器机器状态](server_migration/CURRENT_BRAVA_SERVER.json) · [连接脚本](server_migration/brava_new_server_20261009/connect.sh)。新服务器自身说明为 `/workspace/brava_cfd/README_ZH.md`。
<!-- BRAVA_SERVER_MIGRATION_20261009_END -->

<!-- PRODUCTION_WORKFLOW_CURRENT_STATUS_START -->
**2026-10-04 当前：REFERENCE_PIPE_FIELD_FAIL（CFD_RUNTIME_OR_EXECUTION / TIME_LIMIT），已停止。VALIDATION_ONLY=true。**

单次端面core径向3→6修正PASS（速度0.210666%、Q 0.165549%）；唯一V2全部mesh gates PASS，永久冻结。wall altitude24.9983 µm，中央bulk215.4203 µm/P0=1.173873。唯一CFD在1800s硬限内完成15步，只有第10步快照，没有稳态field；局部u/G/s/curl未评估，oracle/trajectory/dt均未运行。

修复了validation reader的float32序列化严格对应与监督器pidfd退出异常。CFD没有重启，输入和结果未改。监督记录显式恢复，原生exit code与完整内存峰值未知。所有本次native进程已停止，TASK.json移除。不要续算旧CFD、改reference mesh、生成V3或恢复其它路线/BraVa。

入口：[ACTIVE_DEVELOPMENT_SCOPE.md](computation_examples/brava_flow_roi_18mlmin/ACTIVE_DEVELOPMENT_SCOPE.md)；[总报告](computation_examples/brava_flow_roi_18mlmin/production_workflow/PRODUCTION_WORKFLOW_REVIEW_ZH.md)；[图表](computation_examples/brava_flow_roi_18mlmin/production_workflow/OPEN_RESULTS.html)。下文仅历史；原production physics、固定协议与旧1565（141截止、2失败）身份不变。0.5 ms仍无canonical convergence结论。
本轮提交：mesh `9c8a7551a8c9`；CFD停止与接口修复 `82dd025593b0`。
<!-- PRODUCTION_WORKFLOW_CURRENT_STATUS_END -->

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


## 2026-09-30 00:04 UTC：入口离散剖面与出生速度语义补充

M0/M2中/M2细的 INLET 总节点数为 57/113/135，非 WALL 节点数为 35/91/113；非壁入口节点平坦速度分别 0.3481747072/0.2817099407/0.2775147602 m/s。相同总流量与 Flat + Zero_out_perimeter 规则不保证逐点入口速度剖面不变。原注入代码 FiniteSizeAdmission.check 取出生平动速度=当地 FEM 速度、MB 角速度=涡量/2。因此本研究冻结出生时间/位置/半径/取向/ID，沿用同一注入规则，而非另将场依赖的出生速度/角速度冻结。ID440 的 M0/M2中出生位置逐坐标相同，但速度为 0.1287413820/0.2817099407 m/s；其路径变化不能全部归因于下游零速结构消失。两次重放的 injection_admission 源码哈希已核对相同。

此项已集成 M2_MEDIUM_FINE_FIELD_REVIEW_ZH.md/JSON 和总报告；关键速度图标题改为 Frozen birth marks; field-derived birth velocities; pairs OFF。未改任何注入或计算代码。M3 中档仍在计算，继续监控；以实时文件为准。


### 2026-09-30：入口离散剖面可视化与动态报告

新增 `figures/M2_inlet_discrete_profiles.png/.pdf`（300 dpi），使用 M0 / M2 中档 / M2 细档原始入口三角形的 P1 面积平均速率和实际连通关系；统一色标、逐面常色，无科学数据平滑。入口均为平面，实际节点速度共线且非负，因此节点速率均值即 P1 面均值；三档面积积分均核对为 18 mL/min。节点/三角形分别为 57/90、113/202、135/246。图已人工查看，解释为何同一总流量与 Flat/Zero_out_perimeter 规则仍产生不同入口离散剖面及出生采样速度。数值与入口 surface/flow 数组 SHA 归档于 `04_field_audit/M2_MEDIUM_FINE_FIELD_REVIEW.json`。

`render_report.py` 现按实际细档重放完成状态生成 M2 中细档结果段落，出生参数/场依赖初速度的术语已明确。本地报告生成通过，63 个 HTML 链接和 36 份图文件哈希核对通过。45 份脚本及内容寻址归档于 00:26:47 UTC 在服务器核对通过，snapshot SHA `1a7d70d1b2f5cb9c77107cf849fc37f43c8957f29edc330467053ea48b4a1c2d`；report_watch 已重启并成功更新。只修改审计/报告代码，CFD、粒子物理、阈值及原始结果均未修改。

00:30:43 UTC 监控：M3 中档第 72 步；第 70 步检查后连续稳态计数 4/5。仍运行中，尚未宣称 qualified flow。M2 细档重放继续等待共享重计算锁，M3 细档与固定原始边界补充对照保留自动队列。用户要求持续监控，不因后台仍在计算就结束对话。


### 2026-09-30 00:54 UTC：M3 中档完成、M2 细档轨迹与完整 gate

M3_medium_frame180 在第 81 步正常结束，CFD 6720.931323051453 s，独立后处理 PASS、field audit 完成。flow SHA `0783778226a58a9edd7d8575ec742eb938748988a347e3daa17443d8b00150c5`，mesh SHA `7a50be05b2b4f2aa4d4d5d89b45bed867bcdb15a2679080e90501b1c26e98d46`。入口 18 mL/min，出口比例 .334713017247/.332768421965/.332518560788，质量残差 1.76465e-16，壁速严格 0；WSS mean/P95/max 3.50397593/6.53037897/37.63848851 Pa，无严格零速 tetra。1059 条冻结采样记录中 38 条域外。原 388769 中心仍在域内，速度 .030104155584 m/s、div −780.395167625 s^-1、eta .3177322694，壁距 11.125341 µm；原 138584 中心域外，u/G 为 NaN，绝不补零或投影。固定高 eta 球区中 eta≥.5 体积 2.2472362145e-12 m³，但几何变化混杂仍存在。新 M3 两个截面已生成，零单元截面已查看。

M2_fine_frame180_dt0.00050000 已完成：46 EXITED，0 solver failure/censored/supported stationary；O1/O2/O3=16/18/12，307.3815131187439 s。与 M2 中档同冻结出生参数、同 dt、同注入规则与粒子源码，46/46 两次均出流，但完整 medium/fine gate FAIL；逐条仅 ID1565 通过（1/46）。未通过：速度历史 44、最小间隙 39、停留时间 27、末点 26、出口 6、路长 3。改变出口 ID9/168/301/542/708/1463；场分量也 FAIL，不能放行 Phase3。M1 与 M2 已完成的空间 gate 均 FAIL；M3 与固定 PLC 对照尚待完成。

细档 slow-fast 为 7 次：ID301 六次、ID636 一次；中档 16 次，其中 ID301 三次。ID301 中→细 O2→O1、停留 5.9898620254→3.8881474698 s，细档事件窗口 eta 最大值 .05844–.09064，事件位于距原 388769 中心 21.11–23.54 mm，不能都归因于原高 eta 热点；其原选组为 high_near_wall_time。ID636 两档均 O1且保留一次事件，细档 onset .831 s，比中档提前 .0235 s、位置差 62.679785 µm，距原 388769 中心 .176087 mm，速度历史相对误差 .30874617，仍不收敛。窗内峰速可能属于事件前入射高速，新增 JSON 区分 pre/post onset 峰值。没有修改接触或速度科学数组。

新增 `M2_MEDIUM_FINE_TRAJECTORY_REVIEW_ZH.md`、`05_particle_replay/M2_MEDIUM_FINE_TRAJECTORY_REVIEW.json`，原始 ID301/636 速度/间隙/eta 图 `figures/M2_remaining_slow_fast.png/.pdf`，gate/停留/出口转移图 `figures/M2_transport_convergence.png/.pdf`。两图 300 dpi，均已查看；总报告与 OPEN_RESULTS 已纳入。新增 `scripts/report_m2_trajectories.py`，只审计/报告。当前 46 份脚本 snapshot SHA `96481260a9ba82403c2c41c147c17b9781b12ef0cb0ae2eeba549fd00d6dbdaf`，服务器 00:53:45 UTC 全部当前/归档源码核对 PASS。最新本地报告 00:52:15 UTC，79 HTML 链接/46 图哈希核对 PASS。本地 COLLECTION_VERIFIED 已有 5 个候选 CFD +5 个已完成 replay。

当前队列：M3_fine_h075_frame180 于约 00:43:48 UTC 开始 CFD，00:53:58仍在第 1 步非线性2后的矩阵重建；内存曾达42.03 GiB并回收页缓存，OOM=0、无发散。后续 M3 中档/细档 replay 及 fixedplc_validation_followup 等待共享重计算锁/原 campaign 完成。不可重复启动。用户要求持续监控，不因后台计算仍在运行就结束当前对话。


### 2026-09-30 01:37 UTC：M2 同坐标交叉采样已完成

新增只读审计 `scripts/audit_m2_path_field.py`，对全部 46 条中/细档轨迹的共同存活时间区间取原时间戳并集，对记录位置作分段线性重建，交叉采样两档冻结 FEM。每档 256586 个坐标查询；两档域外点均为 0，46 条完整有效。局部单进程分两次加载流场，实际采样耗时 756.1017573 / 816.1113656 s，现已结束，无本地重计算进程需恢复。原始端点 FEM 速度重现最大分量误差 2.22044604925e-16 m/s，域内分类差异 0；向量分解恒等式残差最大5.55111512313e-17 m/s。没有积分新轨迹、移动采样点、平滑或补零。

分解：`u_f(x_f)-u_m(x_m)=[u_f(x_m)-u_m(x_m)]+[u_f(x_f)-u_f(x_m)]`，另用细档路径/中档场作基准复核。全46条时间加权向量 RMS 中位数：总差40.9619942 mm/s，同中档路径场差0.230766220 mm/s，同细档路径场差0.241994398 mm/s。ID301 总差58.972655 mm/s，对应同位置场差0.147123/0.188431；ID636为4.183906与0.282617/0.306076；ID440为23.631719与0.277963/0.303415。**JSON原始单位为m/s，报告与图乘1000显示mm/s。** 路径位置项幅度大，但这是依赖基准且可相消的运动学恒等式，不是可相加的因果比例；路径分离自身仍为响应未收敛表现，不能豁免原 gate，也未分离 G、出生采样、接触及 dt 的贡献。

产物：`04_field_audit/M2_PATH_FIELD_DECOMPOSITION/` 中的 QUERY_MANIFEST、两档 samples、46份 decomposition NPZ、DECOMPOSITION.json；`M2_PATH_FIELD_DECOMPOSITION_ZH.md`；`figures/M2_path_field_decomposition.png/.pdf`（300dpi，已查看）。总报告与 OPEN_RESULTS 已集成，上传服务器后01:37:06逐项核验流场/轨迹/查询/结果哈希 PASS，服务器记录 `REMOTE_TRANSFER_VERIFIED.json`。

M2 轨迹专项报告另补实际近壁模型解释：`assemble_isolated` 的剪切驱动依赖 G，瞬时阻力/接触约束求速度；P1 G逐单元常值、最近壁面特征与接触活动集变化都是存在的敏感机制，目前没有给每个事件做独立因果分摊。冻结下界 `max(2e-9 m, .001*a_ref)` 对ID301/636均为2nm；这是连续介质模型交接下界，不是原子接触/糖萼/壳厚测量。源码与执行identity哈希核对一致，未改任何粒子物理。

当前47份脚本（含归档）服务器核对PASS，snapshot SHA `4ee68ad51bab13edd58290d3c4af109e959579fb8740649b8f4ab225cb3b65c3`，audit源码SHA `ad65a908bce46d3b2f434ee4d388f412476864ce1f71f909b6b873bbe76fe369`。最新本地报告01:35:32 UTC，83个HTML链接/48份图文件哈希核对PASS。report_watch已重启。M3细档01:37:30进度第8步非线性8，内存38.75GiB，无OOM/发散；后续重放与fixedPLC队列不变。继续保持对话活跃并监控。


## 持续监控更新：2026-09-30 01:56 UTC

- 用户要求计算/开发停止前保持当前会话监控；不要因后台进程仍在运行而提前结束。M3_fine_h075_frame180 正在第 12 步，最近保存 step10，连续稳态检查尚为 0/5；无 PETSc 发散、非线性失败或 OOM。后续 M3 重放及 fixedPLC 中/细档仍按隔离重计算锁排队。
- 已完成沿 M2 中、细档同一记录坐标的速度与逐单元 G 对照，共 46×2=92 个比较。复用已校验的 canonical 单元编号，直接读取原 FIELD_ARRAYS，不移动位置、不改场、不重跑轨迹。
- 时间加权相对 RMS 中位数：u 为 0.133066% / 0.140853%，G 为 4.496538% / 4.791847%，分别对应中档/细档路径；G 逐轨迹最大值 21.814579% / 18.493759%。ID440 的 G 为 8.000651% / 7.853748%，ID636 为 7.955350% / 10.434798%。这些沿途统计与原固定区域 gate 的权重和位置不同，不能替换固定区域 G 门槛或解释为力误差百分比。
- 新文件 M2_PATH_GRADIENT_REVIEW_ZH.md、04_field_audit/M2_PATH_FIELD_DECOMPOSITION/PATH_GRADIENTS.json、figures/M2_path_velocity_gradient.png/pdf；已整合进 FINAL_VALIDATION_REPORT_ZH.md 和 OPEN_RESULTS.html。原 46 条轨迹判定保持 M2 FAIL、1/46 PASS；Phase3/4/5 未启动。
- 本地报告校验：87 个 HTML 引用存在、50 个图文件 SHA 正确。48 个当前脚本及执行归档已在服务器逐文件验证；当前快照 SHA256 8e7d3be533d6b1d455264f5f9b621869148886bdb94281a6b6ce57cc1ab3adae。仅重启 report_watch 以加载报告变更，未改动活动求解器或轨迹物理。


## 速度输入时刻审查：2026-09-30 02:04 UTC

- 新增 scripts/audit_event_input_timing.py，只读分析 M2 中/细档全部 23 次筛查事件。v[i] 来自 x[i−1] 的求速输入，作用于此前接受区间；逐状态位移恒等式 max residual 3.49e−17 m。
- M2 细档 7 个 onset 的前后输入均在同一个 tetra，G 没有跨单元跳变；h/a 最小约 2.91986，当前近壁修正为零；v 与起点 FEM u 最大差 2.74e−18 m/s，Δu=GΔx 最大残差 3.81e−17 m/s。onset 是筛查阈值触发，不代表瞬时物理突变；不能据此认证场或步长。
- ID636 窗口最大向量变化位于事件前减速段：t=0.821602382238 s，165.84671→4.381839 mm/s，向量变化162.58314 mm/s；两个实际求速输入 FEM 向量变化162.58428 mm/s。前一步0.5 ms移动82.92335 µm，壁间隙80.629019→0.346824 µm，输入单元改变；当前步经几何截短为0.102382 ms。保持全部接受状态，不把几何安全当积分误差控制。
- 输出 M2_EVENT_INPUT_TIMING_ZH.md、05_particle_replay/EVENT_INPUT_TIMING.json、figures/M2_event_input_timing.png/pdf，整合到总报告/HTML。已目视审查图，本地91个HTML引用和52个图SHA通过。
- 当前49脚本快照 SHA256 004d4a0b5c3f2db745f10aa1ec9ce1788d8a33cb833a2315773efc6d3261141b。服务器部署验证及 report_watch 单独重启；原求解器、物理和门槛均未修改。M3细档当前第15步，无发散/OOM；最近区间平均CPU7.6829核（配额7.68），GPU利用率随阶段波动，显存约7419MiB。


## CUDA 近壁求速分解及检查点：2026-09-30 02:27 UTC

- M3_fine_h075_frame180 已保存 step20（t=0.2 s）：epsilon_mass=9.536743164336735e−7，no-slip PASS；相邻10→20步 E_u=0.0010067926167980944、E_Q=0.000531091094777434，尚高于稳态阈值，连续通过0/5。当前step22，继续原求解。无发散/OOM。
- 额外完成轻量只读 CUDA 求速分解：scripts/audit_m2_velocity_components.py 在 ID301/636 的7791/3374个原接受区间起点，重新调用原GPU wall/resistance/contact内核，没有重新积分、改场或改物理。执行时所有particle源文件/边界/轨迹/G数组SHA校验，0个CPU fallback。重现平动误差max 7.209740646747396e−18 / 1.371419323662774e−18 m/s，原gap逐点完全一致。服务器用时8.604s，host maxRSS1671.992MiB，CUDA张量峰值4220.733MiB；与CFD短时并行前检查了>2GiB剩余cgroup内存，没有触发OOM。
- 在同一完整R下代数分解 V=R^-1 D Vfree + R^-1 b_wall_shear + (Vcontact−Vhydro)。ID301/636正接触乘子区间1399/10，持续0.6959310935288676/0.004397617761977113s；最大接触平动修正0.0002859128641286205/0.005678746352278384 mm/s；最大壁剪切项1.5611152937517072/1.1230889408403073 mm/s。小的瞬时接触修正仍能决定可行位置与长期停留，不能据此删接触或判其不重要。此为当前模型的恒等分解，不是因果比例或物理认证。
- 输出目录04_field_audit/M2_VELOCITY_COMPONENTS，元数据文件特意名为COMPONENTS.json（不能叫SUMMARY.json，避免被field候选glob误读）。已执行源SHA09c45e5c64864469177356d167b13492dc3b712c43e7cba7b5b13c4e0a2bfc09归档；该次输出最初写SUMMARY后仅重命名为COMPONENTS，内容SHA2e77639093e15089f7694897d87b8762475369c79124e2f000658536ea358d25不变。当前helper已改成直接写COMPONENTS；不要重跑覆盖原目录。
- M2_VELOCITY_COMPONENT_REVIEW_ZH.md、REVIEW.json、figures/M2_velocity_components.png/pdf已生成、目视审查并整合总报告/HTML。投影采用当地输入法向/自由流切向，基底变化也可导致投影折点；法向symlog线性区±0.001mm/s仅为显示。平动/角速度闭合残差分开记录，最大1.76162e−18m/s、7.08265e−13rad/s。
- 当前51脚本快照SHAbe5892a78cc5b0528b94f40355ae6b581db452815329559917f11cd012a5e6a9，服务器当前/归档逐文件验证PASS。95个HTML引用、54个图SHA校验通过。仅重启report_watch；原CFD/replay队列和物理代码未变。

2026-09-30 02:40 UTC：brain_report.md 已新增第15节本轮验证工程索引（26个新增入口存在性通过），并将原只读核对范围明确限定为2026-09-29，保留旧路径/生产身份。服务器最新报告02:27:23 UTC通过95链接/54图SHA校验。M3细档当前step27-5，最近保存仍step20；继续监控，不结束活动会话。

2026-09-30 02:50 UTC：M3_fine_h075_frame180 已保存step30（t=0.3s）：mass=9.313227840223519e−10，no-slip PASS；20→30区间E_u=6.46539527783536e−6已低于1e−5，E_Q=1.163560846853065e−6仍高于1e−6，因此连续稳态仍0/5，不能放宽阈值。当前step31，无发散/OOM。监控副本已同步。

2026-09-30 03:10 UTC：M3_fine_h075_frame180 保存step40(t=0.4s)，mass=9.09501627114305e−13、no-slip PASS；30→40区间E_u=1.3411969990950936e−7、E_Q=3.767887099517846e−9，首次在线稳态通过，连续1/5。当前step41-2，无发散/OOM；仍需后续连续通过和独立postprocess/audit，不能提前标记qualified。监控副本已同步。

2026-09-30 03:27 UTC：M3_fine_h075_frame180 保存step50(t=0.5s)，mass=1.0587911840678754e−15、no-slip PASS；40→50区间E_u=3.162921143997744e−9、E_Q=1.6586360945117294e−11，连续在线稳态2/5。当前step52-2，无发散/OOM，平均CPU约7.68核。监控副本同步，后续验证继续排队；仍不提前认定qualified。

2026-09-30 03:38 UTC：M3_fine_h075_frame180 保存step60(t=0.6s)，质量不平衡浮点计算结果0、no-slip PASS；50→60区间E_u=6.835025882349445e−11、E_Q=8.170338637057106e−14，连续在线稳态3/5。当前step61-2，无发散/OOM；仍需两个连续通过区间与独立后处理/局部场审计。

2026-09-30 03:51 UTC：M3_fine_h075_frame180 保存step70(t约0.7s)，mass=5.293955920339377e−16、no-slip PASS；60→70区间E_u=1.2263649719662043e−12、E_Q=2.117582368135751e−15，连续在线稳态4/5。当前step71-2，无发散/OOM；待最后一个连续通过区间，之后仍需postprocess/audit。后续M3重放、fixedPLC中细档队列仍正常等待。


## M3 细档完成与首轮重放：2026-09-30 04:22 UTC

- M3_fine_h075_frame180 CFD 正常完成至 step81/t=0.81s，耗时12013.207763s；后处理与场审计均已完成。flow SHA256 77e8dd5787f716d79f46a830ccd21ce04d4aed4bb75867eb8959c15ead6ffd3c，mesh f8ca5b3350eb69d8bdd617cbf35e1e7b65f8bb0c7acd6b08ebb1b544cfb0e07d。Q=18 mL/min，质量不平衡5.293956e−16，无滑移零速。全局通过不构成局部准确性认证。
- M3 中细档固定点共1059记录，公共1015、非公共44（两档都域外36，仅中档域外2，仅细档域外6）。ALL u/G/最近壁面WSS相对L2为1.832158%/23.677805%/22.988770%；原高η区为3.225192%/57.080631%/43.753141%。场部分FAIL，完整输运gate仍PENDING。没有投影、剔除、补零或改变门槛。
- 原388769中心中→细速度30.104156→42.336765 mm/s，div −780.395168→−199.679852 s⁻¹，η 0.317732→0.0590803。高η球区η≥0.5体积2.2472362e−12→1.28341895e−12 m³，P99 abs(div)582.120015→448.428578。严格零速tet均0，但原138584中心两档均域外，不能声称同点物理修复；仍需固定原边界对照。
- M3_medium重放46条：45 EXITED、1 RIGHT_CENSORED（ID487）、0 SOLVER_FAILURE；O1/O2/O3=14/20/11，耗时472.602775s，共15次slow→fast筛查事件。440出O1/停留1.176582s；1505出O3/2.832163s；636出O3/3.275914s、1事件；301出O3/4.699441s、4事件。ID487观察年龄10.85101168s，不是完整出流停留时间。仅终止成功仍不能证明网格收敛。
- 新增scripts/report_m3_fields.py，对比表由原NPZ独立复算校验，生成04_field_audit/M3_FIELD_REVIEW.json、M3_MEDIUM_FINE_FIELD_REVIEW_ZH.md及figures/M3_field_convergence.png/pdf，已目视检查并整合总报告/HTML。脚本52份，snapshot SHA0148fca16325b6dd342d66050d8dcfd79fd796190012ceb126f43293c19987c2；服务器当前与归档逐项PASS。本地报告04:19:58 UTC，102链接/58图SHA核对PASS。
- 已收集并验证6套新CFD、6套46条replay。resume_campaign正常EXITED。fixedplc_validation_followup先获得重计算锁，当前M1_fixedplc_medium CFD step4，无发散/OOM；replay_followup等待M3细档重放。只重启report_watch加载报告代码，不改活动CFD/replay物理或队列顺序。用户要求保持监控，当前会话仍应继续，不因后台运行结束对话。

2026-09-30 04:30 UTC：M1_fixedplc_medium 首个保存step10/t=0.1s，Q=18mL/min、wall no-slip PASS；mass=0.0009765625000005532，尚高于1e−6；没有可比较的相邻保存区间，连续稳态0/5。CPU最近区间平均7.682核/配额7.68，内存约17.9GiB，无发散/OOM。M3细档replay继续排队。服务器报告04:22:09 UTC通过102链接/58图哈希校验。继续监控，不结束会话。

2026-09-30 04:40 UTC：M1_fixedplc_medium 保存step20/t=0.2s，mass=9.536743166101389e−7，无滑移PASS；10→20区间E_u=0.0005978380476882767、E_Q=0.0005312727440190861，均未通过稳态阈值，连续0/5。当前step21，无发散/OOM；M3细档replay仍按锁排队。六批速度/终止图目视检查可读，无科学数据或显示修改。

2026-09-30 04:47 UTC：M1_fixedplc_medium 保存step30/t=0.3s，mass=9.313229604875494e−10、no-slip PASS；20→30区间E_u=6.992473098982338e−6已通过1e−5，E_Q=1.1735879781046173e−6仍高于1e−6，因此连续稳态0/5。当前step31，无发散/OOM，继续监控与排队验证。

2026-09-30 04:52 UTC：M1_fixedplc_medium 保存step40/t=0.4s，mass=9.087957663249265e−13、no-slip PASS；30→40区间E_u=1.4256495928721642e−7、E_Q=3.784891550631772e−9，首次在线稳态PASS，连续1/5。当前step42，无发散/OOM；仍需后续4个通过区间和独立后处理/场审计。M3细档replay仍排队，保持会话监控。

2026-09-30 04:56 UTC：M1_fixedplc_medium 保存step50/t=0.5s，mass=1.4117215787571673e−15、no-slip PASS；40→50区间E_u=3.1866196118312683e−9、E_Q=1.4976910112735453e−11，连续在线稳态2/5。求解继续，无发散/OOM；尚未做独立postprocess，不能提前标qualified。

2026-09-30 04:59 UTC：M1_fixedplc_medium 保存step60/t=0.6s，mass=5.293955920339377e−16、no-slip PASS；50→60区间E_u=6.823587769328541e−11、E_Q=1.4778960277614096e−14，连续在线稳态3/5。当前step62，无发散/OOM；尚需2个通过区间与独立后处理/审计。

2026-09-30 05:01 UTC：M1_fixedplc_medium 保存step70/t约0.7s，mass=5.293955920339377e−16、no-slip PASS；60→70区间E_u=1.4109426047139342e−12、E_Q=1.6323030754379747e−15，连续在线稳态4/5。尚需最后一个通过区间和独立postprocess/场审计；M3细档replay继续排队，无发散/OOM。


## 固定原始边界中档 CFD 完成：2026-09-30 05:12 UTC

- M1_fixedplc_medium 正常完成至step81，CFD耗时3046.8507647514343s；独立后处理PASS，GPU场审计35.050454s完成。flow SHA34dc701a2abddba75e3aaec9b957e19606746333a61153eefbc8c601ac293d2d；mesh SHA bd55dfecdbcacbeb030bdf4c853c98130a85bd92933c6e7905fe2c61df6beb7f；flow_arrays SHA2ce9f31179fa7e1db779668ef6ec8f7b5421901c4c77d96ea94084c13b0927d3。
- Q18mL/min，O1/O2/O3分流0.333352798730/0.333241934278/0.333405266992，mass5.293956e−16，wall u=0；p范围0.4012090844–2721.25341995Pa，umax0.555804656297m/s，WSS面积均值3.07307300817Pa、P955.71529748269、max11.28195065961，零WSS面0。末区间Δu/p/WSS=4.346787e−15/5.012235e−16/6.817329e−15。
- 全762096tet没有严格零速单元，1059固定记录全部域内。边界全部5角色原三角形坐标集合保留；相对M0采样壁距最大数值差8.24798e−19m。原388769中心：cid640459，u=0.0673137627506m/s，div−835.198383551s⁻¹，η0.583408546，wallgap47.398452042µm。原138584中心（记录778/782）：cid737895，u=0.0174271102756m/s，div−182.600963228，η0.274616758747，wallgap26.208933697µm。原零速位置仍域内且恢复非零速度，可排除移出域这一解释；但高散度点仍明显异常，不构成局部精度通过。
- 高η固定球区η≥0.5体积2.922481573986097e−12m³，P99 absdiv377.593591318、rQ0.059651832984、χ0.737539194676；原零速球区体积4.02524822830543e−12、P99absdiv254.270851183、rQ0.09242695145、χ1.13806810972。M0→中档ALL u/G/WSS差6.432275%/25.432941%/22.784283%，仅基线变化诊断，不是中细档收敛。
- 新增scripts/report_fixedplc_fields.py、FIXED_PLC_FIELD_REVIEW_ZH.md、04_field_audit/FIXED_PLC_FIELD_REVIEW.json、figures/fixed_plc_target_fields.png/pdf，已目视核查与集成总报告/HTML。M0零速处η仍未定义，图中不填0，文字放在完整x范围内避免压轴。
- 当前53脚本，snapshot SHA356f0f329d184156fe7a78ce31f85af7cfaed573bc4c79dae694a93a538b7bf3；服务器当前/归档逐项PASS。仅重启report_watch，本地05:10:33报告106链接/60图SHA通过。收集器已验证7套新CFD、6套replay；M3细档需下一次收集。
- M3_fine_h075_frame180重放已完成：45EXITED、1RIGHT_CENSORED、0FAIL，耗时392.72720337s；正等分析与中细档比较。fixedPLCmedium重放将接入，之后才fixedPLCfine CFD/replay。没有选正式flow或开展Phase3/4/5。继续保持活动会话监控。


## M3 完整输运对照失败：2026-09-30 05:19 UTC

- M3细档重放已收集、哈希验证与分析完成：45EXITED、1RIGHT_CENSORED(ID123)、0FAIL；O1/O2/O3=17/18/10，耗时392.72720337s。中档观察截止为ID487，细档487已O3出流/3.60706760332s；细档123观察年龄11.137332773s仅为停留时间下界。不能把两档45出流/1截止解释为同样结果。
- 中细档完整gate FAIL、0/46逐条PASS。各指标失败数：速度历史46、停留/观察持续时间43、endpoint39、mingap38、出口/删失不一致12、pathlength8、终止分类2。44条在两档均出流，其中10条改出口，ID[5,7,9,168,301,542,636,878,1377,1505]；终止分类改变ID[123,487]。含删失4×4转换矩阵（O1/O2/O3/删失，中档行细档列）[[13,0,1,0],[1,16,3,0],[3,2,5,1],[0,0,1,0]]。
- 事件中档15次，ID220:1/301:4/487:8/636:1/1463:1；细档10次，ID301:2/1377:2/1463:3/1522:3。细档636事件消失，但出口变O2/停留2.35351218374s/minGap30.19470578µm；301变O1/1.89154245423s。440细档O1/1.01080447102s，1505O1/0.986327374014s；所有原数据保留，不把事件消失当同路径问题修复。
- 细档run identity SHAb43731e3d159ef16dea75ea8dd19e80f5206ea251c64e417477d7af673153dac；analysis summary d09737dbc5714e4b5189b3fac8b82049ef1fca17c10cb58be5913cc4dab0b863。当前本地已验证7套新CFD、7套replay。
- 新增scripts/report_m3_trajectories.py、M3_MEDIUM_FINE_TRAJECTORY_REVIEW_ZH.md、05_particle_replay/M3_MEDIUM_FINE_TRAJECTORY_REVIEW.json，以及figures/M3_transport_convergence.png/pdf、M3_trajectory_sensitivity.png/pdf；已目视检查并集成总报告/HTML。散点三角形明确删失时间为下界；区分两档均出流换出口与出流/删失转换；保留速度区间语义、不按事件对齐。
- 当前54脚本snapshot SHAeb8e8e0031e141d6daee4b2b2e2f342a003ad86ba2081dbe29ae2a46b3d4c49a，服务器当前/执行归档逐项PASS。最新本地报告05:18:05UTC，116链接/66图哈希PASS；report_watch重启后需复核服务器最新生成结果。
- replay_followup已正常EXITED（M3两档均结束）。fixedplc_validation_followup正在运行M1_fixedplc_medium 46条重放，05:19时23/46完成；随后还需M1_fixedplc_fine CFD/replay。report_watch RUNNING。M1/M2/M3均FAIL，fixedPLC完整pair仍PENDING；Phase3/4/5未运行。保持会话活跃持续监控，不因队列在后台而结束。


## 固定原始边界中档重放完成：2026-09-30 05:28 UTC

- M1_fixedplc_medium_dt0.00050000完成全部46：40EXITED、6RIGHT_CENSORED、0SOLVER_FAILURE、0SUPPORTED_STATIONARY；O1/O2/O3=16/12/12；耗时887.287274599s；18次slow→fast事件。观察截止ID[9,393,1017,1329,1377,1463]，保留年龄下界，不视为完整停留时间或数值失败。
- 关键轨迹：301 O2/停留10.9888642816s/1事件（绝对出流时刻11.2381090792，仅比11.240815660观察终点早约2.707ms）；393观察截止年龄10.9249286671s、0事件；440 O3/1.42209430397s/minGap33.3457352711µm/0事件；636 O3/1.78806352905s/minGap1.58610944152µm/0事件；1189 O3/6.83404761950s/minGap约2nm/0事件；1329观察截止年龄10.1881442359s、0事件；1505 O3/1.50814242277s/minGap41.9355783854µm/0事件。
- ID1017观察截止时速度约0.26m/s、gap347.898µm，不能把所有观察截止都解释为近壁静止或数值停滞；其他删失也要保留原终止语义。
- fixedplc_validation_followup继续下一项M1_fixedplc_fine CFD；没有重复启动或改变物理参数。收集器正在同步第8套replay，完成后需update_comparisons/assess/render。服务器05:19:00报告116链接/66图SHA核对PASS；当前54脚本快照不变。保持会话监控，不提前结束。

2026-09-30 05:34 UTC：第8套replay已完整收集、哈希验证。最新本地总报告05:32:09UTC，服务器05:28:20UTC，均通过119链接/68图SHA检查。发现collect_completed元数据同步能回写旧RESUME_REPORT，已从独立BRAVA_CONTEXT_HANDOFF恢复缺失的831字符（仅缺05:28中档重放完成段，旧内容确认为新内容前缀）；collector新增RESUME_REPORT_ZH.md排除项，恢复说明继续本地维护后上传。未改科学数据。当前54脚本snapshot SHA628537c4e217a88cb92db424e884f7eba54729257c2b0fa16a09796f1e235af0，服务器当前/归档逐项PASS；只改collector，无须重启report_watch。M1_fixedplc_fine已启动约6分钟，PETSc日志停留各rank矩阵信息、尚无首条NS迭代；CPU平均7.68核持续工作，内存约20.3GiB且缓慢上升，无发散/OOM，继续观察矩阵建立阶段。


## 固定边界出生输入与求速语义：2026-09-30 05:49 UTC

- 新增逐条原始首状态块核对，首块 SHA 与 MANIFEST 一致；46 条出生位置、半径、取向、时间逐项精确一致。M0→M1_fixedplc_medium 最大平动速度向量差5.72326868918e−15m/s，最大角速度向量差39.9915652526rad/s（ID175）；两种单位分开统计。出生速度按各场原规则派生，不能称完整数值初态逐项一致。
- 同时核查已执行 gpu_engine/gpu_batch/gpu_resistance 的源 SHA：球形过阻尼模型每步按当前位置的 u、半涡量、G 与壁几何瞬时求速并覆盖记录 v/omega，旧出生 omega 不作为惯性初值进入下一位置更新。出生角速度差不是后续路径分离的独立原因；没有改变物理、重新积分或分配因果百分比。
- report_fixedplc_fields.py 已生成 figures/fixed_plc_birth_inputs.png/pdf，并扩展 FIXED_PLC_FIELD_REVIEW_ZH.md 和 04_field_audit/FIXED_PLC_FIELD_REVIEW.json；新图已目视核查。细档重放关闭后会按相同规则自动纳入。
- 最新本地总报告05:48:03 UTC，122个HTML链接/70图SHA核对PASS。前一服务器报告05:47:12 UTC也通过122/70；最新语义补充已部署，report_watch仅重启自身待刷新。当前54脚本snapshot SHA cd9b11491b87042d6ad978e7b0a85b4ce5feb7c7f05a1cbb635f3f3921f6c668，服务器当前及归档逐项PASS。
- M1_fixedplc_fine 初始矩阵建立耗时较长，05:45已到step1-10，05:46进入step2；05:48最后NS2-1，线性残差5.603e−11，内存约21.2GiB，无发散或OOM。尚无step10检查点，在线稳态0/5；固定边界完整pair仍PENDING。用户要求保持当前会话监控，不能因后台运行而结束对话。

2026-09-30 06:29 UTC：M1_fixedplc_fine 首个保存检查点 step10/t=0.1s，epsilon_mass=0.0009765624999991414，wall no-slip PASS；质量不平衡仍高于1e−6，且尚无相邻保存区间，因此连续稳态0/5。当前已到step11-3，无发散/病态矩阵警告/非线性超限/OOM，cgroup内存约22.65GiB。前期多次矩阵建立较慢，8个MPI rank运行态、平均CPU接近7.68配额；实际PETSc设置为CUDA Krylov + CPU子域LU，不能仅凭GPU利用率判断没有计算，也不把CPU占用等同并行效率证明。保持原参数继续；固定边界完整pair仍PENDING，Phase3/4/5未启动。服务器最新报告05:49:40UTC，122链接/70图SHA核验PASS，包含出生omega不是惯性初值的说明。


## CPU/GPU 求解成本审查：2026-09-30 06:42 UTC

- 新增 scripts/report_solver_costs.py，只读解析已关闭9个PETSc执行段（7套新CFD，M2中档拆原暂停段/正式续算并单列重启试验）。SOLVER_RESOURCE_REVIEW_ZH.md、SOLVER_RESOURCE_REVIEW.json、figures/solver_resource_costs.png/pdf已生成并目视核查，加入总报告、HTML和brain_report索引。活动M1_fixedplc_fine尚无最终profile，明确排除。
- 解析绑定原profile/execution/PETSC_OPTIONS哈希，校验打印表头、stage时间求和、事件含义及GPU工作计数。MatMult事件GPU %F均100，MatSolve与MatLUFactorNum均0；GPU %F不是利用率、加速比或耗时百分比。事件时间取MPI rank最大累计值且可能嵌套，不相加。原-log_view_gpu_time计时开销保持并明确说明。
- fixedPLC中档已关闭profile：MatLUFactorNum 19次/max733.08s/rank max:min5.1；MatSolve max1579.5s；MatMult max263.54s；PETSc总时间max3040s，外部执行3046.850765s。M3细档LU max4201.6s/max:min2.7，MatSolve max6203.2s。此为当前混合CPU/GPU路径成本，不是局部物理精度或并行效率认证。活动细档日志显示子域LU，最终瓶颈仍待其自身profile，未改求解参数或启动性能试验。
- 当前55脚本snapshot SHA ad0dca31cab6a7b000a8975b23b049e9cb604da83e4af1b026bad6bfa53a003e，服务器当前及归档逐项验证后仅重启report_watch。最新本地总报告06:39:44UTC，126链接/72图SHA及9份profile输入哈希PASS；等待服务器刷新后复核。
- M1_fixedplc_fine 06:40已到step13-7；最近保存step10，mass约9.765625e−4、no-slip PASS、在线稳态0/5。内存约22.7GiB，无发散/非线性失败/OOM；继续活动会话监控，原物理、轨迹和生产身份不变。

2026-09-30 06:43 UTC：服务器06:41:53UTC版报告已通过126链接/72图SHA校验。当前M1_fixedplc_fine已进入step15，健康检查无新增异常，内存约22.7GiB；最近保存仍step10，继续监控。

2026-09-30 06:54 UTC：性能报告措辞改为按实际活动段描述，pending识别改读RESUME_INDEX当前续算段，避免后续恢复时被旧根执行记录误判关闭。无科学输入/结果变化。55脚本当前snapshot SHA490765078a61df9ad8d0a109c4afcabe5365f0224399c59a874e00e4d817e2a4，服务器当前及归档逐项PASS；仅重启report_watch。最新本地报告06:52:49UTC，126链接/72图SHA PASS。细档最新完成NS17-9，最近保存step10，内存约23.1GiB，无健康检查异常；继续监控。

2026-09-30 06:59 UTC：M1_fixedplc_fine保存step20/t=0.2s，epsilon_mass=9.536743162572084e−7已低于1e−6，wall no-slip PASS；10→20步E_u=0.0006045619708760976、E_Q=0.0005327901592196962，仍高于稳态阈值1e−5/1e−6，连续通过0/5。当前step21-2，内存约22.96GiB，无发散/病态矩阵警告/非线性超限/OOM。服务器06:53:56UTC版报告及当前性能helper源码身份已核对，126链接/72图SHA PASS。55脚本snapshot仍490765078a61df9ad8d0a109c4afcabe5365f0224399c59a874e00e4d817e2a4，保持会话监控，未进入Phase3/4/5。


## 历史失败机制与报告链接复核：2026-09-30 07:22 UTC

- 更新P1_FAILURE_MECHANISM_ZH.md：移除已过时的M2中细档/重放尚待完成措辞，改为实际M2完整gate FAIL、1/46通过、6条改出口，并补固定原边界中档同坐标结果。仅更新说明，未改科学数组或门槛。
- 重新读取现有population_v1/analysis/failure_speed_review_20260929/CAUSE_ANALYSIS_ZH.md并核对安全代码。两条历史失败从578387渐近逼近共享内部零速面，G非零特征值−443.471477s⁻¹、约2.255ms离散衰减尺度；末点面距4.1–4.8e−16m，最小重心坐标−4.37e−12/−5.13e−12在约5.38e−12包含容差内。按较小canonical ID归入138584，没有宏观穿面/穿壁证据，不能把失败归因于显式大步穿过面。相邻单元末点求值也仅2.2–2.5e−13m/s。
- 精确区分三项：历史两例触发ROUNDOFF_SCALE_STAGNATION_64_ACCEPTED_SUBSTEPS；CPU另作32实际状态支持复核；代码还存在独立EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS分支，但不能用它改写历史终止。两例严格零速记录82/71个，检查调用时序使其不必正好64；末gap13.630/8.026µm、接触秩0，保留数值失败而非吸附/支持静止。
- 发现服务器缺少本地既有failure_speed_review_20260929和acceleration_diagnosis诊断目录，已原样复制这两个此前不存在的报告目录（约1.4/1.9MiB），没有重新运行诊断或覆盖旧轨迹/CFD。40文件SHA验证一致，5份MD/HTML的43个相对链接全部有效；身份清单HISTORICAL_CAUSE_REPORT_COPY.json。P1说明亦已同步，脚本快照仍55份/490765078a61df9ad8d0a109c4afcabe5365f0224399c59a874e00e4d817e2a4。
- 细档当前最近保存step20，连续稳态0/5；07:18最新NS27-6，此后处于矩阵处理，健康检查无新增异常。继续保持当前会话监控。

2026-09-30 07:24 UTC：M1_fixedplc_fine保存step30/t=0.3s，mass=9.313227840223519e−10、no-slip PASS；20→30区间E_u=6.789556283243992e−6已低于1e−5，但E_Q=1.1761763468993745e−6仍略高于1e−6，连续稳态0/5，未放宽阈值。当前step31-3，无发散/非线性失败/OOM，内存约23.12GiB。上一07:22记录的NS27是较早历史，矩阵处理已结束并继续推进。当前没有活动本地render/collection任务；55脚本快照不变，固定边界细档CFD及后续46条重放仍须完成。

2026-09-30 07:45 UTC：M1_fixedplc_fine保存step40/t=0.4s，mass=9.091486967196157e−13、wall no-slip PASS；30→40区间E_u=1.3432829414410796e−7、E_Q=3.6928835672948594e−9，首次在线稳态PASS，连续1/5。当前step42-1，无发散/病态矩阵警告/非线性超限/OOM，内存约23.28GiB。仍需后续4个连续区间和独立postprocess/field audit，不能提前认定qualified；随后同一46条重放和完整固定边界中细档gate。继续保持活动会话监控，55脚本快照490765078a61df9ad8d0a109c4afcabe5365f0224399c59a874e00e4d817e2a4不变。

2026-09-30 08:02 UTC：M1_fixedplc_fine保存step50/t=0.5s，mass=1.2352563814125215e−15、wall no-slip PASS；40→50区间E_u=2.9470104967084617e−9、E_Q=1.468463962963338e−11，连续在线稳态2/5。最近NS50-2s，内存约23.37GiB，无发散/病态矩阵警告/非线性超限/OOM。仍需3个连续区间与独立后处理/全场审计，固定边界细档46条重放尚未开始，完整pair仍PENDING；保持当前会话监控，不提前结束。脚本和科学输入未改动。


## 固定原边界轨迹终态报告：2026-09-30 08:18 UTC

- 新增 scripts/report_fixedplc_trajectories.py、FIXED_PLC_TRAJECTORY_REVIEW_ZH.md、05_particle_replay/FIXED_PLC_TRAJECTORY_REVIEW.json；三组 figures/fixed_plc_censored_endpoints、fixed_plc_failure_replays、fixed_plc_slow_fast_replays 的 PNG/PDF 已目视核查，加入总报告、HTML、brain_report 索引。只读已关闭重放，细档关闭并完成比较后自动扩展，不提前填入细档结果。
- 中档全部 46 条末接受状态块 SHA 与 MANIFEST 一致，位置与分析末点一致；观察截止记录的原始 v/x/gap/time 逐项精确核对。已完成中档仍为40出流、6观察截止、0失败、0支持静止，18筛查事件。六条删失 ID9/393/1017/1329/1377/1463 末速分别0.177625/0.301134/259.897707/0.516910/0.894015/0.864757 mm/s；ID1017末gap347.898204µm，但全程最小gap约2nm，不能把删失都解释为停滞，也不能用远壁终态排除历史近壁段。ID301中档出流仅早于统一观察终点2.706580758ms，原窗口和门槛均保留。
- 当前56脚本 snapshot SHA2b6c7cdab61d5829454357e2c23393ae265e7fdc689ae78babfdc7818f5a15dd；本地语法与服务器当前/归档哈希逐项通过。仅重启report_watch，CFD不受影响。本地总报告08:15:07UTC，136链接/78图SHA核对PASS；服务器刷新核验待完成。
- 最新08:17监控：细档仍为NS57-2，在线稳态2/5，健康检查无异常；等待step60及后续连续稳态、独立场审计和46条重放。用户要求保持会话监控，不能因为后台任务尚在运行而结束本轮。

2026-09-30 08:20 UTC：M1_fixedplc_fine保存step60/t=0.6s，mass=1.764651973446459e−16、wall no-slip PASS；50→60区间E_u=6.207878511094995e−11、E_Q=1.216286372697972e−13，连续在线稳态3/5。当前NS61-2，无新增求解异常。服务器08:17:49UTC报告已完成136链接/78图SHA核验，新增fixedPLC轨迹helper及46末块/曲线身份也通过。仍需2个连续区间与独立场审计，然后固定46条重放；完整pair PENDING，当前56脚本snapshot 2b6c7cdab61d5829454357e2c23393ae265e7fdc689ae78babfdc7818f5a15dd。保持活动会话监控。

2026-09-30 08:38 UTC：M1_fixedplc_fine保存step70/t=0.7s，mass=3.529303946892918e−16、wall no-slip PASS；60→70区间E_u=1.2469250801689046e−12、E_Q=2.8675594568504962e−15，连续在线稳态4/5。当前NS72-2，无新增求解异常。仍需1个连续区间、独立后处理/场审计和同46条重放；完整pair PENDING。当前56脚本snapshot 2b6c7cdab61d5829454357e2c23393ae265e7fdc689ae78babfdc7818f5a15dd不变。继续保持活动会话监控，不提前结束。


## 固定原边界细档 CFD 完成：2026-09-30 09:02 UTC

- M1_fixedplc_fine 于08:54 UTC正常结束81步/0.81s，执行PASS、12419.276369572s。step80相邻检查E_u=2.747384083001279e−14、E_Q=1.764651973446459e−16，5个连续稳态通过。独立postprocess PASS，最终入口3.000000000000002e−7m³/s=18mL/min、mass打印0、wallu0；出口份额0.3333067539875032/0.3332000184611204/0.3334932275513764；p范围0.469298241689–2727.24882455Pa、umax0.555482321373m/s。最终u/p/WSS变化4.233316802e−15/1.146578582e−15/6.929122183e−15。
- flow SHA b507a8053fc475e8c7728d5bb57040a825eba82c68c6efa03df58b997b455927；mesh SHA d94d41efc45cec9db6f82af542fecdf50560837a59c0cb498aa68d8b7f3fadd0；flow_arrays SHA0b49d72da37e561823031eeeb6b7d5d6dc5c428b1e08abcdccde0d2172f95a8b；WSS SHA6ce78e0b3e154dd5af8c5c1c5a659ffc5f5a6d72adf680b52cc8d453f4c81c32。895647个单元的局部顶点顺序均恢复canonical后求梯度。WSS面积均值3.07660478826Pa、P955.72242383601、max21.1300748585、严格零壁面0。
- 全场审计119.275179s；zero_tetra=[]、1059固定点全部域内、CPU/GPU G差1.33582e−12（不代表物理准确）。原388769中心细档tet99026，u范数0.0858001085324m/s、div+105.313603521s⁻¹、G范数2329.91755075s⁻¹、eta0.045200570933、gap47.398452042µm。原138584中心细档tet881160，speed0.0217377561409、div−32.8983149133、G1332.14411060、eta0.0246957627569、gap26.2089336968µm；仍在域内。
- 高eta球区63668cells，eta≥0.5体积7.56859899273e−13m³（中档2.92248157399e−12），P99absdiv291.708110191/rQ0.029976432593/chi0.424392404；原zero球区23716cells，eta≥0.5体积1.71999785850e−12，P99absdiv213.662427744/rQ0.056900913965/chi0.874613304。固定中细1059共同点ALL u/G/WSS相对L2=2.17688614%/23.85911261%/17.53807549%；高eta4.22512504/44.32914540/31.08804492%；zero2.52025058/34.51213163/23.32111736%。局部场仍明显网格敏感，不能根据异常体积减少宣布收敛。完整pair等细档重放。
- 细档关闭PETSc profile SHA09fbb4d4c2adcc72270c39a5f1caf5d6e4580bbba22482b69fdf4dd241c3536d：19次LU最大累计7961.1s、rank max:min47.2；MatSolve3578.4s、MatMult1595s（嵌套事件不相加）。PETSc总12410s；子域分解负载严重不均衡解释长CPU阶段，不把占满配额当效率认证。原profile配置未改。
- 08:57 UTC后细档46条replay开始，09:01已关闭17/46，仍RUNNING。collection已核验8个新CFD+8个已完成replay；活动第9批不宣称完成。本地compare_fields/update_comparisons/assess已执行，render_report仍进行（PTY66270）；fixedpair PENDING，未进入Phase3/4/5。当前56脚本snapshot仍2b6c7cdab61d5829454357e2c23393ae265e7fdc689ae78babfdc7818f5a15dd。保持活动会话监控。


## 用户再次请求暂停：2026-09-30 09:06 UTC

用户要求“在下一个可停止点暂时停止”，此指令覆盖此前持续开发/监控要求。立即停止新增科学分析与开发。细档CFD已经完成并收回（8套新CFD、8套已分析重放本地核验PASS），细档第9套46条重放仍在运行。实际replay_cohort.py没有中途续算支持，采用自然完成当前批次并封存的可持久化停止边界，避免丢失缓冲状态。

已在服务器09:05:40 UTC向fixedplc_validation_followup单独的父调度PID695948发送SIGSTOP，子replay PID734227继续当前批次；随后report_watch已STOPPED。不能对整个进程组发送停止/杀死信号，直到当前replay已封存。PAUSE_REQUEST_20260930.json记录请求、进程身份及禁止自动恢复。09:05:55监控当前27/46关闭，余19条，绝对活跃时间4.5835–5.796s。此阶段仅安全封存与暂停核验，不启动analyze_replay、完整pair评估、新WSS诊断或新计算。此前准备进一步定位WSS峰值的工作尚未执行，应等用户恢复。

本地09:01:11、服务器08:57:55总报告均通过136链接/78图SHA核验，含8套新CFD、8套已分析replay。M1固定边界细档原中心图、剖面图和10个关闭PETSc执行段成本图已目视核验。本地无活动render/collection作业；PAUSE尚处于等待当前重放封存，不可提前称全部停止。


## 暂停完成：2026-09-30 09:14 UTC

- 服务器09:12:55UTC已全部安全停止：本研究supervisor全部STOPPED/EXITED、autostart/autorestart=false、无遗留研究计算进程、GPU compute列表为空，实例保留。父队列先SIGSTOP，细档原始重放自然完成后再停止父进程，未强杀活动积分。
- 细档第9批原始重放完整关闭：40EXITED、6RIGHT_CENSORED、0SOLVER_FAILURE、0SUPPORTED，875.299645901s。46个MANIFEST/46状态块已复制到本地且与服务器SHA一致。未运行该批analyze_replay，未生成完整fixedPLC中细档gate；不能根据两档40/6相同汇总称轨迹收敛。
- PAUSE_RESUME_ZH.md、PAUSED_STATE.json/PAUSED_STATE_20260930.json、PAUSE_LOCAL_COPY_VERIFIED_20260930.json为当前恢复依据；上次9月29日记录已另存历史。HTML和状态/科学报告添加行政暂停标识，原数据/门槛不变。
- 没有需续算的CFD或半截轨迹；恢复后先分析已关闭细档run，再update_comparisons/assess/render。不得直接重启已完成的fixedplc_validation_followup或覆盖已有replay。用户未明确恢复前，不开展新分析、开发或计算。WSS峰值位置进一步诊断尚未执行，Phase3/4/5未启动。

### 2026-09-30 10:09 UTC：用户明确恢复（覆盖前述暂停状态）

98 个暂停封存文件和 56 个当前脚本/归档本地、服务器 SHA 均通过。已启动独立 supervisor 任务 analyze_fixedplc_fine_resume_20260930，分析已自然关闭的 M1_fixedplc_fine 46 条轨迹；旧 producer 未重启。PAUSED_STATE*.json 保留为历史快照；当前状态看 mesh_particle_validation_v2/CURRENT_EXECUTION_STATE.json 与 PAUSE_RESUME_ZH.md。下一步完整中细档比较并审查 Phase 2G 后续细化。

### 2026-09-30 10:39 UTC：同几何面内/体网格补充验证进行中

本次用户明确恢复后，fixedPLC 细档已完成 GPU 分析，9 批重放、8 套 CFD 收回核验；原固定边界中细档完整 FAIL：0/46 通过、13 条均出流者改出口、4 条终止变化。新增 WSS 峰值同壁面交叉核对，峰值在原高 η 区域，非入口。

继续 Phase 2G：同原始 PLC 分片平面的 50/25 µm 热点面内细分，四端口三角形保持。原生 50 µm 输出含 85 个四壁点及 float32 约 1.85 nm 偏差，FAIL 保留。证明输出边界与声明输入的 float32 舍入像逐点一一对应后，派生网格恢复输入精度，并对四壁点单元涉及的全部内部边作完整邻接扇区一致二分，增加内部自由速度节点；是 CFD 前网格操作，无已有 field 修改。独立完整 QA 保留严格 0.1 nm 距离、原物理点与 all-WALL=0 门槛。

M1_coplanar_medium 已完整 PASS：804079 tetra /158908 nodes，77 条内部边二分消除85个四壁点，父子体积相对差最大5.55e-15，原体域体积差0，全部固定点域内。实际 high-η 区域壁面边长 median42.95µm/壁邻接高度48.04µm。细档 native 原生1063718tet/213845nodes、353个四壁点，当前 mesh_coplanar_refinement_queue 仍在完整原生审计；额外恢复精度试算无负/零体积。活动作业属于本研究独立 supervisor，不重启旧停止任务。

64 个当前脚本快照 e639dec6f8a0231802c90f55108d39ba1097214548cc7a5a970e02af922c9dd6；coplanar_validation_followup.py 已实现但尚未运行，须待两档网格全部完成并生成明确资源/网格审查文件 COPLANAR_CFD_REVIEW.json。新资源/构造计划、尺度报告以及未来所有源版本均保留。276 个受保护原件及冻结副本重新校验 PASS。Phase3/4/5、生产发布均未开始。


### 2026-09-30 11:03 UTC：coplanar 两档 QA 通过，已启动串行 CFD 验证

两档最终完整几何/拓扑 QA PASS：804079 / 1064977 tetra，0 四壁点，父子体积守恒、四端口不变、原分片平面域保留、固定 1059+4 点全部域内。细档通过 288 条真实内部边完整扇区二分消除原生 353 个四壁点；最小内部节点壁距0.533634µm。实际热点壁面边中位数约43→21µm，壁邻接高度约48→24µm。COPLANAR_CFD_REVIEW.json 已记录门槛/mesh SHA、单 CFD/原6h上限、RAM与原物理条件。独立 supervisor 的 coplanar_validation_followup 已启动，中档 CFD RUNNING；report_watch 同时持续更新，只在各输出关闭后做比较。生产场未改。

新增等价壁面诊断：对 M0 原 46 条全部 417292 个接受状态查询，coplanar 中细档无显著法向差；M0 对新表示仅一个显著差异，为 ID393 t1.632s 的双壁近等距位置。固定原 M0 u/G、位置/半径/黏度和原模型，五种壁面表示 GPU 瞬时求解，M0 下一接受区间平移速度重现误差3.03e-20m/s、无fallback。固定原面和coplanar都选择另一法向，接触法向集合仍完全相同；平移向量差2.83152048e-6m/s、角速度向量差约180.975rad/s。进一步分别交换gap/normal以区分舍入影响；不能从一个瞬时状态外推整个轨迹因果或作视觉平滑修正。

68脚本快照754474bb0c0076fcdfe2ff0b64c233993c75629946886fbf04ccefe03818276c已同步，所有已执行版本保留。当前执行状态见 CURRENT_EXECUTION_STATE.json，不用历史 PAUSED_STATE 推断仍暂停。Phase3/4/5尚未开始，旧停止任务不恢复。


### 2026-09-30 11:12 UTC：壁面表示诊断和报告整合已关闭

NORMAL_CONTROLS.json 使用相同输入分开交换 gap/normal，平移向量差分别为1.05340606e-16m/s与2.83152048e-6m/s；完整表示替换2.83152048e-6m/s。单点差异由最近法向选择主导，相对M0该点速率的向量差12.3089%，角速度向量差180.975rad/s。两法向夹角35.752°，完整接触法向集合相同，无fallback，不重新积分，不能外推整个轨迹的因果比例。新WALL_REPRESENTATION_REVIEW_ZH.md与300dpi PNG/PDF已生成。

MESH_CORRESPONDENCE 更新覆盖coplanar两档：M0球区外683415个参考单元中仅118/128个四顶点保持；中档→细档球区外721858个单元仅201个保持。球区外原WALL节点均保持，但体网格非严格嵌套，不能把后续差异全部归因于局部热点。

总报告整合的page变量顺序错误已修复，局部/服务器完整render均成功；本地11:10:25UTC、服务器11:11:11UTC报告均通过151链接/86图件SHA核验，所有新子报告输入/源哈希和10个相对链接也PASS。错误仅在派生报告，无原科学数据变动。68脚本快照6bbbc37f3387fda11102390eaf6207f353a95fb77bd2b1bca36173bf8f8bfa2a，本地/服务器当前源与归档均PASS。当前coplanar中档CFD第5步，未达首个step10保存点，内存约20.8GiB，无发散/病态/非线性失败标记或OOM；8套CFD+9批已分析重放的完成数不变。继续活动监控，不提前结束。


### 2026-09-30 11:34 UTC：coplanar 输出读取技术问题与恢复验证

2026-09-30 11:34 UTC：开发仍活动，未收到用户暂停。原中档运行因输出坐标身份读取问题在原生 step13 安全关闭，保存 result/stFile 010 和013；旧 coplanar_validation_followup EXITED，细档 CFD 与新微泡批次尚未启动。原生输出/检查点不得覆盖。reports/execution.json 不存在，因为原 runner 最后 read 也抛异常；不能编造原 native exit_code 或将原运行标为 PASS。

已核实：输入 Float64，native VTK Points Float32、Velocity/Pressure Float64，max坐标舍入1.878042915e-9m。原读取器阈值约1.07846e-11m，拒绝匹配。158908节点恰为输入float32编码的单射/双射，全部804079个tetra行物理节点一致。原源码VtkData.get_points读double，set_points新建默认vtkPoints未设置double；vtk_xml写原坐标，源码和程序hash保存在01_mesh_capability/NATIVE_OUTPUT_IDENTITY/。

native_output_identity.py 为隔离class，继承原measure积分，仅替换read；精确编码、无碰撞、双射、tetra行身份检查，无距离容差放宽、无插值。canonical_flow_grid 只建立不同身份的派生场，使用实际求解器输入坐标和原节点u/p，保留原native VTK。audit_native_output_identity已PASS：原10/13步可读，mass分别0.0009765625/0.0001220703125，启动未稳态；已完成fixedPLC中档与旧reader的u/p和积分字典逐项相等；ULP坐标破坏、tetra节点变化、舍入碰撞三个反例均拒绝。

当前独立job check_coplanar_output_restart 启动；recover_output_identity.py --stage check 创建segments/output_identity_check_10_to_13，通过coplanar_output_runner.py从原10检查点重复11/12/13；核对u/p及8MPI的Y/A历史，原门槛1e-8。原base runner SHA f27feeaa5fdbd39581117b6bf93ed18b9ccf31412fa65a1af5feebdeaebacac5、物理求解程序43cd7...不改，PETSc/BC/网格不改。验证段的raw CFD状态因非稳态预期FAIL，连续性判据独立，要求exit0且health无fail。

下一步（待完成）：
1. 监测并核验10→13连续性；若失败保留证据并排查，不放宽门槛。
2. 修订研究postprocess_candidate.py的coplanar分支：绑定reader审查和源SHA，使用身份适配器read；原u/p/WSS/质量/稳态阈值保持；canonical派生steady_flow.vtu与原native VTK分开sha，不对原文件改坐标；保存原文件/节点对应/转换身份。不要改B/vendor或生产postprocess_flow.py。
3. 连续性PASS后recover_output_identity.py --stage resume从原13检查点继续，segments/output_identity_resume_0001；已实现，只会求解至稳态并写EFFECTIVE_EXECUTION/RESUME_INDEX，尚未启动。原6h预算扣除原运行保守上界与验证成本，预算不能自动重置。
4. 新串行恢复队列（尚未实现）接上中档续算→独立postprocess→audit_flow→原46重放/分析；然后准备尚不存在的M1_coplanar_fine，通过coplanar_output_runner.py从零执行并相同后处理/重放。旧producer不可重启覆盖原case。使用HEAVY_COMPUTE.lock，添加新producer给report_watch。
5. 更新报告监控：原monitor未显示read_retries，已在本地添加输出读取失败数/最新错误/technical hold/validation_only；这项新修改尚未部署快照。RESTART_CONTINUITY / OUTPUT_IDENTITY_STOPPED / RAW日志均保留。

当前73脚本快照e07f852b08ad00f13d3c97b1b602409ab6928da1666d3f55fe7ad9956f285975服务器已同步；monitor_progress.py后续本地修改未包含于该快照。原8套CFD+9批重放仍为已关闭数量，production、Phase3/4/5未改变。用户要求持续监控，不在后台运行时提前结束。


### 2026-09-30 11:43 UTC：输出身份与连续性通过，中档正式续算已启动

RESTART_CONTINUITY PASS，10→13重放u/p相对L2=3.7339515694e-15 /2.10441913865e-15，最大u/p绝对差7.68829444553e-15m/s /5.00222085975e-12Pa；8rank Y/A最大相对差2.08286188685e-14，原门槛1e-8未改。短段原execution FAIL仅因没有稳态停止条件，exit0/health无fail；其科学身份是重启连续性验证。

canonical VTK写回再读坐标/单元/u/p逐项完全相等，原native不变。研究postprocess_candidate.py仅对两档coplanar走已审查身份适配器，原物理积分/WSS/稳态门槛不变。派生steady_flow使用实际solver输入Float64坐标和原始节点场，记录native与派生场各自SHA，原native保留；没有投影、插值、平滑或求解器变更。

新producer coplanar_identity_followup已启动；中档segments/output_identity_resume_0001从原13继续，剩余CFD墙钟上限19079.52710247s（原6h扣除原运行保守上界与验证成本）。新队列将中档独立后处理/场审计/固定46条重放分析后，再以相同身份读取方式启动尚未开始的fine。旧coplanar_validation_followup不可重启。report_watch已启动，新producer已登记。监控已显式输出read_retry_count/error，原读取问题不可再被健康标记掩盖。

新INTERNAL_P1_ZERO_FACES复核：M0两张内部三WALL点零速面分别邻接138584与578387/145555，10个通过QA候选该计数均0。新fixedPLC全417292个M0记录位置查询：medium/fine法向完全相等、接触数量与开口标记相同，gap最大差7.3319e-18m，不能用M0对新表示的单点法向切换解释固定边界两档大量轨迹差异。详见WALL_REPRESENTATION_REVIEW_ZH.md与INTERNAL_P1_ZERO_FACES_ZH.md。

当前74个脚本，最新快照以CURRENT_STUDY_SCRIPT_SNAPSHOT.json为准。当前后处理及监控/报告修改已本地归档，下一次部署与完整render校验进行中。原8套CFD/9批已分析重放的完成数不变，生产场和Phase3/4/5均不变；持续活动监控。


## 2026-09-30 12:29 UTC：coplanar 中档稳态完成，微泡回放运行中

用户“恢复开发”授权仍生效，持续监控指令仍生效；当前不是暂停点。新队列 coplanar_identity_followup 运行中，report_watch 自动更新报告。M1_coplanar_medium 在 segments/output_identity_resume_0001 的 step81 正常完成，5/5 稳态区间，续算壁钟 2482.025157 s；独立后处理与完整 field audit 已完成。原 reader 中断段与 checkpoint 10→13 连续性验证段保留，不混入合格 CFD 数目。

规范 flow SHA256 `84d99ab2e2eb2084531f4e6a4dd18f0249f60eb7596602a17ff571fa602e4cda`，原生 flow SHA256 `d19431c60297caf914a538086ad804de7ba8eabc4c07a1f95f5d308b084083b8`，mesh SHA256 `cfd5baafbbe8074c655250b1399da2275d172f2804420357c2d9e43781b28c4a`。新增 verify_canonical_final.py 用独立字典序精确编码匹配、编码无碰撞、逐行 tetra 成员和节点 u/p 字节相等检查，服务器 PASS；记录在 reports/CANONICAL_FINAL_IDENTITY_VERIFIED.json，首次临时检查另行保留。规范文件使用真正 solver 输入 Float64 坐标，原生 Float32 输出仍保留；无物理场插值。

全域严格零速 tetra 为0，1059 条固定记录全在域内。原 388769 中心对应新 tetra34821，速度0.06727204861 m/s，div +516.4051726 s⁻¹，G范数1704.419239 s⁻¹，η0.3029801358；原138584中心对应393782，速度0.02092170747 m/s，div −67.86938086，G范数798.2662597，η0.08502098145。原高η球区 η≥0.5 体积2.062699337e−12 m³、P99 abs(div)610.8784339 s⁻¹；原零速球区相应为1.738086977e−12和375.4568453。不能把η降低或零速结构消失当作局部精度认证。入口18mL/min，出口份额0.3332659859/0.3332513332/0.3334826809；WSS均值3.072695514 Pa，最大13.24224440 Pa。

46个冻结初态的中档 replay 已启动，12:29UTC 为10/46终态，尚未完成，不能提前宣称无失败。之后队列才启动 coplanar fine CFD 与相同46条回放。已完成正式验证结果数量为9个新CFD、9个已分析replay；新的中档replay尚不计入完成。ACTIVE_FLOW未改，production mesh未选，Phase3/4/5未启动。

新增 report_coplanar_validation.py 已用闭合中档场生成目标位置图与 PENDING 报告；其完整 replay 分支须待真实数据关闭后核对。报告入口将加入 COPLANAR_VALIDATION_REVIEW_ZH.md 与对应 PNG/PDF，原始科学结果不改。当前76个study脚本和精确归档在服务器全部验证，snapshot SHA256 `8110cdd5fbc2ea70ada848ffae344c82fc0a468529dd9e15ed684a89599318e1`。


## 2026-09-30 12:43 UTC：coplanar 中档重放完成，细档启动

M1_coplanar_medium_dt0.00050000 已关闭并完成原始接受状态SHA审计，耗时856.105162 s：41 EXITED、5 RIGHT_CENSORED、0 SOLVER_FAILURE、0 SUPPORTED_STATIONARY。O1/O2/O3=12/16/13，slow-fast筛查25次；删失ID为4、393、1377、1449、1522。ID440→O2，停留4.640216819s；1505→O2，3.848075123s；636→O3，3.913735675s、最小间隙13.23816035µm、无筛查事件；1189→O3，6.051869297s、1事件（t=5.6225s）；1329→O2，9.673118431s、无事件。ID393观察截止时仍运动，不能称静止。与原M0或其他网格差异不构成完整中细档收敛。

细档 M1_coplanar_fine 已由同一审核过的队列从零启动，12:43仍处于首步初始化，未出现求解告警。当前仅1个昂贵CFD。持续监控指令仍生效，不要结束开发后让服务器任务独自运行。9个新CFD/10个已分析回放均已收集到本地并核验，原正式ACTIVE_FLOW未改。

本次收集发现续算段的SV_MESH绝对服务器链接会被通用元数据同步带入本地，已只将本研究副本中的已知链接转为../../SV_MESH；原路径兼容链接及服务器链接未改。collect_completed.py 已改为完整排除SV_MESH并逐段验证本地目标，9/10结果再次收集验证PASS。

报告已处理原中断段无完整execution.json的真实情况：不填补未知总耗时，中档续算2482.025157s、重启验证273.273859s分列，原段PETSc内部计时另列。已完成的续算不再因旧progress.json误列为活动CFD。新增同几何中档场/轨迹报告已在服务器用真实闭合数据执行，完整中细档门槛仍PENDING。当前76脚本快照及服务器归档校验PASS，SHA256 `b758612fd00f4998498f63269e6375d2013138326fc272761535ac33496eb85d`。


## 2026-09-30 13:09 UTC：细档第4步；中档路径壁面查询已核验

M1_coplanar_fine 仍由 coplanar_identity_followup 单独运行，刚完成第3步进入第4步；尚未到第10步首次输出，没有合格细档场或细档回放。CPU端线性系统准备有数分钟的停顿，实际残差随后继续下降；首步第二次线性求解585迭代达到9.46e-11，未触发失败。内存约30.7GiB，所有求解告警计数为0。用户持续监控和恢复开发指令仍生效。

新增只读 GPU 诊断 audit_coplanar_medium_path_wall_queries 已关闭，11.760177s、峰值CUDA分配48081408字节。对中档闭合轨迹全部376929个接受状态（221465个近壁查询）交叉查询M0/两档coplanar壁面：中细档法向向量最大差0、接触数/最近特征/rim分类差0，最大间隙差8.4568e-18m、最近点差1.0408e-17m。相对M0仍有最大法向向量差0.03949236、最近点差6.6603e-8m，未隐藏。证据在01_mesh_capability/COPLANAR_MEDIUM_PATH_WALL_QUERIES/，全体状态块、查询数组、源码/边界SHA已在本地报告中校验；这是固定已记录位置的查询对照，未重新积分，不外推细档未来路径或完整接触求解。

报告补充明确slow-fast筛查定义：50ms回看，之前最低速率<2mm/s、当前>=10mm/s，50ms事件间隔；0次事件不表示没有低速后加速，ID440/1505曲线上仍有较慢的加速过程。阈值与原始轨迹未改。当前76脚本快照SHA256为d656c1e0b89f2e60db5ae66bc797747b21bf19df30fe24f37c99c511bca86954，本地与服务器当前/归档已验证。最新本地总报告12:55:21UTC，166链接/94图文件哈希PASS；服务器总报告最近12:41:15UTC，后续新helper将在下一次报告生成时集成。


## 2026-09-30 13:54 UTC：细档输出身份验证通过；仍在启动求稳态

细档 M1_coplanar_fine 最新完整原生日志为 step16 / nonlinear10，尚未达到稳态；coplanar_identity_followup 和 report_watch 活动，用户持续监控与恢复开发指令仍生效。首份 result_010 和 stFile_010 已保存并收集本地，FIRST_OUTPUT_CHECK.json 与 FIRST_OUTPUT_COLLECTION_VERIFIED.json 均通过。214133 个节点满足输入 Float64→原生 Float32 的精确编码双射，1064977 个 tetra 行的物理节点集合均相同；原生局部顶点顺序变化不影响身份。最大坐标序列化差 1.8780429154752845e-9 m，未插值或改写原始节点 u/p。step10 质量误差0.0009765625，属于尚未稳态的启动输出，不作为合格细档场。

周期性的长等待已有对应日志证据：实际 SV13Q_BEGIN 的第8、13步均记录 MAX_AGE 重建并正常结束；第1、3步另有 ITERATION_GROWTH 重建。服务器构建源码目录 sv13q_reuse.h 当前 SHA256 db635f95a7e0c187fa9768f1a836b57034bef4aca4319c2c24b0cd8fbddae220，其 adaptive 策略为最多复用5步、健康求解迭代数超过参考1.5倍后请求下次重建。此为当前源码与运行日志的对应核查，不补造编译时未单列的 header 哈希，也不将等待时长直接全部归因于 LU；待本次关闭后的 PETSc profile 才能核对事件耗时。当前冻结 PETSc 选项和可执行程序未改变。13:53:49→13:54:29 的 cgroup CPU 消耗平均约7.683核，配额7.68核；内存约31.10GiB。CPU消耗包含MPI等待，不能等同有效浮点工作或物理精度。

本地 MESH_CONVERGENCE_ASSESSMENT 的待完成清单此前停留在中档回放完成前；已用现有脚本和闭合数据重新生成，去除错误的中档缺失项，原先四组中细档 FAIL 和 coplanar PENDING 不变。正在据此更新本地总报告。仍为9个新CFD、10批已分析回放，细档轨迹尚未启动，ACTIVE_FLOW及Phase3/4/5不变。


## 2026-09-30 14:17 UTC：细档step20已保存；当前step24

用户恢复开发与持续监控仍生效，coplanar_identity_followup / report_watch 活动。细档第20步输出已由既定身份适配器成功读取，读取重试0；epsilon_mass=9.536743166101389e-7，wall no-slip PASS。10→20区间 E_u=0.0006055167280575586、E_Q=0.0005345753589403937，尚未满足既定稳态门槛，连续稳态区间0/5。最新已进入第24步，内存约31.4GiB，求解告警0；第23步MAX_AGE重建后已正常收敛，没有失败重试。细档合格场、细档46条回放仍未完成，9个新CFD/10批已分析回放数量不变。

监控程序新增 latest_linear_solve：只解析完整日志行，显示最新SV13Q_BEGIN的step/logical/attempt/rebuild_reason以及准备中、迭代中、完成状态；当前残差只取该轮BEGIN之后的日志，避免CPU准备期间误显示上一轮残差。已用真实第18/23步的准备和后续迭代记录核验。未改变求解器或科学模型。最新76脚本快照 SHA256 5ea3e20f7057171e7bea1605c0f1a69b354f00a3a473e855bb6517bcbb47a9b1，服务器当前与精确归档均PASS；CURRENT_STUDY_SCRIPT_DEPLOYMENT和CURRENT_EXECUTION_STATE已收回本地。

本地报告已刷新到13:54:38.519977UTC，166链接/94图文件哈希PASS；coplanar待完成清单仅列细档场/审计/回放/配对比较，不再错误列中档回放。最近监控快照在monitoring/LATEST.json，实际时间以文件updated_UTC为准；该本地文件不是实时推送。


## 2026-09-30 14:45 UTC：用户要求安全暂停，细档第37步检查点已保存

最新用户指令是暂停，并汇总已完成任务、结果和后续工作。当前已PAUSED_BY_USER，恢复必须等新的明确指令。coplanar_identity_followup和report_watch均STOPPED，全研究无RUNNING/STARTING/BACKOFF/STOPPING任务，无研究进程、无GPU计算进程。原生CFD通过STOP_SIM在M1_coplanar_fine第37步（0.37s）保存并exit0，health_failures为空。原execution.json保留FAIL/stop=null，因为外部用户停止前未完成连续5区间稳态，不是新的微泡数值失败。

8rank完整Y/A检查点stFile_037.bin SHA256=243b56cd70cb32f40201d18d12b8d71c4990a958fff0c6d75a3032e3b170592d，14531136字节。第10/20/30/37结果和检查点、原生日志、输入及执行身份已收集本地，30文件哈希PASS。服务器证明reports/USER_PAUSE_CHECKPOINT_VERIFIED.json，SHA256=62fcd71b90009ba03dd3b948ccf655a3b92fb0bf624093612bbcc0283aa23b23；本地USER_PAUSE_LOCAL_COLLECTION_VERIFIED.json。原生37步质量误差7.276013e-12、无滑移PASS，但未资格化。最后完整20→30区间Eu6.9514836e-6/EQ1.1882023e-6，EQ仍超门槛。耗时7422.8621943s，原6h预算余14177.1378057s，后续重启验证须从该剩余量扣除，不能自动重置预算。

机器暂停文件PAUSED_STATE_20260930T144554Z.json；CURRENT_EXECUTION_STATE已指向该暂停，旧PAUSED_STATE_20260930.json的09:12暂停历史不覆盖。PAUSE_RESUME_ZH.md已改为完整阶段总结和待办；暂停前报告/HTML/状态保存在pause_snapshots/20260930T144554Z/。仍为10个QA通过候选、9套合格新CFD、10批已分析46条回放；四组已完成中细档比较全部FAIL，同几何细档PENDING。ACTIVE_FLOW/旧1565/微泡物理/Phase3/4/5均未改。

恢复时不能直接重启旧coplanar_identity_followup（会重复中档续算和撞已存在目录）。需要隔离细档原30→37连续性验证，比较u/p和8rankY/A，再从37独立段续算；当前recover_output_identity.py主流程硬编码中档10→13，须适配新入口。之后完成细档稳态/后处理/原生规范场身份/场审计/同一46回放/中细档门槛；空间未稳定时继续有界加密可行性审查。不得越级启动dt/P2/严格无散度阶段或新正式大批量。

当前76脚本快照6500430ebf6a370ae6456ad648e2969a3f7edc96b6fcb593cb1930fbe15870f1，服务器当前与归档核验PASS。最后只做暂停报告渲染与链接/图哈希收尾，不启动新科学计算。


## 2026-09-30 15:40 UTC：用户明确恢复，细档30→37连续性验证已启动

最新指令为“恢复开发，在开发停止之前不要停止 thinking 的信息显示”；14:45暂停要求已被新指令撤销。当前正在隔离supervisor任务check_coplanar_fine_restart_37运行30→37原生检查点重放，report_watch也已恢复。旧coplanar_identity_followup保持STOPPED，不能重新启动。原始37步检查点和30个暂停文件全部哈希一致；本地276项受保护输入再核验PASS。

新增scripts/resume_coplanar_fine.py；recover_output_identity.py的prepare/validate已参数化，主入口仍保留历史中档流程。新细档检查保持u/p和8rank Y/A相对L2<=1e-8，只有PASS才能从37继续。检查与续算合计不能超过原剩余14177.137805700302秒。RESUME_AUTHORIZATION_FINE_37.json记录本轮授权和源身份；CURRENT_EXECUTION_STATE状态已恢复活动，历史PAUSED_STATE_20260930T144554Z.json不改写。

新scripts/fine_user_resume_followup.py已完成并部署，但尚未启动；待检查PASS后按续算、后处理、独立原生/规范场身份核验、全场审计、同46单泡回放、分析顺序推进。保持原生暂停execution.json的FAIL；后续有效执行另存EFFECTIVE_EXECUTION。发布派生流场身份前备份原STUDY_IDENTITY，原始求解输入/结果不修改。不得越级做dt/P2或发布生产网格。

本轮78脚本快照fe0eb1719f8b37aa8c2646e190613efb7e8b7e356be4a3f6701061ab9a10bc54，当前和精确归档服务器核验PASS。15:43监控：正在31步首次预条件器准备，CPU约7.68核，RAM约27.8GiB，GPU约3.2GiB，求解告警0；这些是运行状态，不是连续性通过证据。当前科学完成量仍为9新CFD、10批回放。


## 2026-09-30 16:05 UTC：30→37连续性未过导数相对门槛；已转入20→30有界验证

本轮用户仍明确要求恢复开发，未暂停。check_coplanar_fine_restart_37已于15:57关闭：native exit0、health_failures=[]、实际31至37步完整执行，耗时1017.5058188438416秒。u/p相对L2=3.69920752557423e-15/3.867596768777731e-15；Y各rank也约1e-15，但A相对L2=9.81848e-8至2.06983e-7，超过原1e-8门槛。reports/RESTART_CONTINUITY.json保留FAIL，绝不能说37续算验证通过。

已完成audit_restart_conditioning.py及FINE_RESTART_CONDITIONING_ZH.md/figures/restart_A_conditioning.{png,pdf}。原30步A压力分量范数1582.518降至37步15.7178；重放与原37步速度导数最大差2.488e-12m/s²，压力导数8.871e-9Pa/s；gamma·dt·deltaA/Y范数约1e-14。只支持相对范数条件性解释，不是未来误差上界，也不替代FAIL门槛。当前构建树initialize/Integrator/Array/output/main已核对直接Y/A读取、数组布局与广义alpha校正；不虚构编译时源码映射。

新技术恢复review=RECOVERY_REVIEW_FINE30.json，SHA3775cdde49a134130e436a5e6fc37dddbef42e1566ba183d4899303a17beb7f5。已启动新supervisor任务check_coplanar_fine_restart_30，通过recover_fine_checkpoint30.py从原20检查点重放到原30。保持u/p和各rank Y/A全部<=1e-8；新结果只写reports/RESTART_CONTINUITY_20_TO_30.json，不覆盖原37 FAIL。只在新验证PASS后，从原30（不是37）续算，段名user_pause_resume_from30_0001。原37及失败检查全部保留；原6h扣除原段7422.8621943和第一次验证1017.50581884后剩13159.63198685646秒，第二验证还须扣除。若新验证仍FAIL，不自动重试、改门槛或强行放行。

fine_user_resume_followup脚本已改为等待20→30 PASS，尚未启动；不要用旧fine37入口，其授权钉住旧脚本SHA。postprocess_candidate已按RESUME_INDEX.continuity_validation读取有效独立证明并核对原检查点SHA。recover_output_identity.validate新增输出路径参数和禁止覆盖既有证明。80脚本快照527be58ca626d9715dbc18b1d55fd794ebd95a2638df8763837e47358b879ebd已部署核验PASS。render_report新增诊断报告链接/正文的本地编辑尚未部署（等服务器当前render关闭后再部署，防运行中源身份混杂）。旧M3“等待回放”文字已改为实际完成FAIL，等待本轮render核验。

科学完成量仍9个新CFD/10批回放；没有新生产网格、dt/P2/无散度试验。继续监控并每约一分钟简短中文更新，不要在活动开发结束前final。

16:13补充：render_report的新诊断链接/正文已部署，80脚本当前快照b761d50efd3a385e7540dd57f3aebeb6b621478f49fe01af6f9d343e9ac2aece，服务器当前/归档PASS。报告16:09:36UTC已收回本地，167链接/96图SHA PASS。首次失败验证24文件本地收集PASS。诊断引用的5个当前构建树源文件已复制到01_mesh_capability/NATIVE_RESTART_CONDITIONING_SOURCES，附SHA和不声称编译源码身份的说明。20→30仍在运行，16:21处于第26步MAX_AGE预条件器重建，无求解告警。


## 2026-09-30 16:32 UTC：20→30原门槛通过，细档已从原30步开始续算

check_coplanar_fine_restart_30已PASS：u相对L2=3.7627402492217675e-15，p=5.7214329759145e-15，8rank最大A相对L2=1.0078585437251169e-9，小于未修改1e-8门槛；21至30步原生完整执行、exit0、health_failures=[]。检查段execution.json仍FAIL（不是稳态资格），连续性证明单独PASS。验证耗时1549.606749534607秒。第30→37旧FAIL原样保留，绝不能覆盖其RESTART_CONTINUITY.json。

fine_user_resume_followup已启动并运行recover_fine_checkpoint30.py --stage resume，从原stFile_030.bin恢复到segments/user_pause_resume_from30_0001。原检查点SHA5feed1c46deab62784326f8b543f308867bc3d9428f08500c384d6ab067f0f78，active proof=reports/RESTART_CONTINUITY_20_TO_30.json；不要误用原37失败proof。剩余CFD预算严格为21600−7422.862194299698−1017.5058188438416−1549.606749534607=11610.025237321854秒，从续算启动开始计时。实际XML只改Continue_previous_simulation=true、仍500步，PETSc/BC不变，8MPI+GPU。原37数据保留；30→37重算算入预算。

队列自动执行原5连续区间稳态CFD后，后处理、规范流场身份核验、全场审计、冻结46初态回放与分析；只有每项PASS才继续，不自动生产发布或进入dt/P2。report_watch运行。当前没有合格细档流场或细档微泡结果，仍9新CFD/10批回放。16:36监控处于31步首次预条件器准备，无告警。

新增report_fine_restart_recovery.py将两次检查、原检查点身份与预算整理为FINE_RESTART_RECOVERY_ZH.md/REVIEW.json，已服务器执行PASS，并由render_report自动引用。当前81脚本快照556bb73f90c87ff961c1672ae7355ce4fd78cb921f84669b0dd10d0f244cb490已部署且当前/归档SHA PASS。正在运行的fine_user_resume_followup保留启动时80脚本快照与实际源身份，不追溯改写；新增报告辅助不修改队列引用的科学代码。第二次验证全部文件也已完整收回本地并核对SHA。最新本地总报告仍16:09版，等本轮render关闭后再同步新版并核验链接/图哈希。


## 2026-09-30 16:55 UTC：最新用户指令为暂停，现已安全停止并回收

本条取代前述活动恢复状态。用户16:48:55 UTC要求在下一个可停止点暂停；原生续算段user_pause_resume_from30_0001已完成31至37步并退出0，health_failures=[]。fine_user_resume_followup和report_watch均STOPPED，全部研究队列STOPPED/EXITED，GPU计算进程为空；实例保留。第37步Y/A完整检查点SHA 276f72b5a23ec2203aa6829026c0322017cb929b3959f38e7251d2c1b1c79fe7，26个封闭文件本地回收哈希和独立checkpoint_audit均PASS。原37步文件、首次30→37 FAIL、20→30 PASS全部保留。

本续算耗时1042.891990184784秒，原21600秒总预算现剩10567.13324713707秒，休眠等待不计入。当前仍未稳态合格；原execution.json保留FAIL，不能改称成功或微泡失败。完成量仍9套新CFD、10批46单泡回放，无细档正式场或细档微泡结果，无新生产网格，dt/P2/无散度阶段未开始。

最新暂停记录：mesh_particle_validation_v2/PAUSED_STATE_20260930T165531Z.json；可读总结mesh_particle_validation_v2/PAUSE_RESUME_ZH.md；图形入口OPEN_RESULTS.html。CURRENT_EXECUTION_STATE为PAUSED_BY_USER，等待新的明确恢复指令，禁止自动重启任何队列。恢复前需审查本次37检查点的连续性路径和剩余预算，不直接复用现有段或旧启动器。


## 2026-09-30T17:56:09.883833+00:00：用户再次恢复，复用已验证第30步的新续算已运行

最新指令“恢复开发，在开发停止之前不要停止 thinking 的信息显示”解除16:55暂停。服务器fine_resume_followup_0002已于17:54:47 UTC启动，report_watch运行；旧fine_user_resume_followup及其他暂停队列保持STOPPED。新段为M1_coplanar_fine/segments/user_pause_resume_from30_0002，从字节未改动且20→30原门槛已PASS的原第30步检查点恢复。没有重跑验证，没有放行新旧第37步检查点；两个原生暂停记录及30→37 FAIL全部保留。

新启动预算10567.13324713707秒，原6小时扣除原执行、两次检查和第一续算1042.891990秒。31→37重复计算仍扣预算。技术依据RECOVERY_REVIEW_FINE30_0002.json、授权RESUME_AUTHORIZATION_FINE30_0002.json，正向预检及放宽门槛/漏算已执行段/重置预算三个负向控制均PASS，276个本地保护输入SHA一致。实际XML只改restart=true、仍500步、8MPI+CUDA，PETSc/物理边界/科学模型不改。83脚本快照0414b516cb4341fb3c2b5d93c2374b64a351f7f12670bcfbc070b1fa5685664e服务器当前及归档SHA/AST全PASS。

新入口scripts/resume_validated_fine.py、fine_resume_followup.py串行推进稳态、后处理、精确原生/规范场身份、全场审计、46初态单泡重放和分析。只有逐项成功才前进；不自动发布生产网格或进入dt/P2。report_fine_restart_recovery与render_report已支持多次暂停续算的预算和历史耗时，不隐藏上一续算成本。当前仍9个合格新CFD、10批回放；完整空间gate尚未完成。持续监控脚本monitor_progress.py；暂停时应先冻结新producer，再向当前新段run/STOP_SIM请求原生边界停机，并停止report_watch及后续队列。不要把旧段当成活动段。


## 2026-09-30T19:30:10.428855+00:00：细档稳态与全场审计关闭，46条回放正在运行

M1_coplanar_fine的segments/user_pause_resume_from30_0002已在第81步原生exit0，health_failures=[]，execution PASS；40/50/60/70/80五连续区间满足原稳态门槛。续算耗时4899.31183719635秒，原6小时预算扣除全部原执行、两次验证、两次续算后剩5667.82140994072秒。后处理与独立规范/原生输出身份核验PASS，规范flow SHA d049122a0e77c001a1d70e0b1ef09a88faae7b634b5b0bcb99e461cc55bf168a。原暂停execution FAIL与30→37 FAIL证明仍原样保留，不追溯改写。

全场审计已完成：1064977 tetra，严格零速单元0，1059固定点均在域内；CPU/GPU G最大绝对差9.095e-13 s⁻¹。此数值一致性不等于局部物理精度。中细档固定点u相对L2为1.9433%，G为19.3892%，WSS为20.9880%；部分局部区域u也超过2%，空间场门槛未过。原388769中心细档|div|约86.119 s⁻¹，η约0.037544，相比原场改善但区域梯度/WSS仍不稳定。

181个已关闭CFD/审计文件已本地回收并SHA PASS，见CLOSED_FINE_0002_LOCAL_COLLECTION_VERIFIED.json。当前合格新CFD计10套，已完成分析回放仍10批；第11批M1_coplanar_fine_dt0.00050000正运行，固定46初态/0.5ms/关闭泡间作用，不能提前称其完成或无失败。fine_resume_followup_0002与report_watch继续运行。原正式生产流场、1565条旧协议population和显示参数均未更改；dt/P2阶段未启动。

若用户再次暂停，原生CFD已经完成，不要重复从30步启动，也不要改写已PASS执行；先核对当前回放/分析安全边界。下一步关闭回放和分析、收集轨迹、复核完整空间gate，并按实际资源评估进一步mesh refinement。新增只读资源测算位于01_mesh_capability/COPLANAR_NEXT_REFINEMENT_FEASIBILITY，不生成候选/不运行CFD；25µm表面数组复现完全一致，20/18.75/12.5µm仅计算候选表面规模与非上界资源估计。


## 2026-09-30T19:50:01.761401+00:00：25 µm细档全部关闭，20 µm有界网格实验已启动

M1_coplanar_fine第81步稳态CFD、规范输出身份、全场审计及冻结46单泡回放/分析均完成。细档回放41出流、5观察截止、0数值失败、0有支撑静止；O1/O2/O3=16/17/8，截止ID5/393/534/1329/1377。20次固定阈值slow–fast事件（中档25次）。ID440 O2停留0.963458s，ID1505 O1停留0.845126s，ID636 O3停留1.800663s；ID1189 O3停留9.937649s，ID393/1329观察截止。原旧1565结果和协议身份均不改。

完整50/25µm比较FAIL：全部逐条门槛通过0/46；两档均出流但出口改变12条（123,168,175,341,487,542,650,708,1077,1208,1240,1505），终止类别改变6条（4,5,534,1329,1449,1522）。速度历史门槛46条均未通过；39条终点差、39条持续时间差、27条最小间隙差、11条路径长度差超门槛。固定点u/G/WSS相对L2=1.9433%/19.3892%/20.9880%，无正式网格选择，dt/P2仍未启动。不能把无数值失败称为空间收敛。

376个细档回放/分析文件回收SHA PASS，全部ID及终止记录交叉核验PASS，见CLOSED_FINE_REPLAY_0002_LOCAL_COLLECTION_VERIFIED.json。当前完成量为10套新CFD、11批46轨迹回放（包括同协议M0）。总报告19:39:59UTC版本地171链接/98图SHA PASS，REPORT_MANIFEST SHA d6f5173056b550fd793de532e026d99c9c6c8f4baa16542d7cba22d1754af116。新细档曲线已人工查看，无平滑/时间对齐；完整解释见COPLANAR_VALIDATION_REVIEW_ZH.md。fine_resume_followup_0002与report_watch均已正常EXITED，不能再次启动这些独占输出入口。

依据用户Phase2G继续可行空间细化，19:48:41UTC启动独立coplanar_h020_mesh_only。候选M1_coplanar_fine_h020，对照M1_coplanar_fine；只生成及QA，不自动CFD。表面272904面/136454点，原顶点和端口保留；显式表面准备资源上限30万，旧默认25万不变，所有科学门槛不变。原生TetGen已输出295691点/1412899 tetra，尚待完整拓扑与严格域对应QA和实际CFD资源复核；不可先计作新合格网格或CFD。

依据COPLANAR_H020_MESH_REVIEW.json（SHA49dad3a905716e457dc8c05c27fec4fa4418f9183abcc235b8506fd174982ea3）、COPLANAR_H020_MESH_PLAN_ZH.md；资源预测及25µm输入逐数组复现见COPLANAR_NEXT_REFINEMENT_RESOURCES_ZH.md。prepare_coplanar_plc_probe只新增显式表面资源上限参数，audit_coplanar_transport_mesh新增独立plan入口及源记录；新coplanar_mesh_probe_followup按源/证据锁定串行执行并禁止自动CFD。84脚本快照fd04695f582befc62430fe18de1932bb32e477abc74b69a005d59d6bff50ee3a当前/归档AST及SHA全PASS。旧83快照及旧执行源码保留，粒子/CFD物理未修改。

当前活动任务只有coplanar_h020_mesh_only。如收到新的暂停指令，先核对该网格队列原生生成/QA边界，不要向已关闭CFD写停止标记或从30步重复恢复。网格阶段结束后检查实际质量、规模和内存再决定下一步；20µm未运行CFD，不存在它的流场或轨迹结果。


## 2026-09-30T20:31:04.066971+00:00：20 µm网格通过，新的原始零场CFD正在运行

最新用户指令仍为恢复开发。20 µm候选M1_coplanar_fine_h020已完成网格及严格域对应QA：296162节点、1414875四面体、全WALL零速风险单元0，负/零体积和非流形面0；原分片平面壁面和4端口严格保留。原生Float32 probe的FAIL证据保留，输入精度恢复及471次完整内部边扇分割后最终GATE_MESH PASS。minSICN=0.01165949，比25 µm的0.01811更低；不能只因更细就声称更准。46个封闭网格文件本地SHA核验PASS，见CLOSED_H020_MESH_LOCAL_COLLECTION_VERIFIED.json。

资源复核、85源文件和30证据哈希预检及预算/漏源/错网格/放宽门槛四项负向控制均PASS后，coplanar_h020_cfd_followup于20:18:56.261662 UTC启动。实际case为03_cfd_runs/M1_coplanar_fine_h020（根case，从零场开始，不是旧25 µm续算段）。实际XML/PETSc与正式基线逐字节相同：8MPI+CUDA、P1/P1、dt=0.01、500步、21600秒、5连续稳态区间；压力及微泡模型未变。当前仍首次CPU ASM/LU分解，尚无合格h020流场或轨迹。20:29实际求解器RSS约34.44 GiB，cgroup总42.05 GiB含约7.61 GiB可回收文件缓存，OOM/oom_kill均0，继续监测实际峰值而不把预测当上界。

当前85脚本快照SHA ca522c5f7e9b654306c4d5b49302c85ac2283e41b82ea9ca9dbeb7120d25da33；CFD review SHA f248f5c8092ddcaae934791a489e6f4ae151a2ac5819853ec3bf95020475efb2。见COPLANAR_H020_CFD_REVIEW.json/.md（可读文件为COPLANAR_H020_CFD_REVIEW_ZH.md）、H020_CFD_PREFLIGHT_CONTROLS.json、H020_CFD_STARTUP_LOCAL_VERIFIED.json。启动配置10文件及12个输入哈希已收回本地核验。当前运行producer引用的85脚本不得运行中修改。

队列仅在CFD稳态PASS后顺序执行后处理、独立规范/原生输出身份核验、全场审计、冻结46初态0.5ms回放、分析及空间比较。当前完成科学量仍10套新CFD/11批回放，25 µm已全部结束，50/25 µm完整空间gate仍FAIL（逐轨迹全门槛0/46）。无新生产网格、无dt/P2/无散度试验。20:16:49报告（启动前静态快照）已本地验证171链接/98图SHA PASS，manifest SHA dbd00dfdc05c7ed80afb52cf56ee704df323fa5fa3056f7cdb065ec34e30fb06；其active_CFD=[]只代表生成时刻，实际运行状态见monitoring/LATEST.json与CURRENT_EXECUTION_STATE.json。

只有coplanar_h020_cfd_followup活动，fine_resume_followup_0002/report_watch/coplanar_h020_mesh_only已EXITED，不能重启独占输出旧入口。若收到新的暂停指令，应核对当前阶段并冻结后续推进，向活动h020根case/run/STOP_SIM请求原生安全边界退出，保留检查点与已耗预算；不能用旧25 µm暂停辅助原样操作，不能恢复已关闭25 µm或重置6小时。新的h020检查点须另行审查连续性，旧25 µm证明不可套用。


### 2026-09-30 20:53 UTC 活动监控补充

20 µm CFD已完成第1时间步（14次非线性迭代，末Ri/R0=4.530e-11）并进入第2步；首步观察原文与哈希见H020_FIRST_STEP_OBSERVATION.json。两次CPU LU构建已通过，第二次由ITERATION_GROWTH触发；近时cgroup总内存约38.68GiB，solver进程树RSS约37.42GiB（共享页可能重复计数），OOM/oom_kill=0。未有稳态资格、未有h020流场或轨迹，继续coplanar_h020_cfd_followup，不重启旧任务。

只读原生日志复核：SV13Q age在同一时间步的多次logical求解间保持不变，进入下一时间步才增长；例如25 µm原段logical3至14均step1/age0，logical15为step2/age1。因此不能把max_age=5口头写成“每5次线性求解必重建”；还存在ITERATION_GROWTH等提前触发条件。以上是实际日志观察，不修改任何策略。

待当前85脚本producer完全退出后再部署的本地草稿在pending_development/：report_coplanar_h020_validation.py生成独立25/20 µm比较及6张PNG/PDF；report_fine_restart_recovery.py修正已合格25 µm仍显示“尚需完成”的旧文字并列实际预算。均只做报告，不改科学门槛或结果，尚未部署/运行。两个brava_collect_closed_h020_*.py是对应封闭阶段的独占回收脚本，已做语法检查，未执行；执行前必须核对阶段和manifest不存在。旧25 µm回收脚本不能直接用于新case。


### 2026-09-30 21:22 UTC 首小时资源图已核验

20 µm CFD已完成前3时间步，当前第4步迭代，3次预条件器构建已完成，约38.8GiB总内存且OOM/求解健康告警均0。首小时20:18:56至21:18:56的60份监控快照已封闭保存为monitoring/H020_FIRST_HOUR_STATUS.jsonl（SHA9fcf480688ed8bdf285e9aaba811d30beeca4d6a169e241935a39dcb5559b936），本地只读生成资源PNG/PDF并人工查看，7个文件本地/服务器SHA PASS见monitoring/H020_FIRST_HOUR_RESOURCE_ARTIFACTS_VERIFIED.json。资源脚本位于monitoring/report_h020_first_hour.py，不在运行中85脚本producer目录，不改变其身份。

图文入口monitoring/H020_FIRST_HOUR_RESOURCE_REVIEW_ZH.md。CPU区间加权7.680012核（配额7.68），采样cgroup峰值42.055626GiB（限制42.055664），OOM/kill增量0；包含文件缓存/共享内存，不等于solver私有RSS。GPU是稀疏瞬时样本，不推算平均利用率、时间占比或完整瓶颈。此图只涵盖首小时，不声称完整CFD内存上界、稳态资格或输运改善；主报告98图manifest未因此修改。


21:46 UTC：H020_RECORDED_POSITION_WALL_QUERY_PLAN.json：现有audit_coplanar_wall_queries.py支持直接使用25 µm已关闭46轨迹的全部记录位置，依次查询M0/25/20 µm壁面；本计划冻结源/轨迹/边界哈希，等待活动producer完全退出并取得重计算锁后执行。原50 µm参考查询用时11.76秒、CUDA张量峰值约0.045GiB（不代表CUDA上下文或进程全部显存/主存）。未启动本计划，未增加任何轨迹积分。25/20 µm报告草稿已支持新查询审计路径COPLANAR_FINE_H020_PATH_WALL_QUERIES；当前草稿SHA70131730e1321c29286894cd07f2a96ca8e2161d843ca7c684fce6c925607db6。 当前h020已完成7步，第8步因MAX_AGE启动第4次预条件器构建，暂无异常。


## 2026-09-30 22:11 UTC：20 µm首个保存状态通过完整性核验，CFD继续

第10步已完成并进入第11步。原生stFile_010.bin的8-rank Y/A检查及result_010.vtu有限性、壁面无滑移核验PASS，两个原文件及证明已回收本地SHA PASS，见03_cfd_runs/M1_coplanar_fine_h020/reports/SAVED_STATE_010_VERIFIED.json与SAVED_STATE_010_LOCAL_COLLECTION.json。检查点SHA e1e56de0118089ab487ab5dfdb70edde22f9832ce18f78c404e68c5046c2edf8。这只是保存状态完整性，绝不声称重启连续性或稳态已通过。

第10步time=0.1s：入口18mL/min，出口归一化分流0.3332182392/0.3331620949/0.3336196659；epsilon_mass=0.0009765625，尚大于1e-6稳态物理门槛。压力0.439712…2166.73231Pa，umax0.54997159m/s，壁面最大速度0。尚无可比较的保存间隔，连续稳态区间0/5。输出适配器读取无重试或失败；实际后处理/全场审计/微泡回放尚未开始。

活动任务仍仅coplanar_h020_cfd_followup，实际case仍根case，从零场连续推进。总内存约39GiB、OOM=0；4次预条件器构建均完成。实际源85快照不变；未修改物理、容差、预算或生产场。完成量仍10套新CFD/11批46回放。

后续保存状态回收辅助为/tmp/brava_archive_h020_saved_step.py（副本pending_development/brava_archive_h020_saved_step.py），参数10至500中的10倍数；等待原生进入下一步才读，首个状态interval=null合法。10已执行且输出独占，不要重跑。以后20/30等可独立回收，但这些证明都不是重启连续性认证。


22:47 UTC：20 µm CFD仍活动，当前第16步，CPU约7.68核、cgroup约39.1GiB，health/OOM/read retry均0。新增pending_development/render_report.py草稿，将25/20 µm独立报告、六张PNG/PDF和HTML链接接入主报告；AST已通过，三份待部署report源码本地/服务器SHA一致，身份见H020_REPORT_INTEGRATION_DRAFT.json。草稿尚未部署或执行，活动85源快照ca522c5f7e9b654306c4d5b49302c85ac2283e41b82ea9ca9dbeb7120d25da33保持。brain_report第15节顶部已更新为10套新CFD/11批回放已完成、20 µm正在运行，原历史时间线保留。下一保存状态20步尚未产生，10步不重复回收。


## 2026-09-30 23:08 UTC：用户要求短暂停止汇总后立即恢复（尚在等待保存点）

最新指令为“在下一个可停止点暂时停止，并输出已完成/结果/后续及预计耗时；完成这些后恢复开发并持续显示进度”。这不是等待以后另发恢复指令的长期暂停。本次计划在h020第20步保存并由brava_archive_h020_saved_step.py完整核验/回收后，对当前同一组12个producer/runner/MPI进程做短暂SIGSTOP，生成阶段报告后SIGCONT原进程，不退出求解器、不重启检查点。暂停时间仍计入原6h预算。尚未发出任何真实停止信号，也未写STOP_SIM。若用户另发新的暂停/休眠指令，优先遵守新指令，不沿用自动恢复。

控制工具monitoring/brief_status_hold.py，--stage suspend/resume --study RS --record RS/monitoring/H020_BRIEF_STATUS_HOLD_20260930.json --checkpoint20；默认checkpoint参数值20（实际CLI写作--checkpoint 20）。必须先有SAVED_STATE_020_LOCAL_COLLECTION.json。源SHA7b3d1c1119645e5365f9cffa5d78487c683e959b24e716239c0cf68eb29439d3；缺保存点拒绝操作的负向控制PASS，H020_HOLD_PREFLIGHT_CONTROL.json。工具逐PID验证start_ticks、T状态和CPU不推进，resume验证配置/检查点SHA与同PID，绝不对管理服务发信号。当前producer759659、MPI759669，8原生rank759676..759683；实际工具仍从supervisor和/proc重新查身份。

报告生成工具monitoring/write_status_review.py，要求真实hold记录及第20步证明后才执行，输出根目录PROGRESS_REVIEW_AND_NEXT_STEPS_ZH.md/.json。先暂停时生成，再恢复后刷新状态；尚未执行。已准备完成内容：10套新CFD/11批46回放、50/25µm FAIL0/46、20µm未合格及后续gate；当前阶段估计4–6h，预算触发则额外1–2h或更多；空间通过后dt/整理另2–4h，需P2则规划1–3工作日以上且无实测承诺。ETA证据monitoring/H020_ETA_EVIDENCE_20260930T2304Z/：近期普通非线性迭代中位26s，step8/13重建首条耗时685/689s；按25µm未来NL数及每5步重建，从h020第18步开始到81步约14469s。仅规划情景，不改预算、不保证81步稳态。

23:07 h020已完成第18步4次非线性迭代，预条件器第6次构建完成，内存约39.1GiB，无health/OOM告警。继续监控至20保存点，完成实际短停/汇总/恢复，然后继续原开发，不要因此发最终答复结束仍活动的计算。


## 2026-09-30 23:19 UTC：短暂停止、汇总、自动恢复均已完成；当前继续开发

用户23:00要求的“下个安全点短暂停止并汇总，然后恢复开发”已全部执行。h020第20步stFile_020.bin完整8-rank Y/A及原生流场核验PASS并回收本地；检查点SHA cdbc69781439fc20fe0c1d1c1303ab5796049a3909ca8914ec95d8ed03d9f16a，证明reports/SAVED_STATE_020_VERIFIED.json与SAVED_STATE_020_LOCAL_COLLECTION.json。step20质量误差9.536743159e−7、wall no-slip PASS；10→20区间Eu=0.0006024017628、EQ=0.0005338033864，仍高于1e−5/1e−6，连续稳态0/5。10/20已回收且独占，不重复执行collector。

在完整step20后挂起同一组12进程45.202254秒，23:18:02.622881 UTC挂起，23:18:47.825130 UTC恢复；当前MPI及8rank PID/start_ticks未变，日志已继续增长。没有native退出，没有检查点重启，没有STOP_SIM，没有改写budget或科学配置，短暂停顿计入原21600秒。审计monitoring/H020_BRIEF_STATUS_HOLD_20260930.json与H020_BRIEF_HOLD_RESUMPTION_VERIFIED.json，恢复及3个报告文件本地/服务器SHA一致。后者status=PASS_RESUMED_SAME_PIDS_LOG_ADVANCED。当前状态恢复RUNNING_REVIEWED_COPLANAR_CFD；不是PAUSED，也无需再等恢复指令。

用户要求的完成量、结果、后续和时间估计已通过commentary答复并写入PROGRESS_REVIEW_AND_NEXT_STEPS_ZH.md/.json（9个链接存在，SHA验证），已同步服务器。MD SHA1101e49cdab0f5284761a1c8520ac6b35ca196356e8ea8d1930d627c204ca919。当前阶段估计4–6h，空间通过后全剩余工作约6–10h；若预算触发并需重启验证，额外1–2h或更多；若需P2代表域则规划1–3工作日以上，非保证时限。ETA实测证据留在monitoring/H020_ETA_EVIDENCE_20260930T2304Z/，预算可能不足仍未修改。报告生成器首次使用相对--hold路径报ValueError，随后改用绝对路径运行成功；科学结果不受影响，实际使用时--hold传绝对路径。

23:19监控第21步已完成2次非线性迭代，CPU区间平均包含45秒挂起故短时4.54核，不能误判性能退化；恢复后的独立区间待继续采样。内存39.18GiB、health/read retry=0。继续原producer；完成量仍10套新CFD/11批46回放；h020未稳态、未回放。开发和持续进度指令仍生效。


23:48 UTC监控：h020已完成第26步6次非线性迭代；第23步第7次MAX_AGE预条件器已完成，内存约39.1GiB、CPU约7.68核、health/read retry/OOM均0。用户短暂停止并汇总后恢复的要求已在23:18完成，不要再执行brief_status_hold或等待新的恢复指令。下一个待回收完整保存状态30；10/20已独占回收。

待部署25/20µm报告草稿已在23:32补入本次45.202254秒同进程挂起/继续、预算包含停顿的原始记录，当前pending_development/report_coplanar_h020_validation.py SHA54545916905a8f6f9f9d1ae8286cea975fc7f88189d0d26112a3f871413bac09；旧701317草稿保存在pending_development/previous_drafts。H020_REPORT_INTEGRATION_DRAFT.json已更新，脚本AST通过并同步；尚未部署/执行，活动85脚本仍不变。

只读复核B/scripts/solve_remote.py：policy在启动时读入一次；达到maximum_wall_time_s且stop尚空时写STOP_SIM，等待原生边界保存退出；未见硬kill分支。budget退出没有stop.step，因此保留execution FAIL。当前21600秒未改，未创建额外预算/重启验证/续算计划或作业。实测ETA存在先触预算的风险，若实际发生，保留原native/execution/检查点和全部成本，再按用户已允许的“随DOF记录增加墙钟 allowance”审查独立续算；不能靠改policy文件假装改变活动进程预算，更不能改时钟/内存或追溯改写PASS。原25µm重启证明不能用于h020，旧fine25入口不能原样复用。


2026-10-01 00:07 UTC：单次本地SSH监控查询返回255，原生CFD未中断；00:05直接读取原生日志确认28步已完成、29步迭代，00:06:50常规监控恢复为29步NL5、health=0、OOM=0。没有向求解器或管理服务发信号；尝试清理滞留本地查询前该PID已自行退出。/tmp/brava_monitor_live.py已加连接8秒/存活探测及总体25秒超时，持久副本monitoring/local_helpers/brava_monitor_live.py，SHA48eb3faba10d49be9bec6ef3fa3c8a4938d861b596ab7fd327ecbb38375595d9；旧版归档同目录。活动85服务器脚本不变。记录monitoring/MONITOR_QUERY_RECOVERY_20261001.json。继续至30步保存并回收，不再次执行已完成的短暂hold。


## 2026-10-01 00:13 UTC：第30步检查点已收回；启动收敛图完成

h020第30步已完成、当前第32步NL2，health/read retry/OOM均0，内存39.18GiB，唯一活动producer仍coplanar_h020_cfd_followup。第30步stFile SHA cced807e052650ca93e87d35c1306ae6b07346aa18c3f99a566af54ac44d42b2，size20209024；8rank Y/A、流场有限性、无滑移和本地文件SHA均PASS，local collection 2026-10-01T00:10:20.304084UTC。服务器证明SHA5cd6c5d2122f5117d1049f7c3143c9e68dcab400f3a1f7a721d8bc9a53fb0016。10/20/30均已独占回收，不重跑对应collector；下一个40。此为保存完整性，不是重启连续性或合格流场认证。

20→30：Eu=7.0319738117395925e−6低于1e−5，EQ=1.1857101942204126e−6仍高于1e−6，连续稳态0/5，未改门槛。mass=9.313226075571546e−10、WALLu=0，umax=0.5499245861335055m/s。第30步pressure_range_pa=[0.4384423082164724, 2725.8372827409903]；outlet_fractions={'OUTLET_01': 0.33308154594874734, 'OUTLET_02': 0.33337173066055575, 'OUTLET_03': 0.3335467233906969}。

新增只读启动收敛图：monitoring/H020_STARTUP_STEP30/REVIEW_ZH.md、REVIEW.json、startup_convergence.png/.pdf。本地标准matplotlib生成、已目视检查，绑定三个完整保存状态及其本地回收证明；第10步缺相邻区间，保留缺失不补零，图不是合格稳态场。5文件本地/服务器SHA PASS见monitoring/H020_STARTUP_STEP30_ARTIFACTS_VERIFIED.json，源monitoring/report_h020_startup.py SHA580453e79dc8f7ef76ec5dd1402d06972623abd2734b47a872954cca263af490；PNG SHA272dbe1d5f137625ef3632aaaf0d9c004964e903a07fde4df47553a67155b93c。图在monitoring，不修改主报告98图manifest或活动85源快照。

用户23:00要求的短暂停止/汇总/恢复已完成，当前继续开发，无待确认事项。完成科学量仍10套新CFD/11批回放，h020未稳态；原21600秒预算未更改。当前进度报告PROGRESS_REVIEW_AND_NEXT_STEPS_ZH.md保持23:18阶段快照，最新进度读monitoring/LATEST.json。


## 2026-10-01 00:41 UTC — 20 µm继续运行，预算恢复仅为已验证保护的备用草稿

- 最新用户的“短暂停止→输出完成项/结果/待办/ETA→自动恢复”已于23:18完成（原PID暂停45.202254秒后恢复）。当前开发继续，不再等用户恢复；若收到新的暂停指令则按新指令执行。
- 当前仍只有 coplanar_h020_cfd_followup 运行，原6小时/500步时限和85个源文件不变。00:41监控：原生第38步预条件器第10次重建；37步完成，累计15730秒，CPU约7.679核，内存约37.06 GiB，健康标志全0，OOM/oom_kill=0。memory.max累计事件从746607增到747602，不能把累计历史事件全部归于本次；本次观察间隔新增995次，pressure累计新增约0.017秒。仍0/5稳态区间；10/20/30已收集，下一次40。
- 新增pending_development三个备用草稿（未调度、未启动、未修改活动代码）：recover_h020_wall_budget.py (91b8532ebc826527cded97932e88ba952704e7584e29bad99621d2f4f56f63cc)、h020_wall_budget_followup.py (dd513c9458aa790d10491ff27c02ae85139e0b199d8498fb42f3b18d4ed1cfc8)、review_h020_closed_budget.py (685cae6c74638bef9a7d8d04866aa3c87127034258b57d6201cacf6f0553fff3)。均AST PASS，在真实活动状态下执行保护检查均拒绝：原execution.json未关闭，未创建审查或新段、未启动CFD。
- 方案见pending_development/H020_WALL_BUDGET_RECOVERY_PLAN.json及同名_ZH.md；本地/远端8个草稿/证据SHA一致，证明H020_BUDGET_DRAFT_FILES_VERIFIED.json SHA916da85b22ce2dab74d34bbf19ac26b88d4dba5eb653d846065e42f0cb4076c0。正向预检和原生连续性尚未执行。
- 原用户附件Phase2明确允许随DOF增加记录后提高wall-clock allowance。备用拟增加6小时使h020累计原生额度12小时；这不是用户指定数值，也不是ETA/收敛保证。必须先等原producer完全退出，确认原本段确为预算FAIL、exit0/health[]，保存关闭原始资料、核对原生最终checkpoint，才运行review_h020_closed_budget.py --study <server study>生成独立H020_WALL_BUDGET_REVIEW.json；现在该文件不存在。严禁修改活动policy或绕过其STOP_SIM。
- 具体恢复方案用h020自身20→30重放，1小时上限、原u/p/各rank Y/A相对L2门槛1e-8；仅PASS后从原h020第30步续算，并重新满足5区间稳态。全部原段/验证/续算耗时累计扣除12小时，原FAIL不改写。25 µm的连续性证明不转用。失败不自动重试。
- 闭合审查创建后还须先实际正向preflight并核对资源，再调度pending_development/h020_wall_budget_followup.py --study <study> --review <absolute review>；它依次检查/续算/后处理/独立身份/场审计/冻结46条回放/分析/原有报告。未创建新的supervisor job。若初始段本身PASS，则完全不启用此备用路径。
- 已准备的/tmp/brava_collect_closed_h020_cfd.py仍只接受fresh root PASS；若使用续算路径必须先适配segment和本地SV_MESH兼容链接，不能直接运行旧草稿。报告整合草稿仍要等producer闭合后部署；85脚本当前禁止变更。


## 2026-10-01 00:58 UTC — h020第40步首次通过区间门槛，继续运行

40步已完成独立核验/本地回收（不再重复执行collector40）。原生cp SHA a0bc9fad418af43e2177c38a32775de7ec5dda876352a3911a5495a7b34f5735。30→40 E_u=1.4440512382013205e-7、E_Q=3.6555320171826055e-9；mass=9.093251619169604e-13，壁面u=0，u_max=.549924585365713 m/s。连续通过1/5，仍无合格h020场/轨迹。00:57:57原生进入43步重建第11个PC，原段elapsed16737秒，CPU7.677核/RAM35.05GiB，health全0/OOM0，memory.max累计747602未再增长。原6h政策不变，备用追加预算仍未激活。

新科学图表 monitoring/H020_STARTUP_STEP40/REVIEW_ZH.md + REVIEW.json + startup_convergence.png/.pdf 已运行、看图并双端SHA验证。证明 monitoring/H020_STARTUP_STEP40_ARTIFACTS_VERIFIED.json SHA5a9a35c16ac86b0413c87d302fdb921d6db3d78e2744536d30744c92dd89f795。新通用绘图程序 monitoring/report_h020_startup_series.py SHAf88b1c2e7ab38ee3d4b266e16a2eb6ac40e95e69ed89168556650189b401e5bf，可在各保存状态收集后用 --through-step 50/60/...生成新独立目录（已存在则拒绝覆盖）。原生比较为<=，图标明连续区间数，不单独认证稳态或局部物理精度；不修改原30步图。

时间估计已随实测更新：monitoring/H020_ETA_STEP40/REVIEW_ZH.md及ESTIMATE.json/monitor_snapshot.json，双端证明SHA cfe1f951967eced22a3448ead78087217f2917ab99594538f2f2ed7658f0625b。00:53快照原6h剩5156秒；若假定50..80均通过，到81步不中断尾段约7472秒；因此很可能预算停机。备用h02020→30验证约3090秒、原30→81续算约10677秒，加原段剩余后原生总剩约18923秒；含审计/46条回放/报告，本轮预计仍需5–7小时。假定未来NL次数像已关闭25µm、普通迭代26秒/重建额外661秒，非收敛/时限保证。若仍需P2/P1，全项目预留1–3工作日或更长。当前用户未再暂停。

新的通用关闭CFD回收草稿 pending_development/brava_collect_closed_h020_cfd_any_segment.py，/tmp同名，SHA3aeefd66b49d89f2423780d606be00cca7792188941fe899c85cbcc51e205e11；外层及嵌入远端probe均AST PASS、已同步但未执行。兼容fresh PASS或经过独立验证的预算续算PASS，校验原FAIL保留，回收全部原生段；只把已知远端绝对SV_MESH链接转为本地相对链接。优先使用此新collector；旧root-only版本保留。下一保存状态50，预算到点前继续原producer；不要把备用草稿当已运行或已正向预检。


## 2026-10-01 01:10 UTC — 备用恢复草稿输入保护补充（仍未启用）

pending_development/recover_h020_wall_budget.py更新为SHA6d66846b4504e10ab5eccb3541d4ee39125412caf3f85e2f9a08f287a62c7c42：在每个验证/续算阶段重新核验BUDGET_STOP_PRESERVED内所有原段文件（包括原生分区文件），并拒绝未计入账目的额外分段。其余followup/review_creator源码不变。旧9项草稿/保护证明完整保存于pending_development/previous_drafts/h020_budget_20261001T0108Z，两端均已归档。更新后的三个真实活动状态负向保护均PASS，8项草稿双端SHA一致，新的H020_BUDGET_DRAFT_FILES_VERIFIED.json SHA3adbc16d02651ca5bdbab32bc634e9158bf9ebed6950cfe07316d5a01bec1fd0。85个活动源文件、原6小时政策未变；正向闭合审查、重启连续性、续算均未执行。01:09:34主计算到46步、原生elapsed17431秒、health全0、仍1/5。若以后进入分段，/tmp/brava_h020_memory.py仍读取根case的progress，需要先适配RESUME_INDEX再用于新段资源监控；普通brava_monitor_live.py已经能识别分段。


## 2026-10-01 01:26 UTC — h020第50步已归档，连续稳态2/5

50步collector已执行成功，不重复：cp SHA0d8ef7f82dfcfb7770f8d8fadc218793cf0df98a0d7093e4720752ae36218bab；40→50 E_u=3.256860679840286e-9，E_Q=1.3692905220556472e-11，mass=7.058607893785836e-16，wallu=0，u_max=.5499245853645829。累计2/5，仍未正式合格/未回放。monitoring/H020_STARTUP_STEP50下MD/JSON/PNG/PDF已生成、开图检查并双端SHA验证，证明monitoring/H020_STARTUP_STEP50_ARTIFACTS_VERIFIED.json SHAef871879fa0ab24a854880d8b11ef5427f8237a690dc68b3f26a20728552e2ca。01:25:31监控进入52步，elapsed18391.804秒、CPU7.6805核、RAM39.164GiB、health全0，原6h预算仍有效。

下一保存节点60；10/20/30/40/50都已独立收集。当前用户要求的短暂停止/总结/自动恢复早已完成，开发仍ACTIVE。备用预算恢复代码版本以01:10记录为准，仍未激活；正向review_creator须用/root/particle8_2_runs/env/bin/python（负向测试在缺少closed execution时已提前拒绝，未验证科学依赖），先等原producer EXITED、原预算FAIL关闭且HEAVY锁可用。原runner只写STOP_SIM并等待原生步完成；若6小时点落在耗时约11分钟的PC重建中，实际关闭可能超出时限数分钟，全部真实耗时仍计入累计12小时。不能因deadline到达就假定已保存退出。


## 2026-10-01 01:58 UTC — h020第60步归档，连续稳态3/5，原预算将到期

60步collector已成功执行，不重复：cp SHA9f8b36202e705b29d2b4f5242fc06e5f09c911121508fa09cc5e03e21988207c；50→60 E_u=6.920730656668296e-11，E_Q=5.461597857816791e-14，mass=1.764651973446459e-16，wallu=0，u_max=.5499245853642368。累计3/5，不是正式合格场。10/20/30/40/50/60均已独立回收，下次70（也可能先触发6h预算停止）。01:57:37已进入63步重建第15个PC，elapsed20318.287秒，health全0，原6h剩约1282秒；不能在deadline到达时假定原生已经退出。

monitoring/H020_STARTUP_STEP60的MD/JSON/PNG/PDF已生成、看图并同步验证。最初排版的边缘标签靠近y轴刻度，初版完整保存至H020_STARTUP_STEP60_INITIAL_LAYOUT；新monitoring/report_h020_startup_series_v2.py SHA4042a972e8f8060ac8d1283234c7c5a7a2a73b976f753cd412de50427d596911只调整首末标签对齐，输入哈希及所有数值与初版一致。40/50图和旧程序f88原样保留。9文件（当前5+初版4）双端证明monitoring/H020_STARTUP_STEP60_ARTIFACTS_VERIFIED.json SHAac58814618cd631e1aa2fe87defb1c895c03e94b94fec9683c1afa229b8c082b。后续绘图优先v2。

最近memory.max累计748084，较747602增加482次，发生于58步重建期间；pressure累计增加约6.5ms，OOM/oom_kill仍0。峰值不是局部物理精度证据，CPU/GPU资源策略和源文件保持原样。当前用户没有再次要求暂停；开发ACTIVE。备用预算方案及三个程序仍未正向执行，使用01:10版本driver6d66846b...和followup dd513c...，创建具体闭合审查必须先见原producer EXITED、原执行budget FAIL/exit0/health[]及完整final checkpoint。建议闭合审查生成后先回收原rootcase/审查并逐项校验BUDGET_STOP_PRESERVED文件哈希，再正向preflight/调度新job。不要直接运行只接受CFD PASS的collectors来回收这段预算FAIL；本地原始元数据需要保留，然后再启动恢复。


## 2026-10-01 02:32 UTC — 原h020预算段已关闭，独立20→30连续性验证正在运行

原coplanar_h020_cfd_followup已EXITED。原生最终68步、exit0、health[]、checkpoint PASS；status FAIL仅因为wall time budget exceeded，原6小时停止请求02:19:02.040367 UTC，实际原生wall_time21918.14679002762秒（6h05m18s），等待当前PC/步完成的额外318秒如实计费。原execution SHA8f9d30ca32a002e0d00b217e19e5e9e026228e0ea7e6d6ce6ba8537af0f8e970；68步cp SHAf5d67eaea2ad635088ff18ca615b21e1bf06c9b087d055c8aaacc968879a5bfa。原最后完整稳态保存60步为3/5；未资格化、未回放。原20/30/40/50/60和68数据不覆盖。

已用solverpy运行review_h020_closed_budget.py生成H020_WALL_BUDGET_REVIEW.json SHA362c0507e7fe7a8ac58b4b84653d91000e64ad8555dbe99b878912fd0f6daa67。原47项关闭文件已回收并逐项SHA通过：H020_INITIAL_BUDGET_LOCAL_COLLECTION_VERIFIED.json；原元数据单独归档monitoring/H020_INITIAL_BUDGET_CLOSED_METADATA，包含原身份/policy/input_hashes/execution/progress/保存清单和恢复review。保存清单BUDGET_STOP_PRESERVED SHA ddc37d844e6ee89c315d8f75c988dd9dff751dac85e0ffc15e2891c935ffcb12。原生/源文件/网格均未修改。

正向preflight PASS：H020_WALL_BUDGET_POSITIVE_PREFLIGHT.json SHA7746651eaca38fb4961217b8a855f3ca9ac143dd951d4c3b6e60bf761163a90d。新额度累计43200秒，追加21600，扣原实际21918.14679后在验证前剩21281.85320997238秒。所有验证/续算实际耗时继续累计，原FAIL不改。中文review H020_WALL_BUDGET_REVIEW_ZH.md SHA220676f3f6e837583e9995442d7fc9f5139f27ebce5bfc733f66a37808e92ece。

已通过原scripts/add_job.py调度唯一新任务coplanar_h020_budget_followup_0001，producer PID776396；H020_WALL_BUDGET_FOLLOWUP_JOB.json SHA3fcb2ef12f7ccd6fd74ed921659ce595e97561d8becfd3e05308d656a6df48a2；启动证明H020_WALL_BUDGET_LAUNCH.json SHAdfc37c22009e6102f784f01813c2ec0b8daf699053e1d2a4c25ce4dec757303c。CURRENT=VALIDATING_H020_BUDGET_RESTART，RESUME_INDEX指segments/wall_budget_check_20_to_30_0001，MPI PID776416。先20→30原生重放，上限3600秒，u/p与每rank Y/A原1e-8门槛。仅PASS才由followup自动从原30步继续；失败直接停止，不能自动重试/放宽/转用25µm证明。30→最终续算、审计、46回放均尚未运行。初步资源约7.68核，验证段第一PC构建中，health全0。

代码继续使用driver6d66846b4504e10ab5eccb3541d4ee39125412caf3f85e2f9a08f287a62c7c42、followup dd513c9458aa790d10491ff27c02ae85139e0b199d8498fb42f3b18d4ed1cfc8；活动85脚本保持原快照ca522...，直到新producer也完全结束前仍不可改动。新分段所有文件独立保存，根case仍原预算FAIL。检查段即便原生raw execution FAIL（30步到顶而非稳态停止）也不能单凭该字段认定连续性失败，要读独立RESTART_CONTINUITY_BUDGET_20_TO_30.json；validate会严格校验最终步、全部实际步、health及原1e-8误差。

资源工具/tmp/brava_h020_memory.py已经按RESUME_INDEX选段，当前SHA93ac06843284a86fb1b6bfce4648091ea5808e1a796595edc86120a8ae7d0cc4，原版本和一次首键解析失败版均归档monitoring/local_helpers。一次远端查询成功但本地原解析器假定utc首键导致StopIteration，已改为与JSON键顺序无关并验证9个新原生进程；无CFD中断。记录monitoring/MEMORY_MONITOR_SEGMENT_UPDATE.json SHAfe995fbeb94be1eebc64b73ddf5d6327834ba593b7bde9bcca4dc1d5669fca03。普通monitor_live已经选新段，初次低CPU均值包含前几分钟回收/预检空档，之后7.68正常。不要继续用旧根progress判断新段耗时/进程。

后续：监控新验证与自动续算；验证关闭后独立回收校验，别重跑root已归档10..60 collector来收新段。任何新用户暂停优先按新段实际状态处理。待最终producer闭合后再运行wall query、部署报告整合并做25/20空间审查。pending报告report_coplanar_h020_validation.py最后关于“没有从检查点重启”的句子原指45秒短暂停顿；正式报告必须缩窄该句并补充此后的预算恢复执行链、原FAIL、连续性和累计耗时，避免误称整个h020从未重启。collector使用any_segment版本3aeefd...，原root-only已不适用。


## 2026-10-01 02:39 UTC — 最终h020报告草稿补充预算重启执行链

pending_development/report_coplanar_h020_validation.py当前SHA710cc43f541badecab70616a980a6130ad3f75b75c48903ae1243da8b1f88e17，旧545459...保存在previous_drafts同哈希名源码中；H020_REPORT_INTEGRATION_DRAFT.json已更新。新增原预算FAIL保留、同网格20→30原1e-8连续性、最终合格段与canonical identity、所有执行SHA/耗时账目的交叉核验；缩窄45秒短停未重启的措辞为“仅就那次短停而言”，另列后续预算检查点重启表。尚未部署/完整运行，85活动源码未改。真实远端VALIDATING_RESTART状态下，单独执行报告身份保护函数会拒绝宣称合格且不写报告；证明H020_REPORT_LINEAGE_DRAFT_GUARD.json SHA8e65afb909fc8bbad8260f2e22c310b9ccc061e4861e552e075212bcaa18aa2c。正向真实闭合场报告验证仍待最终CFD/回放结束。02:38:28验证段elapsed586.635秒，仍第21步FIRST_SOLVE PC构建，health全0、CPU7.6806核、RAM35.93GiB。新段自己的稳态计数0不覆盖原段3/5，也不等于连续性失败。


## 2026-10-01 02:58 UTC — 连续性验证第26步，补充闭合原段资源证据

开发仍ACTIVE；23:18短停/进度汇总/自动恢复已履行，无新暂停待处理。唯一活动job仍coplanar_h020_budget_followup_0001，MPI776416，20→30验证段。02:57:38累计1737秒，25步7次非线性迭代完成，26步按MAX_AGE重建第二个PC，尚无新残差；CPU7.682核，health全0，OOM/oom_kill=0，memory.max累计749293未新增。36e2秒验证上限、原1e-8连续性门槛、累计12h账目均未改；原FAIL保留。原运行85源码和活动pending恢复三程序均不修改。

新增只读资源报告01_mesh_capability/COPLANAR_H020_CLOSED_INITIAL_RESOURCES/REVIEW_ZH.md、REVIEW.json、resource_planning.png/.pdf已实际生成、看图并本地/远端SHA核验；证明monitoring/H020_CLOSED_INITIAL_RESOURCE_ARTIFACTS_VERIFIED.json SHA64ea9b5f282998e4637c2b578ee14ef711f5b37b7bb064efd9f7584f64f1f9c0（5文件含源码、5文档链接）。来源monitoring/review_h020_closed_initial_resources.py dbd30c1bb49679b3620407bc7fcf22cc29aa4f0cbfdd47d9882c6ab77dd080c0。20µm原闭合段进程树RSS采样峰值41682.796875 MiB=40.7058553467 GiB（包含runner和后代，5秒采样，共享页可能重复计）；之前37–38GiB是常态而非峰值。25µm全部闭合段同口径最大23.3538131714GiB。节点增长1.383倍、RSS增长1.743倍；原20µm含10%余量审查估计35.53008279GiB，被实测超出14.5673%。独立首小时cgroup采样峰42.055625916/42.055664063GiB，OOM增量0；cgroup与RSS口径不可等同。

18.75µm新两档PLC面数/体网格割线外推326135节点、1542728tetra，两种RSS情景44.8255/47.0463GiB；12.5µm52.1783/58.3627GiB。这些不是实测、上界或所有P1细化不可行证明。未生成新网格/启动新CFD/批准P2或dt；25/20空间判定仍待闭合场+同46回放。若仍失败，下一档需要以当前资源证据重新审查，不能沿用旧线性预测直接批准。

新增pending_development/collect_closed_h020_budget_continuity.py SHAe88755f9fdda5f3923c556187c175298c5faa1628b7755059d4952736ea693d5，外层及远端probe AST PASS、已同步、尚未执行。仅在独立RESTART_CONTINUITY_BUDGET_20_TO_30.json和检查段execution/ledger均闭合后运行；支持真实PASS/FAIL结论且保留raw execution状态，核对原47文件、20/30检查点、完整步骤/8rank误差、原1e-8门槛与账目，创建monitoring/H020_BUDGET_CONTINUITY_CLOSED独立快照，归档整个检查段，不碰后续续算。只补对应segment的本地../../SV_MESH链接。当前不要执行缺闭合数据的collector。

新增pending_development/report_h020_budget_continuity.py SHAcebb86396519c262aa4aec91b9ba1fbd774949b1258e4e3ca39c27efd669a498，AST PASS、已同步、尚未执行。等collector成功后在本地scientific venv用--study运行；生成该独立目录REPORT_ZH.md/REPORT.json/continuity_errors.png/.pdf，不覆盖原始证据，分别列u/p和8rank Y/A相对L2、raw执行状态、独立连续性结论，不称稳态/空间精度通过。生成后需看图、文档链接/来源哈希和双端同步核验。原主报告整合草稿仍710cc43...，等新producer退出再部署。


## 2026-10-01 03:24 UTC — h020原生连续性PASS，已自动进入原30步续算

20→30检查段已关闭，原生exit0/health[]，wall_time3125.663186311722秒（52分05.663秒），在3600秒限额内，checkpoint30完整。独立连续性PASS：u相对L2=3.7038325221178695e-15，p=2.9178908968086167e-15，8rank A最大9.96650696968471e-10，所有u/p/Y/A均<=原1e-8门槛。原生检查段raw status仍FAIL（指定30步结束未稳态资格化），没有改写。proof reports/RESTART_CONTINUITY_BUDGET_20_TO_30.json SHAadb7a7a89e79b57000fd750259a230830ea8c2301d9472ecb0a84b0d656d761a；检查段execution SHAe6e787d3ba1c3d626b7936143641593aab07ec3e7e0cc781417463d03e793afa。原30cp cced807e...，重算30cp34960bc48dd08aa9105cb51659810a63009e303ae986d07f9654fc5c62b7a1dd；检查点非位级相同，但通过原状态误差门槛。

followup已自动从原30步进入segments/wall_budget_resume_from30_0001。producer仍776396、job coplanar_h020_budget_followup_0001；新的MPI PID778658，验证MPI776416已退出。CURRENT RUNNING_H020_REVIEWED_BUDGET_CONTINUATION，RESUME_INDEX status RUNNING，checkpoint30/validation_only=false。03:24:15新段elapsed197.32秒，31步FIRST_SOLVE PC构建，CPU7.680核/RAM30.95GiB，health全0。仍无合格h020场、无h020回放，稳态必须在新段再满足5区间。

原段+验证已实际累计25043.80997633934秒，43200秒总额度在续算启动前剩18156.19002366066秒（约5.0434h）。账目reports/H020_WALL_BUDGET_LEDGER.json当前SHA439aed9e4373821113ce6e70a20a2227d64ec49cc99ea3ec0cabf8a02299408e，后续续算关闭才追加其实际耗时。不得把这个额度当剩余耗时保证，也不得预算重置。

已执行pending_development/collect_closed_h020_budget_continuity.py，30个检查段/原始引用/快照文件本地SHA核验PASS：monitoring/H020_BUDGET_CONTINUITY_CLOSED/MANIFEST.json SHA9488ca510b4cfbc5bbc10b2c47ca6e8f1dbfa197dbd0d2ebdb33a0a4aec40ea9，LOCAL_COLLECTION_VERIFIED SHA 4fb48cdd28b92f0e15b004bae5514a668e9496c6571beb1da2375c52aafb642c。完整验证case已在本地segments/wall_budget_check_20_to_30_0001，SV_MESH指本地../../SV_MESH；不要重跑此collector覆盖独立快照。root原预算FAIL不覆盖。

已执行report_h020_budget_continuity.py生成上述目录REPORT_ZH.md/REPORT.json/continuity_errors.png/.pdf，18个相对L2数据点及门槛使用实值，无平滑。图已目视检查，9文件含source/原manifest/快照/报告双端SHA PASS：monitoring/H020_BUDGET_CONTINUITY_ARTIFACTS_VERIFIED.json SHAb06bd8805f3f55f969e3679f837fe7378575663007144c4afb9cd862476ca36a。图只证明重启连续性，不证明稳态或局部精度。

新增pending_development/archive_h020_continuation_saved_step.py SHAb52724659ed7b91bf1ec8051272dd06d2c18fc3514790ecdb2cb01169ea0f1f7，外层/嵌入probe AST PASS，尚未执行。用于续算40/50/60/70/80完整保存并进入下一步后归档，固定选择新continuation目录且要求index/proof/cp身份；校验checkpoint、复制原生VTU/bin，单独保存段内SAVED_STATE_XXX证明和本地回收证明。不得使用旧root collector归档这些续算状态。完整CFD闭合后仍用any_segment collector 3aeefd...，回放后用原replaycollector。

最终报告草稿report_coplanar_h020_validation.py已加入实测资源报告入口和输入哈希检查，当前SHA9e3580cf3252afc1ee9e50f13f1913577a55bc40179e166e27580fe6e3444279；旧710cc43...已保存previous_drafts。H020_REPORT_INTEGRATION_DRAFT.json已更新，AST PASS、尚未部署/正式执行，旧负向lineage guard记录对应710cc43版（核心身份保护未改变）。活动85源和恢复driver/followup全部不动，直到新producer完全结束。原用户短停/总结/自动恢复已完成，开发ACTIVE。


## 2026-10-01 04:04 UTC — h020续算40步归档，稳态区间1/5

当前仍coplanar_h020_budget_followup_0001 / MPI778658，在segments/wall_budget_resume_from30_0001运行；04:03:34已进入41步MAX_AGE第3次PC构建，elapsed2554.94秒、health全0。续算30→40 E_u=1.444051238026335e-7、E_Q=3.6555320171826055e-9，质量误差9.09501627114305e-13，壁面u=0、u_max=.5499245853657136。新段连续区间1/5，不能与原预算段3/5相加，不是合格稳态场。原生31..35均各5NL、36..40各4NL，与原闭合记录一致；之前总结中35步4NL应以原execution实际5NL为准。

已实际运行pending_development/archive_h020_continuation_saved_step.py 40：段内reports/SAVED_STATE_040_VERIFIED.json SHA1863e9f8934c5945b6f47d6385b2cad19a95bea8f0284be4d7700681e1b48ca7；40cp SHA56c9e018c191e74387914aca0c881153cc5c6a5621a1583f8bd3684df385418d；reports/SAVED_STATE_040_LOCAL_COLLECTION.json SHA450d03f37d8ceaabcf5d4679e9bb36c1e45fafb68bdc0e9c0b4609e724b07a2e。两个原生文件本地SHA已通过，server同步本地核验记录。不重跑collector40。下一保存状态50，使用新continuation脚本而非旧root工具。新段本地目录目前仅40步资料，尚未完整同步活动case，也未创建SV_MESH链接；最终any_segment collector会收齐并处理已知链接。

最近内存约38.2GiB，memory.max累计749293未新增、OOM/oom_kill=0，CPU约7.68核。验证报告及资源报告均已完成并双端核验，最终报告草稿仍9e3580cf...，活动85代码未变。当前用户短停/汇总/自动恢复已完成，开发ACTIVE；下一步继续稳态区间、合格后审计+46回放，再闭合后部署报告/壁面查询和25/20空间审查。


## 2026-10-01 04:40 UTC — h020续算50步归档，稳态2/5，更新进度/ETA报告

同一续算segment和MPI778658继续运行。50步已完整保存并进入51步第5个PC构建；04:40:02新段elapsed4745.62秒，CPU7.679核/RAM34.1GiB，health全0。40→50 E_u=3.2568607014569315e-9、E_Q=1.3692861104257136e-11；mass5.293955920339377e-16，wallu0，u_max=.5499245853645819；新段连续稳态2/5，仍未正式资格化或回放。

已运行archive_h020_continuation_saved_step.py 50，原生文件和证明本地/服务器已回收核验。段内reports/SAVED_STATE_050_VERIFIED.json SHA25253c6dfd4e4c6394c37790c370bf3d7de4b52b444a9e99c4298d59123bc5c4；50cp SHA96e75eecf6e4be749bfe098c4165427210aa2e124e5e9f21cbebeec29144745c；SAVED_STATE_050_LOCAL_COLLECTION.json SHA10f4362dd155cb1cb41ea8c332a430cf78737d03709c49104187d41689fdcc8c。新续算40/50均已收集，不重复collector；下一次60。

新增monitoring/H020_CONTINUATION_STEP50_PROGRESS/PROGRESS_ZH.md、STATUS.json、monitor_snapshot.json，04:38:06快照。实查10套合格新CFD/11批分析完成回放，20µm未计完成。该快照沿用先前原30→81条件模型10677秒，减当前实际续算4629.166秒，原生估计剩6047.834秒约1.68h；含后处理、场审计、46回放、比较/回收/报告预计本轮20µm还需2–3h（从04:38算，不每次重置）。前提60/70/80继续通过并正常停止/审计/回放顺利，非收敛或时限保证。全项目若需P2/P1代表区域验证仍预留1–3工作日或更长，未自动启动P2/dt。原生histor.dat的reported_elapsed可能包含检查点继承的旧计时；预算/ETA使用runner真实wall elapsed与closed ledger，不用native打印绝对elapsed冒充本段耗时。

报告源monitoring/report_h020_continuation_step50_progress.py 1b9726f8358280cc1634f79a0b53b3cbc7b5eac7a029b858a73c11ba749ba1d0；已执行、6链接核验、4文件含源双端SHA PASS，证明monitoring/H020_CONTINUATION_STEP50_PROGRESS_VERIFIED.json SHA09e3c90bd59862fa45317a5509a921d7a80a6f7a7a584b1f35d8d1fc73b7b5d0。不覆盖23:18原用户暂停汇总快照。当前开发ACTIVE；85活动源及恢复三程序原样，报告整合草稿9e3580cf...仍未部署。最近memory.max仍749293、OOM0，无新增事件。


## 2026-10-01 05:08 UTC — h020续算60步归档，稳态3/5

同一segments/wall_budget_resume_from30_0001、MPI778658持续运行；05:07:24进入61步第7次MAX_AGE PC构建，elapsed6388.53秒、health全0。60步已完整保存并独立回收SHA核验，50→60 E_u=6.92073671421056e-11、E_Q=5.43953970814871e-14，mass1.764651973446459e-16，wallu0，u_max=.549924585364235。新段连续稳态3/5，不能合并原预算段的3/5，仍未资格化/回放。

archive_h020_continuation_saved_step.py 60已实际完成，不重复运行：段内SAVED_STATE_060_VERIFIED.json SHA0f32984c64fe337f445462df823d4f10a47b356fd63fef51b5a2400db0929a29；原生60cp SHA44b3e13841779c6dfa54ae46646cb63e56ac58e71da423e8d2ca0447e3325c68；SAVED_STATE_060_LOCAL_COLLECTION.json SHA07c3c966daa9e93dc9e662696b775b592b3f619be159c38dd4bd2c9f1920c386。新续算40/50/60均已收集，下一70；旧root10..60文件独立保留。

最近普通迭代内存约38.5GiB，CPU约7.68核，memory.max仍749293，OOM/oom_kill为0；无新增数值失败。04:38进度/ETA快照保持历史，不每次重置其2–3h估计。全部正式微泡1565及ACTIVE_FLOW未更改，h020仍是独立实验，主报告整合待producer完全结束。开发ACTIVE，无新暂停请求。


## 2026-10-01 05:39 UTC — h020续算70步归档，稳态4/5

同一续算segment、MPI778658持续运行；05:38:51进入71步第9次MAX_AGE PC构建，elapsed8273.42秒、health全0。新段60→70 E_u=1.3919122895068985e-12、E_Q=1.0587911840678754e-15，mass3.529303946892918e-16，wallu0，u_max=.5499245853642302。连续稳态4/5，仍须80步第五区间以及原生正常退出，当前不是正式合格场/没有h020回放。

archive_h020_continuation_saved_step.py 70已实际执行并本地SHA核验：段内SAVED_STATE_070_VERIFIED.json SHAc5daa2f577fdb6ccfd9f18667e69c88bfa30de551acdeb55b789b4576dbb31a9；70cp SHA9be7af45bdca09053506392e63ea92de15cfb5b713b65bf66d342a8fa419f2f7；SAVED_STATE_070_LOCAL_COLLECTION.json SHA9625a0fd9d40a2422c8378b9cfdee01304208c074a6990f25f0515db51bc604a。新续算40/50/60/70已收集，不重复。下一80保存且进入81后收集80；若80触发稳态STOP_SIM仍要等待原生当前步完成，不能先假定退出。

目前CPU约7.68核，普通阶段内存约38.7GiB，memory.max仍749293、OOM0，无新增数值失败。活动源85及恢复程序未改，报告整合草稿9e3580cf...仍待整个producer退出后部署。后续流水线将独立身份核验/场审计/同46条回放/分析/原报告；闭合后再执行实际路径wall query与新报告整合。用户短停/汇总/自动恢复已完成，开发ACTIVE，无待暂停指令。

### 2026-10-01 06:09 UTC：20 µm续算已完成五次稳态检查，等待原生正常退出

实际续算第80步已保存；40/50/60/70/80连续5/5通过，程序于06:08:57 UTC发出five consecutive steady intervals and physical gates正常停止请求。原生求解器已进入第81步预条件器重建，仍须完成该步才能正常退出；此时尚无最终合格场、场审计或20 µm轨迹结果，不能宣称整轮完成。唯一活动生产器仍为coplanar_h020_budget_followup_0001，未启动dt/P2或正式群体批次。

第80步原生检查点/VTU已独立核对并收回本地：E_u=3.0838435797285595e-14，E_Q=2.205814966808074e-16，质量误差3.529303946892918e-16，无滑移壁点最大速度0；O1/O2/O3分流0.3330815453214332/0.3333717340193113/0.3335467206592554。仅表示保存状态和全局稳态检查，不认证局部场精度。

位于03_cfd_runs/M1_coplanar_fine_h020/segments/wall_budget_resume_from30_0001/：检查点080 SHA74e177a3e2f464de57960bbb0626b5125bb1779d8874474948567934bffed770；reports/SAVED_STATE_080_VERIFIED.json SHA df9bdfbf620f874d592ca73e6f013dff56769685fe14d9861c3a2e34e8b1696d；reports/SAVED_STATE_080_LOCAL_COLLECTION.json SHA0e3baded88518cbc41ecdf3b592ba591b8991d86a55db7449ca4b39a98a2905a。收集脚本正常exit0，原数据不变。

用户要求的短暂停止/进展报告/自动恢复早已完成，当前继续开发。待原生PASS退出后，生产器自动执行后处理、canonical身份检查、场审计、冻结46初态回放及比较；仅当整个生产器EXITED后才能部署pending_development中的新报告版本或运行额外重任务。最新报告草稿9e3580cf...已同步但未部署，85项运行源码仍冻结。

### 2026-10-01 06:39 UTC：20 µmCFD与场审计已完成；场门槛FAIL，46条回放进行中

20 µm原生续算最终81步、exit0、health_failures=[]、status PASS，wall_time_s=10816.527707576752；execution SHA47619e36478eec152efe1273976d8f7435511855597f35455af302b02795f465。全部原生段累计35860.33768391609s，12h总预算剩7339.662316083908s；ledger SHAfa458c8b1ef54245e51f785f98fe63822cfe99c945440b955666f616a2702053。原6h失败段仍保留，未改写。最终续算段仍segments/wall_budget_resume_from30_0001，RESUME_INDEX已QUALIFIED_STUDY_FLOW（SHA26f7d7deadfe9af26f6576940554e57c542fabd5f0da54e12014ba39f9891566）。

最终flow SHA cb36025fece24f087f550f4e8753c8c868a1894f67ec444fb498af6b42220925；native SHA be4b2bc24f13140d91863b8fd9d30eaceea45c7d0a79cf9313a69ceedb0870c2；mesh SHA e6fc8966bf25d8d2b2190265ed98dd53cef90c055115bad24b8e317a564f56e3。STUDY_FLOW报告SHA15fad6ae17622ed54bb1083ea7c7f42f25456f345ea3b8738a2f12de6e345ce4，独立canonical身份PASS，证明SHA9eb672b99f6b360918628a67c15aee72e9279e24c1bcc15d6c841f5995c99477。最终质量误差1.764651973446459e-16，最大壁点u0，WSS均值3.0751888565863297Pa、P955.717266737924176Pa、最大27.55816780237174Pa。正式ACTIVE_FLOW未替换。

场审计634.8188254833221s完成，1059固定点全部在域内，zero_tetra=[]、eta_zero_gradient_count0、centroid_zero_speed_count0。CPU/GPU G最大差1.5916157281026244e-12仅为计算一致性。全场max|div|4584.603041364238s^-1。SUMMARY SHA b16a7a22a4aae72961e0c4ae6c798c6218f7689d74c09c6ce73fecf40e5903d2。已用any_segment收集器完成141文件逐一SHA验证：CLOSED_H020_CFD_OUTPUT_MANIFEST.json SHA dcb44eb10747ea6fa93243c3ce77906ec2584c516f716fbecc92121103fc3cb9；CLOSED_H020_CFD_LOCAL_COLLECTION_VERIFIED.json SHA318e298012723275b5f9a8fb54e0e3e99e3be22fa09c98b42f6f8e28fcab152e。无需重跑该排他收集器。

新增独立场比较monitoring/H020_CLOSED_FIXED_POINT_COMPARISON/REVIEW_ZH.md、COMPARISON.json、field_convergence.png/pdf，已执行、查看并双端验证5文件/3链接。源码pending_development/report_h020_closed_fixed_points.py SHA7152095053674db90a1aedee37661fb42a027c4ccaca0cd767584f6a19c37b1b；比较JSON SHA5e117d6cd5b09857b762cd7f0478f8ace4bbb29a66c2357d8cdb2c6a510764b3；双端证明monitoring/H020_CLOSED_FIXED_POINT_ARTIFACTS_VERIFIED.json SHA39241997be3fea0ae7374254830d60708e5b22108c47a07c9671b9404e50713e。此为回放关闭前的场比较快照，不是完整轨迹判定。25/20 ALL：u1.680135%、G16.588714%、WSS22.112642%；六区域u均<2%，G/WSS均超过各5%门槛，场FAIL。high_eta区G25.796676%、WSS37.520462%；zero区G20.226279%、WSS15.439108%。不得选production或直接做正式dt/P2。

原388769中心20 µm：u0.078366828544123m/s，div-64.80922289030491s^-1，G1986.5017187410656s^-1，eta0.032624800813853626；原138584中心u0.02347642529326879m/s，div6.089754709468679s^-1，G914.156253017135s^-1，eta0.006661612486234861。单点改善不代表区域/轨迹收敛。

06:39:01监控：唯一生产器coplanar_h020_budget_followup_0001仍RUNNING，CURRENT_EXECUTION_STATE=REPLAYING_REVIEWED_COPLANAR_COHORT；20 µm冻结46初态回放elapsed376.93s、5526迭代、24/46终止，其余继续，未有整批终止汇总或分析，不得把24条计成出流。GPU瞬时36%，RAM7.56GiB。旧85项源码仍冻结ca522c...。

新增待部署辅助程序pending_development/deploy_closed_h020_reports.py SHA00fcda40b4516de44fc1ac535c1974a869aae188c6e449a1f057f6437f10516b，语法通过并已同步。只允许整个生产器EXITED、所有活跃作业为空、46分析完成、壁面查询完成后部署原三个报告草稿，并先备份旧报告/figures/源码快照到monitoring/PRE_H020_REPORT_INTEGRATION（尚不存在）。实际负向验证H020_REPORT_DEPLOYMENT_ACTIVE_GUARD.json：活动原生持锁时拒绝预检查，85源码/SNAPSHOT不变、无发布输出，不能把它说成部署成功。正向部署尚未执行。

接下来：等待46回放/分析/自动比较/旧报告生产器整体EXITED；收集完整replay；再独占运行原版本audit_coplanar_wall_queries.py（25/20、reference25、输出COPLANAR_FINE_H020_PATH_WALL_QUERIES），然后运行报告部署辅助程序及新版render_report，检查图片/链接/双端哈希。完成25/20完整评估和资源可行性复核后，按原Phase2G决定进一步P1细化/有条件代表域P2，不跳过前置门槛。当前未再暂停，用户要求的短暂停止/报告/自动恢复已经履行。

### 2026-10-01 07:03 UTC：20 µm验证与报告全部闭合；18.75 µm仅网格检查已启动

20 µm冻结46初态回放734.7184898853302s完成，44出流、2观察截止、0数值失败、0有支撑静止；O1/O2/O3=16/15/13，截止ID393、1377，冻结slow-fast筛查24事件（25 µm为20）。完整25/20比较FAIL：0/46全部指标通过；共同出流但出口变化12条[2,4,9,168,440,487,650,1041,1077,1189,1208,1505]；终止类别变化3条[5,534,1329]。速度史46条均未过5%门槛，端点与停留时间各40、最小间隙23、路径长度12未过。历史失败ID440现在O3、residence1.353973870579249s；1505现在O2、1.7897540451226592s，不能把出流当成收敛。旧1565身份不变，也没有新正式群体批次。

376回放文件已收回并逐SHA验证，不能再运行排他收集器：CLOSED_H020_REPLAY_OUTPUT_MANIFEST.json SHAa14672d44b5e96e068ba659ece21d09013f240d7e9d5355db3b7f021bb6e75a7；IDENTITY SHA6c379e3cbc4f98349f15802170e9d0083749eae5fd5ca48e85613a84c99e76a7；analysis/SUMMARY SHA0a17a96bf0e6035565a06fe1764923d244b97e6c0efca20e60ddb86c613e8f18。生产器coplanar_h020_budget_followup_0001于06:47已EXITED，不得重启。

追加GPU只读壁面查询已完成：01_mesh_capability/COPLANAR_FINE_H020_PATH_WALL_QUERIES/REVIEW.json SHA6e819fe82d92ef8f311afe32368e514b61f2608c3efcd0fb90b35742af4cd903。在25 µm全部407106记录位置（242336近壁查询）上，25/20最大间隙差8.51098705401121e-18m、法向向量差0、接触数量变化0、feature/rim变化0。耗时12.764s，无重新积分；只覆盖这些坐标，不证明完整接触求解或所有20 µm路径等价。

已部署三个报告源码，旧85快照与报告/figures保存在monitoring/PRE_H020_REPORT_INTEGRATION；现在CURRENT_STUDY_SCRIPT_SNAPSHOT为86项、SHAdb2266fe8a4033612fe15dd317e41783afe33a81b09857b91df681b1746bfae2。render_report SHA799f09b3f3cd7e175defe898ebe70555c8bd9da6818c5a65f9f2ef2ec0331a91（增加独立场图链接），report_coplanar_h020_validation SHA9e3580cf3252afc1ee9e50f13f1913577a55bc40179e166e27580fe6e3444279，report_fine_restart_recovery SHA44f922737f688296d6f2783ce663ee13ad99c1b5d7a6e1da08790c9e2425f3d6。部署助手已正向执行，勿重跑排他归档。

H020_INTEGRATED_REPORT_EXECUTION.json记录新版报告实际PASS、exit0、58.674601793289185s，ACTIVE_FLOW SHAde168191873967374caea00083ea6c34baec164dd873795a5fef1f6e3c19cb6d不变。报告时刻06:49:55.968724UTC，REPORT_MANIFEST SHA6c755ccf0f3ab033b0d845f02f3a2f3467bd92ff72a279b29da3ead8f92c20f2。pending_development/collect_integrated_h020_report.py SHA4302cacbe24b3d76f3bc9e6b9221cd750e43d91dca0d9a9212d474df0343b8db已实际完成：86源码、106图文件、HTML+h020MD共201链接通过。证明H020_INTEGRATED_REPORT_LOCAL_VERIFIED.json及更新的REPORT_LINKS_AND_HASHES_VERIFIED.json；旧链接证明另存monitoring/H020_REPORT_LINK_PROOF_PREVIOUS.json。三张新h020图已目视核对，显示原始接受状态、无平滑/事件对齐。OPEN_RESULTS.html和COPLANAR_H020_VALIDATION_REVIEW_ZH.md已发给用户。该静态报告的active=[]属于其生成时刻，当前新网格任务见下。

遵照Phase2G，在确认完整25/20 FAIL后继续下一档有界网格可行性检查，尚未开发P2或正式dt。新计划COPLANAR_H01875_MESH_PLAN_ZH.md SHA c307e529580b0294a436eda18e2fd9e831365964d8c7a80023279888272546d9；审查COPLANAR_H01875_MESH_REVIEW.json SHAe4055ebee1ca4a1d0e61b29b2425c4fec0e6f0feedfd07db3dc984af11e5792c；pending_development/h01875_mesh_only_followup.py SHA2710c49491d3ac179e60fc18a114bba1c44c4c6081b2f14be272892502638558。首次未执行草稿5f309...及review ab726...保留于pending_development/previous_drafts/h01875_before_state_guard，当前版本增加启动状态及步骤间PAUSE检查。

新job coplanar_h01875_mesh_only_0001于07:02:41UTC启动，唯一活动job。预检查证明SHA72e99ccaef6e3d92ca6b226e5b95b0bf2e0301b3ee7f1184b32d314d7deb5310、job SHA1b721e403dac7863f2763843d33c20496bca5b251b0e73ba59fe7936762b01e4。CURRENT_EXECUTION_STATE=RUNNING_H01875_MESH_ONLY_RESOURCE_REVIEW。目标probe M1_coplanar_h01875_probe、candidate M1_coplanar_fine_h01875，comparison_medium=h020。表面上限300000→330000、tet上限1900000、native网格900s、prelaunch RAM<30GiB，科学/物理门槛不变。预计mesh+QA约2–20min，非保证时限。CFD_automatic_launch/particle_automatic_launch均false。当前86项源码与新pending驱动冻结，任务结束前勿修改。

07:03观测：输入表面实际154739点/309474面，与已测精确面数一致；原顶点/端口三角形不变，原面面积误差6.83e-15、平面误差2.71e-20m。原生网格已生成337592点/1613894 tetra，原始QA/随后共形拓扑处理仍未闭合，不能算合格候选。旧资源外推18.75为44.83/47.05GiB，超过37GiB CFD规划与42.06GiB容器额度；即使GATE_MESH PASS也不授权CFD，尤其不能用旧线性RSS QA结果绕过最新资源审查。下一步待实际网格QA后重新评估内存，可独立研究更多MPI分区，但目前没有改变MPI/求解器配置或声称内存改善。

原生base runner已核对可从policy.MPI_ranks读取进程数，默认启动command为wrapper+launcher -n ranks+wrapper+exe；而study prepare_candidate_cfd.py的身份metadata固定CPU_MPI_ranks=8，coplanar_output_runner.py候选白名单也仅含50/25/20。未来更多MPI/新候选需要独立准备与身份检查，不能只改数字后直接投入长CFD。目前未启动此类实验。用户没有新的暂停请求，最早短暂停止/报告/自动恢复已完成，开发持续。


### 2026-10-01 07:53 UTC：18.75 µm 网格 QA 闭合；16 MPI 资源短试验已启动，开发 ACTIVE

18.75 µm mesh-only job coplanar_h01875_mesh_only_0001 已于07:40:48.961752UTC结束，H01875_MESH_ONLY_CLOSED=CLOSED_MESH_QA_PASS，实际2288.6561284065247s（38.1min），原2–20min规划偏短。原生337592节点/1613894tet有453 all-WALL且float32墙面精度FAIL，原记录保留；恢复声明PLC输入精度并完成405内部边扇共形细分后，337997节点/1615458tet、all-WALL0、负/零体积及非流形0、1063固定几何探针全部在域内、体积差0、五角色细分PLC三角形坐标集合精确一致、四端口原三角形不变。严格域QA PASS。minSICN最低.020611120735138218，中位.8134400522315466，低质量607。不能将mesh PASS当CFD/空间收敛通过。

72文件（含原始失败、共形候选、输入PLC/ancestry、MPI启动能力记录）已本地/服务器逐SHA回收。collector pending_development/collect_closed_h01875_mesh.py SHA368a21231bfe76ec9e7c29fe3b18113a8792d1ade126c4556203ca164c5f7ebb；CLOSED_H01875_MESH_OUTPUT_MANIFEST SHAa4ac9a4573ebcfb164607299ba254e9ce6fb9ce2e487c1b9ee7d74533b7a5804；CLOSED_H01875_MESH_LOCAL_COLLECTION_VERIFIED SHA406f7b7d1866297b7c75b9e50bf062da31c1341593fb9f9341269879b7db8a63；GATE_MESH SHAd8f1e40cae163f7520daeae77eac81f1ffbf8d07173614a3ca125c7b69a59ac9。收集器已完成，勿重跑排他proof。

实际节点数下8rank外推RSS46.45584956775984/49.55544275313739GiB，高于37GiB规划及42.0556640625GiB容器限额；是规划情景，不是实测或上界。新增01_mesh_capability/COPLANAR_H01875_MESH_RESOURCE_REVIEW/REVIEW_ZH.md、REVIEW.json、mesh_and_resources.png/pdf，图已目视，8文件9链接双端核验，证明H01875_MESH_RESOURCE_REPORT_VERIFIED SHAe1bcf279f80650c2face1cc6acdf154b5c7727553e337a89635b6f469c393c46。报告JSON SHAc9537f94e76b2d6c2c4b80eb53ed0cc65f5664b2dbcb692b1103667119940b70。

原生启动器8rank默认PASS，16默认slot不足FAIL已保留；显式--use-hwthread-cpus可启动16个进程（每个仍affinity0..15，CPU总额度7.68核），只是启动能力，原记录位于01_mesh_capability/MPI_LAUNCH_CAPABILITY_8_16。没有修改全局MPI配置。

为实测更多分区资源和数值一致性，新独立case位于01_mesh_capability/MPI_PARTITION_RESOURCE_PROBE_H020_16/case。使用现有20µm mesh（mesh_arrays f93bff504ef4bd035b58157f9e39e5d73205c082a0cf9f3c934670fcda9c3f2d），零初态到第10步；MPI8→16，启动命令显式硬件线程槽位，总步500→10，其余XML/物性/BC/dt.01/PETSc及非线性容差不变。预算9000s、正常停止宽限900s、采样cgroup39.5GiB软停/41GiB硬停，仅终止本试验进程组。原先未执行的3600s草稿及初次guard保存在pending_development/previous_drafts/mpi16_before_startup_cost_review；修改依据原8rank第10步reported elapsed6660s/137次线性求解，不预设16rank加速。

审查H020_MPI16_RESOURCE_PROBE_REVIEW SHA0c21fb9ec6357cf851fca9ccc2a691a7b8e09169668e5710a87c6e02096acb45；计划H01875_RESOURCE_AND_MPI16_PLAN_ZH.md。驱动pending_development/probe_h020_mpi16_resources.py SHA422e3a9047aff913272065ac9ac75139d2290f20ca33962939612c84bc7f8efb；语法/5处adapter完全逆向还原base源检查通过；实际并发锁负向guard V2 PASS（此时mesh已结束，由test_parent持同一锁，不称仍由mesh持锁）。H020_MPI16_RESOURCE_PREFLIGHT_VERIFIED SHA f81ae286737abaa0ed1eaba81dd6388d31614a77aec08489c5f19780c1ba5fbc 已正向PASS。执行adapter SHA62c844eecf9a5d2b0eb7ebff29c6e28f7563d69073d011255d13c7f0a1678053，11输入/源码本地proof SHA760bd892744e539617239b2789679ed2c15e30733ee345fc3cef59b825028267（输入证明不等于CFD完成）。

job h020_mpi16_resource_probe_0001 于07:48:47.982028UTC启动，PID788441，原生MPI进程组788449。仅此昂贵job活动。CURRENT_EXECUTION_STATE清理旧暂停/续算字段，状态RUNNING_H020_MPI16_RESOURCE_PROBE。运行期间冻结86活动脚本快照db2266fe...及新pending驱动422e3a.../review0c21...，不得改活动源/配置。监控使用monitoring/monitor_mpi16_resource_probe.py，写MPI16_LATEST.json与MPI16_HISTORY.jsonl；旧/tmp/brava_monitor_live.py虽active正确，其CFD摘要仍硬编码h020已完成PASS，不能拿它判断新短试验完成。

原始固定10步execution预计仍FAIL（不满足稳态停止），必须保留。独立RESULT只有在原生正常step10/残差health/无停止或OOM、同canonical u/p/G相对L2与Linf均≤1e-8、采样cgroup<39.5GiB时才PASS_RESOURCE_AND_PARTITION_PROBE。G L2体积权重，u/p同节点原值，不平滑/插值。8rank参考原始step10 VTK SHA1f1acf461a968c7011e500b90b478c1def255c9e4fdc5582d19a347b411c731e；cp10 e1e56de0118089ab487ab5dfdb70edde22f9832ce18f78c404e68c5046c2edf8。试验PASS也不自动授权18.75µm全CFD；仍须资源审查，不能声称稳态/物理/空间精度通过。

07:50首次监控elapsed87.6s，首次PC分解，cgroup32.36GiB/treeRSS26.26GiB，GPU6420MiB，OOM0、无STOP、无保存步；仅早期样本，不能当峰值。预计本短试验加核验2–3h（从07:48启动计，不每次重置）；若后续可行，18.75µm完整CFD/场审计/46回放另6–12h。整个任务仍依赖空间gate，不作全部完成保证；正式dt/P2/P1及生产群体尚未开始。此前用户短暂停止/进度说明/自动恢复已履行，没有待处理暂停请求。


### 2026-10-01 08:12 UTC：用户新指令暂停，等待明确恢复；取消先前自动恢复语义

最新用户要求“下一个可停止点暂时停止，输出已完成/结果/后续/预计时间；完成后我会给指令是否恢复开发”。已执行，**现在 PAUSED，禁止自动恢复或继续开发计算**。上文旧的 ACTIVE/自动恢复描述已被本条覆盖。

16 MPI 首次资源试验于08:03:30UTC结束，881.8395447731018s，cgroup采样41.0255012512207GiB触发预设41GiB硬停止；原生exit1、final_step=null、没有完成第1步/没有可比较输出。无OOM/max事件增量；treeRSS采样35.01164245605469GiB，仅是中止前样本，不是完整峰值上界。RESULT=FAIL_RESOURCE_OR_PARTITION_PROBE；equivalence_metrics=null，不能声称16MPI不一致或一致。execution SHA09e95578b643b20204b58193634c0f444645cc9f8fe8d37c92fc449d3c6631bf，RESULT SHA f4fc74fda621f86b9059abc1b9757cc7a125179e4d8adf108bbf6f559710b4b1。29文件已双端SHA回收：CLOSED_MPI16_RESOURCE_PROBE_MANIFEST SHA3d6bee7b8c76745ed285d32d97a061c9c1d95f60caf8a856ebe49d2353ba4343；CLOSED_MPI16_RESOURCE_PROBE_LOCAL_VERIFIED SHAc7235ef7c661e3553c0f03d20ae215d1f124b056e17351f2c0bcfa31c8aa5ff4。收藏脚本pending_development/collect_closed_mpi16_resource_probe.py SHA3fca03580624dbd98a1d7e89fc6b76a32bd23723915530b3b99697d7610ee0a8已完成，勿重跑排他proof。

期间完成只读缓存初查：pending_development/audit_closed_evidence_page_cache.py，输出01_mesh_capability/MPI16_PAGE_CACHE_AUDIT/RESIDENCY.json。583个已归档证据文件有1.601959228515625GiB驻留页；mincore不读取大文件内容、不请求回收缓存。停前MEMORY_STAT_AFTER_SOFT_STOP记录anon约30.50GiB/file约8.30GiB；这些类别可能交叠，不能任意求和或把1.60GiB直接从cgroup中减掉。文件缓存调查尚不完整，**没有执行posix_fadvise/清缓存，没有资源重试或18.75µm CFD**。只读查询MPI安装得默认mpi_yield_when_idle=false，记录pending_development/MPI_IDLE_PROGRESS_CAPABILITY.json（服务器）；未改变当前MPI配置，未验证该调度参数的数值或性能影响。

新增固定点分组补充：monitoring/COPLANAR_FIXED_POINT_STRATIFICATION/REVIEW_ZH.md、REVIEW.json、fixed_point_strata.png/pdf；源pending_development/report_fixed_sample_strata.py。5文件5链接双端哈希核验、图目视通过，证明monitoring/COPLANAR_FIXED_POINT_STRATIFICATION_VERIFIED SHA faa18074d06d02cf1ed717f161e8a9f5b521901ad54a971d94fed1bd0d4da82f；JSON SHAd0b8961512e254355124f30b70926938e5daadd3d6b8d9bf74e5ad882168681d。25/20最靠壁两个130点组u相对L2均14.223548%（ALL1.680135%）；20场中位壁距.135981/1.359812µm，中位u.100970/1.009701mm/s，最大速度向量差.113110/1.131103mm/s。50/25这两组u16.333698%。各网格130对最近射线点全部同单元；ray fraction是原壁面面心→对顶点比例，不一定法向。球域并集内345点u/G/WSS1.333372/21.975626/27.715555%；并集外714点2.259107/10.974946/6.395081%。这是事后分层诊断，不改变冻结gate、不把采样误差份额解释为因果或全域概率。

实际暂停08:12:32.539850UTC，CURRENT_EXECUTION_STATE=PAUSED_BY_USER_AFTER_MPI16_RESOURCE_PROBE、active_jobs=[]、resume_requires_explicit_user_instruction=true。全部study supervisor无活动任务，重计算锁可取得，旧MPI788449不再存在。暂停快照pause_snapshots/20261001T081232Z/，PAUSE_VERIFIED SHA7666cec50f106c017ad41012169be365158fe818f85d63cb09527b35557dc61a；保存PRE_PAUSE_EXECUTION_STATE、PROGRESS_ZH.md、LOCAL_PAUSE_REPORT_VERIFIED，并追加最终已知原生CFD可执行文件/proc/exe无残留进程和报告远端SHA核验。实例/服务保留。正式ACTIVE_FLOW SHAde168191...、86源码快照db2266fe...、旧1565及动画均未更改。

完整已完成/结果/后续/ETA见pause_snapshots/20261001T081232Z/PROGRESS_ZH.md，包含现有综合图文链接。本轮用户要求已履行后应结束，**等待下一条明确恢复指令**。恢复后先读这份暂停报告及16MPI失败/缓存证据，再审查资源方案；不能直接重启旧失败job或自动启动18.75 CFD。条件规划：资源复核/有界验证2–4h；若可行下档CFD+场审计+46回放6–12h；空间通过后dt约1–2h。顺利情况下剩余约10–18h；若需局部P2/P1另1–3工作日或更久，不保证收敛或绝对完成时间。


<!-- LOCAL_P2_CURRENT_STATUS_START -->
## 2026-10-01 局部 P2 / 近壁审查已按 A 区门槛结束（2026-10-01T13:35:42.769575+00:00）

用户的新指令已替代旧全域 P1 路线：不恢复 18.75 µm CFD、17.5/15 µm 网格、旧微泡批次或旧资源探测。新目录 `computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/`；服务器对应 `/workspace/brava_flow_roi_18mlmin_20260928/p2p1_nearwall_cross_validation/`。新边界使用原“中心线阻力网络 → 出口压力”，未使用校准 CFD；局部入口按网络流量、下游按沿程网络压力。

结论 **LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED**。A 的局部 P1、P2-medium、P2-fine 已完成；P2 两档为 28,268/101,401 tet、130,794/447,679 DOF。内部速度差 1.30%，但 0.1/1/10 µm 近壁速度差 21.76/21.62/20.22%，G 差 21.75/21.47/19.31%，WSS 差均 15.95%，未通过冻结门槛。fine 正常退出，最终非线性残差 1.87e-10、全局流量缺口 9.36e-10；全局守恒不能代替局部精度。异常单元中心细分后有 6 个单侧梯度，η 范围 0.092–0.409，不能择一作改善证明。

近壁五类夹具及 ID1189/637 共 33,198 个保存状态重算完成：**IMPLEMENTATION_MATCHES_DECLARED_MODEL**，最大速度复现误差 3.06e-18 m/s；没有新积分、近壁源码或参数修改。硬接触启用的速度非光滑属于声明模型局限。按 A 区停止规则，B/C 只提取未求解，7 条新轨迹未执行；本轮科学作业已退出，没有自动下一阶段。

[中文报告](computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/FINAL_P2_NEARWALL_REVIEW_ZH.md) · [六张图文入口](computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/OPEN_RESULTS.html) · [机器状态](computation_examples/brava_flow_roi_18mlmin/p2p1_nearwall_cross_validation/CURRENT_EXECUTION_STATE.json)。原正式 ACTIVE_FLOW、旧 1565 条（1422 出流/141 截止/2 失败）和动画身份保持；不能改称新固定灌注协议结果。旧任务暂停记录只属于被替代路线，不是恢复授权。
<!-- LOCAL_P2_CURRENT_STATUS_END -->

<!-- LOCAL_FORMULATION_BOUNDARY_STATUS_START -->
## 2026-10-01 最新任务：Domain A 边界 benchmark 已技术性停止

本段是最新执行状态；旧全域加密、旧资源探测和上一轮局部 P2 任务均不恢复。本轮请求存档为 `computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/provenance/USER_REQUEST.txt`。新目录 `local_formulation_benchmark/`，服务器对应 `/workspace/brava_flow_roi_18mlmin_20260928/local_formulation_benchmark/`。

状态 **TECHNICAL_STOP_BC_B_SOLVER_FAILED**；边界门槛 **NOT_EVALUATED_BC_B_SOLVER_FAILED**，科学主结论 `main_case=null`。仅新跑一次 Domain A 的 BC-B P2-fine，复用原 fine 网格、P2 空间、求解器参数、core 和全部固定探针，BC-A 未重跑。BC-B 按本轮最新提示使用正式全域三维截面数据：入口三分量速度、出口完整应力牵引；BC-A 仍为中心线阻力边界。正式 CFD 未重新校准或修改。

BC-B 在首档延拓 0.25、第 40 次 SNES 迭代后未收敛退出，残差 0.04123218，目标约 6.97e-8。耗时 15.645 分钟；8 MPI，实测使用约 7.62/7.68 CPU 核配额；进程树峰值 10.102 GiB、绝对 cgroup 峰值 16.786 GiB。无收敛 FIELDS，不能计算有效的 BC-A/B core 差异，不能宣称 BOUNDARY_DOMINATED，也不能用此失败证明 DC 优劣。未调参、重试或增加迭代上限。

边界提取/施加 QA 通过；原始 P1 应力内侧单侧导数仍存在离散相容性和物理精度限制。BC-B 入口原始截面流量 17.524174 mL/min 未强制归一化到 18。PETSc 实测不支持 CUDA/HIP/Kokkos Vec/Mat；因无合格 BC-B 字段，GPU 字段比较未运行。不要把预备脚本或历史 GPU 结果写成本轮新计算。

服务器监督程序 `brava_boundary_benchmark` 已 EXITED、autostart=false、autorestart=false，MPI 主进程 796084 已消失。未启动 DC、Nektar++、新轨迹或任何 B/C 计算。已有 `LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED` 与 `IMPLEMENTATION_MATCHES_DECLARED_MODEL` 冻结不变；141 个原 particle 源文件身份检查通过，未重复 near-wall 验证。

[最终判断](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/FINAL_DECISION_ZH.md) · [中文边界报告](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/FINAL_BOUNDARY_REVIEW_ZH.md) · [两张图文入口](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/OPEN_RESULTS.html) · [停止依据](computation_examples/brava_flow_roi_18mlmin/local_formulation_benchmark/WHY_CONTINUE_OR_STOP_ZH.md)。建议下一阶段先做有界的人工截面牵引传递/非线性收敛相容性诊断，仅为建议，未经新指令不得自动执行或恢复旧路线。
<!-- LOCAL_FORMULATION_BOUNDARY_STATUS_END -->


<!-- BOUNDARY_TRANSFER_DIAGNOSIS_STATUS_START -->
## 2026-10-01 最新：Domain A BC-B 有界诊断完成，全部停止

完成时间 2026-10-01T20:23:13.213210+00:00。本段覆盖上文旧任务的 ACTIVE/恢复表述；未经新指令不得自动继续。最新请求存档于 `computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/provenance/USER_REQUEST.txt`。本轮目录 `boundary_transfer_diagnosis/`，服务器 `/workspace/brava_flow_roi_18mlmin_20260928/boundary_transfer_diagnosis/`。

主诊断 **TECHNICALLY_INCONCLUSIVE_AFTER_BOUNDED_DIAGNOSIS**；结果状态 **BOUNDARY_SENSITIVITY_NOT_EVALUABLE**。旧 `local_formulation_benchmark/` 失败案例 48 项冻结文件原样保留。BC-A 未重跑；网格、447679 DOF、8 MPI/1 线程、物理、NS 容差及 0.25/0.5/1 延拓均保持。

源 P1 三截面流量：入口 17.524174 mL/min、出口合计 17.547795 mL/min，缺口 0.134792%，与原 P1 积分散度一致。入口原方法与严格 P2 trace L2 投影相对差 4.06e-16；出口原载荷与 quadrature 直接积分相对差约 2.38e-15。两侧总牵引跳跃在固定压力基准下为 1.02%/0.40%，但切向分量为 24.26%/33.02%；保留不确定性，未证明它是失败主因。固定 BC-A 状态的出口替换代数扰动大于入口替换，不能将残差当物理误差。

全域投影初值残差改善 6.51 倍，未达到预设十倍强相关判据。Stokes 48.25 秒收敛、残差 2.34e-12，仅 DIAGNOSTIC_ONLY。ns_full_projection：NOT_CONVERGED，最后档位 0.5，KSP -3/3000，9.27 分钟；ns_stokes：NOT_CONVERGED，最后档位 0.5，KSP -3/3000，9.43 分钟。旧失败首因分类 NONLINEAR_CONVERGENCE_STALL；本轮详细证据和线性退出原因见 NS_DIAGNOSIS.json，不把失败字段当有效 Boundary Gate 比较。

原正式 CFD、旧 1565 条（1422 出流、141 截止、2 数值失败）、动画及固定灌注协议的结果身份不变。未重跑 near-wall/parity/网格收敛/粒子；`LOCAL_P2_NOT_YET_SPATIALLY_RESOLVED` 与 `IMPLEMENTATION_MATCHES_DECLARED_MODEL` 仍按历史范围保留。没有启动 DC、Nektar++、B/C、extra-fine 或全域加密。Supervisor `brava_boundary_diagnosis` 已退出，autostart=false、autorestart=false；双端交付校验文件位于本轮目录。

[中文诊断报告](computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/FINAL_BOUNDARY_DIAGNOSIS_ZH.md) · [四张图文入口](computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/OPEN_RESULTS.html) · [当前状态](computation_examples/brava_flow_roi_18mlmin/boundary_transfer_diagnosis/CURRENT_EXECUTION_STATE.json)。下一步仅建议保持边界/网格固定，针对失败档位的线性化系统与预条件器做最小有界审查；不自动执行。
<!-- BOUNDARY_TRANSFER_DIAGNOSIS_STATUS_END -->
