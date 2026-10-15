# 当前clot原生软化开发入口

本轮：lps_native_softening_dynamic_commit_and_small_crack_v1 / 20261013T214800Z。

- **结论：P3第244步(122µs)因最大Δd=0.102554>0.10停止；P4 NOT_RUN。**
- 已证明：原生VV真实自由积分、一次提交、自然软化/全失效、当前窗口力/历史/动量一致，能量余额0.144631%。
- 唯一dt/2窗口105–122µs局部改善；没有完整时间步认证或后续自动续算。
- 先读本轮 `verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/FINAL_STATUS.json` 与51项REVIEW报告。
- 服务器环境：`source /workspace/y_clot/env.sh`；项目 `/workspace/y_clot/project`。
- 当前C++：`analysis/dynamic_softening_v1/bond_bpm_peri_variational.h/.cpp`。
- 动态协调：同目录 `fix_lps_variational_dynamic.h/.cpp`（内部style LPS/VARIATIONAL/DYNAMIC）。
- 实际隔离库：`/workspace/y_clot/project/verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/build/liblammps.so`。
- 该库SHA：`866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd`；仅研究使用，未写入env.sh或生产入口。
- 冻结数学：`analysis/softening_feasibility_v1/model.py`，SHA `3eb1396b376cf30e40586c26ba649909ab9382d9eac9607666359258273acf74`。
- 动态驱动：`analysis/dynamic_softening_v1/native.py`、`dynamic_run.py`；读事件/能量用 `postprocess.py`。
- 本轮源快照/补丁/复现步骤：verification本轮的source_snapshot、NATIVE_SOURCE_PATCH.diff、REPRODUCE.txt。
- 核心FP64：main_accepted_trajectory.npz；dt_half_window_trajectory.npz；energy_residual_decomposition.npz。
- 只读证据：DYNAMIC_FORCE_ENERGY_AUDIT、STOP_STATE_TOPOLOGY_GRADIENT_AUDIT、RUN0_AND_TRIAL_STATE_SAFETY。
- 图：results本轮OPEN_RESULTS.html（5组PNG/SVG，P3实际停止点云，无P4假图）。
- 参考341点/12047原键；左右49/47受控，245自由，dx50µm，horizon150.5µm。
- G940Pa，K100000/3Pa，rho1000，s0=.15，主sc=.18；完整m/r/taper不变。
- native BPM f为N/m³，新协调fix在POST_FORCE乘V，原真实质量rhoV保持；导出force/U为N/J。
- 只在真实POST_INTEGRATE提交，run0和试探不提交；accepted token在restart延续。
- 上轮受控运动学基线：verification/lps_variational_softening_native_isolated_v1/20261013T201800Z。
- 上轮库SHA b80bfc0a99f4160d62567c0e2d9ddc611004ad40f79c8e49370293c681560773，保留。
- 冻结新旧二值本构力L2差8.814180703270e-6N保留，不靠经验修正消除。
- 本轮未覆盖：卸载段、活动闭合接触、P4、2MPI动态迁移；宏观强度/Gc/生产伸缩性未认证。
- 生产库 tools/lammps-native-y-wall-shared-v1-build/liblammps.so SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a，不变。
- P0审计28候选，安全删除0；UNCERTAIN及历史科研证据均保留。31001历史文件保护通过。
- 实际WSL文档：/home/lzy/projects/clot_rupture.md；before/after/diff/SHA见本轮verification。
- 下一轮唯一优先：同一341点冻结工况的损伤时间步准入与全时程能量验证；先预注册，不自动运行。
- 不必展开旧PROJECT_ACTIVE_CONTEXT及全部历史目录；按需读取以上少量入口即可。
