# 当前clot原生软化开发入口

本轮：lps_native_softening_full_protocol_timestep_admission_v1 / 20261013T224838Z（2026-10-14 Europe/Zurich）。

- **结论：从无损伤t=0开始，dt=0.25µs在第2562步/640.5µs因Δd=0.11363469>0.10停止。完整700µs FAILED_ADMISSION。**
- 最后通过全部冻结门槛的状态是第2561步/640.25µs；越界现场已保存，不继续积分。
- 已覆盖完整升载与峰值保持，卸载350–640.5µs仅部分覆盖；末保持NOT_RUN。
- 原生历史/力/能量/动量在实际窗口一致，能量相对余额最大0.372426823%；不认证全局时间收敛。
- 232.5µs首次真实闭合接触，1576个步有活动接触；最大软罚力重叠24.6649%dx，物理穿透精度未认证。
- 停止时4834条完全失效、7213条存活；仍是一个341点连通体。新裂面/Gc不可用。
- 先读verification/lps_native_softening_full_protocol_timestep_admission_v1/20261013T224838Z/FINAL_STATUS.json和36项REVIEW报告。
- 新驱动/审计：analysis/full_protocol_timestep_v1/；run_protocol.py有单工况启动锁，禁止对封存目录重跑。
- 冻结工况：同目录verification下FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json及其SHA。
- 核心数据：trajectory.npz、damage_history_sparse.npz、keyframes/、actual_contact_pairs.npz、final_executed.restart。
- 关键失败：NUMERICAL_SAFETY_STOP_RECORD.json、FAILURE_EVENT_DIAGNOSIS.json；键1813，自由节点24/245，独立损伤公式完全吻合。
- 生产环境仍为source /workspace/y_clot/env.sh；本轮只在进程内优先隔离库。
- 复用库：/workspace/y_clot/project/verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/build/liblammps.so
- 库SHA：866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd；本轮没有C++修改/编译。
- 实际C++：analysis/dynamic_softening_v1/bond_bpm_peri_variational.h/.cpp及fix_lps_variational_dynamic.h/.cpp。
- 内部style LPS/VARIATIONAL/DYNAMIC，POST_INTEGRATE一次提交，POST_FORCE将N/m³乘V为N。
- 冻结离线数学：analysis/softening_feasibility_v1/model.py，SHA 3eb1396b376cf30e40586c26ba649909ab9382d9eac9607666359258273acf74。
- 341点/12047键，左右49/47规定点、245自由点；rhoV=1.25e-10kg；G940、K100000/3Pa，s0=.15、sc=.18。
- 每侧60µm，半余弦升载300µs、保持50µs、卸载300µs、末保持50µs；冻结物理加载不变。
- 上轮依据：verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/；旧0.5µs到122µs与105–122µs短窗口均保持。
- 本轮不重扫31001历史文件；按上轮195文件payload及直接依赖做限定SHA核查，删除0。
- 根活跃入口是本轮明确更新的可变文件；其上轮SHA副本保存在本轮ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md；旧manifest不改。
- 生产库SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a，env.sh和已有tracked diff均保持。
- 六组PNG/SVG：results/lps_native_softening_full_protocol_timestep_admission_v1/20261013T224838Z/OPEN_RESULTS.html。
- 实际WSL文档：/home/lzy/projects/clot_rupture.md，before/after/diff/SHA在本轮verification。
- 下一轮唯一优先：同一工况的局部损伤事件分步与力—功—能量一致积分，先设计/预注册；不是盲目全局dt=0.125µs。
- P4缺口试样NOT_RUN；不自动进行V3、CFD、微泡、AI、材料Gc标定或生产部署。
- 后续只读上述相关入口，无需展开全部历史报告。
