# 当前 clot 活跃开发入口

本轮：LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1；20261013T214800Z。
P0已完成；P1–P5仅按本轮用户授权进入，不自动部署。

- 服务器项目（平铺）：`/workspace/y_clot/project`；环境：`source /workspace/y_clot/env.sh`。
- 先读本文件，再读本轮 `verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/FINAL_STATUS.json`（完成后生成）。
- 已接受原生源码：`/workspace/y_clot/project/verification/lps_variational_softening_native_isolated_v1/20261013T201800Z/source/src/BPM/bond_bpm_peri_variational.{h,cpp}`。
- 上轮可读源码/驱动：`analysis/native_softening_v1/`；完整快照在上轮 source_snapshot。
- 冻结离线数学参考：`analysis/softening_feasibility_v1/model.py`（SHA `3eb1396b376cf30e40586c26ba649909ab9382d9eac9607666359258273acf74`）。
- 上轮正式库：`/workspace/y_clot/project/verification/lps_variational_softening_native_isolated_v1/20261013T201800Z/build/liblammps.so`。
- 上轮正式库SHA：`b80bfc0a99f4160d62567c0e2d9ddc611004ad40f79c8e49370293c681560773`。
- 当前新开发目录：`analysis/dynamic_softening_v1/`（P0审计后才新建动力学功能）。
- 当前新隔离源码/库预定路径：`verification/lps_native_softening_dynamic_commit_and_small_crack_v1/20261013T214800Z/source` 和 `build/liblammps.so`。
- 原始试样：`/workspace/y_clot/project/verification/lps_variational_softening_native_isolated_v1/20261013T201800Z/provenance/`，341节点、12047原始键，dx=50µm，horizon=150.5µm。
- 材料：G=940Pa，K=100000/3Pa，rho=1000kg/m³，s0=.15；本轮主动态sc=.18。
- 直接原生依赖：Bond/Atom/Comm/Domain/Update + 原锁定 pair_bpm_peri；将审查 verlet.cpp、fix_nve.cpp 和相关 hook。
- 原生 f 是 N/m³；物理 F=V*f；完整导出U为J。动力学质量/积分约定须本轮实际核验。
- 已验收：固定历史力/能、受控运动学、sorting/2MPI/restart、1dx闭合接触。
- 待本轮认证：真实自由动力学、自动一步一次提交、动力能量、小型非均匀裂纹。
- 未认证：宏观强度、材料Gc、6–8J/m²目标、V3、生产MPI伸缩性。
- 新模型二值失效不等于旧删键LPS；旧历史力L2差8.814180703270e-6N保留。
- 生产库 `tools/lammps-native-y-wall-shared-v1-build/liblammps.so` 保持隔离，禁止替换。
- 历史verification、source_snapshot、implementation_history均受保护；不因old/test/debug命名清理。
- P0安全删除数量：0。受审资产均有现行/复现价值或用途不确定，不人为制造删除。
- 文档：真实WSL `/home/lzy/projects/clot_rupture.md`；更新须before/after/diff/SHA。
- 仅按需读取上述少量入口；历史PROJECT_ACTIVE_CONTEXT中的旧“下一轮”不代表当前指令。
