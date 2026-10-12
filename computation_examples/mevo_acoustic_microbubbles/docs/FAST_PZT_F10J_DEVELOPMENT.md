# F10J — GPU B2 Solver Convergence and Execution Efficiency

The original GPU failure is explained and corrected without changing physics or tolerances. Corrected and optimized fixed-state solvers pass the registered numerical gates. **The optimized GPU is still slower than the matched CPU reference, so production integration is stopped.**

Baseline complete F10I commit: `209b9a77ca8e114e0f3258a3d2a176cbec41c611`; old GPU v2 execution `458dabada891b218d8dc91b256a6de6c3d7d496d`. Independent branch `dev/fast-pzt-f10j-gpu-b2-convergence-20261012`; WSL project `/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10j_20261012/computation_examples/mevo_acoustic_microbubbles`. Remote `/workspace/mevo_acoustic_microbubbles/f10j`; immutable `.execution/<commit>` snapshots, unique Supervisor labels. F10I's 71 manifest entries and all 13 original F5 input files are unchanged and hash-verified. Old compatibility symlinks and unrelated projects were not altered.

## Deliverables and actual execution

`src/mevo_mb/fast_sim/f10j_gpu.py` is an opt-in subclass, not a copy of the solver or production engine integration. `scripts/fast_pzt_f10j/` contains audit, numerical validation, fair timing, offline topology, figures and delivery scripts. Configurations separately preregister diagnosis, correction, optimization and timings. Complete results, negative cases, raw PCG/Armijo traces, logs, input/source hashes and actual execution IDs are in `results/fast_pzt/f10j/`. The final postprocessing commit is not substituted for a scientific execution commit.

- `performance_v1`: `336e28dea7299583f7fec7458d030c0165f771e2`.
- `optimized_v2`: `c986ef9b8b4b70302f6a82f11e1e72e72e8e033a`.
- `guards_v1`: `1b2b7cd4e2df555d340110addb16132125502681`.
- `audit_v1`: `e5688f0998c28a9de2eb5199602349716040d9d0`.
- `process_v1`: `0e8fa5f4cf6722980bd483dd75c0ca094ae42ca8`.
- `corrected_v1`: `8e9d9af5ec9d93d0c368dfc9e47393cf0ef9b950`.

## Layered decisions

|Layer|Decision|
|---|---|
|Mathematical stopping, this fixed state|PASS; original 1e-14 m bound and original true/recursive PCG targets|
|GPU implementation, tested fixtures|PASS; 16 legal small states, 2 expected budget rejections, 1 real32; optimization bitwise equivalent|
|CPU/GPU performance|NO B=1 ADVANTAGE; GPU 22.048 / CPU 11.867 ms|
|Contact / near-field physical calibration|PHYSICS_UNCALIBRATED|
|Group position and shape|Prior F10H/I evidence retained, not newly expanded|
|Strict edge timing|Original 12/29 over 10 us, max130 us FAIL retained|
|Persistent split/merge prediction|PERSISTENT TOPOLOGY VALIDATION NOT RUN|
|AI data|training_ready=false; no production integration|

## Historical topology screening

All analysis uses existing trajectories, fixed IDs, original 2.1R graph and 2 ms dwell. 98 historical ledgers reconstruct consistently. Six failed/terminated files have incomplete reconstruction coverage and are explicitly excluded; other skipped paths include empty schemas, multi-batch data and replicated throughput workloads. Skipped or incomplete data are not labelled as having no events. F10B has no additional trajectory.h5 in this archive; its references are inherited data. Full path/hash, model, accepted-frame scope, topology counts and exclusions are in TOPOLOGY_SCREEN.json. Histories are correlated runs, not independent samples.

F5 main contains 31 persistent aggregation merges before release, then 11 persistent splits and 2 merges after switch-off; one late split remains right-censored. F7 mature histories show no persistent events in examined records; F8 contains historical flow-release changes. F10D–H preparation variants have model-dependent merging. All three actual geometric B2 histories (F10H10/5 us, F10I tightened2ms) remain one connected group and provide no positive persistent-event accuracy test. Frozen-normal F10G B2, PBD M0 and geometric B2 are explicitly different models.

Candidate source: `results/fast_pzt/f5/science_v1/main/trajectory.h5`, stored frame835 at417.5ms, all32 IDs active, finite positions, maximum overlap 0.002467nm. Original M0 first persistent split at418.54ms. Proposed future examination window417.5–421.54ms. PZT/secondary off, original CFD on, original triangular wall; exact initial arrays and IDs in HISTORICAL_F5_CANDIDATE.npz, controls/inputs/hashes in TOPOLOGY_CANDIDATE.json. Minimum stored wall gap is approximately zero, so near-wall coupling and valid B2 restart must be checked first. This is a historical geometry candidate, **not an admitted B2 validation case or B2 truth**. No event is manufactured and no new topology trajectory was run.

## Reproduction and tests

Use the existing environment; remote jobs run committed snapshots through `python.sh` (BLAS threads1). Executed modules: `audit --label audit_v1`, `validate --label corrected_v1`, `validate --label optimized_v2 --optimized`, `performance --label performance_v1`, `cold_process --label process_v1`, `guards --label guards_v1`. Each label is immutable; choose a new label for a separately authorized rerun. Numerical tests run before optimization/performance. Old failure reproduction and GPU merit-kernel vs independent extended-precision tests are actual computations. No pytest that merely mirrors code is substituted for scientific checks.

Offline modules: `analyze_audit`, `topology_screen`, `figures`, `delivery`, `progress_document`, `verify`. Do not rerun scientific histories to regenerate plots. See the dedicated [solver audit](FAST_PZT_F10J_GPU_SOLVER_AUDIT.md) and [performance report](FAST_PZT_F10J_PERFORMANCE.md).

## Required Chinese answers

1. 原GPU卡在基础PCG、第一次方向PCG各自残差，以及最终原保守梯度界；末次Armijo把合格下降步拒绝了。
2. 递推与真实残差基本一致。基础为1.25499e−12／1.25371e−12，目标5.18976e−13；没有发现灾难性递推漂移。
3. 1.850e−14 m不是实际位置误差。旧程序停在尚未接受最后一个下降步的状态；不需要放宽1e−14门槛。
4. 两个接近的总能量计算出现舍入；完整候选被算成增加1.02349e−16，稳定差值实际减少约1.32813e−17，20个候选均被原比较拒绝。
5. 发现数值实现问题：基础及首方向未收敛仍继续、相近总能量比较不可靠。原物理方程及停止目标没有改错或放宽。
6. 修正后32泡通过全部登记的数学、参考和预算检查。保守界7.885585e-16m；五次全步更新。
7. 相对原CPU位置／力差2.605196e-15m／1.319722e-18N；相对独立参考1.155703e-16m／5.851729e-20N；真实残差1.107183e-16m。
8. 中途标量读回35→17→12；包含最终向量证据为47→31→26次数组读回。没有声称消除所有同步。
9. 本轮三组完整成本中位数：旧失败GPU22.789ms、修正版42.256ms、优化版22.048ms。构图和独立进程另报，未冒用历史35.457ms。
10. 没有公平效率优势。CPU11.867ms，优化GPU成本是其1.858倍。
11. 本轮停止GPU B2生产集成，保留独立模块和CPU有限状态研究参考，不继续昂贵长轨迹或大规模重写。
12. 找到F5关声后417.5ms候选，但源自M0且贴壁；需要另行验证B2合法重启，不能称已有B2持续分裂真值。
13. 原团级和逐泡位置／形状指标的短窗数值证据保留；严格连边时间仍失败，物理校准、持续拓扑预测与AI准入未完成。
14. corrected_v1与optimized_v2获固定状态数值准入；原F10I仍失败；负例正确拒绝；性能和生产准入未通过。
15. 四组实际英文PNG/PDF位于 `/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10j_20261012/computation_examples/mevo_acoustic_microbubbles/results/fast_pzt/f10j/figures`；来源与人工阅图记录同目录。
16. 科研进度文档由独立脚本最小补充，备份为 `/home/lzy/projects/microbubble_aggregation.md.bak_F10J_20261012`；before/after、diff、SHA、标题结构检查保存在progress_document/。
17. F10K优先用CPU数值参考准备可复现的短窗持续拓扑与近壁耦合验证；先做重启和预算准入，不把M0事件时刻作为B2拟合目标，也不继续追逐130us瞬时边时差或未经验证的GPU生产吞吐。
