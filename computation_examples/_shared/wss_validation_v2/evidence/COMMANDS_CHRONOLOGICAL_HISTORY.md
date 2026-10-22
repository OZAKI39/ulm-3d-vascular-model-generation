# 当前有效范围与复现入口

**用户于2026-09-27取消细档真实血管CFD。** 后文细档/GPU/MPI变体命令仅为已保留的历史执行记录，禁止再次启动细档或新性能实验。用户随后再次取消全部压力敏感性CFD。当前有效主序列仅为中档CPU8/LU完成及独立检查→两档网格比较/图表→报告。后文任何压力实验命令亦为撤销计划，禁止执行。两档血管敏感性与三档圆管验证须区分。

本地恢复入口仍为`watch_vessel_results.py`及`finish_numerical_sequence.py`；它们已移除细档及压力阶段前置条件，顺序门不再生成/注册/启动任何CFD。先检查在运行进程，避免重复启动。细档目录`reports/user_cancellation.json`保存SKIPPED_BY_USER；`run_solver.py`、初始化入口、等待入口及注册入口都会拒绝该取消算例。

# 历史实际命令及证据（按时间记录）

项目：`/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular`。
新增服务器根：`vast4090:/workspace/wss_validation_v2_20260927T1230Z`。

阶段一：
```bash
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B rotate_visualization/prepare_surface_data.py --case /home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1 --output wss_validation_v2/stage1/recomputed
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B wss_validation_v2/scripts/stage1_checks.py
```

阶段二生成：
```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B wss_validation_v2/scripts/build_pipe.py
```
`build_pipe.py` 要求新目录；不得在既有结果上覆盖重跑。复制整份 v2 到另一个目录并修改根路径或先使用新的 case 名。

服务器脚本及 stage2 通过 `rsync -az` 传入独立目录。私有 supervisord 配置保存在 evidence/supervisord.conf，不修改系统服务或正式求解目录。启动与观察：
```bash
ssh vast4090 '/usr/local/bin/supervisord -c /workspace/wss_validation_v2_20260927T1230Z/supervisord.conf'
ssh vast4090 'supervisorctl -c /workspace/wss_validation_v2_20260927T1230Z/supervisord.conf status'
ssh vast4090 'cat /workspace/wss_validation_v2_20260927T1230Z/logs/sequence_status.json'
```

每个 case 的 reports/launch.json 记录实际 MPI、二进制、PETSc 选项、配置和脚本哈希；reports/progress.json 保留正在运行的稳态判据；run/solver.log 是实际求解器原始日志。run/1-procs 同时保存 VTU 与 stFile 检查点。

若连接中断，任务仍受私有 supervisor 管理。恢复时先读取状态；不要再启动正在运行的任务。已通过的 case 会由 run_sequence.py 跳过；失败 case 保留原目录，复制其 SV_MESH、policy 和初始 run 输入到新命名目录后再运行，禁止覆盖失败日志。若使用 checkpoint 重启，需在独立续算目录中设置 Continue_previous_simulation 和正确 restart 路径，并记录初始化方式；本次默认从零初值计算。

阶段二已运行后处理：
```
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B wss_validation_v2/scripts/analyze_pipe.py --case wss_validation_v2/stage2/pipe_nr4
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B wss_validation_v2/scripts/analyze_pipe.py --case wss_validation_v2/stage2/pipe_nr8
```

阶段二增加线性代数后端独立核对（不改 PDE/网格/BC/收敛门槛）：
```
.venv/bin/python -B wss_validation_v2/scripts/prepare_cpu_backend.py --source wss_validation_v2/stage2/pipe_nr4 --dest wss_validation_v2/stage2/pipe_nr4_cpu
```
仅将 PETSc `-mat_type aijcusparse -vec_type cuda` 换为 `-mat_type aij -vec_type standard`，XML 和网格逐字节一致。实际运行由私有 supervisor 的 cpu_backend_check 管理。原 GPU fine 运行保留。

执行方式核对追加命令：
```
.venv/bin/python -B wss_validation_v2/scripts/prepare_asm_blocks.py --source wss_validation_v2/stage2/pipe_nr16_cpu --dest wss_validation_v2/stage2/pipe_nr16_cpu_asm16
```
CPU/GPU 和 ASM 分块均为独立副本；输入哈希、实际运行命令和 PETSc 生效参数分别留在各例 reports/launch.json 和 run/solver.log。

阶段二最细执行后端核对已实际运行：
```
.venv/bin/python -B wss_validation_v2/scripts/compare_backend_snapshots.py --gpu wss_validation_v2/stage2/pipe_nr16 --cpu wss_validation_v2/stage2/pipe_nr16_cpu_mpi8 --step 5
```
三种场相对差均<1e-15。随后仅在原 GPU 独立验证算例写入 STOP_SIM；它不是收敛验证结果。CPU8 主步长继续；`pipe_mpi8_halfdt` 私有 supervisor 程序等待主步长求解成功后再执行半时间步。监视器新增可选 MPI_ranks 和对应 N-procs 路径，旧运行版本保存为 evidence/run_solver_rank1_verified.py，并与原 launch.json 的 SHA核对一致。

当前统一状态查询（区分主动停止的性能试验和真正求解失败）：
```
ssh vast4090 'python3 /workspace/wss_validation_v2_20260927T1230Z/scripts/status.py'
```
最细有效运行是 `stage2/pipe_nr16_cpu_mpi8`，原 `logs/sequence_status.json` 对应已经正常停止的重复 GPU 队列；不要将该旧队列状态当成当前有效计算状态。半时间步 `stage2/pipe_nr16_cpu_mpi8_halfdt` 在主计算结束前为 NOT_STARTED，supervisor 中的等待程序不等于 CFD 已经运行。

生产入口冻结配置绑定负测试及最后回归：
```
.venv/bin/python -B wss_validation_v2/scripts/frozen_configuration_negative_test.py
.venv/bin/python -B rotate_visualization/prepare_surface_data.py --case /home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1 --output wss_validation_v2/stage1/recomputed
.venv/bin/python -B wss_validation_v2/scripts/stage1_checks.py
```
负测试实际同时改变黏度及相应运动黏度（材料内部仍一致），真实入口根据冻结配置哈希拒绝；测试临时目录自动清理，原始 H0 不写入。

验证脚本配置保持性检查已运行：`xml_case` 逐项比较实际 H0 XML，仅允许总步数、保存间隔和显式时间步变化；材料、边界、方程和非线性设置必须相同。独立临时目录检查输出见 evidence/xml_preservation_check.json。

后续新 CFD 的私有 supervisor 注册入口：
```bash
ssh vast4090 '/root/particle8_2_runs/env/bin/python -B /workspace/wss_validation_v2_20260927T1230Z/scripts/register_remote_program.py --name NAME --case stage3/CASE'
```
可用 `--wait-for stage3/PREVIOUS` 串联同一阶段的算例。它只修改本次独立目录的 supervisord.conf，每次保存原配置；已有算例日志存在则拒绝覆盖。实际注册记录为服务器 logs/registered_NAME.json。

阶段二最细半步长已运行生产后处理、独立日志解析、mesh_metrics.py、collect_results.py 和 pipe_additional_checks.py；与主步长逐场比较见 data/pipe_timestep_comparison.json。阶段二完成标记见 stage2/completion.json。

阶段三已实际执行：
```bash
.venv/bin/python -B wss_validation_v2/scripts/prepare_baseline.py
.venv/bin/python -B wss_validation_v2/scripts/analyze_vessel.py --case wss_validation_v2/stage3/vessel_baseline
.venv/bin/python -B wss_validation_v2/scripts/prepare_baseline_execution.py
.venv/bin/python -B wss_validation_v2/scripts/prepare_vessel_mesh.py --case wss_validation_v2/stage3/vessel_medium --edge-factor .75
.venv/bin/python -B wss_validation_v2/scripts/prepare_vessel_mesh.py --case wss_validation_v2/stage3/vessel_fine --edge-factor .5
```
所有 Python 3.13 命令均使用 PYTHONDONTWRITEBYTECODE=1、OPENBLAS_NUM_THREADS=1；网格生成脚本自行调用官方嵌入 Python。每个新增算例整体通过 rsync -az --relative 同步到服务器同名 stage3 路径。supervisor 程序 vessel_baseline_cpu、vessel_baseline_halfdt、vessel_medium 依次运行；等待进程不等于对应 CFD 已开始。实际状态见每例 reports/launch.json、progress.json、execution.json。

阶段三执行安排更新：GPU8独立试验首轮DIVERGED_BREAKDOWN，进程组已停止、原始日志保留。半步长另用原H0 GPU1配置，实际case为stage3/vessel_baseline_gpu_mpi1_halfdt，与原vessel_baseline比较。CPU8基线继续；待它完成，vessel_medium_after_cpu启动中网格，vessel_fine等待中网格完成。原vessel_baseline_halfdt和vessel_medium只是旧等待程序，已停止并将autostart设为false；CPU8半步长明确未执行。配置变更前文件保存在服务器logs/supervisord_before_halfdt_gpu1.conf。

主三档血管网格实际执行安排（最终生效）：stage3/ORIGINAL_GPU1_MAIN_SEQUENCE.md。两档新网格在CFD启动前通过 select_original_gpu1.py 将执行后端改回原H0 GPU1；reports/unexecuted_cpu8_inputs 保存旧输入，只有policy和PETSC_OPTIONS哈希改变。旧CPU中/细等待程序禁用，vessel_medium_gpu1等待vessel_baseline_gpu_mpi1_halfdt，vessel_fine_gpu1等待vessel_medium。CPU8原网格一致性复算独立继续。严禁将未运行的CPU8中/细输入记成实际结果。

阶段三附加诊断已在基线上实际运行：
```bash
.venv/bin/python -B wss_validation_v2/scripts/diagnose_wall_stencil.py --case wss_validation_v2/stage3/vessel_baseline
.venv/bin/python -B wss_validation_v2/scripts/control_volume_checks.py --case wss_validation_v2/stage3/vessel_baseline
.venv/bin/python -B wss_validation_v2/scripts/compare_backend_snapshots.py --gpu /home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1 --cpu wss_validation_v2/stage3/vessel_baseline_cpu_mpi8 --step 10 --output wss_validation_v2/data/vessel_backend_step10.json
```
watch_vessel_results.py 正在本地监视四个指定新算例完成状态；只有solver执行PASS后才拉取并依次运行正式WSS入口、独立日志/区域/截面分析、mesh_metrics、wall_stencil和CSV汇总。其状态在logs/stage3_analysis_watch.json；原始命令输出在logs/watch_vessel_results.log。它不会启动阶段四，不会把等待中的算例当作已求解。

半步长完成后监视器分类修复的实际命令：
```bash
.venv/bin/python -B wss_validation_v2/scripts/test_log_acceptance.py
.venv/bin/python -B wss_validation_v2/scripts/reassess_recovered_execution.py --case wss_validation_v2/stage3/vessel_baseline_gpu_mpi1_halfdt
# 同一复核在服务器独立副本执行；后续监视器及 flow_parser 已同步
ssh vast4090 'supervisorctl -c /workspace/wss_validation_v2_20260927T1230Z/supervisord.conf start vessel_medium_gpu1'
```
复核只能用于首次 FAIL、正常完成且全部证据满足条件的已有结果；重复运行会拒绝覆盖首次记录。详情见 stage3/MONITOR_RECOVERY_REASSESSMENT.md。中档已实际启动，细档保持同阶段队列；本地 watcher 已重启继续拉取与分析。

中档GPU1实际失败后的当前有效恢复流程：
```bash
.venv/bin/python -B wss_validation_v2/scripts/restore_cpu8_after_gpu_failure.py
# 该脚本在服务器同名独立目录也执行；拒绝覆盖已保留的失败副本。
ssh vast4090 '/root/particle8_2_runs/env/bin/python -B /workspace/wss_validation_v2_20260927T1230Z/scripts/register_remote_program.py --name vessel_medium_cpu8_retry --case stage3/vessel_medium --wait-for stage3/vessel_baseline_cpu_mpi8'
ssh vast4090 '/root/particle8_2_runs/env/bin/python -B /workspace/wss_validation_v2_20260927T1230Z/scripts/register_remote_program.py --name vessel_fine_cpu8 --case stage3/vessel_fine --wait-for stage3/vessel_medium'
```
实际GPU1失败副本及CPU8输入恢复证据见 stage3/CPU8_MAIN_SEQUENCE_AFTER_GPU_FAILURE.md。旧GPU程序不再使用。

后续代数执行检查和当前主配置见`stage3/EXECUTION_UPDATES.md`。已实际运行prepare_subdomain_lu_check.py、prepare_ilu6_check.py、prepare_lu16_check.py、configure_mpi16_launch.py、prepare_fieldsplit_check.py及fix_fieldsplit_configuration.py。所有同场比较使用：
```bash
.venv/bin/python -B wss_validation_v2/scripts/compare_backend_snapshots.py --reference CASE_A --candidate CASE_B --step 5 --output OUTPUT_JSON --change-description '实际执行配置差异'
```
CPU8/LU及CPU16/LU第5步均通过预登记1e-9相对L2门槛。`adopt_subdomain_lu.py`在本地/服务器执行，归档未完成的中档ILU2例，保持网格/XML/物理/dt/容差；目前实际中档supervisor程序为`vessel_medium_LU`。细档尚未启动。分块实验为`vessel_fieldsplit_check`，尚未采用。MPI清理使用与本次case/run完全一致的cwd筛选；新监视器在收到停止信号时会清理该例求解器进程，不能仅假定停止监视器就停止了独立MPI进程。

细档最新运行门：`run_fine_after_medium.py`（服务器supervisor仍名vessel_fine_LU16），先等中档execution PASS和本地独立检查同步的independent_acceptance.json，再调用prepare_initial_guess.py生成明确初值，最后调用run_solver.py真正求解。源码读取/分发链哈希见evidence/initialization_source_identity.json；映射同网格实测见data/initial_guess_identity_test.json。源/目标hash、几何外部投影点和距离、只增加两个初始化XML标签的语义核对均写入细档reports。完成后compare_initial_final.py量化实际CFD对初值的修正。

阶段三/四顺序执行门（已启动；在本地V目录运行）：
```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 ../.venv/bin/python -B scripts/watch_vessel_results.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 ../.venv/bin/python -B scripts/finish_numerical_sequence.py
```
第一条只独立接收/检查阶段三；第二条等待其COMPLETE并再次核对实际解、日志、输入锁后才登记/启动阶段四。恢复时先查进程，禁止重复启动正在运行的watcher或顺序门。两个脚本按已保存完成证据续接，不覆盖求解日志。细档初值登记是此文较早“默认零初值”说明的明确例外。
