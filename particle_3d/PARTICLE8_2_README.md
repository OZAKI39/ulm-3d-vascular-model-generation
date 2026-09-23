# Particle-8.2：多出口诊断与远程并行计算审计

本阶段先做来源审计，再按 A→B→C→D→E→F 顺序推进。所有正式轨迹与 point-tracer 批次必须在 `root@50.115.148.16:4159` 对应的、经 SSH 实测的服务器主机运行。本地 WSL 只开发、测试、至多 10 条 smoke、调度、同步和审核。不会在本地补跑远程大任务。

当前进度与阻塞以 `reports/particle8_2/PARTICLE8_2_VALIDATION.json` 为准，不能把代码已实现视为远程科学验证已通过。此前 P8.1 的环境文件虽记录 WSL 平台，但没有逐轨迹主机回执，所以远程重计算来源记为 `NOT_PROVEN`。

## 模块

- `particle82_provenance.py`：原始资源命令、容器配额、主机身份及远程运行门禁。
- `particle82_batch.py`：Linux fork/COW、按 ID 独立积分、原子文件与提交回执、恢复、稳定排序、RSS/PSS/CPU/swap 记录。
- `particle82_integration.py`：沿用 P8.1 的完整积分顺序、P6.5 推进器与精确查询缓存。仅参数化计算预算 factor=1/4/16；不修改力、尺寸、入口点或接受规则。
- `particle82_tracers.py` / `particle82_point_native.py` / `.cpp`：原 Frozen P1 数据的双精度零半径 RK23 诊断。显式四面体邻接定位，局部误差控制；出口由实际 ODE 末段与原端盖再次分类。早期 VTK 流线实现保留为失败对照，不用于正式微泡。
- `particle82_diagnostics.py`：原 P8.1 每个入口 proposal 的直接 basin 映射、尺寸分组、720 个停止状态字段、独立的新自然出生账本。

零半径诊断有 30 s 时间预算、100000 个接受状态预算及 2 mm 路径长度预算；未解析路径保留在总分母中。局部 cap 末步最多为配置空间步长的 2%，由原速度方向推进，不按目标出口连接路径。其 RK23 数值误差与有限尺寸 P6.5 微泡积分是不同的验证对象。

## 可复现命令

来源审计与精确提交部署：

```bash
.venv/bin/python particle_3d/scripts/audit_particle82_compute.py
.venv/bin/python particle_3d/scripts/deploy_particle82.py --run-id diagnostic_pipeline
```

部署器只读取已提交源码，创建 `/root/particle8_2_runs/<commit>/<run_id>/repo`，不覆盖既有目录。Frozen 的输入清单还引用 `configs`、`logs`、`reports` 中的少量证明文件，必须一并同步并验证散列。SonoVue 外部冻结项目已按其 `SHA256SUMS` 同步至远程只读输入目录；旧模块所需绝对路径通过独立符号链接访问同一输入。

服务器独立环境为 `/root/particle8_2_runs/env`。命令应由 tmux 持有，避免 SSH 断开终止任务。以下命令仅在远程快照的 repo 内运行，`../host_provenance.json` 来自实际 SSH 主机回执：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1 \
  /root/particle8_2_runs/env/bin/python -Bu particle_3d/scripts/run_particle82_points.py \
  --count 100000 --workers 8 --shard-size 256 \
  --host-provenance ../host_provenance.json --output-dir ../point_basin_100000

/root/particle8_2_runs/env/bin/python particle_3d/scripts/benchmark_particle82_server.py \
  --input-ledger /path/to/original/birth_ledger.json --output-dir ../scaling \
  --host-provenance ../host_provenance.json --resume

/root/particle8_2_runs/env/bin/python particle_3d/scripts/run_particle82_batch.py \
  --input-ledger /path/to/ledger.json --output-dir /path/to/run \
  --workers 8 --id-start 1 --id-end 64 --shard-size 1 \
  --config /path/to/config.json --host-provenance ../host_provenance.json --resume
```

配置必须包含 `dt_s`、`guard_factor` 和 `dataset_role`。正式数据角色为 `NATURAL_FLUX_WEIGHTED_DATASET`；benchmark、basin-conditioned、continuation、timestep 数据必须另设角色与目录。示例 workers=8 不代表已经测得最佳 worker 数。

64-ID scaling 在 1/2/4/8/16 workers 上运行，受实际可用 CPU 上限限制。选择完成且 cgroup 峰值内存低于 90% 上限的最高吞吐量。恢复读取速度不能作为轨迹计算吞吐量；完整已有 benchmark 保留原性能记录，未完成试验保留后重新开始完整计时批次。

每条已提交轨迹拥有原始 NPZ、JSON 及独立主机回执。恢复验证源码、配置、Frozen 输入、出生属性和样本散列；写完元数据但未完成回执事务的 ID 属未提交任务，可确定性重做。合并按 stable ID 排序；科学记录散列与 PID/worker/执行时间分离，因此并行完成顺序不会改变科学目录。

## 历史回归与科学边界

P1/P6.5/P7/P8/P8.1 的回归在独立历史 Git 副本中、分别启动 pytest 执行，以保留依赖提交历史并避免旧阶段同名测试模块和绘图全局状态相互影响。不得写入正式 P8.1 输出。

仍采用 Frozen flow、SonoVue 原分布、P7 累积通量时钟、H_D=0.45 feed 参数、C_MB=8.5e12 m⁻³ 和 isotropic V0 假设。计算停止不是生理捕获。真实 RBC passage、完整耦合悬浮液、PK、生产 timestep/neighbor 参数均未建立。本阶段不会启动 Particle-8.3、Particle-7.5 或新 RBC physics。

## 结果导出与审核

`export_particle82_results.py` 从远程自然集合生成同一保存场景、全部接受状态 CSV.gz、出生/终止事件和逐轨迹来源。`plot_particle82.py` 生成 00–09、11，`render_particle82.py` 生成四组 MP4/GIF；`audit_particle82_media.py` 解码所有视频与 GIF 帧并生成 10/storyboard，`audit_particle82_render_determinism.py` 新建两个窗口复现同一活动状态。

永久证据测试通过 `PARTICLE82_RESULTS_ROOT`、`PARTICLE82_POINT_ROOT` 绑定真实远程文件；没有提供证据目录时，这些测试明确 skip，不计作证据已通过。最终完整测试不得 skip。`finalize_particle82_report.py` 只在所有汇总、媒体审核、真实测试记录和代理图片检查完成后生成最终中文报告及验证 JSON。

点路径原始 NPZ 保留完整时间/位置，不做显示抽稀。它们会产生大量文件缓存；如基准前缓存导致 cgroup 内存门槛触发，可使用 `release_particle82_read_cache.py` 对本阶段已落盘的 point NPZ 做只读 `POSIX_FADV_DONTNEED`，记录前后 memory.stat，不删除数据、不改系统 drop_caches。受影响的基准另存后重新完整计时。

`sync_particle82_evidence.py --mode summaries` 同步摘要；`--mode review` 同步自然数据与审核材料；`--mode all` 恢复完整远程原始诊断数组。任何未同步的大文件必须在远程散列清单中说明准确位置，不能声称其已在本地。

### 历史停止轨迹的真实续算

服务器 factor=1 从出生重算与旧 P8.1 的 64-ID 控制并非全部逐字节相同，两个病例的状态数也不同，故这些重算不能冒充历史 continuation。`particle82_checkpoint.py` 从旧 NPZ 的最后已接受位置、速度、角速度、q 和原 radius 恢复；原样保留整个历史前缀。根据末个已接受二分叶节点恢复尚未执行的右子区间端点与深度，保持原 nominal dt 网格、P3 递归算术和 P6.5 trial/接受规则。全部 720 个历史记录的二分区间结构均可恢复。

原 rejected trial 不消耗物理时间，可重新尝试；累计 provider 次数恢复，未接受调用的连续计数从恢复点重新开始。该计算开销差异明确记录，不宣称恢复所有旧调用栈缓存。后缀首步和全部后缀都再次检查 `x_next = x_old + dt*v`，确保从真实保存状态开始推进，禁止接上另一条从出生重算的路径。串并行 benchmark 仍比较同一服务器上的完整重算，不与历史 checkpoint 角色混用。

### 自然批次的计算预算选择

原计算预算固定 natural scheduled=5000、最多 5000、自然批次最多 3 小时、总远程工作最多 8 小时。自然出生账本生成前，用固定 64-ID 基准与旧 720 条续算的实测 wall time，预测候选 guard=1/4/16 的计算耗时：`5000 / baseline_tracks_per_hour + (5000/2200) * continuation_wall_hours`，并加 15% 时间余量，选择预计能满足 3 小时预算的最大已测 factor。此规则在新自然样本生成前记录；不会查看新出口结果来修改预算或重抽样。原最高 guard=16 的诊断仍完整运行，生产 dt/physics 没有改变。该预测不保证运行时间，实际 3 小时限额仍生效。
