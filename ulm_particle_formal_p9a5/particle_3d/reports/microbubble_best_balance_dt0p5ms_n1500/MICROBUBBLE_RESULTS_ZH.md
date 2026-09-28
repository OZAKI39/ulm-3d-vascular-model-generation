# 新流场微泡轨迹：1500 条，dt = 0.5 ms

本批使用 best-feasible-balance-v1 新流场；原有 H0/CORE500 结果保留。
流场 SHA-256：`fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2`。

## 实际结果

- 轨迹总数：1500。出口计数：{'O1': 0, 'O2': 133, 'O3': 1182}。
- 终态：{'COMPLETED': 1315, 'SUPPORTED_STATIONARY': 185}。
- 全部完成文件哈希复核：True；数值检查：True。
- 侵壁、约束、入口逃逸、非有限值、出口重分类检查计数：{'penetration_count': 0, 'handoff_violation_count': 0, 'inlet_escape_count': 0, 'nan_inf_count': 0, 'unclassified_corruption_count': 0}。
- 三个出口均出现微泡：False。

| 出口 | 微泡数 | 占全部 1500 条 | 占实际穿越出口者 | 新流场体积流量占比 |
|---|---:|---:|---:|---:|
| O1 | 0 | 0.000% | 0.000% | 15.515% |
| O2 | 133 | 8.867% | 10.114% | 32.912% |
| O3 | 1182 | 78.800% | 89.886% | 51.573% |

流量占比取同一冻结 FEM 场的官方截面通量；它不是微泡计数应达到的目标。
保持本批粒径分布时，未观测到某个出口穿越，只能说明本批有限样本及当前有限尺寸模型的结果，
不能据此断言该血管没有流量，或任何尺寸的微泡均不可能通过。本文没有完成导致出口分布差异的
独立物理归因研究。`data/point_birth_coverage.json` 中的点示踪诊断单独保存，不混入正式轨迹计数。

## 模型、样本及时间

复用 P9-A.1/P6.5 有限尺寸球形微泡、黏性阻力、近壁作用及接触约束积分，
不是点示踪流线。单向冻结流场，微泡之间以及与 RBC 之间没有相互作用。
刚性静止壁面，动力黏度 0.00345312 Pa·s。SonoVue 原始直径分布条件化至 ≤4 μm，
在新 FEM 入口按正通量采样，有限尺寸单次准入；按来源顺序取首 1500 个有效样本，
不重抽被拒来源事件，不按出口筛选轨迹。`proposal_ledger.jsonl.gz` 保留来源及拒绝记录。
用户已明确选择保持原有粒径分布；若某个出口没有穿越，按零计数报告，
不补充小粒径样本以改变覆盖率。原始选择记录见 `data/user_scope_decision.json`。

名义时间步为 0.0005 s，接触和安全检查可细分接受步；单条轨迹观察窗沿用 3/6/12 s。
接触支持静止与长驻留删失分别报告，不把静止当作出口穿越，也不推断为生理滞留。
出口完成定义仍为微泡**中心**首次穿过官方出口截面，不声称整个有限球体已通过。

连续注入沿用标称浓度假设，末个样本出生时刻为 61165.9 s；
这 1500 条是独立轨迹集合，并非同一瞬间同时注入 1500 个微泡。
动画按轨迹年龄对齐回放；穿越出口的标记到达后消失，接触支持静止的终点保持显示。
后者只是显示已计算终态，没有额外推演其后运动。9.1 像素标记用于辨认，不代表实际微泡直径。动画界面已改为英文，血管及微泡主体显示放大 30%，保持原有配色、背景、视角方向及回放方式。

## CPU / GPU 实际分工

8 个 CPU 工作进程执行原有有限尺寸积分和几何/接触安全判据。
累计 CPU 时间 43160.8 s；批次墙钟时间 5864.1 s；
有效平均占用 7.360 个 CPU 核。
RTX 4090 用 float64 实际计算全体四面体梯度与派生场，并复核保存轨迹的位移、速度、
时间步及数据完整性；CUDA 结果见 `gpu_mesh_validation.json`、`gpu_tracks_validation.json`。
**运动积分仍由 CPU 执行，不能称为 GPU 轨迹积分器。** GPU 检查不是微泡动力学精度证明。

CPU 在原有 handoff 函数的私有依赖中使用编译后的点到三角面片距离计算，
几何定义、边/面比较及并列选择规则不变。500 组真实几何查询的
最近点坐标逐位一致，返回权重数值完全相同；严格字节检查发现 10 组权重有 +0/−0
符号差别。实际加速仅用于 `first_handoff_event`，两处调用均丢弃权重，
因此该差别不进入动力学。不能把最初 `np.array_equal` 测试里的 `kernel_bitwise_pass`
名称理解为所有返回字节完全相同；完整核输出的严格字节检查标记为未通过，
实际使用依赖的检查为通过，均保留在 `native_exact_bytes_verification.json`。
首条完整轨迹的 393 个保存状态及其数组科学哈希逐位一致，
完整试算审计科学哈希一致。生产使用的编译源码、适配器和验证证据哈希进入每条完成记录。
初次浮点累加顺序不一致的失败记录保留在 `native_order_initial_failure.json`；
按本机 NumPy 的实际累加顺序修正后重新通过，未将该初版用于生产。
动画由服务器 NVIDIA EGL/OpenGL 渲染，实际 GPU 型号记录在 `render_manifest.json`。
NVENC 的实际初始化返回 `unsupported device (2)`，因此视频编码使用 CPU 的 libx264；
失败日志保留在 `logs/animation_smoke_retry.log`，不把编码描述为 GPU 编码。

## 数据与复现

- `tracks/mb_XXXXXX/trajectory.npz`：原始接受步；列定义见同目录 `trajectory.json`。
- `tracks/mb_XXXXXX/audit.jsonl.gz`：全部接受/拒绝试算与连续接触证据。
- `tracks/mb_XXXXXX/COMPLETE.json`：输入身份与全部结果文件哈希。
- `data/trajectories_dt0p5ms.npz`：每条轨迹按 0.5 ms 重采样的位置；offsets 切分轨迹，
  末次出口时刻可能不是整步。坐标 m，时间 s，半径 m。只对原有分段线性位置插值。
- `data/trajectory_catalog.csv`、`outlet_summary.csv`：逐轨迹及出口统计。
- `data/contract.json`、`cohort.json`、`preparation.json`：参数、固定样本及源码身份。
- `scripts/campaign.py`：prepare / pilot / production；`scripts/gpu_validate.py`：mesh / tracks；
  `scripts/summarize.py`：汇总；`scripts/render_results.py`：图件及入口页。

复现应在新目录设置 `config.json` 的 source_root、flow_path，并依次执行上述阶段。
既有 H0 源码、原始网格、流场及历史轨迹未覆盖，也未启动任何 CFD。
本批固定 dt 的结果不构成时间步收敛证明；每个出口的微泡比例不必等于体积流量比例。
