# Particle-8.1 / Full-vessel ULM microbubble trajectory simulation

本阶段在完整官方 Frozen 几何与静态 FEM 流场中独立积分 SonoVue 微泡，保存从入口面到首次正式 outlet 中心穿越的完整轨迹。未完成轨迹和入口拒绝保留独立终止原因。结果不等于 RBC-coupled physiological suspension；没有 RBC 输运、MB–MB 耦合、完整 PK、新 CFD、Particle-7.5 或新变形物理。

输出隔离在 `particle_3d/outputs/particle8_1/`。旧 P0–P8 / full3d 文件逐项 SHA256 锁定。Frozen 与 SonoVue 输入沿用原 manifest，所有单位为 SI，几何、半径、墙和 flow 不移动、不调整。

主要模块：

- `particle81_simulation.py`：P7 通量时钟、原 SonoVue inverse CDF、Haar SO(3) 姿态、有限尺寸入口接纳，以及原 P6.5 stepper 的逐接受子步记录。
- `particle81_cache.py`：用既有 P6 private dependency binding 缓存相同位置/半径的 wall gap 与完全相同的 BVH 请求；不修改任何上游模块全局；每个 stepper 的查询适配器共享只读原始网格，仅缓存精确请求的返回值。可选 global distance motion bound 仅作独立等价性实验，正式数据默认关闭。
- `particle81_replay.py`：轨迹 catalog、事件 ledger、只读场景、物理时钟插值与显示帧时间表。
- `particle81_visuals.py` / `particle81_figures.py`：完整 3D vessel、轨迹累积、代表性旅程与十张静态图。

从仓库根目录复现：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python particle_3d/scripts/run_particle81_dataset.py --count 2200 --workers 16
.venv/bin/python particle_3d/scripts/export_particle81_scene.py
.venv/bin/python particle_3d/scripts/audit_particle81_physics.py
.venv/bin/python particle_3d/scripts/audit_particle81_sensitivity.py
.venv/bin/python particle_3d/scripts/benchmark_particle81_cache.py
.venv/bin/python particle_3d/scripts/plot_particle81.py
.venv/bin/python particle_3d/scripts/render_particle81.py
.venv/bin/python particle_3d/scripts/audit_particle81_media.py
.venv/bin/python particle_3d/scripts/verify_particle81.py
```

`--count` 是累计预约数量，不是成功数量。可增加该数延长物理采集；已有 ID、属性、时间前缀不改变。各 worker 的位置随机流由 `[master seed, stable ID, 81]` 确定；并行处理顺序不影响输出物理序列。恢复时核对输入、积分设置和样本文件散列，拒绝不一致续算。

`trajectories/mb_XXXXXX.npz` 保存每一个原推进器接受子步：elapsed time、位置、区间所用速度/角速度、q、原始 wall gap、g_NF、handoff 下限、近场状态和速度评估时刻。绝对物理时钟是 `birth_time_s + elapsed_time_s`；使用 elapsed time 推进，避免长采集时钟损失微小子步精度。初始行速度是自由场诊断；后续行速度属于刚完成的区间，不误标为终点重新求解速度。

同名 JSON 保存原属性、候选入口点、面 ID、barycentric 坐标、随机状态、所有终止原因、出口三角形与完整路径统计。`data/particle8_1_trajectory_catalog.json` 是图和动画的唯一场景索引，`particle8_1_events.json` 是按绝对时钟排序的事件记录。图像源码侧车和逐帧 manifest 绑定同一场景散列，显示代码不抽样新的粒子。

固定设定：H_D=0.45 是保留的 feed 参数，当前独立 MB 积分不包含 RBC，也不强制 tube Hct；C_MB=8.5e12 m⁻³ 是沿用的名义浓度；P7 deterministic cumulative flux accumulator 决定出生时刻。初始姿态采用原 P7 isotropic V0 方法，是模型假设。浓度固定在很长采集窗口内仅是本计算设定，不表示一次实际注射后的 PK 时程。

数值推进仍是 P6.5 球形法向近场阻力、2 nm / 10⁻³a handoff、原 P3 时间二分和连续壁间隙证书。名义 dt=0.25 ms 沿用 P7 验证步长；production timestep、neighbor 参数和非球形 lubrication 未冻结。每只最多 1.5 s 物理驻留、16000 次 provider 调用；另记录精确静止、舍入尺度停滞和无接受状态的拒绝试探保护。最近 64 个接受子步若总共推进不足名义 dt 的 1%，也明确停止为过度细分进度不足。安全停止不等于真实生理捕获，不能计为已穿行轨迹。所有 guard 是计算停止条件，不修改力、位置、半径或已接受子步；这些条件会影响完成比例，不能把完成/停止频率当作生理概率。

动画 01/04 慢放真实到达事件窗口，窗口间空闲时间明确 CUT；02 压缩整个 acquisition clock，历史层只加入已退出的完整路径；03 慢放展示代表性 MB 从 birth 到 outlet 的完整旅程。显示仅做保存状态间线性位置插值和原 P8 SLERP，原始轨迹不抽稀、不平滑，不外推到出口，不旋转物理坐标。相机正交缓慢绕行，原尺寸球体表示当前 MB；历史线按真实最终出口着色。长窗口下经常没有 active MB，不能用同步多泡云团替代真实调度。

静态 localization 图使用原接受位置并按对应物理 dt 加权，避免细分步密集造成假密度。这是 ULM-style occupancy，可视化不包含超声 PSF、测量噪声、检测漏失或声学重建。

最终交付由 `finalize_particle81.py --source-commit <代码提交>` 生成，要求测试、视频解码及实际图片检查记录有效；用户人工验收保持待审。

独立计算分区可使用 `--first-id`，仅用于不重叠 ID 区间；完整 birth ledger 仍保存从 1 起的固定前缀。正式目录只包含 `trajectories/`，诊断参考与版式预览不计入队列。

停止后的微泡状态未知；回放只在有效计算区间内显示 active 轨迹，停止后不外推，也不把计算结束解释成物理排出或死亡。

单 MB journey 在已计算出口端点保持 12 个显示帧（15 fps 下 0.8 s），物理时钟不推进。GIF 每五帧并保留首/中/末关键帧，MP4 与逐帧 JSON 为时间依据。
