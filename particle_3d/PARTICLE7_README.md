# Particle-7：持续入口与粒子生命周期

本阶段实现基础设施，保留 P0–P6.5 的科学模型。参数见 `contracts/PARTICLE7_INLET_POPULATION_V0.json`：入口 RBC 体积流率比例 0.45、微泡浓度 8.5e12/m³、确定性累计通量时刻、暂定各向同性 RBC 姿态、六条独立随机流。旧 RBC 平均漂移关闭。

正式 INLET 的线性法向速度在零通量线处裁剪后精确积分。采样使用线性密度的 Dirichlet 混合表示，与按 q/qmax 接受的拒绝采样等价，但每次抽样都接受。正式出生位置没有向内偏移。原始 SonoVue 模块通过已有只读适配器调用；RBC 通过原 P2 分布批量预取，每只粒子的原始几何和来源保持不变。

`PopulationLifecycle` 按出生和出口事件分段，维护分物种 FIFO、稳定 ID、原始预约、实际入场和退出记录。LAMMPS 仅负责存储、邻居和二进制重启，仍只运行 `run 0`；删除使用 `compress no`。checkpoint 保存 P6 自定义属性模式、二进制哈希、全局状态、所有 RNG 完整状态、下一只 RBC、预取缓冲区、pending 候选记录和 ID 墓碑。

验证分为三种场景：

- 原参数的大规模混合通量账本。1000 个 MB 自然对应约 110 万个 RBC，不能靠调 Q 改变此比例。采用宽 200 µm、长 10 nm 的开口控制区和共同平移速度，有限形状可以跨出开口。每次入场均检查精确支撑平面全程壁间隙与可能相交粒子的 P4 exact pair gap；共同平移保持粒子相对位置。这是生命周期压力测试，不是 RBC 通过血管或悬浮液水动力学验证。
- 实际单 MPI LAMMPS 插入、移动、删除、邻居和约 2358 个混合事件的二进制重启验证。采用同一类开口控制区。
- 真实 Frozen INLET 接纳/拥堵 smoke。接纳的粒子作为静态有限尺寸障碍保留；没有为尚未建立的真实混合 RBC 运动发明传播模型。真实流量控制预约时刻，但这项 smoke 不能解释为真实血管内 43 秒的输运结果。MB 实际尝试，失败保留队列；另设独立单 MB 场景，用原 P6.5 求解器在出生后推进 1 ms，验证无人工偏移的实际入管，结果不计入混合拥堵场景；RBC 使用原 P3 自由扁球及胶囊替代。入场不代表真实 RBC passage 已建立。

控制区 Hct 按中心归属统计整个原始 RBC 体积，因此对于跨越开口的有限形状，它是指定的账本诊断，不是几何交集体积分数，也不能外推生理 tube Hct。长测试和真实接纳测试均不强制管内 Hct 等于 0.45。

复现命令（项目根目录，使用本项目 `.venv`）：

```bash
.venv/bin/python particle_3d/scripts/run_particle7_long.py
.venv/bin/python particle_3d/scripts/run_particle7_validation.py
.venv/bin/python particle_3d/scripts/run_particle7_validation.py --real-only
.venv/bin/python particle_3d/scripts/run_particle7_validation.py --isolated-mb-only
.venv/bin/python particle_3d/scripts/analyze_particle7_ledger.py
.venv/bin/python particle_3d/scripts/plot_particle7.py
.venv/bin/python particle_3d/scripts/verify_particle7.py
```

二进制 checkpoint 目录必须是新目录；可给数值验证脚本传 `--output /tmp/particle7-replay/data` 选择新的输出目录，已有 checkpoint 不会被静默覆盖。完整原始长账本为 gzip 压缩 CSV；其余图均有 CSV/JSON 原始数据和绘图映射。全部原始证据位于 `reports/particle7`。不运行 CFD，不选择生产时间步或邻居设置，不实现 PK、非球形润滑或 Particle-8。
