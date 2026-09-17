# Revised v6 连续微泡灌注实现 Review Summary

> 状态说明（2026-07-15）：本文件记录的是此前“预平衡后开始正式记录”的实现。该初始条件已经按用户最新要求被“正式仿真从空血管开始”取代；当前行为请参阅 `review_summary_empty_lumen_formal_start.md` 和 `REPRODUCTION_SPEC.md`。

审查日期：2026-07-14  
验证环境：`D:\anaconda3\envs\pmp\python.exe`

## 1. 本次实现目标

本次修改将 `generate_microbubble_trajectories.py` 的生产粒子链路从“固定初始数量、粒子离开后按显示数量补充”的模式，改为 `Revised v6.md` 规定的确定性连续灌注模型。生产运行不再提供 `passive_tracer` 模式，只允许有限尺寸、近壁迁移率粒子模型。

入口物理输入改为微泡数浓度 (C_0)。默认配置为 (5\times10^5\) MB/mL；将同一项改为 (1\times10^6\) MB/mL 即可运行较高浓度条件。代码将 MB/mL 乘以 (10^6) 转换成 bubble/m³，再结合有效厚度、已收敛入口速度、有限尺寸可通过区域和既有粒径分布，计算真实注入率。域内粒子数不再是输入参数，而是注入率、停留时间、出口离开和入口等待共同产生的结果。

## 2. 修改了哪些文件，每个文件为什么改

### 入口、配置与主调用链

- `configs/physics_flow_config.yaml`
  - 新增 `inlet_number_concentration_mb_per_ml` 和预平衡安全上限。
  - 默认 (C_0=5\times10^5\) MB/mL，并明确记录 (1\times10^6\) MB/mL 的切换方式。
  - 删除生产 YAML 中固定场景粒子数、初始批量采样、离场补充、随机种子和 `motion_model` 选择项。
- `utils/config.py`
  - 解析并校验正的有限浓度和预平衡上限。
  - 从公开动力学配置对象中删除 `motion_model`；旧配置显式写入 `passive_tracer` 会报错。
  - 旧配置若显式写 `mobility`，仅作为迁移期兼容键接受，不再形成模型分支。
- `generate_microbubble_trajectories.py`
  - 更新入口说明，使其准确描述连续灌注和唯一公开的迁移率模型。
- `utils/particle_advection.py`
  - 删除公开的被动示踪分发；已收敛流场只进入迁移率连续灌注链路。
  - 将有效厚度、入口边界深度和流量容差传给灌注模块。
- `utils/runner.py`
  - 接通新的物理场配置和连续灌注入口。
  - 结果元数据改为记录浓度、注入率、预平衡、真实最小/最大/平均粒子数和入口等待等信息。

### 新增的连续灌注模块

- `utils/particle_inlet_flux.py`（新增）
  - 在每个 root 血管的开放边界下游建立固定注入截面。
  - 计算 (u_n^+(s))、有限半径可通过条件和联合数通量 (j(s,R))。
  - 完成 MB/mL、m³ 和微米制流量之间的显式单位换算，并构造粒径边缘 CDF 与给定粒径下的位置条件 CDF。
- `utils/particle_perfusion_schedule.py`（新增）
  - 用 (t_k=(k+1/2)/\lambda) 生成与时间步无关的等通量事件。
  - 用无随机种子的二维 Halton 序列确定每个永久 ID 的不可变粒径和固定入口位置。
- `utils/particle_equilibration.py`（新增）
  - 同时检查出口数通量与域内人口趋势，不新增独立百分比容差，直接复用已接受流场的流量容差。
  - 用累积停留时间和二分查找避免随预平衡步数增长的重复全历史扫描。
- `utils/particle_perfusion_transport.py`（新增）
  - 从空血管开始预灌注；首次满足双平衡条件的状态直接继承为正式 (t=0)。
  - 在内部时间步中按物理事件时刻拆分推进，不把注入时刻量化到保存帧。
  - 若指定位置发生有限尺寸重叠，粒子保留原 ID、粒径和位置参数，在入口上游等待；不改位置、不重采样、不删除、不强插。
  - 粒子状态槽按永久 ID 单调增加且永不复用。
  - 输出变长帧记录、生命周期时间、预平衡和稳定性诊断。

### 复用、输出与可视化兼容

- `utils/particle_mobility_transport.py`
  - 公开函数在正浓度生产配置下接入新连续灌注驱动，同时复用原有迁移率 RHS、碰撞、Heun/Euler 和连续有限尺寸几何约束。
  - 零浓度固定人口路径仅保留为旧单元测试的私有回归基准，YAML 校验无法进入该路径。
- `utils/particle_progress.py`
  - 新增按物理时间显示的预平衡进度条；正式记录仍按已完成内部子步显示。
- `utils/types.py`
  - 增加计划注入时间、实际准入时间、离开时间和入口等待时间数组。
- `utils/field_io.py`
  - 将上述生命周期数组写入轨迹 NPZ。
- `vis_utils/result_loader.py`
  - 按真实变长帧人口加载数据，并报告真实峰值人口。
  - 为仍存活的永久 ID 保持稳定的临时显示列，避免其他粒子离开后导致轨迹尾迹无故重置；这不会改变或复用仿真状态槽。
- `vis_utils/cli.py`、`vis_utils/plotting.py`
  - 将“目标粒子数”显示改成真实峰值、当前活跃数以及本帧进入/离开数。
- `REPRODUCTION_SPEC.md`
  - 更新当前模型、配置、公式、调用链、变长 NPZ schema、限制条件和常见复现错误。

### 测试

- `test_files/test_particle_continuous_perfusion.py`（新增）
  - 覆盖 Halton 序列、半粒子误差界、双平衡判定、空血管预平衡、正式状态继承、永久 ID、无重叠入口、变长帧、生命周期 NPZ 和稳定显示列。
- `test_files/test_particle_dynamics_config.py`
  - 增加 `passive_tracer` 被拒绝、连续灌注默认参数和配置校验测试。
- `test_files/test_particle_mobility_transport.py`
  - 去掉公开 passive 分发预期，保留原迁移率数值回归测试。

## 3. 核心逻辑变化，按调用链说明

1. `generate_microbubble_trajectories.py` 读取 YAML；生产配置给出 (C_0)、粒径分布、正式帧数和内部时间积分参数。
2. `runner.py` 按原流程生成并验收 CFD 流场；未收敛流场仍不能生成粒子。
3. `particle_advection.py` 预计算速度梯度、黏度扩展、固体壁距离和内法向，然后无条件进入迁移率模型。
4. `particle_inlet_flux.py` 在 root 直段建立固定截面，计算

   \[
   j(s,R)=C_0 b u_n^+(s)p_R(R)I[d_w(s)\ge R],
   \qquad
   \lambda=\int\!\int j(s,R)\,\mathrm ds\,\mathrm dR.
   \]

5. `particle_perfusion_schedule.py` 在完整安全时域内预生成确定事件：

   \[
   t_k=\frac{k+1/2}{\lambda},
   \qquad
   (R_k,s_k)\sim\frac{j(s,R)}{\lambda}
   \]

   其中分布反演由固定二维 Halton 点完成，没有随机数量、随机时刻或随机重采样。
6. `particle_perfusion_transport.py` 从零粒子开始，按内部子步推进；若子步内存在计划事件，则精确拆成“事件前推进—准入判定—事件后推进”。
7. 入口重叠时，事件进入上游等待队列；后续被接受的积分段结束时重试。准入后才成为活动粒子，且状态槽等于永久 ID。
8. `particle_equilibration.py` 在固定保存时间边界检查最近至少一个平均停留时间窗口内的出口通量和人口趋势。两者都满足 `field.flux_tolerance` 后，将当前完整状态定义为正式时间零点。
9. 正式阶段继续同一事件表和同一状态，不重新生成粒子。每帧只写当前真实活动 ID，因此 `frame_offsets` 对应变长切片。
10. `field_io.py` 写入物理场、变长轨迹、永久 ID 注册表和精确生命周期时间；可视化加载器再将活动 ID 映射到稳定的临时显示列。

## 4. 行为变化：改动前 vs 改动后

| 项目 | 改动前 | 改动后 |
| --- | --- | --- |
| 公开粒子模式 | 可在 passive 与 mobility 之间分发 | 仅 mobility；`passive_tracer` 配置直接报错 |
| 域内粒子数 | `n_bubbles` 规定每帧目标数 | 不再规定；约为 \(\lambda\overline T_{res}\)，并自然波动 |
| 初始帧 | 一次性生成同龄粒子 | 从空血管连续灌注，平衡状态继承为正式帧 0 |
| 粒子补充 | 出口离开触发补充 | 出口不触发任何事件；事件只由 \(C_0\) 和入口流量预先决定 |
| 注入时刻 | 与帧/时间步或离场事件相关 | 固定等通量中点时刻，与数值时间步无关 |
| 粒径与位置 | 随机初始/补充采样 | 按联合通量 CDF 和确定性 Halton 序列唯一确定 |
| 入口拥挤 | 可尝试其他位置或最小重叠兜底 | 保持原计划参数在上游等待，绝不强制重叠插入 |
| 状态存储 | 曾有固定人口观察槽概念 | 永久 ID 对应永久状态槽，终止后绝不复用 |
| 帧输出 | 固定宽度记录 | 每帧只保存真实活动粒子，使用 ragged `frame_offsets` |
| 数量调节 | 改 `n_bubbles` | 改物理浓度 (C_0)；(1\times10^6) MB/mL 的注入率约为 (5\times10^5) 条件的两倍 |

## 5. 新增/修改了哪些测试

- 新增确定性 Halton 基数 2/3 精确值测试。
- 新增任意时间内 `|N_in(T)-lambda*T| <= 1/2` 测试。
- 新增“出口通量正确但人口不稳定时不得通过”的双条件平衡测试。
- 新增直通道端到端测试：空域启动、预平衡继承、正式负出生时间、活动人口、无入口重叠、无状态槽复用和进度计数。
- 新增预平衡前已离开粒子的兼容 `birth_frame=-1` 测试。
- 新增四个生命周期时间数组的 NPZ 持久化测试。
- 新增 ragged 帧加载后同一永久 ID 始终使用同一活动期显示列的测试。
- 新增 `passive_tracer` 配置拒绝测试，并保留所有已有迁移率、碰撞、壁面、可视化和 CFD 回归测试。

## 6. 已运行的检查命令及结果

1. 变更文件 Ruff 静态检查：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m ruff check <本次变更的 Python 文件>
   ```

   结果：`All checks passed!`

2. Python 字节码编译：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m compileall -q ulm_microbubble_traj_gen/generate_microbubble_trajectories.py ulm_microbubble_traj_gen/utils ulm_microbubble_traj_gen/vis_utils ulm_microbubble_traj_gen/test_files
   ```

   结果：通过，无语法错误。

3. 新增连续灌注与配置定向测试：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m unittest -v ulm_microbubble_traj_gen.test_files.test_particle_continuous_perfusion ulm_microbubble_traj_gen.test_files.test_particle_dynamics_config
   ```

   结果：12 项通过。

4. 完整回归测试：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -m unittest discover -v ulm_microbubble_traj_gen/test_files
   ```

   结果：88 项通过。PyVista/Trame 测试出现上游 VTK `GetData` 弃用提示，但没有失败。

5. 生产 YAML 加载检查：

   ```powershell
   D:\anaconda3\envs\pmp\python.exe -c "from ulm_microbubble_traj_gen.utils.config import load_config; ..."
   ```

   结果：读取到 `C0=500000 MB/mL`、预平衡上限 `10 s`、Euler、40 个内部子步。

6. 使用已有已收敛真实流场 `results/20260714_190122` 的只读生产参数粒子验证：

   - 注入率：`0.8452573889 bubble/s`；
   - 平均计划间隔：`1.18307159 s`；
   - 预平衡在约 `4.049 s` 达成；
   - 正式 1001 帧真实人口为 `1..2`，平均约 `1.80619`；
   - 入口等待事件：`0`；
   - 轨迹记录：`1808`；
   - 最大内部步位移/网格比：约 `0.0836`，低于 `0.2` 警戒线；
   - 验证未重算 CFD、未写入结果目录。

## 7. 可能的回归风险

1. 当前是二维/2.5D 硬圆盘近似，没有真实三维血管中的面外绕行自由度；高浓度可能产生比真实系统更多的入口等待。
2. 默认当前血管流量下，(5\times10^5) MB/mL 只产生约 `0.845 bubble/s`，因此正式帧看见 1–2 颗微泡是浓度和停留时间的物理结果，不是显示丢失。若提高到 (1\times10^6) MB/mL，注入率线性加倍，但人口仍不是固定值。
3. 几何或流量改变后，10 s 内可能无法满足严格双平衡；此时程序会按设计失败，而不会强行开始记录。
4. 固定截面流量通过插值积分得到。当前真实流场审计中，截面与开放入口参考流量的相对差约 `0.00120`，会写入日志/元数据；它反映截面采样和离散流场差异。
5. 出口时间记录在被接受积分段的末端，等待粒子也在被接受段结束时重试；应继续用 `integration_substeps` 收敛试验评估亚步时间误差。
6. 碰撞搜索仍是确定性 `O(N^2)`；若使用远高于当前要求的浓度，粒子推进时间可能明显增加。
7. 为保护原有数值回归，旧零浓度固定人口内核仍作为私有测试基准存在；它不在生产 YAML 或公开分发中。历史旧 schema 结果仍可被可视化加载器读取。

## 8. 需要人工重点看的 5 个位置

1. `configs/physics_flow_config.yaml`：确认本次实验采用 (5\times10^5) 还是 (1\times10^6) MB/mL，以及 10 s 预平衡安全上限是否符合实验安排。
2. `utils/particle_inlet_flux.py`：重点核对 root 固定截面位置、有效厚度和 (j(s,R)) 的单位换算、有限尺寸可通过判定。
3. `utils/particle_perfusion_schedule.py`：核对等通量中点时刻、Halton 索引起点、粒径边缘 CDF 和条件位置 CDF 的反演。
4. `utils/particle_perfusion_transport.py`：重点看事件拆步、入口等待、双平衡切换到正式 (t=0)、永久 ID 生命周期和正式帧统计。
5. `utils/field_io.py` 与 `vis_utils/result_loader.py`：确认 ragged `frame_offsets`、四个生命周期时间数组和稳定显示列符合下游分析预期。

## 9. 是否有任何超出需求范围的改动

没有修改 CFD 求解方程、血管生成器、Goldman 迁移率系数、近壁剪切载荷、碰撞力公式、连续壁面几何约束、PyVista/VTK CFD 场景或历史 `results` 数据。

为使浓度驱动的变长帧可被现有动画正确使用，本次包含了必要的可视化加载兼容修改；为避免预平衡在 40 子步配置下产生近似二次历史扫描，也包含了等价的监视器性能优化。这两项都直接服务于 Revised v6 连续灌注主链路，不属于无关功能扩展。没有安装新的第三方包。
