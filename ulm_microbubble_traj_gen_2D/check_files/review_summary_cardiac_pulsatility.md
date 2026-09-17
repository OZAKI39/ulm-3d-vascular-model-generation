# Cardiac Pulsatility Implementation Review Summary

## 1. 本次实现目标

本次实现让 `generate_microbubble_trajectories.py` 当前实际使用的
`utils` 管线支持心跳脉动。接受后的稳态 CFD 场被解释为心动周期
平均参考场；程序在每个粒子内部积分阶段，根据根入口到当前位置的血管
拓扑路径距离、脉搏传播速度和当前物理时间，连续计算局部脉动倍率。

同时修正旧 `utils` 心跳实现中已确认的问题：传播延迟不再错误地乘以
“每周期采样点数”，而是由物理延迟除以采样步长得到；相位采用周期线性
插值，不再使用 `floor + round + clip` 的跳变索引。

实现仍属于“稳态平均场 + 延迟相位调制”的准稳态运动学模型，不是包含
血管顺应性、瞬态压力传播和波反射的完整脉动 CFD。ECG 形状只作为正值
血流倍率代理，不解释为真实测量的流量或压力波形。

## 2. 修改了哪些文件，每个文件为什么改

### 新增文件

- `ulm_microbubble_traj_gen/cardiac_waveform.py`
  - 保存一个紧凑的 ECG 形状单周期代理。
  - 负责周期线性插值、时间导数、周期积分和周期均值归一化。
  - 允许新旧轨迹管线共享一致的波形时间语义。
- `ulm_microbubble_traj_gen/utils/cardiac_pulsatility.py`
  - 根据血管父子拓扑建立根入口路径距离场和物理传播延迟场。
  - 在任意粒子位置、任意内部阶段时间计算局部倍率及空间梯度。
  - 路径距离参数按 vessel ID 向量化映射，避免逐血管重复扫描整张网格。
- `ulm_microbubble_traj_gen/test_files/test_cardiac_pulsatility.py`
  - 覆盖波形、均值、周期、积分、传播延迟、父子血管连续性、旧延迟修复、
    脉动注入反演和 mobility RHS 调制。
- `ulm_microbubble_traj_gen/check_files/review_summary_cardiac_pulsatility.md`
  - 本文件。

### 配置与入口调用链

- `ulm_microbubble_traj_gen/utils/config.py`
  - 新增 `CardiacPulsatilityConfig` 及 BPM、传播速度、初始相位、波形采样数
    等严格校验。
  - 缺少该配置段时默认关闭，保证程序化旧测试配置兼容。
- `ulm_microbubble_traj_gen/configs/physics_flow_config.yaml`
  - 生产配置默认开启心跳；当前为 300 BPM、25000 um/s、2048 点/周期，
    并启用周期平均流量保持。
- `ulm_microbubble_traj_gen/utils/particle_advection.py`
  - 在稳态水动力场准备完成后建立空间心跳模型，并传给唯一公开 mobility
    推进器。
- `ulm_microbubble_traj_gen/utils/runner.py`
  - 将心跳配置接入主调用链，并把关键心跳信息写入 domain metadata。

### 粒子动力学与连续灌注

- `ulm_microbubble_traj_gen/utils/particle_mobility_transport.py`
  - RHS 新增明确的物理时间。
  - 同步调制局部流体速度、完整速度梯度、涡量/旋转输入和轨迹剪切力。
  - 使用 `grad(a*u0) = a*grad(u0) + u0 tensor grad(a)` 保持局部梯度与
    延迟相位场一致。
  - 碰撞力不直接乘心跳倍率，仍通过局部 mobility 转换为附加运动。
- `ulm_microbubble_traj_gen/utils/particle_perfusion_schedule.py`
  - 保留关闭心跳时原有的固定平均注入率路径。
  - 开启心跳时使用 `lambda(t)=lambda_mean*a_inlet(t)`，反演累计数通量
    生成与时间步无关的确定性事件时间。
- `ulm_microbubble_traj_gen/utils/particle_perfusion_transport.py`
  - Euler 使用子段起始时间；Heun 分别使用当前时间和预测结束时间。
  - 保存帧也使用该帧的真实物理时间重新计算 RHS。
  - 打印 BPM、周期、周期数、倍率范围、最大延迟和瞬时注入率范围。
  - 保存逐记录倍率、单周期波形、路径距离、延迟和详细元数据。
- `ulm_microbubble_traj_gen/utils/transport_trajectory.py`
  - 修正旧图网络输运的延迟量纲和整数相位跳变；其他分叉与轨迹逻辑未改。

### 输出、测试和规范

- `ulm_microbubble_traj_gen/utils/types.py`
  - 为轨迹对象增加可选心跳数组，不改变既有必需字段。
- `ulm_microbubble_traj_gen/utils/field_io.py`
  - 心跳开启时把新增数组写入轨迹 NPZ；关闭时不写这些可选键。
- `ulm_microbubble_traj_gen/test_files/test_particle_dynamics_config.py`
  - 增加心跳默认值、显式加载和非法配置测试。
- `ulm_microbubble_traj_gen/test_files/test_particle_continuous_perfusion.py`
  - 增加从空血管开始的完整脉动灌注、Heun 子步和 NPZ 落盘测试。
- `ulm_microbubble_traj_gen/REPRODUCTION_SPEC.md`
  - 更新当前作用域、配置、公式、调用链、输出字段、假设、常见错误和明确
    未实现的完整瞬态 CFD 范围。

`generate_microbubble_trajectories.py` 本身没有加入物理代码；它原本就通过
`utils.runner` 启动主流程。本次在下游模块化接入，避免把配置、波形、
拓扑传播、积分和 I/O 堆到命令行入口文件。

## 3. 核心逻辑变化，按调用链说明

1. `generate_microbubble_trajectories.py` 加载 `physics_flow_config.yaml`。
2. `config.py` 生成并验证心跳配置。
3. `runner.py` 仍先求解、验收稳态 CFD 场；该场现在被解释为周期平均参考场。
4. `particle_advection.py` 建立稳态速度梯度等水动力场，同时建立心跳模型。
5. `cardiac_pulsatility.py`：
   - 递归计算每条血管近端相对根入口的累计 X-Z 路径长度；
   - 加上网格点在当前直血管段上的投影距离；
   - 除以传播速度得到秒单位延迟；
   - 对延迟求空间梯度，并将场延拓到近壁插值需要的相邻网格。
6. `cardiac_waveform.py` 生成一个周期的正值代理并在默认配置下归一化为
   周期均值 1。任意正负时间都通过周期坐标连续求值。
7. `particle_perfusion_schedule.py` 对入口累计数通量进行确定性二分反演，
   得到脉动事件时刻；Halton 半径/位置、永久 ID 顺序不变。
8. `particle_perfusion_transport.py` 从空血管的正式时间零开始，在每个事件
   分段和内部子步向 RHS 传递精确时间。
9. `particle_mobility_transport.py` 在粒子位置插值延迟，构造局部脉动速度、
   梯度和剪切，再进入原有自由流/近壁 mobility、碰撞和几何约束。
10. `field_io.py` 保存轨迹及心跳重建数据，`runner.py` 保存可读元数据。

## 4. 行为变化：改动前 vs 改动后

| 项目 | 改动前 | 改动后 |
| --- | --- | --- |
| 主 `utils` 背景流 | 所有粒子时刻使用同一稳态场 | 稳态场作为周期平均参考，在每个内部阶段按局部延迟相位调制 |
| 传播延迟 | 主流程没有；旧 `utils` 的换算量纲错误 | `tau=根路径距离/传播速度`，直接保留为秒 |
| 波形取值 | 旧流程用整数索引、取整和裁剪 | 单周期、严格周期化、连续线性插值 |
| 周期平均流量 | 旧代理不保证均值为 1 | 生产默认归一化为周期均值 1 |
| 空间位置 | 旧流程使用微泡累计行程 | 主流程使用根入口到 Eulerian 空间位置的血管拓扑路径 |
| 速度梯度/旋转 | 稳态 | 与局部倍率及其空间梯度同步变化 |
| 轨迹 WSS | 静态代理采样 | 静态平均代理乘当前局部脉动倍率 |
| 注入事件 | 固定平均间隔 | 反演积分脉动数通量，高流量阶段更密、低流量阶段更疏 |
| Euler/Heun 时间 | RHS 无显式心跳时间 | Euler 用当前时间；Heun 预测 RHS 用结束时间 |
| 初始状态 | 空血管 | 仍为空血管，不恢复预平衡 |
| ID/槽位 | 永久 ID，不复用槽位 | 不变 |
| 心跳关闭 | 稳态连续灌注 | 倍率恒等于 1，保留原固定速率计划和原 perfusion ID |
| 输出 | 无主流程心跳重建数据 | 新增逐记录倍率、单周期波形、路径距离和延迟场 |
| 静态 CFD HTML | 稳态场 | 仍显示周期平均参考场，不伪装成瞬时心动相位 |

## 5. 新增/修改了哪些测试

新增 `test_cardiac_pulsatility.py`，覆盖：

- 300 BPM 对应 0.2 s 周期；
- 波形为正、周期连续、周期均值为 1；
- 周期积分与正负时间取模；
- 脉动累计注入通量的事件时间反演；
- 旧 `utils` 中 1000 um / 25000 um/s = 0.04 s = 40 个 1 ms 样本，
  不再得到旧错误的 8 个样本；
- 父血管到子血管的根路径距离连续；
- 1000 um 路径对应精确 0.04 s 延迟；
- mobility RHS 中流体速度、粒子速度、WSS 和记录倍率同步变化。

修改配置测试，覆盖默认关闭兼容、显式开启和非法 BPM/传播速度/相位/采样数。

修改连续灌注测试，覆盖心跳开启、空血管起点、Heun 内部阶段、永久 ID、
新 perfusion ID，以及新增 NPZ 键完整落盘。

原有连续灌注、mobility、碰撞、壁面、可视化、流场和几何测试全部继续运行。

## 6. 已运行的检查命令及结果

### 编译检查

```powershell
D:\anaconda3\envs\pmp\python.exe -m compileall -q ulm_microbubble_traj_gen
```

结果：通过，无语法错误。

### 专项测试

```powershell
D:\anaconda3\envs\pmp\python.exe -m unittest ulm_microbubble_traj_gen.test_files.test_cardiac_pulsatility ulm_microbubble_traj_gen.test_files.test_particle_continuous_perfusion ulm_microbubble_traj_gen.test_files.test_particle_dynamics_config
```

结果：21 项全部通过。

### 全量回归

```powershell
D:\anaconda3\envs\pmp\python.exe -m unittest discover -s ulm_microbubble_traj_gen\test_files
```

结果：97 项全部通过，用时约 7.4 s。PyVista/trame 仅打印已有的 VTK 9.6
弃用警告，没有测试失败。

环境中未安装 `pytest`，所以第一次 `python -m pytest ...` 返回
`No module named pytest`。没有安装新包；测试文件基于标准库 `unittest`，
后续全部验证使用 `unittest` 完成。

### 主入口 quick-test

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\generate_microbubble_trajectories.py --quick-test
```

结果：

- 流场求解与物理验收通过；
- 心跳配置正确加载：300 BPM、0.2 s/周期；
- 脉动灌注计划、Numba mobility 推进、NPZ 和 domain metadata 保存成功；
- 输出目录：`ulm_microbubble_traj_gen/results/20260715_005607`；
- 新增数组实测形状：
  - `record_cardiac_multiplier`: `(42,)`
  - `cardiac_waveform_time_s`: `(2048,)`
  - `cardiac_waveform_multiplier`: `(2048,)`
  - `cardiac_path_distance_um`: `(371, 319)`
  - `cardiac_delay_s`: `(371, 319)`
- 波形数组算术均值为 1.0；最大传播延迟约 0.14568 s。

该命令最终在原有 VTK 粗网格完整流线验收处退出：quick-test 的 8 um 网格
不能为全部出口构建根到出口完整流线。诊断文件已正常写入上述目录。此失败
发生在心跳轨迹和结果文件保存之后，不属于本次心跳实现；全量 PyVista 单元
测试仍通过。quick-test 同时打印内部位移比警告，这是原有粗网格/2 子步测试
配置的分辨率警告，不是生产配置 40 子步的结果。

## 7. 可能的回归风险

1. **模型边界**：空间延迟调制是准稳态运动学近似，不是严格的顺应性瞬态
   不可压缩 CFD。不能把新增结果解释为压力波或完整 Womersley 解。
2. **代理波形**：默认曲线来自 ECG 形状，不是主体特异的实测流量曲线；
   当前归一化后范围约为 0.59 到 1.78。
3. **分叉相位场**：路径距离按 raster vessel ID 归属；极粗网格或复杂交汇核心
   可能产生较陡的局部相位梯度，应在生产分辨率下检查。
4. **峰值速度**：心跳峰值会提高瞬时位移。生产运行仍需关注内部步长/网格/
   半径警告并做子步收敛检查。
5. **部分周期人口**：仿真时长不是整数心动周期时，初始相位会影响总注入事件
   数和每帧人口，这是物理计划的预期变化。
6. **文件大小**：轨迹 NPZ 新增两个网格场，最大网格下会增加压缩文件体积。
7. **可视化语义**：现有静态 CFD HTML 仍是平均场；若误认为某个瞬时相位，
   会造成解释错误。
8. **多根入口**：多个根入口目前共享同一入口相位零点；若以后需要不同根入口
   的测量相位，需扩展为逐入口波形。

## 8. 需要人工重点看的 5 个位置

1. `ulm_microbubble_traj_gen/cardiac_waveform.py`
   - 单周期生成、均值归一化、周期插值和积分原函数是否符合希望的波形语义。
2. `ulm_microbubble_traj_gen/utils/cardiac_pulsatility.py`
   - 父子拓扑路径、网格 vessel ID 归属、延迟场和延迟梯度是否符合血管模型。
3. `ulm_microbubble_traj_gen/utils/particle_mobility_transport.py` 的
   `_evaluate_rhs`
   - 动态速度梯度公式、WSS 倍率以及碰撞力不直接缩放的边界。
4. `ulm_microbubble_traj_gen/utils/particle_perfusion_schedule.py`
   - 累计数通量与事件时间二分反演，特别是半事件计数规则。
5. `ulm_microbubble_traj_gen/utils/particle_perfusion_transport.py`
   - Euler/Heun 时间传递、空血管起点、逐记录倍率和 metadata 输出。

## 9. 是否有任何超出需求范围的改动

没有无关功能改动。

为避免主管线修好而旧 `utils` 继续保留已知量纲错误，额外对
`utils/transport_trajectory.py` 做了同范围的窄修复；未修改其分叉、Poiseuille、
帧组装或其他轨迹逻辑。同步更新 `REPRODUCTION_SPEC.md`、新增测试和生成本
review summary 都属于实现可复现性与交付检查范围。

没有修改现有可视化样式、PhiFlow 稳态求解器、Goldman mobility 系数、碰撞
公式、微泡直径分布、空血管初始条件、永久 ID 规则或数组槽位复用策略；也没有
安装或升级任何环境功能包。
