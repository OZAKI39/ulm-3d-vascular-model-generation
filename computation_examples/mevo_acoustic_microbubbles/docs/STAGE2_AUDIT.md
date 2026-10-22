# Stage 2A：现有工程审计

审计日期：2026-10-01。基线提交 `fc3fac0e6fe7f0a57bd90812800714a97fff94a5`，
保留为 `stage1-demo-baseline`；开发分支 `stage2-validated-physics`。
开始修改前，本地和服务器工作区均干净、提交一致。
已逐文件校验 Stage-1 DELIVERY_MANIFEST：18,440 文件、2,214,355,317 字节，零缺失/哈希差异。
现有结果不删除、不覆盖；新计算写入 `results/stage2/`。

## 结构、证据和决策

| 代码 | 当前模型/职责 | 已有证据 | 决策 |
|---|---|---|---|
| control/pzt.py | 独立disc声源的幅值、相位、开关、频率 | 控制向量单测 | 保留API；不是压电FEM |
| acoustics/kwave_backend.py | k-wave-python统一API、CuPy GPU、CW末段录制和复数解调 | 实际CPU/GPU运行；1/2 MHz四源验证 | 保留求解管线；C++ ABI修复单独限时评估 |
| acoustics/basis.py | complex64基底、复梯度、组合、HDF5、有界插值 | 线性复场测试，basis/direct误差9.78e-8/1.56e-7 | 保留数据格式、相位约定、接口 |
| acoustics/convergence.py | 0.1/0.08 mm同坐标比较 | 压力L2差1.40%/3.43%；力差1.23%/3.72% | 保留回归证据；不能称网格无关 |
| bubbles/response.py | 固定f0/阻尼的线性振子 | 单位/符号测试 | 仅unit/smoke/debug；不得作为SonoVue结论 |
| MarmottantBubbleModel | 尚未实现的占位接口 | 明确抛NotImplementedError | 保留兼容性；不冒称已有非线性壳模型 |
| bubbles/forces.py | harmonic primary；截断/正则化secondary | 符号、缩放、互易测试 | primary保留；原pair实现标为LegacyRegularizedBjerknes、demo only |
| bubbles/tracker.py | 准稳态Stokes、RK2/RK4、位移保护子步 | 精确Poiseuille测试、重复迁移、两档步长 | 注入response/validity接口；保留旧算例回归路径 |
| bubbles/contact.py | 位置投影和NNLS接触速度约束 | 30泡接触链、全帧几何检查 | 仅几何约束；不是润滑、黏附或融合 |
| bubbles/population.py | 固定种子、单分散/lognormal/样本 | 可复现、无初始重叠 | 保留；粒径分布和壳参数分布分离 |
| bubbles/metrics.py | 几何连通cluster、最近邻及pair分布 | 保存全部接受子步 | 保留定义；接触直径饱和不能验证轨迹 |
| flow/fields.py | ZeroFlow/Poiseuille/静态VTU点数据插值 | synthetic四面体VTU、域外拒绝 | 保留FlowField.velocity(x,t)，不启动患者算例 |
| coupling/geometry.py | 有限圆柱和球半径几何投影 | mask/投影测试 | 保留demo边界，新增wall model接口 |
| coupling/experiment.py | 控制量、场和粒子耦合 | 实际combined_field保存、source archive | 接入response table；不在轨迹步内积分MHz ODE |
| io/artifacts.py、vtk_export.py | 独立run、源码快照、HDF5/VTP/VTR/PVD | 10个完整粒子run全帧检查 | 保留；Stage-2单独记录第三方版本/validity |
| visualization/plots.py | 原值PNG/PDF，300 dpi | 图像已核验 | 保留；本轮优先文献benchmark |
| tests/、configs/、scripts/ | 16项测试、两频率配置和有限supervisor任务 | 16 passed，ruff passed | 保留回归并增加分层标记 |

## 必须解决的物理和数值问题

1. 壳弹性/黏度、气体参数及半径依赖无SonoVue文献映射；不能延用固定2.5 MHz振子作生理模型。
2. 32泡/20 ms结果虽形成几何cluster，但dt对照末位置RMS差5.72 µm，平均位移差6.13%。
   最近邻约2.5 µm受接触约束支配，不能证明个体轨迹收敛。
3. 原secondary分母使用max(d, R_i+R_j+regularization)，并忽略互相修正的径向响应、迟滞/近场。
   Stage-2严格路径必须记录/中止超适用范围；不得沿用匿名正则化解释近场。
4. 100泡开发探针成本高；旧失败与中断证据保留。性能瓶颈还包括接触投影，不能只计pair力。
5. 迁移dt/2最初被几何子步抵消；补充实际减半步长已有证据。声场细化仍使横向位移变化7.96%。
6. 源幅值不是泡位置局部PNP；response table查询必须使用局部压力、绝对相位一致、严禁静默外插。
7. 数学澄清：对于严格单频的压力梯度，完整周期的非线性V(t)平均力只依赖V的基频投影。
   即便径向运动非线性，`-<V grad p>`与精确基频公式仍应相等；差异可以来自多谐波压力、
   暂态、非整数窗口或数值积分误差，不能为了“展示非线性”人为制造差别。
8. 初期声学run只有开发提交与参数，并非完整未提交源码快照；后期粒子run有源码zip和逐文件hash。
9. C++/CUDA曾缺libhdf5_serial_hl.so.100；安装版0.6.2统一C++接口忽略record_start_index。
   即使修复动态库，benchmark也必须对齐实际采样时刻/窗口，不能直接把全时程当末4周期。

## 分阶段修改与闸门

按用户指定顺序：2A审计封存 → 2B APECSS原版编译/示例 → 2C PROTEUS复用审计 →
2D有出处的SonoVue参数 → 2E非线性径向验证/并行响应表 → 2F力与平移验证 →
2G文献两泡验证 → 2H Lazarus三档压力 → 2I NumPy/CuPy扩展性 → 2J受validity约束的圆管演示。

径向验证失败不进入clustering；两泡验证失败不开展大规模aggregation。
依赖/论文不可获取时保存错误和缺失项，继续不依赖该阻塞的工作，不把未完成项写成通过。
每个新模型先查作者/官方实现，再考虑C/MATLAB包装，最后才是有方程号的最小reference实现。
APECSS上游不直接改写；必要patch另存。生产response采用离线表，CPU多进程且按cgroup限额选worker。

Stage-1原配置和结果作为legacy regression，后续物理默认另设Stage-2配置；保留现有声场缓存身份。
本轮不修改BraVa正式CFD、1565条旧轨迹、原灌注协议、旧路径兼容链接或已停止任务。
