# F10R：仿真控制序列

状态：`EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED`，`training_ready=false`。

本轮只表达离线仿真操作；没有真实仪器通信或电压输出。用户暂时搁置AFG1062和功放电脑控制，故不实现SCPI、PyVISA、USB、串口或固件，也不根据假设填入仪器时延和增益。

## 几何与物理范围

沿用F10M-R3/F10Q：PDMS外形32.5×9×9 mm，中心直通道直径1.6 mm、X两端封闭、双顶部1.5 mm孔的轴位于x=±10.85 mm。运动腔体充水、背景流严格为零；无声流、PDMS反射或弹性传播。四枚7 mm圆片在下表位置朝内发声。声场仍是原2.5 MHz均匀水圆片代理，复压力和一致梯度缓存为3/1.5 µm，原1.78 µm泡径半径的APECSS响应不外推。历史32泡4×4×2接触初态是人为合法演示初态，不是本轮声学聚合结果。

| 原始固定顺序 | 圆片中心(mm) | 朝内单位法向 | 未来计划通道 |
|---|---|---|---|
| Y− | (0, −4.5, 0) | (0, 1, 0) | 仪器B / CH1 |
| Y+ | (0, 4.5, 0) | (0, −1, 0) | 仪器B / CH2 |
| X− | (−16.25, 0, 0) | (1, 0, 0) | 仪器A / CH1 |
| X+ | (16.25, 0, 0) | (−1, 0, 0) | 仪器A / CH2 |

上表只是未来映射，不代表硬件已经连接或同步。

## 输入与时序

入口`src/mevo_mb/control/scheduled_pzt.py`提供不可变的`Segment`与`Schedule`。每段包含唯一编号、起始时间、持续时间、固定2.5 MHz频率、四路启停/幅值/相位、模型身份、幅值单位与时间模型。JSON往返保留全部浮点精度；原F10Q控制直接读取INPUT.json。

幅值只允许0–1500 **Pa-equivalent**，上限来自F10Q探索范围。它既不是Vpp，也不是团心实际声压，更不是功放安全限值。没有任何电压换算。实际泡处声压每步仍严格进入原APECSS检查，非法压力拒绝而非裁剪。

- `ONE_HOT_PZT`：每段最多一枚开启，允许OFF。
- `AXIS_PAIR_CONTROL`：每个同轴对共用幅值和启停，保留各对相位差；`axis_pair`辅助函数以负侧相位为零，δ=φplus−φminus。
- `FOUR_CHANNEL_INDEPENDENT`：保留全部独立幅相与启停。

三者共用原Sequential QP、F2、APECSS、GPU场及壁面、阈值和数据格式。没有修改物理模型以获得不同模式的成功。

`IDEAL_SCHEDULED_SEGMENTS`表示理想离线事件。采用公共连续声学时钟`Re[P exp(+iωt)]`，相位偏置在事件处改变，不重置振荡器时间。积分先走到事件；若事件不对齐原时间网格，缩短该段最后一步，在事件时间启用下一段。测试覆盖0.7503 ms非网格切换，实际科研两段切换为1 ms。没有假设四条USB命令可在1 µs内同步生效。

每次事件保存命令时间、模拟生效时间、前后全精度控制、启停、原始粒子位置和逐泡复声压。轨迹另保存实际起始时间/步长/段编号；不会因更换控制重复推进或漏掉时间。

## 本轮冻结的控制与预算

配置`configs/straight_tube_exp1/f10r/protocol.json`、`PREREGISTRATION.json`在任何新轨迹前保存目标、初态、阈值与搜索范围。X=(0.2,0,0) µm、Y=(0,0.2,0) µm，绝对目标容差均0.075 µm；不随结果改动目标。

`REGISTERED_SCHEDULES.json`保留原X、OFF、X+相位−6/+6/+10°完整轨迹计划、两个1 ms的−3/+3°X控制、X−400再X+400 Pa-equivalent的单片序列。其他±1/±3/−10°只做场诊断，不伪造轨迹。原始、绕Z旋转±5°三个固定ID初态均先检查合法性，没有投影。

`FIELD_SELECTED_SCHEDULES.json`保存24个有限X成对场候选筛选；`CONTROL_REGISTRATION.json`记录一次幅值解析细化，`RUN_SCHEDULES.json`在任何新轨迹前冻结结果。共同X幅1036.2016581041792 Pa-equivalent、δX=210°，Y关；初始实际最大声压约1200 Pa。未重新求k-Wave。

四个Y候选先做0.5 ms：单Y+、同轴Y对、原Y控制的X幅值缩至0.25或0.5倍。只有基础质量通过、保持一个几何连通分量、|Rg变化|≤5%、方向正确且有至少1 nm的2 ms线性预估信号才进入完整验证；该预估仅用于筛选，不代替实际轨迹。选形态变化最小者并在完整计算前保存`FULL_SELECTION.json`。

单次最多45 s，远端总预算480 s、最多32条，全部科学数值原则上≤600 s。完整窗口2 ms、dt=1 µs、半步0.5 µs。新X对、新X分段、新低压Y及相位+10°做半步/细缓存；原X和OFF复用已归档F10Q对照，在新接口全轨迹精确重放之后才准许这样比较。未对旋转初态或±6°完整重复数值准入，标为NOT_RUN。

## 验收与证据

原17.8 nm重叠预算、178 nm修正预算、硬球1 nm实际几何门槛、0.01 nm真实互补性、原方程残差和严格零壁修正原样保留。Sequential QP最多五次来自F10Q已归档修订，不是F10R放宽。

基础数值、COM、形状、图拓扑、固定ID、目标达成分别统计。实际目标误差可在单一分辨率下观察，但未经完整对照不称为数值准入；几何连通不是物理融合。全部Y结论标`HOMOGENEOUS_WATER_PROXY_SPECIFIC`，原因和实测代理场诊断见Y报告。

按顺序的复核入口为`prepare`（仅空目录预登记）、`diagnostics`、`review_peaks`、`select`、独立Supervisor的`remote`、只读回传`archive`、`analyze`、`verify`、`figures`。历史结果不覆盖。已有成果直接查阅`analysis/FINAL_GATE.json`、`analysis/COSTS.json`、`local/INDEPENDENT_VERIFICATION.json`和`figures/FIGURE_PROVENANCE.json`。
