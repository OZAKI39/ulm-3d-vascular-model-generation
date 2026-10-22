# F10L：无流直管与三枚9 mm PZT实验环境

F10L工程和本轮必要数值检查通过，已完成三路新声场、两条真实短轨迹、四组英文PNG/PDF及远程回传。当前是实验装置的可核对演示，仍为 **EXPERIMENTAL_GEOMETRY_PENDING / MATERIAL_ACOUSTICS_UNCALIBRATED / SOURCE_PRESSURE_UNCALIBRATED / PHYSICS_UNCALIBRATED，training_ready=false**，不是已校准的真实实验预测。

## 基线与独立工作区

完整F10K HEAD：850fd2df3225351f4fe9564f39cbe72a38662b30；原工作区干净，九项原归档检查通过。最新公开同步分支sync/mevo-acoustic-f10k-20261013的15c42e0e422ff5d1521e5cc8c12ede2688b3de59是整理导出；与完整F10K的722个源码、配置和脚本文件逐字节一致。没有用整理快照覆盖完整历史。

本地工作区：/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10l_20261013/computation_examples/mevo_acoustic_microbubbles

分支：dev/fast-pzt-f10l-straight-tube-three-pzt-20261013。服务器：/workspace/mevo_acoustic_microbubbles/f10l_straight_tube；每次执行使用.execution/commit不可变源码，复用项目独立.venv，Supervisor管理且有超时。执行快照见results/straight_tube_exp1/f10l/execution/。

## 装置

暂定内径1.6 mm；演示长度20 mm、壁厚0.1 mm、外径1.8 mm。实验参数文件中真实壁厚、长度、材料仍为空。两端默认为开放，越过端口停止并保留有效前缀，绝不反射；也支持有明确有限泡半径的封闭端。局部坐标沿X，整体中心/旋转、端部方向、源中心与法向均可配置。

| PZT | 发声面中心/mm | 朝向 | 直径/mm | 最近管外壁或端面间隙/mm |
|---|---|---|---:|---:|
| Side Y+ | (0, +3.9, 0) | (0, −1, 0) | 9 | 3 |
| Side Y− | (0, −3.9, 0) | (0, +1, 0) | 9 | 3 |
| End X+ | (+13, 0, 0) | (−1, 0, 0) | 9 | 3 |

只有这三枚；没有Z±源。侧面间隙7.8 mm，侧面到端部面最小距离8.5 mm，有限圆片未相交。注射装置尺寸未知，无法断言不存在硬件干涉。

## 声学计算

从新9 mm三圆片重新计算，没有复用原四枚6 mm底稿。均匀水c=1500 m/s、ρ=1000 kg/m³、f=2.5 MHz均为演示值。数值源是透明加性压力源，不是含真实反射、背衬和电压响应的实体PZT。机械管壁没有被冒充为声学材料层；多材料声速/密度/衰减、弹性剪切及标定保留接口。

由装置及BLI支持推导网格[346, 160, 160]，h=0.1 mm；PML在外，实际扩展为[378, 192, 192]。每波长6格，不称网格收敛。CFL=0.25，每周期24点，4个记录周期，正号时间相量Re[P exp(+iωt)]。实际BLI支持未到物理域边界。三次求解以及源准备/保存共373.71 s。各路原始时间、显存和稳态差在acoustic/basis_*/metadata.json；稳态门槛为0.01。F6一致梯度对实际复声压的有限差分相对差为3.096e-10。

基场由每次1000 Pa数值输入归一为Pa/source-Pa。ON每路输入300 Pa、同相；泡处实际峰值为1111.838–1412.282 Pa，RMS为峰值/√2。这些不是驱动电压、真实功率或水听器标定结果。

采样整卡峰值显存为15100 MiB（约14.75 GiB），低于16 GiB登记上限，但高于12.90 GiB预估；不是精确分配上界。两条短轨迹验证了合法性及代码一致性，未进行时间步收敛扫描，不宣称积分误差或局部物理精度已全面收敛。

首次启动仅进行不必要的自动二进制下载，650 s进程上限结束，GPU求解次数为0。恢复时只给已安装k-wave-python 0.6.2添加显式Python-only导入守卫，原文件、改后文件与SHA存于environment/；未改任何数值求解文件，未伪造已下载二进制。正式三路使用Python/CuPy。环境安装增加k-wave-python 0.6.2和CuPy CUDA12 13.6.0等依赖，其版本约束将项目matplotlib由3.11.2调整为3.10.9；没有更换CUDA驱动。完整运行环境见dynamics/ENVIRONMENT.json。

## 完整逐泡对照

保留原M0：一次/二次Bjerknes、F2探索性近场、Stokes过阻尼、局部润滑、确定性有序求和及PBD。新增解析圆柱壁面，不读取旧Y管CFD或三角壁面。原残余重叠0.01R、单轮接触和墙面修正0.1R预算未放宽。F3只增加显式物理参数开关，默认冻结条件仍保留；远壁10步与旧F6 M0坐标逐位一致。

初态为32个固定ID、等半径1.78 µm泡，seed=20261013，云团半径40 µm，中心(50,50,50) µm，初始最小泡心距离7.753 µm。半径/频率匹配原APECSS；实际表范围为1.779–1.781 µm、2.499–2.501 MHz、非零0.1 Pa至40000 Pa。任务文字中的下界不能改写原SI单位。零声压走精确静止分支，其他越界拒绝。

| 指标 | OFF | THREE_PZT_ON |
|---|---:|---:|
| 有效物理时间 | 1 ms | 1 ms |
| 实际步数/步长 | 1000 / 1 µs | 1000 / 1 µs |
| 最大逐泡位移/nm | 0 | 458.044 |
| 质心总位移/nm | 0 | 12.8097 |
| Rg起末/µm | 30.102220 → 30.102220 | 30.102220 → 29.940736 |
| 长短比起末 | 1.166419 → 1.166419 | 1.166419 → 1.167637 |
| 连通分量起末 | [32, 32] | [32, 32] |
| 最大接受重叠/nm | 0 | 0 |
| 最大壁面修正/nm | 0 | 0 |
| 积分/保存墙钟/s | 0.342503 | 0.366437 |
| 含准备和观测的单案例/s | 0.705974 | 0.570525 |

两条均合法且逐泡保留。OFF所有保存位置完全不变，背景速度/梯度严格为零。ON局部力和压力使用同一GPU算子在真实保存帧重算；第一帧CPU独立参考相对差为[1.0752521371204333e-16, 4.62446888426884e-16, 4.433760454490621e-16]。短窗口并不证明聚团、到达目标或形态控制。没有拉长动画、合并成球或改动坐标来制造效果。

## 验证、图件和归档

共21项必要工程检查通过，包括坐标旋转/端部翻转、反向或Z侧源拒绝、未知墙厚、4096位置CPU/GPU严格零流、有限半径、小墙面投影、超0.1R拒绝、开放端停止、封闭端约束、相量符号、原M0等价和错误控制维数拒绝。压力/力和场域越界拒绝继续生效；数值通过不等于近场和局部水动力物理准确。

四组图位于results/straight_tube_exp1/f10l/figures/：

1. Fig_F10L_1_apparatus_geometry：真比例3D及XY投影，红管、灰蓝PZT、发声法向和初态位置。
2. Fig_F10L_2_geometry_cross_sections：XY/XZ/YZ、内外径/厚度/间隙、PML与ROI、真实泡半径约束。
3. Fig_F10L_3_acoustic_basis：真实三路及登记合成场的幅值、相位和F6梯度方向。
4. Fig_F10L_4_no_flow_bubble_dynamics：真实固定ID球体、质心、Rg、长短比、分量、接触和墙面质量。

PNG和PDF均实际打开检查；FIGURE_PROVENANCE.json记录输入/图件/脚本SHA，VISUAL_REVIEW.json记录逐图审阅。未平滑或补造科学数据。已独立核对55个远程文件共19.804 MB的SHA-256。原F10K代码及结果未改；首次失败日志和不可变执行快照保留。

## 下一步与限制

真实管长/内外径、墙厚材料、内外介质、有效发声面/间距基准、端口与注射装置、频率/电压/相位/声压标定、微泡浓度/半径分布仍待实验室提供。先确认这些输入和匹配响应表，再预登记可分辨短距离目标、无声/固定控制/反馈对照及到达误差；随后研究受声波尺度、对称性和近场物理限制的形态指标。三路复声压线性叠加不等于任意三维导航或任意微米形状。

未模拟/校准项、模型接口、复现命令分别见STRAIGHT_TUBE_EXPERIMENT_MODEL.md、STRAIGHT_TUBE_LAB_PARAMETERS.md与下节。没有重跑F1–F10K、旧CFD、频率/压力/数量扫参、AI训练或新接触模型。

## 复现入口

```bash
export PYTHONPATH=src:.
# WSL：预检；无Warp时checks只执行CPU检查
python -m scripts.straight_tube_exp1.preflight
python -m scripts.straight_tube_exp1.checks
# 服务器：独立空结果目录中执行一次，不能覆盖已完成声场/轨迹
python -m scripts.straight_tube_exp1.acoustics
# GPU检查需要上一步的新basis.h5
python -m scripts.straight_tube_exp1.checks
python -m scripts.straight_tube_exp1.dynamics
# WSL归档、绘图、审阅后生成门槛和报告
python -m scripts.straight_tube_exp1.archive
python -m scripts.straight_tube_exp1.figures
python -m scripts.straight_tube_exp1.delivery
python -m scripts.straight_tube_exp1.verify
```

科学文件按独立结果目录只创建一次；源码先提交，再由remote模块创建新标签受限任务。
科研进度文件仅补一个F10L中文段落；备份为microbubble_aggregation.md.bak_F10L_20261013，before/after/diff/SHA与标题结构检查见progress_document/CHECK.json。
