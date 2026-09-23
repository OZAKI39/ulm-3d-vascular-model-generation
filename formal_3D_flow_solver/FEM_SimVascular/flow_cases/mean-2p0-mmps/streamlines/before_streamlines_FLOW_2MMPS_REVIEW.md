# 新入口平均速度 2.0 mm/s：三维血流场审核报告

本轮已从零初始场重新完成三维 FEM 求解，新的入口积分流量、总体守恒、稳态变化量和导出文件均通过检查。已生成黑底静态图、两段旋转 MP4 与实际视频关键帧总览。结果保存在独立 `mean-2p0-mmps` case；旧冻结场和原测试均保留。

## 为什么重算，以及实际改动

旧入口平均速度约 0.352841 mm/s，本轮改为 2.0 mm/s。原入口面积为 7.756795804427715e-12 m²，因此目标流量为 **1.551359160885543e-14 m³/s = 15.513591609 pL/s**。新解由原生求解器计算，未把旧速度场乘倍率作为结果。

本工程直接把**带符号流量**写入 XML。源码 `baf_ini.cpp` 先把周边节点清零，再对入口分布按面积积分归一化；`set_bc.cpp` 使用 `Value × 归一化分布 × 外法向` 得到速度。因此流入值必须是 **−1.551359160885543e-14 m³/s**。验证器将入口积分取负后作为正的流入量；JSON/CSV 同时保留带符号积分。`Flat + Zero_out_perimeter` 的内部节点速度可以高于截面平均值，这不改变规定的平均流量。

配置与旧 XML 逐项比较，唯一变化是 `INLET/Value`。入口仍为 `Dir / Steady / Flat / Impose_flux=true / Zero_out_perimeter=true`；壁面无滑移；三个出口为零牵引 `Neu / Steady / 0`。没有指定出口分流比例，也没有添加阻力或 RCR。密度为 1056 kg/m³，动力黏度为 0.00345312 Pa·s。几何、70,363 个节点及 371,402 个四面体均保持原值。

零牵引约束的是出口总应力，不能理解为把每个出口节点的压力强制设成零；记录中的出口平均压力可以非零。

## 计算在哪里完成

使用远程服务器 `f7c62a262077`，CPU 为 AMD Ryzen 7 7800X3D（8 核、16 线程），物理内存约 61 GiB；容器 CPU 配额为 7.68 核、内存上限约 42.06 GiB。显卡为 NVIDIA RTX 4090（24 GiB）。实际使用已验证的 CPU/GPU 混合求解管线，1 个 MPI 进程、1 个 OpenMP 线程。完整硬件原始记录见 [host.json](reports/host.json)。

- 远程结果目录：`/workspace/flow_mean_2p0_mmps_20260922/flow_cases/mean-2p0-mmps`
- WSL 本地结果目录：`/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps`
- 原生求解耗时：28.99 分钟；退出码：0。
- 求解器 SHA-256：`0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7`。
- 完整命令、PETSc 选项、残差与 GPU 采样：[execution.json](reports/execution.json)；原始日志：[solver.log](run/solver.log)。

## 时间步、稳态与数值稳定性

时间步取 `min(0.5 h10 / Umean, 0.05 Dh² / ν)`。新入口速度对应的对流限制为 6.891799e-05 s，黏性限制为 1.504060e-07 s，因此本轮仍采用 **1.504060091688221e-07 s**。这是重新计算限制后的选择。参考雷诺数约为 0.001918。

每 10 步保存一次场。稳态要求连续 5 个保存区间的体积加权速度相对变化不超过 1e−5、边界流量变化除以目标流量不超过 1e−6，同时入口误差、质量不平衡均不超过 1e−6，速度和压力有限，壁面保持零速。第 **70 步**首次满足完整条件；写入原生停止文件后，在第 **71 步**正常结束并保存完整重启状态。最终物理时间为 1.067882665099e-05 s。

全部接受的线性校正及非线性时间步通过残差检查。预条件器陈旧导致的可恢复尝试次数为 0，原日志与接受历史均保留。最终场相对首次达标场的速度变化为 2.104e-15，流量变化为 3.051e-16。入口启动阶段的守恒误差随推进下降，未将尚未稳定的早期场冻结。具体数据见 [steady_intervals.csv](reports/steady_intervals.csv)。

## 入口与出口积分结果

实际入口带符号积分为 **-1.551359160440232e-14 m³/s**；其绝对值对应的截面平均速度为 **1.999999999426 mm/s**。入口相对目标误差为 **2.870457e-10**。此微小差异在既定 1e−6 相对容差以内。

| 出口 | 实际流量（m³/s） | 实际流量（pL/s） | 占出口总量 |
|---|---:|---:|---:|
| OUTLET_01 | 6.605558022269e-16 | 0.660556 | 4.257917% |
| OUTLET_02 | 1.321835981977e-14 | 13.218360 | 85.205026% |
| OUTLET_03 | 1.634675982406e-15 | 1.634676 | 10.537057% |

出口总流量为 **1.551359160440232e-14 m³/s**。入口与出口总量差除以目标流量为 **0.000000e+00**，通过守恒检查。分流比例来自求解结果和独立面积积分，未作为边界条件指定。可复查 [boundary_flows.csv](reports/boundary_flows.csv) 与 [physics_validation.json](reports/physics_validation.json)。

## 速度、压力及新旧结果比较

新场速度范围为 **0.000000–7.097925 mm/s**；代表性的体积平均速率约为 **0.875943 mm/s**。该平均采用四面体四点积分计算速度向量模值，属于数值积分近似；节点算术均值和体积 RMS 也保留在 JSON 中。壁面最大速度为 0.000e+00 m/s。原始压力范围为 -1.266852–2461.935400 Pa，未人为平移压力。

旧场最大速度为 1.252216 mm/s，新场约为它的 **5.668291 倍**；入口流量增加至 **5.668280 倍**。新场与按流量倍率缩放的旧场之间仍有 6.161507e-06 的体积 L2 相对差异，该数值仅用于比较。此低雷诺数下速度量级近似随流量增加，但正式导出始终使用本次重新求得的场。

## 如何阅读静态图和动画

所有三维图使用黑底、浅色少量文字以及右侧独立颜色条，单位统一为 mm/s，色标固定为 0–7.5 mm/s。这样不同视角和近景可以直接比较，颜色条不会遮挡血管。Figure 01 严格显示真实外表面速度：无滑移壁面为零，入口和出口截面可非零。

总体图、四视角图、分叉近景及视频使用透明的真实外形，加一层位于管腔内部的速度采样曲面。曲面由节点到外表面的有符号距离约 −0.35 μm 的等值面生成，在其上插值 FEM 速度向量后取模；**彩色内部曲面不是壁面速度，也不是截面平均速度或沿视线最大值**。该处理仅服务显示，原始网格和导出流场未改变。采样曲面上的极值可能小于整个三维场的最大值，统一色标仍覆盖整个场。

两段视频均为 18 秒的平稳 360° 相机环绕，物理数据不旋转、不随帧变化。18 秒是展示时长，**不是血流时间演化长度**。完整视图逐帧检查了模型投影边界；近景有意聚焦最大速度节点附近的分叉，远端血管可离开近景画面。Figure 04 的六个画面直接从最终 MP4 解码提取，便于快速人工审核。

验收过程中，画面差异检查曾发现初版视频重复使用缓存画面，尽管相机参数已改变。已加入每帧显式渲染并重新导出；最终两段视频均通过完整解码与旋转检查。初版 MP4、失败日志及修正说明保留在 [render_attempt_01](reports/render_attempt_01/)，并新增静止视频拒绝测试。色标与文字在独立的右侧二维区域绘制，避免透明几何遮挡标注。

- [Figure_00_overview.png](figures/Figure_00_overview.png)
- [Figure_01_surface_speed.png](figures/Figure_01_surface_speed.png)
- [Figure_02_multiview.png](figures/Figure_02_multiview.png)
- [Figure_03_flow_validation.png](figures/Figure_03_flow_validation.png)
- [Figure_04_storyboard.png](figures/Figure_04_storyboard.png)
- [Figure_05_branch_detail.png](figures/Figure_05_branch_detail.png)

- [Animation_01_full_vessel.mp4](animations/Animation_01_full_vessel.mp4)：1920×1080，24 帧/秒，18 秒，432 帧均已解码。
- [Animation_02_branch_detail.mp4](animations/Animation_02_branch_detail.mp4)：1920×1080，24 帧/秒，18 秒，432 帧均已解码。

![总体概览](figures/Figure_00_overview.png)

![关键帧总览](figures/Figure_04_storyboard.png)

## 新冻结场与后续适用范围

新候选冻结场位于 [frozen_flow](frozen_flow/)。`steady_flow_mean_2p0_mmps.vtu` 与最终原生 VTU 逐字节一致，`flow_arrays_si.npz` 保存坐标、四面体、边界标签、速度和压力，单位分别为 m、m/s、Pa；另外保留完整原生重启文件。导出后重新加载，检查了全部点坐标、四面体连接、场值、文件散列和流量积分。NPZ 使用原正向四面体次序，VTU 保留求解器的局部顶点次序，两者单元连接等价。

该结果适合作为后续粒子或微泡模拟的**新背景流场候选**，后续调用应显式选择新路径和清单，不能继续误用旧基线。本轮没有运行新的粒子计算，也没有覆盖任何旧冻结数据；原始 checkout 与隔离副本的 23 个冻结文件均已核对散列。

本轮验证说明在既有刚性壁面、牛顿流体、当前网格和出口模型下求解稳定且守恒。尚未开展本新工况的独立网格加密、时间步减半对照或实验比较，不能把当前 PASS 解释为这些研究已经完成。后续若需要讨论小尺度血细胞效应、壁面变形或不同出口负载，仍需另行建立并验证相应模型。

## 测试、证据与复现

本轮正式检查 **30 项全部通过**，其中包含原冻结文件回归检查及新增边界条件、目标流量、出口守恒、冻结导出、图像检查、视频解码证据检查；错误符号、错误流量、改动出口、空白图和截断视频均有拒绝测试。原有测试未删除。测试记录：[final_tests.xml](reports/final_tests.xml)；显示证据：[artifact_validation.json](reports/artifact_validation.json)；总记录：[FLOW_2MMPS_VALIDATION.json](FLOW_2MMPS_VALIDATION.json)。

新增代码与永久测试：

- [scripts/flow_2mmps/audit_outputs.py](../../scripts/flow_2mmps/audit_outputs.py)
- [scripts/flow_2mmps/finalize.py](../../scripts/flow_2mmps/finalize.py)
- [scripts/flow_2mmps/postprocess_remote.py](../../scripts/flow_2mmps/postprocess_remote.py)
- [scripts/flow_2mmps/prepare.py](../../scripts/flow_2mmps/prepare.py)
- [scripts/flow_2mmps/render.py](../../scripts/flow_2mmps/render.py)
- [scripts/flow_2mmps/solve_remote.py](../../scripts/flow_2mmps/solve_remote.py)
- [scripts/flow_2mmps/validate.py](../../scripts/flow_2mmps/validate.py)
- [tests/flow_2mmps/test_deliverables.py](../../tests/flow_2mmps/test_deliverables.py)
- [tests/flow_2mmps/test_physics_guards.py](../../tests/flow_2mmps/test_physics_guards.py)

复查已有结果时，在项目目录依次运行 `validate.py`、`audit_outputs.py` 和 `pytest tests/flow_2mmps`。如需重新渲染，运行 `render.py`。相关 Python 包版本见 [software_versions.json](reports/software_versions.json)。重算应在新的独立目录运行 `prepare.py` 和 `solve_remote.py`；准备与求解脚本遇到已有 run/log 会拒绝覆盖。远程求解脚本使用已记录的固定求解器安装路径，迁移服务器需显式调整该安装位置并核对散列。

```text
AUTOMATED_CHECKS = PASS
FLOW_SOLVE = PASS
INFLOW_TARGET_MATCH = PASS
FLOW_CONSERVATION = PASS
VISUALIZATION_OUTPUTS = PASS
ANIMATION_DECODE = PASS
STAGE_RESULT = PASS
```
