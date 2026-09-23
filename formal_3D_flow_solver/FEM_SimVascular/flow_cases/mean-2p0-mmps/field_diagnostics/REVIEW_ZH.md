# 完整血管流场：残差、压力与壁面剪切应力可视化

本次使用入口平均速度 **2.0 mm/s** 的已完成 FEM 流场，生成残差图、压力场及壁面剪切应力旋转动画，并同步放大原有流线动画。原求解结果、网格和原有流线均保持不变，没有重新运行 CFD 或积分微泡轨迹。

## 可直接使用的结果

- [统一浏览页面](OPEN_RESULTS.html)
- [9 页 PPT，内嵌三段视频](Flow_Diagnostics_Presentation.pptx)：非线性双面板主图、逐步首末配对图、残差诊断总图；压力静态图与视频；WSS 静态图与视频；流线静态图与视频。视频点击播放。
- [残差曲线 PNG](figures/residual_convergence.png)、[矢量 PDF](figures/residual_convergence.pdf)、[可编辑 SVG](figures/residual_convergence.svg)。单幅图另存为 `residual_nonlinear`、`residual_linear_termination`、`residual_gmres_histories` 的 PNG/PDF/SVG。
- [压力场旋转视频](animations/pressure_full_vessel_enlarged.mp4)、[4K 静态图](figures/pressure_overview_4k.png)。
- [壁面剪切应力旋转视频](animations/wss_full_vessel_enlarged.mp4)、[4K 静态图](figures/wss_overview_4k.png)。
- [放大的流线旋转视频](animations/streamlines_full_vessel_enlarged.mp4)、[4K 静态图](figures/streamlines_overview_4k.png)。

## 绘图与放大

残差图采用白色背景、英文标注、细线和简洁坐标轴，PNG 为 320 dpi，同时提供矢量 PDF/SVG。非线性主图拆成“每步初始残差的全程演化”和“本步下降数量级”两个面板，辅助配对图每条竖线仅对应同一步。全程参考量、实际非线性容差以及第 4、23 步的停止原因均已明确标注。低残差平台不解释为已证明的物理精度或浮点下限。详见 [残差专项说明与源码核验](RESIDUAL_REVIEW_ZH.md)。

三维图沿用此前流线图的黑底、turbo 配色、英文标题和右侧独立色标。本版取消随视角变化的自动居中与缩放。**相机位置、焦点、相机朝上方向和缩放比例全部固定**；只对显示对象应用绕固定轴的刚体旋转。旋转轴通过血管几何中心附近，并沿原始几何的主要延伸方向；轴在计算开始时确定一次，显示时保持竖直。三类场使用同一旋转轴、中心、角速度与相机。

轴的 SI 坐标系单位方向为 `[-0.5580255480746352, -0.6281610770700111, 0.5422371703878719]`，固定中心为 `[127.68544759998213, 87.31737167766876, 119.98094592122365]` µm，正交相机半高度固定为 **58.714130 µm**，角速度 **20°/s**。利用整个 360° 旋转的几何包络一次性确定尽可能大的固定尺度，完整血管约占三维视窗高度的 94.5%。没有加粗或修改血管几何。某些视角分支在投影中相互遮挡属于真实三维投影，旋转后可分辨。

上一个自动缩放版本及对应代码已备份到 `../field_diagnostics_revisions/before_fixed_axis/`，本次主目录、网页和 PPT 均指向新的固定轴版本。

每段视频 **1920×1080、24 fps、432 帧、18 秒**。这里的 18 秒是固定轴展示旋转的播放时间，不是瞬态血流演化时间。实际冻结结果为第 71 步，时间步长 1.504060091688e-07 s，结果时刻 1.067882665099e-05 s。

## 残差的含义与读取结果

数据来自 `../run/solver.log`，共 **71 个时间步、167 次线性求解、89,212 行 PETSc 迭代监测记录**。重新读取原日志后，与已保存的求解历史逐条匹配。

非线性主图 A 只画每步第一条残差的全程演化，B 单独画 `log10(本步第一条/本步最后一条)`，以数量级表示步内下降量。辅助图用竖线配对同一步首末记录，避免把各步末值连成一条迭代历史。

全程固定参考为第 1 步第 1 条初始 KSP 范数，`R_ref = 5.714650454390e-03`。新图从更高精度的 KSP 记录重构比值，全部 167 条结果均与原 NS 行定义及打印值一致（允许原 NS 文本的四位有效数字舍入）。最后一条记录约为 **3.105e-14**；其含义是最后一次修正关联的已记录残差，不是最终更新后额外评估的精确误差。

实际非线性停止条件是至少 2 次迭代后，全程参考或本步参考任一个比值达到 **1e-10**。图中的虚线表示全程参考这一项，不代表已确定的浮点下限。第 4、23 步分别以 3、2 次迭代满足条件，较前一步少一次；原始日志确认两处均无预条件器重建，线性求解均正常收敛。低残差不解释为 CFD 解的物理误差，网格与时间步精度仍需独立研究。

完整定义、源码位置、停止规则逐条验证及跳变证据见 [残差专项说明](RESIDUAL_REVIEW_ZH.md)、[审计 JSON](NONLINEAR_RESIDUAL_AUDIT.json)、[逐步 CSV](data/nonlinear_step_audit.csv)。

四面板总图 C 为每次线性求解结束时的 **真实相对残差**，D 为第 1 步第 1 次、第 30 步第 1 次、第 71 步第 2 次求解的实际 GMRES 迭代历史。87 次按相对容差 **1e-10** 终止，80 次按绝对容差 **1e-24** 终止。由于后期右端向量很小，达到绝对容差后，相对残差可以高于 1e-10；这不是求解失败。

所有线性求解的真实最终残差均小于 `max(rtol × 初始残差范数, atol)`；最大比值为 **0.999887610**。日志没有失败线性求解、未识别 NS 行或恢复重试。线性残差属于 PETSc 内部缩放线性系统，不能解释为压力误差（Pa），残差收敛也不等于离散误差或网格收敛已经证明。

原始数据：[逐线性求解 CSV](data/residual_linear_solves.csv)、[全部真实残差迭代 CSV](data/residual_true_monitor.csv)、[逐时间步 CSV](data/residual_time_steps.csv)。

## 压力场

显示完整血管外表面及进出口截面上的原始 P1 节点压力；它是体压力场在边界上的取值。保留现有压力参考，未重新归零，也未把出口零牵引条件误解释成所有出口节点压力严格等于零。

原场压力范围 **-1.266852–2461.935400 Pa**。色标固定为 **−5–2500 Pa**，所有帧一致，未裁去微小负压。压力数据：[SI 坐标 VTP](data/pressure_surface_si.vtp)。

## 壁面剪切应力

原冻结结果只有速度和压力，未直接保存 WSS。本次从原 **P1 四面体速度梯度**、**0.00345312 Pa·s** 的动力黏度和真实 WALL 面法向计算切向黏性牵引的幅值。该定义与本项目 `svMultiPhysics` 的 `bpost` 一致；源码保存的向量采用相反符号约定，幅值不受影响。只处理 **45,221 个 WALL 三角面**，入口和出口截面不作为血管壁。

每个壁面三角面都匹配到唯一相邻四面体；法向朝向流体域外侧。128 个单元还通过独立 4×4 仿射拟合复核，速度梯度最大绝对差为 **3.083e-10 s⁻¹**。WSS 与法向的最大点积为 **1.099e-14 Pa**。

| 量 | 数值（Pa） |
|---|---:|
| 原始三角面最小 WSS | 0.189134 |
| 原始三角面最大 WSS | 51.141688 |
| 壁面面积加权平均 WSS | 8.246773 |
| 节点显示场最大 WSS | 46.515593 |

为避免 P1 单元梯度的片状颜色影响阅读，动画显示的是相邻三角面 **WSS 幅值按面积加权得到的节点场**，随后在面内插值。没有平滑几何，原始面值完整保存在 CSV/VTP；节点显示最大值与原始面最大值不同，论文统计应使用原始面值。固定色标 **0–55 Pa** 覆盖原始面值全范围。

这是已有解的 WSS 后处理，尚未新增网格收敛或独立 WSS 实验验证。当前结果是冻结稳态壁面剪切应力，不是脉动周期平均 WSS，也不包含 OSI。

原始数据：[逐壁面 WSS CSV](data/wall_wss_facets.csv)、[含原始面值、节点显示值及法向的 VTP](data/wall_wss_si.vtp)。

## 检查与复现

本项目 `tests/flow_2mmps` **42 项永久测试通过**。诊断部分包括 8 项 WSS 科学测试、3 项固定轴检查，本次增加 4 项残差语义测试，覆盖归一化不变性、两种参考的 OR 停止规则、最少迭代与上限区别以及第 4、23 步的日志核验。三段固定轴 MP4 沿用已逐帧解码验证的版本，视频哈希未改变。代理已查看新的残差图，科学解释的最终人工审核由用户完成。

PPTX 结构、9 页数量和 3 段内嵌视频的字节哈希通过检查；当前未进行原生 PowerPoint 界面试播。

数据计算、三维绘图、打包代码分别为 `scripts/flow_2mmps/compute_field_diagnostics.py`、`render_field_diagnostics.py`、`finalize_field_diagnostics.py`；仅重绘残差可运行 `scripts/flow_2mmps/residual_figures.py`。在项目根目录使用 `/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python` 运行。

验证记录：[计算检查](COMPUTE_VALIDATION.json)、[媒体检查](MEDIA_VALIDATION.json)、[源文件哈希](SOURCE_LOCK.json)、[永久测试](tests.xml)、[视觉检查](VISUAL_INSPECTION.json)、[总验证](VALIDATION.json)。本次产物的 SHA-256 清单见 `SHA256SUMS.txt`。

定义依据：[本地 svMultiPhysics WSS 后处理源码](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/post.cpp)，[PETSc 真实残差监测](https://petsc.org/release/manualpages/KSP/KSPMonitorTrueResidual/)，[PETSc 线性求解与停止条件](https://petsc.org/release/manual/ksp/)。
