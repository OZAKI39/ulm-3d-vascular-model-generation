# Step 3 ParaView 人工检查（PENDING）

自动检查：PASS。Step2 的人工通过不能代替本轮流场验收。请先读 STEP3_REPORT.md，特别是回流与 closure。

## 打开与准备

1. File → Open，打开 /home/lzy/projects/compre_output/step3/20260913_002844/output/flow.pvd，点击 Apply、Reset Camera。
2. 选中 flow，Filters → Threshold，Scalars 选 FluidMask，Lower/Upper 均设 1，Apply。所有物理场只在这个结果上查看。
3. 着色选 VelocityMagnitude_m_s；确认单位 m/s，Rescale to Data Range。第 0 帧应静止。
4. 切到最后一帧（1000 步，2.0366811e-06 s），最大速度约 0.000602064 m/s。比较帧时固定颜色范围，避免自动缩放误导。
5. 从原 flow 建另一条 Threshold，PortLabel 取 1–4，按 PortLabel 着色。1=inlet，2=outlet_01，3=outlet_02，4=outlet_03。叠加 Step2 的 cap_triangles.vtp、mapped_centers_normals.vtp；这些冻结文件只读。端口 marker 的 PhysicalFieldValid=0，物理场零值仅占位。
6. 打开 output/measurement_inlet.vtp 与三个 measurement_outlet_*.vtp 看积分面。需箭头时，对 FluidMask=1 的结果做 Cell Data to Point Data，再 Glyph，Vectors 选 Velocity_m_s，限制箭头数量并缩小 Scale Factor。这些插值仅用于显示。

## A–J 检查

| 项目 | 操作与判断 | 人工结果 |
|---|---|---|
| A 管腔 | 叠加冻结 STL，旋转并 Slice，FluidMask=1 应只在管腔内；允许已记录的原生 0.001 LU inflate 容差 | PENDING |
| B 管壁 | 放大壁与 cap-wall 接缝，无明显穿墙箭头或高速泄漏；近壁中心速度非零本身不是壁 BC 失败 | PENDING |
| C inlet | label 1 的整体箭头进入血管，与 CSV 正 Qin 对照 | PENDING |
| D outlet_01 | label 2 与 CSV signed 流量一致；正表压启动可向内流，记录帧号与方向 | PENDING |
| E outlet_02 | label 3 的较高正表压启动向内推动，应结合压力、密度及 CSV 解读 | PENDING |
| F outlet_03 | label 4 保留负表压，检查实际向外流与 CSV 一致 | PENDING |
| G 分叉连续性 | Slice/Glyph 检查主干、分叉无明显断层；有限传播区不可误当稳态断流 | PENDING |
| H 孤立高速 voxel | 用颜色找到极值，放大邻格与 Slice，检查是否单格无几何原因的异常 | PENDING |
| I checkerboard/blow-up | 依次看 0、10、100、1000 的速度、DensityLU、GaugePressure_Pa，无非物理棋盘格或爆涨 | PENDING |
| J 数量级 | 核对单位 m/s 和时间 s；目标入口均速约 0.00035 m/s，与自动最大速度比较 | PENDING |

记录检查人、日期、ParaView 版本、A–J 各项 PASS/FAIL/UNVERIFIED，及异常帧号、端口、坐标。人工异常不能因自动通过而忽略。明确人工确认后才可改 Step3=PASS。

本检查不能升级为稳态、生理、实验、生产 CFD、最终网格或 RBC-ready 验证。若对高 closure 或回流无法作出可信解释，请记录 UNVERIFIED/FAIL，先审查数值合同，不调整物理目标来消除现象。
