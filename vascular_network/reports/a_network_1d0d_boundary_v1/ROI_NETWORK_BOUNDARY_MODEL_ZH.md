# ROI 外部网络模型：目前只完成拓扑审计

结构根 2410 已精确找到，但现有记录只把它作为假定入口；没有找到将它标注为水力源点的数据集证据。外围端点的水力身份也未确认。按本轮源点 STOP 条件，未执行 A 的 operating-point 求解，也未创建或运行新 3D case。

已用原始半径构建全部 9785 条边的 SI 阻力：

`R_e = (8 μ L / π) (r0²+r0·r1+r1²) / (3 r0³ r1³)`，μ=0.00345312 Pa·s。

这是 r(s) 线性变化的精确积分，r0=r1 时连续退化为 `8μL/(πr⁴)`，不使用平均半径。半径正、长度正、没有零长边；没有删环或强制树化。原 SWC 一节点单 parent 的格式无法表达任意多父连接，观测 0 环不证明生理网络没有吻合环。

删除 ROI 内部 108 段真实边后，保留三个精确切点及 O3。包含四个端口的外部组件互不连接，外部 degree 分别为 **1、1、1、0**，见 [external_network_topology.json](data/external_network_topology.json)。盒内未选中血管仍属于外部网络，人工 extension 从未加入真实 A 图。

这说明拓扑上没有三个 outlet 之间的外部连接。要把 O1/O2 化为有限 R，仍需要确定各自远端如何接到参考压力。O3 外部没有边：若将其视为零压 reservoir，就是固定压力约束/R=0 的理想化；若视为未解析的盲端，就是 no-flow，不能输出有限 R3。两种含义不可混淆。

通用 sparse solver 与 Schur 工具已实现并通过解析测试。完整图只建立 scipy CSR matrix，没有 dense N×N。数值解、条件数、节点守恒残差目前均未计算，见 network_solver_audit.json。未把拓扑无连接画成已验证的 conductance heatmap。

合法的 exterior reduction 应保留 ports，内部 sparse 消元：`Y=G_BB−G_BI G_II⁻¹ G_IB`。源点若是非零固定压力，ROI 端口关系一般是 **q_B=Y p_B+f_source**；只有把源点也保留为边界变量，或正确处理齐次参考后，才能写齐次映射。通用测试已覆盖 affine forcing、对称性、被动性、full/reduced 等价和奇异 Z 不可逆情形。

本案例 Y、Z、R1/R2/R3、coupling magnitude 均为 null。CSV 只有 schema，无伪造数字。这不是生理 ground truth；即使后续补齐边界定义，结果也只是 steady Newtonian、geometry-informed 的边界模型。
