# Manufactured Poiseuille P0/P1基线与几何结论

当前仍无合格particle field winner。以下保留P0/P1准确度、P2原生失败及高保真PLC的最终证据；MMG/Omega_h结论仅按历史索引查阅。

沿用冻结R1=596.779230 µm、R2=626.573561 µm、L=12R、Umean=0.1 m/s、rho=1056 kg/m³、mu=0.00345312 Pa·s及production svMultiPhysics二进制/求解容差。P0/P1/P2是网格策略，不是有限元阶数。旧PLC修复、strict containment及旧P0/P1/P2科学结果冻结，本轮不重算旧34测试或旧1565轨迹。

## 解析几何与边界定义

第一套表面设计采用512周向节点、30轴向段、cap 16径向环；closed/manifold/conforming PLC直接进入官方native TetGen体网格路径，P0/P1不重网格表面。没有换mesher。两档半径使用完全相同归一化规则。

对已安装solver作有限范围源码核查后，未找到明确可用的全速度Dirichlet pressure gauge/nullspace输入。因此采用用户允许的fallback：入口直接评价物理坐标解析速度，wall no-slip，出口 t=(-pI+mu(G+Gᵀ))n，p(z)=(-8muU/R²)(z−L)，n=+ez。此traction同时含横向黏性项与法向压力；不是标量零traction。源码路径、行号和SHA见[能力证据](manufactured_pipe/MANUFACTURED_BC_CAPABILITY.json)，实际输入见各case的BOUNDARY_FIDELITY.json。

连续解析圆管问题有Poiseuille解；三角平壁仍是有界几何近似，不能声称解析圆管解同时也是多面体域的严格解。CFD前独立验收有限面最大径向偏差≤0.05 µm、cap面积误差≤0.05%、速度trace L2≤1%、Q误差≤0.5%、rim零速度。**无nodal flux normalization。** 稳态守恒监视以实际离散入口Q为基准；科学标定Q仍对解析Q比较，不混用二者。

## 固定探针与验收

从正式server SonoVue source实际实例导出 a_min=0.375 µm，near-wall h_lower=max(2 nm,0.001a)=2 nm，所以 **y_min_particle=0.377 µm**。粒径合同SHA与fixed-window合同匹配，未改变粒径分布。hard probes为0.377/1/2/5/10 µm；bulk为0.05/0.1/0.25/0.5R。每组固定3截面×32方位=96点，严格重心坐标包含性复用，不删除、平移或外推域外点。[粒径/间隙证据](manufactured_pipe/PARTICLE_PROBE_CONTRACT.json)。

Particle velocity≤5%，particle gradient/WSS≤8%，内部Q≤5%，dp/dz≤8%，bulk velocity≤5%。0.1 µm保留为WALL_STRESS_DIAGNOSTIC，单独失败不否决生产；div与bulk跨截面变化仅诊断。每组相对L2是整体范数，不代表逐点最大误差保证。

顺序固定：P0失败才P1，P1失败才P2；首个PASS后停止R1加密，原规则只作一次R2。P2如表面被破坏则GEOMETRY_FAIL，不运行CFD。所有后续BraVa/轨迹阶段依赖前置门槛。

## R1_P0：表面与边界

中央区 z/L=0.4–0.6，有限三角形内部最大绝对径向偏差 **0.011256 µm**，门槛0.05 µm。密集射线内缩统计 median/P95/P99/样本max：0.008322 / 0.011138 / 0.011149 / 0.011151 µm。分位数是采样估计；验收使用中央区裁切后的每个finite triangle精确径向极值界，不只检查顶点。

cap面积 INLET/OUTLET = 1.118835837e-06 / 1.118835837e-06 m²；面积相对误差 0.0025% / 0.0025%。wall面积 2.685256627e-05 m²（误差 0.0006%），体积 8.012376048e-09 m³（误差 0.0025%）。

入口离散解析profile相对L2误差 **0.1268%**，Q误差 **0.1345%**，中心节点速度 0.2 m/s，rim为0，缩放因子1。两端cap均检查解析速度trace表示能力；实际出口施加完整Cauchy traction，不施加出口速度Dirichlet。详情：[表面和trace QA](manufactured_pipe/R1_P0/mesh/FIDELITY_QA.json)。

### R1_P0：原始CFD精度

**FAIL_ENGINEERING_CALIBRATION**。失败门槛：particle_velocity, WSS。37,625节点 / 184,712 tetra；CFD 1032.46 s，采样RSS峰值 6.956 GiB。求解收敛和局部精度分别验收。

| 壁距 / 角色 | 域内 / 96 | CFD速度误差 | CFD梯度误差 | 解析节点插值速度误差 | CFD−插值速度范数 / 解析范数 |
|---|---:|---:|---:|---:|---:|
| 0.1 µm / WALL_STRESS_DIAGNOSTIC | 96 | 8.9851% | 6.5153% | 7.9806% | 2.0508% |
| 0.377 µm / PARTICLE_RELEVANT | 96 | 7.1799% | 6.4410% | 6.1615% | 2.0908% |
| 1 µm / PARTICLE_RELEVANT | 96 | 6.7039% | 6.3075% | 5.7349% | 2.0494% |
| 2 µm / PARTICLE_RELEVANT | 96 | 6.5153% | 6.2365% | 5.5816% | 2.0267% |
| 5 µm / PARTICLE_RELEVANT | 96 | 6.2596% | 5.8229% | 5.3624% | 1.9905% |
| 10 µm / PARTICLE_RELEVANT | 96 | 5.8810% | 5.0293% | 4.9901% | 1.9994% |
| 29.84 µm / BULK | 96 | 4.4697% | 2.5527% | 3.6134% | 1.9104% |
| 59.68 µm / BULK | 96 | 2.6797% | 3.6642% | 2.1108% | 1.5554% |
| 149.2 µm / BULK | 96 | 2.1497% | 12.5746% | 3.2534% | 2.2602% |
| 298.4 µm / BULK | 96 | 3.2463% | 12.7426% | 4.0191% | 1.3370% |


WSS相对L2误差 **14.2619%**；三个内部截面Q误差 **2.1965% / 2.0066% / 2.1095%**；dp/dz拟合值 -7872.46 Pa/m，解析值 -7756.65 Pa/m，误差 **1.4930%**。实际入口速度与输入最大差 1.388e-17 m/s。

中央区按tetra体积加权 div(u) RMS = 22.2281 s⁻¹，除以U/R为 0.132653；仅诊断。WSS直接来自相邻tetra原始对称速度梯度，在实际壁面上取切向traction，无节点平均或显示平滑。

插值分解为逐点向量恒等式 u_h−u_exact=(u_h−Iu_exact)+(Iu_exact−u_exact)，**范数不相加**。Iu_exact只用于诊断，不替换CFD；CFD−插值项也不是纯求解器误差，仍包含离散方程、剩余几何和边界影响。不能仅靠一套网格把所有误差严格因果分离。原始证据：[标定JSON](manufactured_pipe/R1_P0/cfd/analysis/CALIBRATION.json)。

WSS补充诊断：同一固定wall targets上的解析节点插值WSS误差 19.3373%；中段全部wall facets面积加权CFD WSS误差 17.8226%。后者只作诊断，不替换原固定targets验收值。P1梯度在单元间不连续，边上最近facet选择有单侧trace含义；不做平均以改变验收。

最终wall-adjacent tetra法向高 min/P5/median/P95/max = 0.2995 / 13.49 / 75.61 / 126.3 / 146 µm。中心轴D/h（h为所在tetra六边均值）min/P5/median/P95 = 4.85 / 4.863 / 5.633 / 6.378。这是局部尺度代理，不是字面数出的跨径单元数；最终单元类型与中间wedge证据分别报告。

## R1_P1：表面与边界

中央区 z/L=0.4–0.6，有限三角形内部最大绝对径向偏差 **0.011256 µm**，门槛0.05 µm。密集射线内缩统计 median/P95/P99/样本max：0.008322 / 0.011138 / 0.011149 / 0.011151 µm。分位数是采样估计；验收使用中央区裁切后的每个finite triangle精确径向极值界，不只检查顶点。

cap面积 INLET/OUTLET = 1.118835837e-06 / 1.118835837e-06 m²；面积相对误差 0.0025% / 0.0025%。wall面积 2.685256627e-05 m²（误差 0.0006%），体积 8.012376048e-09 m³（误差 0.0025%）。

入口离散解析profile相对L2误差 **0.1268%**，Q误差 **0.1345%**，中心节点速度 0.2 m/s，rim为0，缩放因子1。两端cap均检查解析速度trace表示能力；实际出口施加完整Cauchy traction，不施加出口速度Dirichlet。详情：[表面和trace QA](manufactured_pipe/R1_P1/mesh/FIDELITY_QA.json)。

### R1_P1：原始CFD精度

**FAIL_ENGINEERING_CALIBRATION**。失败门槛：particle_velocity, particle_gradient, WSS, bulk_velocity。48,136节点 / 255,770 tetra；CFD 1950.05 s，采样RSS峰值 10.051 GiB。求解收敛和局部精度分别验收。

| 壁距 / 角色 | 域内 / 96 | CFD速度误差 | CFD梯度误差 | 解析节点插值速度误差 | CFD−插值速度范数 / 解析范数 |
|---|---:|---:|---:|---:|---:|
| 0.1 µm / WALL_STRESS_DIAGNOSTIC | 96 | 10.7670% | 8.2986% | 8.1671% | 3.0551% |
| 0.377 µm / PARTICLE_RELEVANT | 96 | 8.9494% | 8.2566% | 6.3393% | 3.1183% |
| 1 µm / PARTICLE_RELEVANT | 96 | 8.5192% | 8.2411% | 5.9039% | 3.1399% |
| 2 µm / PARTICLE_RELEVANT | 96 | 8.3788% | 8.1775% | 5.7479% | 3.1633% |
| 5 µm / PARTICLE_RELEVANT | 96 | 8.1377% | 7.7211% | 5.4987% | 3.1853% |
| 10 µm / PARTICLE_RELEVANT | 96 | 7.7551% | 6.9687% | 5.1274% | 3.1937% |
| 29.84 µm / BULK | 96 | 6.3118% | 4.2144% | 3.7994% | 2.9735% |
| 59.68 µm / BULK | 96 | 4.2156% | 2.4543% | 2.1499% | 2.3485% |
| 149.2 µm / BULK | 96 | 0.9464% | 5.9747% | 1.0504% | 0.4956% |
| 298.4 µm / BULK | 96 | 0.3963% | 6.5263% | 0.6864% | 0.4309% |


WSS相对L2误差 **16.0894%**；三个内部截面Q误差 **0.7489% / 0.7293% / 0.7009%**；dp/dz拟合值 -7777.5 Pa/m，解析值 -7756.65 Pa/m，误差 **0.2688%**。实际入口速度与输入最大差 1.388e-17 m/s。

中央区按tetra体积加权 div(u) RMS = 11.7233 s⁻¹，除以U/R为 0.0699622；仅诊断。WSS直接来自相邻tetra原始对称速度梯度，在实际壁面上取切向traction，无节点平均或显示平滑。

插值分解为逐点向量恒等式 u_h−u_exact=(u_h−Iu_exact)+(Iu_exact−u_exact)，**范数不相加**。Iu_exact只用于诊断，不替换CFD；CFD−插值项也不是纯求解器误差，仍包含离散方程、剩余几何和边界影响。不能仅靠一套网格把所有误差严格因果分离。原始证据：[标定JSON](manufactured_pipe/R1_P1/cfd/analysis/CALIBRATION.json)。

WSS补充诊断：同一固定wall targets上的解析节点插值WSS误差 21.2600%；中段全部wall facets面积加权CFD WSS误差 20.1836%。后者只作诊断，不替换原固定targets验收值。P1梯度在单元间不连续，边上最近facet选择有单侧trace含义；不做平均以改变验收。

最终wall-adjacent tetra法向高 min/P5/median/P95/max = 0.2995 / 10.58 / 81.53 / 127.4 / 146.6 µm。中心轴D/h（h为所在tetra六边均值）min/P5/median/P95 = 11.46 / 12.12 / 13.18 / 15.39。这是局部尺度代理，不是字面数出的跨径单元数；最终单元类型与中间wedge证据分别报告。

## R1_P2：原生表面拓扑失败，未运行CFD

复用与P0/P1完全相同SHA的512边形PLC、h_bulk/R=1/6及请求4层、edge_size_fraction=0.5、decreasing ratio=0.8。仅采用一套有源码/API依据的替代原生表面设置：surface_mesh_flag=True、use_mmg=False，仍由native TetGen生成体网格。此设置切换到原生VMTK表面重网格/投影路径，未改写已冻结adapter或改用其他volume mesher。

依据：[官方2023-05 native源码](https://github.com/SimVascular/SimVascular/blob/f6179e1f07e36cde191cb5e15945ab5eee310113/Code/Source/sv/Mesh/TetGenMeshObject/sv_TetGenMeshObject.cxx)。其中BL调用位于surface flag分支内；默认MMG使用与bulk尺寸绑定的size band。VMTK路径带原surface投影，但这不保证finite facets或拓扑合格。版本、库SHA与行号见[接口核查](manufactured_pipe/NATIVE_LAYER_CAPABILITY.json)，不是对安装binary与tag逐字重建一致性的宣称。

**实际在surface检查处终止：80条free edges、230条non-manifold edges。** Native操作7.30 s后失败；没有volume.vtu，没有本轮wedge intermediate，没有可导出的失败surface，因此不能编造该输出的radial deviation/cap/Q或CFD误差。输入PLC的拓扑/保真PASS不等于remesh输出PASS。证据：生成记录（历史引用 `manufactured_pipe/R1_P2/mesh/GENERATION.json`；结论见本报告和最终摘要，旧路径已清理）、原生日志（历史引用 `monitoring/manufactured_mesh_R1_P2_0001/job.log`；结论见本报告和最终摘要，旧路径已清理）。可选MITK/DICOM的libicu警告也出现在成功P0任务，日志中的直接终止原因是surface拓扑，不把它误诊为求解器或GPU故障。

旧P2曾生成的4层wedge及其厚度证据仍冻结，不能移作本轮P2成果。本轮P0/P1均未满足particle gates，P2无法提供合法待验收volume，无合格R1 winner；按停止规则不追加尝试，不运行R2或BraVa/粒子阶段。

## 本轮五个问题的实际结论

- Q1 / R1_P0：particle velocity相对L2为5.8810%–7.1799%，WSS为14.2619%，均未达到全部门槛。
- Q1 / R1_P1：particle velocity相对L2为7.7551%–8.9494%，WSS为16.0894%，均未达到全部门槛。
- Q2：P1中心轴D/h中位由5.63提高到13.18，中心深处速度更准确；但wall-adjacent tetra法向高中位由75.61升至81.53 µm，近壁速度/WSS未改善。成本从1032.46升至1950.05 s。单纯降低bulk尺寸不能保证需要的近壁分辨率；这两档结果不是严格嵌套网格收敛率。
- Q3：本轮P2在原生表面阶段失败，未获得CFD证据，无法量化其精度收益或CFD成本。
- Q4：没有可认证的最低成本winner；R2未运行。
- Q5：没有新BraVa场或8/46轨迹结果，旧zero-field/slow-fast问题尚未验证改善。

## BraVa尺度与尚未回答的问题

已有89个截面的D/h统计复用原数据：min=5.8383、P5=6.0856、median=7.1573、P95=8.1889。h是所在tetra六边均值，不能声称等同字面跨径单元数。目标P0名义D/h=5，P1/P2名义D/h=12；新BraVa实际值要等合格candidate才能给出。[局部尺度记录](manufactured_pipe/BRAVA_LOCAL_SCALE_QA.json)。

CPU用于native meshing与8-rank CFD；未迁移直管求解至GPU。资源守护记录进程树RSS、OOM及保守nonreclaimable上界；active file cache另列，不因可回收文件缓存误杀任务，不清缓存或改系统参数。粒子GPU阶段仅在门槛通过后启动。

旧1565结果仍为458/527/437出流、141观察截止、2数值失败；不是新固定灌注协议结果。固定协议仍为灌注1.240815660 s+观察10 s、dt=0.5 ms、浓度4e9 m⁻³。现有正式BraVa场/出口压力保持，未为外观修改科学数据。

[总报告](PRODUCTION_WORKFLOW_REVIEW_ZH.md) · [图文入口](OPEN_RESULTS.html) · [当前策略](PRODUCTION_MESH_POLICY.json)。
