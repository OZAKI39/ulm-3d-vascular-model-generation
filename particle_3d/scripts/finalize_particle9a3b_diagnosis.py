"""Write evidence-grounded Chinese diagnosis; never updates production policy."""
from pathlib import Path
import csv
import json

REPO=Path(__file__).resolve().parents[2]
REPORT=REPO/'particle_3d/reports/particle9a3b_flowfield_conservation'
FEM=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular')


def main():
    read=lambda name:json.loads((REPORT/'data'/name).read_text())
    s=read('diagnosis_summary.json');stats=read('divergence_statistics.json');wall=read('wall_flux_audit.json')
    inlet=read('inlet_flux_audit.json');glob=read('global_gauss_balance.json');native=read('native_to_export_audit.json')
    rows=list(csv.DictReader((REPORT/'data/slab_gauss_balance.csv').open()))
    t='\n'.join('| '+ ' | '.join([f"{float(r['section_s_m'])*1e6:.6f}",f"{float(r['Q_section_m3_s']):.12e}",
      f"{float(r['Q_loss_relative_inlet'])*100:.6f}%",f"{float(r['Q_loss_m3_s']):.9e}",
      f"{float(r['volume_integral_divergence_m3_s']):.9e}",f"{float(r['closure_residual_relative_inlet']):.3e}"])+ ' |' for r in rows)
    st='\n'.join('| '+ ' | '.join([label,str(stats[key]['count']),f"{stats[key]['min']:.6f}",f"{stats[key]['max']:.6f}",
          f"{stats[key]['RMS']:.6f}",f"{stats[key]['volume_weighted_mean']:.9g}"])+ ' |'
          for key,label in [('domain','全域'),('root_inlet_to_first_junction','root（中心线归属）'),('wall_adjacent','WALL 面相邻'),('interior','非 WALL 面相邻'),('exact_inlet_to_last_legal_section_slab','入口→54.996557 µm 精确 slab')])
    text=f'''# 一句话结论

**NATIVE_FIELD_LOCAL_CONSERVATION_LIMITATION_FOUND。流量的局部不守恒首次出现于原生 svMultiPhysics 的已解析 P1 速度场；checkpoint、原生 VTU、frozen VTU 和 Particle 使用的是同一个速度场，导出或 Particle 读取没有新增损失。**

6 个真实内部截面的流量损失为 **0.290603%–3.601346%**，由封闭上游体积的负散度积分解释，最大 Gauss 残差 / Q_in 为 **{s['max_gauss_residual_relative']:.3e}**。没有提高 tolerance，没有修正/缩放速度，没有运行微泡。

# 为什么入口守恒但内部截面不守恒

把整个网格想成许多小水管单元。入口与所有出口合起来，流量差几乎为零；但某些小四面体像有数值上的微小水源，另一些像有水汇。全管的水源和水汇可以抵消，截到半路时却未必抵消。因此“总入口等于总出口”不等于“中间每一刀都等于入口”。这里不是物理血液被吸收，而是离散速度场没有逐单元零散度。

本次独立重算全局 inlet/outlet 相对差 **{glob['global_inlet_outlet_relative_error']:.3e}**；全域 ∫div(u)dV = **{glob['volume_integral_divergence_m3_s']:.6e} m³/s**，确实接近零。这与内部百分量级误差并不矛盾。

# 原生 CFD 自己有没有这个问题

**有，准确说是原生离散求解所保存的 resolved nodal velocity 已有这个问题。** 原生为 TET4、连续线性速度 P1 / 线性压力 P1，VMS 稳定化；不是隐藏的高阶速度被导出降阶。

依据实际 solver.xml、与运行二进制构建源目录对应且逐文件 SHA 相同的源码，以及 checkpoint 数值。原生 checkpoint 中 70,363 节点 × 3 速度 DOF 与原生 VTU 的最大差、RMS 差、相对 L2 差均为 **0**；压力也逐值相同。运行二进制 SHA 与原算例 execution.json 一致。完整源码定位见 [NATIVE_VELOCITY_SPACE_AUDIT_ZH.md](NATIVE_VELOCITY_SPACE_AUDIT_ZH.md)。

VMS 连续性残差包含测试函数加权的 div(u) 和稳定化细尺度项。它没有把每个单元的 div(u)=0 当作单独约束。源码中的残差细尺度 `up` 不是一个已保存的高阶输运速度，不能虚构一条“更守恒的 native 高阶曲线”。本轮证明了损失所在层和数学机制；没有凭一张网格把误差幅度进一步归因到某个网格参数、稳定化常数或入口过渡。

# 导出 VTU 有没有增加问题

**没有可测新增误差。** 原生输出直接复制当前解的节点速度到 `vtkDoubleArray`；`scripts/flow_2mmps/validate.py:97` 再用 `shutil.copyfile` 保存 frozen VTU。三个 VTU 文件（native、frozen、Particle reference）的 SHA-256 均为 `129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d`。

没有 interpolation、projection、averaging、速度 float32 转换或网格重采样。坐标在 VTK 中为 float32；求解器输入网格本身也为 float32，读入 double 后再写出，实际数值变化为 0。原网格 NPZ 的 float64 坐标同样只是这些数值的精确提升，不能冒充未量化原坐标。**cannot test authoritative precision loss**：没有另一套未量化且真正用于本次求解的坐标；不能据此评估重新用更高精度网格求解会怎样。当前求解输入→导出→Particle 的坐标变化为 0；六截面 native/Particle 重算的最大通量差为 3.1554e-30 m³/s（Q_in 的 2.034e-16，局部顶点排序引起运算次序差），不是坐标量化新增损失。

见 `data/native_to_export_audit.json`、`data/precision_audit.json`。

# Particle 读取有没有增加问题

**没有。** reference 仍是原 VTU 的逐字节副本；速度保持 float64 bitwise equal；坐标 float32→float64 后 numeric equal，最大差为 0（dtype 改变，所以字节比较不是 equal）。canonical tetra 相对 flow 文件只交换局部前两个顶点 `[1,0,2,3]`，每行四个节点集合与物理单元相同，没有换单元或错配节点。

`FrozenFEMField` 在同一个四面体上用重心权重做 P1 插值。永久测试直接抽查真实 tetra 内点与原 VTU 节点速度加权值相同；独立边方程求导与 Particle 梯度的最大散度差为 1.0914e-11 s⁻¹，属浮点运算顺序差。未发现明确读取、映射或重建软件 bug。见 `data/export_to_particle_audit.json`。

# P1 divergence 是什么情况

每个四面体内部速度为线性函数，所以 div(u) 是常数。本轮解析求解节点差分的线性方程，再取速度梯度的迹；没有 finite difference，没有平滑。

| 区域 | tetra 数 | min (s⁻¹) | max (s⁻¹) | RMS (s⁻¹) | 体积加权均值 (s⁻¹) |
| --- | ---: | ---: | ---: | ---: | ---: |
{st}

root 的描述性统计按 tetra 质心到 authoritative root / daughter 中心线的最近归属划分，并要求 root arclength 在首分叉前；这不是分叉交界面的精确裁切定义。另列入口至最后合法截面的**真实裁切 slab**统计，体积权重使用真实截取体积。WALL-adjacent 指至少具有一个真实 WALL 三角面，其余组可能接触 cap 或仅接触 WALL 顶点，不能读作“离所有边界很远”。两种 root 统计的定义均明确保存在 JSON，Gauss 检查使用精确 slab，完全不依赖中心线归属近似。

全域及分组的 mean、median、P01/P05/P95/P99、普通/体积加权 RMS 均在 `data/divergence_statistics.json`；371,402 个 tetra 的质心、体积、散度、root arclength、wall/root 标签在 `data/tetra_divergence.csv`。

![Root divergence](figures/01_root_divergence.png)

图 01 以约 1 µm 分箱显示分位带，并单独显示体积加权均值，避免散点堆叠。既有正也有负的局部散度；root 整体具有净负贡献，不能说每个位置都负。

# 流量少掉的水去了哪里

**由离散速度场中的净数值水汇解释，不是 WALL 漏水。** 使用外法向、Q_in 定义为正的入流：

`Q_section + Q_wall − Q_in = D = ∫Ω div(u)dV`

所以 `Q_in − Q_section − Q_wall = −D`。不能把正的损失与同号的负散度积分直接写成等号。本轮 Q_wall=0，于是流量损失正好是 `−D`。

真实 slab 是四面体网格与截面上游半空间的交，按共享面连通性保留 inlet 连通分量。完整 tetra 用原体积，截穿 tetra 用凸多面体几何的行列式体积；D 为这些体积乘单元常散度之和。边界分别来自真实 INLET、真实 WALL 和切面。6 个 slab 均完整包含入口、不含其他出口；另外切到的远处血管不会混入。没有圆柱或中心线管体替代。

| s (µm) | Q_section (m³/s) | 损失 / Q_in | Q_in−Q_section (m³/s) | D (m³/s) | 闭合残差 / Q_in |
| --- | ---: | ---: | ---: | ---: | ---: |
{t}

最大 VTK / 独立 NumPy 截面通量差 / Q_in = **{s['max_two_integrators_difference']:.3e}**。Gauss 比较采用 signed flux，另保存 positive flux；这 6 个截面两者在舍入精度内一致（最大差 3.1554e-30 m³/s）。没有把正向积分误当净通量用于 Gauss。

![Representation comparison](figures/02_flux_by_representation.png)

图 02 三条曲线重合。数据 CSV 另保存原生 VTU 一列，它也重合；native 曲线来自真实 checkpoint 数值，不是假定。

![Gauss balance](figures/03_slab_gauss_balance.png)

图 03 使用 `−∫div(u)dV` 与正损失比较，单位 pL/s；WALL 修正单独为零。

# WALL 有没有漏水

**没有。** {wall['node_count']:,} 个 WALL 节点的三个速度分量全部精确等于 0。{wall['triangle_count']:,} 个 WALL 三角面的最大法向速度、最大绝对三角形流量、所有绝对三角形流量之和、总 signed leakage 均为 **0**。真实 slab 中裁切后的 WALL 部分同样为 0；并非只靠正负抵消。

![Wall audit](figures/04_wall_leakage_audit.png)

入口本身也复算一致：面积 {inlet['inlet_area_m2']:.12e} m²，平均法向速度 {inlet['mean_normal_velocity_m_s']*1000:.12f} mm/s，最大节点速度 {inlet['max_speed_m_s']*1000:.9f} mm/s，Q_in={inlet['Q_in_m3_s']:.15e} m³/s。入向单位法向为 {inlet['inward_normal']}。INLET.vtp 只保存 GlobalNodeID，没有另一套独立速度；各层都通过这些 ID 取同一节点速度。

# 0.0001% tolerance 原来是干什么的

`mass_limit=1e-6` 是相对量，即 0.0001%。它用于 **global inlet-vs-target** 以及 **sum(outlets)-inlet** 两个 boundary flux gate，都以 Q_target 归一化。Stage SV1.3 `prepare.py:35` 明确设置，2.0 mm/s 算例从 Stage Q production_policy 继承。`postprocess.py:55-57` 给出两个误差定义，`sv13.py:49-52` 将其用于稳态停止门槛。不是任意内部平面的误差保证。

P9-A.3 按上一轮要求把原策略沿用到 internal-section gate，因此检测到了当前阻断；这属于验收用途扩展，不能据此宣称原生求解未满足其已有全局验收。完整引用和 SHA 见 `data/mass_tolerance_provenance.json`。

# 目前是否合理要求内部截面也达到这个值

可以把它保留为**尚未达成的严格输运目标**，但没有证据说当前原生 P1 场必然能满足它。**保留原全局 1e-6，不改现有内部 gate，本轮不推荐放宽值：NO_TOLERANCE_RECOMMENDATION_YET。**

当前 mesh h10=0.275672 µm；最小检查偏移仅 0.068918 µm，仍损失 0.290603%。节点/导出误差为 0、积分闭合为约 1e-15，说明百分量级差并非积分噪声。仅有一个网格和一种离散方案，没有网格收敛序列，也没有允许多少输运浓度/粒径通量偏差的预算，因此不能从观察到 0.3% 就倒推出 1% 合理。加密候选位置不是网格收敛验证。

# 真正应该修哪一层

需要面对**速度场的局部守恒性质**。简单重复导出、直接读原生 checkpoint 或替换 Particle 的同阶 P1 插值都不会改变结果。若保持已有 CFD，最小研究方向是隔离的、受约束的保守输运场重建原型；它会改变 Particle 用的速度，因此必须作为新表示独立审核，不能混称“无损重新导出”。

数学上连续 P1 已属于 H(div)；只标记 H(div) 或换容器不会自动零散度。下一方案必须同时施加局部零散度/面通量闭合，并检验 WALL 全矢量 no-slip、梯度与原端口 flux。详见 [NEXT_FIX_OPTIONS_ZH.md](NEXT_FIX_OPTIONS_ZH.md)。本轮没有实施任何方案，也不认定必须立刻重跑 CFD；是否需要新 CFD 算例取决于保守重建能否满足这些约束及误差预算。

![Pipeline](figures/05_velocity_representation_pipeline.png)

# 2D→3D interior injection 是否还能继续

**BLOCKED_BY_NATIVE_FIELD_CONSERVATION。** 在保持当前 gate 和当前 field 的前提下，不能继续 joint sampler、2000 events、30 smoke 或正式 500。后续若独立保守场通过审核，才可进入 YES_AFTER_FLOWFIELD_FIX；本轮没有修复，不能提前写成该状态。

# 验证、复现和证据范围

永久测试目录为 `particle_3d/tests/particle9a3b_flowfield_conservation/`，8 项测试验证常速截面、线性零/非零散度 Gauss、Particle 解析梯度、零 WALL 节点、入口复算、原生/导出/Particle 节点与采样一致性、两种截面积分器。**本地 8/8 通过（3.246 s），服务器 8/8 通过（1.547 s），无 skip。**

此外在同一个真实 371,402 tetra 网格、约 30 µm 的实际裁切体积上替换三种解析控制场，复用已求得的真实几何体积。常速、线性 div-free 和 div=30 s⁻¹ 三组的 Gauss 相对残差分别约 6.16e-15、4.10e-15、4.02e-16。控制场可能有非零 WALL 通量，均显式纳入闭合；没有把 div-free 错解成任何切面通量都相同。

服务器目录：`{s['source_repo']}`。实际 **{s['workers']} workers**（6 个不同 PID），OMP/OpenBLAS/MKL/NumExpr 全部为 1；计算 runtime **{s['runtime_s']:.6f} s**，包含散度、6 slabs、各表示对照、真实网格解析控制和数据写出，不含开发、SSH 传输及画图。环境和 PID 在 `data/diagnosis_summary.json`，复现见 [REPRODUCE.md](REPRODUCE.md)。测试实际结果、保护校验和交付时长另见 `data/delivery_verification.json`，不把计算用时误称为整个任务用时。

源码与数据均为新增文件；**26,691 个既有文件的 SHA 全部未变**（包括旧 P9-A.3 全部53个文件），原 HEAD、index 均未变。`git diff --stat` 仍为原有12个文件、194行新增/13行删除，不包含本轮未跟踪新增文件；本轮独立统计在 `logs/new_files_stat.txt`。未 commit/push。根因结论限于当前固定算例，不把它外推成真实血流存在压缩性，也不声称已经有经过验证的保守修复。

# 本轮没有做什么

没有改 CFD physics、边界条件或 solver；没有改 P9-A.1、P6.5、handoff、contact solver、微泡分布或 Method B/C；没有改 tolerance；没有 rescale internal velocity；没有继续 joint sampler；没有生成 2000 bubble events，没有运行 30 smoke，没有生成新的正式 500。旧 `particle9a3_interior_inlet` 保持只读证据。
'''
    (REPORT/'FLOWFIELD_MASS_CONSERVATION_DIAGNOSIS_ZH.md').write_text(text)
    source=FEM/'vendor/svMultiPhysics_stage_q/Code/Source/solver'
    audit=f'''# 原生 velocity space 审计

**确定结论：本次求解的 resolved velocity 是四节点四面体连续 nodal P1；pressure 同阶 P1；VMS 稳定化。原生速度没有先采用高阶再被导出降为 P1。**

## 本次二进制与源码的联系

算例 execution.json 的二进制 SHA-256 是 `0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7`，本轮服务器重算一致。二进制位置为 `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics`；对应 CMakeCache 的 SV_SOURCE_DIR / CMAKE_HOME_DIRECTORY 指向同级 `external/svMultiPhysics-reuse/Code`。

本地 vendor 12 个所引用源文件与该服务器 source 的 SHA 全部相同，关键文件也与保存的 Stage Q `source_patch.json.after` 相同。见 `data/native_build_source_verification.json`、两份 native_source_hashes JSON。固定源码行与完整 SHA 摘录保存在 [reference/native_source_excerpts.txt](reference/native_source_excerpts.txt)，没有用在线最新版替代本次运行源。

## 判断链

| 证据 | 本次实际内容 | 推论 |
| --- | --- | --- |
| [solver.xml]({FEM}/flow_cases/mean-2p0-mmps/run/solver.xml) | `<Add_equation type="fluid">`；没有 Use_taylor_hood_type_basis；固定网格、无 mesh scale 参数 | 检查源码缺省路径 |
| 实际 mesh-complete.mesh.vtu / result_071.vtu | 70,363 点，371,402 cells；全部 VTK cell type 10，每单元4节点 | TET4，没有 TET10 的边中点 DOF |
| [read_files.cpp:1742]({source}/read_files.cpp:1742) | THflag 默认 false，未定义时不改；1764 行 nFs=1 | 没启用 Taylor–Hood；仅 TRI6/TET10 才允许该分支 |
| [fs.cpp:300]({source}/fs.cpp:300) | 第一函数空间继承 mesh 的 eType/eNoN/N/Nx | 使用 TET4 网格基函数 |
| [fs.cpp:85]({source}/fs.cpp:85) | lStab 分支把速度和压力两个局部空间都设为 lM.fs[0] | pressure 也为同一 P1 空间 |
| [nn.cpp:87]({source}/nn.cpp:87) | TET4 映射 Tetra4；174–184 给局部节点排序 | 没有隐含额外节点 |
| [BasisTraits.h:115]({source}/FE/Basis/BasisTraits.h:115) | complete_lagrange_alias_order 中 Tetra4 返回1、Tetra10 返回2 | 程序分支直接确定 polynomial order=1 |
| [LagrangeBasis.cpp:30]({source}/FE/Basis/LagrangeBasis.cpp:30) | named linear elements 固定 order=1，Tetra10 属 quadratic aliases | 原生 P1，有明确源证据 |
| [fluid.cpp:496]({source}/fluid.cpp:496) | nFs==1 时 vmsStab=true | 同阶空间使用 VMS |
| [fluid.cpp:1534]({source}/fluid.cpp:1534) / [1727]({source}/fluid.cpp:1727) | divU 为三方向导数和，连续性残差 `Nq(a)*divU - upNx` | 弱连续性，不是逐 tetra divU=0 约束 |

这里“P1/P1”不表示原连续 Navier–Stokes 方程允许压缩，而是其离散近似不提供逐单元强零散度。残差稳定化细尺度 `up` 在积分点计算；没有保存的高阶 nodal field 可直接替换当前速度。

## Checkpoint → 原生 VTU → frozen

[output.cpp:200]({source}/output.cpp:200) 取 current velocity；296–299 写 header、cplBC、Y_n、A_n。[Array.h:346]({source}/Array.h:346) 使用 `row + col*nrows_`。

本例 header `<8i3d` 共56 bytes，前8整数 `(1,1,1,70363,0,4,0,71)`；文件长度 4,503,288 = 56 + 64×70,363 bytes。无额外耦合状态/位移块；Y_n 为每节点4个 double（u_x,u_y,u_z,p），随后同尺寸 A_n。按源格式读取的全体 Y_n 与原生 VTU 速度、压力逐值完全相同；不是根据“看着像速度”猜偏移。

[vtk_xml.cpp:923]({source}/vtk_xml.cpp:923) 使用同一 current solution；1100–1107 的 outGrp_Y 逐节点复制 lY。[VtkData.cpp:166]({source}/VtkData.cpp:166) 写 vtkDoubleArray，没有投影或平均；240–261 的 vtkPoints 实际输出 float32 坐标。输入 mesh 坐标本来就 float32，`vtk_xml_parser.cpp:524` 读入 double，输入/输出全点坐标本轮证明数值相同。

[validate.py:97]({FEM}/scripts/flow_2mmps/validate.py:97) 直接 copyfile；102 行检查 hash。Particle [migrate_particle9a_flow.py:56]({REPO}/particle_3d/scripts/migrate_particle9a_flow.py:56) 再 copy2。所有节点速度、坐标、单元行均审计，只有 canonical 局部顶点排序变化 `[1,0,2,3]`。

## 哪些是确定的，哪些仍未知

P1 阶次、double 原生速度、checkpoint 转换链、无新增数值误差、native 六截面已有损失，均由配置/源/数值共同确定。VMS 弱连续性能够容许此类局部散度，实际 Gauss 数值证明当前损失来自该已解析场。

本轮没有网格收敛序列，没有求解器当前内部残差的重新装配，也没有独立未量化的 authoritative solver coordinates。因此**不声称**已经证明误差幅度主要由网格粗细、某个 stabilization 参数或坐标初始量化单独造成；这不妨碍定位到 native field 的局部守恒限制。
'''
    (REPORT/'NATIVE_VELOCITY_SPACE_AUDIT_ZH.md').write_text(audit)


if __name__=='__main__':main()
