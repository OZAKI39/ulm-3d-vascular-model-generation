# 原生 velocity space 审计

**确定结论：本次求解的 resolved velocity 是四节点四面体连续 nodal P1；pressure 同阶 P1；VMS 稳定化。原生速度没有先采用高阶再被导出降为 P1。**

## 本次二进制与源码的联系

算例 execution.json 的二进制 SHA-256 是 `0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7`，本轮服务器重算一致。二进制位置为 `/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics`；对应 CMakeCache 的 SV_SOURCE_DIR / CMAKE_HOME_DIRECTORY 指向同级 `external/svMultiPhysics-reuse/Code`。

本地 vendor 12 个所引用源文件与该服务器 source 的 SHA 全部相同，关键文件也与保存的 Stage Q `source_patch.json.after` 相同。见 `data/native_build_source_verification.json`、两份 native_source_hashes JSON。固定源码行与完整 SHA 摘录保存在 [reference/native_source_excerpts.txt](reference/native_source_excerpts.txt)，没有用在线最新版替代本次运行源。

## 判断链

| 证据 | 本次实际内容 | 推论 |
| --- | --- | --- |
| [solver.xml](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps/run/solver.xml) | `<Add_equation type="fluid">`；没有 Use_taylor_hood_type_basis；固定网格、无 mesh scale 参数 | 检查源码缺省路径 |
| 实际 mesh-complete.mesh.vtu / result_071.vtu | 70,363 点，371,402 cells；全部 VTK cell type 10，每单元4节点 | TET4，没有 TET10 的边中点 DOF |
| [read_files.cpp:1742](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/read_files.cpp:1742) | THflag 默认 false，未定义时不改；1764 行 nFs=1 | 没启用 Taylor–Hood；仅 TRI6/TET10 才允许该分支 |
| [fs.cpp:300](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/fs.cpp:300) | 第一函数空间继承 mesh 的 eType/eNoN/N/Nx | 使用 TET4 网格基函数 |
| [fs.cpp:85](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/fs.cpp:85) | lStab 分支把速度和压力两个局部空间都设为 lM.fs[0] | pressure 也为同一 P1 空间 |
| [nn.cpp:87](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/nn.cpp:87) | TET4 映射 Tetra4；174–184 给局部节点排序 | 没有隐含额外节点 |
| [BasisTraits.h:115](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/FE/Basis/BasisTraits.h:115) | complete_lagrange_alias_order 中 Tetra4 返回1、Tetra10 返回2 | 程序分支直接确定 polynomial order=1 |
| [LagrangeBasis.cpp:30](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/FE/Basis/LagrangeBasis.cpp:30) | named linear elements 固定 order=1，Tetra10 属 quadratic aliases | 原生 P1，有明确源证据 |
| [fluid.cpp:496](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/fluid.cpp:496) | nFs==1 时 vmsStab=true | 同阶空间使用 VMS |
| [fluid.cpp:1534](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/fluid.cpp:1534) / [1727](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/fluid.cpp:1727) | divU 为三方向导数和，连续性残差 `Nq(a)*divU - upNx` | 弱连续性，不是逐 tetra divU=0 约束 |

这里“P1/P1”不表示原连续 Navier–Stokes 方程允许压缩，而是其离散近似不提供逐单元强零散度。残差稳定化细尺度 `up` 在积分点计算；没有保存的高阶 nodal field 可直接替换当前速度。

## Checkpoint → 原生 VTU → frozen

[output.cpp:200](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/output.cpp:200) 取 current velocity；296–299 写 header、cplBC、Y_n、A_n。[Array.h:346](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/Array.h:346) 使用 `row + col*nrows_`。

本例 header `<8i3d` 共56 bytes，前8整数 `(1,1,1,70363,0,4,0,71)`；文件长度 4,503,288 = 56 + 64×70,363 bytes。无额外耦合状态/位移块；Y_n 为每节点4个 double（u_x,u_y,u_z,p），随后同尺寸 A_n。按源格式读取的全体 Y_n 与原生 VTU 速度、压力逐值完全相同；不是根据“看着像速度”猜偏移。

[vtk_xml.cpp:923](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/vtk_xml.cpp:923) 使用同一 current solution；1100–1107 的 outGrp_Y 逐节点复制 lY。[VtkData.cpp:166](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/vendor/svMultiPhysics_stage_q/Code/Source/solver/VtkData.cpp:166) 写 vtkDoubleArray，没有投影或平均；240–261 的 vtkPoints 实际输出 float32 坐标。输入 mesh 坐标本来就 float32，`vtk_xml_parser.cpp:524` 读入 double，输入/输出全点坐标本轮证明数值相同。

[validate.py:97](/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/scripts/flow_2mmps/validate.py:97) 直接 copyfile；102 行检查 hash。Particle [migrate_particle9a_flow.py:56](/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/scripts/migrate_particle9a_flow.py:56) 再 copy2。所有节点速度、坐标、单元行均审计，只有 canonical 局部顶点排序变化 `[1,0,2,3]`。

## 哪些是确定的，哪些仍未知

P1 阶次、double 原生速度、checkpoint 转换链、无新增数值误差、native 六截面已有损失，均由配置/源/数值共同确定。VMS 弱连续性能够容许此类局部散度，实际 Gauss 数值证明当前损失来自该已解析场。

本轮没有网格收敛序列，没有求解器当前内部残差的重新装配，也没有独立未量化的 authoritative solver coordinates。因此**不声称**已经证明误差幅度主要由网格粗细、某个 stabilization 参数或坐标初始量化单独造成；这不妨碍定位到 native field 的局部守恒限制。
