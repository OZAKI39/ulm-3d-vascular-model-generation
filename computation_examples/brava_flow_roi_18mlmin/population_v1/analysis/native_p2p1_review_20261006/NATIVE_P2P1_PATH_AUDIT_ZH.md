# 原生 Taylor–Hood 路径审计

2026-10-06。分类：**NATIVE_P2P1_PATH_CONFIRMED**。已实际执行原生基函数、混合函数空间、局部三维 fluid assembly 和压力边节点处理；完整 P2/P1 GPU CFD、VTU/restart 往返运行仍为 NOT_RUN。此分类证明调用链存在且单元级可执行，不是 production 认证。

权威是现有 binary 的 retained exact source，commit `c3f0bb892b765b718f61069ecd9726dbc6d177fd`。重新核对 2431 个源码文件、生产 binary SHA `43cd7d7b9ae9bde91275982225c6bf084859141cdaf2ec1030e1b83c3b94691d` 与指定 PETSc。未下载公共 master 替代，未修改生产源码或重新构建求解器。完整路径、SHA、源码带行号摘录见 [机器审计](NATIVE_P2P1_PATH_AUDIT.json)。

| 环节 | 实际证据 |
|---|---|
| 输入与开关 | `vtk_xml_parser.cpp:51/67` 读取 VTK_QUADRATIC_TETRA；`read_files.cpp:1741–1776` 将 Taylor–Hood 开关映射到 `nFs=2`，三维只允许 TET10 |
| 节点顺序 | `nn.cpp:150–202`：reference corners 为 `(1,0,0),(0,1,0),(0,0,1),(0,0,0)`；六边为 01/12/20/03/13/23；FE permutation 为 `[1,2,3,0,5,9,8,4,6,7]` |
| 定向与所有权 | `read_msh.cpp:657–704` 同步交换角点及中点；`distribute.cpp:209–245` 分配全部 `eNoN` 节点。输入负责保证每条共享边只有一个全局节点，solver 不会把 TET4 自动升阶 |
| 速度与压力 | `fs.cpp:73–197,320–383`：速度10节点/15点积分，压力4节点/4点积分；两空间在同一 reference 坐标上求值 |
| 边界 | `read_msh.cpp:447–520`、`set_bc.cpp:995–1144`：WALL face.gN 的全部节点均施加强 Dirichlet，包含 TRI6 边中点；仅端点在 WALL 上的内部边不因此自动固定 |
| 3D fluid | `eq_assem.cpp` → `fluid::construct_fluid` → `fluid_3d_m` / `fluid_3d_c`，分别使用 eNoNw=10、eNoNq=4 |
| 压力边处理 | `fs::thood_val_rc`、`Integrator.cpp:140` 在求解前使边节点压力行成为单位行；corrector 从四角点恢复用于序列化的边压力 |
| PETSc/GPU | `PetscLinearAlgebra::assemble` → `FsilsLinearAlgebra` → `lhsa_ns::do_assem` → Val/R；`petsc_impl.cpp` 使用通用节点 local-to-global 映射。历史生产 options 为 aijcusparse/cuda；本轮没有重新跑 CUDA benchmark |
| 输出/重启 | `VtkData.cpp:737` 输出 QUADRATIC_TETRA；`vtk_xml.cpp` 保留全10节点；`output.cpp`/`initialize.cpp` 保存/恢复全部节点状态。源码链已读，不声称做了完整高阶 restart round-trip |

单元级 independent unknowns 为 `3×10+4=34`，存储为 `4×10=40` 个槽位，多出的6个压力槽位不是独立高阶压力 DOF。冻结 V2 全域升阶后的对应数字为 814201 个独立未知量、1037320 个存储槽位；这只是表示数据规模，不能作为实际求解内存或耗时测量。

[原生 fixture](NATIVE_UNIT_FIXTURE.json) 链接现有 production object files，仅替换测试 main；没有重写 Taylor–Hood assembly 或执行线性求解。4106 个原生 basis 样本、10/4 mixed spaces、边节点动量贡献、压力—边速度耦合、压力边单位行均通过。局部残差数组4×10，切线数组16×10×10；线性速度场连续性 residual 与解析值最大差2.776e-17。

retained tests 中存在 FE basis 单测和 Stokes manufactured P2/P1 输入。未找到可据以声称已通过三维 TET10 GPU fluid 的现成运行证据，详见 [索引](NATIVE_REGRESSION_INVENTORY.json)。两次初始 harness 错误均已修正并保留监督记录：枚举名；错误使用只服务 reduced pressure spaces 的 quadrature overload。正确调用 `nn::select_ele → fs::init_fs_msh` 后通过。它们不是生产 CFD 失败。

原生开关有两项必须披露的语义：`nFs=2` 自动令 `vmsFlag=false`，连续性 PSPG/tauM 关闭，动量 tauC/tauB 关闭而部分 tauM/up 项仍存在；因此即使参数全部不动，不能把完整离散方程差异说成“只增加速度多项式阶次”。其次，速度采用10节点几何，压力梯度采用4角点几何。后者在 D2 的实际 `construct_fluid` 残差中暴露出曲面映射/积分权重不一致，详见 [主报告](NATIVE_P2P1_DEVELOPMENT_REVIEW_ZH.md)。未绕过该门槛启动 CFD。
