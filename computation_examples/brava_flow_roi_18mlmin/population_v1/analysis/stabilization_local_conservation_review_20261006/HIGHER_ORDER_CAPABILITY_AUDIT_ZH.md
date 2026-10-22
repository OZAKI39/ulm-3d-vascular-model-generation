# svMultiPhysics higher-order capability 静态审计

分类：**HIGHER_ORDER_FLUID_READY_TO_PROTOTYPE**。含义是已有可复用的 TET10/P1 **流体实现链**，不是 production 准入或本次 GPU P2 实测通过。新 P2 CFD=0。

| 能力 | 实际 production 源码证据 | 结论/限制 |
|---|---|---|
| Tetra10 输入 | vtk_xml_parser.cpp:51/67；read_msh.cpp:657–700 | 读取 VTK_QUADRATIC_TETRA，并处理二次节点交换 |
| 不同速度/压力空间 | read_files.cpp:1741–1776；fs.cpp:321–373 | Use_taylor_hood_type_basis=true + TET10 → nFs=2；速度10节点，压力4顶点 |
| fluid 装配 | fluid.cpp:480–750、1443、1768 | 独立 eNoNw/eNoNq、两轮积分、Hessian、速度/压力块；不仅是一个 Tetra10 类 |
| 边中点压力 | fs.cpp:394–474；Integrator.cpp:140、1088起 | 求解前处理无物理 pressure DOF 的边节点，corrector 由顶点插值 pressure |
| PETSc/GPU | PetscLinearAlgebra.cpp:74起 | 以可变 num_elem_nodes 交给 FSILS 装配，Val 经 TH 处理后传入 PETSc；未见 TET4 硬编码阻断；CUDA/P2 实际运行仍 NOT_RUN |
| VTK 导出 | vtk_xml.cpp:1498–1539；VtkData.cpp:737 | 输出全部 eNoN、全部节点数据，10节点转 VTK quadratic tetra；不能将可视化 pressure 中点值误作独立 pressure DOF |
| restart | output.cpp:250–347；实际 solution arrays | 保存全节点 Yn/An/Dn 等；有基础路径，未验证 P2 restart roundtrip、分区重启或 GPU 稳定性 |
| 任意阶 Lagrange | FE/Basis/BasisFactory.cpp、NodeOrderingConventions.cpp | 有 C0 scalar 任意阶基础；当前 reader 的 TH 仅允许 TRI6/TET10，fluid 数组 MAX_SIZE=27；不宣称任意阶 fluid 即插即用 |
| 当前 particle sampler | particle_3d/audit.py:43、field.py:40–83；原 GPU tetra kernel | 明确只接受 celltype10、四顶点及 P1 常梯度；不能保留/采样 TET10 的边中点速度 DOF，不能静默降阶 |

静态链条可直接复用，但以后原型必须验证速度/压力约束、节点顺序、导出、高阶 G/curl/shear 与 GPU solve；不要手写替代 CFD solver。P2/P1 也不是逐单元强守恒方法，C0 P2 的 G 跨单元仍可能跳变；“更高阶”本身不保证修复1189或822。

本次查询 upstream 默认分支固定为 `d18b32339c95b839dc797e572ba46eac948c1af7`（commit date 2026-09-29T07:03:42Z）。已读取其 read_files.cpp、fs.cpp、fluid.cpp、vtk_xml.cpp 和 Integrator.cpp；同样保留 TH reader、mixed fluid assembly、edge-pressure/corrector。文件 SHA 见 [UPSTREAM_AUDIT_IDENTITY.json](UPSTREAM_AUDIT_IDENTITY.json)。[固定 upstream 代码](https://github.com/SimVascular/svMultiPhysics/tree/d18b32339c95b839dc797e572ba46eac948c1af7/Code/Source/solver)和[官方文档](https://simvascular.github.io/documentation/multi_physics.html)仅补充能力判断；生产算法归属仍由实际 c3f0bb… 源树决定。下载失败的 upstream 文件不冒充已经获得 SHA。

**唯一建议的下一步方向：复用现有 svMultiPhysics Taylor–Hood P2/P1 做有限、独立的解析场及高阶采样验证原型。** 依据是已知 P1 结构限制与现成源码路径，尚非“修复必然有效”的证据。本轮不启动、不生成新 mesh、不改变历史 BC；未来任何新 BraVa 边界设计遵从中心线阻力估算方法。
