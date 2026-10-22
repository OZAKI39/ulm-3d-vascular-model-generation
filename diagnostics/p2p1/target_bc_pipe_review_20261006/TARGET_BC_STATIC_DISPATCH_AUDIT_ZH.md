# 目标 BC 的原生静态执行路径

仅使用当前 production retained exact source；带行号、SHA 的证据见 TARGET_BC_STATIC_CALL_GRAPH.json。原生 parser 默认 Weakly_applied=false（Parameters.cpp:592），与速度/压力阶数无关。

| 目标 | XML | 原生有效路径 | weak |
|---|---|---|---|
| Pipe INLET | Dir / General | set_bc_dir → set_bc_dir_l → all_fun::igbc | 无 |
| Pipe WALL | Dir / Steady / 0 | set_bc_dir → set_bc_dir_l → face.gN逐节点赋值 | 无 |
| Pipe OUTLET_01 | Traction + VTP | set_bc_neu → set_bc_trac_l → nn::gnnb / prescribed vector integration | 无 |
| BraVa INLET | Dir / Steady / impose_flux | bc_ini归一化 → set_bc_dir → set_bc_dir_l | 无 |
| BraVa WALL | Dir / Steady / 0 | set_bc_dir → set_bc_dir_l | 无 |
| BraVa O1/O2/O3 | Neu / Steady | set_bc_neu → set_bc_neu_l → b_assem_neu_bc → b_fluid | 无 |

1. Dir 在3D mixed P2/P1下不会自动转 weak。read_bc:753–756 只读取显式 Weakly_applied；set_bc_dir:998 排除 weak，set_bc_dir_w:1230 仅接受 weak。
2. General 的指定空间值由 read_temp_spat_values 装入 gm，set_bc_dir_l 的 igbc 读取；仍由 strong 路径逐节点施加。
3. Impose_flux 仅置 bType_flx；bc_ini:455 对 profile 的曲面积分归一化。它不修改 weakDir，也不改变边界装配路径。
4. Traction + 文件在 read_bc:150–179 设置 trac/gen/flat，明确 weakDir=false 后提前返回；set_bc_neu 分发至 set_bc_trac_l，不经过 set_bc_dir_wl。
5. 当前三个稳态 Neu 设置 Neu/bfs/std/flat，分发至 set_bc_neu_l，再以给定压力负号构造 hg，进入 b_assem_neu_bc。不是对求解的 P1 pressure 做边界反解；原 Backflow_stabilization_coefficient=0 保持。
6. Zero_out_perimeter 对非general profile在 bc_ini:440 将与其它face共享节点的 s 清零，然后进行所需flux归一化。**General BC 在bc_ini:293提前返回，故 pipe该XML flag不额外处理gm节点值**；边缘值依赖提供的解析数据及wall强约束。不能笼统声称该flag总会处理General数据，更不能说它切换weak。
7. 强Dirichlet按face.gN遍历全部节点，不按P1角点数截断。TRI6的三个边中点因此同样约束；eq.dof=nsd+1时最后pressure槽被排除。P2重放检查该语义。
8. 这三个目标路径不需要P1 pressure物理inverse mapping。Dir直接赋值；traction和Neu以面几何积分规定载荷。弱路径的第二次TET4反解仍是独立库缺陷，本轮不修改未使用的路径。

Native mixed nFs=2还会使用原生Taylor–Hood稳定化开关（部分VMS项关闭），不能把条件性P2 CFD宣称为所有稳定化项完全同原P1。体装配共同映射补丁继承上一轮已验证对象；边界调用不应被标记为执行了volume common-geometry helper。

运行证据另见 TARGET_BC_RUNTIME_DISPATCH.json：必须与此图一致才可建立路径。完整原输入的P1 setup和两种阶数下的局部dispatch是不同证据，均不等于BraVa P2全域计算或新场认证。
