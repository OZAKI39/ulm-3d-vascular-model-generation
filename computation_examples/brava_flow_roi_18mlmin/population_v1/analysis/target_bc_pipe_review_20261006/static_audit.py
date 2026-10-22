from common import *

def main():
 targets={'Parameters.cpp':[(530,595)],'read_files.cpp':[(129,181),(220,239),(499,519),(685,776),(1750,1806)],'baf_ini.cpp':[(275,321),(438,467)],'set_bc.cpp':[(984,1042),(1070,1100),(1158,1237),(1440,1459),(1472,1517),(1650,1690),(1926,2032)],'eq_assem.cpp':[(32,169)],'consts.h':[]}
 excerpts=[]
 for filename,ranges in targets.items():
  p=SOURCE/filename;lines=p.read_text().splitlines()
  for lo,hi in ranges:excerpts.append(dict(path=str(p.relative_to(RETAINED)),sha256=sha(p),lines=[lo,hi],text='\n'.join(f'{i+1}: {lines[i]}' for i in range(lo-1,min(hi,len(lines))))))
 branches={
 'Dir':{'path':['Parameters::read_xml','BoundaryConditionParameters::set_values','read_files_ns::read_eq','read_files_ns::read_bc','baf_ini_ns::bc_ini','set_bc::set_bc_dir','set_bc::set_bc_dir_l'],'condition':'bType_Dir set, weakDir=false; no mixed-order automatic override','pressure_inverse':False},
 'Traction':{'path':['Parameters::read_xml','read_files_ns::read_bc','read_trac_bcff','baf_ini_ns::bc_ini (general early return)','set_bc::set_bc_neu','set_bc::set_bc_trac_l','nn::gnnb (full TRI6)'],'condition':'bType_trac + bType_gen + bType_flat; explicit early weakDir=false return','pressure_inverse':False},
 'Neu':{'path':['Parameters::read_xml','read_files_ns::read_bc','baf_ini_ns::bc_ini','set_bc::set_bc_neu','set_bc::set_bc_neu_l','eq_assem::b_assem_neu_bc','fluid::b_fluid'],'condition':'bType_Neu + bType_bfs; pressure scalar prescribed, not a pressure-field inverse','pressure_inverse':False},
 'Weak':{'path':['read_files_ns::read_bc','set_bc::set_bc_dir_w','set_bc::set_bc_dir_wl','nn::get_nnx(TET10)','nn::get_nnx(TET4)'],'condition':'Dir AND fluid/FSI AND explicit Weakly_applied=true; absent from target XML','pressure_inverse':True}}
 put('TARGET_BC_STATIC_CALL_GRAPH.json',dict(retained_source_sha_manifest=sha(IDENTITY/'EXACT_DEPLOYED_SOURCE_MANIFEST.json'),branches=branches,excerpts=excerpts))
 (HERE/'TARGET_BC_STATIC_DISPATCH_AUDIT_ZH.md').write_text('''# 目标 BC 的原生静态执行路径

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
''')
if __name__=='__main__':main()
