#!/usr/bin/env python3
"""Create the Chinese review report from retained evidence, then print STOP status."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import sha256, timestamp, write_json
from fem3d.mesh_input import read_contract


def main():
    reports=ROOT/"reports/stage01"
    def read(path):
        return json.loads((ROOT/path).read_text())
    contract,lock=read_contract(ROOT)
    adapter=read("reports/stage01/source_adapter.json")
    reference=read("reports/stage01/reference_integrity.json")
    preservation=read("reports/stage01/stage00_preservation.json")
    comparison=read("reports/stage01/mesh_comparison.json")
    visuals=read("reports/stage01/visualization_manifest.json")
    tests=ET.parse(reports/"pytest_results.xml").getroot().find("testsuite").attrib
    failed=int(tests["failures"])+int(tests["errors"])
    skipped=int(tests["skipped"])
    passed=int(tests["tests"])-failed-skipped
    qc={p:read(f"outputs/stage01/{p}/qc/geometry_qc.json") for p in ("coarse","medium")}
    meshmeta={p:read(f"outputs/stage01/{p}/metadata/meshing.json") for p in qc}
    reloads={(p,r):read(f"outputs/stage01/{p}/qc/reload_r{r}.json") for p in qc for r in (1,2)}
    gates={"source_and_contract_sha256":True,"reference_projects_unchanged":reference["status"]=="PASS",
        "stage00_artifacts_preserved":preservation["status"]=="PASS",
        "one_SI_adapter_conversion":adapter["conversion_factor"]==1e-6 and adapter["solver_units"]=="m",
        "exact_welding":adapter["welding"]["maximum_merge_distance_m"]==0 and adapter["welding"]["welding_tolerance_m"]==0,
        "all_profile_geometry_hard_gates":all(v["hard_gate_status"]=="PASS" for v in qc.values()),
        "all_four_fresh_process_reloads":all(v["status"]=="PASS" for v in reloads.values()),
        "five_review_images":len(visuals["images"])==5 and all(sha256(reports/name)==entry["sha256"] for name,entry in visuals["images"].items()),
        "pytest_no_failures":failed==0}
    status="CONDITIONAL PASS" if all(gates.values()) else "FAIL"
    final={"timestamp":timestamp(),"status":status,"hard_gates":gates,"human_review":"PENDING",
        "conditions":["User reviews boundary tags, four caps, tetra cutaway and worst elements", "User accepts the documented nondegenerate low-quality cells for subsequent development, or defines a further authorized mesh study"],
        "tests":{"passed":passed,"failed":failed,"skipped":skipped},"stage2_started":False,
        "reference_integrity":reference["per_root"],"stage00_preservation":preservation["status"]}
    write_json(reports/"stage01_status.json",final)
    def table(headers,rows):
        return "\n".join(["| "+" | ".join(headers)+" |","|"+"|".join(["---"]*len(headers))+"|"]+["| "+" | ".join(map(str,row))+" |" for row in rows])
    def fmt(v):
        return f"{v:.8e}"
    weld=adapter["welding"]
    profile_rows=[]
    for p,row in comparison["profiles"].items():
        quality=row["quality"]["gmsh_min_sicn"]
        profile_rows.append([row["mesh_profile"],f'{row["vertex_count"]:,}',f'{row["tetrahedron_count"]:,}',f'{quality["minimum"]:.6f}',f'{quality["P5"]:.6f}',f'{quality["median"]:.6f}',f'{row["wall_time_s"]:.3f} s',f'{row["peak_rss_kib"]/1024:.1f} MiB'])
    port_rows=[]
    for p in qc:
        for name,port in qc[p]["ports"].items():
            port_rows.append([p,name,port["entity_id"],fmt(port["source_area_m2"]),fmt(port["area_m2"]),f'{port["relative_error"]:.3e}',f'{port["normal_dot_product"]:.15f}'])
    port_detail=[]
    for name in qc["coarse"]["ports"]:
        vals=[qc[p]["ports"][name] for p in qc]
        port_detail.append([name,fmt(max(v["absolute_error_m2"] for v in vals)),fmt(max(v["centroid_displacement_m"] for v in vals)),fmt(vals[0]["source_max_plane_deviation_m"]),fmt(max(v["max_plane_deviation_m"] for v in vals))])
    quality_rows=[]
    for p in qc:
        for key,label in (("gmsh_min_sicn","形状质量（高为好）"),("max_to_min_edge_ratio","最长/最短边（低为好）"),("volume_m3","体积 m³")):
            quality_rows.append([p,label]+[fmt(qc[p]["quality"][key][q]) for q in ("minimum","P1","P5","median","P95","maximum")])
    low_rows=[]
    worst_text=[]
    for p in qc:
        quality=qc[p]["quality"]
        locations=quality["low_quality_nearest_boundary_counts"]
        low_rows.append([p,quality["advisory_count"],f'{quality["advisory_fraction"]*100:.4f}%',locations["INLET"],locations["OUTLET_01"],locations["OUTLET_02"],locations["OUTLET_03"],locations["WALL"]])
        worst=quality["worst_elements"][0]
        coords=", ".join(f"{v*1e6:.6f}" for v in worst["centroid_m"])
        worst_text.append(f"- {p} 最差单元：数组编号 {worst['cell_index']}，Gmsh element {worst['gmsh_element_id']}；中心 ({coords}) μm；最近边界为 {worst['nearest_boundary_patch']}，到该边界三角形中心约 {worst['distance_to_nearest_boundary_triangle_center_m']*1e6:.3f} μm。")
    reload_rows=[]
    for (p,r),row in reloads.items():
        reload_rows.append([p,r,row["cell_count"],row["topology"]["exterior_facet_count"],row["status"],"逐 tetra / 标签一致"])
    text=f"""# Stage 1 — 3D FEM 体网格与边界标签

## 这一步想解决什么？

Stage 0 已确认并冻结输入表面，用户已确认 Stage 0 技术与人工审核通过。Stage 1 只检查能否在这层表面内部建立可信的四面体体网格，并保留 wall、inlet 和三个 outlet 的含义。两个开发网格均完成几何检查和 DOLFINx 保存、关闭、新进程重载。此阶段没有创建速度/压力空间或求解血流；Stage 2 未开始。

当前结论：**STAGE 1 STATUS: {status}**。硬性几何与标签检查通过，但低质量非退化单元及五张审核图仍需用户人工判断。

## 输入是什么？

唯一权威契约：[Stage 0 source_contract.json](../stage00/source_contract.json)。SHA256：`{lock['source_contract_sha256']}`。

冻结 VTP：`{contract['geometry_path']}`。SHA256：`{contract['geometry_sha256']}`。输入单位 **μm**，只在 VTP → SI adapter 中执行一次 `x_m = x_um × 1e-6`。Gmsh、QC、MSH 和 XDMF/HDF5 坐标均为 **m**，面积为 m²，体积为 m³；图像只把显示副本改为 μm。

配套 `_m.stl` 的 SHA256 为 `{contract['meter_geometry_sha256']}`，只用于源文件一致性核对。STL 不保留 `CellEntityIds`，不能作为唯一输入；端口角色来自已冻结契约，未按位置或轴方向猜测。

{table(['entity / facet tag','含义','三角形数'],[[1,'WALL',67071],[4,'INLET',56],[3,'OUTLET_01',44],[5,'OUTLET_02',42],[2,'OUTLET_03',49]])}

五类边界合计 **67,262** 个三角形。体单元只有一个 physical group：**FLUID = 100**，与 facet tag 的含义分开记录。

## 怎么生成体网格？

先合并数值完全相同的重复顶点记录，再把五类已带三角网格的 discrete surface 组成一个闭合 surface loop 和一个 volume。Gmsh 4.15.2 使用三维 Delaunay 路径（Algorithm3D=1）填充内部，`Mesh.MeshOnlyEmpty=1` 保留已有表面网格；允许 Gmsh 三维内部单元优化，未进行表面重网格、平滑、修补、降采样、cap 修改或 extension 裁剪。操作后的全部边界顶点及带标签三角形都逐个核对。[Gmsh 官方说明](https://gmsh.info/doc/texinfo/) 定义了 discrete entity、MeshOnlyEmpty 和质量指标。

{table(['顶点合并审计项','实测'],[['输入点记录',weld['points_before_welding']],['精确合并后的独立坐标',weld['unique_points_after_welding']],['合并的重复记录',weld['merged_record_count']],['最大合并距离 / 容差','0 m / 0 m'],['未被任一三角形引用的独立点记录',weld['unreferenced_unique_point_records']],['实际边界顶点',weld['mesher_boundary_vertex_count']],['合并前 / 后三角形数','67262 / 67262'],['移动的保留顶点','0']])}

40 个未引用点只是文件内多余记录，没有属于表面三角形的几何；只从 mesher 的点表中排除它们，不删除或改变任何三角形。合并顺序确定；重复执行结果一致；人工构造的相近但不相同点不会被合并。

两个 profile 的内部目标边长分别为 **1.0e-6 m** 和 **5.0e-7 m**；`optional_interior_size_control` 本次为 null，没有额外尺寸场。冻结表面边长中位数为 {adapter['boundary_edge_lengths_m']['median']*1e6:.4f} μm，已经很细，因此 coarse 不能把边界变粗。线程数固定为 1；开发限额为 150 万 tetra、300 s、8 GiB 虚拟地址空间，不自动提升资源限额。

远端只负责生成体网格、几何 QC 和 DOLFINx 转换/重载。WSL 始终是源码与报告的唯一来源；通过 Stage 0 的 probe/sync/run/fetch 脚本加 `--stage 1` 调用，每次运行前校验同步源码 SHA256。远端 Ryzen 7 7800X3D、约 61 GiB RAM，存在 RTX 4090，**GPU used = false**；Gmsh 只用单进程单线程。DOLFINx 0.11.0 的 XDMF/HDF5 保存方式依据[官方 Gmsh 示例](https://raw.githubusercontent.com/FEniCS/dolfinx/v0.11.0/python/demo/demo_gmsh.py)。没有进行 MPI 性能基准。

coarse 的独立测试与五张预览图生成后，才由 [coarse_gate.json](coarse_gate.json) 允许 medium 开始。完整运行命令、版本、输入/配置 SHA、主机、资源、stdout/stderr、退出码保存在 `logs/stage01/` 及取回的远端日志中；每个网格的 `metadata/meshing.json` 记录本次 meshing 的独立耗时与进程峰值 RSS。

## 网格是否真的填满血管内部？

两张网格均只有 **1 个通过共享面连通的流体体积**。四面体总和约为 **{qc['coarse']['volume_closure']['tetra_volume_m3']:.12e} m³**；与冻结闭合表面的有向体积比较，coarse 相对差 {qc['coarse']['volume_closure']['relative_error']:.3e}，medium 相对差 {qc['medium']['volume_closure']['relative_error']:.3e}。外表面没有额外内腔、缺失面或非流形面；完整外表面等于原始 67,262 个带标签三角形。

![内部体网格切面](tetrahedral_cutaway.png)

图中显示 medium 入口向内 4 μm 附近的真实体单元局部：左侧切去一半体积，右侧从切面正面看四面体与平面的交线。彩色多边形覆盖内部截面，未用一层空心表面冒充体网格。局部图用于理解结构；全体网格是否完整由全部单元连接关系、外表面对应和总体积核对共同检查。另保留 [coarse 同类切面](coarse_review/tetrahedral_cutaway.png)。

## 边界标签是否保留下来了？

是。两个 profile 均有 WALL 67,071、INLET 56、OUTLET_01 44、OUTLET_02 42、OUTLET_03 49 个 exterior facets。**未标记 = 0，重复标记 = 0，标在内部面上的边界标签 = 0**。每个边界三角形的外法向依据其唯一邻接 tetra 的内部位置确定，再与契约比较，未按 xyz 符号猜测。

![最终边界标签，两个视角](boundary_tags.png)

![四个端口与法向](port_closeups.png)

下表面积由最终体网格的 oriented exterior facets 直接积分。对线性三角形，叉积计算面积就是精确的常数面积积分，不需要流场函数空间。

{table(['profile','boundary','entity id','source area (m²)','FEM area (m²)','relative difference','normal dot product'],port_rows)}

{table(['端口','最大绝对面积误差 (m²)','最大中心位移 (m)','源 cap 平面偏差 (m)','最终最大平面偏差 (m)'],port_detail)}

两种网格的最大边界位移均为 **0 m**，全部边界坐标位级相等，带标签三角形 connectivity 完全对应。表面总面积为 {qc['coarse']['boundary_fidelity']['final_surface_area_m2']:.15e} m²；三角形实际几何 bounds 为 `{qc['coarse']['boundary_fidelity']['source_bounds_m']}` m，与最终网格边界相同。这里的 bounds 只统计真正被三角形引用的顶点，40 个未引用点不属于表面。比较没有只依赖 bounding box。

## 四面体有没有坏单元？

两个 profile 的 **负体积、零体积、非有限坐标/体积、非有限质量指标均为 0**，有向 Gmsh tetra 没有 inverted 单元。DOLFINx 可重新排列 tetra 的局部顶点编号，因此重载时比较每个 tetra 的完整几何并验证正的几何体积，不把合法的编号奇偶排列误判成倒置。

采用 Gmsh `minSICN`（有符号逆条件数）描述形状：**越高越好，规则四面体为 1，接近 0 表示很扁或细长**。最长/最短边比作为补充；体积本身不直接等同于质量。[指标定义见 Gmsh 官方 API](https://gmsh.info/doc/texinfo/#gmsh_002fmodel_002fmesh_002fgetElementQualities)。

![四面体质量分布](mesh_quality_distribution.png)

{table(['profile','指标','minimum','P1','P5','median','P95','maximum'],quality_rows)}

把 `minSICN < 0.1` 作为本项目的 **ADVISORY / MANUAL_REVIEW 显示触发值**，不是论文级合格线或收敛标准。低质量单元仍全部保留并报告：

{table(['profile','低质量单元数','占比','最近 INLET','最近 OUTLET_01','最近 OUTLET_02','最近 OUTLET_03','最近 WALL'],low_rows)}

多数低质量单元最接近端口 cap，图中最差的若干单元也集中在端口附近；这是与冻结 cap 扇形细长三角形相邻的现象。位置归类使用最近边界三角形中心，只用于诊断，不参与标签赋值。其余靠近 wall 的单元尚不能自动区分 junction 与窄血管，不宣称已经完成这种解剖位置判读。

{chr(10).join(worst_text)}

![最差单元及位置](worst_elements.png)

右图给出最差单元附近的真实形状，左图的红点仅为方便定位而放大。coarse 的[独立最差单元图](coarse_review/worst_elements.png)也保留。这些正体积 sliver 可能使后续高阶 FEM 的数值条件和局部误差变差；Stage 1 没有通过实际求解评估影响，也没有为改善统计数字移动冻结表面。

## coarse 和 medium 有什么区别？

{table(['profile','vertices','tetra','minimum quality','P5 quality','median quality','meshing runtime','peak RSS'],profile_rows)}

medium 四面体数量为 coarse 的 {comparison['tetrahedron_ratio_medium_to_coarse']:.3f} 倍，两个网格边界完全相同。medium 的质量中位数更高，但最差值略低，说明整体加密不能自动消除冻结 cap 邻近的坏形状。耗时只指生成脚本中 Gmsh 建模、网格生成及导出，不包括 SSH、环境探测、后续 QC 或渲染；RSS 是单进程高水位，不是 MPI 内存总和。

目前有两个可用于后续开发审查的网格；**未做网格收敛分析，也不能称为网格已收敛**。正式产物分别为 [coarse XDMF](../../outputs/stage01/coarse/mesh/fluid.xdmf) 和 [medium XDMF](../../outputs/stage01/medium/mesh/fluid.xdmf)，各自必须与同目录 `fluid.h5` 一起保存。MSH、原始 NumPy 网格、QC 与 metadata 同时保留。机器可读对照见 [mesh_comparison.json](mesh_comparison.json)。

## 做了哪些自动测试？

完整本地 pytest：**{passed} passed，{failed} failed，{skipped} skipped**，见 [pytest_results.xml](pytest_results.xml)。新增的七个 Stage 1 测试文件永久保留，覆盖冻结 SHA/VTP 标签、精确合并及近点攻击样例、正体积与单连通体、外表面完整分区、全部三角形几何/标签一致、端口面积/中心/法向/平面、保存重载以及图像内容。测试会从实际返回的网格重新计算 QC，而非只读取 PASS 字符串；缺标签、重复标签、未知标签、顶点位移、标签替换和丢面均有失败样例。

唯一 skip 是 WSL 未安装 DOLFINx 时的旧 Stage 0 环境 smoke；它不是 Stage 1 reload 的替代。Stage 1 的四次远端重载均在保存进程结束后的新 Python/MPI 进程中运行，检查 cell tags、facet tags、全体 tetra 几何和端口积分，owned cells/facets 只统计一次以避免 ghost 重计：

{table(['profile','MPI ranks','tetra','boundary facets','结果','验证'],reload_rows)}

初次 coarse 本地回归曾把总面积浮点求和要求为位级相等，出现一次约 1 ulp 的 reduction 差异；测试已改为 `rtol=1e-14, atol=0`，边界坐标和三角形对应仍要求精确相等。原始失败 XML 和日志保留。Stage 0 的 reference 测试现把临时结果写到 pytest 临时目录，避免重跑测试覆盖历史报告；检查内容不变。

完整 pytest 之后，又使用 Stage 0 的全量 SHA256、大小、权限、mtime 与目录/符号链接机制重新检查两个只读工程：

{table(['参考工程','modified','deleted','added'],[[root,val['modified'],val['deleted'],val['added']] for root,val in reference['per_root'].items()])}

共核对 {reference['entries']:,} 个参考条目；没有执行参考工程代码。Stage 0 的 {preservation['entries']} 个报告、输入、输出和日志文件也与 Stage 1 开始快照完全一致。证据见 [reference_integrity.json](reference_integrity.json) 和 [stage00_preservation.json](stage00_preservation.json)。

## 人工还应该看什么？

- **boundary_tags.png**：与 [Stage 0 原始总览](../stage00/source_geometry_overview.png) 对照，确认一个入口、三个出口、壁面及向外箭头的实际位置。
- **port_closeups.png**：四个 cap 的细长扇形三角形仍是原始几何；检查有无裂口、明显翘曲和错误标签。
- **tetrahedral_cutaway.png**：检查体积内部切面是否填满；对照 coarse 预览与 medium 正式图。
- **mesh_quality_distribution.png** 与 **worst_elements.png**：检查最差单元的 cap 邻近位置、扁薄程度，并决定这些单元是否可用于下一步开发。

当前人工审核状态为 **PENDING**。Codex 的图像检查只检查渲染是否可读，不替代用户的几何和网格质量审核。

## 还存在什么问题？

coarse 有 {qc['coarse']['quality']['advisory_count']} 个、medium 有 {qc['medium']['quality']['advisory_count']} 个触发低质量提示的单元，最小体积均约 {qc['coarse']['quality']['volume_m3']['minimum']:.3e} m³；严格正体积不等于适合任意求解器。冻结 cap 的长细三角形不能在本阶段擅自修改，进一步改善若涉及源表面必须另行确定授权范围。后续离散误差、求解稳定性和血流正确性尚未评估。

DOLFINx 0.11.0 的几何 dofmap 访问在重载日志中产生兼容性弃用提示，四次验证均正常完成；当前路径锁定此环境，未宣称兼容未来版本。另一个已明确的限制是 optional interior size control 预留为 null，本次仅验证常数内部目标尺寸。

## 是否可以进入 Stage 2？

**STAGE 1 STATUS: {status}**

所有本阶段硬性检查通过。条件是用户完成上述图像审核，并接受已公开的非退化低质量单元用于后续开发，或先指定进一步的网格改进范围。在这之前不能把本报告当作无条件的网格质量批准。

本次工作到 Stage 1 为止，**没有开始 Stage 2**。没有实现或运行 Stokes/Navier–Stokes、速度/压力空间、Real 约束、流量乘子、血流场、WSS、流线、RBC 或微泡流程。
"""
    (reports/"REPORT.md").write_text(text)
    lines=["Stage 1 completed.",""]
    for p in qc:
        lines.extend([f"{p}:",f"    tetrahedra = {qc[p]['tetrahedron_count']}",f"    vertices = {qc[p]['vertex_count']}",f"    boundary facets = {qc[p]['topology']['exterior_facet_count']}",f"    status = {qc[p]['status']}",""])
    lines.extend(["tests:",f"    passed = {passed}",f"    failed = {failed}",f"    skipped = {skipped}","", "report:","    reports/stage01/REPORT.md","","human review required:"])
    lines.extend("    "+name for name in visuals["images"])
    lines.extend(["",f"STAGE 1 STATUS: {status}","Stage 2 not started."])
    summary="\n".join(lines)+"\n"
    (reports/"terminal_summary.txt").write_text(summary)
    print(summary)
    if status=="FAIL":
        raise SystemExit(1)


if __name__=="__main__":
    main()
