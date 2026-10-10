"""Package only saved geometry / zero-step / display evidence; no solver calls."""
from pathlib import Path
import ast
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
import numpy as np
import pyvista as pv
from PIL import Image
import imageio.v2 as imageio

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pd_clot.lammps_fragmentation.y_bifurcation_sources import OUT,RESULTS,digest
from pd_clot.lammps_fragmentation.config import dump


def read(name):return json.loads((OUT/name).read_text())


def main():
    info=read('GEOMETRY_SUMMARY.json');fit=read('GEOMETRY_AND_WALL_FIT_CHECK.json');clear=read('LUMEN_INTERSECTION_AND_CLEARANCE.json')
    flow=read('ACTUAL_Y_ASSET_AUDIT.json');reg=read('SOURCE_ASSET_REGISTER.json');run=read('LAMMPS_RUN0_CHECK.json')
    for f in ['GEOMETRY_AND_WALL_FIT_CHECK.json','LUMEN_INTERSECTION_AND_CLEARANCE.json','ACTUAL_Y_ASSET_AUDIT.json','LAMMPS_RUN0_CHECK.json',
              'STRICT_GEOMETRY_AUDIT.json','FLOW_COORDINATE_CHECK.json','BACKGROUND_FLOW_OVERLAY_CHECK.json','Y_CFD_SOURCE_UNCHANGED.json','HISTORY_PROTECTION_FINAL.json']:
        assert read(f)['status']=='PASS',f
    assert run['time_steps_advanced']==0 and run['clot_nodes']==info['clot_nodes']
    assert 'run 0\n' in (OUT/'in.y_clot_geometry_check.lammps').read_text()
    surface=np.load(OUT/'SURFACE_NODES_AND_NORMALS.npz');total=surface['A_i_m2'].sum()
    assert np.isclose(total,info['exposed_area_m2'],rtol=1e-12)
    m=pv.read(RESULTS/'paraview/Y_CLOT_SCENE.vtm');assert m.n_blocks==6
    assert m['clot_points'].n_points==info['clot_nodes'] and m['anchors'].n_points==info['anchor_nodes']
    assert m['bonds'].n_cells==info['bond_count']
    raw=pv.read(RESULTS/'paraview/background_flow.vtu')
    source=next(f for f in reg['files'] if f['relative_path']==reg['states']['simvascular/steady']['snapshot'])
    native=pv.read(source['local_readonly_copy'])
    assert np.array_equal(raw.points,native.points) and np.array_equal(raw['Velocity'],native['Velocity']) and np.array_equal(raw['Pressure'],native['Pressure'])
    data=RESULTS/'FIGURE_DATA'
    for name in ['CURVED_MURAL_CLOT_GEOMETRY.npz','CURVED_WALL_ANCHORS.npz','BOND_GROUPS.npz','SURFACE_NODES_AND_NORMALS.npz']:
        shutil.copy2(OUT/name,data/name)
    files={str(p.relative_to(ROOT)):dict(sha256=digest(p),bytes=p.stat().st_size) for p in (RESULTS/'paraview').rglob('*') if p.is_file()}
    dump(data/'SOURCE_IDENTITIES.json',dict(geometry_id='CURVED_MURAL_DOME_V1',source_geometry_id=flow['geometry_id'],source_mesh_id=flow['mesh_id'],source_flow_id=flow['flow_id'],
         dataset_files=files,geometry_generated_from='actual server WALL triangles clipped to one frozen footprint; radial inward height profile',visual_only_smoothing='normal interpolation, no point displacement',CLOT_PRESENT_CFD_SOLVED='NOT_RUN'))
    reader=imageio.get_reader(RESULTS/'y_clot_static_camera_orbit.mp4');count=reader.count_frames();assert count==72
    sheet=Image.new('RGB',(1280,960),'white')
    for k,frame in enumerate([0,18,36,54]):sheet.paste(Image.fromarray(reader.get_data(frame)).resize((640,480)),((k%2)*640,(k//2)*480))
    reader.close();sheet.save(OUT/'MOVIE_QA_CONTACT_SHEET.png')
    gif=Image.open(RESULTS/'y_clot_static_camera_orbit.gif');assert gif.n_frames==72
    html=(RESULTS/'interactive_y_clot_scene.html').read_text();assert len(html)>100000 and '<html' in html.lower()
    qa=read('VISUALIZATION_QA.json')
    # Human/model visual inspection is recorded separately after opening these PNGs.
    review=read('VISUAL_INSPECTION.json');assert review['status']=='PASS'
    qa.update(status='PASS',visual_review=review,paraview_scene_reloaded=True,raw_background_arrays_bitwise_equal=True,
              mp4_decoded_frames=count,gif_frames=gif.n_frames,interactive_html='SELF_CONTAINED_EXPORTED;_BROWSER_INTERACTION_NOT_AUTOMATED',
              html_sha256=digest(RESULTS/'interactive_y_clot_scene.html'))
    dump(OUT/'VISUALIZATION_QA.json',qa)
    statuses=dict(REAL_Y_VESSEL_MESH='VERIFIED',UNOBSTRUCTED_STEADY_CFD_IMPORTED='VERIFIED',Y_CFD_SOURCE_UNCHANGED='PASS',
                  MOTHER_CURVED_WALL_PLACEMENT='PASS',MURAL_CLOT_SHAPE='READY',LAMMPS_CUBIC_POINTS='READY',CURVED_WALL_INTERFACE_BONDS='READY',
                  LAMMPS_DATA_RUN0='PASS',GEOMETRY_WALL_CLEARANCE='PASS',BACKGROUND_CFD_OVERLAY='PASS',SURFACE_GEOMETRY_LOAD_SCHEMA='READY',
                  CLOT_PRESENT_CFD_SOLVED='NOT_RUN',CLOT_SURFACE_TRACTION_VALIDATED='NOT_ESTABLISHED',CURVED_INTERFACE_STRENGTH_VALIDATED='NOT_ESTABLISHED',
                  CLOT_FRAGMENTATION_DYNAMICS='NOT_RUN',BUBBLE_STREAMING_CONNECTED='NOT_YET',PHASE_FIELD_B1_B2='TERMINATED',PD2_FRACTURE='STOPPED',
                  BIOLOGICAL_QUANTITATIVE_VALIDATION='NOT_ESTABLISHED',PATIENT_SPECIFIC_CLOT='NOT_ESTABLISHED',PRIMARY_USE='EXPLORATORY_GEOMETRY_AND_BACKGROUND_FLOW_INTEGRATION')
    final=dict(task='Y_BIFURCATION_CLOT_PLACEMENT_AND_VISUAL_INTEGRATION_V1',status='READY_FOR_HUMAN_REVIEW',**statuses,
               shape_count=1,clot_nodes=info['clot_nodes'],wall_anchors=info['anchor_nodes'],bonds=info['bond_count'],
               geometry=info,run0=run,own_jobs_remaining=[],automatic_followup=False,git_commit_or_push=False,
               run0_invocations=2,run0_scope='Earlier draft and corrected final geometry; both zero time steps. Earlier evidence retained.',
               surface_mapping_node_count=len(surface['surface_node_id']))
    dump(OUT/'FINAL_STATUS.json',final)
    perf=dict(geometry=read('GEOMETRY_PERFORMANCE.json'),final_run0=run,
              rendering=read('VISUALIZATION_PERFORMANCE.json'),interactive_export=read('INTERACTIVE_EXPORT_PERFORMANCE.json'),
              historical_hash_check=read('HISTORY_PROTECTION_FINAL.json'),
              source_transfer_timing='NOT_MEASURED',development_attempts='Preserved in development_attempts; interrupted first geometry preprocessing RSS observed, full resource peak not established',
              scientific_dynamic_steps=0,new_CFD_solves=0,new_bubble_solves=0,resource_scope='Actual geometry / initialization / rendering only; no benchmark campaign')
    dump(OUT/'PERFORMANCE_REAL_TASK.json',perf)
    src=OUT/'source_snapshot';src.mkdir(exist_ok=True)
    for p in (ROOT/'pd_clot/lammps_fragmentation').glob('y_*.py'):
        ast.parse(p.read_text());shutil.copy2(p,src/p.name)
    for name in ['interface_model.py','calibrate_failure.py','lammps_writer.py','config.py']:
        shutil.copy2(ROOT/'pd_clot/lammps_fragmentation'/name,src/name)
    dump(OUT/'SOURCE_CODE_IDENTITY.json',{p.name:digest(p) for p in src.glob('*.py')})
    original_volume=info['old_hemisphere_analytic_volume_m3']
    config_id=next(f['sha256'] for f in reg['files'] if f['relative_path']=='configs/stage4/bifurcation.yaml')
    report=f'''# Y 分叉附壁血栓装配与可视化审查

结果：READY_FOR_HUMAN_REVIEW。本轮只构造一个 CURVED_MURAL_DOME_V1，完成实际 Y 网格导入、曲壁血栓、立方格材料点与三层壁锚点、新连接、零步初始化及背景 CFD 展示。没有启动含血栓 CFD、微泡或断裂动力学。CLOT_PRESENT_CFD_SOLVED = NOT_RUN。

1. 原始文件和身份

已通过 SSH 在服务器 `/workspace/mevo_acoustic_microbubbles` 逐级读取 config、geometry/mesh/steady/import STATE 及其实际引用文件，下载 15 项只读快照并逐项核对 SHA256。15 项也均与本地平行项目镜像相同。服务器所查项目根下没有当前 clot_clearance；新开发和零步检查在已有本地血栓工作树 `{ROOT}` 执行，没有在服务器伪造血栓仓库或修改平行项目。

- 配置：`/workspace/mevo_acoustic_microbubbles/configs/stage4/bifurcation.yaml`；SHA256 `{config_id}`。
- 原四面体网格：`/workspace/mevo_acoustic_microbubbles/{reg['states']['simvascular/mesh']['mesh']}`；SHA256 `{flow['mesh_id']}`。
- 原稳态场：`/workspace/mevo_acoustic_microbubbles/{reg['states']['simvascular/steady']['snapshot']}`；SHA256 `{flow['flow_id']}`。
- WALL 与三个端口位于网格同目录 `mesh-surfaces/`；每个原文件的路径、大小、mtime、完整 hash 和只读快照路径见 SOURCE_ASSET_REGISTER.json / SOURCE_CONFIG.json。

ACTUAL_Y_ASSET_AUDIT.json 的 {len(flow['rows'])} 项手交信息均有 ACTUAL_FILE_VERIFIED 记录；配置来源和实际网格测量分别标注。母管长度 14 mm、半径 1.3 mm；两子管长度 12 mm、半径 1 mm、偏转 ±30°。真实体网格 {native.n_points:,} 节点、{native.n_cells:,} 四面体。密度 1000 kg/m³、黏度 0.001 Pa·s、入口抛物线通量归一化和两个零参考自然牵引出口同时从实际 solver.xml 核对。原场第 121 步、dt=0.01 s，时间 1.21 s。此身份与 BraVa 四端口 18 mL/min 案例分开。

2. 放置位置、壁面拟合与间隙

中心线投影固定 (-6,0,0) mm，名义脚点 (-6,0,1.3) mm，实际 WALL 径向交点 (-6,0,{fit['actual_footpoint_m'][2]*1000:.9f}) mm。脚点修正 {fit['nominal_footpoint_correction_m']*1e6:.6f} µm，低于 dx=50 µm，不移动轴向位置。通过裁剪/线性细分实际 WALL 三角形获得连续足迹，每个足迹点沿局部径向向内使用指定 h=H(1−s²)^(3/2) 高度，未另造 Y 管。

足迹是一个连通面片；面积 {info['attachment_area_m2']*1e6:.6f} mm²。附着顶点对原壁面最大/中位距离为 {fit['attachment_vertex_gap_max_m']*1e6:.9f}/{fit['attachment_vertex_gap_median_m']*1e6:.9f} µm。最大表面正向偏差 {fit['maximum_surface_penetration_m']*1e6:.9f} µm，小于预定 0.002 µm 数值容差；这是对同一离散表面的重合精度，不是生物壁面的纳米精度。材料中心离壁最小 {fit['particle_center_wall_distance_min_m']*1e6:.3f} µm，与体素/积分点离散有关，不能把中心距误报为连续实体悬空。

最终全部 {info['clot_nodes']:,} 材料点及 {fit['voxel_support_samples_checked']:,} 个体素中心/角点/面中点/边中点通过严格四面体归属；{info['anchor_nodes']:,} 个锚点中心在原管腔外。表面、三端口与对侧壁的距离记录在 LUMEN_INTERSECTION_AND_CLEARANCE.json；入口最近顶点距离约 {clear['minimum_distance_to_port_surface_vertices_m']['INLET']*1000:.3f} mm，未触及端口或分叉区域。锚点是固定壁代理点，未声称其名义体积构成真实血管壁。

3. 新形状与旧平面半球

| 指标 | 曲壁附壁隆起 | 旧平面半球 |
|---|---:|---:|
| 名义 a、b、H / 半球 R | 0.5625 mm | 0.5625 mm |
| 连续实体体积 | {info['continuous_mesh_volume_m3']*1e9:.6f} mm³（实际闭合网格） | {original_volume*1e9:.6f} mm³（2πR³/3） |
| LAMMPS 积分体积 | {info['clot_volume_m3']*1e9:.6f} mm³ | {info['old_hemisphere_nominal_volume_m3']*1e9:.6f} mm³（历史保存结果） |
| 材料点数 | {info['clot_nodes']} | {info['old_hemisphere_clot_nodes']} |
| 附着形态 | 圆柱实际曲壁上的连续椭圆足迹 | 平面圆盘 |

新连续体积约为旧解析半球的 {info['continuous_mesh_volume_m3']/original_volume*100:.2f}%，不能称等体积变形。足迹轴向跨度 {info['footprint_axial_extent_m']*1000:.6f} mm、名义圆柱展开周向跨度 {info['footprint_circumferential_arc_extent_m']*1000:.6f} mm；保存网格最大高度 {info['maximum_mesh_height_m']*1000:.6f} mm。外露面积 {info['exposed_area_m2']*1e6:.6f} mm²，总表面积 {info['total_surface_area_m2']*1e6:.6f} mm²。材料中心 bbox、连续实体 bbox 与质心见 GEOMETRY_SUMMARY.json。

连续实体中心截面 x=−6 mm 面积 {clear['center_section_clot_area_m2']*1e6:.6f} mm²，除以同位置实际管腔截面 {clear['center_section_lumen_area_m2']*1e6:.6f} mm² 得 {clear['center_section_area_occupancy']*100:.4f}%。这是几何截面占据率，不能当水力阻力或流量变化。另列材料柱的 dx² 占据率 {clear['discrete_center_column_area_occupancy']*100:.4f}%。

连续边缘平滑降至零；有限 dx 无法填满全部薄缘。保守掩码去除 {info['conservative_wall_voxel_removals']} 个不满足壁面支持点条件的候选体素，另去除 {info['removed_lattice_whiskers']} 个单邻居细须。积分体积相对连续网格相差 {(1-info['clot_volume_m3']/info['continuous_mesh_volume_m3'])*100:.2f}%，没有进行分辨率收敛或擅自改 dx。医学实体视图显示生成材料点所用的同一连续参考网格，求解审查视图显示实际保留点和连接，二者差异明确保留。

4. LAMMPS 装配和零步读取

- 材料点 {info['clot_nodes']:,}；固定锚点 {info['anchor_nodes']:,}，三层分别 {info['wall_layers_counts']['1']}/{info['wall_layers_counts']['2']}/{info['wall_layers_counts']['3']}。
- BULK {info['bulk_bonds']:,}；INTERFACE {info['interface_bonds']:,}；WALL {info['wall_bonds']:,}；总连接 {info['bond_count']:,}。
- 三层按真实壁面外侧 0–50、50–100、100–150 µm 法向距离壳层选点，所有坐标保持全局 dx=50 µm 立方格；不移动格点贴墙。每点体积 dx³=1.25e−13 m³，每点质量 1.25e−10 kg。血栓总质量 {info['clot_mass_kg']:.9g} kg，排除锚点。
- δ=150.5 µm；无重复点、重复连接、孤立点；实际最大连接长度 150 µm。10 个候选 BULK 跨自由面连接被剔除；界面连接必须穿过实际 attachment patch 且不穿过 exposed surface。全部 WALL 连接另经精确线段/管壁相交检查，穿腔数为 0。
- 复用原 write_typed_data、LPS prelude、SI 内力密度到节点力转换与现有变体，未复制第二套求解器。变体二进制 SHA256 `{run['executable_sha256']}` 未改变。
- 最终 `read_data` + `run 0` 成功，{run['atoms']:,} 个 atom、{run['bonds']:,} 个 bond 与 NPZ 逐项对应，只有 step=0、time=0，坐标不变，速度/损伤/初始力为零。无时间积分 fix 或外加载荷。

开发记录保留：初步 VTK FindCell 曾接受一个略越壁的支持点，严格重心坐标检查拒绝后增强到每体素 27 点及 1 nm 保守裕量，最终去除该体素及一个细须点。先前草稿的 run 0 和图保存在 development_attempts，正式数据以本目录顶层 FINAL_STATUS 为准。总计两次零步初始化、零次动力学；没有多位置或多形态扫参。

5. 冻结 CFD 展示能说明什么

流线与箭头从原 VTU 速度数组生成，显示原管腔内的流向、速度分布，以及新血栓相对于原流场的位置。56 条显示流线中 4 条在第一次碰到新固体时截断，全部保存线段中点均在新固体之外。没有构造绕流、涡旋或假微泡。所有背景图与交互视图明确标记 UNOBSTRUCTED BACKGROUND CFD — NOT CLOTTED FLOW FIELD。

不能据此得到血栓造成的压降、分流变化、近血栓壁面剪切、表面准确 traction、微泡聚集或血栓破裂。background_flow.vtu 保留原无血栓场；新实体内部的原场只代表历史数值，不能作为现态流体。对可见流线做 mask 不会产生新的物理解。

6. 原流量和分流是否保留

从实际端口网格及原速度重新积分，入口约 {-flow['boundary_outward_flux_m3_s']['INLET']*6e7:.8f} mL/min；A/B 约 {flow['boundary_outward_flux_m3_s']['OUTLET_A']*6e7:.8f}/{flow['boundary_outward_flux_m3_s']['OUTLET_B']*6e7:.8f} mL/min，仍为原近各半分流。全域原速度范围 {flow['velocity_range_m_s']} m/s，压力范围 {flow['pressure_range_Pa']} Pa。原字段和导出背景字段的坐标、Velocity、Pressure 数组逐项相等。这些全部是原无血栓结果。

7. 展示和审查来源

四张 PNG/PDF、72 帧/6 s 镜头环绕 MP4/GIF、交互 HTML 和 ParaView VTM 均离线读取同一份保存坐标。仅在显示时统一 m→mm，各方向同尺度；没有单独放大或移动血栓。光滑效果仅插值表面显示法向，未移动实体点或 LAMMPS 点。局部视图从管腔侧显示隆起，并提供真实 y=0 剖面。

自动检查了 PNG 非空、PDF、MP4/GIF 帧、ParaView 数据重载、源哈希与图形数组对应；四张主图和动画四视角拼图人工/模型目视检查通过。HTML 已导出内嵌场景，未进行浏览器交互自动化；ParaView 数据已实际重新读入。首次交互导出缺插件，只在本任务 visual_dependencies 隔离安装 trame-pyvista / nest_asyncio2，未改变平行项目环境。

8. 未来载荷 provider 契约

SURFACE_NODES_AND_NORMALS.npz 给出 {len(surface['surface_node_id'])} 个初始外露面材料节点的全局坐标、外露法向、实际表面积权重 A_i、到真实壁面的距离、投影后的轴向切向、周向切向、壁到管腔方向、原局部流向及 geometry/mesh/field 身份。只将外露面三角形面积分配一次到最近外露材料节点，A_i 总和 {total:.12g} m²；不加入附着面面积，不把全部体积节点当受载表面。薄缘节点可有两种拓扑邻近标签，面积集合仍互不重复。

provider 必须输出 load_source、surface_node_id、position_xyz_m、time_s/时间序列、source_geometry_id、source_flow_id、coordinate_frame、source_units、validity_domain，以及 force_xyz_N，或 traction_xyz_Pa 和 nodal_area_m2。后者逐节点用 F_i=t_i A_i 转换。未来沿真实曲壁切平面投影方向，不能继续用旧平底半球的原点径向。当前数据为初始映射，不支持断裂新增表面自动更新。

BACKGROUND_Y_FLOW = REQUIRES_CLOT_PRESENT_CFD_FOR_TRACTION；ANALYTIC_STREAMING 仅 schema；BUBBLE_STREAMING / CFD_STREAMING = NOT_YET。本轮没有把任何这些载荷施加给 LAMMPS。

9. 材料、接触和界面适用范围

保留 Liu 同源血栓未损伤小应变映射 K=33333.333333333336 Pa、G=940 Pa，工程密度 1000 kg/m³，bulk s00=0.15、interface s00=0.140624859375、α=0，以及已存在的 MODIFIED_1P0DX 接触变体。1.0dx 是有限工况下的工程接触选择，不是生物接触距离。界面是同文献猪皮剥离 proxy；曲壁法向、邻域和连接组织均变了，旧平面宏观强度不能自动迁移为本场景认证。

CURVED_INTERFACE_STRENGTH_VALIDATED = NOT_ESTABLISHED；材料映射、曲壁界面、接触及像素/网格展示都不构成临床、患者特异或实验精度验证。固定锚层不代表可变形血管壁或 FSI。

10. 真正的含血栓 traction 需要什么

下一阶段若获授权，应以原管腔减去血栓实体形成新的流体域，为暴露血栓表面建立独立 boundary tag，重建相容体网格及近表面解析度；移除附着区域的流体空间并核对端口和流固边界拓扑。可把原入口流量与出口牵引作为清楚声明的初始 BC，但必须重新求解新的含血栓速度/压力，再由总应力 σ=−pI+μ(∇u+∇uᵀ) 及统一法向约定映射表面牵引，并检查表面积、总力、单位和局部误差。若血栓会变形或断裂，另需定义几何更新与耦合策略。

原 CFD 的求解域包括现在的固体区域，且未在新 exposed surface 上施加流固无滑移边界，因此不能把旧场压力/梯度直接命名为新血栓上的真实 traction。以上仅列出必要后续，不在本轮执行。

交付与停止

- results 目录 OPEN_RESULTS.html：四图、视频、交互 HTML 和 ParaView 入口。
- verification 顶层 FINAL_STATUS.json、输入 data/deck、run0 日志、全部几何/流场/source 检查及 evidence.zip。
- 历史 {read('HISTORY_PROTECTION_FINAL.json')['protected_files']} 个文件、{read('HISTORY_PROTECTION_FINAL.json')['protected_bytes']} bytes 哈希保持；15 项服务器和本地镜像原资产保持。没有 commit/push 或未完成本轮作业。
- 现阶段主用途：EXPLORATORY_GEOMETRY_AND_BACKGROUND_FLOW_INTEGRATION。停止并等待人工审查，不自动执行 CFD、bubble 或 fracture。

实现参考：[LAMMPS BPM/PERI](https://docs.lammps.org/bond_bpm_peri.html)、[read_data](https://docs.lammps.org/read_data.html)、[PyVista UnstructuredGrid](https://docs.pyvista.org/api/core/_autosummary/pyvista.unstructuredgrid)。官方 BPM/PERI 的等间距立方格、节点体积和读取规则用于输入装配；工程接触变体的身份单独保存。
'''
    name='REVIEW_Y_BIFURCATION_CLOT_PLACEMENT_AND_VISUAL_INTEGRATION_V1.md'
    (OUT/name).write_text(report);shutil.copy2(OUT/name,RESULTS/name)
    cards=[('y_bifurcation_clot_overview','真实 Y 管与血栓位置'),('curved_mural_clot_wall_fit','曲壁贴合与真实剖面'),('original_background_flow_with_clot_overlay','原始无血栓背景流叠加'),('lammps_curved_interface_audit','LAMMPS 原始点与连接')]
    html='''<!doctype html><html lang="zh"><meta charset="utf-8"><title>Y 血管附壁血栓审查</title><style>body{max-width:1150px;margin:35px auto;padding:0 22px;font:17px/1.65 system-ui;color:#243347;background:#fbfcfe}img,video{width:100%}a{color:#006e87}section{margin:32px 0;padding:20px;background:white;border:1px solid #d9e1e8;border-radius:14px}strong{color:#7e1428}</style><h1>Y 形血管 · 曲壁附壁血栓</h1><p><strong>原始无血栓背景 CFD，未重算含血栓流场。无断裂动力学。</strong></p>'''
    html+=f'<p>材料点 {info["clot_nodes"]:,} · 固定锚点 {info["anchor_nodes"]:,} · 三类连接 {info["bond_count"]:,} · run 0 通过</p><p><a href="{name}">中文审查报告</a> · <a href="interactive_y_clot_scene.html">交互式三维场景</a> · <a href="paraview/Y_CLOT_SCENE.vtm">ParaView 场景</a></p>'
    for f,title in cards:html+=f'<section><h2>{title}</h2><a href="{f}.pdf">PDF</a><img src="{f}.png" alt="{title}"></section>'
    html+='<section><h2>静态镜头环绕（6 秒）</h2><video controls loop muted src="y_clot_static_camera_orbit.mp4"></video><p>仅相机运动；几何、流场和血栓不随时间演化。</p><a href="y_clot_static_camera_orbit.gif">GIF</a></section></html>'
    (RESULTS/'OPEN_RESULTS.html').write_text(html)
    (RESULTS/'paraview/README.md').write_text('Open Y_CLOT_SCENE.vtm together with its Y_CLOT_SCENE/ directory. Coordinates are SI metres. Hide bonds / points for solid view; show them for audit. Background streamlines stop at the clot. background_flow.vtu is the ORIGINAL unobstructed field, not recomputed. No surface traction is certified.\n')
    (OUT/'REPRODUCE.md').write_text(f'''# Reproduction boundaries

Root: {ROOT}
Python: /home/lzy/projects/computation_examples/_shared/python_environment/bin/python
The original external project is read-only. source_assets holds verified server snapshots.
Y_CLOT_REQUEST.json parameterizes position, footprint, height and physical constants. SOURCE_CONFIG.json gives external paths and identity transform. One request was executed.

Geometry entry: `python -m pd_clot.lammps_fragmentation.y_scene_task geometry --configuration verification/{OUT.name}/Y_CLOT_REQUEST.json`.
Initialization entry: `python -m pd_clot.lammps_fragmentation.y_scene_task run0`; completed run0 is guarded against accidental repetition. Use a fresh output task directory for a separately authorized rerun; preserve this evidence.
Offline figure entry: `PYTHONPATH=verification/{OUT.name}/visual_dependencies:. python -m pd_clot.lammps_fragmentation.y_scene_visualization`.
The figure entry reads saved files, traces the stored vector field only and never launches CFD or dynamics.
Installed export dependencies are isolated in visual_dependencies; shared / parallel environment unchanged.
The lightweight ZIP excludes native source VTU/VTP copies, full original background_flow.vtu and superseded draft renders; exact source paths/hashes remain in the register. Full local ParaView output is retained.
''')
    end='\n'.join(subprocess.check_output(['git',*args],cwd=ROOT,text=True) for args in [['branch','--show-current'],['rev-parse','HEAD'],['status','--short']])
    (OUT/'ENDING_GIT_STATE.md').write_text(end)
    dump(OUT/'DELIVERY_CHECK.json',dict(status='PASS',required_statuses=statuses,all_arrays_finite=True,paraview_reloaded=True,
          surfaces_and_groups_correspond=True,only_run_zero=True,physics_source_unchanged=True,source_hashes_unchanged=True,
          images_visually_inspected=True,movie_frames_decoded=72,html_exported=True,browser_automation='NOT_RUN',source_code_syntax='PASS'))
    paths=[]
    for parent in [OUT,RESULTS]:
        for p in parent.rglob('*'):
            if not p.is_file():continue
            rel=p.relative_to(parent)
            if p.name.startswith('evidence.zip') or p.name in ['HASH_MANIFEST.json','FINISH_DELIVERY.log','PACKAGE_VERIFICATION.json'] or '__pycache__' in p.parts:continue
            if 'source_assets' in rel.parts or 'visual_dependencies' in rel.parts or 'rendered_draft' in rel.parts:continue
            if p.name=='background_flow.vtu':continue
            paths.append(p)
    manifest={str(p.relative_to(ROOT)):dict(sha256=digest(p),bytes=p.stat().st_size) for p in paths}
    dump(OUT/'HASH_MANIFEST.json',manifest);paths.append(OUT/'HASH_MANIFEST.json')
    with zipfile.ZipFile(OUT/'evidence.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in paths:z.write(p,str(p.relative_to(ROOT)))
    with zipfile.ZipFile(OUT/'evidence.zip') as z:assert z.testzip() is None
    (OUT/'evidence.zip.sha256').write_text(digest(OUT/'evidence.zip')+'  evidence.zip\n')
    print('DELIVERY_READY',len(paths),(OUT/'evidence.zip').stat().st_size,flush=True)


if __name__=='__main__':main()
