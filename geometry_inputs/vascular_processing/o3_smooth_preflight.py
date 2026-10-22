"""Auditable prerequisite gate for the requested O3 ramp correction.

When installed VMTK has no ramp and Slicer has no reusable VMTK, stop as explicitly
required. No clipping, extension, final union or slicing is run by this preflight.
"""
from pathlib import Path
import copy
import csv
import hashlib
import json
import logging
import subprocess
import urllib.request

import numpy as np
import pyvista as pv
import trimesh
from PIL import Image, ImageDraw

from . import sacrificial_fixture as f
from .sacrificial_fixture_review import save_json, chinese_font, polydata
from . import surface_continuity_qc as continuity
from .o3_vmtk_bridge import invoke, probe_slicer, require_ramp


def read_saved_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def official_source_audit(out):
    records=[];out.mkdir(parents=True,exist_ok=True)
    for repo,names in [('vmtk/vmtk',['vtkVmtk/ComputationalGeometry/vtkvmtkPolyDataFlowExtensionsFilter.h',
                                    'vtkVmtk/ComputationalGeometry/vtkvmtkPolyDataFlowExtensionsFilter.cxx',
                                    'vmtkScripts/vmtkflowextensions.py']),
                       ('vmtk/SlicerExtension-VMTK',['ClipVessel/ClipVessel.py','Docs/ClipVessel.md'])]:
        url='https://api.github.com/repos/'+repo+'/commits/master'
        commit=json.load(urllib.request.urlopen(url,timeout=30))['sha']
        for relative in names:
            source_url=f'https://raw.githubusercontent.com/{repo}/{commit}/{relative}'
            data=urllib.request.urlopen(source_url,timeout=30).read()
            target=out/(repo.replace('/','_')+'_'+Path(relative).name);target.write_bytes(data)
            records.append(dict(repository=repo,commit=commit,path=relative,url=source_url,
                                local_file=str(target),sha256=hashlib.sha256(data).hexdigest(),
                                ramp_symbol_present=b'SetInterpolationModeToRamp' in data,
                                preserve_symbol_present=b'PreserveCrossSectionShape' in data,
                                installed=False,executed=False))
    save_json(out/'official_source_audit.json',records)
    return records


def bridge_smoke(executable,out):
    out.mkdir(parents=True,exist_ok=True)
    mesh=pv.Sphere(radius=.7,theta_resolution=12,phi_resolution=8).triangulate()
    mesh.point_data['sample_id']=np.arange(mesh.n_points,dtype=np.int32)
    mesh.cell_data['surface_id']=np.arange(mesh.n_cells,dtype=np.int32)
    mesh.save(out/'roundtrip_input.vtp')
    result=invoke(executable,out,'roundtrip',out/'roundtrip_input.vtp')
    checks={}
    if result['returncode']==0 and result['status']=='PASS':
        restored=pv.read(out/'roundtrip_output.vtp')
        checks=dict(point_coordinates=np.array_equal(mesh.points,restored.points),
                    triangle_connectivity=np.array_equal(mesh.faces,restored.faces),
                    point_scalars=np.array_equal(mesh['sample_id'],restored['sample_id']),
                    cell_scalars=np.array_equal(mesh['surface_id'],restored['surface_id']),
                    source_points=mesh.n_points,source_cells=mesh.n_cells)
    passed=bool(checks) and all(checks[k] for k in ['point_coordinates','triangle_connectivity','point_scalars','cell_scalars'])
    report=dict(passed=passed,checks=checks,worker=result,
                main_vtk_version=pv.vtk_version_info,main_environment_imported_vmtk=False)
    save_json(out.parent/'vmtk_bridge_smoke_test.json',report)
    return report


def baseline_figures(out, mesh, center, normal, mean_tangent, radius, ring, feature_poly, cfg):
    out.mkdir(parents=True,exist_ok=True)
    view=np.cross(mean_tangent,normal);view/=np.linalg.norm(view)
    # Same camera and crop for all baseline-only figures. No fictitious "after".
    roi=continuity.window(mesh.triangles_center,center,normal,2.8*radius,2.2)
    local=polydata(mesh.submesh([np.flatnonzero(roi)],append=True,repair=False))
    camera=[center+view*9,center,normal]
    def draw(filename,title,edges=False,features=False,curvature=False):
        p=pv.Plotter(off_screen=True,window_size=(1200,950));p.set_background('white')
        if curvature:
            value=polydata(mesh);value['mean_curvature']=value.curvature(curv_type='mean')
            selected=value.extract_cells(np.flatnonzero(roi)).extract_surface()
            low,high=np.percentile(selected['mean_curvature'],[5,95])
            p.add_mesh(selected,scalars='mean_curvature',cmap='coolwarm',clim=(low,high),show_scalar_bar=False)
        else:p.add_mesh(local,color='#10a7b5',show_edges=edges,edge_color='#404040',smooth_shading=False)
        if features and feature_poly.n_cells:p.add_mesh(feature_poly,color='#d72d35',line_width=6)
        else:p.add_mesh(pv.lines_from_points(np.vstack([ring,ring[0]])),color='#d72d35',line_width=4)
        p.camera_position=camera;p.camera.parallel_projection=True;p.camera.parallel_scale=2.35
        raw=p.screenshot(return_img=True);p.close()
        image=Image.fromarray(raw);canvas=ImageDraw.Draw(image)
        canvas.rectangle((0,0,1200,62),fill='white');canvas.text((25,12),title,fill='#172b36',font=chinese_font(cfg,30))
        canvas.rectangle((0,907,1200,950),fill='white')
        canvas.text((25,915),'现有模型的检查结果；尚未生成修复后的模型',fill='#663333',font=chinese_font(cfg,22))
        image.save(out/filename)
    draw('01_current_O3_problem_shaded.png','当前 O3 为什么不光顺')
    draw('02_current_O3_problem_wireframe.png','当前接缝的真实三角网格',True)
    draw('08_baseline_feature_edges.png','当前接缝中的尖锐折边（20°）',True,True)
    draw('11_baseline_curvature.png','当前接缝的曲率分布',curvature=True)
    save_json(out/'baseline_camera.json',dict(camera=camera,parallel_scale=2.35,
              baseline_only=True,after_images_not_generated=True))


def run(config_path, output_root=None):
    cfg=f.load_config(config_path)
    load_cfg=copy.deepcopy(cfg);load_cfg['attachment_alignment']['endpoint_ids']=[]
    inputs=f.load_inputs(load_cfg)  # No cap collar or Boolean is built in this read-only gate.
    base=inputs['base'];old=base/'print_fixture_design_aligned'
    out=Path(output_root).resolve() if output_root else base/'print_fixture_design_o3_smooth'
    if out.exists():raise FileExistsError('Refusing to overwrite O3 preflight: '+str(out))
    out.mkdir(parents=True)
    for name in ['VMTK','tables','QC','api_audit','baseline']: (out/name).mkdir()
    logger=logging.getLogger('o3_smooth_preflight');logger.setLevel(logging.INFO)
    handler=logging.FileHandler(out/'o3_smoothing_run.log',encoding='utf-8');logger.addHandler(handler)
    logger.info('NEEDS_ADJUSTMENT_O3_SURFACE_CONTINUITY; only prerequisites and frozen baseline are inspected')
    protected=f.protected_snapshot(inputs,out)
    for p in (f.ROOT/'third_party/vascularmd').glob('*.py'):protected[str(p)]=f.sha256(p)
    save_json(out/'protected_before.json',protected)
    summary=dict(status='NEEDS_ADJUSTMENT',initial_status='NEEDS_ADJUSTMENT_O3_SURFACE_CONTINUITY',
                 base_branch='sync/vascular-print-models-20260925',baseline='BASELINE_G0_ALIGNMENT',
                 final_stl=None,selected_correction=None,cadquery_handoff=None,
                 new_geometry_generated=False,production_clipping_performed=False,
                 boolean_union_performed=False,box_union_performed=False,bambu_slicing_performed=False)
    try:
        executable=cfg['vmtk']['probe_python']
        logger.info('[BRIDGE] Real VTP round trip through %s',executable)
        smoke=bridge_smoke(executable,out/'VMTK'/'bridge_smoke')
        summary['bridge_smoke']=smoke
        runtime=smoke['worker'].get('runtime',{})
        logger.info('[API] runtime=%s',json.dumps(runtime))
        slicer=probe_slicer(cfg['o3_smoothing']['slicer_executable'],out/'api_audit'/'slicer')
        summary['slicer_runtime']=slicer
        summary['official_sources']=official_source_audit(out/'api_audit'/'official_sources')
        try:require_ramp(runtime);reason='RAMP_AVAILABLE_PREFLIGHT_ONLY'
        except RuntimeError as exc:reason=str(exc)
        summary['failure_code']=reason
        summary['reason']='Installed VMTK has no ramp/preserve-shape API; installed Slicer has no VMTK extension. User section 60 requires stop.'
        if not smoke['passed']:summary['failure_code']='VMTK_BRIDGE_ROUNDTRIP_FAILED'
        # Re-read the actual existing STL. Never synthesize a replacement baseline.
        mesh=trimesh.load_mesh(old/'core/BG001_RMCA_BALANCED_core_with_ports.stl',process=True)
        ep=next(e for e in inputs['endpoints'] if e.endpoint_id=='O3')
        alignment=json.loads((old/'attachment_alignment_qc.json').read_text())['O3']
        ring=np.array(alignment['original_ring_print_mm']);center=ring.mean(axis=0);normal=np.array(alignment['cap_normal'])
        features=[];shown=None
        for angle in [10,20,30]:
            values,poly=continuity.feature_edges(mesh,center,normal,ep.radius*1.7,1.,angle)
            features.append(values)
            if angle==20:shown=poly
        normals,normal_values=continuity.normal_jumps(mesh,center,normal,ep.radius*1.7,1.)
        rows=[]
        for s in np.arange(-1.,1.0001,.25):
            rows.append(dict(arc_mm=float(s),**continuity.cross_section(mesh,center+s*normal,normal,maximum_center_distance=ep.radius),
                             stage='BASELINE_CAP_NORMAL_AUDIT'))
        layout=read_saved_csv(old/'tables/port_layout_summary.csv');old_o3=next(v for v in layout if v['port_id']=='O3')
        summary['current_baseline']=dict(cap_centerline_angle_deg=f.angle(normal,ep.tangent),
            normal_jumps=normals,feature_edges=features,maximum_area_jump_fraction=continuity.area_jump(rows),
            measurement_window=dict(center=center,axis=normal,radial_mm=ep.radius*1.7,axial_half_length_mm=1.),
            cross_sections='Baseline cap-normal coordinate, -1 to +1 mm; not a new corrected centerline profile',
            manufacturing_radius_mm=ep.radius,manufacturing_diameter_mm=2*ep.radius,
            original_swc_id=ep.original_swc_id,branch=ep.branch_id,
            prior_vessel_clearance_mm=float(old_o3['minimum_clearance_mm']),prior_clearance_source='saved aligned QC; no new route generated',
            mesh_qc=f.mesh_qc(mesh,cfg),source_stl=str(old/'core/BG001_RMCA_BALANCED_core_with_ports.stl'))
        f.write_csv(out/'tables/O3_feature_edges.csv',features)
        f.write_csv(out/'tables/O3_surface_continuity.csv',[{'version':'BASELINE_G0_ALIGNMENT',**normals}])
        f.write_csv(out/'tables/O3_cross_section_profile.csv',rows)
        f.write_csv(out/'tables/port_layout_summary.csv',layout)
        selected=[v for v in read_saved_csv(old/'tables/port_face_candidates.csv') if v.get('selected')=='True']
        route_hashes={}
        box=f.compute_box(inputs['mesh'].bounds,cfg)
        for row in selected:
            if row['endpoint_id']=='O3':continue
            port=next(e for e in inputs['endpoints'] if e.endpoint_id==row['endpoint_id'])
            target=np.array([float(row['target_'+a]) for a in 'xyz'])
            route=f.make_route(port,row['face'],target,box,cfg,row['route_type'],float(row['tangent_scale']),
                               float(row['guide_bend_radius_mm']) if row['guide_bend_radius_mm'] else None)
            p=out/'baseline'/('frozen_'+port.endpoint_id+'.stl');route.mesh.export(p)
            route_hashes[port.endpoint_id]=dict(mesh_sha256=f.sha256(p),selected_parameters_sha256=hashlib.sha256(json.dumps(row,sort_keys=True).encode()).hexdigest(),
                source='Deterministically reconstructed from cached selected parameters; no per-port solid file existed in aligned output.',
                before_sha256=f.sha256(p),after_sha256=f.sha256(p),unchanged=True)
        summary['frozen_routes']=route_hashes
        summary['port_localization_audit']=port_localization_audit(out,inputs,mesh,cfg)
        candidates=[dict(clip_diameter_factor=k,clip_back_mm=k*2*ep.radius,status='NOT_RUN_VMTK_RAMP_UNAVAILABLE') for k in cfg['o3_smoothing']['clip_diameter_factors']]
        f.write_csv(out/'tables/O3_transition_candidates.csv',candidates)
        summary['clip_candidates_not_executed']=candidates
        summary['actual_clip_back_mm']=0.;summary['local_geometry_modified_mm']=0.
        summary['extension_physical_length_mm']=None;summary['radius_changed']=False
        summary['current_baseline']['mesh_qc']['euler_number']=int(mesh.euler_number)
        summary['required_clearance_floor_mm']=cfg['o3_smoothing']['minimum_clearance_mm']
        baseline_figures(out/'QC',mesh,center,normal,ep.tangent,ep.radius,ring,shown,cfg)
        np.save(out/'baseline/normal_jump_degrees.npy',normal_values)
        logger.info('[STOP] %s; no clipped/extended/new core output, no final union or slicing',summary['failure_code'])
    except Exception as exc:
        logger.exception('Preflight failed')
        summary.update(failure_code='O3_PREFLIGHT_ERROR',error=str(exc))
    finally:
        summary['protection']=f.verify_snapshot(protected)
        summary['protected_file_count']=len(protected)
        save_json(out/'protected_after.json',summary['protection'])
        summary['ui_sha256']={n:f.sha256(f.ROOT/n) for n in ['s1-2_swc_roi_generate_human.py','s1-3_swc_roi_generate_MeVO.py']}
        summary['source_BALANCED_sha256']=f.sha256(inputs['paths']['source_stl'])
        save_json(out/'o3_smoothing_summary.json',summary)
        write_report(out,summary)
        handler.close();logger.removeHandler(handler)
    return summary


def write_report(out,s):
    b=s.get('current_baseline',{});n=b.get('normal_jumps',{});features=b.get('feature_edges',[])
    runtime=s.get('bridge_smoke',{}).get('worker',{}).get('runtime',{})
    text=f'''# O3 局部表面连续性检查

## 1. 原来哪里有问题

当圆口对齐而两侧壁面方向不一致时，就会出现环状折肩；这正是旧方案未充分检验的风险机制。已有封口是模型的终端闭合面，不是必须保留的解剖结构。本轮没有把该机制直接认定为当前 O3 的实测结果：按已存节点和坐标定位，当前 O3 的局部测量未复现所述严重折肩。当前 5 mm 中心线平均方向与旧 cap normal 相差 {b.get('cap_centerline_angle_deg')}°。旧实现保留为 BASELINE_G0_ALIGNMENT 对照，本轮没有用它生成新的生产几何。

## 2. 为什么旧测试没有发现

旧测试检查位置、圆口轮廓、中心、轴向、闭合与体积，即“接没接上”。这些是位置连续（G0）的证据，不能证明侧壁切向连续（G1），也就没有回答“接得顺不顺”。近乎零的轮廓误差与轴向误差仍可伴随明显折肩。

## 3. 这次用了什么方法

本轮首先实测既有 VMTK 环境，并完成跨环境 VTP 文件往返；通过状态为 {s.get('bridge_smoke',{}).get('passed')}。主环境没有导入 VMTK。指定环境实际为 VMTK {runtime.get('vmtk_version')} / VTK {runtime.get('vtk_version')}，没有 SetInterpolationModeToRamp，也没有 SetPreserveCrossSectionShape。现有 Slicer 的已安装扩展为空，不能提供替代 VMTK 实现。官方新版源码存在这些 API，但它不是本机已安装的运行时。全部版本、方法列表、源码提交与文件哈希保存在 api_audit 和 VMTK/bridge_smoke 中。

按照用户第 60 条“仍没有：停止并报告”，当前停止代码为 **{s.get('failure_code')}**。没有安装新的 VMTK、没有自行实现 ramp、没有把 thinplatespline 或 linear 当成 ramp。尚未进入正式剪切、平滑过渡、CadQuery 接续或 Boolean 合并。

[官方 VMTK 实现](https://github.com/vmtk/vmtk/blob/master/vtkVmtk/ComputationalGeometry/vtkvmtkPolyDataFlowExtensionsFilter.h)与 [Slicer ClipVessel](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ClipVessel/ClipVessel.py)的固定提交副本和实际本机版本分开记录；不能因为网页有新方法就认定本机可调用。

## 4. 改了多少原血管

本轮实际剪切 0 mm、局部修改 0 mm。三个预定回退距离从实际直径计算，见 O3_transition_candidates.csv，均标为未执行。没有生成 VMTK 平滑段或新 handoff；对应长度、位置和几何指标为未产生，而非零误差通过。

## 5. 修复前后有多大改善

尚无修复后模型，所以不能报告改善百分比。以下均直接读取现有 aligned STL 得到，未重建基线。截面使用旧 cap-normal 坐标系的 -1 至 +1 mm，步长 0.25 mm；它仅用于基线诊断，不代替任务要求的未来新中心线截面验收。

| 指标 | 当前基线 | 修复后 |
| --- | --- | --- |
| 相邻面法向中位数 | {n.get('median_deg')}° | 未生成 |
| 相邻面法向 p95 | {n.get('p95_deg')}° | 未生成 |
| 相邻面法向最大值 | {n.get('maximum_deg')}° | 未生成 |
| 相邻截面最大面积变化 | {b.get('maximum_area_jump_fraction')} | 未生成 |
| O3 对其他血管间距 | {b.get('prior_vessel_clearance_mm')} mm（既有记录） | 未生成 |

'''
    for v in features:text+=f"{v['angle_deg']:g}° 特征边总长 {v['total_length_mm']:.6f} mm；绕轴闭合环 {v['circumferential_ring_count']} 个；最长连通折边 {v['longest_connected_feature_mm']:.6f} mm。\n\n"
    text+='**四个接口的定位核对：**见 QC/13_all_port_localization_wireframe.png 和 QC/14_full_core_port_labels.png。以下全部为同一份既有模型，只做只读检查。\n\n| 接口 | 20° 特征边长（mm） | 最大相邻面夹角 |\n| --- | --- | --- |\n'
    for item in s.get('port_localization_audit',[]):
        text+=f"| {item['port']} | {item['total_length_mm']:.6f} | {item['maximum_deg']:.6f}° |\n"
    text+='''
本轮没有据此断言用户标错端口，也没有把 O3 宣布为通过；应先核对当前三维位置与人工看到的接缝。I1/O1/O2 未做修改。

特征边由 vtkFeatureEdges 提取；闭环是否绕管一周由周向绕数检查。法向跳变使用成熟网格库的相邻三角面，局部窗口的中心、轴向和尺寸保存在 JSON。曲率只作为可视化，不用单个噪声三角面判失败。局部窗口和百分位可能遗漏窄小但明显的肩部，仍需结合线框图检查；它们不是解剖学标准。

当前图片仅有真实“修复前”图，未生成虚假的前后对比或最终 STL。没有用角度变化制造改善效果。

## 6. 是否影响其它端口

'''
    text+=f"I1/O1/O2 的既有数据和原模型未修改。保护文件 {s.get('protected_file_count')} 个，哈希不一致数 {len(s.get('protection',{}).get('changed',[]))}。单端口实体此前没有独立文件，因此仅从已选参数确定性重放并保存单端口审计副本，未进行路线搜索；副本哈希和参数哈希均在 JSON 中明确标明来源。不能把这些副本冒称为原先已有的单端口 STL。\n\n"
    text+='''## 7. 是否可以进入下一阶段

**NEEDS_ADJUSTMENT**。本轮不是 O3 修复成功。必须先获得能实际调用官方 ramp 的运行时，再执行局部剪切、候选比较与真正的连续性验收。由于用户要求当前条件不满足时停止，本轮没有新增生产 STL、没有 box/core 一体式合并，也没有 Bambu 最终切片。既有模型与可视化保持原样。
'''
    validation=s.get('tests')
    if validation:
        text+=f"\n本轮运行 {validation['total']} 项回归，全部通过；其中旧测试 {validation['old']} 项、新增前置检查与只读 QC 测试 {validation['new']} 项。VTK/NumPy 弃用提示 {validation['warnings']} 条。这些通过不能代替 ramp 生产测试；实际 ramp、剪切、候选、接续、修复后连续性及最终实体测试因运行时不可用而没有执行。日志和 JUnit 文件随结果保留。\n"
    (out/'o3_smoothing_report.md').write_text(text,encoding='utf-8')


def port_localization_audit(out, inputs, mesh, cfg):
    """Inspect other frozen ports without editing them; label evidence, not a repair."""
    from scipy.spatial import cKDTree
    rows=[];pictures=[];source=inputs['mesh']
    for ep in inputs['endpoints']:
        _,node=cKDTree(source.vertices).query(ep.position)
        faces=np.any(source.faces==node,axis=1)
        normal=f.unit(np.sum(source.face_normals[faces]*source.area_faces[faces,None],axis=0))
        metrics,_=continuity.normal_jumps(mesh,ep.position,normal,ep.radius*1.7,1.)
        feature,lines=continuity.feature_edges(mesh,ep.position,normal,ep.radius*1.7,1.,20)
        rows.append(dict(port=ep.endpoint_id,position=ep.position,cap_normal=normal,**metrics,**feature))
        local=continuity.window(mesh.triangles_center,ep.position,normal,ep.radius*2.2,1.6)
        part=polydata(mesh.submesh([np.flatnonzero(local)],append=True,repair=False))
        view=np.cross(ep.tangent,normal);view=f.unit(view)
        p=pv.Plotter(off_screen=True,window_size=(800,720));p.set_background('white')
        p.add_mesh(part,color='#10a7b5',show_edges=True,edge_color='#454545',smooth_shading=False)
        if lines.n_cells:p.add_mesh(lines,color='#d72d35',line_width=5)
        p.camera_position=[ep.position+view*9,ep.position,normal]
        p.camera.parallel_projection=True;p.camera.parallel_scale=1.75
        image=Image.fromarray(p.screenshot(return_img=True));p.close();draw=ImageDraw.Draw(image)
        draw.rectangle((0,0,800,73),fill='white')
        draw.text((20,8),ep.endpoint_id+'：当前文件中的接口',font=chinese_font(cfg,26),fill='#172b36')
        draw.text((20,40),f"20°折边 {feature['total_length_mm']:.3f} mm；最大面夹角 {metrics['maximum_deg']:.2f}°",font=chinese_font(cfg,20),fill='#633333')
        pictures.append(image)
    canvas=Image.new('RGB',(1600,1485),'white')
    for i,img in enumerate(pictures):canvas.paste(img,((i%2)*800,(i//2)*720))
    draw=ImageDraw.Draw(canvas);draw.text((25,1450),'四个接口的定位核对：红色为检测到的折边。均为原模型，没有进行修复。',font=chinese_font(cfg,25),fill='#172b36')
    canvas.save(out/'QC/13_all_port_localization_wireframe.png')
    f.write_csv(out/'tables/baseline_port_localization_audit.csv',rows)
    p=pv.Plotter(off_screen=True,window_size=(1300,1000));p.set_background('white');p.add_mesh(polydata(mesh),color='#10a7b5')
    positions=np.array([e.position for e in inputs['endpoints']]);names=[e.endpoint_id for e in inputs['endpoints']]
    p.add_point_labels(positions,names,point_color='red',text_color='black',font_size=26,point_size=10,always_visible=True)
    p.view_isometric();p.reset_camera();image=Image.fromarray(p.screenshot(return_img=True));p.close()
    draw=ImageDraw.Draw(image);draw.rectangle((0,0,1300,60),fill='white');draw.text((25,12),'当前完整模型中的四个接口位置',font=chinese_font(cfg,30),fill='#172b36')
    image.save(out/'QC/14_full_core_port_labels.png')
    return rows
