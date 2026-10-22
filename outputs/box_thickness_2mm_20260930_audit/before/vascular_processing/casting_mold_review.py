"""Static visual evidence and a nine-section Chinese casting-mold report."""
from pathlib import Path
import json
import numpy as np
from PIL import Image
import pyvista as pv
import trimesh
from .sacrificial_fixture_review import polydata,chinese_font,save_json
from .print_frame_review import caption

NAMES=['01_accepted_vascular_core','02_intact_five_wall_box','03_core_inside_box_before_union',
 '04_port_wall_overlap','05_final_one_body_abs_mold','06_four_wall_cross_sections','07_pdms_volume_preview',
 '08_pdms_ligament_overview','09_O3_pdms_gap','10_core_to_wall_clearance','11_internal_support_risk',
 '12_top_support_access','13_orientation_top3','14_selected_print_orientation','15_bambu_slice_preview']


def scene(core,box=None,bounds=None,view=None,opacity=.2,crop=None,highlight=None,lines=None,points=None,
          solid=None,bed=False,scale=None,support=None):
    plot=pv.Plotter(off_screen=True,window_size=(1200,850));plot.set_background('#f5f8fb')
    def surface(m):
        data=polydata(m) if isinstance(m,trimesh.Trimesh) else m
        return data.clip_box(bounds=np.asarray(crop).T.ravel(),invert=False) if crop is not None else data
    if core is not None:plot.add_mesh(surface(core),color='#0797a7')
    if box is not None:plot.add_mesh(surface(box),color='#91a1b3',opacity=opacity,show_edges=True,edge_color='#6b7e8d')
    if solid is not None:plot.add_mesh(surface(solid),color='#c394df',opacity=.18)
    if highlight is not None and len(highlight.faces):plot.add_mesh(surface(highlight),color='#f29227')
    if lines:
        for p,q,label in lines:
            plot.add_mesh(pv.Line(p,q),color='#d33b46',line_width=5)
            plot.add_point_labels(np.array([(np.asarray(p)+q)/2]),[label],font_size=20,text_color='#ba2432',point_size=0,always_visible=True,shape_opacity=.8)
    if points:
        for point,passed in points:plot.add_mesh(pv.Sphere(radius=.8,center=point),color='#32ac63' if passed else '#e44848')
    if support is not None and len(support):
        selected=support[::max(1,int(np.ceil(len(support)/80000)))];n=len(selected)
        plot.add_mesh(pv.PolyData(selected.reshape(-1,3),lines=np.column_stack([np.full(n,2),np.arange(2*n).reshape(-1,2)]).ravel()),color='#e79722',line_width=1.1)
    if bounds is None:
        bounds=core.bounds if core is not None else box.bounds
    bounds=np.array(bounds);center=bounds.mean(0);span=float(np.max(bounds[1]-bounds[0]))
    if bed:plot.add_mesh(pv.Plane(center=[center[0],center[1],0],direction=[0,0,1],i_size=span*1.15,j_size=span*1.15),color='#dce3e9',opacity=.35)
    view=np.array([1.6,-2.3,2.6] if view is None else view,dtype=float)
    up=[0,1,0] if abs(view[2])/np.linalg.norm(view)>.99 else [0,0,1]
    plot.camera_position=[center+view*span,center,up];plot.enable_parallel_projection();plot.camera.parallel_scale=scale or span*.65
    plot.add_axes();image=plot.screenshot(return_img=True);plot.close();return image


def montage(panels,columns):
    width=panels[0].width;height=max(p.height for p in panels);image=Image.new('RGB',(columns*width,int(np.ceil(len(panels)/columns))*height),'white')
    for k,p in enumerate(panels):image.paste(p,((k%columns)*width,(k//columns)*height))
    return image


def printing_view(row):
    # View from the actual opening side while remaining above the build plate.
    opening=np.asarray(row['transform_4x4'])[:3,2]
    return 3*opening+np.array([.6,-.6,2.])


def wall_cross_sections(path,inputs,union,box,cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import Rectangle
    font=FontProperties(fname=cfg['visualization']['chinese_font'])
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for ax,port in zip(axes.ravel(),inputs['ports']):
        axis='XY'.index(port['face'][-1]);other=1-axis;origin=port['target'];normal=np.zeros(3);normal[other]=1
        side=0 if port['face'][0]=='-' else 1;edge=box['inner'][side,axis];wall_edge=box['outer'][side,axis];z=origin[2]
        ax.add_patch(Rectangle((min(edge,wall_edge),z-3),abs(wall_edge-edge),6,color='#bdc8d2',alpha=.7))
        for mesh,color,linewidth in [(inputs['core'],'#e58b23',2.),(union,'#08799b',1.3)]:
            section=mesh.section(plane_origin=origin,plane_normal=normal)
            if section is not None:
                for loop in section.discrete:ax.plot(loop[:,axis],loop[:,2],color=color,lw=linewidth)
        ax.set_xlim(min(edge,wall_edge)-3,max(edge,wall_edge)+3);ax.set_ylim(z-3,z+3);ax.set_aspect('equal')
        ax.set_title(port['port_id']+'：端口贯穿实体侧壁',fontproperties=font);ax.set_xlabel('XY 壁法向坐标 / mm',fontproperties=font);ax.set_ylabel('Z / mm');ax.grid(alpha=.18)
    fig.suptitle('灰色为墙体区域，橙色为血管截线，蓝色为最终整体的真实截线',fontproperties=font,fontsize=14)
    fig.savefig(path,dpi=150);plt.close(fig)


def render(out,inputs,box,union,contacts,partition,pdms,wall_rows,top,selected,summary,cfg):
    out=Path(out);folder=out/'QC';core=inputs['core'];bounds=union.bounds;saved=[]
    def save(index,title,caption_lines,**kwargs):
        im=scene(kwargs.pop('core',core),bounds=kwargs.pop('bounds',bounds),**kwargs);path=folder/(NAMES[index-1]+'.png')
        caption(im,title,caption_lines,cfg).save(path);saved.append(str(path))
    save(1,'已验收血管芯：保持原样',['四个接口、半径、路径和衔接均未修改。'])
    save(2,'四侧封闭、底部封闭、顶部开放',['新生成无装配孔盒体；尺寸与已审核旧盒一致。'],core=box['mesh'])
    save(3,'合并前：血管穿过完整侧壁',['灰色半透明盒体与青色已验收血管；端口穿墙处随后作实体合并。'],box=box['mesh'],opacity=.16)
    panels=[]
    for port in inputs['ports']:
        contact=contacts[port['port_id']];center=contact.bounds.mean(0);local=np.array([center-5,center+5])
        panels.append(caption(scene(core,box['mesh'],bounds=local,crop=local,opacity=.16,highlight=contact),
            port['port_id']+'：端口与侧壁实体重叠',[f"橙色为实际交叠实体，体积 {abs(contact.volume):.4f} mm³。"],cfg))
    path=folder/(NAMES[3]+'.png');montage(panels,2).save(path);saved.append(str(path))
    save(5,'一体式 ABS 牺牲倒模件',['此图显示最终 Boolean 网格：一个主颜色、一个连通实体。'],core=union)
    path=folder/(NAMES[5]+'.png');wall_cross_sections(path,inputs,union,box,cfg);saved.append(str(path))
    future=pv.read(out/'diagnostic/future_PDMS_volume.vtp')
    save(7,'未来 PDMS 块与血管通道',['紫色为内腔减去血管芯后的 PDMS 空间；青色位置在 ABS 去除后成为管腔。',
        '这是诊断预览，不是制造文件。'],solid=future)
    thin=pdms['top10'][:2];gaplines=[(r['first_point_mm'],r['second_point_mm'],f"{r['distance_mm']:.4f} mm") for r in thin]
    save(8,'PDMS 最薄位置',['红线为非连接血管表面之间的 FCL 三角面距离；已排除真实分叉连接区。'],lines=gaplines)
    pair=pdms['o3_pair'];a=np.array(pair['first_point_mm']);b=np.array(pair['second_point_mm']);local=np.array([np.minimum(a,b)-4,np.maximum(a,b)+4])
    save(9,'O3 附近最薄 PDMS：保持验收间距',[f"最薄距离 {pair['distance_mm']:.6f} mm；相对验收值变化 {pdms['o3_deviation_mm']:.6f} mm。"],bounds=local,crop=local,lines=[(a,b,f"{pair['distance_mm']:.4f} mm")])
    row=next(r for r in wall_rows if not r['legal_port_crossing_excluded'])
    save(10,'血管到内壁的 PDMS 厚度',[f"除合法穿墙带外，最近位置 {wall_rows[0]['distance_mm']:.2f} mm，位于排除带边界。",
        f"图示不涉及穿墙排除带的位置：{row['distance_mm']:.3f} mm。"],box=box['mesh'],lines=[(row['vascular_point_mm'],row['wall_point_mm'],f"{row['distance_mm']:.3f} mm")])
    risk=union.submesh([selected['support']['risk_face_ids']],append=True,repair=False)
    save(11,'选定打印姿态下的内部支撑风险',['橙色为按真实工艺阈值识别的悬垂表面；图像以倒模姿态显示。',
        '包含血管及盒内壁；橙色表面积不等于真实支撑体积。'],box=box['mesh'],highlight=risk)
    rows=selected['support']['rows'];points=[(r['representative_surface_mm'],r['status']!='TRAPPED_SUPPORT_RISK') for r in rows if not r['vascular_branch'].startswith('BOX')]
    save(12,'只能从顶部进入的支撑清理通道',['俯视：绿色为找到顶部通路的采样邻域，红色为未找到通路。',
        '3 mm 探针检查是几何代理，不保证实际支撑可取出。'],box=box['mesh'],view=[0,0,3],points=points)
    panels=[];scale=max(max(r['bbox_extents_mm']) for r in top[:3])*.65
    for row in top[:3]:
        m=union.copy();m.apply_transform(row['transform_4x4'])
        panels.append(caption(scene(m,bounds=[row['bbox_min_mm'],row['bbox_max_mm']],bed=True,scale=scale,view=printing_view(row)),f"打印候选 {row['rank']}",
            [f"内部风险表面积 {row['internal_support_risk_area_mm2']:.1f} mm²；阻挡风险 {row['trapped_risk_count']}",f"高度 {row['z_height_mm']:.2f} mm；顶部开口随模型一起旋转。"],cfg))
    path=folder/(NAMES[12]+'.png');montage(panels,len(panels)).save(path);saved.append(str(path))
    oriented=union.copy();oriented.apply_transform(selected['transform_4x4'])
    save(14,'选定打印姿态',['+X 侧壁贴打印板，开口朝侧方；打印后须转回开口朝上，再倒入 PDMS。'],core=oriented,bounds=oriented.bounds,bed=True,view=printing_view(selected))
    actual=summary.get('bambu',{}).get('selected');path=folder/(NAMES[14]+'.png')
    if actual and actual.get('native_preview'):
        im=Image.open(actual['native_preview']).convert('RGB');im=im.resize((1000,1000))
        caption(im,'Bambu 原生切片文件预览',['直接读取本次 3MF 内 Bambu 生成的预览图；没有用几何渲染替代。',
            f"实际切片成功；时间 {actual['estimated_print_time_seconds']:.0f} s，材料 {actual['filament_used_g']:.2f} g。"],cfg).save(path)
        arrays=np.load(actual['actual_paths_file']);support=arrays['segments'][arrays['kind']==1]
        placed_core=core.copy();placed_core.apply_transform(selected['transform_4x4'])
        placed_box=box['mesh'].copy();placed_box.apply_transform(selected['transform_4x4'])
        caption(scene(placed_core,box=placed_box,opacity=.08,bounds=oriented.bounds,support=support), '实际 G-code 支撑轨迹（补充图）',
            ['为看清盒内支撑，本诊断图将墙体透明显示；生产墙体完整封闭。',
             '橙色来自实际挤出轨迹，按归档喷头偏移恢复模型坐标；不是独立支撑实体。'],cfg).save(folder/'16_actual_support_toolpaths.png')
    else:
        im=Image.new('RGB',(1200,850),'#f5f8fb');caption(im,'Bambu preview not available from current CLI',
            ['当前没有可用的真实 Bambu 预览，不以模拟图替代。'],cfg).save(path)
    saved.append(str(path))
    if summary['failures']:
        caption(scene(core,box['mesh'],bounds=bounds,view=[0,0,3],points=points),'当前结果需要调整',summary['failures'][:3],cfg).save(folder/'failure_evidence.png')
    return saved


def failure_image(out,inputs,box,union,error,cfg):
    mesh=union if union is not None else inputs['core']
    caption(scene(mesh),'倒模件检查未完成',[error[:90]],cfg).save(Path(out)/'QC/failure_evidence.png')


def write_report(out,s):
    out=Path(out);box=s.get('box',{});union=s.get('union',{});pd=s.get('pdms',{});a=s.get('support',{});b=s.get('bambu',{});o=s.get('orientation',{});selected=o.get('selected',{});actual=b.get('selected') or {}
    dims=lambda x:' × '.join(f'{n:.2f}' for n in x) if x is not None else '未完成'
    parts=['# BG001 RMCA 一体式 ABS 牺牲倒模件评估\n']
    def add(n,title,text):parts.append(f'## {n}. {title}\n\n{text}\n')
    add(1,'这次做了什么','将已验收的血管芯与新生成的完整盒体作实体合并。盒体四侧及底部封闭，顶部完整开放；没有装配间隙孔、顶盖或顶部横梁。旧开放框架作为历史探索保留，本轮生产目标改为 OPEN_TOP_FIVE_WALL_CASTING_BOX。没有重新生成血管、调整半径、重做端口路径或接口光顺。')
    add(2,'最终是不是一个 ABS 整体',f"{'是' if union.get('connected_components')==1 and union.get('passed') else '尚未确认'}。连通分量 {union.get('connected_components')}；封闭检查 {union.get('watertight')}，流形检查 {union.get('manifold')}，退化三角形 {union.get('degenerate_faces')}。合并后的血管材料损失 {union.get('core_material_loss_mm3')} mm³，盒体材料损失 {union.get('box_material_loss_mm3')} mm³。原始血管 SHA256 为 `{s['input']['accepted_core_hash']}`。")
    add(3,'盒子结构',f"BOX_DIMENSIONS_REUSED：沿用已审核尺寸，本轮没有重新优化或扩大盒体。\n\n| 项目 | 数值 |\n|---|---|\n| 内部尺寸 | {dims(box.get('inner_dimensions_mm'))} mm |\n| 外部尺寸 | {dims(box.get('outer_dimensions_mm'))} mm |\n| 侧壁厚度 | {box.get('wall_thickness_mm')} mm |\n| 底厚 | {box.get('bottom_thickness_mm')} mm |\n| 顶部开口 | {box.get('top_opening_width_mm',0):.2f} × {box.get('top_opening_length_mm',0):.2f} mm |\n| 开口面积 | {box.get('top_opening_area_mm2',0):.2f} mm² |\n\n上行射线与内腔实体交集检查确认顶部无封板。单独盒体 STL 和 STEP 位于 reference/，最终制造件位于 final/。")
    table='| 端口 | 原固定侧壁 | 重叠体积 (mm³) | 轴向交叠 (mm) | 实体连接 |\n|---|---|---:|---:|---|\n'
    for r in s.get('ports',[]):table+=f"| {r['port_id']} | {r['wall']} | {r['intersection_volume_mm3']:.6f} | {r['intersection_axial_span_mm']:.3f} | {r['status']} |\n"
    add(4,'四个端口',table+'\n要求交叠长度至少达到壁厚的 90%，并且交叠体积为正。端口在完整侧壁中形成实体连接，没有悬在装配孔里。四个血管衔接口的局部几何另存于 tables/frozen_interface_qc.csv。')
    pair=pd.get('minimum_pair',{})
    add(5,'PDMS 最薄区域',f"全局最薄非连接区域位于 {pair.get('first')} 与 {pair.get('second')} 之间，距离 {pd.get('minimum_ligament_mm',0):.6f} mm，几何分级为 {pair.get('geometry_grade')}。O3 间距 {pd.get('o3_gap_mm',0):.6f} mm，相对验收值变化 {pd.get('o3_deviation_mm',0):.6f} mm。距离由 FCL 对最终交付三角面计算，排除同分支及真实分叉、端口衔接的局部连接区；没有把连接处的零距离当作 PDMS 薄区。\n\n对内壁厚度，四个合法穿墙区域先排除 {pd.get('wall_crossing_exemption_mm')} mm 的局部轴向带。带外最小值 {pd.get('minimum_core_wall_mm')} mm 位于排除带边缘，属于该定义的边界值；中央原生血管到壁面的最小厚度为 {pd.get('minimum_native_core_wall_mm')} mm。未来 PDMS 占据内腔减去血管芯的空间，诊断体积 {pd.get('future_volume_mm3',0)/1000:.2f} mL。诊断 VTP 只用于理解外形和管腔关系。")
    quantity=actual.get('support_quantity',{});baseline=s.get('baseline',{}).get('support',{})
    add(6,'内部打印支撑',f"需要内部支撑。竖直基准先进行了关闭支撑、开启支撑的实际切片；其顶部通路代理识别到 {baseline.get('trapped_risk_count')} 个阻挡风险邻域，因此继续搜索姿态。选定姿态中，{a.get('candidate_regions')} 个潜在支撑邻域有 {a.get('top_accessible_count')} 个找到顶部通路，其中 {a.get('probe_accessible_count')} 个通过直径 {a.get('probe_diameter_mm')} mm 探针检查；剩余阻挡风险 {a.get('trapped_risk_count')} 个。顶部通路始终在倒模坐标下计算，侧壁从未被当作移除出口。\n\n实际支撑类型取自当前 Bambu 工艺，为 {actual.get('actual_project_settings',{}).get('support_type','未确认')}。根据真实挤出路径在内腔中的分布，盒内支撑估计使用耗材 {quantity.get('estimated_internal_support_mass_g',0):.2f} g、对应挤出体积 {quantity.get('estimated_internal_support_volume_mm3',0):.2f} mm³。这是按 G-code E 值和实际耗材直径、密度计算的估计，不是独立支撑实体体积。\n\nProbe accessibility is a geometry proxy, not a guarantee that real support can be removed. 当前没有独立支撑网格，状态为 BAMBU_SUPPORT_GEOMETRY_UNAVAILABLE；实际拆除仍需要人工审核和实物试验。")
    add(7,'Bambu 实际切片与打印姿态',f"实际状态：{b.get('status')}。版本 {b.get('version')}，打印机 {b.get('printer')}，喷嘴 {b.get('nozzle_mm')} mm，ABS 配置 {b.get('filament')}，工艺 {b.get('process')}。选定切片预计耗时 {b.get('print_time_seconds')} s、总耗材 {b.get('filament_usage_g')} g。\n\n共检查 {o.get('candidate_count')} 个确定性姿态，采用显式字典序比较阻挡风险、内部支撑风险面积、顶部可达性、外部支撑、总支撑量、高度和占地；真实切片后用实际支撑挤出量细化同等几何候选的排序。最终姿态 {selected.get('method')}，整体包围盒 {dims(selected.get('bbox_extents_mm'))} mm，高度 {selected.get('z_height_mm',0):.2f} mm。打印后按照 casting_restore_transform.json 恢复顶部朝上再灌注。\n\n四个端口、主要血管分支、四壁与底板均使用真实非支撑挤出路径作采样保留检查，结果见 tables/slicer_feature_preservation.csv；3MF 内部网格也单独核对刚性变换。真实 Bambu 原生预览保存为 QC/15_bambu_slice_preview.png，补充支撑路径图来自实际 G-code。没有发送、上传或启动打印。工艺阈值来源及未输出文本的真实 --help 探测均保留在 bambu/。")
    comparison='\n\n实际开启支撑的切片比较（关闭支撑的结果同样保存在 CSV 中）：\n\n| 候选 | 打印分钟 | 总耗材 g | 盒内支撑估计 g | 特征保留检查 |\n|---|---:|---:|---:|---|\n'
    for record in b.get('results',[]):
        if record.get('support_mode')!='ON' or record.get('status')!='BAMBU_SLICED':continue
        comparison+=f"| {record['orientation_id']} | {record['estimated_print_time_seconds']/60:.1f} | {record['filament_used_g']:.2f} | {record['support_quantity']['estimated_internal_support_mass_g']:.2f} | {'通过' if record['feature_preservation']['passed'] else '失败'} |\n"
    if actual:
        frame=actual.get('gcode_coordinate_frame',{})
        comparison+=f"\n选定候选为 {selected.get('candidate_id')}。+X 侧壁贴打印板，开口在打印时朝侧方；脱离打印板后，将整个模具转回开口朝上，底板在下，再进行灌注。正向与逆向矩阵分别保存在 print_transform_final.json 和 casting_restore_transform.json，均不含缩放。\n\n路径检查采用归档配置中的喷头偏移 {frame.get('offset_added_mm')} mm 恢复模型坐标；未改变 G-code 或模型。坐标规则依据 [Bambu GCode 实现](https://github.com/bambulab/BambuStudio/blob/master/src/libslic3r/GCode.cpp)。启动和擦拭路径、支撑路径均不计入血管保留检查。该检查以 3 mm 空间分箱、0.60 mm 最近路径容差验证结构存在，并非完整的挤出珠体积重建。\n"
    off_warnings=sum(bool(r.get('unsupported_floating_thin_warning_lines')) for r in b.get('results',[]) if r.get('support_mode')=='OFF')
    comparison+=f"\n关闭支撑的切片中有 {off_warnings} 个候选收到 Bambu 的悬空结构提示，不能仅凭导出成功就认定可打印。选定的开启支撑切片有 {len(actual.get('unsupported_floating_thin_warning_lines',[]))} 条悬空、无支撑或薄结构警告。\n"
    if s.get('coordinate_audit_correction'):
        comparison+='\n初次检查因未恢复喷头偏移而误报缺失；初次结果保留在 diagnostic/initial_coordinate_audit/。校正坐标后复查了相同的 10 份真实切片，归档 SHA256 未变。\n'
    parts[-1]+=comparison
    add(8,'当前最重要的风险','1. 顶部可达性仅是几何代理；实际树状支撑可能需要分段切除，必须确认不会牵拉或折断细血管。\n2. O3 附近约 1.69 mm 的 PDMS 薄区虽然通过既定回归门槛，仍需实物验证。\n3. 侧放显著增加盒内壁支撑。当前按“清理通路优先”的既定顺序选定姿态，并不是最省料的姿态。\n4. 尚未验证 ABS 实际打印、PDMS 灌注、丙酮去除和流动实验。')
    add(9,'最终状态',f"**{s['status']}**。这是供人工审阅的结果，不是制造验证完成。\n\n受保护文件 {s.get('protected_file_count')} 个，哈希不匹配 {len(s.get('hash_mismatches',[]))} 个。硬失败项 {len(s.get('failures',[]))} 个。尚未取得独立支撑实体，部分非选定切片的首层耗时元数据异常已记录但未用于决策；完整警告见 JSON。永久测试在 tests/test_s1_6_abs_casting_mold.py，测试结果保存在 test_results.xml。生成入口为 s1-6_abs_casting_mold.py；复运行时使用 --output-root 指向新目录，以保护所有历史结果。")
    if s.get('tests'):
        t=s['tests'];parts[-1]+=f"\n实际测试：{t['passed']} 项通过，{t['failed']} 项失败；其中包含小型合成模型与真实 BG001 结果回归。验证记录见 validation_results.json。\n"
    (out/'abs_casting_mold_report.md').write_text('\n'.join(parts),encoding='utf-8')
