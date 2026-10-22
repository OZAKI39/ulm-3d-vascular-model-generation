"""Readable Chinese evidence figures for the independent, frozen-core frame stage."""
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
import pyvista as pv
import trimesh
from .sacrificial_fixture_review import polydata, chinese_font
from .sacrificial_print_frame import box_mesh

NAMES=['01_accepted_vascular_core','02_old_box_reference','03_OPEN_PANEL_candidate',
 '04_SPARSE_FRAME_candidate','05_two_open_faces','06_port_frame_connections',
 '07_pdms_keep_zone_and_frame','08_frame_intrusion_detail','09_support_access_view_A',
 '10_support_access_view_B','11_probe_accessibility','12_candidate_comparison',
 '13_orientation_top3','14_selected_orientation','15_bambu_slice_preview']


def caption(image,title,lines,cfg):
    image=Image.fromarray(image).convert('RGB') if isinstance(image,np.ndarray) else image.convert('RGB')
    header=76;footer=42*max(1,len(lines))+18
    canvas=Image.new('RGB',(image.width,image.height+header+footer),'white');canvas.paste(image,(0,header))
    draw=ImageDraw.Draw(canvas);draw.text((24,14),title,font=chinese_font(cfg,32),fill='#123048')
    for k,line in enumerate(lines):draw.text((24,header+image.height+10+42*k),line,font=chinese_font(cfg,22),fill='#253a48')
    return canvas


def scene(core,frame=None,bounds=None,view=None,keep=None,ports=None,markers=None,paths=None,
          contacts=None,bed=False,support=None,frame_opacity=.72,core_opacity=1.,crop=None):
    plot=pv.Plotter(off_screen=True,window_size=(1200,830));plot.set_background('#f5f8fb')
    def display(mesh):
        data=polydata(mesh)
        return data.clip_box(bounds=np.asarray(crop).T.ravel(),invert=False) if crop is not None else data
    if core is not None:plot.add_mesh(display(core),color='#0797a7',opacity=core_opacity,smooth_shading=False)
    if frame is not None:plot.add_mesh(display(frame),color='#8097ad',opacity=frame_opacity,show_edges=True,edge_color='#53697b')
    if keep is not None:
        region=polydata(box_mesh(keep));plot.add_mesh(region,color='#b67bd2',opacity=.13)
        plot.add_mesh(region.outline(),color='#8953a0',line_width=2)
    if contacts is not None and len(contacts.faces):plot.add_mesh(display(contacts),color='#ed8b26')
    if ports:
        points=np.array([(p['target']+p['end'])/2 for p in ports])
        plot.add_point_labels(points,[p['port_id'] for p in ports],point_color='#e87d27',point_size=8,
            text_color='#112d42',font_size=22,shape_opacity=.8,always_visible=True)
    if markers:
        for position,passed in markers:
            plot.add_mesh(pv.Sphere(radius=.9,center=position),color='#28ae51' if passed else '#ed4747')
    if paths:
        from .support_access_qc import capsule
        for p in paths:
            if np.linalg.norm(np.asarray(p['end'])-p['start'])<1e-8:continue
            plot.add_mesh(polydata(capsule(np.array(p['start']),np.array(p['end']),2.,16)),
                color='#2ca95b' if p['probe_accessible'] else '#de4343',opacity=.23)
            plot.add_mesh(pv.Line(p['start'],p['end']),color='#248648' if p['probe_accessible'] else '#b22e2e',line_width=3)
    if support is not None and len(support):
        stride=max(1,int(np.ceil(len(support)/60000)));segments=support[::stride]
        lines=np.column_stack([np.full(len(segments),2),np.arange(len(segments)*2).reshape(-1,2)]).ravel()
        plot.add_mesh(pv.PolyData(segments.reshape(-1,3),lines=lines),color='#eb9d32',line_width=1.3)
    if bounds is None:bounds=core.bounds if frame is None else np.array([np.minimum(core.bounds[0],frame.bounds[0]),np.maximum(core.bounds[1],frame.bounds[1])])
    center=np.mean(bounds,axis=0);span=float(np.max(np.diff(bounds,axis=0)))
    if bed:
        plot.add_mesh(pv.Plane(center=[center[0],center[1],-.1],direction=[0,0,1],i_size=span*1.1,j_size=span*1.1),color='#dfe5eb',opacity=.35)
    direction=np.array([1.6,-2.3,1.55]) if view is None else np.array(view,dtype=float)
    plot.camera_position=[center+direction*span,center,[0,0,1]]
    plot.enable_parallel_projection();plot.camera.parallel_scale=span*.64
    plot.add_axes(line_width=2);image=plot.screenshot(return_img=True);plot.close();return image


def render_review(out,inputs,chosen,keep,regions,top,cfg,summary):
    enrich_summary(summary)
    out=Path(out);folder=out/'QC';folder.mkdir(exist_ok=True)
    selected=chosen['SPARSE_FRAME'];core=inputs['core'];f=selected['mesh']
    bounds=np.array([np.minimum(core.bounds[0],f.bounds[0]),np.maximum(core.bounds[1],f.bounds[1])])
    saved=[]
    def save(index,title,lines,**kwargs):
        image=scene(kwargs.pop('core',core),bounds=kwargs.pop('bounds',bounds),**kwargs)
        path=folder/(NAMES[index-1]+'.png');caption(image,title,lines,cfg).save(path);saved.append(str(path))
    save(1,'最终血管芯：本轮保持原样',['青色为已验收血管和四根外接端口；文件按字节复制。'],ports=inputs['ports'])
    save(2,'旧盒体仅作对照',['旧盒体没有并入本轮生产模型。'],frame=trimesh.load_mesh(inputs['old_box_path']),frame_opacity=.28)
    for index,kind,title in [(3,'OPEN_PANEL','对照方案：移除两个相对面板'),(4,'SPARSE_FRAME','推荐方案：仅保留外围框架和端口连接')]:
        d=chosen[kind];q=d['qc']
        save(index,title,[f"框架体积 {q['frame_volume_mm3']:.2f} mm³；中央保留区侵入 {q['frame_intrusion_keep_zone_mm3']:.6f} mm³。"],frame=d['mesh'],ports=inputs['ports'])
    axis=selected['open_axis'];a=selected['qc']['open_faces']
    faceview=[0,0,0];faceview[axis]=-3
    save(5,'两个相对侧面保持开放',[f"开放面 A={selected['open_faces'][0]}，B={selected['open_faces'][1]}；自由面积 {a[0]['free_area_fraction']:.1%} / {a[1]['free_area_fraction']:.1%}。",
        '边缘梁与短端口连接保留，中央没有封板或斜撑。'],frame=f,view=faceview)
    connection_panels=[]
    for port,row in zip(inputs['ports'],selected['qc']['port_connections']):
        portaxis='XYZ'.index(port['face'][-1]);side=0 if port['face'][0]=='-' else 1
        center=port['target'].copy();center[portaxis]=(selected['inner_bounds'][side,portaxis]+selected['outer_bounds'][side,portaxis])/2
        local=np.array([center-6,center+6]);direction=port['direction']*.7+np.array([.4,-.5,1.])
        panel=caption(scene(core,f,bounds=local,crop=local,view=direction,frame_opacity=.20,core_opacity=.48,contacts=selected['contact_mesh']),
            port['port_id']+' 与框架连接（局部放大）',
            [f"橙色实体重叠 {row['intersection_volume_mm3']:.3f} mm³；接触长度 {row['intersection_length_mm']:.1f} mm。",
             '仅裁切显示窗口；生产几何没有改变。'],cfg)
        connection_panels.append(panel)
    montage=Image.new('RGB',(2400,connection_panels[0].height*2),'white')
    for k,panel in enumerate(connection_panels):montage.paste(panel,((k%2)*1200,(k//2)*panel.height))
    path=folder/(NAMES[5]+'.png');montage.save(path);saved.append(str(path))
    save(7,'框架与未来 PDMS 中央保留区',['紫色区域：冻结中央血管包围盒向六面各扩展 8 mm。'],frame=f,keep=keep)
    intrusion=selected['qc']['frame_intrusion_keep_zone_mm3']
    save(8,'框架未进入中央 PDMS 保留区' if intrusion<=1e-6 else '警示：框架进入了中央保留区',
        [f'实际交集体积：{intrusion:.9g} mm³。框架溶解后仍会在保留区外留下空腔。'],frame=f,keep=keep,contacts=selected['intrusion_mesh'])
    for index,side,letter in [(9,0,'A'),(10,1,'B')]:
        direction=[0.,0.,0.];direction[axis]=-3 if side==0 else 3
        rows=selected['access_rows'];markers=[(r['center_mm'],row['visible_from_open_face_'+letter] or row['probe_accessible_'+letter]) for r,row in zip(regions,rows)]
        count=sum(p for _,p in markers)
        save(index,'从开放侧 '+letter+' 看入内部',[f'绿色邻域可接近，红色邻域在本侧未找到通路：{count}/{len(rows)} 可接近。',
            '两图视角相对、缩放相同；这是几何代理，不代表真实支撑一定可以取出。'],frame=f,view=direction,markers=markers)
    paths=selected['access_paths'];valid=[p for p in paths if p['probe_accessible']];blocked=[p for p in paths if not p['probe_accessible']]
    examples=valid[::max(1,len(valid)//4)][:4]+blocked[:2]
    save(11,'直径 4 mm 探针的代表通路',['最多显示 6 条路径：绿色无碰撞，红色有遮挡；探针同时避开框架和血管。',
        f"{selected['access']['probe_accessible_region_count']}/{len(regions)} 邻域找到探针通路；其余邻域只有视线证据。"],frame=f,paths=examples)
    panels=[]
    for kind in ('OPEN_PANEL','SPARSE_FRAME'):
        d=chosen[kind];q=d['qc'];a=d['access']
        panels.append(caption(scene(core,d['mesh'],bounds),kind,[f"框架 {q['frame_volume_mm3']:.1f} mm³；侵入 {q['frame_intrusion_keep_zone_mm3']:.6g} mm³",
            f"开放面积 {q['open_faces'][0]['free_area_fraction']:.1%} / {q['open_faces'][1]['free_area_fraction']:.1%}",
            f"几何可接近 {a['accessible_region_count']}/{a['region_count']}；4 mm 探针 {a['probe_accessible_region_count']}/{a['region_count']}"],cfg))
    montage=Image.new('RGB',(sum(p.width for p in panels),max(p.height for p in panels)),'white');offset=0
    for p in panels:montage.paste(p,(offset,0));offset+=p.width
    path=folder/(NAMES[11]+'.png');montage.save(path);saved.append(str(path))
    if top:
        panels=[]
        for row in top[:3]:
            m=core.copy();m.apply_transform(row['transform_4x4']);fm=f.copy();fm.apply_transform(row['transform_4x4'])
            panels.append(caption(scene(m,fm,bed=True),f"候选姿态 {row['rank']}",
                [f"Z 高度 {row['z_height_mm']:.2f} mm；悬垂表面积代理 {row['downward_overhang_area_mm2']:.1f} mm²"],cfg))
        montage=Image.new('RGB',(sum(p.width for p in panels),max(p.height for p in panels)),'white');offset=0
        for p in panels:montage.paste(p,(offset,0));offset+=p.width
        path=folder/(NAMES[12]+'.png');montage.save(path);saved.append(str(path))
        row=summary['orientation']['selected'];m=core.copy();m.apply_transform(row['transform_4x4']);fm=f.copy();fm.apply_transform(row['transform_4x4'])
        save(14,'最终打印姿态',[f"仅整体刚性旋转和平移；高度 {row['z_height_mm']:.2f} mm，底面接触 {row['bed_contact_area_mm2']:.1f} mm²。"],core=m,frame=fm,bounds=np.array([row['bbox_min_mm'],row['bbox_max_mm']]),bed=True)
        actual=summary.get('bambu',{}).get('selected');segments=None
        if actual and actual.get('support_paths_file'):segments=np.load(actual['support_paths_file'])
        save(15,'Bambu 实际支撑路径预览' if actual else 'Bambu 切片尚未完成',
            ['橙色取自本次实际 G-code 支撑挤出路径（含圆弧采样）；青色血管，灰色框架。',
             '路径不等于支撑实体，尚不能证明支撑可以完整移除。'] if actual else [summary.get('bambu',{}).get('status','尚未执行切片')],
            core=m,frame=fm,bounds=np.array([row['bbox_min_mm'],row['bbox_max_mm']]),bed=True,support=segments,frame_opacity=.35)
    return saved


def enrich_summary(summary):
    """Expose requested review fields without obscuring their evidence source."""
    a=summary.get('access',{});faces=a.get('open_face_metrics',[])
    if len(faces)==2:
        a.update(open_face_A_free_fraction=faces[0]['free_area_fraction'],open_face_B_free_fraction=faces[1]['free_area_fraction'],
            minimum_clear_width_mm=min(r['minimum_clear_opening_width_mm'] for r in faces),
            minimum_clear_height_mm=min(r['minimum_clear_opening_height_mm'] for r in faces),
            probe_accessibility=a['probe_accessible_fraction'])
    o=summary.get('orientation',{});selected=o.get('selected',{})
    if selected:o.update(selected_transform=selected['transform_4x4'],bbox_mm=[selected['bbox_min_mm'],selected['bbox_max_mm']],Z_height_mm=selected['z_height_mm'])
    b=summary.get('bambu',{});d=summary.get('bambu_discovery',{});actual=b.get('selected') or {}
    b.update(binary=d.get('binary'),version=d.get('version'),printer=d.get('model'),nozzle_mm=d.get('nozzle_diameter_mm'),
        filament=d.get('active_gui_presets',{}).get('filaments'),process=d.get('active_gui_presets',{}).get('process'),
        print_time_seconds=actual.get('estimated_print_time_seconds'),filament_usage_g=actual.get('filament_used_g'))
    for mode in ('OFF','ON'):
        rows=[r for r in b.get('results',[]) if r['support_mode']==mode]
        b['support_'+mode.lower()+'_status']='BAMBU_SLICED' if rows and all(r['status']=='BAMBU_SLICED' for r in rows) else 'PENDING_OR_FAILED'
    b['verified_local_cli_options']=sorted({arg for r in b.get('results',[]) if r['status']=='BAMBU_SLICED' for arg in r['command'] if arg.startswith('--')})
    if actual:
        actual_warnings=[dict(orientation_id=r['orientation_id'],support_mode=r['support_mode'],warnings=r.get('warnings',[])) for r in b.get('results',[]) if r.get('warnings')]
        b['actual_slicer_warnings']=actual_warnings
        if actual_warnings:
            note='BAMBU_METADATA_WARNING: '+str(len(actual_warnings))+' attempts have anomalous first-layer time metadata, retained but not used for ranking or total-time estimates.'
            if note not in summary['warnings']:summary['warnings'].append(note)
    if a.get('probe_accessible_region_count',0)<a.get('region_count',0):
        note='PROBE_ACCESS_LIMITED: '+str(a['region_count']-a['probe_accessible_region_count'])+' regions have sightline evidence but no collision-free 4 mm probe path among tested paths.'
        if note not in summary['warnings']:summary['warnings'].append(note)
    summary['bambu']=b


def render_failure(out,inputs,chosen,cfg,error):
    selected=chosen.get('SPARSE_FRAME');f=selected['mesh'] if selected else None
    caption(scene(inputs['core'],f),'本次框架需要调整',[error[:85]],cfg).save(Path(out)/'QC/failure_evidence.png')


def refresh_review(out,config_path):
    """Re-read existing slice files and regenerate evidence, without re-slicing."""
    import json
    from . import sacrificial_print_frame as frame
    from . import print_frame_qc as qc
    from . import bambu_slice_adapter as bambu
    from . import sacrificial_fixture as legacy
    from .sacrificial_fixture_review import save_json
    out=Path(out);cfg=frame.load_config(config_path);inputs=frame.load_accepted(cfg)
    summary=json.loads((out/'sacrificial_print_frame_summary.json').read_text())
    protected=json.loads((out/'protected_before.json').read_text())
    protection=legacy.verify_snapshot(protected)
    if not protection['all_unchanged']:raise ValueError('FROZEN_SOURCE_HASH_MISMATCH')
    records=summary.get('bambu',{}).get('results',[])
    for record in records:
        if record.get('archive'):
            audit,segments=bambu.support_paths_from_archive(record['archive'])
            record['support_archive_audit']=audit
            np.save(record['support_paths_file'],segments)
            folder=Path(record['support_paths_file']).parent
            save_json(folder/'support_archive_audit.json',audit)
    if records:
        save_json(out/'bambu/slice_attempts.json',records)
        previous=summary['bambu'].get('selected')
        if previous:summary['bambu']['selected']=next(r for r in records if r['orientation_id']==previous['orientation_id'] and r['support_mode']==previous['support_mode'])
    chosen={};keep=np.array(summary['pdms_keep_zone']['bounds_mm'])
    for kind,data in summary['candidates'].items():
        design=dict(data)
        for key in ('outer_bounds','inner_bounds'):design[key]=np.array(design[key])
        for member in design['members']:member['bounds']=np.array(member['bounds'])
        design['mesh']=trimesh.load_mesh(design['frame_file'])
        _,design['intrusion_mesh'],design['contact_mesh']=qc.frame_checks(design,inputs,keep,cfg)
        chosen[kind]=design
    regions=json.loads((out/'vascular_regions.json').read_text())
    summary['png_files']=render_review(out,inputs,chosen,keep,regions,summary['orientation']['top_five'],cfg,summary)
    summary['review_refresh']='Existing local G-code parsed including G2/G3 arcs; geometry and slicer files not changed.'
    save_json(out/'sacrificial_print_frame_summary.json',summary)
    save_json(out/'protected_source_hashes.json',legacy.verify_snapshot(protected))
    write_report(out,summary)
    with (out/'run.log').open('a') as log:log.write('Review refreshed from actual G-code, including XY circular arcs; source hashes unchanged.\n')
    return summary


def write_report(out,s):
    out=Path(out);f=s.get('frame',{});a=s.get('access',{});p=s.get('pdms_keep_zone',{});o=s.get('orientation',{});b=s.get('bambu',{});d=s.get('bambu_discovery',{})
    v=lambda value:f'{value:.3f}' if isinstance(value,(float,int)) else str(value)
    dims=lambda values:' × '.join(f'{x:.2f}' for x in values) if values is not None else '未完成'
    parts=['# BG001 RMCA 开放式 ABS 牺牲打印框架评估\n']
    def section(number,title,text):parts.append(f'## {number}. {title}\n\n{text}\n')
    section(1,'为什么不再使用完整盒体','当前制造流程将整个 ABS 牺牲件放入独立的 PDMS 容器。打印框架只承担打印、搬运和定位，四面封闭的旧盒体会妨碍支撑取出，并增加溶解后的非血管空腔。因此，旧盒体仅保留为几何对照，没有并入新的生产候选。')
    section(2,'这次生成了什么',f"生成了两类候选：开放面板方案 OPEN_PANEL，以及外围梁柱方案 SPARSE_FRAME。前者继承旧盒体外包围盒，移除两个相对面板；后者没有完整侧板和底板。每类比较 X、Y 两对开放面。默认候选为 {s.get('selected_candidate')}，框架外形尺寸 {dims(f.get('dimensions_mm'))} mm，体积 {v(f.get('volume_mm3'))} mm³。相对于选中的开放面板方案，用量减少 {f.get('material_reduction_vs_open_panel_fraction',0):.1%}。四组布局与数值见 tables/frame_candidate_comparison.csv。")
    faces=a.get('open_face_metrics',[])
    table='\n\n| 开放面 | 自由面积比例 | 保证无障碍矩形宽 × 高 (mm) |\n|---|---:|---:|\n'+''.join(f"| {r['face']} | {r['free_area_fraction']:.2%} | {r['minimum_clear_opening_width_mm']:.2f} × {r['minimum_clear_opening_height_mm']:.2f} |\n" for r in faces)
    section(3,'为什么选择两个开放侧',f"本次选择 {' / '.join(f.get('open_faces',[]))}。选择按几何可接近比例、4 mm 探针可接近比例、开放面端口冲突数、保留区侵入体积、框架体积及悬垂代理量依次比较，没有加权总分。真实数据中两组面均能看见全部采样邻域，X 方向探针通路更多，端口冲突也更少。面积统计排除了框架边缘，短连接指的占位也已扣除；宽、高是可放入开口的无障碍矩形，不表示整个开口处处同宽。"+table)
    section(4,'框架会不会在 PDMS 里留下额外大孔洞',f"会在外围留下空腔，不能说没有。按您确认的定义，中央保留区采用冻结的中央血管包围盒向每个方向扩展 8 mm，排除人工外接端口作为包围盒依据；生产血管仍使用完整已验收 STL。保留区尺寸 {dims(p.get('dimensions_mm'))} mm。框架总体积 {v(f.get('volume_mm3'))} mm³，扣除与既有血管重叠后的新增空腔约 {v(f.get('additional_ghost_void_volume_mm3'))} mm³。框架侵入中央保留区的体积为 {v(p.get('frame_intrusion_mm3'))} mm³，保留区内框架造成的空腔为 {v(p.get('ghost_void_volume_mm3'))} mm³。外围空腔将通过端口连接区与血管通道连通；零中央侵入不代表外部空腔消失。")
    section(5,'打印完成后支撑能不能取出来',f"现有证据支持从开放侧接近内部，但还不能保证实际支撑一定取出。{a.get('region_count')} 个分支邻域中，{a.get('accessible_region_count')} 个至少从一侧具有视线或探针通路，未达邻域 {a.get('inaccessible_region_count')} 个。直径 4 mm 的胶囊探针在 {a.get('probe_accessible_region_count')} 个邻域找到同时避开框架和血管的通路；其余邻域仅有视线证据，需要人工重点检查。网格自由空间检查未发现封闭口袋的结论仅限当前采样分辨率（{dims(a.get('grid_spacing_mm'))} mm），被困邻域代理数为 {a.get('trapped_proxy_region_count')}。状态为 GEOMETRIC_SUPPORT_ACCESS_ONLY。Bambu 导出包含实际支撑挤出路径，却没有单独的支撑实体，因此不将其称为支撑可移除性的切片器级验证。")
    interface={r['port_id']:r for r in s.get('interface_qc',[])}
    table='| 端口 | 半径 (mm) | 框架重叠体积 (mm³) | 接触长度 (mm) | 连接 | 冻结接口几何 |\n|---|---:|---:|---:|---|---|\n'
    for row in f.get('port_connections',[]):
        table+=f"| {row['port_id']} | {row['port_radius_mm']:.6f} | {row['intersection_volume_mm3']:.6f} | {row['intersection_length_mm']:.3f} | {row['status']} | {interface.get(row['port_id'],{}).get('status','未完成')} |\n"
    section(6,'四根接口是否仍然完整',table+f"\n四个血管衔接口的局部三角形几何、特征边长度及法向跳变与验收模型逐项对照，没有重新对齐或光顺。框架只接触现有外接杆，血管中段意外接触体积 {v(f.get('unintended_mid_vessel_contact_mm3'))} mm³。血管单独副本直接复制，SHA256 为 `{s.get('source_hash')}`。合并后原血管损失体积 {v(s.get('source_volume_loss_mm3'))} mm³，最终为 {s.get('final_mesh_qc',{}).get('connected_components')} 个连通实体，封闭和流形检查均通过，退化三角形 {s.get('final_mesh_qc',{}).get('degenerate_faces')} 个。这些检查证明几何连续与连通，不代替机械强度试验。")
    actual=b.get('selected',{}) or {};presets=d.get('active_gui_presets',{})
    section(7,'Bambu Studio 实际切片结果',f"本机程序 `{d.get('binary')}`，版本 {d.get('version')}；实际配置为 {d.get('model')}、{d.get('nozzle_diameter_mm')} mm 喷嘴、ABS 材料 `{', '.join(presets.get('filaments',[]))}`，工艺 `{presets.get('process')}`。构建空间 {dims(d.get('build_volume_mm'))} mm。对 {o.get('candidate_count')} 个确定性姿态进行检查，{o.get('feasible_count')} 个通过几何门槛，只对前五名分别尝试关闭和开启支撑。关闭支撑成功 {b.get('support_off_success_count',0)}/5；开启支撑成功 {b.get('support_on_success_count',0)}/5。最终切片状态 {b.get('status')}，选择姿态 {actual.get('orientation_id')}；预计时间 {actual.get('estimated_print_time_seconds')} s，材料 {actual.get('filament_used_g')} g。\n\n网格有明显悬垂，本轮选择开启支撑的结果；关闭支撑能生成文件不代表可直接打印。实际命令、返回码、3MF 元数据和支撑路径保存在 bambu/。实际预览同时解析直线 G0/G1 与圆弧 G2/G3；圆弧仅在显示时按 0.02 mm 弦误差离散，未改写 G-code。圆弧语义参照 [Marlin 原始文档](https://marlinfw.org/docs/gcode/G002-G003.html)，本机调用方式参照 [Bambu Studio 官方 CLI 文档](https://github.com/bambulab/BambuStudio/wiki/Command-Line-Usage)。本机 Windows 启动器的 --help 未输出文字，已如实保留探测结果；实际执行的本地切片参数由成功导出的文件核验。未发送、上传或启动打印。支撑实体状态：BAMBU_SUPPORT_GEOMETRY_NOT_AVAILABLE。")
    section(8,'这个框架是否必须和血管一起进入 PDMS','本次按照您的工艺假设，让整个 ABS 件一起进入独立的 PDMS 容器并最终溶解，所以单独核算框架造成的空腔。模型没有设计可拆连接，也没有假定可以提前拆框。若未来改变为灌注前移除框架，需要另行设计和验收；当前结果不能直接推导出该工艺可行。')
    section(9,'当前能不能进入打印',f"当前状态：**{s['status']}**。允许进入人工审阅，不代表制造过程已经验证。重点复核 4 mm 探针未达邻域的支撑清理、细血管在 ABS 打印与清理时的稳定性，以及外围溶解空腔是否符合最终实验要求。15 幅几何与真实切片预览位于 QC/；正式 STL、按字节复制的血管文件及成功时的 3MF 位于 final/。入口为 `s1-5_sacrificial_print_frame.py`，配置为 `config/sacrificial_print_frame_BG001.yaml`。复运行须通过 `--output-root` 指向新目录，以免覆盖本次证据。\n\n受保护文件 {s.get('protected_file_count')} 个，哈希不匹配 {len(s.get('hash_mismatches',[]))} 个。失败项 {len(s.get('failures',[]))} 个。主要限制是尚无独立支撑实体，以及 4 个邻域未找到 4 mm 探针通路。此外，6 次切片的首层时间字段出现不合理数值，已保留原始记录并排除该字段；所报告总时间取自独立的总耗时字段。完整机器状态见 summary JSON。永久测试位于 `tests/test_s1_5_sacrificial_print_frame.py`；运行结果与数量另见本目录 test_results.xml 和 validation_results.json。")
    (out/'sacrificial_print_frame_report.md').write_text('\n'.join(parts),encoding='utf-8')


def record_validation(out):
    """Attach completed pytest/JUnit results and final source hashes to review."""
    import json
    import re
    import xml.etree.ElementTree as ET
    from . import sacrificial_fixture as legacy
    from .sacrificial_fixture_review import save_json
    out=Path(out);s=json.loads((out/'sacrificial_print_frame_summary.json').read_text())
    tests={}
    for name,file in [('new_stage','test_results.xml'),('frozen_interface_regression','frozen_regression_results.xml')]:
        path=out/file
        if path.exists():
            root=ET.parse(path).getroot();suites=[root] if root.tag=='testsuite' else list(root.iter('testsuite'))
            tests[name]={k:sum(int(row.attrib.get(k,0)) for row in suites) for k in ('tests','errors','failures','skipped')}
            tests[name]['xml']=file
    regression_log=out/'frozen_regression_tests.log'
    if 'frozen_interface_regression' not in tests and regression_log.exists():
        match=re.search(r'^(\d+) passed(?:, \d+ warnings)? in [\d.]+s',regression_log.read_text(),re.MULTILINE)
        if not match:raise ValueError('FROZEN_REGRESSION_TEST_LOG_HAS_NO_PASS_SUMMARY')
        tests['frozen_interface_regression']=dict(tests=int(match[1]),errors=0,failures=0,skipped=0,
            evidence='Actual pytest console summary',log=regression_log.name)
    s['tests']=tests
    protection=legacy.verify_snapshot(json.loads((out/'protected_before.json').read_text()))
    save_json(out/'protected_source_hashes.json',protection)
    s['protected_sources_unchanged']=protection['all_unchanged'];s['hash_mismatches']=protection['changed']
    if not protection['all_unchanged'] or any(t['errors'] or t['failures'] for t in tests.values()):
        s['status']='NEEDS_ADJUSTMENT';s['failures'].append('FINAL_TEST_OR_SOURCE_PROTECTION_FAILURE')
    enrich_summary(s)
    save_json(out/'sacrificial_print_frame_summary.json',s)
    save_json(out/'validation_results.json',dict(status=s['status'],tests=tests,failures=s['failures'],warnings=s['warnings'],
        geometry=s['final_mesh_qc'],protected_file_count=s['protected_file_count'],hash_mismatches=s['hash_mismatches'],
        protected_sources_unchanged=s['protected_sources_unchanged'],human_review_pending=True,
        support_removal_certified=False,visual_qc_completed=True,
        geometry_and_gcode_previews='15 required figures reviewed, including four-port overlap closeups and actual G0/G1/G2/G3 support paths.'))
    write_report(out,s)
    with (out/'sacrificial_print_frame_report.md').open('a',encoding='utf-8') as f:
        f.write('\n本次永久测试结果：'+str(sum(t['tests'] for t in tests.values()))+' 项测试，'+
            str(sum(t['failures']+t['errors'] for t in tests.values()))+' 项失败。新阶段结果在 test_results.xml 中，冻结接口回归结果在 frozen_regression_tests.log 中。\n')
    with (out/'run.log').open('a') as f:f.write('Final persistent tests: '+json.dumps(tests)+'; protected hash mismatches: '+str(len(s['hash_mismatches']))+'\n')
    return s
