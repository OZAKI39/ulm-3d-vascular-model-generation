"""RFUS review: measured coverage, costs, removal evidence and explicit failures."""
from __future__ import annotations
import json,pickle,shutil,zipfile
from pathlib import Path
import numpy as np
import trimesh
import shapely
from shapely.geometry import GeometryCollection
from .support_inputs import save,table,verify,numeric
from .support_v3_1_review import polygons_mesh,pvmesh
from .support_v3_geometry import union


def image(path,title,meshes=(),points=(),lines=(),note='',camera='iso'):
    """RFUS-only review renderer; historical plots and interactive UI stay intact."""
    import pyvista as pv
    plot=pv.Plotter(off_screen=True,window_size=(1450,1050));plot.set_background('white')
    for mesh,color,opacity in meshes:
        if mesh is not None and len(mesh.faces):plot.add_mesh(pvmesh(mesh),color=color,opacity=opacity)
    for xyz,color,size in points:
        if len(xyz):plot.add_points(np.asarray(xyz),color=color,point_size=size,render_points_as_spheres=True)
    for item in lines:
        xyz,color=item[:2];width=item[2] if len(item)>2 else 3
        if len(xyz)>1:plot.add_lines(np.asarray(xyz),color=color,width=width,connected=True)
    plot.add_text(title,font_size=17,color='#17212b',position='upper_left')
    plot.add_text(note,font_size=11,color='#333333',position='lower_left');plot.add_axes()
    if camera in ('top','xy'):plot.view_xy()
    elif camera=='xz':plot.view_xz()
    elif camera=='yz':plot.view_yz()
    else:plot.view_isometric()
    plot.reset_camera();plot.camera.zoom(.87);path.parent.mkdir(parents=True,exist_ok=True);plot.screenshot(str(path));plot.close()


def native_self_support(parsed, select, footprints, box_footprints, profile):
    """Native supports may use 0.2 mm layers while the model uses 0.1 mm."""
    rows=[];moves=parsed['extrusion'][select];widths=parsed['extrusion_widths'][select]
    heights=parsed['extrusion_heights'][select];tops=np.round(moves[:,:,2].mean(1),5)
    groups=np.unique(np.round(np.c_[tops,heights],5),axis=0)
    for z,h in groups:
        if z<=profile['first_layer_height']+.001:continue
        ix=(np.abs(tops-z)<1e-5)&(np.abs(heights-h)<1e-5)
        current=shapely.union_all(shapely.buffer(shapely.linestrings(moves[ix,:,:2]),widths[ix]/2,cap_style='flat'))
        previous=shapely.union_all([p for source in [footprints,box_footprints] for zz,p in source.items() if abs(zz-(z-h))<.001])
        missing=current.difference(previous.buffer(min(profile['line_width']/2,h)))
        for p in shapely.get_parts(missing):
            if p.area>.005:rows.append(dict(z_mm=float(z),actual_layer_height_mm=float(h),area_mm2=float(p.area)))
    return dict(unsupported_region_count=len(rows),rows=rows,method='Archived per-move LAYER_HEIGHT; same actual-layer growth and 0.005 mm2 area criterion')


def comparative(c,records):
    from .support_rfus_travel import events,actual_sections
    from .support_rfus_validation import underbody_coverage
    from .support_v3_validation import actual_self_support
    out=c['out']/'comparison';out.mkdir(exist_ok=True);rows=[]
    for name,folder,archive,native in [
        ('BAMBU_ORGANIC',c['base']/'slicer/bambu_native_current',c['base']/'slicer/bambu_native_current/sliced.3mf',True),
        ('V3_SPARSE',c['base']/'v3_1/iteration_00',c['base']/'v3_1/iteration_00/X1C_project.3mf',False)]:
        cache=out/(name+'.json')
        if cache.exists() and json.loads(cache.read_text()).get('comparison_metric_version')==2 and json.loads(cache.read_text()).get('target_fingerprint')==c.get('corridor_digest'):rows.append(json.loads(cache.read_text()));continue
        parsed=events(archive);boxfp=actual_sections(parsed,c['box'],c['profile']);vfp=actual_sections(parsed,c['vascular'],c['profile'])
        if native:
            select=np.array(['support' in tag.lower() for tag in parsed['features']]);segments=parsed['extrusion'][select];widths=parsed['extrusion_widths'][select];z=np.round(segments[:,:,2].mean(1),5)
            fp={float(level):shapely.union_all(shapely.buffer(shapely.linestrings(segments[z==level,:,:2]),widths[z==level]/2,cap_style='flat')) for level in np.unique(z)}
            filament=float(parsed['extrusion_filament_mm'][select].sum());baseline=json.loads((folder/'slice_result.json').read_text());seconds=baseline.get('estimated_print_time_seconds');segment_count=None;trapped=None;access=None;volume=None
        else:
            mesh=trimesh.load_mesh(folder/'custom_support_v3.stl');fp,filament=actual_sections(parsed,mesh,c['profile'],with_filament=True);old=json.loads((folder/'candidate.json').read_text());seconds=old['slice'].get('estimated_print_time_seconds');segment_count=old['removal'].get('fragment_count');trapped=old['removal']['trapped_pieces'];access=None;volume=float(mesh.volume)
        model={z:boxfp.get(z,GeometryCollection()).union(vfp.get(z,GeometryCollection())) for z in boxfp.keys()|vfp.keys()};cov=underbody_coverage(c['corridors'],fp,model,c['profile']);selfcheck=native_self_support(parsed,select,fp,boxfp,c['profile']) if native else actual_self_support(fp,boxfp,c['profile'])
        with zipfile.ZipFile(archive) as z:settings=json.loads(z.read('Metadata/project_settings.config'))
        diameter=numeric(settings['filament_diameter']);actualvolume=filament*np.pi*(diameter/2)**2
        r=dict(name=name,centerline_underbody_percent=cov['percent'],max_uncovered_arc_mm=cov['max_uncovered_underbody_arc_mm'],support_self_unsupported=selfcheck['unsupported_region_count'],actual_support_extrusion_volume_mm3=actualvolume,geometry_support_volume_mm3=volume,segment_count=segment_count,trapped_segments=trapped,accessible_handles=access,estimated_print_time_seconds=seconds,archive=str(archive),metric='Same complete 1990-sample corridor target; archived variable extrusion width; actual E volume',unavailable_metrics='Native Organic has no designed RFUS segment or handle identity; removal not inferred')
        r.update(comparison_metric_version=2,self_support_method=selfcheck['method'],target_fingerprint=c.get('corridor_digest'))
        save(cache,r);table(out/(name+'_underbody.csv'),cov['rows']);rows.append(r)
    for r in records:
        s=r['slice'];rows.append(dict(name=r['candidate_id'],centerline_underbody_percent=s['underbody']['percent'],max_uncovered_arc_mm=s['underbody']['max_uncovered_underbody_arc_mm'],support_self_unsupported=s['self_printability']['unsupported_region_count'],actual_support_extrusion_volume_mm3=s.get('actual_support_extrusion_volume_mm3'),geometry_support_volume_mm3=r['support_volume_mm3'],segment_count=r['segment_count'],trapped_segments=r['removal']['trapped_segments'],accessible_handles=r['removal']['accessible_handles'],estimated_print_time_seconds=s['estimated_print_time_seconds'],failed_gates=r['failed_gates']))
    table(c['out']/'method_comparison.csv',rows);return rows


def scenes(c,record):
    import pyvista as pv
    folder=Path(record['folder']);g=folder/'geometry';qc=c['out']/'QC_ROBUST';qc.mkdir(exist_ok=True)
    with (folder/'segment_geometry.pkl').open('rb') as f:segments=pickle.load(f)
    v=c['vascular'];box=c['box'];support=trimesh.load_mesh(g/'custom_full_underbody_support.stl')
    meshes={k:trimesh.load_mesh(g/f'support_{k}.stl') for k in ['roofs','curtains','handles','necks']}
    base=[(box,'#cccccc',.08),(v,'#963f46',.30)];blue='#548bb1';orange='#d78b35'
    note='RFUS | O000 | X1 Carbon / 0.2 mm | software diagnostics; physical parameters uncalibrated'
    image(qc/'01_frozen_input.png','Frozen vascular + ports + 2 mm box',base,note=note)
    lines=[(b['lower'][np.isfinite(b['lower']).all(1)],'#126284') for b in c['corridors']]
    image(qc/'02_all_underbody_corridors.png','All branches and ports: full lower-envelope targets',base,lines=lines,note='No FreeFloating or sparse risk mask selects support')
    image(qc/'03_full_roof_ridges.png','Narrow ridges | actual profile line multiples',base+[(meshes['roofs'],orange,1)],note='Minimum gap 0.20 mm + layer quantization; no U-shaped or circumferential saddle')
    image(qc/'04_full_curtains.png',f"{record['segment_count']} separate curtain segments",base+[(meshes['curtains'],blue,1)],note=note)
    seams=np.array([s['position'] for s in record['seams']]).reshape(-1,3)
    image(qc/'05_segment_seams.png','Full-height separation seams',[(support,blue,.7)],points=[(seams,orange,8)],note='Short roof breaks remain in the actual coverage audit; no automatic bridge credit')
    bif=np.asarray(c['topology']['bifurcations']).reshape(-1,3)
    image(qc/'06_bifurcation_segments.png','Parent and daughters retain underbody targets',base+[(support,blue,.6)],points=[(bif,orange,11)],note='Shorter segments near bifurcations; all unresolved routes are recorded')
    portparts=union([s['mesh'] for s in segments if s['is_port']]);image(qc/'07_external_port_support.png','Port curtains: inner bottom / external plate roots',base+[(portparts,blue,1)],note=note)
    image(qc/'08_removal_handles.png',f"{record['handle_count']} individual removal handles",base+[(meshes['curtains'],blue,.25),(meshes['handles'],orange,1)],note='Each grip belongs to one segment; strength and manual force uncalibrated')
    image(qc/'09_base_necks.png','Cuttable base necks | explicit isolated cut planes',[(support,blue,.2),(meshes['necks'],orange,1)],note='Not a calibrated breakaway-force claim')
    image(qc/'10_support_only_top.png','Support only: top view',[(support,blue,1)],camera='xy',note=note)
    image(qc/'11_support_only_side.png','Support only: side view',[(support,blue,1)],camera='xz',note=note)
    part=pvmesh(support).clip(normal='x',origin=support.bounds.mean(0)).triangulate()
    clipped=trimesh.Trimesh(part.points,np.asarray(part.faces).reshape(-1,4)[:,1:],process=False)
    image(qc/'12_support_only_section.png','Support section at mid-X',[(clipped,blue,1)],note='Clipped review view; original support geometry is unchanged')
    image(qc/'13_support_plus_vascular.png','Support and frozen vascular model',base+[(support,blue,.65)],note=note)
    u=json.loads((folder/'underbody_coverage.json').read_text());good=np.array([r['xyz'] for r in u['rows'] if r['supported']]).reshape(-1,3);bad=np.array([r['xyz'] for r in u['rows'] if not r['supported']]).reshape(-1,3)
    image(qc/'14_centerline_underbody_coverage.png',f"Actual underbody {u['percent']:.4f}% | longest gap {u['max_uncovered_underbody_arc_mm']:.3f} mm",[(v,'#bbbbbb',.18)],points=[(good,'#16884a',5),(bad,'#d72232',8)],note='Green: actual G-code supported | Red: actual unsupported | full arc sampling')
    with (folder/'actual_layer_geometry.pkl').open('rb') as f:hit,miss=pickle.load(f)
    image(qc/'15_actual_toolpath_coverage.png',f"Required layer area: {record['slice']['required_layer_coverage']['percent']:.4f}%",[(v,'#cccccc',.15),(polygons_mesh(hit),'#16884a',1),(polygons_mesh(miss),'#d72232',1)],note='Green/red from actual per-move-width footprints, not nominal support STL')
    sf=record['slice']['self_printability'];points=np.array([[*r['representative_xy'],r['z_mm']] for r in sf['rows']]).reshape(-1,3)
    image(qc/'16_support_self_printability.png',f"Actual support self-unsupported regions: {sf['unsupported_region_count']}",[(support,blue,.6)],points=[(points,'red',8)],note='Previous actual support and approved box extrusion; no secondary support tree')
    pathlines=[]
    for s in record['removal']['segments']:
        if s['removal_path']:pathlines.append((np.array([np.array(m)[:3,3] for m in s['removal_path']['transforms']]),'#258a44'))
    image(qc/'17_removal_order.png','Swept removal routes after individual base cuts',base+[(support,blue,.25)],lines=pathlines,note=f"Found sequence: {len(record['removal']['sequence'])}; trapped: {record['removal']['trapped_segments']}; remaining segments are obstacles")
    lookup={r['segment_id']:r for r in record['removal']['segments']};hg=[];hb=[]
    for s in segments:
        if s.get('handle_target'):(hg if lookup[s['segment_id']]['handle_accessible'] else hb).append(s['handle_target'])
    image(qc/'18_handle_accessibility.png','3 mm tool probe: individual handle access',base+[(meshes['handles'],orange,.8)],points=[(np.array(hg).reshape(-1,3),'green',10),(np.array(hb).reshape(-1,3),'red',10)],note='Geometric tool proxy; accessible handle may replace a deep base-cut strategy')
    travel=record['travel'];risk=sorted(travel['rows'],key=lambda r:r['nozzle_tip_min_clearance_mm'])[:250]
    image(qc/'19_travel_near_miss.png',f"Actual travel near-miss proxy: {travel['near_miss_segments']} moves",base+[(meshes['handles'],orange,.7)],lines=[(np.array([r['start'],r['end']]),'#d72232') for r in risk],note='Worst 250 displayed; all moves in CSV. Only previously extruded material participates')
    plot=pv.Plotter(off_screen=True,window_size=(1800,1250));plot.set_background('white');plot.add_mesh(pvmesh(v),color='#a85659',opacity=.22);plot.add_mesh(pvmesh(box),color='#cccccc',opacity=.08)
    palette=['#4b8dba','#65a887','#aa83b6','#d4974b','#7d9cb9','#b17983'];centers=[];labels=[]
    for i,s in enumerate(segments):
        plot.add_mesh(pvmesh(s['mesh']),color=palette[i%len(palette)],opacity=.9)
        centers.append(s.get('handle_target') or s['mesh'].bounds.mean(0));r=lookup[s['segment_id']];order=r.get('sequence_index');labels.append(f"{i+1} / {order+1 if order is not None else '?'}")
    table(qc/'segment_display_index.csv',[dict(display_number=i+1,segment_id=s['segment_id'],branch_id=s['branch_id'],removal_order=(lookup[s['segment_id']]['sequence_index']+1) if 'sequence_index' in lookup[s['segment_id']] else None,status=lookup[s['segment_id']]['status']) for i,s in enumerate(segments)])
    plot.add_point_labels(np.asarray(centers),labels,font_size=12,point_size=2,shape_opacity=.5,always_visible=False)
    plot.add_text(f"{record['candidate_id']} | {record['status']}\nSegment / removal order (? unresolved); no physical validation",font_size=17,color='black');plot.add_axes();plot.view_isometric();plot.reset_camera();plot.screenshot(str(qc/'20_final_human_review.png'));plot.close()


def generate(c,records,runtime):
    records=[r for r in records if r.get('evaluation_complete')]
    if not records:raise ValueError('NO_ACTUAL_RFUS_RESULTS')
    from .support_rfus_removal import topology_audit
    from .support_rfus_validation import gates,refresh_underbody
    for r in records:
        refresh_underbody(c,r)
        if 'segment_topology' not in r:topology_audit(r)
        r['failed_gates']=gates(r)
        if not r['segment_topology']['all_independent'] and 'SEGMENT_NOT_INDEPENDENT' not in r['failed_gates']:r['failed_gates'].append('SEGMENT_NOT_INDEPENDENT')
        r['status']='NEEDS_ADJUSTMENT' if r['failed_gates'] else 'READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED'
        save(Path(r['folder'])/'candidate.json',r)
    good=[r for r in records if not r['failed_gates']]
    # Coverage first; failed candidates remain diagnostics, never production.
    selected=min(good or records,key=lambda r:(-r['slice']['underbody']['percent'],r['removal']['trapped_segments'],r['slice']['self_printability']['unsupported_region_count'],r['support_volume_mm3']))
    out=c['out'];folder=Path(selected['folder']);comparisons=comparative(c,records);scenes(c,selected)
    protection=verify(json.loads((out/'validation/protected_before.json').read_text()));save(out/'validation/source_protection.json',protection)
    for name in ['curtain_segments.csv','roof_ridges.csv','removal_handles.csv','base_necks.csv','underbody_coverage.csv','support_self_printability.csv','segment_removal.csv','travel_near_miss.csv','bambu_validation.json']:
        shutil.copy2(folder/name,out/name)
    shutil.copytree(folder/'geometry',out/'geometry',dirs_exist_ok=True)
    production=None
    if good:
        production=out/'production_selected';production.mkdir(exist_ok=True);shutil.copytree(folder/'geometry',production/'geometry',dirs_exist_ok=True);shutil.copy2(folder/'X1C_diagnostic_full_support.3mf',production/'X1C_robust_full_support.3mf')
    status='READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED' if good else 'NEEDS_ADJUSTMENT'
    summary=dict(mode='ROBUST_FULL_UNDERBODY_SUPPORT',status=status,output=str(out),diagnostic_selected=selected['candidate_id'],production_selected=str(production) if production else None,selected_folder=str(folder),candidates=[{k:r[k] for k in ['candidate_id','status','folder','segment_count','handle_count','support_volume_mm3','failed_gates']} for r in records],metrics=selected['slice'],removal=selected['removal'],travel={k:v for k,v in selected['travel'].items() if k!='rows'},comparison=comparisons,source_protection=protection,scale_factor=1.,orientation='O000',historical_research_preserved=True,physical_status=selected['physical_status'],printer_commands_sent=False,runtime_seconds=runtime)
    save(out/'robust_full_support_summary.json',summary)
    s=selected['slice'];rem=selected['removal'];u=s['underbody'];v3=next(r for r in comparisons if r['name']=='V3_SPARSE');native=next(r for r in comparisons if r['name']=='BAMBU_ORGANIC')
    answers=[
        '1. 已停止追求最少支撑。本轮以所有分支、端口的完整下托路径为目标，旧的稀疏优化与姿态搜索只保留作研究对照。',
        f"2. 是。全部 {len(c['corridors'])} 条原始分支及端口路径、{u['sample_count']} 个固定弧长样本均进入目标集合。几何上已有盒体的位置也须有真实刀路承托证据；路由失败不会从分母删去。",
        f"3. 当前选出的{'审阅' if good else '诊断'}件是 {selected['candidate_id']}，包含 {selected['segment_count']} 段帘状支撑。",
        f"4. 每段平均弧长 {selected['mean_segment_length_mm']:.3f} mm；高曲率及分叉附近采用更短分段，局部避让可进一步缩短。",
        f"5. 该候选主帘厚度为 {selected['curtain_lines']} 条线宽。读取的配置线宽为 {c['profile']['line_width']:.3f} mm；同时比较 2/3 线帘和 1/2 线底颈。",
        f"6. 中心线下侧实际刀路覆盖为 {u['percent']:.6f}%；完整需支撑层多边形覆盖为 {s['required_layer_coverage']['percent']:.6f}%。两个指标分母不同，不相互替代。",
        f"7. 最长未承托弧长为 {u['max_uncovered_underbody_arc_mm']:.6f} mm，未覆盖样本为 {u['unsupported_samples']} 个。",
        f"8. 导出支撑与冻结血管/端口相交体积为 {selected['vascular_port_intersection_mm3']:.9g} mm³；与盒体的非法相交体积为 {selected['illegal_box_intersection_mm3']:.9g} mm³。所有偏移及双帘尝试保存在 route_attempts.csv。",
        '9. 没有生成环抱式鞍座或分叉包壳。屋脊沿下表面，仅为窄条；分叉仍有未完成承托的位置时会显示在红色失败图上。',
        f"10. 支撑自身仍有 {s['self_printability']['unsupported_region_count']} 个实际欠支撑层区域。这是打印门槛，不能因 STL 看起来连续而忽略。",
        f"11. 生成 {selected['handle_count']} 个独立抓取柄，当前顺序中 {rem['accessible_handles']} 个通过 3 mm 工具代理检查。",
        f"12. {'是' if rem['inaccessible_segments']==0 else '否'}。仍有 {rem['inaccessible_segments']} 个内部支撑段没有通过柄/底颈任一操作通路检查。柄可达允许替代深处底颈可达，但不代表真实剥离力已验证。",
        f"13. 仍有 {rem['trapped_segments']} 个受困支撑段；已寻找 {len(rem['sequence'])} 段的拆除顺序。模拟包括血管、端口、盒壁和未拆支撑，先断开设计底颈再检查整段与柄的连续扫掠。",
        f"14. Bambu Auto Support 已关闭：{s['auto_support_off']}；实际原生支撑刀路线段数为 {s['native_support_toolpaths']}。",
        f"15. 真实 X1 Carbon 整模切片成功：{s['slicer_success']}；支撑加盒体基础的独立切片成功：{s['support_only_slicer_success']}。切片成功不等同于制造合格。",
        '16. 本方法把承托目标扩展到整个弧长，结构和拆除单元更明确；是否比 V3.1 更可靠须看完整刀路与拆除结果。'+('全部软件门槛已通过。' if good else '本轮仍未通过全部门槛，不能宣称已经比旧方法可靠。'),
        '17. 相比 Organic，分段、底颈与独立柄使预期操作方式更明确；没有实物拆除试验，不能声称已经更容易拆。原生树没有本方法的段/柄编号，缺少的对照拆除指标以空值记录，没有编造。',
        '18. 0.20 mm 屋脊间隙、底颈切割/剥离行为、帘强度、柄强度和实际操作力均未经过实物校准。CUTTABLE BASE NECK IS NOT PHYSICALLY CALIBRATED。',
        '19. '+('可进入人工审阅，仍需校准物理参数。' if good else '当前不建议进入完整模型实物试打。状态为 NEEDS_ADJUSTMENT；STL 与 3MF 作为诊断结果保留。')]
    report='# BG001 RMCA 全下托帘状支撑：实现与真实验证\n\n'+'\n\n'.join(answers)+'\n\n## 实际比较与代价\n\n'
    report+='所有比较均重新读取历史真实 G-code，并按本轮相同的完整弧长样本和逐段 LINE_WIDTH 测量。没有将 V3 的旧稀疏层覆盖百分比与 RFUS 的完整中心线百分比混为一谈。\n\n'
    report+='| 方法 | 全弧长下托覆盖 | 最长缺口 mm | 自身欠支撑层区 | 实际支撑挤出体积 mm³ | 预计总时长 h |\n|---|---:|---:|---:|---:|---:|\n'
    for r in comparisons:report+=f"| {r['name']} | {r['centerline_underbody_percent']:.4f}% | {r['max_uncovered_arc_mm']:.4f} | {r['support_self_unsupported']} | {r['actual_support_extrusion_volume_mm3'] or 0:.2f} | {(r['estimated_print_time_seconds'] or 0)/3600:.2f} |\n"
    for name,old in [('V3 稀疏肋',v3),('Bambu Organic',native)]:
        dv=s.get('actual_support_extrusion_volume_mm3',0)-old['actual_support_extrusion_volume_mm3'];dt=(s['estimated_print_time_seconds']-old['estimated_print_time_seconds'])/60
        report+=f"\n相对 {name}，当前诊断件实际归属于支撑的挤出体积变化 {dv:+.2f} mm³，总切片预计时间变化 {dt:+.2f} min。挤出量来自实际 E 值及线段归属，不代表称重实测材料。\n"
    report+='\n## 方法及证据边界\n\n中心线采样间距不超过真实线宽；所有下表面均来自冻结 STL 的射线交点，未用半径猜测替换真实表面。新支撑独立生成，完整模型仅作固定 O000 的刚体平移，缩放为 1。分段缝与共享投影区域仍参加实际覆盖测量。\n\n'
    report+='端口端面上的共面射线采用一个 float32 坐标间距以内的数值偏移重试，交点仍取自冻结 STL。若采样缓存因此更新，会对同一归档 G-code 重新测量，并保留更新前结果；不会把无交点简单记为已承托。\n\n'
    report+=f"屋脊顶面按 {c['profile']['layer_height']:.2f} mm 层高向下对齐；设计间隙以下限 {c['profile']['roof_gap_mm']:.2f} mm 加一个层高以内的离散化表示。该层高量化在刀路覆盖审计中单列，不作为额外覆盖面积。\n\n"
    report+='实际刀路使用归档的逐段 LINE_WIDTH，包含变线宽壁。覆盖容差只来自 float32 坐标精度；分层高度量化单独记录，未为获得 100% 而放宽面积容差。末端和壁内点仍保留在分母。\n\n'
    report+=f"移动近失检查发现 {selected['travel']['near_miss_segments']} 个风险移动段。检查真实非挤出移动，仅计入当时已经挤出的材料；使用喷嘴尖与可配置圆柱代理，并排除一条线宽内的离开/落点接触。这是 GEOMETRIC TRAVEL PROXY，不是完整喷头碰撞证明。\n\n"
    report+='运行中修复了半线宽分段缝条件与闭合帘体端面三角化，并把实际挤出宽度与恒定配置宽度区分开。开发中间结果与日志保留在 history/ 和 validation/，不替代正式候选的检查结果。\n\n'
    report+=f"支撑 STL 导出检查错误：{selected.get('support_export_error') or '无'}；合并 STL 导出检查错误：{selected.get('assembly_export_error') or '无'}。若有错误，文件仅供定位问题，切片成功也不能消除此硬性失败项。\n\n"
    report+=f"断开底颈后，有 {selected['segment_topology']['disconnected_or_invalid_segments']} 段不是单一闭合材料连通体，另有 {len(selected['segment_topology']['touching_segment_pairs'])} 对段存在名义接触。其整体运动路径即使通过，也不能据此声称可用一个柄拿走整段；此项单独阻止生产通过。\n\n"
    report+=f"受保护文件 {protection['checked_files']} 个，SHA256 差异 {len(protection['changed'])} 个。未创建 s1-9，未修改既有可视化，未上传、发送或启动打印。\n\n失败门槛：{', '.join(selected['failed_gates']) or '无'}。\n"
    (out/'robust_full_underbody_report.md').write_text(report)
    (out/'geometry/README.md').write_text(f"# 本轮模型状态\n\n{status}。本目录来自 {selected['candidate_id']}；存在 STL 不表示已经可打印。真实切片和完整失败项见 ../robust_full_underbody_report.md。\n")
    return summary
