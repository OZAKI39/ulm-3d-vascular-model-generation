"""Evidence-led Chinese report, 20 QC figures and protected V2 comparison."""
from __future__ import annotations
import csv,json,pickle,shutil
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save,table,digest
from .support_v3_pipeline import gates,rank

STEMS=['v2_best_reference','v3_orientation_candidates','v3_critical_intervals','long_tubular_intervals','rib_segmentation','rib_windows',
       'roof_ridges','rib_self_printability','support_only_top','support_only_side','support_only_section','support_plus_vascular',
       'cut_points','tool_access','fragment_removal_sequence','actual_toolpath_coverage','uncovered_regions','v2_vs_v3','bambu_warning_localization','final_candidate_review']


def compare(out,best,profile):
    v2=out.parent/'v2';old=json.loads((v2/'candidates/O000/candidate.json').read_text())
    cache=out/'validation/v2_same_gate_actual_coverage.json'
    if not cache.exists() or 'contact_area_mm2' not in json.loads(cache.read_text()):
        from .support_qc import parse_actual_paths
        from .support_v3_validation import extrusion_sections,coverage,contact_burden
        source=json.loads((v2/'candidates/O000/bambu/slice_result.json').read_text())
        paths=parse_actual_paths(source['archive']);segments=paths['segments'][paths['kind']==0]
        mesh=trimesh.load_mesh(v2/'candidates/O000/custom_support.stl');footprints,kept=extrusion_sections(segments,mesh,profile)
        # Match orientation O000 and the bridgeless, strict polygon criterion.
        with (out/'candidates/O000/frozen_detection_cache.pkl').open('rb') as f:_,field,_,_,_=pickle.load(f)
        actual=coverage(field,footprints,profile)[0]
        metric=contact_burden(kept,trimesh.load_mesh(v2/'candidates/O000/vascular_with_ports.stl'),profile)
        save(cache,{**{k:v for k,v in actual.items() if k!='rows'},'contact_area_mm2':metric['contact_area_mm2']})
        table(out/'validation/v2_same_gate_actual_coverage.csv',actual['rows'])
    oldcov=json.loads(cache.read_text());oldprofile=json.loads((v2/'profile_runtime.json').read_text())
    same_profile=all(oldprofile.get(k)==profile.get(k) for k in ['nozzle_diameter','line_width','layer_height','first_layer_height','process_name','printer','material','speeds','temperatures'])
    rows=[]
    for name,r,cov in [('V2_O000',old,oldcov),('V3_'+best['candidate_id'],best,best.get('slice',{}).get('actual_vascular_coverage',{}))]:
        s=r.get('slice',{});rem=r.get('removal',{})
        rows.append(dict(method=name,actual_coverage_percent=cov.get('percent'),coverage_pass=cov.get('all_required_covered'),same_profile=same_profile,
            coverage_gate='NO_BRIDGE_CREDIT_ACTUAL_GCODE_REQUIRED_POLYGONS_1e-5_mm2',
            self_unsupported=r.get('self_buildability',{}).get('unsupported_region_count',r.get('self_buildability',{}).get('unsupported_regions')),
            self_metric='FROZEN_GEOMETRY_LAYER',actual_self_unsupported_regions=s.get('actual_self_buildability',{}).get('unsupported_region_count'),
            near_vessel_burden_mm2=oldcov.get('contact_area_mm2') if name.startswith('V2') else s.get('contact_area_mm2'),
            burden_definition='downward_triangle_area_near_actual_isolated_support_paths',
            burden_percent_improvement=None,burden_comparison_reason='No percentage: both methods must first pass complete actual coverage',
            support_volume_mm3=r.get('support_volume_mm3'),root_count=r.get('root_count'),component_count=r.get('component_count',r.get('support_component_count')),
            cut_fragments=rem.get('fragment_count'),trapped_fragments=rem.get('trapped_pieces'),unresolved_warning_count=s.get('unresolved_warning_count'),print_time_seconds=s.get('estimated_print_time_seconds')))
    table(out/'tables/v2_vs_v3.csv',rows);return rows


def graphics(out,best,results,comparison):
    import pyvista as pv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .sacrificial_fixture_review import polydata
    from .support_inputs import transform
    folder=Path(best['folder']);qc=out/'QC';qc.mkdir(exist_ok=True)
    model=trimesh.load_mesh(folder/'mold_without_support.stl');v=trimesh.load_mesh(folder/'vascular_with_ports.stl');box=trimesh.load_mesh(folder/'box.stl');support=trimesh.load_mesh(folder/'custom_support_v3.stl')
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    with (folder/'frozen_detection_cache.pkl').open('rb') as f:det,field,_,intervals,bs=pickle.load(f)
    def shot(n,title,items,points=(),segments=(),view='iso',note=''):
        p=pv.Plotter(off_screen=True,window_size=(1400,1000));p.set_background('white')
        for mesh,color,alpha in items:
            if mesh is not None and len(mesh.faces):p.add_mesh(polydata(mesh),color=color,opacity=alpha,smooth_shading=False)
        for xyz,color,size in points:
            if len(xyz):p.add_points(np.asarray(xyz),color=color,point_size=size,render_points_as_spheres=True)
        for xyz,color,width in segments:
            if len(xyz)>1:p.add_lines(np.asarray(xyz),color=color,width=width,connected=True)
        p.add_text(title,font_size=17,color='black',position='upper_left');p.add_text(note,font_size=11,color='#333333',position='lower_left')
        p.add_axes();p.view_xy() if view=='top' else p.view_xz() if view=='side' else p.view_isometric()
        p.reset_camera();p.camera.zoom(.9);p.screenshot(str(qc/f'{n:02}_{STEMS[n-1]}.png'));p.close()
    def lines_from_polys(items,color):
        lines=[]
        for z,p in items:
            for q in (p.geoms if hasattr(p,'geoms') else [p]):
                if q.geom_type=='Polygon' and q.area>1e-7:
                    xy=np.array(q.exterior.coords);lines.append((np.c_[xy,np.full(len(xy),z)],color,3))
        return lines
    oldfolder=out.parent/'v2/candidates/O000';old=trimesh.load_mesh(oldfolder/'custom_support.stl');oldv=trimesh.load_mesh(oldfolder/'vascular_with_ports.stl');oldrecord=json.loads((oldfolder/'candidate.json').read_text())
    oldcov=oldrecord['slice']['actual_vascular_coverage']['percent'];oldself=oldrecord['self_buildability']['unsupported_regions'];oldtrapped=oldrecord['removal']['trapped_pieces']
    shot(1,'V2 O000: reference only',[(old,'#b79673',1),(oldv,'#747b84',.25)],note=f'{oldcov:.2f}% historical coverage; {oldself} self deficits; {oldtrapped} trapped pieces. Production remains empty.')
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,r in zip(axes.flat,results):
        path=Path(r['folder'])/'custom_support_v3.stl'
        if path.exists():
            mesh=trimesh.load_mesh(path);q=mesh.vertices[::max(1,len(mesh.vertices)//15000)];ax.scatter(q[:,0],q[:,2],s=.1,color='#b67f39')
        ax.set_title(r['candidate_id']+' / rigid scale 1');ax.set_aspect('equal');ax.set_xlabel('print X (mm)');ax.set_ylabel('print Z (mm)')
    fig.suptitle('8 bounded candidates; geometry failures remain diagnostic');fig.tight_layout();fig.savefig(qc/'02_v3_orientation_candidates.png',dpi=150);plt.close(fig)
    intervals_lines=[];long_lines=[]
    for x in intervals:
        b=bs[x['branch_index']];q=b['xyz'][x['start_index']:x['end_index']+1];intervals_lines.append((q,'#c64b3b',5))
        if x['classification']=='LONG_TUBULAR_SPAN':long_lines.append((q,'#246aac',6))
    shot(3,'Critical support intervals',[(v,'#cccccc',.4)],segments=intervals_lines,note=f'{len(intervals)} intervals; port/vessel semantics and complete layer polygons in CSV.')
    shot(4,'Long tubular intervals',[(v,'#cccccc',.4)],segments=long_lines,note='SWC arc association; the frozen detector is unchanged.')
    palette=['#e09b3d','#3f8f9c','#8260b2','#85a44a']
    shot(5,'Segmented vertical ribs',[(r['ribs'],palette[i%4],1) for i,r in enumerate(ribs)],note='Continuous narrow crown over arched windows; cuts separate the crown.')
    windows=[w['xyz'] for r in ribs for w in r['windows']]
    shot(6,'Removal windows: 3 mm clear reference',[(support,'#e0bd83',.65)],points=[(windows,'#009ec2',11)],note='Window width is a geometric probe reference, not a verified tool.')
    shot(7,'Narrow roof ridges',[(r['roofs'],'#bd4661',1) for r in ribs]+[(v,'#a9afb3',.2)],note='Lower envelope ridge, finite layerwise flare; no U-shaped wrapping support.')
    selfpts=[[*r['representative_xy'],r['z_mm']] for r in best.get('self_buildability',{}).get('rows',[])]
    actual_selfpts=[[*r['representative_xy'],r['z_mm']] for r in best.get('slice',{}).get('actual_self_buildability',{}).get('rows',[])]
    shot(8,f'Self-printability: {len(selfpts)} mesh / {len(actual_selfpts)} actual extrusion deficits',[(support,'#d8ac69',1)],points=[(selfpts,'#a333bb',10),(actual_selfpts,'red',9)],note='Purple=STL layers; red=actual extrusion layers. Bridge credit = 0; region tolerance 0.005 mm2.')
    shot(9,'Support only / top',[(support,'#bd8842',1)],view='top')
    shot(10,'Support only / side',[(support,'#bd8842',1)],view='side')
    # Local standalone curtain is explicitly labelled: it makes all five
    # construction regions legible instead of hiding them behind the full box.
    sample=max(ribs,key=lambda r:len(r['windows']))
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(sample['path_xy'],axis=0),axis=1))]
    ss=np.linspace(.001,arc[-1]-.001,240);zz=np.linspace(sample['base'].min()+.001,sample['top'].max()-.001,180)
    xx=np.interp(ss,arc,sample['path_xy'][:,0]);yy=np.interp(ss,arc,sample['path_xy'][:,1])
    q=np.c_[np.tile(xx,len(zz)),np.tile(yy,len(zz)),np.repeat(zz,len(ss))];labels=np.zeros(len(q),int)
    for k,key in enumerate(['ribs','roofs','necks','feet'],1):labels[sample[key].contains(q)]=k
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    fig,ax=plt.subplots(figsize=(12,6));colors=['white','#d79d45','#bb4264','#7846b8','#54777b']
    ax.pcolormesh(ss,zz,labels.reshape(len(zz),len(ss)),cmap=ListedColormap(colors),vmin=0,vmax=4,shading='nearest')
    ax.set_aspect('equal');ax.set_xlabel('Arc along actual curved rib (mm)');ax.set_ylabel('Print height (mm)')
    ax.set_title('Actual longitudinal centre section: '+sample['rib_id'])
    ax.legend(handles=[Patch(color=c,label=n) for c,n in zip(colors[1:],['Rib','Roof ridge','Cut neck','Root'])],loc='upper left',bbox_to_anchor=(1.01,1))
    for window in sample['windows']:ax.annotate('3 mm window',(window['arc_mm'],window['clear_shoulder_z_mm']-.7),ha='center',fontsize=9)
    fig.tight_layout();fig.savefig(qc/'11_support_only_section.png',dpi=170);plt.close(fig)
    shot(12,'Support plus frozen vascular core',[(support,'#d6a45b',.85),(v,'#50768f',1),(box,'#9fa5a8',.12)])
    cutrows=list(csv.DictReader((folder/'cut_points.csv').open()));cp=[json.loads(r['origin']) for r in cutrows]
    shot(13,'Geometric cut planes / explicit notches',[(support,'#d6a45b',.65)],points=[(cp,'#b83d6e',8)],note='Finite-plane audit checks other material, not just the isolated neck primitive.')
    tools=json.loads((folder/'tool_access.json').read_text()) if (folder/'tool_access.json').exists() else []
    paths=[]
    for r in tools:
        attempt=r if r['passed'] else r['attempts'][0]
        paths.append((transform(np.array([attempt['start'],attempt['end']]),best['orientation']['print_transform']),'#24894a' if r['passed'] else '#cc3939',2))
    shot(14,'3 mm probe access from casting top',[(support,'#d6a45b',.4),(v,'#667d91',.5)],segments=paths,note='Green=proxy path clear; red=blocked. Actual pliers/scissors NOT validated.')
    removal=best.get('removal',{});fragment_items=[];motions=[]
    for i,r in enumerate(removal.get('fragments',[])):
        mesh=trimesh.load_mesh(folder/'removal_fragments'/f"{r['fragment_id']}.stl");mesh.apply_transform(best['orientation']['print_transform'])
        fragment_items.append((mesh,'#ce3939' if r['status']=='TRAPPED_FRAGMENT' else palette[i%4],.8))
        if r.get('path'):motions.append((transform(np.array(r['path']['transforms'])[:,:3,3],best['orientation']['print_transform']),'#288b4d',1))
    shot(15,'Ordered removal: remaining fragments are obstacles',fragment_items,segments=motions,note=f"{removal.get('trapped_pieces')} trapped; casting-top exit checked on the whole piece.")
    covpath=folder/'bambu/actual_coverage_geometry.pkl';supported=[];missing=[]
    if covpath.exists():
        with covpath.open('rb') as f:supported,missing=pickle.load(f)
    shot(16,'ACTUAL TOOLPATH COVERAGE',[(v,'#aaaaaa',.2)],segments=lines_from_polys(supported,'#249148')+lines_from_polys(missing,'#d82929'),
         note='GREEN = actual extrusion support; RED = uncovered required polygons. No nominal roof colouring.')
    shot(17,'Uncovered actual layer polygons',[(v,'#bbbbbb',.2)],segments=lines_from_polys(missing,'#d82929'),note='Any visible red region above tolerance blocks production.')
    fig,axes=plt.subplots(1,4,figsize=(18,5));metrics=['actual_coverage_percent','self_unsupported','actual_self_unsupported_regions','trapped_fragments']
    for ax,k in zip(axes,metrics):
        vals=[r.get(k) or 0 for r in comparison];ax.bar([r['method'] for r in comparison],vals,color=['#999999','#bf934f']);ax.set_title(k)
        ax.bar_label(ax.containers[0],labels=['N/A' if r.get(k) is None else f'{r[k]:.2f}' for r in comparison])
    fig.suptitle('Matched coverage gate; NO contact-area reduction percentage claim');fig.tight_layout();fig.savefig(qc/'18_v2_vs_v3.png',dpi=150);plt.close(fig)
    warningpath=folder/'bambu/warning_localization.csv';warnings=list(csv.DictReader(warningpath.open())) if warningpath.exists() else []
    bad=[json.loads(r['xyz']) for r in warnings if r.get('xyz') and r['classification']!='INTENTIONAL_AIR_GAP']
    shot(19,'Unresolved warning localization',[(model,'#bbbbbb',.25)],points=[(bad,'#d32929',6)],note='Independent geometry / G-code correlation; native Bambu warnings do not provide these coordinates.')
    shot(20,'DIAGNOSTIC REVIEW / '+best['candidate_id'],[(box,'#bcc5c8',.2),(v,'#477286',1),(support,'#d29b45',1)],note='NEEDS_ADJUSTMENT; no physical coupon. See failed gates before considering a full print.')


def generate(out,results,profile,inputs,protection,runtime):
    valid=[r for r in results if r.get('planned_coverage')];best=min(valid,key=rank) if valid else None
    if not best:
        summary=dict(status='NEEDS_ADJUSTMENT',production_selected=None,results=results,source_protection=protection)
        save(out/'vascular_support_v3_summary.json',summary);return summary
    eligible=[r for r in valid if not gates(r)];selected=min(eligible,key=rank) if eligible else None
    status='READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED' if selected else 'NEEDS_ADJUSTMENT' if any(r.get('slice',{}).get('slicer_success') for r in valid) else 'SLICER_VALIDATION_PENDING'
    main_ids=json.loads((out/'top4_candidate_ids.json').read_text());sliced_ids=[r['candidate_id'] for r in results if r.get('slice',{}).get('slicer_success')]
    summary=dict(status=status,output=str(out),architecture='FDM_SEGMENTED_RIB_SUPPORT_V3',production_selected=selected,diagnostic_selected=best,
        source_protection=protection,orientation_count=len(results),top4_candidate_ids=main_ids,top4_actual_slices=sum(n in sliced_ids for n in main_ids),
        all_actual_sliced_candidates=len(sliced_ids),extra_diagnostic_slice_ids=[n for n in sliced_ids if n not in main_ids],
        physical_status='NO PHYSICAL COUPON AVAILABLE',unvalidated=['ROOF_GAP','BRIDGE_CAPABILITY','BREAK_FORCE','ABS_RIB_STRENGTH'],runtime_seconds=runtime)
    save(out/'vascular_support_v3_summary.json',summary)
    save(out/'sensitivity/line_selection_final.json',dict(diagnostic_candidate=best['candidate_id'],
        candidates=[dict(candidate_id=r['candidate_id'],geometry_parameters=r.get('geometry_parameters'),actual_upgrade=r.get('actual_thin_wall_upgrade'),
                         critical_features_preserved=r.get('slice',{}).get('all_features_preserved')) for r in results],physical_strength_validated=False))
    if selected:
        dst=out/'production_selected';dst.mkdir(exist_ok=True)
        for name in ['mold_without_support.stl','custom_support_v3.stl','mold_with_custom_support_v3.stl','print_transform_v3.json','casting_restore_transform_v3.json','support_removal_plan.json','X1C_project.3mf']:shutil.copy2(Path(selected['folder'])/name,dst/name)
    elif (out/'production_selected').exists():
        raise RuntimeError('STALE_PRODUCTION_SELECTION_REQUIRES_EXPLICIT_ARCHIVAL')
    folder=Path(best['folder']);geometry=out/'geometry';geometry.mkdir(exist_ok=True)
    table(out/'tables/v3_orientation_candidates.csv',[dict(orientation_id=r['candidate_id'],gravity_in_casting=r['orientation']['gravity_in_casting'],
        scale_factor=r.get('scale_factor'),planned_coverage_percent=r.get('planned_coverage',{}).get('percent'),actual_coverage_percent=r.get('slice',{}).get('actual_vascular_coverage',{}).get('percent'),
        geometric_self_unsupported=r.get('self_buildability',{}).get('unsupported_region_count'),actual_self_unsupported=r.get('slice',{}).get('actual_self_buildability',{}).get('unsupported_region_count'),
        trapped_fragments=r.get('removal',{}).get('trapped_pieces'),real_slice_success=r.get('slice',{}).get('slicer_success'),failed_gates=gates(r)) for r in results])
    for name in ['custom_support_v3.stl','mold_with_custom_support_v3.stl','mold_without_support.stl','support_ribs.stl','support_roofs.stl','support_pillars.stl','support_cut_necks.stl','support_exterior.stl','print_transform_v3.json','casting_restore_transform_v3.json']:
        if (folder/name).exists():shutil.copy2(folder/name,geometry/name)
    mappings={'v3_critical_intervals.csv':'critical_intervals.csv','v3_rib_segments.csv':'rib_segments.csv','v3_cut_points.csv':'cut_points.csv',
              'v3_fragments.csv':'fragments.csv','v3_removal_paths.csv':'removal_paths.csv','v3_actual_coverage.csv':'bambu/actual_coverage.csv','v3_bambu_warnings.csv':'bambu/warning_localization.csv'}
    for dst,src in mappings.items():
        if (folder/src).exists():shutil.copy2(folder/src,out/'tables'/dst)
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:rib_data=pickle.load(f)
    segments=[]
    for rib in rib_data:
        length=float(np.linalg.norm(np.diff(rib['path_xy'],axis=0),axis=1).sum());start=0.
        stops=[(w['arc_mm']-w['clear_width_mm']/2,w['arc_mm']+w['clear_width_mm']/2) for w in rib['windows']]+[(length,length)]
        for index,(stop,following) in enumerate(stops):
            segments.append(dict(segment_id=f"{rib['rib_id']}_LEG{index:02}",parent_rib_id=rib['rib_id'],interval_id=rib['interval_id'],branch_id=rib['branch_id'],
                projected_arc_start_mm=start,projected_arc_end_mm=stop,leg_length_mm=max(0.,stop-start),width_mm=rib['width_mm'],roof_width_mm=rib['roof_width_mm'],
                is_port=rib['is_port'],definition='Vertical leg between clear tool windows; crown cut and removal fragments audited separately'))
            start=following
    table(out/'tables/v3_rib_segments.csv',segments)
    table(out/'tables/v3_roof_segments.csv',[dict(roof_id=r['rib_id'],width_mm=r['roof_width_mm'],interval_id=r['interval_id']) for r in json.loads((folder/'rib_segments.json').read_text())])
    comparison=compare(out,best,profile);s=best.get('slice',{});rem=best.get('removal',{});cov=s.get('actual_vascular_coverage',{});failed=gates(best)
    answers=[f"最终姿态：{best['candidate_id']}，打印 +Z 在浇注坐标中为 {best['orientation']['gravity_in_casting']}；当前是诊断候选，生产选择{'已建立' if selected else '为空'}。",
        '选择理由：先检查完整覆盖、自打印、碰撞、拆除及实际切片门槛；失败候选优先按覆盖排序，保留盒体欠支撑惩罚。本轮完整规划 8 个姿态，并未重复对 128 个姿态作昂贵建模。',
        f"必须支撑区间：{best['critical_intervals']} 个，含血管与端口；每个区间保留弧长、曲率、半径、重力关系及对应层多边形。",
        f"实际 Bambu 刀路覆盖：{cov.get('percent',0):.4f}%；分母为冻结检测器的需支撑层多边形面积之和。",
        f"漏支撑：{cov.get('uncovered_region_count')} 个层区域未通过，每区域数值容差 0.00001 mm²；不以接近 100% 代替完整覆盖。",
        f"薄肋段数：{best['rib_segments']}，工具窗 {best['window_count']}。",
        f"窄屋脊数：{best['roof_ridges']}，不采用包围血管的 U 形托架。",
        f"辅助柱数：{best['pillar_count']}；主要结构是分段薄肋。",
        '未校准桥：主候选采用 NO_UNCALIBRATED_BRIDGE；工具窗上方为逐层闭合拱，不使用 4/6/8 mm 自由桥。',
        f"支撑自身欠支撑区域：STL 几何逐层检查 {best['self_buildability']['unsupported_region_count']} 个，实际挤出足迹逐层检查 {s.get('actual_self_buildability',{}).get('unsupported_region_count')} 个。任一检查非零都不通过；没有桥接豁免。",
        f"可动碎片：{rem.get('fragment_count')} 块，其中内部 {rem.get('internal_fragment_count')} 块；另有明确设计的盒面根部残留 {rem.get('retained_root_feet')} 块。",
        f"受困碎片：{rem.get('trapped_pieces')} 块。尚未找到路径表示本次有限路径搜索未通过，不冒充数学上的不可拆除证明。",
        f"顶部 3 mm 工具探针：{rem.get('tool_access_failed')} 个切割点未通过；Tool access is a geometric proxy only。",
        f"Bambu Auto Support：{'OFF' if s.get('auto_support_off') else '未确认'}；原生支撑刀路 {s.get('native_support_toolpaths')} 条。",
        f"真实 X1C 切片：整件 {s.get('slicer_success')}，支撑加合法盒体基础 {s.get('support_only_slicer_success')}；切片成功与生产合格分别判定。",
        f"相对 V2 的改进：由大量离散顶垫/分支转为沿区间的薄肋。V2 的几何逐层自身欠支撑为 42 个，本次几何检查为 {best['self_buildability']['unsupported_region_count']} 个；但新增实际刀路检查仍有 {s.get('actual_self_buildability',{}).get('unsupported_region_count')} 个，不能宣布自打印目标已经实现。V2 12 个受困件对应本次 {rem.get('trapped_pieces')} 个。两者未同时通过完整覆盖，禁止声称接触面积降低百分比。",
        '未实验验证的参数：屋脊间隙、桥接能力、断裂/切割力、ABS 薄肋强度。NO PHYSICAL COUPON AVAILABLE；上述参数全部 NOT PHYSICALLY VALIDATED。']
    text=['# BG001 RMCA 分段薄肋支撑 V3 实现与验证','',f'状态：**{status}**。','']+[f'{i+1}. {a}' for i,a in enumerate(answers)]
    text+=['','## 方法与结果解释','',
        f"完整规划姿态 {len(results)} 个；当前主切片组为 {', '.join(main_ids)}，全部真实切片成功的候选共 {len(sliced_ids)} 个。额外诊断切片 {', '.join(n for n in sliced_ids if n not in main_ids) or '无'} 仍保留，并未覆盖或删除。",
        '本次保留 V1/V2 作为历史和基准，不改变已冻结的血管半径、端口、2 mm 盒体或可视化界面。检测沿用三角面、逐层截面、FreeFloating 和 SWC 弧长融合。新的规划单位为连续危险区间；根部与工具窗将其分成短薄肋，顶部以有限斜率过渡到窄屋脊。其几何假设是竖直薄壁比长细斜杆更少依赖悬空逐层生长，但这不是 ABS 实验强度结论。',
        '切割截面以有限工具平面分别核算已设计窄颈与其余支撑材料。根部下半颈可作为明确记录的盒面残留；不得仅因碎片很小就排除。拆除恢复到冻结浇注坐标，完整碎片、血管、盒体和未拆碎片都参加 FCL 扫掠检查。工具探针是独立检查。',
        '真实覆盖从 Bambu G-code 普通模型挤出线提取，在对应高度与显式支撑截面相交后，按实际配置线宽形成挤出足迹，再投影到危险多边形。血管本体和盒体刀路不会被整体视为支撑。QC 16 中绿色代表实际承托，红色代表缺失。',
        'V2 历史 88.00% 使用了短桥判据。本轮在 V3 validation 目录重新计算 V2 的无桥、严格容差覆盖，不覆盖其原始结果。v2_vs_v3.csv 同时保留历史指标的定义。只有双方同一覆盖门槛通过后，才允许进一步声明支撑负担改善比例。',
        '',f"失败门槛：{', '.join(failed) if failed else '无软件门槛失败'}。",'',
        f"保护复核：{protection['checked_files']} 个上游、历史输出和冻结检测文件 SHA256 {'全部不变' if protection['all_unchanged'] else '异常'}。V2 production_selected 保持原值。",'',
        '当前建议：'+('可进入人工几何复核；物理参数仍须审慎试验。' if selected else '不建议立即进行下一次完整模型试打。先修复所列软件门槛；当前 STL/3MF 仅供诊断审阅。'),
        '未上传、发送或启动打印机。']
    testfile=out/'test_results.xml'
    if testfile.exists():
        import xml.etree.ElementTree as ET
        suites=list(ET.parse(testfile).getroot().iter('testsuite'));count=sum(int(x.get('tests',0)) for x in suites);failures=sum(int(x.get('failures',0))+int(x.get('errors',0)) for x in suites)
        text+=['',f'永久测试：{count} 项，失败或错误 {failures} 项；原有 V1/V2 测试保留，新增 13 类几何和回归测试。']
    (out/'vascular_support_v3_report.md').write_text('\n'.join(text)+'\n')
    graphics(out,best,results,comparison)
    return summary
