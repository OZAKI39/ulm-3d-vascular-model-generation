"""Readable Chinese report and non-interactive review figures for every outcome."""
from pathlib import Path
import csv,json,shutil
import numpy as np
import trimesh
from .support_aware_review import scene,annotate
from .support_inputs import save,table

def markers(rows,color='#e6a02e',radius=.35):
    shapes=[]
    for r in rows:
        p=np.asarray(r['xyz'] if isinstance(r,dict) else r);m=trimesh.creation.icosphere(subdivisions=1,radius=radius);m.apply_translation(p);shapes.append(m)
    return [(trimesh.util.concatenate(shapes),color,1.)] if shapes else []

def path_mesh(path):
    """Render real extrusion segments as lines via PyVista, not invented solids."""
    import pyvista as pv
    segments=np.load(path)['segments'];data=pv.PolyData(segments.reshape(-1,3));data.lines=np.column_stack([np.full(len(segments),2),np.arange(len(segments)*2).reshape(-1,2)]).ravel();return data

def slicer_view(core,box,folder,path,title):
    from .sacrificial_fixture_review import polydata
    import pyvista as pv
    from PIL import Image
    p=pv.Plotter(off_screen=True,window_size=[1100,800]);p.set_background('white');p.add_mesh(polydata(core),color='#3b9969');p.add_mesh(polydata(box),color='#849aaa',opacity=.13)
    f=Path(folder)/'actual_support_paths.npz'
    if f.exists():p.add_mesh(path_mesh(f),color='#cb9242',line_width=1)
    p.view_isometric();p.camera.zoom(1.1);p.show(screenshot=str(path),auto_close=True)
    annotate(Image.open(path).convert('RGB'),title,['黄色：实际切片支撑挤出刀路；绿色：冻结血管']).save(path)

def arc_coverage_view(out,det,guided,path):
    """Show every sampled branch, including safe uncontacted intervals."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    from PIL import Image
    from .support_contact_planner import influence
    branches=json.loads((out/'branches.json').read_text());cfg=json.loads((out/'resolved_config.json').read_text())
    actual=Path(guided['folder'])/'actual_contact_coverage.csv' if guided else None
    ids={r['point_id'] for r in csv.DictReader(actual.open()) if r['actual_support_underneath']=='True'} if actual and actual.exists() else set()
    touched=[p for p in det['raw'] if p['point_id'] in ids]
    missing=set((guided.get('coverage') or {}).get('uncovered_point_ids',[])) if guided else {p['point_id'] for p in det['raw'] if p['risk_level']=='CRITICAL'}
    fig,ax=plt.subplots(figsize=(11,6.5));labels=[];rows=[]
    for j,b in enumerate(branches):
        raw={p['sample_index']:p for p in det['raw'] if p['branch_id']==b['branch_id']};hits=[p for p in touched if p['branch_id']==b['branch_id']]
        s=np.array(b['arc']);colors=[]
        for i,arc in enumerate(s):
            p=raw.get(i);supported=any(abs(q['arc_position']-arc)<=influence(q,cfg) for q in hits)
            status='SUPPORTED' if supported else 'UNSUPPORTED_SAFE'
            if p and p['risk_level']=='CRITICAL':status='UNSUPPORTED_CRITICAL' if p['point_id'] in missing else 'SUPPORTED'
            colors.append({'SUPPORTED':'#399d63','UNSUPPORTED_SAFE':'#aeb9c3','UNSUPPORTED_CRITICAL':'#ce4f4f'}[status])
            rows.append(dict(branch_id=b['branch_id'],arc_position_mm=float(arc),status=status))
        segments=[[(s[i],j),(s[i+1],j)] for i in range(len(s)-1)]
        ax.add_collection(LineCollection(segments,colors=colors[:-1],linewidths=10))
        labels.append(b['branch_id'])
    ax.set(xlim=(0,max(max(b['arc']) for b in branches)+1),ylim=(-.7,len(branches)-.3),xlabel='Centerline arc length (mm)',yticks=range(len(branches)),yticklabels=labels)
    ax.invert_yaxis();ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.2)
    ax.legend(handles=[Line2D([0],[0],color=c,lw=7,label=l) for c,l in [('#399d63','Supported arc'),('#aeb9c3','Unsupported, noncritical'),('#ce4f4f','Unsupported, critical')]],loc='upper center',bbox_to_anchor=(.5,1.13),ncol=3)
    fig.tight_layout();fig.savefig(path,dpi=150);plt.close(fig)
    annotate(Image.open(path).convert('RGB'),'逐分支弧长覆盖：基于实际切片支撑位置',['绿色：接触位置的工程影响范围；灰色：非危险未支撑段；红色：危险遗漏','人工端口与解剖分支均列出；绿色不表示整段表面与支撑连续接触']).save(path)
    table(out/'tables/actual_branch_arc_coverage.csv',rows)
    spans=[]
    for region in det['regions']:
        candidates=[p for p in det['raw'] if p['branch_id']==region['branch'] and region['arc_start']<=p['arc_position']<=region['arc_end'] and p['risk_level']=='CRITICAL']
        spans.append({**region,'coverage_status':'UNCOVERED' if any(p['point_id'] in missing for p in candidates) else 'COVERED' if candidates else 'NONCRITICAL'})
    table(out/'tables/unsupported_spans.csv',spans)

def access_view(mesh,box,guided,path):
    import pyvista as pv
    from PIL import Image
    from .sacrificial_fixture_review import polydata
    p=pv.Plotter(off_screen=True,window_size=[1100,800]);p.set_background('white')
    p.add_mesh(polydata(mesh),color='#3b9969');p.add_mesh(polydata(box),color='#849aaa',opacity=.10)
    record=Path(guided['folder'])/'support_accessibility.json' if guided else None
    rows=json.loads(record.read_text())['rows'] if record and record.exists() else []
    if guided and (Path(guided['folder'])/'actual_support_paths.npz').exists():p.add_mesh(path_mesh(Path(guided['folder'])/'actual_support_paths.npz'),color='#cb9242',line_width=1)
    for row in rows:
        route=row.get('probe_path')
        if route:
            line=pv.Line(route['start'],route['end']);p.add_mesh(line.tube(radius=1.5),color='#48afc7',opacity=.28)
            p.add_mesh(line,color='#246782',line_width=2)
    p.view_isometric();p.camera.zoom(1.05);p.show(screenshot=str(path),auto_close=True)
    annotate(Image.open(path).convert('RGB'),'真实支撑分量的顶部 3 mm 探针路径',[f"{sum(r['probe_3mm_reachable'] for r in rows)}/{len(rows)} 个采样连通分量可接近；蓝色透明管为通过碰撞检查的探针路径",'可接近只表示能够触达抓取位置，不保证完整支撑可以无损抽出']).save(path)

def figures(out,mesh,box,det,contacts,topology,records,selected):
    qc=out/'QC';qc.mkdir(exist_ok=True);base=[(mesh,'#3b9969',1),(box,'#8ca4b5',.12)]
    scene(base,qc/'01_input_vascular_model.png','冻结输入：血管、端口与 2 mm 盒体')
    folder=records['BAMBU_NATIVE_BASELINE']['folder'];slicer_view(mesh,box,folder,qc/'02_native_bambu_support_problem.png','Bambu 原生支撑：实际刀路')
    data=np.load(out/'support_requirement_surface.npz');risk=data['required'].astype(int)
    from scipy.spatial import cKDTree
    critical=[p['xyz'] for p in det['raw'] if p['risk_level']=='CRITICAL']
    if critical:risk[(cKDTree(critical).query(mesh.triangles_center)[0]<.7)&data['required']]=2
    scene([(mesh,'#3b9969',1)],qc/'03_support_requirement_surface.png','支撑需求：绿色安全／橙色注意／红色危险',risk=risk)
    points=np.load(out/'detector_points.npz')
    scene(base+markers(points['freefloating']),qc/'04_freefloating_detection.png','FreeFloating 检测点（MIT CPU 移植）')
    scene(base+markers(critical,'#cc4545'),qc/'05_swc_long_span_risk.png',f"危险长悬空段：{det['critical_region_count']} 段，最长 {det['longest_dangerous_span_mm']:.2f} mm")
    scene(base+markers(det['raw'],'#e4a035',.18),qc/'06_all_candidate_support_points.png',f"筛选前：{len(det['raw'])} 个候选位置")
    scene(base+markers(contacts,'#d84d43',.4),qc/'07_pruned_support_points.png',f'筛选后：{len(contacts)} 个计划接触点')
    keep=[]
    for p in topology['bifurcations']:
        m=trimesh.creation.icosphere(subdivisions=2,radius=2.);m.apply_translation(p);keep.append((m,'#b743b3',.18))
    guided=records.get('GUIDED_ORCA')
    maskfile=Path(guided['folder'])/'contact_face_mask.npy' if guided else None
    face_risk=np.zeros(len(mesh.faces),int)
    if maskfile and maskfile.exists():
        touched=np.load(maskfile);face_risk[touched]=1
        if topology['bifurcations']:face_risk[touched & (cKDTree(topology['bifurcations']).query(mesh.triangles_center)[0]<2.)]=2
    scene([(mesh,'#3b9969',1)]+keep,qc/'08_bifurcation_keepout.png','实际接触面与分叉 2 mm 避让区（从下方观察）',risk=face_risk,
        lines=['橙色：实际接触面；红色：进入分叉邻域的接触面；紫色：避让区域',f"邻域内实际接触区间：{guided.get('support_near_bifurcation_count') if guided else '尚未核验'}；属于明确保留的覆盖例外"])
    arc_coverage_view(out,det,guided,qc/'09_support_coverage.png')
    if guided:slicer_view(mesh,box,guided['folder'],qc/'10_guided_orca_support.png','Guided Orca：实际支撑刀路')
    else:scene(base,qc/'10_guided_orca_support.png','Guided Orca：切片结果尚未可用')
    custom=out/'custom_clever_support.stl';scene(base+([(trimesh.load_mesh(custom),'#d69738',1)] if custom.exists() else []),qc/'11_custom_clever_tree.png','自定义树支撑：研究候选，不作为生产推荐')
    bridges=sorted((out/'research').glob('bridge*.stl'))
    scene(base+([(trimesh.load_mesh(bridges[0]),'#d69738',1)] if bridges else []),qc/'12_bridge_pillar_candidate.png','立柱与短桥：实验方案')
    access_view(mesh,box,guided,qc/'13_support_accessibility.png')
    from PIL import Image
    a=Image.open(qc/'02_native_bambu_support_problem.png');b=Image.open(qc/'10_guided_orca_support.png');im=Image.new('RGB',(a.width+b.width,max(a.height,b.height)),'white');im.paste(a,(0,0));im.paste(b,(a.width,0));im.save(qc/'14_baseline_vs_guided.png')
    chosen=records.get(selected or '',guided or records['BAMBU_NATIVE_BASELINE'])
    slicer_view(mesh,box,chosen['folder'],qc/'15_selected_support_plan.png','人工审核候选：'+(selected or '没有通过全部硬门槛的方案'))
    native=Path(chosen['folder'])/'native_preview.png'
    if native.exists():shutil.copy2(native,qc/'16_real_slicer_preview.png')
    else:slicer_view(mesh,box,chosen['folder'],qc/'16_real_slicer_preview.png','实际切片刀路预览；引擎未导出原生缩略图')

def report(out,summary):
    d=summary['detection'];comparison=summary['baseline_comparison'];base=next((r for r in comparison if r['name']=='BAMBU_NATIVE_BASELINE'),{});guided=next((r for r in comparison if r['name']=='GUIDED_ORCA'),{})
    def fmt(v):return '尚未取得有效结果' if v is None else f'{v:.2f}' if isinstance(v,float) else str(v)
    lines=['# BG001 RMCA 血管专用支撑规划报告','',f"当前状态：**{summary['status']}**。本轮输入为最新 2 mm 五面盒与既有血管、端口；上游模型 SHA256 复核结果：{summary['source_protection']['all_unchanged']}。",'',
        '## 为什么原支撑需要调整','',
        '细血管的长悬空段可能在自动支撑中被遗漏；增加自动支撑密度又可能让树枝沿血管连续包裹。本轮将接触位置与支撑结构分开处理：先确定危险位置，再让成熟切片器生成支撑。真实 ABS 打印和机械拆除仍是后续实验。','',
        '## 检测与接触位置','',
        f"检测到 **{d['critical_regions']} 段危险悬空区域**，最长 **{d['longest_span_mm']:.2f} mm**。{d['raw_points']} 个候选位置被筛选为 {d['selected_points']} 个计划接触点。5 mm 是本次工程阈值，并非通用物理定律；4、5、6 mm 的敏感性结果另存为表。",'',
        f"Bambu 原生切片约有 **{fmt(base.get('support_tip_count'))} 个弧长接触区间**；Guided Orca 实际约有 **{fmt(guided.get('support_tip_count'))} 个**。此处是从挤出刀路投影到血管中心线后估计的断续接触区间，不等于切片器内部树叶节点数；连续包裹和网格接触片数量另外记录，避免混淆。",'',
        '按这一共同估计口径，本次**未达成接触数量下降的目标**。筛选后的计划点数不能替代实际切片的接触区间数；本次改进主要体现在减少连续接触面积，不能宣称所有拆除指标均已改善。','',
        '## 必要位置是否仍有遗漏','',
        f"Guided Orca 是否仍有危险段没有覆盖：**{'NO' if guided.get('coverage_percent')==100 else 'YES'}**。实际覆盖率 {fmt(guided.get('coverage_percent'))}%，最大未覆盖危险弧长 {fmt(guided.get('maximum_unsupported_mm'))} mm。计划覆盖与刀路覆盖分开计算；仅计划通过不能作为生产通过。",'',
        '## 接触与拆除负担','',
        f"Guided Orca 的接触面积相对 Bambu 变化：减少 {fmt(guided.get('contact_area_reduction_percent'))}%。该面积由下表面三角形与实际支撑挤出轨迹的邻近关系估计；所有方案使用同一估计方法。",'',
        f"分叉 2 mm 避让区附近接触区间数：{fmt(guided.get('bifurcation_contacts'))}。这表示靠近分叉，不自动等于支撑与分叉实体相交。预计被困支撑数：{fmt(guided.get('trapped_support_count'))}。探针可达性是离散几何估计，能抓到树干不等于一定能完整抽出。",'',
        '分叉完全避让尚未实现。计划中进入避让区的点逐项标记为覆盖例外：移除后会漏掉危险要求，保留理由见 support_contacts.csv。真实切片接触范围可能大于计划点，因此另行核验实际接触面；QC/08 从下方着色显示实际接触，不能把“碰撞检查通过”解读成“分叉没有接触”。','',
        '## 同条件比较','',
        '|方案|接触区间数（估计）|接触面积 mm²|危险覆盖率 %|被困风险数|切片时间估计 s|硬门槛失败|','|---|---:|---:|---:|---:|---:|---|']
    for r in comparison:lines.append('|'+ '|'.join([r['name'],fmt(r.get('support_tip_count')),fmt(r.get('contact_area_mm2')),fmt(r.get('coverage_percent')),fmt(r.get('trapped_support_count')),fmt(r.get('estimated_print_time_seconds')),', '.join(r.get('hard_gate_failures',[])) or '无'])+'|')
    lines+=['','表中的 SUPPORT_COLLISION_NOT_CLEARED 表示保守挤出包络的间距下界未通过，不等同于已经证明发生实体碰撞；失败候选均未被推荐。接触面积与支撑质量均为刀路估计。自定义树的名义接触圆盘面积不可与这些实际刀路面积直接比较，所以其实际接触面积留空。','',
        'Orca 使用其 X1C 0.2 mm 官方配置结构承接本机保存的 Vessel 工艺值。层高、首层、喷嘴、线宽、顶部和 XY 间隙均逐项保持；不同引擎的启动宏使用各自官方实现，转换明细见 profiles/orca_translation_audit.json。','',
        '## 下一次打印选择与剩余风险','',f"通过全部软件硬门槛后的排序推荐：**{summary.get('selected_candidate') or '暂无；仍需调整'}**。Guided Orca 项目文件可人工审核，STL 中的 enforcer/blocker 是控制区域，不是直接打印实体。自定义 Clever Tree 和短桥均为研究候选，其失败结果也保留。",'',
        '最大剩余风险是软件中可支撑、可触达的细血管，在 ABS 冷却收缩和真实拆支撑受力时仍可能断裂。没有执行真实打印、上传或启动打印，也不声明制造验证完成。','',
        f"本次推荐参数：{json.dumps(summary['support'].get('chosen_parameters'),ensure_ascii=False)}。具体候选目录：`{summary['support'].get('artifact_directory')}`。支撑估计用料 {fmt(guided.get('estimated_support_mass_g'))} g；整件打印时间估计 {fmt(guided.get('estimated_print_time_seconds'))} s，并非本次计算耗时。",'',
        '## 方法与复现','',
        '四种检测证据分别来自表面法向、按实际层高切片的上下层承托、FreeFloating 图覆盖与局部低点、SWC/已冻结端口路径的弧长分析。Greedy 覆盖选择沿同一分支弧长计算，弯曲度和局部半径会缩短影响范围；不是固定欧氏球。各指标分别输出，候选先过硬门槛，再按接触面积等字典序排序。','',
        '`./run_python.sh s1-8_vascular_support_planner.py --stage all`','',
        '输入来源见 source_discovery.json；检测、接触、覆盖见 tables/ 和 JSON；12 个以内 Orca 搜索候选见 slicer/；阶段图见 QC/；输入哈希见 protected_before.json 和 source_protection.json；永久测试位于 tests/test_s1_8_vascular_support_planner.py。']
    (out/'vascular_support_report.md').write_text('\n'.join(lines)+'\n')
