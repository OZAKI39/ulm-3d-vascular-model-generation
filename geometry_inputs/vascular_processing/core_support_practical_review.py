"""Final practical evidence, with unsuccessful cuts kept explicitly visible."""
from pathlib import Path
import json,pickle,shutil
import numpy as np
import trimesh
from .support_inputs import save

NAMES=['01_baseline_support','02_true_missing_vs_gap_quality','03_true_missing_repairs','04_final_structural_support_coverage','05_macro_support_printability','06_support_components_before_split','07_strategic_cut_seams','08_support_components_after_split','09_major_component_contact_map','10_last_contact_pieces','11_bifurcation_contact_pieces','12_wrap_risk_map','13_sequential_removal_step_1','14_sequential_removal_step_2','15_sequential_removal_step_final','16_final_support_only','17_final_support_plus_vascular','18_final_Bambu_preview']


def render(base,out,folder,summary,structural,macro,removal):
    import pyvista as pv
    import shapely
    out=Path(out);folder=Path(folder);qc=out/'QC_FINAL';qc.mkdir(exist_ok=True)
    g=pickle.load((folder/'geometry.pkl').open('rb'));original=pickle.load((base/'candidates/O002_R25_C2_LAYER_BASE/geometry.pkl').open('rb'));v=g['vascular'];support=g['support'];colors=['#4682b4','#ff9933','#55aa66','#dd6688','#9365b8','#aa8844','#55bbbb','#999966','#cc7755','#7777dd','#bbaa66','#cc88cc','#66aacc','#aacc88']
    def poly(m):return pv.PolyData(m.vertices,np.c_[np.full(len(m.faces),3),m.faces].ravel())
    def shot(index,subtitle='',meshes=(),points=(),labels=(),core=True,lines=()):
        plot=pv.Plotter(off_screen=True,window_size=(1500,1050));plot.set_background('#f6f8fa')
        if core:plot.add_mesh(poly(v),color='#b75549',opacity=.55,smooth_shading=True)
        for mesh,color,alpha in meshes:
            if len(mesh.faces):plot.add_mesh(poly(mesh),color=color,opacity=alpha)
        for xyz,color,size in points:
            if len(xyz):plot.add_points(np.asarray(xyz),color=color,point_size=size,render_points_as_spheres=True)
        if labels:plot.add_point_labels(np.array([q for q,label in labels]),[label for q,label in labels],font_size=11,shape_opacity=.3,always_visible=True)
        if lines:plot.add_lines(np.asarray(lines).reshape(-1,3),color='#333333',width=3)
        plot.add_text(NAMES[index]+'\n'+subtitle+'\n'+summary['status'],font_size=11,color='#223344');plot.add_axes();plot.show_bounds(xtitle='X (mm)',ytitle='Y (mm)',ztitle='Z (mm)',fmt='%.1f',n_xlabels=3,n_ylabels=3,n_zlabels=3,font_size=10);plot.view_isometric();plot.reset_camera();plot.camera.zoom(.9);plot.screenshot(qc/(NAMES[index]+'.png'));plot.close()
    scaffold=[(support,'#77aabb',.45)]
    shot(0,'Existing O002 baseline: original support',[(original['support'],'#77aabb',.65)])
    baseline=json.loads((out/'baseline_audits/O002_R25_C2_LAYER_BASE/final_true_missing_support.json').read_text())
    def classified(rows):
        return [([r['xyz'] for r in rows if r['true_missing_support']],'#ed2020',15),([r['xyz'] for r in rows if r['new_class']=='SUPPORT_TOO_CLOSE'],'#efba22',12),([r['xyz'] for r in rows if r['numerical_exception']],'#9944dd',15),([r['xyz'] for r in rows if r['new_class']=='ACCEPTABLE'],'#36a65a',5)]
    shot(1,'Baseline: RED real missing / YELLOW close gap / PURPLE numerical',scaffold,classified(baseline['rows']))
    repairs=[]
    for file in sorted((out/'candidates').glob('PRACTICAL_REPAIR_R*/true_missing_repairs.json')):repairs+=json.loads(file.read_text())['changes']
    shot(2,'Local repair targets only; red point remains unresolved',scaffold,[([r['xyz'] for r in repairs],'#2388dd',12),([r['xyz'] for r in structural['rows'] if r['true_missing_support']],'#ed2020',18)])
    shot(3,f"Actual extrusion: {structural['true_missing_support_count']} true missing; {structural['gap_quality_counts']['TOO_CLOSE']} too-close",scaffold,classified(structural['rows']),labels=[(r['xyz'],r['sample_id']) for r in structural['rows'] if r['true_missing_support'] or r['numerical_exception']])
    macropts=[];micropts=[]
    for r in macro['rows']:
        if r.get('wkt'):
            q=shapely.from_wkt(r['wkt']).representative_point();(macropts if r['status']=='MACRO_FAIL' else micropts).append([q.x,q.y,r['z_mm']])
    shot(4,f"Actual support: {macro['macro_fail_count']} macro failures / {macro['micro_deficit_count']} attached micro deficits",scaffold,[(macropts,'#ff2222',13),(micropts,'#e6aa22',9)],core=False)
    shot(5,f"Before post-print cuts: {removal['initial_components']} material component",[(support,colors[0],.9)],core=False)
    cutlabels=[(r['bulk_cut_grip_point'],r['cut_id']) for r in removal['cuts'] if r['effective']]
    shot(6,'Strategic finite blade envelopes; access still audited separately',scaffold,labels=cutlabels)
    directory=folder/'practical_removal';meshes={r['component_id']:trimesh.load(directory/'major_components'/(r['component_id']+'.stl'),force='mesh') for r in removal['major_components']}
    major=[(meshes[r['component_id']],colors[i%len(colors)],.9) for i,r in enumerate(removal['major_components'])]
    labels=[(meshes[r['component_id']].centroid,r['component_id']) for r in removal['major_components']]
    shot(7,f"Virtual cuts produce {len(major)} solids; not all are safely removable",major,labels=labels,core=False)
    shot(8,'Labels show contact arc, not validated break force',major,labels=[(meshes[r['component_id']].centroid,f"{r['component_id']}: {r['contact_arc_mm']:.1f} mm") for r in removal['major_components']])
    residual=[];bif=[]
    for i,r in enumerate(removal['final_contact_pieces']):
        mesh=trimesh.load(directory/'contact_residues'/(r['piece_id']+'.stl'),force='mesh');item=(mesh,colors[i%len(colors)],1.);residual.append(item)
        if r['near_bifurcation']:bif.append(item)
    shot(9,'Roof residue only: long connected pieces remain; bulk cut NOT certified',residual,core=False)
    shot(10,'Bifurcation residue: smaller-piece priority, no forced lateral pull',bif)
    wr=json.loads((directory/'wrap_audit.json').read_text());bs={b['branch_id']:b for b in g['branches']};pts=[]
    for low,high,color in [(0,120,'#44aa66'),(120,150,'#eebb33'),(150,180,'#ee8833'),(180,1e6,'#dd2222')]:pts.append(([bs[r['branch_id']]['xyz'][r['sample_index']] for r in wr['rows'] if low<=r['angular_coverage_deg']<high],color,7))
    shot(11,'Wrap: green <120 / yellow 120-150 / orange 150-180 / red >180',scaffold,pts)
    steps=removal['sequence']['steps'];states=[min(1,len(steps)),min(2,len(steps)),len(steps)]
    for idx,n in zip([12,13,14],states):
        removed={r['component_id'] for r in steps[:n] if r['action']=='REMOVE_SHORT_COMPONENT'};visible=[(m,c,.65) for r,(m,c,_) in zip(removal['major_components'],major) if r['component_id'] not in removed];done=steps[n-1] if n else None
        pts=[];lines=[];lab=[]
        if done:
            point=np.array(done.get('point',done.get('grip_point')));d=np.array(done.get('access_direction',done.get('direction')));pts=[([point],'#222222',16)];lines=[[point,point+6*d]];lab=[(point,done.get('cut_id',done.get('component_id'))+' '+done['action'])]
        title=f'Diagnostic ordered actions: step {n}/{len(steps)}'
        if idx==14:title+=f"; BLOCKED {len(removal['sequence']['remaining_components'])} components"
        # Major solids shown separately for readability. Uncut kerfs are
        # restored in this figure so future cuts are not visually empty air.
        cutdone={r['cut_id'] for r in steps[:n] if r['action']=='CUT_STRATEGIC_CONNECTOR'}
        visible += [(trimesh.load(p,force='mesh'),'#555555',1.) for p in (directory/'cut_kerfs').glob('*.stl') if p.stem not in cutdone]
        shot(idx,title,visible,pts,lab,lines=lines)
    shot(15,'Selected diagnostic support, before post-print cuts',[(support,'#6699aa',1.)],core=False)
    shot(16,'Frozen vascular core + separate custom support',scaffold)
    shutil.copy2(folder/'bambu/native_preview.png',qc/(NAMES[17]+'.png'))
    save(out/'reports/QC_FINAL_manifest.json',dict(figures=[str(qc/(name+'.png')) for name in NAMES],all_exist=all((qc/(name+'.png')).exists() for name in NAMES),native_preview_source=str(folder/'bambu/native_preview.png'),diagnostic=summary['status']=='NEEDS_ADJUSTMENT'))


def report(out,s,structural,removal):
    out=Path(out);a=s['practical_audit'];m=s['macro'];seq=removal['sequence'];missing=[r for r in structural['rows'] if r['true_missing_support']]
    answers=[('本轮是否继续追求数学意义上的 100%？','没有。结构承托与间隙质量已分开。TOO_CLOSE 只影响拆卸和表面风险，不再当作缺承托；过大的间隙仍不能提供实际承托。'),
    ('真正缺承托还剩多少？',f"{structural['true_missing_support_count']} 个，逐点证据见 final_true_missing_support.json。"),
    ('最长真正裸悬空弧长？',f"{structural['max_true_unsupported_arc_mm']:.6f} mm。这是用相邻已承托采样点界定的保守区间，不是假定该长度已经获得桥接验证。"),
    ('多少只是支撑过近？',f"{structural['gap_quality_counts']['TOO_CLOSE']} 个；合适 {structural['gap_quality_counts']['ACCEPTABLE']} 个；过远 {structural['gap_quality_counts']['TOO_FAR']} 个。"),
    ('有没有浮空血管区域？',f"实际血管挤出连通区域起始检查中，未承托起始区域为 {a['floating_vascular_regions_unresolved']} 个；中心线下表面检查仍独立保留，不能互相抵消。"),
    ('有没有浮空支撑岛？',f"{m['floating_support_islands']} 个。内部空腔的负向壳没有被误当作浮空实体。"),
    ('有没有宏观支撑失效？',f"{m['macro_fail_count']} 项；另有 {m['micro_deficit_count']} 项附着的微小逐层亏缺。判据结合实际连通、上一层承托深度和连续传播，不以所有多边形差为零作为要求。"),
    ('Bambu 是否还有未解决关键警告？',f"{a['warnings']['critical_unresolved']} 条。对象级 floating 警告结合实际 G-code、真实缺点与宏观检查定位；有真实缺点时没有直接改名为非关键警告。"),
    ('原先巨大支撑有几个组件？',f"{removal['initial_components']} 个实体连通组件。"),
    ('虚拟切割后有几个主要组件？',f"{removal['major_component_count']} 个。这里只证明材料图确实分开，不能等同于全部工具可达或安全取出。"),
    ('最大组件接触弧长？',f"{removal['maximum_component_contact_arc_mm']:.6f} mm；由各分支实际邻近接触采样的弧长跨度求和。"),
    ('最大有效力臂？',f"建议几何夹持点到最近接触点为 {removal['max_effective_lever_arm_mm']:.6f} mm；同时到最远接触点可达 {removal['max_farthest_contact_distance_mm']:.6f} mm。未通过工具接近检查的夹持点只是待验证点，不能用较小的最近距离证明整块长梁已安全；没有输出力或应力。"),
    ('最后接触小块最长多少？',f"诊断性顶面残片最长尺寸 {removal['max_last_contact_piece_length_mm']:.6f} mm，以实际网格有向包围盒的最长边计算；接触弧长另列。bulk/roof 分离没有自动视作已执行剪切，仍连接成长梁的部分明确不合格。"),
    ('有没有闭合套笼？',f"实际切片截面检查发现 {removal['closed_cage_count']} 个。最大局部包覆角 {removal['maximum_wrap_angle_deg']:.3f}°；120° 不再是拒绝门槛。"),
    ('有没有需要弯曲血管才能拆的组件？','流程禁止弯曲血管。未找到刚性安全路径的组件保持阻塞，没有假设血管变形后即可通过。'),
    ('有没有分叉需要横向强拉？','没有授权这种动作。分叉附近采用更短接触残片的相对目标，取出方向必须朝外并避开轴向拖拽；未通过时保留为失败。'),
    ('预计多少次主要剪切？',f"规划 {removal['primary_cut_count']} 次战略剪切，顺序可达性实际通过 {removal['actual_executed_virtual_cut_count']} 次。其余未被计入已完成；后续残片剪切次数尚不能可靠确定。"),
    ('是否所有剪切必须一开始可达？','不需要。每次移除后更新实际障碍物；没有把未来尚未剪掉的材料提前当作空隙。'),
    ('是否给出完整顺序？',f"完整性：{seq['complete']}。已输出 {len(seq['steps'])} 个几何检查通过的局部动作；仍有 {len(seq['remaining_components'])} 个组件、{len(seq['pending_cuts'])} 个剪口阻塞。参见 SAFE_SUPPORT_REMOVAL_GUIDE.md，失败时不是生产操作许可。"),
    ('现在是否值得进入第一次完整真实打印试验？','可以进入 REAL PRINT TRIAL，不能称为制造验证。' if s['production_frozen'] else '当前不应把这一整件当作通过的软件生产候选。剩余真实缺点与未闭合拆卸计划仍阻止验收；本轮已在规定次数内停止，没有用改变验收口径掩盖它们。')]
    lines=['# BG001 RMCA 血管芯支撑最终实用验收报告','',f"状态：**{s['status']}**。当前对象：`{s['selected_candidate']}`。",'', '## 验收问题与结果','']
    for i,(q,value) in enumerate(answers,1):lines += [f'{i}. **{q}** {value}','']
    lines += ['## 真实运行与解释','',f"本轮完成 {s['repair_round_count']} 轮真实局部修补、{s['topology_round_count']} 轮拓扑切片；没有新增姿态搜索。已有三个候选由当前数据重新分类，原 O002 基线有 5 个真实缺点与 17 个过近点。第二备选虽缺点较少，但存在宏观支撑失败，未在全部要求上严格优于基线。详情见 reports/baseline_selection.json。",'',f"保留历史严格指标用于对照：中心线/端口合并 {s['legacy_metrics']['coverage']['centerline_underbody_percent']:.6f}%，逐层区域 {s['legacy_metrics']['layer']['percent']:.6f}%。它们没有作为新 READY 的数学 100% 门槛。",'',f"真实 Bambu 版本 {s['bambu']['version']}，X1 Carbon / 0.2 mm 喷嘴 / ABS / 0.10 mm Standard Vessel；Auto Support OFF：{s['bambu']['auto_support_off']}。每个 3MF 含独立血管和显式支撑两部分。C3 在当前候选中的截面与刀路联合验证：{removal['thin_web_preserved']}；保留不等于实物断开力已经验证。",'','## 未解决位置','']
    for r in missing:lines += [f"- `{r['sample_id']}`，{r['new_class']}，坐标 {np.round(r['xyz'],6).tolist()} mm，可用高度范围内最近实际挤出横向距离 {r['nearest_support_path']['lateral_distance_mm']:.6f} mm。较低层见证与最近刀路分别记录，不能用较低幕墙冒充顶面。"]
    lines += ['','## 每轮实测比较','','|阶段|轮次|真实缺承托|宏观失效|C3 联合验证|主要组件|最大接触弧长 (mm)|顺序完整|','|---|---|---|---|---|---|---|---|']
    for r in s['round_history']:
        a=r['practical'];lines.append(f"|{r['stage']}|{r['round']}|{a['structural']['true_missing_support_count']}|{a['macro']['macro_fail_count']}|{r.get('c3_preserved','—')}|{r.get('major_components','—')}|{r.get('max_contact_arc_mm','—')}|{r.get('sequence_complete','—')}|")
    lines += ['',f"实际原生切片调用 {s.get('actual_native_slice_attempt_count','见清单')} 次，清单见 validation/native_run_inventory.json。C3 联合验证同时要求实际刀路见证存在、薄颈截面没有被相邻原幕墙增厚；二者分别记录在候选的 practical_removal/removal_audit.json。",'']
    lines += ['',f"刀路见证点保留：{removal['actual_extrusion_thin_link_witnesses_preserved']}；独立 C3 薄颈截面与刀路的联合验证：{removal['thin_web_preserved']}。顶面残片由原构造顶面及实际局部新增材料提取，排除幕墙平移差产生的边缘薄壳。",'','## 输出与复核','',f"所有新结果集中在 `{out}`。旧 closure_v2 及其测试基准保留。",'','- `selected/frozen_vascular_core.stl`：与上游冻结源逐字节相同。','- `selected/vascular_core_print_frame.stl`：仅刚体变换后的打印坐标副本，与支撑对齐。','- `selected/custom_support.stl`：当前支撑。','- `selected/final_project.3mf`：可导入 Bambu 的两组件工程。','- `selected/final_sliced.3mf` 与 `selected/plate_1.gcode`：真实切片档案与刀路。','- `tables/`：重新分类、宏观审计、剪口、主要组件、最终接触残片。','- `QC_FINAL/`：18 张检查图；红色缺承托与黄色过近严格区分，阻塞的拆卸步骤明确标注。','- `production_candidate_freeze.json`：源/支撑/3MF/profile/代码 SHA、版本、变换、时间及真实冻结状态。失败时 production_frozen=false。','',f"源保护检查：{s['source_protection']['checked_files']} 个文件，改变 {len(s['source_protection']['changed'])} 个。",'','```bash','cd /home/lzy/projects/vascular_printing','./run_python.sh s1-11_vascular_core_support.py --stage review','```','','复核复用完成的真实切片；NEEDS_ADJUSTMENT 返回码为 2，不代表程序崩溃。默认入口使用本轮有界流程，不再执行旧 30+ 次循环。已完成的最多三轮不会因重新运行而从零开始。','', '## 实物不确定性','','REAL ABS SUPPORT REMOVAL FORCE NOT VALIDATED。REAL VESSEL DAMAGE RISK NOT VALIDATED。SURFACE QUALITY NOT VALIDATED。所有力臂、包覆与扫掠都是几何代理；未给出未经试验的 N、MPa 或断裂力。','']
    tests=out/'validation/test_results.json'
    if tests.exists():lines += ['## 回归测试','',json.dumps(json.loads(tests.read_text()),ensure_ascii=False),'']
    (out/'vascular_core_support_final_practical_report.md').write_text('\n'.join(lines))
    guide=['# 支撑顺序拆除检查指南','',f"状态：**{s['status']}**。",'','当前未完整通过时，以下表格只用于核对已计算动作，不能按顺序强行完成剩余拆卸。整棵血管不得作为把手，禁止扭动大组件、沿血管轴线拖拽或横拉分叉。先剪远离血管的连接，再取下 bulk，最后只处理短接触残片；未验证的 bulk/roof 分离不会被写成成功步骤。','','坐标为打印坐标系，单位 mm。接近方向指从目标向工具所在空间的向量；取出方向为组件移动向量。检查图 13–15 显示前两步及最后可完成状态。','','|步骤|动作及真实 ID|夹持/剪切位置 (mm)|接近或取出方向|附近血管|禁止动作|','|---|---|---|---|---|---|']
    for step in seq['steps']:
        guide.append(f"|{step['step']}|{step['action']} {step.get('cut_id',step.get('component_id'))}|{np.round(step.get('point',step.get('grip_point')),3).tolist()}|{np.round(step.get('access_direction',step.get('direction')),3).tolist()}|{step.get('nearby_branch',step.get('nearby_branches'))}|向血管内推、轴向拉、整体扭动、横拉分叉|")
    guide += ['','## 必须停止的位置','']
    for r in seq['blocked']:guide += [f"- **{r['component_id']}**：{r['reason']}；尚未可达的相关剪口：{r['pending_cuts']}。禁止把此组件整体当作杠杆拉下。"]
    guide += ['','尚未完成的剪口：'+', '.join(seq['pending_cuts'])+'。','', '本指南没有假设 ABS 的实际剪断力、剥离力或细管承载能力已验证。','']
    (out/'SAFE_SUPPORT_REMOVAL_GUIDE.md').write_text('\n'.join(guide))
