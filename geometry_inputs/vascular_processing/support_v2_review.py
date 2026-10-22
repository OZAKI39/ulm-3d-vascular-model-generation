"""Evidence-led V2 review: same-scale views, real toolpaths and explicit failures."""
from pathlib import Path
import json,shutil,collections
import numpy as np
import trimesh
from .support_inputs import save,table
from .support_aware_review import annotate
from .sacrificial_fixture_review import polydata
from .support_v2_primitives import primitive_gallery

COLORS={'ROOF_PAD':'#d44b52','VERTICAL_PILLAR':'#d59a32','SELF_SUPPORTING_BRANCH':'#428cab','Y_MERGE':'#3266a1','LONGITUDINAL_RAIL':'#9361a6','SHORT_BRIDGE':'#6e4691','BREAKAWAY_FOOT':'#665d50','BREAKAWAY_SEGMENT_NECK':'#e968b0','BUTTRESS':'#3eac9b'}
STEMS=['v1_bad_cone_support','orientation_search_overview','top_orientation_candidates','fdm_unsupported_map','freefloating_evidence','support_required_intervals','dense_roof_candidates','selected_roof_pads','fdm_primitive_library','merge_graph','shared_trunks','rail_supports','bridge_scaffolds','root_domains','breakaway_segmentation','support_only_final','support_plus_vascular','removal_sequence','v1_vs_v2','bambu_auto_support_off_preview']

def read(path,default=None):return json.loads(path.read_text()) if path.exists() else default

def comparison(out,results):
    v1=read(out.parent/'legacy_explicit_support_v1/vascular_support_summary.json',{})
    old=v1.get('selected',{});s=old.get('slice',{})
    old_graph=read(out.parent/'legacy_explicit_support_v1/custom_candidates'/old.get('candidate_id','')/'support_graph.json',{})
    rows=[dict(method='legacy_explicit_support_v1',candidate=old.get('candidate_id'),orientation='CASTING_UPRIGHT',critical_coverage_percent=old.get('planned_coverage',{}).get('percent'),coverage_definition='V1_5MM_ARC_SAMPLES_NOT_COMPARABLE_TO_V2_POLYGON',roof_count=old.get('planned_tip_count'),merge_nodes=len(old_graph['merge_nodes']) if 'merge_nodes' in old_graph else None,root_count=len(old_graph['roots']) if 'roots' in old_graph else None,components=old.get('removal',{}).get('component_count'),rail_count=0,bridge_count=0,support_volume_mm3=old.get('support_volume_mm3'),trapped_pieces=old.get('removal',{}).get('trapped_support_count'),bambu_warnings=s.get('critical_warnings'),near_vessel_actual_area_mm2=s.get('contact_area_mm2'))]
    for r in results:
        sl=r.get('slice',{});cov=sl.get('actual_vascular_coverage',r.get('vascular_coverage',{}))
        rows.append(dict(method='FDM_NATIVE_VASCULAR_SUPPORT_V2',candidate=r['candidate_id'],orientation=r['orientation']['gravity_in_casting'],critical_coverage_percent=cov.get('percent'),coverage_definition='V2_COMPLETE_LAYER_POLYGONS_ACTUAL_ROOFS' if sl.get('slicer_success') else 'V2_COMPLETE_LAYER_POLYGONS_PLANNED_ROOFS',roof_count=r.get('selected_roof_count'),merge_nodes=r.get('merge_count'),root_count=r.get('root_count'),components=r.get('support_component_count'),rail_count=r.get('rail_count'),bridge_count=r.get('bridge_count'),support_volume_mm3=r.get('support_volume_mm3'),trapped_pieces=r.get('removal',{}).get('trapped_pieces'),bambu_warnings=sl.get('warnings'),near_vessel_actual_area_mm2=sl.get('contact_area_mm2'),failure=r.get('hard_gate_failures',r.get('error'))))
    table(out/'V1_vs_V2_topology.csv',rows);table(out/'v1_vs_v2.csv',rows)
    return rows

def report(out,summary,results):
    best=summary.get('selected');comparison(out,results)
    if not best:
        (out/'vascular_support_v2_report.md').write_text('# FDM 血管支撑 V2\n\n状态：NEEDS_ADJUSTMENT。没有完成可导出的候选；全部失败原因见 planning_results.json。\n');return
    s=best.get('slice',{});r=best.get('removal',{});cov=s.get('actual_vascular_coverage',best.get('vascular_coverage',{}));k=best.get('primitives_by_kind',{});m=best.get('mechanical_proxy',{})
    direction=np.array(best['orientation']['gravity_in_casting']);changed=not np.allclose(direction,[0,0,1]);g=best.get('hard_gate_failures',[])
    area=s.get('contact_area_mm2');old=read(out.parent/'legacy_explicit_support_v1/vascular_support_summary.json',{});native=old.get('native_benchmark',{})
    # The same actual downward-face metric used by V1 is retained. Do not
    # substitute nominal pad area or silently equate different coverage gates.
    native_area=old.get('baseline_contact_area_mm2')
    if native_area is None:
        selected=old.get('selected',{});v1area=selected.get('slice',{}).get('contact_area_mm2');decline=selected.get('contact_area_reduction_percent')
        if v1area is not None and decline is not None and decline<100:native_area=v1area/(1-decline/100)
    answers=[
      f"最终打印姿态改变了吗？当前诊断候选为 {best['candidate_id']}，{'改变了打印重力方向' if changed else '仍为浇注开口朝上'}；浇注坐标中的打印 +Z 为 {np.round(direction,5).tolist()}。没有通过全部门槛时，正式生产候选保持空值。",
      f"为什么选择此姿态？它是 128 个粗筛姿态、12 个完整规划候选中，先检查硬门槛，再按危险区域覆盖和近血管负担排序的诊断候选。未通过者不会因只生成很少支撑、少失败一个条件而获得不合理优势。粗筛盒体悬空负担为 {best['orientation']['box_overhang_burden_mm2']:.2f} mm²；尚未证明它已经比所有对照更容易实物打印。",
      f"危险血管区间多少？重新检测得到 {best.get('critical_intervals')} 个融合区间；主判据包含 {cov.get('region_count')} 个需支撑层区域，端口采用相同标准。旧版 19 个区间没有直接复用。",
      f"是否全部覆盖？{'是' if cov.get('all_required_covered') else '否'}。{'真实刀路' if s.get('slicer_success') else '设计顶垫'}覆盖面积为 {cov.get('percent',0):.2f}%，仍有 {cov.get('uncovered_region_count')} 个层区域未充分承托。百分比是各层面积之和，非三维表面积。",
      f"用了多少顶垫？从 {best.get('raw_roof_count')} 个密集候选中，实际连接 {best.get('selected_roof_count')} 个短顶垫。删除冗余项以完整多边形覆盖不变为条件，没有固定数量上限。",
      f"有多少共享合并？{best.get('merge_count')} 个 Y 合并；逐项修复和拒绝原因见 merge_attempts.csv。",
      f"有多少根部？{best.get('root_count')} 个实际生成的根部，包含构建板及合法盒面根；非构建板根具有可折断窄颈。",
      f"用了多少柱、轨道和桥？柱段 {k.get('VERTICAL_PILLAR',0)} 个、短轨道 {best.get('rail_count',0)} 段；其中承担跨越功能的桥跨 {best.get('bridge_count',0)} 个。桥跨数与轨道数存在重叠，不相加。",
      f"有没有退化回锥柱林？没有使用长锥形默认构件。退化检测结果：{best.get('topology_warning') or '未触发'}；这不等于支撑负担已经足够低。",
      f"支撑自身还需要支撑吗？{'逐层筛查未发现额外需求' if best.get('self_buildability',{}).get('passed') else '仍发现欠支撑需求'}，未通过区域 {best.get('self_buildability',{}).get('unsupported_regions')} 个。",
      f"有没有被困碎片？几何路径检查发现 {r.get('trapped_pieces','尚未计算')} 个；{r.get('fragment_count','尚未计算')} 个待取出碎片按浇注开口检查，并把尚未取出的碎片作为障碍。",
      f"关闭自动支撑后 Bambu 真实切片成功了吗？{s.get('slicer_success',False)}；归档配置确认 Auto Support OFF：{s.get('auto_support_off',False)}；独立支撑及合法盒基础切片：{s.get('support_only_slicer_success',False)}。模型和细构件的保留另见刀路检查，不以切片成功代替。",
      f"相比 V1 解决了什么？从固定姿态和整图退回独立柱，改为逐姿态检测、量化矩形顶垫/柱/合并/轨道、盒外落根和分段拆除检查；本候选有 {best.get('merge_count')} 个共享合并。仍未通过的条件全部保留，不能宣称已解决覆盖与拆除。",
      f"相比 Bambu Organic 近血管负担如何？本候选实际刀路的近血管面面积为 {f'{area:.3f} mm²' if area is not None else '尚无结果'}，原生基准为 {f'{native_area:.3f} mm²' if native_area is not None else '见保留的原生基准表'}。由于覆盖率和姿态不同，面积下降不能证明等覆盖条件下更优。",
      f"下一步最值得打印什么？{'先人工检查 '+best['candidate_id']+'，再打印 V2 校准试片' if not g else '先打印并测量 V2 校准试片；当前完整候选不推荐作为正式打印件'}。0.20 mm 间隙及 2 mm 桥跨仍分别标记 PHYSICAL_GAP_UNVALIDATED、BRIDGE_LIMIT_UNCALIBRATED。"
    ]
    lines=['# BG001 RMCA：FDM 姿态与支撑联合规划 V2','',f"软件状态：**{summary['status']}**。这不是制造、材料强度或实物拆除验证。",'']+[f'{i+1}. {a}' for i,a in enumerate(answers)]
    lines+=['','## 方法与判定依据','',
        '首先以 Fibonacci 重力方向与确定性滚转生成 128 个刚体姿态，保留 12 个不同重力方向进行完整检测。血管和 2 mm 盒体仅旋转平移，比例严格为 1。逐姿态重新运行已有四类检测器，并以逐层欠支撑多边形为主判据。5 mm 弧长仅保留为旧版比较量。',
        '短顶垫沿局部血管方向生成，宽厚由真实线宽和层高量化。短桥必须有两端基础；软件桥接上限来自匹配工艺的实物校准文件，否则保守采用 2 mm，并保存 2/4/6/8 mm 敏感性结果。顶垫间隙保存 0.10/0.15/0.20/0.25 mm 候选，未以软件指标宣称某一间隙最佳。',
        '生成的支撑单独作逐层分析，只允许合法盒面作为基础；血管不能用于掩盖支撑自身的悬空。STL 单精度接缝的数值清理仅针对新支撑，最大容差见图文件，未平滑或缩放原血管。',
        '机械细长比按整根连续主干的无侧向支撑高度与最窄断颈计算；增加断颈不等于增加侧向约束，不能把每一小段分别计算后宣称主干更稳定。该指标仅用于保守排序，不代替 ABS 强度试验。',
        '分段移除先模拟切断设计窄颈，再逆变换回浇注姿态。内部碎片必须经过真实盒顶开口；尝试平移、倾斜和组合路径，逐段扫掠碰撞，已经移出的碎片才从障碍中删除。虚拟断口厚度为 0.05/0.06 mm，真实断口间隙、手工具空间和折断力未验证。另检查相邻未削弱的支撑是否填充了设计断面，不能把虚拟切割当作可以轻易折断的证明。有限路径中失败被保守标记为被困，不能据此证明所有可能路径均不存在。',
        '', '## 本机工艺与切片证据','',f"工艺：{summary['profile']['printer']} / {summary['profile']['process_name']}；喷嘴 {summary['profile']['nozzle_diameter']} mm，线宽 {summary['profile']['line_width']} mm，层高 {summary['profile']['layer_height']} mm。速度、温度、桥参数及完整工艺见 profile_runtime.json 与 profiles/。",'']
    for r0 in results:
        if r0['candidate_id'] in summary.get('top4',[]):
            a=r0.get('slice',{});lines.append(f"- {r0['candidate_id']}：Bambu 成功={a.get('slicer_success',False)}；Auto Support OFF={a.get('auto_support_off',False)}；全部构件保留={a.get('all_features_preserved',False)}；未解释警告对应区域={a.get('unresolved_warning_count','待查')}；被困碎片={r0.get('removal',{}).get('trapped_pieces','待查')}。")
            failed_exports=[row['file'] for row in r0.get('semantic_part_export_audit',[]) if not row['closed_positive_volume']]
            if failed_exports:lines.append(f"  辅助拆分文件 {', '.join(failed_exports)} 未通过 STL 精度检查，只作诊断；完整支撑和整件 STL 的核验单独记录。碎片 STL 精度失败数为 {r0.get('removal',{}).get('precision_failed_fragments',0)}，不能作为独立打印件使用。")
    lines+=['','Bambu 原生悬空警告没有逐项坐标时，floating_warning_localization.csv 是独立几何层与真实 G-code 的对应证据，不冒充原生诊断位置。只有实际顶垫刀路覆盖完整的区域才允许解释为有意间隙；其余保留为未解释风险。QC/20 是实际挤出轨迹重建，不是软件界面截图。', '', '## 当前失败与边界','']
    meanings={'UNSUPPORTED_VASCULAR_OR_PORT_REGIONS':'设计顶垫尚未完整承托血管或端口','SUPPORT_REQUIRES_SUPPORT':'支撑自身仍有欠支撑区域','TRAPPED_OR_UNASSESSED_REMOVAL_SEGMENTS':'仍有未找到安全拆除路径的碎片','BREAK_SECTION_CROSSES_UNWEAKENED_SUPPORT':'预定断面穿过了相邻未削弱的支撑，不能假定可以轻易折断','ACTUAL_ROOF_COVERAGE_INCOMPLETE':'实际切片刀路未完整覆盖需要承托的区域','UNRESOLVED_FLOATING_WARNING':'切片悬空警告仍有未被实际承托证据解释的区域'}
    lines += ['- '+meanings.get(x,x)+f'（`{x}`）' for x in g] or ['全部软件硬门槛通过，仍需实物试片。']
    fracture=r.get('break_design_audit',{});exit_check=r.get('opening_exit_audit',{})
    lines+=['',f"诊断候选的断面检查：{'通过' if fracture.get('passed') else '未通过'}。虚拟切割在设计窄颈之外还穿过了 {fracture.get('cut_volume_outside_designed_necks_mm3',0):.3f} mm³ 支撑材料，不能据此宣称真实折断可行。已找到路径的 {exit_check.get('checked_fragments',0)} 个碎片，完整离开浇注开口的终点检查{'通过' if exit_check.get('passed') else '未通过或未完成'}。发生 STL 接缝精度失败的碎片为 {r.get('precision_failed_fragments',0)} 个；这类碎片仅用带索引的原始布尔网格继续碰撞诊断，不算 STL 验证通过。"]
    for candidate in results:
        if candidate['candidate_id'] in summary['top4'] and candidate.get('contact_semantic_backfill'):
            lines+=['',f"{candidate['candidate_id']} 的早期缓存存在分叉距离缺值，已用实际 SWC 分叉坐标补齐，记录见 contact_semantic_backfill.json。该补全没有改变几何或重新选择顶垫；缺值在原轮次未参与分叉惩罚，不追认原选择为使用过完整惩罚。"]
    coupon=read(out/'coupon/bambu_validation.json',{})
    lines+=['',f"V2 试片：Bambu 实际切片成功={coupon.get('slicer_success',False)}，自动支撑关闭={coupon.get('auto_support_off',False)}，{len(coupon.get('feature_preservation',{}).get('rows',[]))} 个检查特征全部保留={coupon.get('feature_preservation',{}).get('passed',False)}，短跨实际挤出证据={coupon.get('bridge_toolpaths',{}).get('passed',False)}。部分短跨由外围线完成，并非原生 Bridge 标签；两端均与上一层真实刀路对应才予确认。以上不是 ABS 物理桥接校准。"]
    lines+=['',f"机械代理：最长柱 {m.get('max_pillar_height_mm',0):.2f} mm，最大细长比 {m.get('max_slenderness',0):.2f}；它不构成 ABS 强度结论。",'', '## 输出与复现','',
        '入口仍为 s1-8，没有创建 s1-9。默认和 --mode custom 均运行 V2；--mode benchmark 保留原生对照。',
        '```bash','cd /home/lzy/projects/vascular_printing','./run_python.sh s1-8_vascular_support_planner.py --stage all','./run_python.sh s1-8_vascular_support_planner.py --stage report','./run_python.sh -m pytest tests/test_s1_8_vascular_support_planner.py -q','```','',
        'candidates/<姿态编号>/ 保存六类 STL、图、层覆盖、合并日志、真实 Bambu/Orca 归档及拆除顺序；top4_candidate_ids.json 指向保留的四个候选。print_transform_v2.json 与 casting_restore_transform_v2.json 对应诊断候选，生产选项以摘要中 production_selected 为准。coupon/ 保存试片及实际切片。',
        f"上游保护：检查了 {summary.get('source_protection',{}).get('checked_files',0)} 个文件，SHA256 {'全部不变' if summary.get('source_protection',{}).get('all_unchanged') else '存在变化，详见 source_protection.json'}。未发送、上传或启动任何打印任务。V1 快照及旧测试、原生基准完整保留在父目录 legacy_explicit_support_v1 和 benchmark；V1 禁止进入 production。"]
    tests=out/'test_results.xml'
    if tests.exists():
        import xml.etree.ElementTree as ET
        suites=list(ET.parse(tests).getroot().iter('testsuite'));lines+=['',f"永久测试：{sum(int(x.get('tests',0)) for x in suites)} 项，失败/错误 {sum(int(x.get('failures',0))+int(x.get('errors',0)) for x in suites)} 项。"]
    (out/'vascular_support_v2_report.md').write_text('\n'.join(lines)+'\n')

def generate(out,summary,results):
    report(out,summary,results)
    best=summary.get('selected')
    if not best:return
    folder=Path(best['folder']);qc=out/'QC';qc.mkdir(exist_ok=True)
    for name in ['critical_intervals.csv','roof_candidates.csv','selected_roofs.csv','support_graph_v2.json','support_primitives.csv','support_roots.csv','support_segments.csv','removal_sequence.json','support_removal_plan.json','support_removal_plan.md','bambu_validation.json']:
        if (folder/name).exists():shutil.copy2(folder/name,out/name)
    import pyvista as pv
    from PIL import Image
    mold=trimesh.load_mesh(folder/'mold_without_support.stl');vascular=trimesh.load_mesh(folder/'vascular_with_ports.stl');box=trimesh.load_mesh(folder/'box.stl');support=trimesh.load_mesh(folder/'custom_support.stl')
    graph=read(folder/'support_graph_v2.json');center=mold.bounds.mean(0);extent=float(mold.extents.max());camera_scale=extent*.62
    def shot(name,title,items,notes=(),points=None,lines=(),view='iso',local_center=None,scale=None):
        p=pv.Plotter(off_screen=True,window_size=[1100,850]);p.set_background('white')
        for item,color,opacity in items:
            if item is None:continue
            if isinstance(item,trimesh.Trimesh):
                if not len(item.faces):continue
                item=polydata(item)
            p.add_mesh(item,color=color,opacity=opacity,smooth_shading=False)
        if points is not None and len(points):p.add_mesh(pv.PolyData(np.array(points)),color='#cf3939',point_size=6,render_points_as_spheres=True)
        for a,b,color in lines:p.add_mesh(pv.Line(a,b),color=color,line_width=2)
        direction=np.array({'iso':[1,-1.25,1.1],'front':[0,-1,0],'side':[1,0,0],'top':[0,0,1]}[view],float);direction/=np.linalg.norm(direction);c=center if local_center is None else np.array(local_center);span=camera_scale if scale is None else scale
        p.camera_position=[c+direction*extent*3,c,[0,1,0] if view=='top' else [0,0,1]];p.enable_parallel_projection();p.camera.parallel_scale=span
        image=Image.fromarray(p.screenshot(return_img=True));p.close();annotate(image,title,notes).save(qc/name)
    ghost=[(vascular,'#4f9b72',.22),(box,'#8e9cad',.07)];base=[(vascular,'#4f9b72',1),(box,'#8e9cad',.08)]
    # Read only our own local planning cache, for exact primitive grouping.
    import pickle
    with (folder/'scaffold_local_cache.pkl').open('rb') as f:_,_,pieces,_=pickle.load(f)
    def colored(kinds=None):
        result=[]
        for kind,color in COLORS.items():
            meshes=[m for m,r in zip(pieces,graph['primitives']) if r['kind']==kind and (kinds is None or kind in kinds)]
            if meshes:result.append((trimesh.util.concatenate(meshes),color,1))
        return result
    legacy=read(out.parent/'legacy_explicit_support_v1/vascular_support_summary.json',{});old=legacy.get('selected',{});oldpath=Path(old.get('folder',''))/'custom_support.stl';v1=None
    if oldpath.exists():
        from .support_inputs import discover
        inputs=discover();v1=trimesh.load_mesh(oldpath);v1.apply_transform(np.array(best['orientation']['print_transform'])@np.linalg.inv(inputs['print_transform']))
    shot('01_'+STEMS[0]+'.png','V1：独立锥柱退化，仅作为历史对照',ghost+([(v1,'#c78d43',1)] if v1 is not None else []),notes=['同尺度、同观察方向；不是新生产候选'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    orientations=read(out/'orientation_candidates.json',[]);fig=plt.figure(figsize=(11,8));ax=fig.add_subplot(projection='3d');xyz=np.array([r['gravity_in_casting'] for r in orientations]);c=ax.scatter(*xyz.T,c=[r['coarse_rank_value'] for r in orientations],cmap='viridis');ax.set(xlabel='Casting X',ylabel='Casting Y',zlabel='Casting Z',title='128 rigid orientations / 32 gravity directions');fig.colorbar(c,ax=ax,label='Coarse FDM burden (surrogate)');fig.savefig(qc/'02_orientation_search_overview.png',dpi=140);plt.close(fig)
    thumbs=[]
    for r in results:
        path=Path(r['folder'])/'mold_without_support.stl'
        if not path.exists():continue
        mesh=trimesh.load_mesh(path);name=f"orientation_{r['candidate_id']}.png";shot(name,r['candidate_id'],[(mesh,'#65a282',1)],notes=[f"层覆盖 {r.get('vascular_coverage',{}).get('percent',0):.2f}%"],local_center=mesh.bounds.mean(0))
        thumb=Image.open(qc/name);thumb.thumbnail((330,255));cell=Image.new('RGB',(330,255),'white');cell.paste(thumb,((330-thumb.width)//2,(255-thumb.height)//2));thumbs.append(cell)
    canvas=Image.new('RGB',(1320,255*max(1,int(np.ceil(len(thumbs)/4)))),'white')
    for i,im in enumerate(thumbs):canvas.paste(im,((i%4)*330,(i//4)*255))
    canvas.save(qc/'03_top_orientation_candidates.png')
    rows=read(folder/'vascular_layer_regions.json',[]);shot('04_fdm_unsupported_map.png','逐层 FDM 欠支撑区域',ghost,points=[[*r['representative_xy'],r['z_mm']] for r in rows if r['classification']=='SUPPORT_REQUIRED'],notes=['点显示区域位置；完整多边形参与数值覆盖，不以点替代面积'])
    det=read(folder/'detection.json',{});raw=det.get('raw',[]);npz=folder/'detector_points.npz';ff=np.load(npz)['freefloating'] if npz.exists() else [p['xyz'] for p in raw if 'FREEFLOATING' in p.get('detectors_triggered',[])]
    shot('05_freefloating_evidence.png','FreeFloating：必须承托的证据',ghost,points=ff)
    shot('06_support_required_intervals.png','重新检测的危险血管区间',ghost,points=[p['xyz'] for p in raw if p['risk_level']=='CRITICAL'])
    dense=read(folder/'roof_candidates_pre_graph.json',[]);shot('07_dense_roof_candidates.png',f'密集顶垫候选：{len(dense)}',ghost,points=[r['xyz'] for r in dense])
    shot('08_selected_roof_pads.png','最终已连接短顶垫',ghost+colored({'ROOF_PAD'}))
    gallery=primitive_gallery(summary['profile']);items=[]
    for i,(kind,mesh) in enumerate(gallery.items()):mesh.apply_translation([(i%3)*16,(i//3)*16,0]);items.append((mesh,COLORS[kind],1))
    shot('09_fdm_primitive_library.png','FDM 构件库：量化截面与短顶垫',items,notes=list(COLORS),local_center=[16,16,4],scale=32)
    shot('10_merge_graph.png',f"共享合并图：{graph['merge_count']} 个 Y 节点",ghost,lines=[(e['start'],e['end'],'#3e84ae') for e in graph['edges']],points=[n['xyz'] for n in graph['nodes']])
    shot('11_shared_trunks.png','共享主干与载荷加粗',ghost+colored({'VERTICAL_PILLAR','SELF_SUPPORTING_BRANCH'}))
    shot('12_rail_supports.png',f"分段短轨道：{graph['rail_count']} 段",ghost+colored({'LONGITUDINAL_RAIL'}),notes=['零段时保留空结果，不用示例冒充实际结构'])
    shot('13_bridge_scaffolds.png',f"两端承托的桥跨：{graph['bridge_count']} 个",ghost+colored({'LONGITUDINAL_RAIL','SHORT_BRIDGE','VERTICAL_PILLAR'}),notes=['2 mm 保守软件上限；物理桥接能力未校准'])
    shot('14_root_domains.png','构建板与合法盒面的根部',[(box,'#8e9cad',.13)]+colored({'BREAKAWAY_FOOT'}),points=[r['xyz'] for r in graph['roots']],notes=[str(dict(collections.Counter(r['domain'] for r in graph['roots'])))])
    shot('15_breakaway_segmentation.png','设计断颈与拆分位置',colored(),points=[n['xyz'] for n in graph['segment_necks']])
    legend=['红：顶垫；金：立柱；蓝：合并臂；紫：短轨道/桥','粉：分段颈；灰：底脚；青：局部渐扩托台']
    shot('16_support_only_final.png','支撑实体：各构件按功能着色',colored(),notes=legend)
    for view in ['front','side','top']:
        shot(f'16_support_only_{view}.png',f'支撑正交视图：{view}',colored(),notes=legend,view=view)
    clipped=polydata(support).clip(normal=[0,1,0],origin=center)
    shot('16_support_section.png','支撑剖切：保留半侧几何',[(clipped,'#d49a3c',1)],view='front',notes=['剖切只用于检查；未修改导出 STL'])
    shot('17_support_plus_vascular.png','冻结血管与 V2 支撑',base+colored(),notes=[summary['status']])
    removal=best.get('removal',{});fragitems=[];mat=np.array(best['orientation']['print_transform'])
    for row in removal.get('segments',[]):
        fp=folder/'removal_fragments'/f"fragment_{int(row['segment_id'][1:]):04}.stl"
        if fp.exists():
            piece=trimesh.load_mesh(fp);piece.apply_transform(mat);color='#a64fbc' if row.get('fragment_stl_precision_failure') else '#d84242' if row['status']=='TRAPPED_REMOVAL_FRAGMENT' else '#777777' if row['status']=='APPROVED_FOOT_STUB_REMAINS_ON_BOX' else '#63a779';fragitems.append((piece,color,1))
    shot('18_removal_sequence.png','分段移除：红色表示未找到安全路径',ghost+fragitems,notes=[f"被困 {removal.get('trapped_pieces','待查')}；顺序见 support_removal_plan.md",'绿色：找到路径；灰色：允许保留的脚垫；紫色：碎片 STL 精度未通过','路径在浇注坐标检查；此图仅恢复到共同打印视角'])
    left=Image.open(qc/'01_v1_bad_cone_support.png');right=Image.open(qc/'17_support_plus_vascular.png');canvas=Image.new('RGB',(left.width+right.width,max(left.height,right.height)),'white');canvas.paste(left,(0,0));canvas.paste(right,(left.width,0));canvas.save(qc/'19_v1_vs_v2.png')
    actual=folder/'bambu/actual_custom_paths.npz';extra=[]
    if actual.exists():
        points=np.load(actual)['points'];extra=[(pv.PolyData(points[::max(1,len(points)//300000)]),'#c78d36',1)]
    shot('20_bambu_auto_support_off_preview.png','Bambu 真实自定义支撑挤出轨迹',base+extra,notes=[f"Auto Support OFF：{best.get('slice',{}).get('auto_support_off',False)}；实际切片成功：{best.get('slice',{}).get('slicer_success',False)}",'真实 G-code 轨迹重建；不是 Bambu 原生界面截图'])
    save(out/'visual_qc_manifest.json',dict(candidate=best['candidate_id'],main_images=[str(qc/f'{i+1:02}_{stem}.png') for i,stem in enumerate(STEMS)],common_camera_center=center.tolist(),common_parallel_scale_mm=camera_scale,geometry_views_same_scale=True,source_geometry_unchanged=True,toolpath_preview_is_native_gui_screenshot=False))
