"""Same-scale explicit-support QC and a candid Chinese manufacturing report."""
from pathlib import Path
import json,shutil
import numpy as np
import trimesh
from .support_aware_review import annotate
from .sacrificial_fixture_review import polydata
from .support_review import markers,path_mesh
from .support_inputs import save

def generate(out,mesh,box,det,contacts,graph,best,summary):
    import pyvista as pv
    from PIL import Image
    qc=out/'QC';qc.mkdir(exist_ok=True);folder=Path(best['folder']);support=trimesh.load_mesh(folder/'custom_support.stl')
    lo=np.minimum(mesh.bounds[0],box.bounds[0]);hi=np.maximum(mesh.bounds[1],box.bounds[1]);center=(lo+hi)/2;extent=float(max(hi-lo))
    def shot(name,title,items,notes=(),lines=None,points=None):
        p=pv.Plotter(off_screen=True,window_size=[1100,800]);p.set_background('white')
        for item,color,opacity in items:
            if isinstance(item,trimesh.Trimesh):item=polydata(item)
            p.add_mesh(item,color=color,opacity=opacity,smooth_shading=False)
        for a,b,color,width in lines or []:p.add_mesh(pv.Line(a,b),color=color,line_width=width)
        if points is not None and len(points):p.add_mesh(pv.PolyData(np.asarray(points)),color='#cb4545',point_size=7,render_points_as_spheres=True)
        direction=np.array([1.,-1.25,1.1]);direction/=np.linalg.norm(direction);p.camera_position=[center+direction*extent*2,center,[0,0,1]];p.enable_parallel_projection();p.camera.parallel_scale=extent*.66
        image=Image.fromarray(p.screenshot(return_img=True));p.close();annotate(image,title,notes).save(qc/name)
    base=[(mesh,'#399d63',1),(box,'#99aabb',.12)];ghost=[(mesh,'#399d63',.25),(box,'#99aabb',.08)]
    native=out/'slicer/bambu_native_current/actual_support_paths.npz';native_items=[(path_mesh(native),'#cb9242',1)] if native.exists() else []
    shot('01_frozen_vascular_input.png','冻结输入：血管、端口、2 mm 盒体；比例 1.0',base)
    shot('02_bambu_native_support.png','Bambu 自动支撑：只作比较基准',base+native_items)
    raw=np.load(out/'detector_points.npz')
    shot('03_freefloating_raw_detection.png',f"FreeFloating 原始检测：{len(raw['freefloating'])} 点",ghost,points=raw['freefloating'])
    shot('04_layerwise_fdm_detection.png',f"逐层新增承托不足：{len(raw['layer'])} 个区域样本",ghost,points=raw['layer'],notes=['区域完整多边形参与计算；这里仅显示代表位置'])
    long=[p['xyz'] for p in det['raw'] if 'SWC_LONG_SPAN' in p['detectors_triggered'] and p['risk_level']=='CRITICAL']
    shot('05_swc_long_span_detection.png',f"长悬空段：最长 {det['longest_dangerous_span_mm']:.2f} mm",ghost,points=long)
    shot('06_critical_intervals.png',f"四类证据合并：{len(det['critical_intervals'])} 个必须支撑区间",ghost,points=[p['xyz'] for p in det['raw'] if p['risk_level']=='CRITICAL'])
    shot('07_raw_tip_candidates.png',f"筛选前：{len(det['raw'])} 个下表面候选位置",ghost,points=[p['xyz'] for p in det['raw']])
    shot('08_pruned_tip_candidates.png',f"筛选后：{len(contacts)} 个计划接触位置",ghost,points=[p['xyz'] for p in contacts],notes=['计划接触点不等于已落地、可切片的显式支撑尖端'])
    topology=json.loads((out/'topology.json').read_text());keep=[]
    for xyz in topology['bifurcations']:
        p=trimesh.creation.icosphere(subdivisions=2,radius=2);p.apply_translation(xyz);keep.append((p,'#bd4b91',.3))
    near=[t['xyz'] for t in graph['tips'] if t['root_connected'] and t['distance_to_bifurcation']<2]
    shot('09_bifurcation_keepout.png',f"分叉邻域内已生成接触：{len(near)} 个",ghost+keep,points=near,notes=['红点为避让例外；具体原因逐项保存在接触表'])
    edges=[(graph['nodes'][e['start']]['xyz'],graph['nodes'][e['end']]['xyz'],'#bc7e27',2) for e in graph['edges']]
    shot('10_clever_support_graph.png',f"向下合并支撑图：{len(graph['merge_nodes'])} 个合并节点",ghost,lines=edges)
    load=[(a,b,c,max(2.,w)) for (a,b,c,_),w in zip(edges,[e['diameter_mm']*2 for e in graph['edges']])]
    shot('11_load_aware_trunks.png','承担更多分支的树干加粗；模型血管不变',ghost,lines=load)
    shot('12_stability_repairs.png',f"稳定性修复：{len(graph['stability_repairs'])} 项，未增加血管接触点",base+[(support,'#d89b3c',1)],notes=['倾斜合并失败时改用竖直支柱；研究短桥仍须实物试验'])
    shot('13_breakaway_roots.png',f"盒内底板可折断根部：{len(graph['roots'])} 个",[(box,'#99aabb',.1),(support,'#d89b3c',.45)],points=[r['estimated_break_location'] for r in graph['roots']],notes=[f"颈宽 {best['neck_width_mm']:.2f} mm；红点为预定折断位置；不代表折断力已验证"])
    colors={'taper':'#d34e4e','branch':'#498aac','pillar':'#d49a36','bridge':'#9461a5'}
    parts=[]
    from .support_custom_graph import edge_geometry
    from .support_breakaway import root_geometry
    for e in graph['edges']:
        m=edge_geometry(graph,e);parts.append((m,colors.get(e['kind'],'#d49a36'),1))
    for r in graph['roots']:
        m,_,_=root_geometry(r['xyz'][:2],r['xyz'][2],r['trunk_diameter_mm'],r['neck_width_mm'],summary['profile'],r['root_id'],r['trunk_load']);parts.append((m,'#81529c',1))
    kinds={kind:sum(e['kind']==kind for e in graph['edges']) for kind in ['branch','pillar','bridge']}
    shot('14_custom_support_only.png',f"显式支撑：{best['planned_tip_count']} 个尖端，{best['removal']['component_count']} 个组件",parts,notes=['红：尖端；蓝：合并分支；金：树干／立柱；紫：可折断底脚',f"本候选：合并节点 {len(graph['merge_nodes'])}；立柱 {kinds['pillar']}；短桥 {kinds['bridge']}"])
    missed=set((best.get('actual_coverage') or best['planned_coverage'])['uncovered_point_ids'])
    shot('15_custom_support_plus_vascular.png','冻结血管、端口与显式支撑',base+[(support,'#d49a36',1)],points=[p['xyz'] for p in det['raw'] if p['point_id'] in missed],notes=[f"关键弧长采样覆盖 {best['planned_coverage']['percent']:.2f}%；红点为缺少支撑的位置"])
    routes=[];ghosts=[]
    for row in best['removal']['rows']:
        if row['status']=='TRAPPED_SUPPORT_COMPONENT':
            path=folder/f"broken_component_{row['component_id']:03}.stl"
            if path.exists():ghosts.append((trimesh.load_mesh(path),'#cb4545',1))
        route=row.get('tilt_path')
        if route:
            poses=[np.array(x) for x in route['transforms']];routes.extend((a[:3,3],b[:3,3],'#2a89a5',3) for a,b in zip(poses[:-1],poses[1:]))
    shot('16_component_removal_paths.png',f"折颈后整组件移除：被困风险 {best['removal']['trapped_support_count']} 个",base+[(support,'#d49a36',.55)]+ghosts,lines=routes,
        notes=['逐组件测试上移、±10°/±20°倾斜、有限横移；FCL 距离界检查整个扫掠路径','红色组件：预设路径中未找到取出路线；仍需验证实际拆除顺序'])
    comp=summary['comparison'];bambu=next((r for r in comp if r['method']=='BAMBU_NATIVE_BASELINE'),{});custom=comp[-1]
    shot('_native_compare.png','Bambu 自动支撑（基准）',base+native_items,notes=[f"接触区间估计 {bambu.get('tip_count')}；面积 {bambu.get('contact_area_mm2',0):.2f} mm²",f"危险覆盖 {bambu.get('coverage_percent',0):.2f}%；体积 {bambu.get('support_volume_mm3',0):.2f} mm³"])
    area=custom.get('contact_area_mm2');shot('_custom_compare.png','显式支撑（失败项完整保留）',base+[(support,'#d49a36',1)],notes=[f"几何尖端 {custom['tip_count']}；实际面积 {area:.2f} mm²" if area is not None else f"几何尖端 {custom['tip_count']}；实际面积尚未切片核验",f"危险覆盖 {custom['coverage_percent']:.2f}%；体积 {custom['support_volume_mm3']:.2f} mm³"])
    a=Image.open(qc/'_native_compare.png');b=Image.open(qc/'_custom_compare.png');im=Image.new('RGB',(a.width+b.width,max(a.height,b.height)),'white');im.paste(a,(0,0));im.paste(b,(a.width,0));im.save(qc/'17_native_vs_custom_contact.png')
    actual=folder/'bambu_assembly/actual_custom_extrusion.npz';extra=[]
    if actual.exists() and best['slice'].get('slicer_success'):
        points=np.load(actual)['points'];extra=[(pv.PolyData(points),'#d49a36',1)]
    shot('18_real_bambu_auto_support_off_preview.png','Bambu 实际刀路：自动支撑关闭',base+extra,notes=[f"真实切片成功：{best['slice'].get('slicer_success',False)}；自动支撑关闭：{best['slice'].get('auto_support_off',False)}",'黄色为实际挤出轨迹采样，不是切片器原生截图'])
    report(out,det,contacts,graph,best,summary)
    stems=['frozen_vascular_input','bambu_native_support','freefloating_raw_detection','layerwise_fdm_detection','swc_long_span_detection','critical_intervals','raw_tip_candidates','pruned_tip_candidates','bifurcation_keepout','clever_support_graph','load_aware_trunks','stability_repairs','breakaway_roots','custom_support_only','custom_support_plus_vascular','component_removal_paths','native_vs_custom_contact','real_bambu_auto_support_off_preview']
    save(out/'visual_qc_manifest.json',dict(main_images=[str(qc/f'{i+1:02}_{name}.png') for i,name in enumerate(stems)],camera_center=center,parallel_scale=extent*.66,orientation_direction=[1,-1.25,1.1],common_scale=True))

def report(out,det,contacts,graph,best,summary):
    actual=best['slice'];cov=best.get('actual_coverage') or best['planned_coverage'];rooted=[p for p in graph['tips'] if p['root_connected']];counts=len(rooted);raw=len(det['raw']);reduction=100*(1-counts/raw)
    failures=best['hard_gate_failures'];area=actual.get('contact_area_mm2');pct=best.get('contact_area_reduction_percent',-1)
    answers=[
      '原 Bambu 自动树在可支撑下表面连续布置接触，导致多处细血管被包裹；减少密度又可能遗漏长悬空段。',
      f"四类证据合并为 {len(det['critical_intervals'])} 个必须支撑的区间，最长危险长段 {det['longest_dangerous_span_mm']:.2f} mm，包含中央血管及冻结端口。",
      f"FreeFloating CPU 移植保留 {det['freefloating_raw_count']} 个原始点，其中 {det['freefloating_mandatory_mapped']} 个通过下表面/逐层证据映射进入必须支撑要求；没有宣称跑通原始 OpenGL 应用。",
      '逐层 FDM 检测计算完整欠支撑多边形，可以识别仅看三角形角度或区域中心点遗漏的长而窄悬空区；结果见 layer_unsupported.csv。',
      f'下表面原始候选位置共 {raw} 个。',
      f'弧长筛选从 {raw} 个候选缩减到 {len(contacts)} 个计划点，减少 {100*(1-len(contacts)/raw):.2f}%；真正连接合法底脚的显式尖端为 {counts} 个。计划点中未连接的部分属于失败，不算优化收益。',
      f"支撑图中 {len(graph['merge_nodes'])} 个合并节点、{len(graph['roots'])} 个根部；折断根颈后有 {best['removal']['component_count']} 个独立组件。",
      f"分叉 2 mm 邻域内显式接触 {best['bifurcation_contacts']} 个；例外位置与原因在 selected_support_tips.csv 和 QC/09。",
      f"整体扫掠移除检查发现 {best['removal']['trapped_support_count']} 个被困组件；只记录预设路径内的几何结果，不代表实物拆除成功。",
      f'实际同口径刀路接触面积较 Bambu 基准下降 {pct:.2f}%，但当前覆盖不同，因此不能宣称已取得满足相同覆盖条件的优化收益。' if area is not None else '实际接触负担尚未取得可比切片结果，不能用名义尖端圆盘面积代替。',
      f"支撑自身逐层可打印检查：{'通过' if best['self_printability']['passed'] else '未通过'}；关闭自动支撑的独立支撑/合法底板切片：{best.get('support_only_slice',{}).get('slicer_success',False)}。切片成功本身不保证能实际打印。",
      f"Bambu X1 Carbon 0.2 mm 实际整件切片：{actual.get('slicer_success',False)}；Auto Support OFF：{actual.get('auto_support_off',False)}；危险覆盖 {cov['percent']:.2f}%。"
    ]
    lines=['# BG001 RMCA 显式血管支撑实现与验证','',f"状态：**{summary['status']}**。当前正式算法为 CUSTOM_LITERATURE_AWARE_SUPPORT；原生 Bambu、Orca、Guided Orca 全部只作 benchmark。",'']
    lines += [f'{i+1}. {a}' for i,a in enumerate(answers)]
    domains=summary.get('critical_coverage_by_domain',{})
    lines += ['','## 实际覆盖与比较口径','',
        f"完整冻结模型的关键弧长采样覆盖率为 {cov['percent']:.2f}%，已覆盖 {cov['covered_requirements']} / {cov['critical_requirements']} 个必须承托的采样位置；最长仍未覆盖的连续弧段为 {cov['max_critical_unsupported_span']:.2f} mm。100% 是软件弧长覆盖门槛，不是材料强度结论。"]
    lines.append('')
    for label,name in [('anatomical_vascular','中央解剖血管'),('artificial_ports','人工端口')]:
        d=domains.get(label)
        if d:lines.append(f"- {name}：覆盖 {d['coverage_percent']:.2f}%，{d['uncovered_samples']} / {d['required_samples']} 个危险采样仍未覆盖。")
    lines += ['', '基准的 tip count 是从真实刀路估计的弧长接触区间数，显式方案是可定位的实体尖端数，两者不能直接作为同一定义的数量降幅。接触面积均用真实血管三角形及相同刀路距离规则计算；基准的覆盖已按本轮融合检测重新计算。基准的拆除项仍是旧探针检查，不等同于本轮整组件扫掠移除。',
        '本轮移除测试的障碍包括冻结血管、端口和盒体。各支撑组件分别测试，暂未把其它未拆组件作为障碍，也没有验证钳子操作空间或真实折断力；实际拆除顺序仍需人工与试片核验。',
        f"本轮检测、点筛选及三候选几何/移除检查累计 {summary.get('core_geometry_search_seconds_including_geometry_and_removal_qc',0):.2f} 秒，不包括外部切片。4/5/6 mm 阈值敏感性见 tables/span_sensitivity.csv；5 mm 是本工艺的待验证工程参数，不是通用 FDM 定律。"]
    warnings=actual.get('critical_warnings',[])
    if warnings:lines += ['', '真实装配切片保留了以下警告（未按“切片成功”忽略）：','', *['- '+str(w) for w in warnings], '', '显式尖端的有意空气间隙也可能触发悬空区域警告，不能仅凭警告定位具体缺陷；危险段覆盖另由独立检测判定。']
    profile=summary['profile']
    lines += ['', '## 工艺与真实验证记录','',
        f"读取本机 {profile['printer']}、ABS Vessel、{profile['process_name']}：喷嘴 {profile['nozzle_diameter']:.2f} mm，层高 {profile['layer_height']:.2f} mm，线宽 {profile['line_width']:.2f} mm。目标顶部间隙 {profile['support_top_gap']:.2f} mm，非尖端支撑与血管的目标距离 {profile['support_xy_gap']:.2f} mm。数值来自实际保存的工艺配置。",
        f"诊断候选支撑体积 {best['support_volume_mm3']:.2f} mm³；非法血管/端口交集 {best['union_audit']['vascular_and_ports_overlap_volume_mm3']:.6f} mm³；非法侧壁/底板交集 {best['union_audit']['illegal_wall_or_bottom_overlap_volume_mm3']:.6f} mm³；原模型损失体积 {best['union_audit']['frozen_mold_removed_volume_mm3']:.6f} mm³。仅批准脚垫与盒底相连。",
        f"实际装配刀路保留 {actual.get('actual_preserved_tip_count','尚未验证')} 个尖端；原生自动支撑刀路数 {actual.get('native_support_toolpath_count','尚未验证')}。本候选最终合并节点 {len(graph['merge_nodes'])} 个：未通过自身打印检查的合并已退回立柱，未将未实现的树合并效果当作成果。"]
    tests=out/'test_results.xml'
    if tests.exists():
        import xml.etree.ElementTree as ET
        suites=list(ET.parse(tests).getroot().iter('testsuite'))
        total=sum(int(s.get('tests',0)) for s in suites);failed=sum(int(s.get('failures',0))+int(s.get('errors',0)) for s in suites)
        lines.append(f"永久测试：{total} 项，失败或错误 {failed} 项；详细结果为 test_results.xml 与 test_run.log。")
    coupon_result=out/'coupon/slice_validation.json'
    if coupon_result.exists():
        cr=json.loads(coupon_result.read_text())
        lines.append(f"27 格试片实际切片成功：{cr.get('slicer_success',False)}；自动支撑关闭：{cr.get('auto_support_off',False)}；全部尖端和窄颈保留：{cr.get('all_coupon_features_preserved',False)}。物理打印、折断和拆除尚未执行。")
    lines += ['','## 当前没有通过的条件','',*['- '+v for v in failures], '',
        '全部失败指标均参与最终状态，不能用接触面积下降掩盖覆盖不足。当前文件是可复核的诊断候选；只有所有门槛通过才会写入 production_selected。',
        '', '## 为什么盒外端口可能与新根部约束冲突','',
        '支撑根部只允许位于盒内底板。冻结盒体具有完整侧壁，部分端口在盒外且低于盒顶。向下生长且自身可打印的支撑无法从盒内越过实心侧壁去承托这些位置；在盒外落床、以侧壁为根或删去这些危险要求都违反本轮约束。代码逐项记录 EXTERIOR_TIP_NO_INNER_BOTTOM_ROUTE，并保持其覆盖失败，不移动任何血管或改盒。',
        '', '## 可折断底脚与下一步试验','',
        f"每个底脚由受控圆形脚垫、窄颈和向上渐扩的树干组成，只有脚垫与盒内底板有意相交。本候选颈宽 {best['neck_width_mm']:.2f} mm；所有接触面积、颈高、载荷和预定断裂位置见 support_roots.csv。根颈以实际线宽的整数倍派生，未计算或宣称经过验证的折断力。",
        '下一步先做 coupon，暂不推荐打印当前未通过门槛的完整模型。27 格试片覆盖三种颈宽、三种尖端宽度以及 0.15/0.20/0.25 mm 几何间隙；模拟底板、渐细支撑与两端有支柱的小直径水平 ABS 圆杆。详见 coupon/coupon_test_plan.md。',
        '', '## 这次实际利用了哪些开源成果','',
        'FreeFloating 的 MIT 有向图覆盖和局部低点算法用于检测，保存未筛选点及每点证据；原应用编译依赖问题保留了日志。',
        'Clever Support 的作者下载包已从官方服务器取得并检验。其内容只有 14 个 OBJ 示例而非程序源码；论文第 5、6 节的向下锥体合并、支撑图与几何稳定性流程已实际阅读。不能把此包描述为已经复用的官方源代码。',
        'meshy 的 MIT 实现提供了由高到低合并、碰撞回退、渐细连接和 sqrt(load) 加粗规则的可审计参考；版权及源文件保留。当前适配禁止原算法可能使用的血管/模型表面根，增加受控底脚、有限路径拆除和切片特征保留检查。',
        'Orca 和 Bambu 的自动树只作对照；Bambu 在正式验证中关闭自动支撑，只将已经生成的显式实体转换为 X1C 刀路。未复制 AGPL 大段源码。短桥是论文思想的独立几何候选，仍须实物试验。',
        '', '## 输入、输出与复现','',
        '`./run_python.sh s1-8_vascular_support_planner.py --stage all`',
        '', '显式 STL 与 VTP 位于 geometry/；可折断试片位于 coupon/；各候选参数、几何检查、真实切片与日志位于 custom_candidates/；旧结果完整保留于 benchmark/guided_v1/，不参与生产选择。',
        f"冻结输入 SHA256：{summary['source_protection']['checked_files']} 个文件，变化 {len(summary['source_protection']['changed'])} 个。scale_factor = 1.0；未上传或启动打印。",'',
        '最大剩余风险包括必要位置无法由合法根部连接、整棵支撑不能取出，以及 ABS 收缩、细颈稳定性和真实折断损伤；这些风险必须逐项解决并通过试片与实物实验，不能由“切片完成”代替。']
    (out/'vascular_support_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
