"""Evidence figures and Chinese closure report; never paint failures green."""
from pathlib import Path
import json,pickle,shutil
import numpy as np
from .support_inputs import save

NAMES=['01_baseline_uncovered_points','02_coverage_repairs','03_100pct_actual_coverage','04_support_self_closure','05_support_segmentation','06_real_components_after_cut','07_weak_links','08_last_contact_pieces','09_bifurcation_support_detail','10_wrap_angle_map','11_initial_tool_access','12_sequential_tool_access','13_removal_load_paths','14_safe_removal_sequence','15_final_support_only','16_final_support_plus_vascular','17_real_Bambu_preview']


def render(base,out,g,record,e,removal,iterations):
    import pyvista as pv
    import trimesh,shapely
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out);qc=out/'QC_SAFE';folder=Path(record['folder']);v=g['vascular'];support=g['support'];c=json.loads((folder/'coverage.json').read_text())
    def mesh(m):return pv.PolyData(m.vertices,np.c_[np.full(len(m.faces),3),m.faces].ravel())
    def shot(i,title,core=True,scaffold=False,points=None,colors=None,labels=None,lines=None,parts=None):
        p=pv.Plotter(off_screen=True,window_size=(1400,1000));p.set_background('#f6f8fa')
        if core:p.add_mesh(mesh(v),color='#ad5b50',opacity=.25 if i in [0,2] else .7,smooth_shading=True)
        if scaffold:p.add_mesh(mesh(support),color='#739caf',opacity=.45)
        if parts:
            palette=['#3182bd','#31a354','#756bb1','#de2d26','#fd8d3c']
            for j,part in enumerate(parts):p.add_mesh(mesh(part),color=palette[j%len(palette)],opacity=.75)
        if points is not None and len(points):
            if colors is None:p.add_points(np.asarray(points),color='#ff2020',point_size=13,render_points_as_spheres=True)
            else:p.add_points(np.asarray(points),scalars=np.asarray(colors),rgb=True,point_size=9,render_points_as_spheres=True)
            if i==2 and colors is not None:
                bad=np.asarray(colors)[:,0]>200
                if bad.any():p.add_points(np.asarray(points)[bad],color='#ff2020',point_size=14,render_points_as_spheres=True)
        if labels:p.add_point_labels(np.asarray([x[0] for x in labels]),[x[1] for x in labels],font_size=11,always_visible=True,shape_opacity=.2)
        if lines:p.add_lines(np.asarray(lines).reshape(-1,3),color='#df6b13',width=3)
        p.add_text(NAMES[i]+'\n'+title,font_size=12,color='#283c48');p.add_axes();p.show_bounds(grid='back',xtitle='X (mm)',ytitle='Y (mm)',ztitle='Z (mm)',n_xlabels=3,n_ylabels=3,n_zlabels=3,font_size=11,fmt='%.1f');p.view_isometric();p.reset_camera();p.screenshot(qc/(NAMES[i]+'.png'));p.close()
    baseline=json.loads((out/'reports/baseline/coverage.json').read_text());points=[r['lower_surface'] for r in baseline['rows'] if not r['supported'] and r['lower_surface'] is not None]
    # Baseline coordinates are transformed into the chosen diagnostic frame.
    bs=json.loads((out/'reports/baseline/vascular_core_support_summary.json').read_text());source_M=np.asarray(bs['orientation']['transform']);M=np.asarray(record['orientation']['transform']);T=M@np.linalg.inv(source_M)
    points=np.asarray(points)@T[:3,:3].T+T[:3,3] if points else []
    shot(0,f"Baseline: {baseline['unsupported_samples']} failed samples",points=points)
    repair=folder/'local_repairs.json';ops=json.loads(repair.read_text())['operations'] if repair.exists() else []
    shot(1,f'{len(ops)} measured local operations in selected diagnostic',scaffold=True,points=[r['point'] for r in ops])
    rows=[r for r in c['rows'] if r['lower_surface'] is not None];colors=[[35,170,90] if r['supported'] else [255,20,20] for r in rows]
    shot(2,f"ACTUAL {c['centerline_underbody_percent']:.6f}% | FAIL {c['unsupported_samples']}/{c['sample_count']} | no production selection",points=[r['lower_surface'] for r in rows],colors=colors)
    sc=json.loads((folder/'support_self.json').read_text());q=[[*shapely.from_wkt(r['wkt']).representative_point().coords[0],r['z_mm']] for r in sc['rows']]
    shot(3,f"Actual self-unsupported regions: {sc['self_unsupported_regions']}",core=False,scaffold=True,points=q)
    shot(4,f"Planned segments: {len(g['segments'])}",scaffold=True,labels=[(np.r_[s['xy'][len(s['xy'])//2],s['top'][len(s['top'])//2]],s['segment_id']) for s in g['segments']])
    parts=[trimesh.load_mesh(p) for p in sorted((out/'virtual_components').glob('G*.stl'))]
    shot(5,f"ACTUAL virtual-cut components: {len(parts)}; only certified cuts simulated",parts=parts)
    shot(6,'Weak-link geometry optimization BLOCKED until coverage closure',scaffold=True)
    shot(7,'Current remaining contact groups; long pieces are NOT safe final pieces',parts=parts,labels=[(r['contact_point'],r['segment_id']+f" {r['last_contact_piece_length_mm']:.1f}mm") for r in removal['rows']])
    forks=np.unique(np.asarray([b['xyz'][0] for b in g['branches'] if not b['is_port']][1:]).reshape(-1,3),axis=0)
    fork_labels=[(q,f'BIF{i+1} / release NOT certified') for i,q in enumerate(forks)]
    fork_labels += [(cut['point'],cut['cut_id']) for cut in removal['cut_checks'] if len(forks) and np.min(np.linalg.norm(forks-np.asarray(cut['point']),axis=1))<3]
    shot(8,'All bifurcation locations and nearby proposed cuts; no approved release direction',scaffold=True,points=forks,labels=fork_labels)
    shot(9,'Local angular wrap proxy on remaining contact groups',scaffold=True,labels=[(r['contact_point'],r['segment_id']+f" {r['wrap_angle_deg']:.1f}deg") for r in removal['rows']])
    cuts={r['cut_id']:r for r in removal['cut_checks']};initial=removal['sequence']['initial_accessible']
    shot(10,f'Initially reachable cuts: {len(initial)}/{len(cuts)}',scaffold=True,points=[cuts[k]['point'] for k in initial])
    fig,ax=plt.subplots(figsize=(11,5));rounds=removal['sequence']['rounds'];ax.bar([r['round'] for r in rounds],[len(r['new_accessible_cuts']) for r in rounds],label='Newly reachable cuts');ax.set_xlabel('Sequential round');ax.set_ylabel('Count');ax.legend();ax.set_title('Obstacles are removed only after successful cut and swept release');fig.tight_layout();fig.savefig(qc/(NAMES[11]+'.png'),dpi=140);plt.close(fig)
    load_lines=[];load_labels=[]
    for r in removal['rows']:
        force=np.asarray(r['force_point']);contact=np.asarray(r['contact_point']);d=np.asarray(r['removal_direction']);d=d/np.linalg.norm(d)
        foot=force+d*np.dot(contact-force,d)
        load_lines.extend([[force,contact],[contact,foot],[force,force+3*d]])
        load_labels.extend([(force,r['segment_id']+f" lever {r['lever_arm_to_vessel_mm']:.2f}mm"),(contact,'CONTACT'),(foot,'PERPENDICULAR ARM')])
    load_labels.extend([(x['center'],x['link_id']) for x in g.get('weak_links',[])])
    shot(12,'Force / contact / perpendicular lever arm / proposed direction; no stress prediction',scaffold=True,lines=load_lines,labels=load_labels)
    shot(13,'No approved production sequence' if not removal['sequence']['complete'] else 'Geometric sequence only; physical forces uncalibrated',scaffold=True,labels=[(r['contact_point'],r['segment_id']+(' BLOCKED' if r['sequence_step'] is None else f" #{r['sequence_step']}")) for r in removal['rows']])
    shot(14,'Custom support only',core=False,scaffold=True);shot(15,'Frozen core + separate support | diagnostic',scaffold=True)
    shutil.copy2(folder/'bambu/native_preview.png',qc/(NAMES[16]+'.png'))
    save(out/'reports/QC_SAFE_manifest.json',dict(candidate=record['candidate_id'],figures=[str(qc/(n+'.png')) for n in NAMES],all_exist=all((qc/(n+'.png')).exists() for n in NAMES),green_only_claim=False))


def report(base,out,s,r,e,removal):
    out=Path(out);c=s['coverage'];l=s['layer_coverage'];selfcheck=s['support_self'];status=s['status'];rows=removal['rows']
    near=[r for r in rows if r['distance_to_bifurcation_mm'] is not None and r['distance_to_bifurcation_mm']<2]
    text=['# s1-11 实际承托闭合与安全拆除检查','',f'状态：**{status}**。生产候选：'+('无。' if s['production_selected'] is None else s['production_selected'])+' 当前文件仅用于诊断。','',
          f"1. 中心线实际承托是否 100%：{'是' if c['centerline_underbody_percent']==100 else '否'}，{c['centerline_underbody_percent']:.8f}%。",
          f"2. 延长端口是否 100%：{'是' if c['port_underbody_percent']==100 else '否'}，{c['port_underbody_percent']:.8f}%。",
          f"3. 完整逐层区域是否 100%：{'是' if l['percent']==100 and l['passed'] else '否'}，{l['percent']:.8f}%。",
          f"4. 最长未通过承托检查的弧长：{c['max_uncovered_arc_mm']:.8f} mm；失败采样点 {c['unsupported_samples']} 个。过近顶面也作为失败，不能用较低幕墙层代替。",
          f"5. 支撑自身未通过上一层承托检查的区域：{selfcheck['self_unsupported_regions']} 个。",'',
          '## 局部修补与停止条件','',f"本轮读取既有报告、逐点覆盖、多边形、拆除及空移记录；保留原姿态搜索结果，只处理 O002 与既有两个候选姿态。完成 {s['iteration_count']} 次有真实切片指标的修补。每条路线最多 8 次；仅在连续两次缺点数和缺口面积均不改善时提前停止。运行异常单独记录，不算作收敛。",'',
          '所有变化仅发生在支撑的实际缺口附近，尺寸采用当前线宽或层高的倍数，网格顶部的半层偏置来自切片平面位置；没有放宽覆盖阈值或改动血管。缺口数据和每次变化见 `support_closure_iterations.csv` 及各候选的 `local_repairs.json`。', '',
          f"本次诊断展示选择 `{s['diagnostic_selected']}`，排序首先比较三类覆盖率中最低的一项。它只用于呈现失败证据，不是已验收的生产候选；局部迭代是否接受仍要求缺点数与缺口面积均不退步。",'',
          '覆盖未达到三个 100% 且支撑自身未闭合时，不开始修改拆除几何。下面仅为当前诊断模型的只读虚拟剪切、载荷路径与顺序检查；没有将其冒充生产拆卸方案。','',
          f"6. 计划支撑段数：{removal['planned_segment_count']}。",
          f"7. 有限虚拟剪切后的实际组件数：{removal['actual_component_count']}。",
          f"8. 最大组件接触弧长：{removal['max_removal_component_arc_span_mm']:.3f} mm；包围盒各轴上限：{removal['max_removal_component_bbox']} mm。若仍为巨大连通组件，则不满足安全拆卸要求。",
          f"9. 候选剪切 {removal['cut_count']} 次，其中几何安全检查通过 {removal['safe_cut_count']} 次；不是已批准操作次数。",
          '10. 不要求所有剪口初始可达。每次剪切前保留尚未切除的材料，每次移除前保留其它组件，随后更新障碍物。','',
          '## 顺序与血管载荷风险','']
    for rr in removal['sequence']['rounds']:text += [f"- 第 {rr['round']} 轮：新增可达剪口 {rr['new_accessible_cuts']}；移除 {rr['removed']}；剩余 {rr['remaining']}。"]
    text += ['',
             '11. 不允许通过血管弯曲来获得移除路径。未找到刚体路径表示“未能证明可拆”，不等于已证明必须弯曲血管。',
             f"12. 检出的闭合截面套环：{removal['closed_cages']} 个。",
             f"13. 最大局部包覆角：{removal['max_wrap_angle']:.3f}°。",
             f"14. 超过 120° 的组件：{removal['wrap_over_120']} 个。",
             f"15. 超过 150° 的组件：{removal['wrap_over_150']} 个，默认失败；本轮未启用人工审核例外。",
             f"16. 分叉附近当前接触组件最大弧长：{max((r['last_contact_piece_length_mm'] for r in near),default=0):.3f} mm。",
             f"17. 最大几何力臂：{removal['max_lever_arm_mm']:.3f} mm。",
             '18. Weak-link ratio：尚未优化或通过真实切片确认，保留为空。实际血管截面由冻结 STL 截取；不使用半径推算，也不将面积比解释为断裂力比。','',
             '## 验证边界','',f"源保护核对 {s['source_protection']['checked_files']} 个文件，改变 {len(s['source_protection']['changed'])} 个。真实 Bambu 切片成功：{e['slice']['slicer_success']}；自动支撑关闭：{e['auto_support_off']}。未解决关键警告：{e['unresolved_warnings']}。",'',
             '软件只能检查几何、实际挤出和碰撞代理，不能保证 ABS 一定先断或血管一定不会断。当前不建议完整实物打印。拆除原则仍为先剪断、再移除长幕墙、最后处理短接触块；未认证的位置不能直接横向掰动。','',
             '## 未通过的硬门槛','',*['- '+x for x in s['failures']],'',
             '## 文件入口','', '`selected/`：血管 STL、自定义支撑 STL、仅拼接的审阅 STL 和原生 `X1C_project.3mf`。血管与支撑保持两个 normal part。', '',
             '`QC_SAFE/`：17 张图。图 03 文件名沿用任务要求，但只要失败点存在就显示红点及真实百分比，不能画成全绿。','',
             '`safe_removal_sequence.json`：生产方案为空；诊断检查与未完成原因单独记录。', '',
             '```bash','cd /home/lzy/projects/vascular_printing','./run_python.sh s1-11_vascular_core_support.py --stage closure','```','']
    baseline=json.loads((out/'reports/baseline/vascular_core_support_summary.json').read_text())
    bc=json.loads((out/'reports/baseline/coverage.json').read_text());bl=json.loads((out/'reports/baseline/layer_coverage.json').read_text())
    text += ['## 基线对照与真实运行证据','',f"原基线：`{baseline['selected_candidate']}`。各项独立列出，不能以逐层覆盖改善替代中心线及支撑自身合格。",'',
             '|指标|本轮开始时的基线|最终诊断候选|','|---|---:|---:|',
             f"|中心线覆盖 / %|{bc['centerline_underbody_percent']:.8f}|{c['centerline_underbody_percent']:.8f}|",
             f"|端口覆盖 / %|{bc['port_underbody_percent']:.8f}|{c['port_underbody_percent']:.8f}|",
             f"|逐层覆盖 / %|{bl['percent']:.8f}|{l['percent']:.8f}|",
             f"|失败采样点|{bc['unsupported_samples']}|{c['unsupported_samples']}|",
             f"|逐层缺口面积累加 / mm²|{bl['uncovered_area_sum_mm2']:.8f}|{l['uncovered_area_sum_mm2']:.8f}|",'',
             '停止记录（运行异常与真实停滞分开）：','', '|路线|轮次|原因|','|---|---:|---|']
    for row in s['stop_conditions']:text.append(f"|{row['baseline']}|{row['iteration']}|{row['reason']}"+('：'+row['error'] if row.get('error') else '')+'|')
    import csv
    attempts=list(csv.DictReader((out/'support_closure_iterations.csv').open()))
    accepted=sum(row['accepted'].lower()=='true' for row in attempts)
    actual_c=json.loads((Path(r['folder'])/'coverage.json').read_text())
    bad=[row for row in actual_c['rows'] if not row['supported']]
    too_close=sum(bool(row.get('too_close_support_layer_tops')) for row in bad)
    text += ['',f"局部修补实际接受 {accepted} 次，未通过回退条件的试验全部保留供追踪。当前诊断的 {len(bad)} 个失败点中，{too_close} 个存在过近的顶面，其余 {len(bad)-too_close} 个缺少有效顶面承托；不能把这两类都描述成完全裸悬空。",'',
             f"支撑自身悬空区域：原基线 {baseline['support']['self_unsupported_regions']} 个，当前诊断 {selfcheck['self_unsupported_regions']} 个。任何单项改善均不能抵消这一未通过状态。"]
    tests=out/'validation/test_results.json'
    if tests.exists():
        t=json.loads(tests.read_text());text += ['',f"永久测试：{t['tests']} 项，失败 {t['failures']}，错误 {t['errors']}，跳过 {t['skipped']}。测试通过只说明检查逻辑按预期执行，不代表制造门槛通过。"]
    for folder in ['actual_thin_web','actual_thin_web_C3']:
        file=out/'validation'/folder/'thin_web_result.json'
        if file.exists():
            a=json.loads(file.read_text());text += ['',f"真实合成连接片 `{folder}`：连接片保留判据 {'PASS' if a['actual_preservation']['passed'] else 'FAIL'}。详细逐层证据见 `validation/{folder}/`；连接片保留不能替代完整覆盖及支撑自身检查。"]
    travel=s['travel']
    text += ['', 'Bambu 的对象级悬浮警告通过候选目录中的 `warning_localization.json` 对应到实际失败采样点、逐层缺口及支撑自身悬空区域；这些区域未闭合前均保留 UNRESOLVED。', '',
             f"实际刀路空移审计：{travel['travel_segments']} 段空移，{travel['near_miss_segments']} 段触发近碰代理。结果使用逐时刻已打印几何及喷头代理，属于几何风险筛查，不是实测碰撞或已校准的失效概率。明细见所选候选的 `travel_near_miss.json`。", '']
    (out/'vascular_core_support_safe_removal_report.md').write_text('\n'.join(text))
    guide=['# SAFE SUPPORT REMOVAL GUIDE','',f'状态：**{status}**。','', '当前没有通过全部门槛的 selected production geometry，因此本文件不提供可直接照做的生产剪切步骤。以下按真实组件编号列出阻碍，避免把未认证方案当作操作指南。','','|组件 / 血管分支|风险|拟夹持位置|拟剪切位置|禁止施力|','|---|---|---|---|---|']
    for r in rows:guide.append(f"|{r['segment_id']} / {r['contact_branch']}|接触长 {r['contact_length_mm']:.2f} mm；力臂 {r['lever_arm_to_vessel_mm']:.2f} mm；{r['wrap_status']}|{np.round(r['force_point'],2).tolist()}（仅支撑侧，未批准）|{r['planned_cut_ids']}（须逐一认证）|沿血管轴强拖、分叉横拉、抓长幕墙直接掰离|")
    guide+=['','操作顺序必须先 A 剪断主体、B 取走长幕墙、C 暴露内部剪口、D 去除短接触片。当前缺少完整合法顺序的组件不得通过弯折血管强行取出。']
    (out/'SAFE_SUPPORT_REMOVAL_GUIDE.md').write_text('\n'.join(guide)+'\n')
