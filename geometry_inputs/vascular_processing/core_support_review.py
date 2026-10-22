"""Evidence figures and Chinese reporting, with failed gates kept visible."""
from __future__ import annotations
import json,pickle,shutil
from pathlib import Path
import numpy as np
from .support_inputs import save

NAMES=['01_frozen_vascular_core','02_orientation_candidates','03_selected_orientation','04_all_lower_envelope_targets','05_continuous_support_roofs','06_vertical_curtains','07_projection_collision_repairs','08_bifurcation_support','09_port_extension_support','10_cut_notches','11_support_only_front','12_support_only_side','13_support_only_top','14_support_plus_core','15_wrap_angle_audit','16_actual_Gcode_underbody_coverage','17_uncovered_arc_map','18_support_self_printability','19_removal_access','20_removal_sequence','21_travel_near_miss','22_final_Bambu_preview']


def render(out,record,evaluation,removal,travel,ranking):
    import pyvista as pv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import shapely
    import trimesh
    out=Path(out);folder=Path(record['folder']);qc=out/'QC';qc.mkdir(exist_ok=True);g=pickle.load((folder/'geometry.pkl').open('rb'));vascular=g['vascular'];support=g['support'];coverage=json.loads((folder/'coverage.json').read_text());wrap=json.loads((folder/'wrap_audit.json').read_text());selfcheck=json.loads((folder/'support_self.json').read_text())
    def poly(m):return pv.PolyData(m.vertices,np.c_[np.full(len(m.faces),3),m.faces].ravel())
    def shot(name,core=True,scaffold=False,part=None,points=None,colors=None,view='iso',lines=None,labels=None):
        p=pv.Plotter(off_screen=True,window_size=(1400,1000));p.set_background('#f6f8fa')
        if core:p.add_mesh(poly(vascular),color='#ad5b50',opacity=.20 if name in NAMES[15:17] else .82,smooth_shading=True)
        if scaffold:p.add_mesh(poly(support),color='#739caf',opacity=.55)
        if part is not None and len(part.faces):p.add_mesh(poly(part),color='#3a92aa',opacity=.85)
        if points is not None and len(points):
            if colors is None:p.add_points(np.asarray(points),color='#204f77',point_size=8,render_points_as_spheres=True)
            else:p.add_points(np.asarray(points),scalars=np.asarray(colors),rgb=True,point_size=8,render_points_as_spheres=True)
        if lines is not None and len(lines):p.add_lines(np.asarray(lines).reshape(-1,3),color='#ca3434',width=2)
        if labels:p.add_point_labels(np.asarray([x[0] for x in labels]),[x[1] for x in labels],font_size=12,shape_opacity=.3,always_visible=True)
        p.add_text(name+'\n'+record['candidate_id']+' | software inspection only',font_size=12,color='#283c48');p.add_axes();p.show_bounds(grid='back',location='outer',xtitle='X (mm)',ytitle='Y (mm)',ztitle='Z (mm)',n_xlabels=3,n_ylabels=3,n_zlabels=3,fmt='%.1f',font_size=12,padding=.08)
        if name==NAMES[15]:
            failed_points=[r['lower_surface'] for r in coverage['rows'] if not r['supported'] and r['lower_surface'] is not None]
            if failed_points:p.add_points(np.asarray(failed_points),color='#ff2020',point_size=14,render_points_as_spheres=True)
            p.add_text('GREEN: correct-gap actual roof | RED: missing roof or invalid gap\n'+str(coverage['unsupported_samples'])+' / '+str(coverage['sample_count'])+' samples FAIL',position='lower_left',font_size=11,color='#283c48')
        if view=='front':p.view_xz()
        elif view=='side':p.view_yz()
        elif view=='top':p.view_xy()
        else:p.view_isometric()
        p.reset_camera();p.screenshot(qc/(name+'.png'));p.close()
    shot(NAMES[0]);fig,ax=plt.subplots(figsize=(12,5));ax.bar(np.arange(len(ranking)),[r['estimated_curtain_volume_mm3'] for r in ranking]);ax.set_xticks(range(len(ranking)),[r['orientation_id'] for r in ranking],rotation=90,fontsize=7);ax.set_ylabel('Coarse support volume proxy (mm3)');ax.set_title('48 deterministic rigid orientations');fig.tight_layout();fig.savefig(qc/(NAMES[1]+'.png'),dpi=150);plt.close(fig)
    shot(NAMES[2]);lower=np.vstack([b['lower'][b['valid']] for b in g['branches']]);shot(NAMES[3],points=lower);shot(NAMES[4],part=g['roofs']);shot(NAMES[5],part=g['curtains'])
    repairs=[s for s in g['segments'] if s['route']!='DIRECT'];q=np.vstack([np.c_[s['xy'],s['top']] for s in repairs]) if repairs else np.empty((0,3));shot(NAMES[6],scaffold=True,points=q)
    forks=np.array([b['xyz'][0] for b in g['branches'] if not b['is_port']]);shot(NAMES[7],scaffold=True,points=forks);ports=np.vstack([b['lower'] for b in g['branches'] if b['is_port']]);shot(NAMES[8],scaffold=True,points=ports)
    labels=[(n['point'],n['notch_id']) for n in g['notches']];shot(NAMES[9],scaffold=True,labels=labels)
    for name,view in zip(NAMES[10:13],['front','side','top']):shot(name,core=False,scaffold=True,view=view)
    shot(NAMES[13],scaffold=True);fig,ax=plt.subplots(figsize=(10,5));ax.hist([r['angular_coverage_deg'] for r in wrap['rows']],bins=35);ax.axvline(120,color='red',label='Hard limit 120 deg');ax.set_xlabel('Near-surface support angular extent (deg)');ax.set_ylabel('Cross-section samples');ax.legend();fig.tight_layout();fig.savefig(qc/(NAMES[14]+'.png'),dpi=150);plt.close(fig)
    rows=[r for r in coverage['rows'] if r['lower_surface'] is not None];points=np.array([r['lower_surface'] for r in rows]);colors=np.array([[35,170,90] if r['supported'] else [225,30,30] for r in rows]);shot(NAMES[15],points=points,colors=colors)
    bad=[r['lower_surface'] for r in rows if not r['supported']];shot(NAMES[16],points=np.asarray(bad),colors=np.tile([225,30,30],(len(bad),1)))
    selfpoints=[[*np.asarray(shapely.from_wkt(r['wkt']).representative_point().coords[0]),r['z_mm']] for r in selfcheck['rows']];shot(NAMES[17],core=False,scaffold=True,points=selfpoints,colors=np.tile([225,30,30],(len(selfpoints),1)))
    accesslines=[[a['point'],np.asarray(a['point'])+np.asarray(a['direction'])*6] for a in removal['access'] if a['passed']];accesslabels=[(a['point'],a['notch_id']+(' OK' if a['passed'] else ' BLOCKED')) for a in removal['access']];shot(NAMES[18],scaffold=True,labels=accesslabels,lines=accesslines)
    fraglabels=[(r['centroid'],r['fragment_id']+(' #'+str(r['suggested_order']) if r['passed'] else ' BLOCKED')) for r in removal['fragments']];shot(NAMES[19],scaffold=True,labels=fraglabels)
    shot(NAMES[20],scaffold=True,lines=[[r['start'],r['end']] for r in travel['rows']]);shutil.copy2(folder/'bambu/native_preview.png',qc/(NAMES[21]+'.png'));save(out/'reports/QC_manifest.json',dict(figures=[str(qc/(n+'.png')) for n in NAMES],source_candidate=record['candidate_id'],all_exist=all((qc/(n+'.png')).exists() for n in NAMES)))


def report(out,s,r,e,removal):
    out=Path(out);c=s['coverage'];sp=s['support'];b=s['bambu'];rem=s['removal'];failed=s['failures'];good=not failed;status=s['status'];folder=Path(r['folder']);M=np.asarray(s['orientation']['transform']);normal=M[:3,:3].T@np.array([0.,0.,1.])
    answers=[('本轮范围','只新增 s1-11 的血管打印支撑及其检查；s1-10 外模、插片、PDMS、溶剂通路、裁切系统均为只读。'),('血管是否改变',f"源 STL SHA256 为 `{s['input']['vascular_hash']}`，开始/结束 {s['source_protection']['mismatch_count']} 项不一致。打印副本只做刚体变换，缩放为 1。"),('当前姿态',f"{s['orientation']['selected']}；原始坐标中的打印向上方向为 {np.round(normal,5).tolist()}。"),('选择原因','先比较打印范围、支撑体积、高度和投影相交，再按真实承托和硬门槛筛选。当前是未通过门槛的诊断候选，不是生产选择。' if not good else '通过全部软件硬门槛后，比较旅行风险、支撑体积与预计时间。'),('支撑范围',f"{c['branch_count']} 条解剖分支、{c['port_count']} 条延长端口，共 {c['sample_count']} 个真实 STL 下表面射线样本。"),('是否实现 100%',f"中心线下侧 {c['centerline_underbody_percent']:.8f}%；端口 {c['port_underbody_percent']:.8f}%；逐层区域 {s['layer_coverage']['percent']:.8f}%。"),('最长未支撑弧长',f"{c['max_uncovered_arc_mm']:.8f} mm；未覆盖样本 {c['unsupported_samples']} 个。"),('支撑是否穿过血管',f"实际输出 STL 的交叠体积 {sp['intersection_mm3']:.10g} mm³；门槛是否通过：{s['gates']['NO_VASCULAR_INTERSECTION']}。"),('支撑自身可打印性',f"实际刀路仍有 {sp['self_unsupported_regions']} 个不满足上一层承托检查的区域。"),('Auto Support','已关闭，支撑为独立 normal part 的显式 ABS 几何。'),('支撑分段',f"计划 {rem['segment_count']} 段，虚拟剪切后 {rem['fragment_count']} 个实体片段。"),('分段长度',f"平均 {rem['mean_segment_length_mm']:.3f} mm，中位数 {rem['median_segment_length_mm']:.3f} mm，最长 {rem['max_segment_length_mm']:.3f} mm。"),('剪口可达性',f"{rem['accessible_notches']}/{rem['cut_notches']} 个剪口通过 3 mm 圆柱探针接近检查。"),('环抱风险',f"{rem['wrap_traps']} 个局部截面超过 120° 门槛。"),('最大包覆角',f"{rem['maximum_wrap_angle_deg']:.3f}°；95 分位 {rem['wrap_angle_p95']:.3f}°。"),('移除方向',f"{rem['removable_segments']}/{rem['fragment_count']} 个虚拟剪切片段通过连续平移或小角度转动扫掠检查；其余片段保持为障碍物。"),('当前 gap',f"{sp['roof_gap']:.2f} mm，标记为 PHYSICAL_REMOVAL_GAP_UNCALIBRATED。"),('gap 选择依据','0.20 和 0.25 mm 均进行真实比较。覆盖下降的候选不因易拆假设而获得通过；全部未过门槛时只输出诊断最优项。'),('真实切片',f"本机 Bambu X1 Carbon / 0.2 mm / ABS / 0.10mm Standard Vessel：{b['slicer_success']}。"),('floating 警告',str(e['slice']['warning_audit']['all'])),('警告是否解释完整',f"未解决警告 {len(b['unresolved_warnings'])} 条。只有中心线、逐层承托和支撑自身都通过，才允许解释为有意保留的 Z gap。"),('预计时间与材料',f"{b['estimated_time']/3600:.2f} 小时、{b['filament']:.2f} g ABS，均为切片估算。"),('是否进入完整实物打印','当前不建议；先解决下列软件硬门槛。' if not good else '可进入人工复核后的首次工艺验证打印，尚未证明物理可靠。'),('仍需实物验证的参数','ABS 剪除及剥离力、表面质量、细分支损伤、gap、剪口/探针与真实工具的一致性、支撑稳定性。')]
    lines=['# 独立 ABS 血管芯支撑研究与运行报告','',f'当前状态：**{status}**。','']
    for i,(title,answer) in enumerate(answers,1):lines.extend([f'{i}. **{title}**：{answer}',''])
    lines+=['## 方法与边界','','冻结的血管包含完整树与四端口。中心线只提供弧长及局部方向，下表面来自真实 STL 的向下射线；采样不粗于当前 0.22 mm 线宽。支撑采用连续窄顶面和竖直幕墙，剪口主要位于幕墙体。比较 2/3 线宽、0.20/0.25 mm 间隙及局部逐层补撑。没有把支撑与血管做布尔合并。','','实际 G-code 解析包含坐标模式、挤出量、圆弧和每段线宽；只将落在显式支撑截面中的真实正挤出计为支撑。同一 XY 位置若存在更高、使间隙过小的支撑层，不允许用较低的幕墙层冒充顶面承托。逐层区域也只统计暴露的支撑顶面。误差界由 STL 浮点精度及 G-code XY 输出小数位推导，间隙另记录实际层高量化。逐层门槛不以 99% 代替 100%，也不给未经验证的长桥接豁免。','','原始中间日志曾将幕墙内部的较低层计为合格承托，因此出现过 100% 结果。该口径已由最近暴露顶面和更高层排除检查取代；最终结论以当前 coverage.json、layer_coverage.json、candidate_comparison.json 和本报告为准。历史日志保留以便追溯。','','## 未通过门槛','',*['- '+x for x in failed],'','## 运行证据','','- `tables/candidate_comparison.json`：所有真实切片候选的覆盖、间隙、线宽及支撑自身比较。','- `candidates/<ID>/coverage.json`：每一个下表面样本的通过/缺失证据。','- `candidates/<ID>/layer_coverage.json`：完整逐层需求和缺失多边形。','- `candidates/<ID>/removal.json`、`wrap_audit.json`、`travel_near_miss.json`：可达性、包覆角、扫掠和实际旅行风险。','- `bambu/profiles/`：当前用户预设、源路径、哈希、引擎版本；只在本次导出副本关闭 Auto Support。','- `validation/source_protection.json`：s1-10 全部冻结文件与核心输入的结束检查。','- `QC/16_actual_Gcode_underbody_coverage.png`：绿点具有正确间隙的实际顶面承托；红点为顶面缺失或间隙不合格。','- `QC/17_uncovered_arc_map.png`：承托检查未通过的位置（包括间隙过小），未通过时不会清空红点。','','## 入口与输出','','```bash','cd /home/lzy/projects/vascular_printing','./run_python.sh s1-11_vascular_core_support.py','```','','`--stage prepare/geometry/slice/review` 支持分步复查。输出根目录为本报告所在目录。`selected/` 在状态 NEEDS_ADJUSTMENT 时仅存诊断候选；3MF 保留血管和支撑两个独立 part，可以用于 Bambu 人工检查，不能把“切片成功”当作生产验收。没有上传或启动打印机。','']
    counts=s['orientation'];errors=json.loads((out/'reports/generation_errors.json').read_text())
    lines+=['## 本次实际运行范围','',f"共搜索 {counts['candidate_count']} 个确定性三维姿态，前 6 个姿态的 18 套基础组合中，{counts['full_candidate_count']} 套完成 STL 输出、{len(errors)} 套因 STL 精度与闭合检查未通过而保留为失败记录；没有用无效网格代替成功结果。包括加宽顶面和补撑试验在内，共 {counts['sliced_candidate_count']} 套真实血管候选完成原生切片与实际覆盖检查。",'']
    horizontal_path=out/'validation/horizontal_actual_result.json'
    if horizontal_path.exists():
        horizontal=json.loads(horizontal_path.read_text())
        lines += [f"30 mm 合成水平管的实际 G-code 中心线下侧覆盖为 {horizontal['coverage']['centerline_underbody_percent']:.2f}%，支撑自身不合格区域 {horizontal['support_self']['self_unsupported_regions']} 个。该合成模型使用管外独立见证柱，避免有意间隙导致全局空层被切片器拒绝；见证柱不计入承托覆盖，也不加入真实血管模型。其逐层区域覆盖为 {horizontal['layer_coverage']['percent']:.4f}%，因此这一项回归通过不代表所有制造门槛通过。",'']
    lines += ['拆除各分项记录见 `tables/removal_scorecard.json`。检查图中的红点同时表示顶面缺失和间隙不合格，具体原因保存在逐样本记录。','']
    test_path=out/'validation/test_results.json'
    if test_path.exists():
        tests=json.loads(test_path.read_text())
        lines += [f"程序回归：{tests['tests']} 项，失败 {tests['failures']} 项，错误 {tests['errors']} 项。测试包含失败检出能力；测试通过不能覆盖本报告列出的真实制造失败。",'']
    (out/'vascular_core_support_report.md').write_text('\n'.join(lines))
    guide=['# 血管支撑拆除参考顺序','',f'状态：**{status}**。以下是几何代理生成的参考，未经过实际剪除力或损伤验证。','']
    if failed:guide+=['当前候选存在未通过的承托或拆除检查，不建议据此直接打印或强行拆除。','']
    guide+=['1. 打印完成并冷却后，从打印板整体取下血管与支撑，不要先拉血管。','2. 对照 `QC/10_cut_notches.png` 和 `QC/19_removal_access.png`，先确认编号剪口及开放接近方向。','3. 只剪可达且不会伤及血管的幕墙剪口；探针未通过的位置不要强行操作。','4. 按下表优先从外围片段开始，逐步向中心处理；支撑顶面向下或侧向剥离。','5. 分叉附近最后处理。表中没有成功路径的片段仍需重新设计。','','|建议顺序|片段|几何移除方向|','|---|---|---|']
    for row in sorted(removal['fragments'],key=lambda r:r['suggested_order'] if r['suggested_order'] is not None else 1e9):guide.append(f"|{row['suggested_order'] or '未通过'}|{row['fragment_id']}|{np.round(row['path']['direction'],3).tolist() if row['path'] else '无已认证路径'}|")
    (out/'vascular_support_removal_guide.md').write_text('\n'.join(guide)+'\n')
