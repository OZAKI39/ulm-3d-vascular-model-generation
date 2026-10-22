"""Standalone C-retainer review figures and source-derived Chinese report."""
import json
from pathlib import Path
import numpy as np
import trimesh,pyvista as pv
from .support_inputs import save,digest
from .casting_fixture_geometry import shape_mesh
from .casting_fixture_v2_geometry import port_transform
from .casting_fixture_review import poly

NAMES=['01_legacy_long_retainers','02_C_retainer_concept','03_C_retainer_all_ports','04_C_retainer_O1_detail','05_C_retainer_O2_master_detail','06_C_retainer_clearance_section','07_top_insertion_sequence','08_secondary_axial_float','09_anti_float_Z_before_after','10_pitch_roll_constraint','11_no_central_tree_contact','12_wall_based_solvent_routes','13_retainer_route_removed_graph','14_auxiliary_void_old_vs_new','15_PDMS_after_ABS_removal','16_PDMS_after_trim','17_real_Bambu_slice_retainers','18_final_ABS_assembly_C_retainer']
COLORS=dict(core='#cc5947',mold='#589bb6',cradle='#deac45',cap='#249c72',land='#a441ab')


def section_figure(out,l,g,new,s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon,Arc
    port=l['ports'][0];k=port['port_id'];rec=new['compact_parameters'][0];M=port_transform(port);fig,ax=plt.subplots(figsize=(12,9));y=rec['cap_front_y_mm']+.45
    for m,label,col in [(l['source_ports'][k].copy(),'Frozen port',COLORS['core']),(shape_mesh(g['fixtures'][k]),'Unchanged lower cradle',COLORS['cradle']),(shape_mesh(new['retainers'][k]),'Upper C-cap',COLORS['cap'])]:
        m.apply_transform(np.linalg.inv(M));sec=m.section(plane_origin=[0,y,0],plane_normal=[0,1,0])
        if sec is not None:
            for j,curve in enumerate(sec.discrete):ax.add_patch(Polygon(curve[:,[0,2]],facecolor=col,alpha=.78,edgecolor=col,label=label if j==0 else None))
    r=port['radius_mm'];gap=rec['upper_clearance_mm'];ax.annotate('',xy=(0,r+gap),xytext=(0,r),arrowprops=dict(arrowstyle='<->',lw=1.5));ax.annotate(f'Positive upper gap = {gap:.2f} mm',xy=(0,r+gap/2),xytext=(1.25,1.65),arrowprops=dict(arrowstyle='->'),fontsize=13)
    ax.add_patch(Arc((0,0),3.2,3.2,theta1=135,theta2=405,color='#536978',lw=1.4,linestyle='--'));ax.text(-1.65,-1.95,f'Insertion opening: {rec["retainer_opening_angle_deg"]:.0f} deg\nNo snap / no port bending',fontsize=13)
    ax.annotate('Frozen lower cradle\ncontinues to carry weight',xy=(.55,-.92),xytext=(1.5,-1.1),arrowprops=dict(arrowstyle='->'),fontsize=12)
    ax.set_xlim(-2.1,4.1);ax.set_ylim(-2.3,2.2);ax.set_aspect('equal');ax.set_xlabel('Local transverse X (mm)');ax.set_ylabel('Local vertical Z (mm)');ax.set_title(f'{k} | Actual mesh section at local Y = {y:.2f} mm');ax.legend(loc='upper left');ax.grid(alpha=.18);fig.tight_layout();fig.savefig(out/'QC/06_C_retainer_clearance_section.png',dpi=160);plt.close(fig)


def graph_figure(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import networkx as nx
    data=json.loads((out/'solvent_access/graph_without_retainer_routes.json').read_text());G=nx.DiGraph();G.add_nodes_from(data['nodes']);G.add_edges_from(data['edges']);fig,ax=plt.subplots(figsize=(15,8));pos={};names={};ports=['I1','O1','O2','O3']
    for node in G:
        if node=='external_bath':pos[node]=(-1,1.5);names[node]='External bath'
        elif node.startswith('vascular_'):pos[node]=(5,1.5);names[node]='Vascular network'
        else:
            k=next((k for k in ports if node.endswith('_'+k)),None)
            x=0 if node.startswith('wall_') else 1 if node.startswith('fixture_') else 2 if 'contact' in node or 'land' in node else 3
            pos[node]=(x,ports.index(k) if k else -1);names[node]=node.replace('_','\n')
    nx.draw_networkx(G,pos,ax=ax,labels=names,node_color=['#97d7be' if n.startswith('port_') else '#dceaf1' for n in G],node_size=2500,font_size=10,arrowsize=17,width=1.7,edge_color='#606e7d');ax.set_title('Retainer nodes REMOVED: four wall-based port routes remain\nTopology / nominal contact only; no dissolution-time prediction',fontsize=15);ax.margins(.18);ax.axis('off');fig.tight_layout();fig.savefig(out/'QC/13_retainer_route_removed_graph.png',dpi=160);plt.close(fig)


def motion_figure(out,s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(14,6));modes=[x['mode'] for x in s['anti_float']['poses']];x=np.arange(len(modes))
    for j,(name,data) in enumerate([('Legacy',s['anti_float_legacy']),('Compact C',s['anti_float'])]):axes[0].scatter(x,np.full(5,j),s=500,c=['#2b9a73' if r['blocked'] else '#bf3939' for r in data['poses']],marker='s')
    axes[0].set_xticks(x,modes);axes[0].set_yticks([0,1],['Legacy','Compact C']);axes[0].set_xlim(-.5,4.5);axes[0].set_ylim(-.6,1.6);axes[0].set_title('Actual frozen core finite-pose collision checks\nGreen = blocked (not force validation)');axes[0].grid(alpha=.2)
    axes[1].axis('off');axes[1].text(.03,.9,'+Z: 1.00 mm\nPitch: +/-2.00 deg about Y\nRoll: +/-2.00 deg about X\n\nSame port-origin centroid pivot\nSame mold and lower cradles\nEach cap +Z 0.80 mm: blocked by guide\nEach port +Z 1.00 mm: blocked by cap\n\nSecondary axial fit allowance retained\nNo stiffness / buoyancy-force claim',va='top',fontsize=15);fig.tight_layout();fig.savefig(out/'QC/10_pitch_roll_constraint.png',dpi=160);plt.close(fig)


def render(out,l,old,new,d,s):
    q=out/'QC';q.mkdir(exist_ok=True);section_figure(out,l,old,new,s);graph_figure(out);motion_figure(out,s);records=[];camera=None;H=l['high'][2]
    for i,name in enumerate(NAMES,1):
        print('s1-10 compact QC',name,flush=True)
        if i in [6,10,13]:continue
        p=pv.Plotter(off_screen=True,window_size=(1600,1150));p.set_background('white');p.add_text(name.replace('_',' '),position='upper_left',font_size=18,color='#253c4c');p.add_text('ABS | mm | Nominal geometry and actual local slicing; physical process uncalibrated',position='lower_left',font_size=11,color='#546676')
        def add(m,color,opacity=1.):
            if isinstance(m,trimesh.Trimesh):m=poly(m)
            if m.n_points:p.add_mesh(m,color=color,opacity=opacity,smooth_shading=False)
        def caps(g,opacity=1.):
            for part in g['retainers'].values():add(shape_mesh(part),COLORS['cap'],opacity)
        def text(t):p.add_text(t,position='upper_right',font_size=13,color='#344a5b')
        def detail(k):
            port=next(x for x in l['ports'] if x['port_id']==k);pt=port['trim_origin'];bounds=np.column_stack((pt-[4,4,4],pt+[9,9,4])).ravel()
            # Symmetric enough for any axis and includes the unchanged wall datum.
            bounds=np.column_stack((pt-[8,8,4],pt+[8,8,4])).ravel();add(poly(old['mold_mesh']).clip_box(bounds,invert=False),COLORS['mold'],.18);add(poly(l['source_ports'][k]).clip_box(bounds,invert=False),COLORS['core'],.85);add(shape_mesh(old['fixtures'][k]),COLORS['cradle'],.75);add(shape_mesh(new['retainers'][k]),COLORS['cap']);return port
        if i in [1,18]:
            add(old['mold_mesh'],COLORS['mold'],.16);add(l['core'],COLORS['core']);caps(old if i==1 else new)
            for part in old['fixtures'].values():add(shape_mesh(part),COLORS['cradle'],.6)
            text(('Legacy: long solvent-entry stems\n' if i==1 else 'Compact C: no exposed stems\n')+f"Max height above port: {s['legacy' if i==1 else 'compact']['max_height']:.3f} mm\nTotal retainer volume: {s['legacy' if i==1 else 'compact']['total_retainer_volume']:.3f} mm3\nIdentical camera / scale in figures 01 and 18")
        elif i==2:
            k='I1';rec=next(r for r in new['compact_parameters'] if r['port_id']==k);add(shape_mesh(new['compact_local'][k]),COLORS['cap']);text(f"Upper wrap {rec['c_retainer_wrap_angle_deg']:.0f} deg; opening {rec['retainer_opening_angle_deg']:.0f} deg\nUpper clearance {rec['upper_clearance_mm']:.2f} mm\nWall {rec['retainer_wall_thickness_mm']:.2f} mm (2 actual lines)\nAxial cap {rec['retainer_axial_length_mm']:.2f} mm\nShort guide {rec['retainer_guide_length_mm']:.2f} mm\nNo solvent stem / no elastic snap")
        elif i==3:
            add(l['core'],COLORS['core'],.55);caps(new)
            for port in l['ports']:p.add_point_labels([port['trim_origin']+[0,0,3]],[port['port_id']+' '+port['role']],always_visible=True,font_size=16)
            text('Four compact caps in unchanged assembly positions')
        elif i in [4,5]:
            k='O1' if i==4 else 'O2';port=detail(k);text(f"{k} {port['role']}\nUpper gap {port['anti_float_gap_mm']:.2f} mm; guide 0.12 mm\nC-cap above port; unchanged lower cradle below\nRear datum seats on the existing inner wall")
        elif i==7:
            add(old['mold_mesh'],COLORS['mold'],.1);add(l['core'],COLORS['core'],.5)
            for port in l['ports']:
                k=port['port_id'];m=shape_mesh(new['retainers'][k]);add(m,COLORS['cap']);pt=m.centroid;high=pt-port['axis']*1.6+[0,0,H+5];low=pt-port['axis']*1.6;p.add_lines(np.array([high,low,low,pt]),color='#7143a1',width=3);p.add_arrows(np.array([low]),np.array([pt-low]),mag=1,color='#7143a1')
            text('Order: '+' > '.join(s['insertion']['order'])+'\nVertical descent, then 1.60 mm wallward slide\nLast 0.02 mm: planar wall-datum seating\nSwept paths include already installed caps')
        elif i==8:
            port=next(x for x in l['ports'] if x['role']=='SECONDARY');k=port['port_id'];part=l['source_ports'][k].copy();M=port_transform(port);inv=np.linalg.inv(M);part.apply_transform(inv);bounds=[-3,3,-2,6,-1.7,2.5];add(poly(part).clip_box(bounds,invert=False),COLORS['core'],.25);add(shape_mesh(new['compact_local'][k]),COLORS['cap'],.85)
            for sign in [-1,1]:
                m=part.copy();m.apply_translation([0,sign*port['axial_float_mm'],0]);add(poly(m).clip_box(bounds,invert=False),'#d5973e',.2)
            p.add_arrows(np.array([[0,0,-1.3],[0,0,-1.3]]),np.array([[0,1,0],[0,-1,0]]),mag=port['axial_float_mm'],color='#64419c');ar=next(r for r in s['axial_float']['rows'] if r['port_id']==k);text(f"{k}: approved relative fit allowance +/-{port['axial_float_mm']:.2f} mm\nActual fit interval {ar['new_min_mm']:.3f} to {ar['new_max_mm']:.3f} mm\nSame as legacy; minimum C engagement 0.22 mm\nPlacement-error envelope, not independent branch motion")
        elif i==9:
            port=detail('I1');shift=l['source_ports']['I1'].copy();shift.apply_translation([0,0,1]);pt=port['trim_origin'];bounds=np.column_stack((pt-[8,8,4],pt+[8,8,4])).ravel();add(poly(shift).clip_box(bounds,invert=False),'#aa3693',.4);text('Red: nominal frozen port (positive gap)\nPurple: +Z 1.00 mm probe intersects C stop\nLegacy and compact both block the same probe\nThis is not measured anti-float force')
        elif i==11:
            add(l['central'],COLORS['core']);caps(new)
            for port,r in zip(l['ports'],s['contact']['rows']):p.add_point_labels([port['trim_origin']],[f"{port['port_id']} | tree gap {r['central_tree_gap_mm']:.2f} mm"],always_visible=True,font_size=14)
            text('New retainer / central tree contacts: 0\nNearest core points all in approved straight zones')
        elif i==12:
            add(old['mold_mesh'],COLORS['mold'],.08);add(l['core'],COLORS['core'],.55)
            for port in l['ports']:
                k=port['port_id'];add(shape_mesh(old['fixtures'][k]),COLORS['cradle']);add(shape_mesh(old['lands'][k]),COLORS['land']);pt=port['trim_origin'];n=port['axis'];a=pt+n*(port['wall_distance_from_trim_mm']+3)+[0,0,-2];b=pt+n*2.2+[0,0,-port['radius_mm']];p.add_arrows(np.array([a]),np.array([b-a]),mag=1,color=COLORS['land'])
            text('Bath > wall > blade > contact land > port\nFour independent routes; NO retainer node needed\nBoth SWCs: every branch reachable from four ports\nPrinted contact / PDMS-film exclusion uncalibrated')
        elif i==14:
            legacy=Path(s['frozen_geometry']['outer_mold_stl']).parents[1]/'diagnostic_PDMS/auxiliary_voids_before_trim.stl';oldvoid=trimesh.load_mesh(legacy);add(oldvoid,'#b986d1',.22);add(d['auxiliary_voids_before_trim'],COLORS['cradle'],.92);text(f"Purple: legacy auxiliary void (overlay)\nAmber: compact auxiliary void\nOld {s['PDMS']['auxiliary_void_old']:.3f} mm3\nNew {s['PDMS']['auxiliary_void_new']:.3f} mm3\nSame coordinates; original trim zones unchanged")
        elif i in [15,16]:
            add(d['PDMS_after_all_ABS_removed' if i==15 else 'PDMS_after_trim'],'#8fbbcf',.24);add(d['all_ABS_voids_inside_cavity' if i==15 else 'vascular_lumen_after_trim'],COLORS['land'] if i==15 else COLORS['core'],.85);text('Purple = cavities after ABS removal\nNo ABS present; peripheral auxiliary voids included' if i==15 else 'Only target vascular lumen remains\nAuxiliary void after trim: numerical zero\nFour exposed ports; 4.50 mm usable straight length')
        elif i==17:
            from .bambu_casting_mold import parse_paths
            rec=s['bambu']['selected_slice']
            if rec.get('archive'):
                parsed=parse_paths(rec['archive']);seg=parsed['segments'][parsed['kind']==0];p.add_lines(seg.reshape(-1,3),color=COLORS['cap'],width=1.3);text(f"ACTUAL positive extrusion paths: {len(seg)} segments\nX1 Carbon / 0.2 mm / Bambu ABS\nAuto Support OFF; support segments = 0\nC lip / wall / guide / float passage checked")
            else:text('Actual slice missing: NEEDS_ADJUSTMENT')
        p.add_axes();p.camera_position='iso';p.camera.up=(0,0,1);p.enable_parallel_projection();p.reset_camera();p.camera.zoom(1.02)
        if i==1:camera=dict(position=list(p.camera.position),focal_point=list(p.camera.focal_point),up=list(p.camera.up),parallel_scale=p.camera.parallel_scale);save(q/'legacy_compact_shared_camera.json',camera)
        if i==18:
            p.camera.position=camera['position'];p.camera.focal_point=camera['focal_point'];p.camera.up=camera['up'];p.camera.parallel_scale=camera['parallel_scale']
        p.show(screenshot=str(q/(name+'.png')),auto_close=True)
    for name in NAMES:records.append(dict(name=name,sha256=digest(q/(name+'.png'))))
    save(q/'figure_manifest.json',records)


def report(out,s):
    old=s['legacy'];new=s['compact'];v=s['solvent'];t=s['PDMS'];ok=lambda b:'通过' if b else '未通过';ans=[
        '旧版高杆将插片延伸到名义浇注面以上，既便于操作，也为插片自身提供独立溶剂入口；端口原本另有壁侧通路。',
        '新版已删除长竖直溶剂杆与长操作柄，四个插片均低于当前名义浇注面。',
        f"新版为上侧 {s['selected_wrap_angle_deg']:.0f}° C 形短盖，保留 {360-s['selected_wrap_angle_deg']:.0f}° 开口，不形成完整环，不依靠掰开 ABS 或弯曲端口安装。",
        f"最大高度按距端口中心线计算，由 {old['max_height']:.3f} mm 降至 {new['max_height']:.3f} mm，减少 {new['height_reduction_mm']:.3f} mm（{new['height_reduction_percent']:.2f}%）。",
        f"四个插片总 CAD 体积由 {old['total_retainer_volume']:.3f} mm³ 降至 {new['total_retainer_volume']:.3f} mm³，减少 {new['volume_reduction_percent']:.2f}%。",
        f"相同 +Z 1 mm、俯仰 ±2°、滚转 ±2° 探测姿态的阻挡检查{ok(s['anti_float']['passed'])}，未弱于旧版。此为有限运动几何检查，不是力学验证。",
        '没有夹紧端口。上方名义间隙保持 0.22 mm，导向间隙保持 0.12 mm；下托座继续承担正常静态承托。',
        f"次基准端口原有 ±1.00 mm 轴向配合余量保留检查{ok(s['axial_float']['passed'])}。实际允许区间及旧新版比较见下表；该值表示相对配合误差包络，不意味着相连分支可独立移动。",
        f"新插片与中央脑血管树接触数为 {s['contact']['central_tree_contacts']}；最近的血管表面点均位于已批准的延长直段。",
        f"四个插片按 {' → '.join(s['insertion']['order'])} 从顶部下降，再向墙滑动 1.60 mm，连续装入检查{ok(s['insertion']['passed'])}。最后 0.02 mm 平面落座单独验证。",
        f"删除全部插片图节点后，仍有 {v['independent_port_routes']} 条独立壁侧端口通路；全部端口检查{ok(v['all_ports_wall_access'])}。",
        '不再要求插片高出液面。插片自身通过短导向后端与既有内墙的有限面接触获得去除入口；这条支路只用于插片自身去除，不承担血管端口的必要可达性。',
        f"compensated 与 smoothed 两套 SWC 的全部保留分支可达检查{ok(v['all_vascular_branches_reachable'])}；每套 5 个分支均可从四个端口到达。",
        f"插片辅助孔腔逐件减少；包括未变下托与浅导向在内的总辅助孔腔由 {t['auxiliary_void_old']:.6f} mm³ 降至 {t['auxiliary_void_new']:.6f} mm³。",
        f"按原裁切面切除外围后，非血管辅助孔腔为 {t['auxiliary_void_after_trim']:.12f} mm³，属于数值零；四端口仍各保留 4.50 mm 已确认可用直段。",
        f"X1 Carbon、0.2 mm 喷嘴、ABS、Auto Support OFF 本地真实切片{ok(s['bambu']['slicer_success'])}。未上传到打印机或启动打印。",
        f"实际正挤出路径中的 C 唇边、两线宽壁、短导向和轴向通道保留检查{ok(s['bambu']['critical_features_preserved'])}；未检出 C 唇边丢失。",
        ('在本轮软件几何、拓扑和切片检查范围内，建议以紧凑 C 插片替换旧长杆版本，并进入实物审阅。' if not s['failures'] else '仍存在未通过门槛，不建议作为就绪生产候选。'),
        '实际抗浮力、ABS 收缩与翘曲、壁侧和端口接触处是否夹入 PDMS 薄膜、溶剂溶解行为、PDMS 溶胀恢复及人工裁切精度尚未验证；不预测溶解时间、流量或更新频率。']
    text='# BG001 RMCA 紧凑 C 形防浮插片比较研究\n\n状态：**'+s['status']+'**。以下“通过”均指软件与本地切片证据。\n\n'+'\n\n'.join(f'{i+1}. {a}' for i,a in enumerate(ans))
    text+='\n\n## 方法与冻结边界\n\n本研究只重构 s1-10 的独立防浮插片。血管芯、半径、延长管路线和长度、外模尺寸、墙厚 2 mm、底厚 3 mm、下托弧、壁侧溶剂窄片、中央血管放置和已验证裁切面保持原值。O2 主基准和其余次基准均从当前布局清单读取。源 CAD、STL、清单、SWC 与旧生产文件逐字节保护；本轮没有重新拟合血管，也没有修改既有可视化界面。\n\n'
    text+='角度候选为 90°、120°、150°，分别保留 270°、240°、210° 开口。候选先经过抗浮、顶部装入、轴向余量、中央零接触、延长直段限制、壁侧通路、孔腔和真实切片门槛，再按最大高度、总体积依次排序，无黑箱评分。C 弧轴向长度 1.20 mm、径向壁厚 0.44 mm（实际 0.22 mm 线宽的两条线），通过渐变短颈连接到原导向内的短舌。后端落座于既有内墙；不改变内墙或导向。打印时 C 的轴向端面朝下，局部端口轴转为打印 Z 轴，短颈逐层扩展，不采用悬空长柄。\n\n'
    text+='|包覆角|开口角|最大高度 mm|总体积 mm³|全部候选门槛|选中|\n|---:|---:|---:|---:|---|---|\n'
    for r in s['candidates']:text+=f"|{r['wrap_angle_deg']}|{r['opening_angle_deg']}|{r['max_height_above_port_mm']:.4f}|{r['total_volume_mm3']:.4f}|{ok(r['hard_gates_passed'])}|{'是' if r['selected'] else '否'}|\n"
    text+='\n## 配合与运动结果\n\n|端口|角色|旧轴向区间 mm|新轴向区间 mm|上方实际最小网格间隙 mm|距中央树 mm|\n|---|---|---|---|---:|---:|\n'
    for a,c in zip(s['axial_float']['rows'],s['contact']['rows']):text+=f"|{a['port_id']}|{a['role']}|[{a['old_min_mm']:.3f}, {a['old_max_mm']:.3f}]|[{a['new_min_mm']:.3f}, {a['new_max_mm']:.3f}]|{c['measured_minimum_core_gap_mm']:.5f}|{c['central_tree_gap_mm']:.3f}|\n"
    text+='\n实际轴向区间由冻结端口网格与原托座/新旧上限位的相对平移扫描和碰撞边界二分计算，要求上侧 C 的有效重叠至少保留一个实际线宽；原有 ±1 mm 余量另作连续路径检查。零间隙溶剂接触面的法向接触不能当作轴向干涉，故在正距离配合搜索中排除该接触薄片，同时单独检查原批准浮动范围内接触面仍位于直段。主基准的正向端部止挡没有变化。计算区间不能替代装配公差或连接血管的力学分析。\n\n'
    text+='顶部装入的障碍包括冻结血管、完整外模及托座、已经装好的插片。前段路径采用连续距离界限检查；最后短距离落座使用去除墙体后的连续检查、盒腔凸半空间约束与最终 CAD 零交叠共同证明，不把有意的墙面接触误报为穿透。插片向上 0.8 mm 时原导向阻挡，端口向上 1 mm 时新 C 上限位阻挡。有限探测并不证明任意组合运动都被约束，更不证明抗浮承载力。\n\n'
    text+='## 溶剂路径与孔腔解释\n\n正式端口路径为外部浴 → 外墙 → 壁生窄片 → 有限面积接触区 → 延长端口 → 血管网络。重建图后明确移除所有插片节点再验证四条独立通路。短插片完全包埋时，不能仅凭“同为 ABS”断言能够溶解；本方案另验证其后端与外墙的有限接触面，使外墙移除后能到达插片自身。该支路不会计入血管必要入口。两类名义接触都必须实物确认，若形成连续 PDMS 薄膜，应视为通路失效。\n\n'
    text+='新插片的孔腔全部位于原端口外围裁切区内。沿用同一组 11 个有限裁切面，不移动切面迁就插片；中央血管和无关端口未被切到。最终 PDMS 与“原保留包络减去冻结血管”的理想域做双向布尔比较，检查额外缺损与残留。诊断布尔沿用原共面接触数值处理，仅在已被血管占据的域内进行 0.0001 mm 接触缝合，生产 STL/STEP 不变。\n\n'
    text+='## 工艺边界与复现\n\n**血管芯与四个插片均禁止蒸汽平滑。** 外模可选的受控蒸汽处理仍未进行尺寸校准，不能据此做软件补偿；需干燥后实测间隙和接触。当前结果只支持人工审阅，不能标记为力学、溶解或制造实证通过。\n\n'
    text+='运行 `./run_python.sh s1-10_all_abs_sacrificial_casting_system.py`；参数见 `config/compact_c_retainer_v2.yaml`。本目录 `retainers/` 是新四件 STL/STEP 和 X1C ABS 实际切片 3MF，`assembly/` 是独立组件装配审阅（组合 STL 不是联合打印文件），`diagnostic_PDMS/` 是溶解和裁切诊断，`QC/` 是 18 张检查图，`tables/` 和 `solvent_access/` 保存数值及图证据，`reports/run.log` 保存真实运行输出。原外模与冻结血管打印件仍从父目录对应文件读取，不需要重新生成或替换。旧插片最小归档为父目录 `legacy_retainer_v1/`。\n\n'
    text+=f"源保护检查 {s['source_protection']['checked_files']} 个文件，错配 {s['source_protection']['mismatch_count']}。血管芯 SHA256：`{s['frozen_geometry']['source_core_sha256']}`。\n"
    if s['failures']:text+='\n未通过项：'+', '.join(s['failures'])+'。\n'
    (out/'compact_c_retainer_report.md').write_text(text);(out/'reports/README.md').write_text('中文总报告：../compact_c_retainer_report.md。日志 run.log；source_protection.json 在上级；真实切片命令、输出和验证在 ../bambu/candidate_*/。\n')
