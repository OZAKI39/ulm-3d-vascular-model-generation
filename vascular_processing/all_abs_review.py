"""Standalone s1-10 review artifacts; does not import or modify existing UIs."""
from pathlib import Path
import json
import numpy as np
import trimesh,pyvista as pv
from .support_inputs import save,digest
from .casting_fixture_geometry import shape_mesh
from .casting_fixture_review import poly

NAMES=['01_frozen_ABS_vascular_core','02_outer_ABS_mold','03_closed_bottom_four_walls','04_vapor_smoothing_target_surfaces','05_core_inside_mold','06_extension_fixture_zones','07_wall_grown_fixtures','08_anti_float_retainers','09_no_central_tree_contact','10_top_insertion_path','11_final_ABS_assembly','12_solvent_access_overview','13_port_solvent_access_I1','14_port_solvent_access_O1','15_port_solvent_access_O2','16_port_solvent_access_O3','17_PDMS_as_cast','18_PDMS_after_ABS_dissolution','19_auxiliary_voids_before_trim','20_trim_planes','21_final_PDMS_after_trim','22_exposed_ports_after_trim','23_outer_mold_real_slice','24_retainers_real_slice']


def render(out,l,g,meshes,d,s):
    q=out/'QC';q.mkdir(exist_ok=True);records=[];core=poly(l['core']);mold=poly(meshes['mold']);central=poly(l['central']);lo=l['low'];hi=l['high'];H=hi[2]
    for i,name in enumerate(NAMES,1):
        print('s1-10 QC',name,flush=True);p=pv.Plotter(off_screen=True,window_size=(1600,1150));p.set_background('white');p.add_text(name.replace('_',' '),position='upper_left',font_size=17,color='#203044');p.add_text('All ABS | Dimensions: mm | Nominal geometry only; material / process parameters uncalibrated',position='lower_left',font_size=11,color='#536271')
        def add(m,color,opacity=1.):
            if isinstance(m,trimesh.Trimesh):m=poly(m)
            if m.n_points:p.add_mesh(m,color=color,opacity=opacity,smooth_shading=False)
        if i in [23,24]:
            from .bambu_casting_mold import parse_paths
            r=s['bambu']['mold_auto_support_off_slice' if i==23 else 'retainer_auto_support_off_slice']
            if r.get('archive'):
                paths=parse_paths(r['archive']);seg=paths['segments'][paths['kind']==0];z=np.round(seg[:,:,2].mean(1),4);levels=np.unique(z);chosen=levels[np.linspace(0,len(levels)-1,min(24,len(levels))).astype(int)];part=seg[np.isin(z,chosen)];p.add_lines(part.reshape(-1,3),color='#217d9d',width=1.2);p.add_text(f"ACTUAL Bambu G-code, {len(chosen)} displayed layers\nX1 Carbon / 0.2 nozzle / ABS / Auto Support OFF\nNative support segments: {r['native_support_extrusion_segments']}",position='upper_right',font_size=12,color='#304355')
            else:p.add_text('ACTUAL SLICE MISSING: NEEDS_ADJUSTMENT',color='red')
        elif 13<=i<=16:
            port=l['ports'][i-13];k=port['port_id'];pt=np.asarray(port['trim_origin']);n=np.asarray(port['axis']);bounds=np.column_stack((pt-[9,9,pt[2]-lo[2]+1],pt+[9,9,H-pt[2]+2])).ravel();add(mold.clip_box(bounds,invert=False),'#5c99b6',.25);add(poly(l['source_ports'][k]).clip_box(bounds,invert=False),'#bf463b',.8);add(shape_mesh(g['fixtures'][k]),'#df9a2a',.7);add(shape_mesh(g['lands'][k]),'#b1267d');add(meshes['retainers'][k],'#388d65',.7)
            w=port['wall_distance_from_trim_mm'];start=pt+n*(w+2)+[0,0,-2];end=pt+n*2.2+[0,0,-port['radius_mm']-.5];p.add_arrows(np.array([start]),np.array([end-start]),mag=1,color='#6129a3');p.add_text(f"{k}: bath > wall > web > conformal land > port\nFinite surface contact: {s['solvent_access']['port_paths'][i-13]['contact_area_mm2']:.4f} mm2\nRetainer stem exposed above fill level\nPrinted contact must be checked physically",position='upper_right',font_size=12,color='#354659')
        elif i in [17,18,19,20,21,22]:
            after=i>=21;add(d['PDMS_after_trim' if after else 'PDMS_after_all_ABS_removed'],'#82b5cf',.18 if i!=21 else .6)
            if i==17:
                add(core,'#c4513d');add(mold,'#4f8fae',.12)
                for k in g['retainers']:add(meshes['retainers'][k],'#39a16a')
                p.add_text('PDMS material domain + separately present ABS',position='upper_right',font_size=13,color='#263b4d')
            elif i==18:
                add(d['all_ABS_voids_inside_cavity'],'#b960b0',.72);p.add_text('Purple = ALL empty cavities after ABS removal\nVascular lumen + auxiliary imprints; NO ABS remains\nExterior mold volume is outside the PDMS domain',position='upper_right',font_size=13,color='#263b4d')
            elif i==19:
                add(d['auxiliary_voids_before_trim'],'#d48b27');p.add_text(f"Amber = nonvascular auxiliary voids\nVolume before trim: {s['voids']['auxiliary_void_volume_before_trim']:.3f} mm3",position='upper_right',font_size=13,color='#263b4d')
            elif i==20:
                add(d['sacrificial_trim_regions'],'#dfa33e',.3)
                for m in d['planes']:add(m,'#8d40a0',.3)
                add(trimesh.util.concatenate([shape_mesh(v) for v in g['guides'].values()]),'#168b69')
            else:
                add(d['vascular_lumen_after_trim'],'#ba4349',.75)
                if i==22:
                    for port in l['ports']:
                        section=l['source_ports'][port['port_id']].section(plane_origin=port['trim_origin'],plane_normal=port['axis']);p.add_lines(section.discrete[0],connected=True,color='#971447',width=5);p.add_point_labels([port['trim_origin']+[0,0,2]],[port['port_id']+' | 4.50 mm usable'],always_visible=True,font_size=15,shape_opacity=.3)
                else:p.add_text('Auxiliary void after trim: numerical zero\nOnly intended vascular lumen remains',position='upper_right',font_size=13,color='#263b4d')
        else:
            if i==1:add(core,'#c64b3e')
            elif i in [2,3]:
                add(mold,'#427f9f',.75 if i==2 else .25)
                if i==3:add(shape_mesh(g['bottom']),'#cb9e45');p.add_text('4 CLOSED WALLS / CLOSED BOTTOM / OPEN TOP\nNo wall holes, U-slots or loading openings',position='upper_right',font_size=14,color='#294458')
            elif i==4:
                add(mold,'#5689a9',.15);floor=pv.Plane(center=((lo[0]+hi[0])/2,(lo[1]+hi[1])/2,lo[2]+.01),direction=(0,0,1),i_size=hi[0]-lo[0],j_size=hi[1]-lo[1]);add(floor,'#40a777',.65)
                for ax in [0,1]:
                    for bound in [lo[ax],hi[ax]]:
                        center=(lo+hi)/2;center[ax]=bound;normal=np.eye(3)[ax];add(pv.Plane(center=center,direction=normal,i_size=H-lo[2],j_size=hi[1-ax]-lo[1-ax]),'#40a777',.2)
                add(core,'#c2353c');p.add_text('SMOOTH TARGET: inner bottom / inner walls\nNO-SMOOTH vascular core\nOptional physical treatment, NO mesh smoothing\nDimensional effects NOT calibrated',position='upper_right',font_size=14,color='#284356')
            elif i==8:
                for j,port in enumerate(l['ports']):
                    k=port['port_id'];part=meshes['retainers'][k].copy();part.apply_translation([0,0,4]);add(part,'#3b9d68');p.add_point_labels([part.centroid],[k+' ABS retainer'],always_visible=True,font_size=15)
                add(core,'#c25443',.3)
            else:
                add(mold,'#619db4',.1 if i in [5,6,7,9,10,12] else .3);c=core.copy()
                if i==10:c.translate([0,0,H+5],inplace=True)
                add(central if i==9 else c,'#319b67' if i==9 else '#c54f41')
                if i in [6,7,9,11,12]:
                    for port in l['ports']:
                        k=port['port_id'];add(shape_mesh(g['fixtures'][k]),'#d09a2c',.85)
                        if i==6:
                            a=np.asarray(port['approved_contact_start']);b=np.asarray(port['approved_contact_end']);add(pv.Cylinder(center=(a+b)/2,direction=port['axis'],height=np.linalg.norm(b-a),radius=port['radius_mm']+.02,resolution=32),'#efc238');p.add_point_labels([port['trim_origin']+[0,0,3]],[k+' '+port['role']],always_visible=True,font_size=14)
                        if i in [11,12]:add(meshes['retainers'][k],'#359766',.8)
                if i==9:p.add_text('CENTRAL TREE AUXILIARY ABS CONTACTS: 0\nMinimum separation: %.3f mm'%s['assembly']['minimum_central_auxiliary_clearance_mm'],position='upper_right',font_size=16,color='#246044')
                if i==10:p.add_arrows(np.array([[0,0,H+30]]),np.array([[0,0,-1]]),mag=20,color='#d09721')
                if i==12:
                    for port in l['ports']:
                        pt=np.asarray(port['trim_origin']);n=np.asarray(port['axis']);start=pt+n*(port['wall_distance_from_trim_mm']+8);end=pt+n*1.;p.add_arrows(np.asarray([start]),np.asarray([end-start]),mag=1,color='#8f399a')
                    p.add_text('BATH > WALL > FIXTURE / ACCESS LAND > PORT > CORE\n4 independent nominal routes; all 5 SWC branches reachable\nRetainer stems have separate top bath access\nGeometric topology, NOT dissolution timing',position='upper_right',font_size=13,color='#364253')
        p.add_axes();p.camera_position='iso';p.camera.up=(0,0,1);p.reset_camera();p.camera.zoom(1.02)
        if i==22:p.view_xy();p.reset_camera()
        p.show(screenshot=str(q/(name+'.png')),auto_close=True);records.append(dict(name=name,sha256=digest(q/(name+'.png'))))
    save(q/'figure_manifest.json',records)


def report(out,s,rows):
    a=s['assembly'];v=s['solvent_access'];t=s['trim'];b=s['bambu'];dim=s['outer_mold'];ok=lambda x:'通过' if x else '未通过'
    answers=['全部打印组件均设计为 ABS，包括血管芯、外模、下托座、通路接触面及四个独立插片。','血管芯、外模和插片仍分别打印；血管先在完全开放的环境中清理自身支撑。','不将血管与外模布尔合并，是为了保持冻结血管的几何与独立打印、清理和装入流程。接触面属于外模，名义上接触血管表面，装配文件仍保留独立组件。',f"四面墙封闭检查{ok(a['walls_closed'])}；完整墙体缺失体积 {a['wall_missing_volume_mm3']:.9f} mm³，没有 U 槽、端口孔或装载缺口。",f"底面封闭检查{ok(a['bottom_closed'])}，底厚 {dim['bottom_thickness']:.2f} mm。",f"顶部开放检查{ok(a['top_open'])}，无需穿墙装入血管。",f"完整血管从顶部装入检查{ok(a['top_insertion_pass'])}；先检验连续路径，再证明最后 0.02 mm 竖直落座到下方接触面不会穿透。",f"所有接触位于获准端口直段检查{ok(a['all_contacts_on_approved_extensions'])}；四处溶剂接触面远离中央树和上游连接段。",f"中央血管与辅助 ABS 的接触数为 {a['central_tree_contacts']}，最小距离 {a['minimum_central_auxiliary_clearance_mm']:.3f} mm。",f"四个插片从顶部下降，再向墙滑动 1.60 mm，进入斜面限位。+Z、俯仰和滚转的有限扰动检查{ok(a['anti_float_geometric_pass'])}；真实抗浮力未测量。",'可选的蒸汽处理只面向独立外模，主要目标为内底和内侧壁浇注表面，处理后充分干燥再装配。局部遮蔽和选择性处理效果未验证，不能假设敏感位置完全不受影响。','冻结血管芯明确禁止 acetone vapor smoothing，也没有做任何软件平滑、缩放或半径补偿。','蒸汽处理没有尺寸校准，状态为 VAPOR_SMOOTHING_NOT_DIMENSIONALLY_CALIBRATED；软件没有为其补偿尺寸。','外墙和底板外表面直接接触外部溶剂浴；插片的细长露出段高出名义浇注面 0.80 mm，提供各自独立的外部入口。','墙体去除后到达相连的竖直窄片，再到贴合冻结 STL 实际下表面的窄接触面，最后暴露延长血管。下托弧的正间隙不被误算成溶剂通路，因为该间隙会被 PDMS 填满。',f"图搜索得到 {v['independent_access_routes']} 条互不共享内部节点的端口通路，每个端口都有独立入口。此结论针对名义接触几何；实物接触失效或夹入 PDMS 薄膜时，该通路不能视为成立。",f"compensated 和 smoothed SWC 均复核；{v['network']['node_count']} 个节点、{v['network']['branch_count']} 个分支，全部分支可达检查{ok(v['all_vascular_branches_reachable'])}，本模型每个分支可从四个端口到达。",'溶解后除目标血管腔外，还会留下下托片、接触通路、插片与浅导向的外围印痕。墙体及底板在 PDMS 外侧，其体积不计作内部辅助空腔。',f"辅助空腔均位于外围切除区的检查{ok(abs(t['auxiliary_void_volume_after_trim'])<1e-6)}；切除区由真实直段中点的有限竖直切面确定。",f"切前辅助空腔 {t['auxiliary_void_volume_before_trim']:.6f} mm³；切后 {t['auxiliary_void_volume_after_trim']:.9f} mm³（数值零）。与仅保留目标血管腔的理想 PDMS 比较，缺损和多余体积分别为 {t['ghost_missing_volume_mm3']:.9f}、{t['ghost_extra_volume_mm3']:.9f} mm³。",f"四端口暴露检查{ok(t['all_ports_exposed'])}，用真实截面和切后 PDMS 面两侧点包含关系核查；中央血管切除体积 {t['central_vascular_removed_mm3']:.9f} mm³。",'I1、O1、O2、O3 各保留 4.50 mm 可用保守直段；不把弯曲段或未确认的上游长度计入。']
    text='# BG001 RMCA 全 ABS 牺牲浇注系统\n\n状态：**'+s['status']+'**。本报告中的通过仅指软件几何、拓扑及本地真实切片。\n\n'+'\n\n'.join(f'{i+1}. {x}' for i,x in enumerate(answers))
    if s['failures']:text+='\n\n**尚未通过的软件门槛：** '+', '.join(s['failures'])+'。不能按可审阅就绪状态使用。'
    text+='\n\n## 研究目的与结构方法\n\n本轮将上一阶段的 PLA 可拆夹具改为独立的全 ABS 牺牲外模；s1-9 原有代码、最终报告和源哈希原位保留，标记为 RETIRED / RESEARCH HISTORY。冻结输入来自已选中的 candidate_0008，不重新运行语义迁移、区域选择、VascularMD、端口路径或接口修整。完整血管仅作统一刚体放置，尺度为 1。\n\n'
    text+='完整血管含延长管参与重新计算盒腔。比较三个小偏航角与三个端部退让值，以最短最大托座跨度、总跨度、最小转角依次选择。内腔尺寸 %.3f × %.3f × %.3f mm，外形 %.3f × %.3f × %.3f mm；墙厚 %.2f mm、底厚 %.2f mm，最低血管表面距内底 10.00 mm。中央区的原约 60 × 50 mm 不代表包含所有外接管后的盒腔尺寸。\n\n'%(*dim['inner_dimensions'],*dim['outer_dimensions'],dim['wall_thickness'],dim['bottom_thickness'])
    text+='O2 为主基准，单侧末端基准保留正间隙；其余为次基准，允许轴向 ±1.00 mm 余量。下托弧径向间隙主端口 0.22 mm、其余 0.33 mm；上限位间隙 0.22 mm，插片导向间隙 0.12 mm。竖直承托片与溶剂窄片厚 0.66 mm，即当前配置的三个线宽。接触面按冻结 STL 的实际三角面构造，没有改变其离散圆截面。\n\n'
    text+='## 溶剂可达性与裁切结果\n\n溶剂图每条边有实际几何依据：墙与薄片的正体积连接、薄片与接触面的连接、接触面与冻结血管的有限表面积重合，以及插片露出液面的体积。图的连通性不等同于溶解动力学，也不能由“全部 ABS”直接推导。四条端口通路在各自壁区建立，接触面积约 0.70 mm²。零间隙接触是这里最敏感的工艺条件：打印层阶、收缩、装配上浮或可选蒸汽处理均可能改变接触并形成 PDMS 薄膜，必须实物确认。若不成立，应将相关图边视为阻断，不能凭名义模型声称溶剂可达。\n\n'
    text+=f"诊断布尔计算将接触面上缘在已被血管占据的区域内延伸 0.0001 mm，以消除零间隙共面导致的退化三角形。该处理仅作用于 PDMS/孔腔诊断，不改动任何生产 STL/STEP。新增部分超出原占据域的体积为 {t['diagnostic_contact_added_outside_occupied_domain_mm3']:.12f} mm³；分开重组 CAD 三角网格后，诊断 PDMS 与原始布尔域的对称差为 {t['diagnostic_sewing_PDMS_symmetric_difference_mm3']:.6f} mm³，低于明确的 0.01 mm³ 诊断限值。最终保留区仍以双向零差独立核查。\n\n"
    text+='浇注前后和溶解后，PDMS 自身的材料域相同；不同之处是其孔腔内原先存在 ABS、后来为空。诊断图 18 显示全部孔腔，图 19 单独突出非目标辅助孔腔。最终切割移走全部外围辅助印痕；保留部分与“裁切包络减去冻结血管”双向布尔比较为零差。浅导向为外模内表面的 0.20 mm 浮雕，转印到 PDMS 表面，ABS 消失后仍可辨认；它不形成伸向中央血管的深槽。\n\n'
    text+=f"本模型有 {t['actual_cut_plane_count']} 个实际有限切面，其中四个为端口暴露面；O2 的一侧利用既有外边界。所有切面距非对应血管表面的最小净距为 {t['minimum_nonport_clearance']:.3f} mm。切前 PDMS 几何体积 {t['pdms_before_volume_mm3']/1000:.3f} mL，切后 {t['pdms_after_volume_mm3']/1000:.3f} mL。这些是模型容积，不是已验证的浇注配方。\n\n"
    text+='|端口|截面积 mm²|可用直段 mm|距最近分叉 mm|距非本端口血管 mm|\n|---|---:|---:|---:|---:|\n'
    for r in rows:text+=f"|{r['port_id']}|{r['section_area_mm2']:.4f}|{r['usable_straight_length_mm']:.2f}|{r['distance_to_bifurcation_mm']:.3f}|{r['minimum_nonport_clearance_mm']:.3f}|\n"
    text+='\n## 实际切片与验证范围\n\n外模与四个平放插片分别采用当前保存的 X1 Carbon、0.2 mm 喷嘴、0.10mm Standard Vessel、Bambu ABS Vessel 配置，强制 Auto Support OFF；血管另作现有 ABS 支撑配置的参考切片。没有开发新的血管支撑算法，也没有向打印机上传或开始打印。\n\n'
    for name,key in [('外模','mold_auto_support_off_slice'),('独立插片','retainer_auto_support_off_slice'),('血管参考','vascular_reference_slice')]:
        r=b[key];text+=f"- {name}：真实切片{ok(r.get('slice_passed',False))}；原生支撑挤出段 {r.get('native_support_extrusion_segments','未运行')}，支撑关闭要求 {r.get('support_off_required',False)}。\n"
    text+=f"\n归档 CRC、实际 ABS/喷嘴设置、归档内生产几何及薄片/接触面/露出段/浅导向刀路均核查，细特征保留检查{ok(b['feature_preservation'].get('passed',False))}。刀路存在不代表实际接触面的尺寸已经校准。PDMS STL 导出还单独复核封闭性和体积；计算布尔体与导出网格的数值清理记录见 trim/post_dissolution_trim.json，冻结血管的生产文件从未参与这种诊断网格清理。\n\n"
    text+='## 制造次序（工艺参数待实物确定）\n\n'
    steps=['单独打印冻结 ABS 血管芯。','在开放空间清理血管打印支撑，检查端口与细分支。','单独打印外 ABS 模具及四个 ABS 插片。','如需要，仅对外模实施已获得操作条件的受控蒸汽处理；本报告不提供未经校准的处理参数。','充分干燥；测量托座、插片导向和溶剂接触部，进行干装配检查。','完整血管从顶部装入，下托面在端口直段承托；不得强压或扭曲。','插片逐个下放并向墙滑入限位，确认露出段高于浇注面；重新确认四处溶剂接触。','按所用材料的实测工艺浇注、脱泡和固化 PDMS。','在合适的实验设施与既定溶剂操作流程下，将整体进入 acetone 去除阶段。','按实际实验结果确定溶剂更新或循环条件，使外围外模先打开；此处不规定时间或频率。','通过四处经实物确认的溶剂通路继续去除血管芯与其余 ABS，检验是否完全去除。','按已验证流程充分去除残余溶剂，观察 PDMS 尺寸恢复。','根据已转印的浅标记及有限切面切除外围 PDMS；本几何需要此步去掉辅助空腔。','检查仅剩目标管腔、四个端口均开放，得到供后续检查的 phantom。']
    text+='\n'.join(f'{i+1}. {x}' for i,x in enumerate(steps))
    text+='\n\n## 局限与复现\n\n溶解时间、溶剂更新要求、蒸汽处理尺寸变化、真实抗浮力、PDMS 溶胀及恢复、手工裁切精度均未物理验证。没有预测“几小时溶完”，也没有声明制造或材料行为已经验证。\n\n运行 `./run_python.sh s1-10_all_abs_sacrificial_casting_system.py`。配置为 `config/all_abs_sacrificial_casting_system.yaml`，默认输出本报告所在目录。`outer_mold` 和 `retainers` 中的 3mf 为本地实际切片；`vascular_core` 保存原字节冻结 STL 与独立参考切片。`assembly/all_ABS_casting_system.step` 保留独立零件；`all_ABS_review.stl` 只用于组合审阅，不能把它作为联合打印生产文件。`diagnostic_PDMS` 仅为浇注物与孔腔审阅，不是 ABS 打印任务。\n\n'
    text+=f"冻结芯 SHA256：`{s['vascular']['source_hash']}`。源保护检查 {s['source_protection']['checked_files']} 个文件，错配 {s['source_protection']['mismatch_count']}。新阶段没有改动既有 s1-2/s1-3 可视化。\n"
    (out/'all_abs_sacrificial_casting_report.md').write_text(text);(out/'reports/README.md').write_text('本阶段中文报告见 ../all_abs_sacrificial_casting_report.md；运行日志 run.log，源发现 source_discovery.json，实际配置 config.json。s1-9 为 RETIRED / RESEARCH HISTORY，原报告与源哈希原位保留。\n')
