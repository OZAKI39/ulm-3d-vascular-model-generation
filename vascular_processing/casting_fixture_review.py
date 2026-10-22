"""Standalone BCF review figures. Existing s1-2/s1-3 UI is never imported/edited."""
import json
from pathlib import Path
import numpy as np
import pyvista as pv
import trimesh
from .support_inputs import save,digest

NAMES=['01_frozen_vascular_core','02_bottomless_frame','03_frame_on_dish','04_no_PLA_bottom_section','05_top_loading_slots','06_vascular_loading_path','07_ports_in_saddles','08_master_datum_port','09_secondary_floating_ports','10_closure_inserts','11_anti_float_constraints','12_fixture_contact_map','13_no_central_vascular_contact','14_bottom_clearance','15_wall_clearance','16_bottom_flange_and_seal','17_port_seal_regions','18_final_assembled_fixture','19_PDMS_cavity_preview','20_exploded_assembly']


def poly(mesh):return pv.PolyData(mesh.vertices,np.c_[np.full(len(mesh.faces),3),mesh.faces].reshape(-1))

def render(out,layout,geometry,meshes,summary):
    out=Path(out);q=out/'QC';q.mkdir(exist_ok=True);H=layout['high'][2];ports=layout['ports'];lo=layout['low'];hi=layout['high'];frame=poly(meshes['frame']);core=poly(layout['core']);central=poly(layout['central']);records=[]
    floor=pv.Plane(center=(0,0,0),direction=(0,0,1),i_size=2*(abs(lo[0])+12),j_size=2*(abs(lo[1])+12))
    titles=['Independent frozen ABS vascular core','Bottomless PLA frame: four walls, no floor','Fixture on reference smooth dish floor','SECTION: PDMS bottom is dish surface, not PLA','Top-open U slots and lower saddles','Core loading: vertical downward translation','Straight port sections in wall-adjacent saddles','Master datum: '+layout['master_port']+' + external end stop','Secondary ports: axial floating channels','Separate removable closures and locking pins','Anti-float: closures positively retained by external pins','Capture zones only on extended ports','Central cerebral tree: NO CONTACT','Core surface to dish floor = 10.00 mm','Central tree to PLA: %.2f mm minimum'%summary['assembly']['central_min_clearance_mm'],'Outward flange and continuous underside seal channel','External port seal grooves: liquid tightness NOT tested','Final assembly: ABS core + separate PLA parts','PDMS casting cavity: dish directly forms its bottom','Exploded: dish / frame / core / closures / lock pins']
    for i,name in enumerate(NAMES):
        print('BCF QC',name,flush=True);p=pv.Plotter(off_screen=True,window_size=(1500,1100));p.set_background('white');p.add_text(titles[i],position='upper_left',font_size=16,color='#152433',font='arial')
        p.add_text('Dimensions: mm | Geometric review only | No physical seal/force validation',position='lower_left',font_size=10,color='#455666',font='arial')
        is_exploded=i==19;draw_floor=i in [2,3,5,6,10,12,13,14,15,17,18,19]
        if draw_floor:
            f=floor.copy();f.translate([0,0,-10 if is_exploded else 0],inplace=True);p.add_mesh(f,color='#d4dbe0',opacity=.7);p.add_point_labels([[0,lo[1]-8,-10 if is_exploded else 0]],['PETRI DISH FLOOR\nreference plane; actual dish not measured'],font_size=15,text_color='#26333f',point_size=0,shape_opacity=0,always_visible=True)
        if i not in (0,9):
            f=frame.copy()
            if i in (3,18):f=f.clip(normal=(1,1,0),origin=(0,0,0),invert=True)
            if i==15:f=f.clip(normal=(0,0,1),origin=(0,0,3),invert=True)
            p.add_mesh(f,color='#2888b7',opacity=.20 if i in (6,7,8,10,11,12,13,14,17,18,19) else .9,smooth_shading=True)
        if i not in (1,4,9,15,16):
            c=core.copy()
            if i in (5,19):c.translate([0,0,18 if i==5 else 30],inplace=True)
            p.add_mesh(c,color='#cf5548' if i!=12 else '#27a475',smooth_shading=True)
        if i in (6,7,8,9,10,11,16,17,18,19):
            for j,port in enumerate(ports):
                k=port['port_id'];show=(i!=7 or k==layout['master_port']) and (i!=8 or k!=layout['master_port'])
                if not show:continue
                for kind,color in [('closures','#46aa74'),('pins','#985ac8')]:
                    m=poly(meshes[kind][k]);shift=np.zeros(3)
                    if i==9:shift=np.array([0,0,12])
                    if is_exploded:shift=np.array([0,0,60 if kind=='closures' else 70])
                    m.translate(shift,inplace=True);p.add_mesh(m,color=color,opacity=.95,smooth_shading=True)
                c=np.array(port['wall_intersection']);label=k+(' : MASTER' if k==layout['master_port'] else f' : float +/-{port["axial_float_mm"]:.1f}')
                p.add_point_labels([c+[0,0,64 if is_exploded else 4]], [label],font_size=16,text_color='#18374b',point_size=4,point_color='#d83032',shape_opacity=.6,always_visible=True)
        if i in (4,5,10,19):
            for port in ports:
                c=np.array(port['wall_intersection']);p.add_arrows(np.array([c+[0,0,14]]),np.array([[0,0,-1.]]),mag=7,color='#bb7714')
        if i in (8,11):
            for port in ports:
                c=np.array(port['capture_center']);n=np.array(port['axis'])
                if i==8 and port['role']=='MASTER_DATUM':continue
                if i==8:p.add_arrows(np.array([c-n*2,c+n*2]),np.array([n,-n]),mag=3,color='#ac6919')
                else:
                    patch=pv.Cylinder(center=c,direction=n,radius=port['radius_mm']+.04,height=2,resolution=48);p.add_mesh(patch,color='#ebbc18');p.add_point_labels([c+[0,0,3]],[port['port_id']+' capture zone'],shape_opacity=0,always_visible=True,font_size=15)
        if i==12:p.add_point_labels([central.center],['NO CONTACT\npositive central clearance'],font_size=22,shape_opacity=.7,text_color='#076347',always_visible=True)
        if i==13:
            point=layout['core'].vertices[np.argmin(layout['core'].vertices[:,2])];p.add_lines(np.array([[point[0],point[1],0],point]),color='#be7427',width=4);p.add_point_labels([point/2],['10.00 mm'],font_size=20,always_visible=True)
        if i==15:
            loop=geometry['seal_loop'].copy();loop[:,2]=.35;p.add_lines(loop,connected=True,color='#bd6333',width=5);p.add_text('Blue = PLA, orange = outer seal path\nCavity bottom remains entirely open',position='upper_right',font_size=13,color='#243846')
        if i==16:
            for port in ports:
                c=np.array(port['wall_intersection'])+np.array(port['axis'])*2.02
                ring=pv.Cylinder(center=c,direction=port['axis'],radius=port['radius_mm']+.65,height=.06,resolution=64);p.add_mesh(ring,color='#d09228');p.add_point_labels([c],[port['port_id']+' exterior seal'],font_size=15,always_visible=True)
        if i==18:
            box=pv.Box(bounds=[lo[0],hi[0],lo[1],hi[1],0,H]);p.add_mesh(box,color='#8ed4ef',opacity=.08)
        if i==3:p.add_text('OPEN BOTTOM\nNo PLA floor / grid / crossbar',position='upper_right',font_size=17,color='#075079')
        p.add_axes();p.camera_position='iso';p.camera.up=(0,0,1);p.reset_camera();p.camera.zoom(1.05 if i!=19 else .95);p.show(screenshot=str(q/(name+'.png')),auto_close=True)
        records.append(dict(name=name,sha256=digest(q/(name+'.png'))))
    save(q/'figure_manifest.json',records)


def report(out,s):
    out=Path(out);f=s['frame'];a=s['assembly'];p=s['ports'];z=s['anti_float'];c=s['cleanup'];b=s['bambu'];status=s['status'];core=s['vascular']
    steps=['检查培养皿平坦底面尺寸，确认能容纳完整法兰；本轮尚未取得实际测量。','单独打印 ABS 血管芯；在开放空间完整清除其支撑并检查细血管。','分别用 PLA 打印无底框架、各端口闭合插片和锁销。','先干装检查端口余量、插片滑动及锁销；不要强压细血管。','在培养皿光滑底面定位框架，按外侧法兰密封路径建立临时密封。','将清理好的完整血管芯从顶部竖直放入，各直端口进入下托座。',f'先安装 {p["master_port"]} 主基准插片和锁销，再装其余插片和锁销。','检查外侧端口密封区域，先进行实际漏液和定位试验；本轮未验证胶黏剂或密封性能。','试验通过后倒入 PDMS 并按实际材料工艺固化；观察血管芯是否移动。','移除锁销、插片和 PLA 框架，取出 PDMS；之后才以既定 acetone 流程去除内部 ABS。']
    text=f'''# 独立 ABS 血管芯与无底 PLA 浇注夹具

状态：**{status}**。这是几何、装配和本地切片审阅，不是制造、漏液或抗浮力实验认证。

## 为什么分开制造

过去血管、盒底与复杂支撑共处于狭小空间，清理支撑容易受阻。本轮把 ABS 血管芯独立打印，先在开放环境中完成清理，再装入另行打印的 PLA 框架。血管的中心线、半径、端口连接和已有可视化均未重建。

## 底面和接触

PLA 框架没有底板：**{'NO PLA BOTTOM' if f['bottom_open'] else '检查未通过'}**。四侧壁下缘与向外法兰共面，未来 PDMS 底面直接接触培养皿光滑底面。没有内部底网、横梁或中央夹点。

共有 {a['nominal_capture_sites']} 处端口捕获位置，全部位于冻结的延长直段。名义装配保留间隙，所以名义硬接触数是 {a['nominal_hard_contact_count']}；这与“4 个承载/限位位置”不是同一统计。中央脑血管是否接触夹具：**{'YES' if a['central_vascular_contact_count'] else 'NO'}**。中央表面到 PLA 的最小三角面距离为 {a['central_min_clearance_mm']:.3f} mm。内腔非法 PLA 侵入体积 {a['illegal_fixture_intrusion_volume']:.9f} mm³；底面净距 {a['core_to_dish_clearance_mm']:.3f} mm。

## 装配与基准

完整血管芯纯竖直向下装入：{'通过' if a['top_loading_path']['passed'] else '未通过'}。FCL 表面距离加区间运动上界核查整条路径，不仅检查起终点。各闭合插片及锁销插入：{'通过' if a['closure_insertion_pass'] else '未通过'}。

主基准为 **{p['master_port']}**。各端口保守可用直段长度接近，按最长直段、较大直径、分叉间距依次选择；主基准额外设外侧端面参考挡块。其余端口采用轴向贯通捕获槽，保守可用浮动为 {p['axial_float']} mm（值表示 ± 范围）。径向装配余量为 {p['clearances']} mm；来自本机线宽的候选搜索，尚未用实物收缩数据校准。

## 抗浮与闭合

下托座和上插片形成端口包围限位。插片由外侧横向锁销连接侧壁，不能仅凭“插片从上方滑入”就认为它不会被浮力抬起。当前 +Z、pitch、roll 几何约束为 {z['geometric_Z_lock']}、{z['pitch_lock']}、{z['roll_lock']}；正 Z 自由行程上界约 {z['limits']['Z']['positive']:.3f} mm。小位姿线性约束与实际三角网格的大位移碰撞分别保存；不把 secondary 轴向浮动误当成刚性锁死。

U 槽由可拆插片关闭，装配缝和端口周围具有外侧密封槽；几何接近连续并不表示 PLA 配合天然液密。法兰底部另有单一闭合周边密封路径，所有密封特征均在内腔外。

## 培养皿与打印

内腔尺寸为 {f['inner_dimensions']} mm；墙厚 {f['wall_thickness']:.2f} mm，向外法兰宽 {f['flange_width']:.2f} mm。完整法兰的最小平坦圆形容纳直径为 **{f['minimum_required_dish_flat_diameter']:.3f} mm**。它取矩形法兰对角线，两个对角顶点给出同一必要且充分下界；没有假定用户培养皿是 90、100 或 150 mm。当前培养皿状态：**{s['dish']['status']}**。

本机配置是 X1 Carbon / 0.2 mm 喷嘴 / 已安装 PLA profile，未更改 GUI 预设。框架 Auto Support OFF 切片：{b.get('frame_slice',{}).get('support_free_slice_passed',False)}；插片和锁销独立摆放后 Auto Support OFF 切片：{b.get('closure_slice',{}).get('support_free_slice_passed',False)}。同时检查归档设置、自动支撑刀路和切片警告。独立 ABS 芯仅使用已有自动支撑配置做参考切片，不研发新支撑。

## 尚未实验验证

PHYSICAL_SEAL_NOT_VALIDATED；PHYSICAL_LEAK_TEST_NOT_DONE；ANTI_FLOAT_FORCE_NOT_PHYSICALLY_VALIDATED。ABS/PLA 收缩、公差、PLA 插片及锁销强度、浇注浮力、夹具取出和真实密封均待试验。未推荐或认证任何具体胶黏剂；未上传、发送或启动打印机。

## 下一次实际制造顺序

'''+ '\n'.join(f'{i+1}. {v}' for i,v in enumerate(steps))+f'''

## 清理与复现

删除 {c['deleted_file_count']} 个退休支撑输出及缓存，释放 {c['deleted_bytes']:,} 字节。清理前 {c['before_file_count']} 文件、{c['before_bytes']:,} 字节；清理后原库存 {c['after_inventory_file_count']} 文件、{c['after_inventory_bytes']:,} 字节（不含本轮新交付物）。没有 Git 历史可恢复，因此源码和原测试保留，s1-8 仅为研究历史。旧 UI 哈希及迁移别名冲突的明确维护记录在 research_retrospective/retained_test_maintenance，不改旧报告或生产输入。

运行：`./run_python.sh s1-9_vascular_casting_fixture.py`。所有输出位于本报告父目录；STL/STEP 在 frame、closures、assembly，血管原样副本在 vascular_core，真实切片在 bambu。assembly/full_fixture_review.stl 是多零件审阅文件，不是供单材料合并打印的文件。

冻结血管 SHA256：`{core['source_hash']}`。上游几何/数据保护与永久测试结果见 reports；失败门槛：{s['failures']}。
'''
    (out/'vascular_casting_fixture_report.md').write_text(text)
