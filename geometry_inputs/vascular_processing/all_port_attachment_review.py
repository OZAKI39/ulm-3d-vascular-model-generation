"""Fixed-camera Chinese before/after evidence for all four attachments."""
from pathlib import Path
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw
from . import sacrificial_fixture as f, surface_continuity_qc as qc
from .sacrificial_fixture_review import chinese_font, polydata, save_json, core_face_labels

COLOR='#10a7b5'


def title(image, text, cfg, second=None):
    draw=ImageDraw.Draw(image);height=82 if second else 54
    draw.rectangle((0,0,image.width,height),fill='white')
    draw.text((18,8),text,font=chinese_font(cfg,26),fill='#172b36')
    if second:draw.text((18,45),second,font=chinese_font(cfg,19),fill='#505050')
    return image


def local_picture(mesh,e,cfg,camera,wire=False,features=False):
    local=qc.window(mesh.triangles_center,e.position,e.tangent,e.radius*2.2,1.65)
    part=mesh.submesh([np.flatnonzero(local)],append=True,repair=False)
    p=pv.Plotter(off_screen=True,window_size=(760,720));p.set_background('white')
    p.add_mesh(polydata(part),color=COLOR,show_edges=wire,edge_color='#404040',smooth_shading=False)
    if features:
        settings=cfg['surface_continuity']
        _,lines=qc.feature_edges(mesh,e.position,e.tangent,e.radius*settings['radial_radius_fraction'],settings['seam_half_width_mm'],20)
        if lines.n_cells:p.add_mesh(lines,color='#d72d35',line_width=5)
    p.camera_position=camera['position'];p.camera.parallel_projection=True
    p.camera.parallel_scale=camera['parallel_scale'];p.camera.clipping_range=camera['clipping_range']
    image=Image.fromarray(p.screenshot(return_img=True))
    observed=dict(position=[list(v) for v in p.camera_position],parallel_scale=float(p.camera.parallel_scale),
        clipping_range=list(p.camera.clipping_range))
    np.testing.assert_allclose(observed['position'],camera['position'],atol=1e-10)
    np.testing.assert_allclose(observed['clipping_range'],camera['clipping_range'],atol=1e-10)
    if observed['parallel_scale'] != camera['parallel_scale']:
        raise ValueError('REVIEW_CAMERA_ZOOM_CHANGED')
    camera.setdefault('rendered_camera_checks',[]).append(observed)
    p.close()
    return image


def render_review(out,baseline,core,inputs,routes,cfg,summary):
    out=Path(out);qdir=out/'QC';qdir.mkdir(exist_ok=True)
    cameras={};files=[];before_pics=[];feature_pairs=[]
    for j,e in enumerate(inputs['endpoints']):
        mean=np.array(e.attachment_qc['reference_swc_mean_tangent'])
        view=np.cross(mean,e.tangent)
        if np.linalg.norm(view)<1e-8:view=qc.axes(e.tangent)[0]
        view=f.unit(view)
        camera=dict(position=[(e.position+view*9).tolist(),e.position.tolist(),e.tangent.tolist()],
            parallel_scale=1.85,clipping_range=[.1,30.],mesh_color=COLOR,
            window_size=[760,720],crop_radial_mm=e.radius*2.2,crop_axial_half_length_mm=1.65)
        cameras[e.endpoint_id]=camera
        for wire,offset,mode in [(False,2,'shaded'),(True,3,'wireframe')]:
            panels=[]
            for label,mesh in [('修复前',baseline),('真实封口对齐后',core)]:
                img=local_picture(mesh,e,cfg,camera,wire)
                title(img,e.endpoint_id+' · '+label,cfg,'同一相机、缩放、裁剪范围与颜色')
                panels.append(img)
            combined=Image.new('RGB',(1520,720),'white')
            combined.paste(panels[0],(0,0));combined.paste(panels[1],(760,0))
            name=f'{offset+j*2:02}_{e.endpoint_id}_before_after_{mode}.png';combined.save(qdir/name);files.append(name)
            if wire:before_pics.append(panels[0])
        feature_panels=[]
        for label,mesh,key in [('前',baseline,'baseline_continuity'),('后',core,'new_continuity')]:
            metrics=summary['per_port'][e.endpoint_id][key]['features'][1]
            img=local_picture(mesh,e,cfg,camera,True,True)
            title(img,e.endpoint_id+' · '+label+'：20° 折边',cfg,f"红线总长 {metrics['total_length_mm']:.4f} mm")
            feature_panels.append(img)
        feature_pairs.append(feature_panels)
    canvas=Image.new('RGB',(1520,1440),'white')
    for j,img in enumerate(before_pics):canvas.paste(img,((j%2)*760,(j//2)*720))
    name='01_all_ports_baseline.png';canvas.save(qdir/name);files.append(name)
    canvas=Image.new('RGB',(1520,2880),'white')
    for j,pair in enumerate(feature_pairs):
        for k,img in enumerate(pair):canvas.paste(img,(k*760,j*720))
    name='10_feature_edges_before_after.png';canvas.save(qdir/name);files.append(name)
    save_json(out/'camera_settings.json',dict(before_after_identical=True,per_port=cameras,
        shading='Flat geometric shading; normals are not visually smoothed'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    font=FontProperties(fname=cfg['visualization']['chinese_font'])
    names=[e.endpoint_id for e in inputs['endpoints']];x=np.arange(4)
    old=[summary['per_port'][p]['baseline_continuity']['normal_jumps']['p95_deg'] for p in names]
    new=[summary['per_port'][p]['new_continuity']['normal_jumps']['p95_deg'] for p in names]
    fig,ax=plt.subplots(figsize=(10,5))
    ax.bar(x-.18,old,.36,label='修复前',color='#8e9aa3');ax.bar(x+.18,new,.36,label='对齐后',color=COLOR)
    ax.set_xticks(x,names);ax.set_ylabel('相邻面夹角 p95（度）',fontproperties=font)
    ax.set_title('接口附近的面夹角：相同范围、相同算法',fontproperties=font)
    ax.set_ylim(0,max(20,max(old+new)*1.2));ax.legend(prop=font);ax.grid(axis='y',alpha=.2)
    fig.tight_layout();name='11_normal_jump_comparison.png';fig.savefig(qdir/name,dpi=160);plt.close(fig);files.append(name)
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,name in zip(axes.flat,names):
        for label,key,color in [('修复前','baseline_continuity','#8e9aa3'),('对齐后','new_continuity',COLOR)]:
            rows=summary['per_port'][name][key]['cross_sections']
            ax.plot([r['arc_mm'] for r in rows],[r['area_mm2'] for r in rows],label=label,color=color)
        ax.set_title(name);ax.set_xlabel('沿封口法向的距离（mm）',fontproperties=font)
        ax.set_ylabel('截面积（mm²）',fontproperties=font);ax.legend(prop=font);ax.grid(alpha=.2)
    fig.suptitle('封口前后截面积连续性',fontproperties=font);fig.tight_layout()
    name='12_cross_section_comparison.png';fig.savefig(qdir/name,dpi=160);plt.close(fig);files.append(name)
    positions=np.array([e.position for e in inputs['endpoints']])
    for colored,filename in [(False,'13_final_full_core.png'),(True,'14_all_ports_final_labels.png')]:
        p=pv.Plotter(off_screen=True,window_size=(1500,1100));p.set_background('white')
        mesh=polydata(core)
        if colored:
            mesh.cell_data['role']=core_face_labels(core,inputs['mesh'],routes)
            p.add_mesh(mesh,scalars='role',cmap=[COLOR,'#23ac62','#f08c28'],clim=[0,2],show_scalar_bar=False)
        else:p.add_mesh(mesh,color=COLOR,smooth_shading=False)
        p.add_point_labels(positions,names,text_color='black',point_color='#d72d35',font_size=24,
            point_size=8,always_visible=True,shape_opacity=.6)
        p.view_isometric();p.reset_camera()
        img=Image.fromarray(p.screenshot(return_img=True));p.close()
        title(img,'完整血管与四个真实封口对齐端口',cfg,'原血管保持不变；图中未包含盒子')
        img.save(qdir/filename);files.append(filename)
    return dict(png_files=sorted(files),same_camera_before_after=True,camera_settings=str(out/'camera_settings.json'))


def write_report(out,s):
    rows=[]
    for port,d in s.get('per_port',{}).items():
        if 'new_continuity' not in d:continue
        old,new=d['baseline_continuity'],d['new_continuity'];c=d['cap']
        rows.append(f"| {port} | {old['features'][1]['total_length_mm']:.4f} → {new['features'][1]['total_length_mm']:.4f} | "
            f"{old['normal_jumps']['p95_deg']:.4f} → {new['normal_jumps']['p95_deg']:.4f} | "
            f"{new['features'][1]['circumferential_ring_count']} | {100*new['maximum_area_jump_fraction']:.3f}% | {d['clearance_mm']:.4f} | {d['status']} |")
    details=s.get('per_port',{})
    improvements={p:d['baseline_continuity']['features'][1]['relevant_circumferential_length_mm']-d['new_continuity']['features'][1]['relevant_circumferential_length_mm']
                  for p,d in details.items() if 'new_continuity' in d and p!='O3'}
    best=max(improvements,key=improvements.get) if improvements else '尚未确定'
    o3=details.get('O3',{});o3passed=o3.get('gate',{}).get('passed',False)
    text=f'''# 四个血管接口的真实封口对齐与连续性复核

## 1. 之前真正的问题在哪里

实际读取旧 core STL 后可见，主要硬接缝位于 I1、O1、O2。其 20° 折边总长分别约为 3.7955、3.8222、4.1462 mm，接近各自一圈周长。O3 已经采用真实封口对齐，旧模型的 20° 折边长度为零，本轮将其作为良好对照。图片中的位置标记用于避免把邻近分支误认作 O3。

## 2. 为什么会出现问题

旧 I1/O1/O2 用约 5 mm 中心线的平均方向连接人工圆管，并将整段等半径圆管向原血管内部延伸。弯曲血管的真实截面朝向与平均方向并不相同，伸入段可能露出原表面，形成斜接或环状台阶。方向差本身不是失败条件；O3 虽相差约 42°，真实封口对齐后仍可连续衔接。

## 3. 这次改了什么

四个接口复用同一封口匹配方法，通过 fitted SWC、节点映射、manifest 和冻结 STL 的中心顶点确认身份，读取 24 点真实边界及面积加权法向。生产轮廓使用实际多边形，半径不缩放。埋入段为长度优先 0.20 mm、内端半径比例 0.65 的短锥台，只有完全通过原有 10⁻⁶ mm³ 外露体积限制才接受；允许的长度回退仅为 0.15、0.10 mm。

I1 的 STL 封口中心位于轮廓点的轴向误差范围之外，平均平面因此产生微小外露薄片。按照允许的真实封口平面投影，将同一法向的平面定位到已确认中心顶点的轴向位置；最大点投影约 0.00000632 mm，小于 0.0001 mm 上限。没有修改冻结 STL，也没有提高容差。O3 不触发这一条件，其封口轮廓、法向、埋入段和路线保持原数值。所有平面、投影量和尝试结果均保存在逐接口 JSON 中。

I1/O1/O2 仅重算连接到原目标前 3 mm 位置的第一弯段。其两端位置与方向固定，用单次局部五次曲线求解满足弯曲半径；没有运行盒壁分配或目标搜索。四个原目标点、外部直管、盒子和端口半径保持不变。新 core 仅由原血管和四个本轮端口依次合并。

布尔计算会在封口附近留下相距数微米的重复顶点，造成极薄的折叠三角形。对本轮合并结果，仅在接口 ±1 mm 范围内用原有 0.00001 mm 数值清理容差合并这些重复点，并优先保留原血管顶点；不平均移动点、不补洞、不重铺全模型网格。清理前后的顶点编号、最大位移、面数和体积均留有记录，原始 STL 只读。该步骤用于去除布尔数值碎片，不改变血管半径或用显示平滑掩盖问题。

## 4. 四个接口逐一结果

| 接口 | 20° 折边长度：前 → 后（mm） | 面夹角 p95：前 → 后（°） | 新闭合锐环 | 最大相邻截面积变化 | 血管间距（mm） | 状态 |
|---|---:|---:|---:|---:|---:|---|
'''+ '\n'.join(rows)+f'''

接口分析位于 4 mm 外层观察范围内，统一使用法向 ±1 mm、径向 1.7 倍半径的接缝带，避免混入旁支和远处弯段。这与上一轮诊断范围一致。20° 折边还按“距封口平面不超过 0.75 mm，边方向与局部圆周切向夹角不超过 45°”筛选，得到锐边长度占周长比例；比例允许超过 1，因为两个台阶可能形成两圈。24 边形的纵向棱边不计作环状接缝。

p95 使用接缝带内全部相邻面，不能单独代表环状硬肩：旧模型大量平滑三角形也可能使 p95 较低。报告保留真实结果；若没有下降 50%，只能依据约定的绝对值 ≤20° 判据及锐环、圆周比例、实际图片联合判断，不能声称达到 50% 改善。最大夹角仅作辅助。截面为沿封口法向 −0.5 至 +1.0 mm、每 0.1 mm 的切片，arc_mm 是有符号法向距离，不冒充弯曲中心线弧长。

## 5. 哪个改善最大

以环向锐边实际减少长度衡量，改善最大的是 {best}。详细的前后值见 tables/all_port_continuity_summary.csv；02–09 为同相机、同缩放、同裁剪、同颜色的实体及线框图，10 为统一 20° 红色折边图。表面采用真实几何着色，没有通过显示平滑隐藏接缝。

## 6. O3 有没有退化

{'O3 通过回归条件：新 p95 不超过旧值加 1°，20° 折边仍为零，间距不低于 1.6 mm。' if o3passed else 'O3 尚未通过完整回归条件，不能认定没有退化。'} 其封口接法和路线保持原样；整体合并可能改变少量三角形分割，因此仍对最终交付 STL 重新测量。

## 7. 是否还需要新版 VMTK

`{s.get('vmtk_status','NOT_EVALUATED')}`

本轮没有安装、升级或调用生产 VMTK ramp。只有对齐后仍有闭合锐环或 p95 超过 30°，且人工图像确认明显硬接缝，才建议升级处理；普通纵向多边形棱线不构成升级理由。

## 8. 是否可以继续 box/core union

`{s['status']}`

这表示当前自动验证结果，仍须人工查看 02–14 号图后决定下一阶段。未执行盒子与血管的一体合并，也未执行切片。受保护文件共 {s.get('protected_file_count',0)} 个，哈希保持不变：{s.get('protected_sources_unchanged',False)}。

最终模型：`core/BG001_RMCA_BALANCED_core_with_ports_all_aligned.stl`。完整机器结果在 all_port_attachment_summary.json；永久测试与真实运行日志分别见 test_results.xml、run.log。运行入口为 `s1-4_sacrificial_box_and_ports.py --all-ports-aligned`，输出目录已存在时拒绝覆盖，可用 `--output-root` 指定另一独立审核目录。
'''
    if s.get('failures'):text+='\n未通过项目：\n\n'+'\n'.join('- '+v for v in s['failures'])+'\n'
    (Path(out)/'all_port_attachment_report.md').write_text(text,encoding='utf-8')
