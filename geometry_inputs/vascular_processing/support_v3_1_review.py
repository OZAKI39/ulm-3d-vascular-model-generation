"""Evidence figures and Chinese outcome report for the two sequential studies."""
from __future__ import annotations
import json,pickle,shutil
from pathlib import Path
import numpy as np
import trimesh
import shapely
from .support_inputs import save,table,transform


def pvmesh(mesh):
    import pyvista as pv
    return pv.PolyData(mesh.vertices,np.c_[np.full(len(mesh.faces),3),mesh.faces].ravel())


def polygons_mesh(items):
    vertices=[];faces=[];offset=0
    for z,p in items:
        for poly in shapely.get_parts(p):
            if poly.geom_type!='Polygon' or poly.area<1e-10:continue
            xy,tri=trimesh.creation.triangulate_polygon(poly)
            vertices.extend(np.c_[xy,np.full(len(xy),z)]);faces.extend(tri+offset);offset+=len(xy)
    return trimesh.Trimesh(vertices=np.asarray(vertices).reshape(-1,3),faces=np.asarray(faces,dtype=int).reshape(-1,3),process=False)


def image(path,title,meshes=(),points=(),lines=(),note='',camera='iso'):
    import pyvista as pv
    plot=pv.Plotter(off_screen=True,window_size=(1450,1050));plot.set_background('white')
    for mesh,color,opacity in meshes:
        if mesh is not None and len(mesh.faces):plot.add_mesh(pvmesh(mesh),color=color,opacity=opacity)
    for xyz,color,size in points:
        if len(xyz):plot.add_points(np.asarray(xyz),color=color,point_size=size,render_points_as_spheres=True)
    for xyz,color,width in lines:
        if len(xyz)>1:plot.add_lines(np.asarray(xyz),color=color,width=width,connected=True)
    plot.add_text(title,font_size=17,color='#17212b',position='upper_left')
    plot.add_text(note,font_size=11,color='#333333',position='lower_left');plot.add_axes()
    if camera=='top':plot.view_xy()
    else:plot.view_isometric()
    plot.reset_camera();plot.camera.zoom(.87);path.parent.mkdir(parents=True,exist_ok=True);plot.screenshot(str(path));plot.close()


def orientation_images(context,qc):
    import pyvista as pv
    base,out,inputs,profile,discovery,mold,v,box,branches,top,rows,baseline,metrics=context
    bottom=box.submesh([np.flatnonzero(box.triangles_center[:,2]<=inputs['inner_bounds'][0][2]+1e-4)],append=True)
    def render(names,path,columns=2):
        nr=int(np.ceil(len(names)/columns));plot=pv.Plotter(off_screen=True,shape=(nr,columns),window_size=(1800,max(600,nr*480)))
        for i,name in enumerate(names):
            row=next(r for r in rows if r['orientation_id']==name);m=np.asarray(row['print_transform'])
            plot.subplot(i//columns,i%columns);plot.set_background('white')
            for original,color,opacity in [(box,'#cad0d5',.15),(bottom,'#d99537',.3),(v,'#aa2634',1.)]:
                mesh=original.copy();mesh.apply_transform(m);plot.add_mesh(pvmesh(mesh),color=color,opacity=opacity)
            vertices=transform(mold.vertices,m);lo,hi=vertices.min(0),vertices.max(0)
            plate=pv.Plane(center=[*(lo[:2]+hi[:2])/2,-.1],direction=[0,0,1],i_size=hi[0]-lo[0]+8,j_size=hi[1]-lo[1]+8)
            plot.add_mesh(plate,color='#bccbd5',opacity=.25)
            opening=np.array(inputs['inner_bounds']).mean(0);opening[2]=inputs['inner_bounds'][1][2]
            point=transform([opening],m);normal=np.asarray(row['opening_normal_print']).reshape(1,3)
            plot.add_arrows(point,normal,mag=12,color='#1268ad')
            # Duplicate the same vector outside the translucent box so an
            # inverted ceiling can never hide the required direction glyph.
            glyph=np.array([[lo[0]-15,hi[1]+15,(lo[2]+hi[2])/2]])
            plot.add_arrows(glyph,normal,mag=15,color='#1268ad')
            plot.add_text('Red: vessels | Amber: original bottom\nBlue: opening normal | Grey: plate',font_size=8,color='#333333',position='lower_right')
            plot.add_text(f"{name} | down-angle {row['angle_to_print_down_deg']:.2f} deg",font_size=12,color='black')
            plot.add_axes();plot.view_isometric();plot.reset_camera();plot.camera.zoom(.8)
        plot.screenshot(str(path));plot.close()
    render([r['orientation_id'] for r in rows],qc/'01_orientation_grid.png',4)
    render(['INV_0'],qc/'02_INV0.png',1)
    for n,d in [(3,15),(4,30),(5,45)]:render([r['orientation_id'] for r in rows if r['orientation_id'].endswith(str(d))],qc/f'{n:02}_INV{d}_group.png',2)


def actual_scene(folder):
    with (folder/'semantic_coverage_geometry.pkl').open('rb') as f:geo=pickle.load(f)
    hits=[];miss=[]
    for h,m in geo.values():hits+=h;miss+=m
    return [(trimesh.load_mesh(folder/'vascular_with_ports.stl'),'#9ba3aa',.15),
            (polygons_mesh(hits),'#198754',1.),(polygons_mesh(miss),'#d52232',1.)]


def compare_actual(left,right,path):
    import pyvista as pv
    plot=pv.Plotter(off_screen=True,shape=(1,2),window_size=(2100,1050))
    for i,folder in enumerate((left,right)):
        r=json.loads((folder/'candidate.json').read_text());plot.subplot(0,i);plot.set_background('white')
        for mesh,color,alpha in actual_scene(folder):
            if len(mesh.faces):plot.add_mesh(pvmesh(mesh),color=color,opacity=alpha)
        cov=r['slice']['actual_vascular_coverage'];plot.add_text(f"{'V3 / O000 baseline' if i==0 else 'Final V3.1 diagnostic'}\nActual vascular + ports: {cov['percent']:.4f}%\nUncovered layer regions: {cov['uncovered_region_count']}",font_size=16,color='black')
        plot.add_text('Green: actual G-code supported | Red: actual unsupported\nRequired layer polygons; no bridge credit',font_size=12,color='black',position='lower_left')
        plot.add_axes();plot.view_isometric();plot.reset_camera();plot.camera.zoom(.85)
    plot.screenshot(str(path));plot.close()


def phase_a_qc(context,records,best):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    base,out,inputs,profile,discovery,mold,v,box,branches,top,rows,baseline,metrics=context
    qc=out/'QC_INVERTED';qc.mkdir(exist_ok=True);orientation_images(context,qc)
    fig,ax=plt.subplots(figsize=(12,6));names=[r['orientation_id'] for r in metrics]
    ax.bar(names,[r['box_bottom_unsupported_burden_mm2'] for r in metrics],label='Bottom slab (ceiling)',color='#d99537')
    ax.bar(names,[r['box_wall_unsupported_burden_mm2'] for r in metrics],bottom=[r['box_bottom_unsupported_burden_mm2'] for r in metrics],label='Walls',color='#5f84ac')
    ax.set_ylabel('Coarse downward face area (mm²)');ax.tick_params(axis='x',rotation=60);ax.legend();fig.tight_layout();fig.savefig(qc/'06_box_ceiling_burden.png',dpi=180);plt.close(fig)
    surrogate=next(r for r in metrics if r['orientation_id']==best['candidate_id']);folder=Path(best['folder'])
    model=trimesh.load_mesh(folder/'mold_without_support.stl');vascular=trimesh.load_mesh(folder/'vascular_with_ports.stl');support=trimesh.load_mesh(folder/'custom_support_v3.stl')
    vis=surrogate['visibility'];xyz=np.asarray(vis['points_print']);mask=np.array(vis['visible_mask'])
    image(qc/'07_cavity_buildplate_visibility.png',f"{best['candidate_id']} | opening visibility",[(model,'#999999',.12)],[(xyz[mask],'green',7),(xyz[~mask],'red',7)],note='Coarse rays must cross the real opening; not release proof')
    xyz=np.asarray(surrogate['direct_targets_print']);mask=np.array(surrogate['direct_clear_mask'])
    image(qc/'08_direct_buildplate_support.png',f"{best['candidate_id']} | direct plate opportunity",[(model,'#999999',.1)],[(xyz[mask],'green',3),(xyz[~mask],'red',4)],note='Area-weighted surface-ray surrogate; final coverage uses actual G-code')
    image(qc/'09_best_inverted_support.png',f"Best diagnostic inverted: {best['candidate_id']}",[(vascular,'#c23e47',.8),(model,'#cccccc',.15),(support,'#408cbd',.7)],note='Not a production selection; failed gates remain in full_candidate_metrics.csv')
    release=best['plate_release'];attempt=release.get('path') or release['attempts'][0];poses=np.array(attempt['transforms']);local=model.copy();local.apply_translation(-np.array(attempt['moving_local_pivot']))
    first,last=poses[0],poses[-1]
    from scipy.spatial.transform import Rotation,Slerp
    mid=np.eye(4);mid[:3,:3]=Slerp([0,1],Rotation.from_matrix(np.array([first[:3,:3],last[:3,:3]])))(.5).as_matrix();mid[:3,3]=(first[:3,3]+last[:3,3])/2
    if len(poses)>2:mid=poses[len(poses)//2] # an actual tested tilt/upward keyframe
    fixed=trimesh.load_mesh(folder/'plate_rooted_support.stl')
    for number,label,pose in [(10,'start',first),(11,'mid',mid),(12,'end',last)]:
        moved=local.copy();moved.apply_transform(pose)
        image(qc/f'{number:02}_plate_release_{label}.png',f"{best['candidate_id']} | plate release {label}",[(moved,'#aaaaaa',.4),(fixed,'#2977ad',1.)],note=f"Swept result: {release['status']} | a failed path is a diagnostic attempt, not a feasible motion")
    image(qc/'13_actual_coverage_inverted.png',f"Actual G-code coverage | {best['candidate_id']}",actual_scene(folder),note='Green supported; red unsupported; vascular, ports and box-required polygons included')
    control=next(r for r in records if r['candidate_id']=='UPRIGHT_O000')
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    names=['VASCULAR','PORT','BOX_BOTTOM','BOX_WALLS']
    for i,r in enumerate((control,best)):
        axes[i].bar(names,[r['semantic_coverage'][n]['percent'] for n in names],color=['#bb4545','#b39151','#6999bb','#69a8a0']);axes[i].set_ylim(0,105);axes[i].set_title(r['candidate_id']);axes[i].tick_params(axis='x',rotation=25);axes[i].set_ylabel('Actual required-area coverage (%)')
    fig.suptitle('100% with zero required area is vacuous; see region counts in CSV');fig.tight_layout();fig.savefig(qc/'14_O000_vs_best_inverted.png',dpi=180);plt.close(fig)


def phase_b_qc(base,b):
    if not b['executed']:return
    out=base/'v3_1';qc=out/'QC_V3_1';qc.mkdir(exist_ok=True)
    folders=[Path(x['folder']) for x in b['iterations']];first,last=folders[0],folders[-1];one=folders[min(1,len(folders)-1)];two=folders[min(2,len(folders)-1)]
    for folder,n,title in [(first,'01_baseline_uncovered','Baseline actual coverage'),(one,'04_iteration1_actual_coverage','Iteration 1 actual coverage'),(last,'09_final_actual_coverage','Final actual coverage')]:
        image(qc/(n+'.png'),title,actual_scene(folder),note='Actual Bambu support extrusion footprints; green supported / red missing')
    for folder in folders:image(folder/'QC/actual_coverage.png',folder.name,actual_scene(folder),note='Actual G-code evidence; not nominal support geometry')
    v=trimesh.load_mesh(last/'vascular_with_ports.stl');support=trimesh.load_mesh(last/'custom_support_v3.stl')
    grains=json.loads((one/'residual_grains.json').read_text()) if (one/'residual_grains.json').exists() else []
    pts=np.array([np.mean(g['centroid_path'],axis=0) for g in grains]).reshape(-1,3)
    image(qc/'02_residual_grains.png',f'Branch/layer residual grains: {len(grains)}',[(v,'#aaaaaa',.3)],[(pts,'#c62736',9)],note='Clustering uses component, SWC branch, arc adjacency, XY proximity and tangent similarity')
    a=trimesh.load_mesh(first/'custom_support_v3.stl');c=trimesh.load_mesh(one/'custom_support_v3.stl')
    delta=trimesh.boolean.difference([c,a],engine='manifold');image(qc/'03_iteration1_repairs.png','Iteration 1 accepted local additions',[(a,'#aeb5bc',.2),(delta,'#e28d28',1.),(v,'#bb5555',.35)],note=f'Accepted added volume: {abs(float(delta.volume)):.5f} mm3; rejected edits are logged')
    d=trimesh.load_mesh(two/'custom_support_v3.stl');delta2=trimesh.boolean.difference([d,c],engine='manifold');image(qc/'05_iteration2_delta.png','Iteration 2 support change',[(c,'#aaaaaa',.2),(delta2,'#e28d28',1.)],note=f'Added volume: {abs(float(delta2.volume)):.5f} mm3; no whole-field dilation')
    for folder,name in [(first,'06_support_self_residual'),(last,'10_final_support_self_check')]:
        r=json.loads((folder/'candidate.json').read_text());sr=r['slice']['actual_self_buildability']['rows']
        meshes=[(trimesh.load_mesh(folder/'custom_support_v3.stl'),'#aaaaaa',.25),(polygons_mesh([(x['z_mm'],shapely.from_wkt(x['polygon_wkt'])) for x in sr]),'#d52232',1.)]
        image(qc/(name+'.png'),f"Actual support-self residual: {len(sr)} regions",meshes,note='Previous actual support/box footprint + min(width/2, layer height) reach')
    final=json.loads((last/'candidate.json').read_text());rem=final['removal'];meshes=[(v,'#bbbbbb',.35)]
    restore=np.asarray(final['orientation']['casting_restore_transform']);inv=np.linalg.inv(restore)
    for fr in rem['fragments']:
        path=last/'removal_fragments'/f"{fr['fragment_id']}.stl"
        if fr['status']=='TRAPPED_FRAGMENT' and path.exists():
            m=trimesh.load_mesh(path);m.apply_transform(inv);meshes.append((m,'#d52232',1.))
    image(qc/'07_trapped_fragment_repairs.png',f"Remaining trapped fragments: {rem['trapped_pieces']}",meshes,note='All attempted shorten/split/window/notch/root edits are recorded separately')
    tools=json.loads((last/'tool_access.json').read_text());good=[];bad=[]
    for t in tools:
        end=t.get('target') or (t['attempts'][0]['end'] if t.get('attempts') else None)
        if end is not None:(good if t['passed'] else bad).append(transform([end],inv)[0])
    image(qc/'08_cut_access_repairs.png',f"3 mm tool proxy | inaccessible {len(bad)}",[(v,'#aaaaaa',.3),(support,'#8aacc6',.3)],[(good,'green',9),(bad,'red',11)],note='Cut relocation attempts preserve coverage; a geometric tool proxy is not physical validation')
    lines=[]
    for fr in rem['fragments']:
        if fr.get('path'):
            poses=np.array(fr['path']['transforms']);lines.append((transform(poses[:,:3,3],inv),'#228452',2))
    image(qc/'11_final_removal_paths.png','Final ordered whole-fragment removal paths',[(v,'#bb5555',.35),(support,'#aaaaaa',.25)],lines=lines,note=f"Found {len(rem['sequence'])} paths; trapped {rem['trapped_pieces']}. Paths are FCL-swept, not endpoint-only.")
    preview=last/'bambu/native_preview.png'
    if preview.exists():shutil.copy2(preview,qc/'12_final_bambu_preview.png')
    else:raise FileNotFoundError('ACTUAL_BAMBU_PREVIEW_MISSING: '+str(preview))
    compare_actual(first,last,out/'actual_coverage_before_after.png')


def metric(r):
    s=r['slice'];rem=r['removal'];return dict(candidate_id=r['candidate_id'],
        direct_buildplate_support_fraction=r.get('direct_buildplate_support_fraction'),
        **{name.lower()+'_coverage_percent':r['semantic_coverage'][name]['percent'] for name in ('VASCULAR','PORT','BOX_BOTTOM','BOX_WALLS')},
        **{name.lower()+'_required_regions':r['semantic_coverage'][name]['region_count'] for name in ('VASCULAR','PORT','BOX_BOTTOM','BOX_WALLS')},
        support_self_unsupported=s['actual_self_buildability']['unsupported_region_count'],
        plate_release=r['plate_release']['status'],trapped_fragments=rem['trapped_pieces'],inaccessible_cuts=rem['tool_access_failed'],
        auto_support_off=s['auto_support_off'],slicer_success=s['slicer_success'],support_only_slice=s.get('support_only_slicer_success'),
        actual_first_layer=r.get('actual_first_layer'),internal_root_statistics=r.get('internal_root_statistics'),failed_gates=r['failed_gates'])


def generate(context,records,winner,b,protection,runtime):
    base,out,inputs,profile,discovery,mold,v,box,branches,top,rows,baseline,metrics=context
    inverted=[r for r in records if r['candidate_id']!='UPRIGHT_O000']
    from .support_v3_1_diagnostics import internal_root_statistics
    for record in records:
        internal_root_statistics(record,inputs);save(Path(record['folder'])/'candidate.json',record)
    best=winner or min(inverted,key=lambda r:(len(r['failed_gates']),-r['slice']['actual_vascular_coverage']['percent']))
    public=[metric(r) for r in records];table(out/'full_candidate_metrics.csv',public);table(out/'bambu_results.csv',public)
    table(out/'direct_buildplate_support.csv',[dict(candidate_id=r['candidate_id'],fraction=r.get('direct_buildplate_support_fraction'),method='Frozen REQUIRED polygon representative rays, area weighted; coarse accessibility only') for r in records])
    table(out/'plate_release_results.csv',[dict(candidate_id=r['candidate_id'],status=r['plate_release']['status'],level=r['plate_release'].get('level'),attempt_count=len(r['plate_release']['attempts']),trapped_nonplate=r['removal']['trapped_pieces']) for r in records])
    table(out/'inverted_vs_O000.csv',[metric(next(r for r in records if r['candidate_id']=='UPRIGHT_O000')),metric(best)])
    phase_a_qc(context,records,best);phase_b_qc(base,b)
    status='READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED' if winner or b.get('winner') else 'NEEDS_ADJUSTMENT'
    a=dict(candidates=[r['orientation_id'] for r in rows],exact_inverted=next(r for r in rows if r['orientation_id']=='INV_0'),
        top_candidates=json.loads((out/'full_candidate_ids.json').read_text()),best_diagnostic=best['candidate_id'],
        direct_buildplate_support_fraction={r['candidate_id']:r['direct_buildplate_support_fraction'] for r in records},
        box_ceiling_burden={r['orientation_id']:r['box_bottom_unsupported_burden_mm2'] for r in metrics},
        actual_coverage={r['candidate_id']:r['semantic_coverage'] for r in records},
        support_self_failures={r['candidate_id']:r['slice']['actual_self_buildability']['unsupported_region_count'] for r in records},
        plate_release={r['candidate_id']:r['plate_release']['status'] for r in records},
        trapped_fragments={r['candidate_id']:r['removal']['trapped_pieces'] for r in records},
        slicer_result={r['candidate_id']:r['slice']['slicer_success'] for r in records},
        winner=winner['candidate_id'] if winner else None,status='INVERTED_STUDY_PASS' if winner else 'INVERTED_STUDY_FAILED',
        failed_gates={r['candidate_id']:r['failed_gates'] for r in records})
    summary=dict(status=status,phase_a=a,phase_b=b,source_protection=protection,output=str(base/'v3_1' if b['executed'] else out),
        production_selected=str(out/'production_selected') if winner else str(base/'v3_1/production_selected') if b.get('winner') else None,
        scale_factor=1.,actual_profile=profile,runtime_seconds=runtime,printer_commands_sent=False)
    dest=base/'v3_1' if b['executed'] else out;save(dest/'support_orientation_and_closure_summary.json',summary)
    exact=next(r for r in records if r['candidate_id']=='INV_0');ceiling=next(r for r in metrics if r['orientation_id']=='INV_0')['box_bottom_unsupported_burden_mm2']
    answers=[
        f"1. 确实测试了精确倒置 INV_0。矩阵计算的开口向量为 {a['exact_inverted']['opening_normal_print']}，与打印向下方向夹角 {a['exact_inverted']['angle_to_print_down_deg']:.8f}°。",
        '2. 已生成 casting-local X/Y 两轴的 ±15°、±30°、±45°，共 13 个倒置姿态及 1 个正放对照；按规定筛选 Top 6，并保证 INV_0 完整建模与切片。',
        f"3. 当前最优倒置诊断姿态为 {best['candidate_id']}；正式胜出者为 {a['winner'] or '无'}。未通过全部门槛的排名仅供诊断。",
        f"4. 该候选实际与内腔相交的支撑实体共 {best['internal_root_statistics']['internal_components']} 个，其中 {best['internal_root_statistics']['internal_plate_components']} 个直接落在打印板。另按冻结需支撑层区域面积加权的直接落板机会为 {100*best['direct_buildplate_support_fraction']:.2f}%；后者是代表点射线机会，不是已成功覆盖比例。",
        f"5. INV_0 原盒底形成约 {ceiling:.2f} mm² 的朝下顶棚粗估面积；实际盒底覆盖率 {exact['semantic_coverage']['BOX_BOTTOM']['percent']:.4f}%，盒壁覆盖率 {exact['semantic_coverage']['BOX_WALLS']['percent']:.4f}%。",
        f"6. 最优诊断姿态的完整模型抬离结果为 {best['plate_release']['status']}。所有倒置候选均检查固定支撑与血管、端口及五面盒体的连续扫掠；非落板支撑另行检查切割与取出。",
        f"7. 倒置组状态为 {a['status']}；所有实际切片候选的逐项失败原因已写入 full_candidate_metrics.csv。"]
    if not winner:
        answers += ['8. 倒置未通过，是因为实际层覆盖、支撑自身逐层承托或拆除等软件门槛仍存在缺口；没有以射线可见率或名义 STL 覆盖代替刀路证据。',
                    '9. 精确倒置增加大面积盒底顶棚负担；45°倾斜降低此粗估负担，但盒壁、血管/端口残差和拆除风险并未同步消失。完整分类见四类 coverage 列及 plate-release CSV。',
                    '10. 因此自动回到 V3 的 O000，继续已有分段薄肋架构；未继续搜索新的姿态。']
    if b['executed']:
        answers += [f"11. 本次重新切片基线有 {b['residual_regions_before']} 个未完整覆盖的层区域，实际覆盖率 {b['coverage_before']:.6f}%。",
            f"12. 实际执行 {b['iteration_count']} 轮局部修复并重新切片，停止原因是 {b['stop_reason']}。",
            f"13. 最终血管与端口合计实际覆盖率为 {b['coverage_after']:.6f}%；四类解剖/盒体指标分别保存在每轮 actual_coverage.csv。",
            f"14. 仍有 {b['residual_regions_after']} 个未完整覆盖层区域。红色残差由实际 Bambu 挤出足迹差集得到，见 actual_coverage_before_after.png。",
            f"15. 支撑自身实际欠支撑区域为 {b['support_self_after']} 个，基线为 {b['support_self_before']} 个。",
            f"16. 受困碎片由 {b['trapped_before']} 个变为 {b['trapped_after']} 个；每个基线失败碎片均记录了肋缩短、分段、窗口、切口和根域修改的尝试及拒绝原因。",
            f"17. 不可达切割点由 {b['inaccessible_cut_before']} 个变为 {b['inaccessible_cut_after']} 个。切口上移须同时通过真实薄颈截面、3 mm 工具通道和整块碎片取出检查。",
            '18. 本轮实际 X1 Carbon 切片使用保存的 0.2 mm 喷嘴、ABS Vessel、0.10mm Standard Vessel 配置，Auto Support OFF；实际原生支撑路径另行核对为零。',
            f"19. 最终真实 X1C 本地切片成功：{b['slicer_result']}。3MF、实际 G-code、命令及 stdout/stderr 随每轮保存；从未上传、发送或启动打印。",
            '20. '+('可以进入人工审阅，但物理参数仍未校准。' if status.startswith('READY') else '目前不建议进入完整模型实物试打；NEEDS_ADJUSTMENT，保留残差失败地图，按迭代上限/停滞规则停止继续补料。')]
    report='# BG001 RMCA 倒置研究与刀路闭环修复 V3.1\n\n'+'\n\n'.join(answers)+'\n\n## 方法、边界与可复核性\n\n'
    report+='冻结血管、半径、端口和 2 mm 盒体，只作单位缩放系数为 1 的刚体摆放；修改对象仅为显式支撑。V1/V2/V3 历史结果均保留。本轮两个新输出目录分别记录倒置实验与 O000 闭环。\n\n'
    report+='覆盖依据为 R[n] 与实际 Bambu custom-support 挤出足迹 C[n] 的 GEOS 差集 U[n]。残差按结构类别、SWC 分支、弧长、相邻层、XY 接近及切向一致性聚类。宽度改变量采用真实线宽的整数倍；局部候选先做碰撞、有限切面、工具和碎片通路筛查，接受后才重新切片。保守筛查失败表示当前有限搜索未找到可靠通路，不是数学上的不可拆除证明。\n\n'
    localization=base/'v3_1/validation/baseline_fragment_collision_localization.json'
    if localization.exists():
        for row in json.loads(localization.read_text()):
            hits=[x for x in row['collision_localization'] if 'FROZEN_VASCULAR_AND_PORTS' in x['colliders']]
            if hits:report+=f"基线碎片 {row['fragment_id']}：向上移动 {hits[0]['upward_mm']:.2f} mm 的定位样本已与冻结血管/端口碰撞。完整有限路径是否通过仍以 FCL 扫掠结果为准。\n\n"
    if (base/'v3_1/history/export_guard_failed_iteration_01').exists():
        report+='网格导出调试时曾出现因新接缝精度失败而连带移除旧肋的错误；该次真实切片完整保存在 history/export_guard_failed_iteration_01，未计入正式闭环轮次。现已改为只撤回导致导出失败的新修改，并在每轮核对基线肋标识集合。修正后的第一轮重审此前接受的局部修改及全部支撑自身残差；此前拒绝的替代方案仍留存，后续轮次重新计算全部实际残差。\n\n'
    report+='数值面积容差取 min(历史 1e-5 mm² 上限，4 × float32 坐标 ULP × 实际线宽)，没有为获得完整覆盖而放宽；名义自打印与实际挤出足迹自打印分别保留。首层接触面积、岛数和 brim 几何连通仅是几何代理，不代表 ABS 附着强度通过验证。\n\n'
    report+='二维差集复用 [Shapely/GEOS 官方接口](https://shapely.readthedocs.io/en/stable/reference/shapely.difference.html)，未重写多边形布尔算法。依赖与许可证见 third_party/support_sources.json。\n\n'
    report+=f"源保护：{protection['checked_files']} 个文件，SHA256 差异 {len(protection['changed'])} 个。屋脊间隙、短桥上限、断裂力及 ABS 肋强度仍分别标记为未校准；没有物理 coupon 不是单独的软件失败理由。\n"
    (dest/'orientation_then_v3_1_report.md').write_text(report)
    (out/'inverted_study_report.md').write_text('# 开口倒置实验\n\n'+'\n\n'.join(answers[:10] if not winner else answers)+'\n\n完整后续报告：../v3_1/orientation_then_v3_1_report.md\n')
    return summary
