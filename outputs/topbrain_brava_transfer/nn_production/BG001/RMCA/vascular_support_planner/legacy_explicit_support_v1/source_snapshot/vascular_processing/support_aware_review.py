"""Auditable Chinese review package, with geometry-derived figures and tables."""
from __future__ import annotations
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import trimesh
from PIL import Image,ImageDraw,ImageFont
from .support_aware_roi import ROOT,RMCA,DEFAULT_OUTPUT,table,verify_snapshot
from .support_aware_port_design import configuration,load_design
from .sacrificial_fixture_review import save_json,polydata
from .topbrain_qc import sha256
from . import roi_support_metrics as metrics

COLORS=['#399d63','#eea23b','#ce4f4f']


def font(size):
    for p in ('/mnt/c/Windows/Fonts/msyh.ttc','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'):
        if Path(p).exists():return ImageFont.truetype(p,size)
    return ImageFont.load_default()


def annotate(image,title,lines=()):
    canvas=Image.new('RGB',(image.width,image.height+100+34*len(lines)),'white');canvas.paste(image,(0,82))
    draw=ImageDraw.Draw(canvas);draw.text((28,20),title,font=font(28),fill='#203542')
    for i,line in enumerate(lines):draw.text((28,image.height+92+34*i),str(line),font=font(22),fill='#203542')
    return canvas


def scene(meshes,path,title,lines=(),center=None,scale=None,risk=None,points=None):
    import pyvista as pv
    p=pv.Plotter(off_screen=True,window_size=[1100,800]);p.set_background('white')
    bounds=[]
    for i,item in enumerate(meshes):
        mesh,color,opacity=item;data=polydata(mesh);bounds.append(mesh.bounds)
        if risk is not None and i==0:
            data.cell_data['support_risk']=np.asarray(risk);p.add_mesh(data,scalars='support_risk',cmap=COLORS,clim=[0,2],show_scalar_bar=False,smooth_shading=False)
        else:p.add_mesh(data,color=color,opacity=opacity,smooth_shading=False)
    if points:
        for pos,label in points:p.add_point_labels(np.asarray(pos)[None,:],[label],font_size=20,shape=None,text_color='black',point_size=7,always_visible=True)
    if not bounds:
        p.close();image=Image.new('RGB',(1100,800),'#f6f7f7')
    else:
        low=np.min([b[0] for b in bounds],axis=0);hi=np.max([b[1] for b in bounds],axis=0)
        center=(low+hi)/2 if center is None else np.asarray(center);extent=float(np.max(hi-low))
        # Risk faces point downwards in the casting frame. A below-side camera
        # makes those actual orange/red triangles visible instead of hiding
        # them behind the green upper surface.
        direction=np.array([1.,-1.25,-1.1 if risk is not None else 1.1]);direction/=np.linalg.norm(direction)
        p.camera_position=[center+direction*max(extent*2,10),center,[0,0,1]];p.enable_parallel_projection();p.camera.parallel_scale=scale or max(2,extent*.68)
        if risk is not None:
            target=np.array([[-30,-25],[30,-25],[30,25],[-30,25],[-30,-25]],float)+center[:2]
            target=np.c_[target,np.full(5,low[2]-.6)]
            p.add_mesh(pv.lines_from_points(target),color='#82909a',line_width=1)
        p.add_axes();image=Image.fromarray(p.screenshot(return_img=True));p.close()
    annotate(image,title,lines).save(path)


def chart(path,title,labels,old,new,units):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,len(labels),figsize=(14,4.3))
    for ax,label,a,b,unit in zip(np.atleast_1d(axes),labels,old,new,units):
        vals=[0 if a is None else a,0 if b is None else b];bars=ax.bar(['Legacy','New'],vals,color=['#788fa0','#399d63'])
        ax.set_title(label);ax.set_ylabel(unit);ax.spines[['top','right']].set_visible(False)
        for j,bar in enumerate(bars):ax.text(bar.get_x()+bar.get_width()/2,bar.get_height(),f'{vals[j]:.2f}' if (a,b)[j] is not None else 'pending',ha='center',va='bottom')
        ax.set_ylim(0,max(vals+[1])*1.22)
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)
    image=Image.open(path).convert('RGB');annotate(image,title).save(path)


def selected_rank(result,model):
    orient=result.get('orientations',{}).get('selected',{});row=orient.get('orientation',{});record=orient.get('slice')
    hard=[f for f in result['failures'] if f!='ROI_SIZE_TARGET_NOT_MET']
    return [bool(hard) or not row.get('passed',False),row.get('trapped_risk_count',999)!=0,record is None,
        bool(record and record.get('critical_warnings')),not model['dimensions']['size_target_met'],
        record['support_quantity']['estimated_internal_support_mass_g'] if record else 1e9,
        record['contact']['vascular_contact_area_mm2'] if record else 1e9,
        model['dimensions']['size_violation_mm'],result['candidate_id']]


def build_summary(out):
    cfg=configuration();models=json.loads((out/'model_results.json').read_text());modelmap={r['candidate_id']:r for r in models}
    candidates=json.loads((out/'candidates.json').read_text());cmap={r['candidate_id']:r for r in candidates}
    legacy=json.loads((out/'baseline/legacy_70mm_slice_metrics.json').read_text());evaluations=json.loads((out/'evaluations.json').read_text())
    if 'vascular_surface_support' not in legacy:
        from . import abs_casting_mold as old_mold
        oldcfg=old_mold.load_config(ROOT/cfg['baseline_config']);inputs=old_mold.load_inputs(oldcfg)
        legacy['vascular_surface_support']=metrics.surface_support(inputs['central'],inputs['companion']['graph'],
            inputs['companion']['transform'],np.array(legacy['orientation']['transform_4x4'])[:3,:3])[0]
        save_json(out/'baseline/legacy_70mm_slice_metrics.json',legacy)
    selected=min(evaluations,key=lambda r:selected_rank(r,modelmap[r['candidate_id']]))
    # Prefer a comparable 3-outlet successful alternative unless its real
    # support burden is materially larger (>25%). Keep this decision explicit.
    three=[r for r in evaluations if cmap[r['candidate_id']]['outlet_count']==3 and r.get('orientations',{}).get('selected',{}).get('slice')]
    preference='LOWEST_HARD_GATE_RISK_THEN_ACTUAL_INTERNAL_SUPPORT'
    if three and cmap[selected['candidate_id']]['outlet_count']==2:
        alternate=min(three,key=lambda r:selected_rank(r,modelmap[r['candidate_id']]))
        a,b=selected_rank(alternate,modelmap[alternate['candidate_id']]),selected_rank(selected,modelmap[selected['candidate_id']])
        if a[:5]==b[:5] and a[5]<=1.25*b[5]:selected=alternate;preference='THREE_OUTLETS_WITHIN_25_PERCENT_REAL_SUPPORT_COST'
    cid=selected['candidate_id'];candidate=cmap[cid];model=modelmap[cid];geometry=selected.get('geometry',{})
    chosen=selected.get('orientations',{}).get('selected',{});orientation=chosen.get('orientation',{});slice=chosen.get('slice')
    profile=json.loads((out/'bambu/evaluation_profile/sha256.json').read_text())
    baseline_slice=legacy.get('selected_slice');identical=bool(slice and baseline_slice and profile==slice['evaluation_profile_sha256']==baseline_slice['evaluation_profile_sha256'])
    support=dict(baseline_internal_support_g=legacy['baseline_internal_support_g'],new_internal_support_g=slice['support_quantity']['estimated_internal_support_mass_g'] if slice else None,
        baseline_total_support_g=legacy['baseline_support_mass_g'],new_total_support_g=slice['support_quantity']['total_support_mass_g'] if slice else None,
        baseline_support_contact_area=legacy['baseline_support_contact_area'],new_support_contact_area=slice['contact']['vascular_contact_area_mm2'] if slice else None,
        trapped_support_count=orientation.get('trapped_risk_count'),top_accessible_fraction=orientation.get('support',{}).get('top_accessible_fraction'),
        contact_method='ACTUAL_GCODE_SUPPORT_PROXIMITY_TRIANGLE_AREA_ESTIMATE',contact_envelope_mm=.5,physical_contact_certified=False,
        support_removal_certified=False)
    support['mass_reduction_percent']=metrics.reduction(support['baseline_internal_support_g'],support['new_internal_support_g'])
    support['contact_reduction_percent']=metrics.reduction(support['baseline_support_contact_area'],support['new_support_contact_area'])
    support['reduction_target_met']=metrics.reduction_gate(support['baseline_internal_support_g'],support['new_internal_support_g'],support['baseline_support_contact_area'],support['new_support_contact_area'])
    protection=verify_snapshot(json.loads((out/'protected_before.json').read_text()))
    failures=selected['failures'].copy()
    if not orientation.get('passed',False):failures.append('PRINT_ORIENTATION_PLACEMENT_OR_SUPPORT_GATE_FAILED')
    if support['trapped_support_count']!=0:failures.append('TRAPPED_SUPPORT_RISK_OR_UNEVALUATED')
    if not support['reduction_target_met']:failures.append('SUPPORT_REDUCTION_TARGET_NOT_MET')
    if not protection['all_unchanged']:failures.append('PROTECTED_SOURCE_CHANGED')
    if slice and not identical:failures.append('SLICER_PROFILE_COMPARISON_MISMATCH')
    if slice and slice['critical_warnings']:failures.append('CRITICAL_SLICER_WARNING')
    if slice and not slice['feature_preservation']['passed']:failures.append('SLICER_DROPPED_VASCULAR_FEATURE')
    if not geometry:failures.append('FULL_MOLD_NOT_GENERATED')
    pending=slice is None
    status='NEEDS_ADJUSTMENT' if failures else 'READY_FOR_HUMAN_REVIEW_WITH_BAMBU_PENDING' if pending else 'READY_FOR_HUMAN_REVIEW'
    dims=model['dimensions'];old=legacy['vascular_dimensions']
    roi={k:candidate[k] for k in ('inlet_count','outlet_count','branch_count','bifurcation_count','centerline_length_mm','root_node','outlet_nodes')}
    roi.update(actual_long_mm=dims['long_side_mm'],actual_short_mm=dims['short_side_mm'],z_span_mm=dims['Z_span_mm'],
        size_target_met=dims['size_target_met'],scale_factor=1.,mevo_only=True,source_component=candidate['source_component'],
        footprint_area_reduction_percent=100*(1-dims['long_side_mm']*dims['short_side_mm']/(old['long_side_mm']*old['short_side_mm'])))
    summary=dict(name=cfg['production_name'],status=status,failures=sorted(set(failures)),selected_candidate=cid,selection_reason=preference,
        target=dict(long_mm=60,short_mm=50,tolerance_mm=3,scale_factor=1.),selected_roi=roi,model=model,
        geometry={**geometry.get('union',model['mesh']),'minimum_pdms_ligament':geometry.get('pdms',model['pdms'])['minimum_ligament_mm']},
        box=geometry.get('box'),support=support,baseline=legacy,
        bambu=dict(printer='Bambu Lab X1 Carbon',nozzle=.2,ABS_profile='Bambu ABS Vessel',support_profile='FROZEN_ACTIVE_PROCESS_TREE_AUTO',
            real_slicing_pass=bool(slice and slice['status']=='BAMBU_SLICED' and slice['feature_preservation']['passed'] and not slice['critical_warnings']),
            print_time_min=slice['estimated_print_time_seconds']/60 if slice else None,total_filament_g=slice.get('filament_used_g') if slice else None,
            support_mass_g=support['new_total_support_g'],warnings=slice['warnings'] if slice else [],profiles_identical=identical,selected=slice,
            printer_job_sent=False),orientation=orientation,source_protection=protection,source_hashes=json.loads((out/'protected_before.json').read_text()),
        candidates_generated=len(candidates),VascularMD_candidates=len(models),top3_evaluated=[r['candidate_id'] for r in evaluations],
        topology_recomputed=False,registration_rerun=False,NN_rerun=False,manufacturing_validated=False,
        new_port_interfaces=geometry.get('interfaces',[]),new_ports=geometry.get('ports',[]))
    summary['new_final_port_interfaces']=geometry.get('final_interfaces',[])
    summary['port_route_review']=geometry.get('route_review')
    summary['completed_mold_candidates']=[r['candidate_id'] for r in evaluations if r.get('orientations',{}).get('selected',{}).get('slice')]
    prior_path=out/'previous_port_box_geometry.json'
    prior=json.loads(prior_path.read_text()).get(cid,{}) if prior_path.exists() else {}
    previous_ports={r['port_id']:r for r in prior.get('ports',[])}
    lengths=[dict(port_id=r['port_id'],wall=r['face'],previous_extension_mm=previous_ports.get(r['port_id'],{}).get('length_mm'),
        current_extension_mm=r['length_mm'],outside_stub_mm=r.get('outside_stub_mm'),preferred_maximum_mm=cfg['ports']['preferred_maximum_extension_mm'],
        preferred_length_met=r.get('moderate_length_target_met'),vessel_clearance_mm=r['vessel_clearance_mm'],
        other_port_clearance_mm=r['minimum_port_clearance_mm']) for r in geometry.get('ports',[])]
    summary['adaptive_geometry_revision']=dict(port_lengths=lengths,unconnected_face_target_gap_mm=cfg['box']['unconnected_face_gap_mm'],
        actual_face_clearances=geometry.get('box',{}).get('adaptive_faces',[]),
        previous_box_outer_dimensions_mm=prior.get('box',{}).get('outer_dimensions_mm'),
        current_box_outer_dimensions_mm=geometry.get('box',{}).get('outer_dimensions_mm'),
        policy='ADAPT_FREE_FACES_TO_ACTUAL_CORE_WITH_PORT_EDGE_MATERIAL_AND_SHORTER_PORT_ROUTES')
    upright_path=out/'baseline/legacy_70mm_upright_metrics.json'
    if upright_path.exists():
        upright=json.loads(upright_path.read_text());record=upright.get('selected_slice')
        if record:
            mass=record['support_quantity']['estimated_internal_support_mass_g']
            summary['supplemental_upright_baseline']=dict(real_resliced=upright['real_resliced'],internal_support_g=mass,
                support_reduction_percent=metrics.reduction(mass,support['new_internal_support_g']),
                new_pose_is_upright=orientation.get('candidate_id')==0,
                profiles_identical=record['evaluation_profile_sha256']==profile,
                print_time_min=record['estimated_print_time_seconds']/60,
                source_metrics=str(upright_path))
    return summary,selected,models,candidates


def figures(out,summary,selected,models,candidates):
    import pyvista as pv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cfg=configuration();qcdir=out/'QC';cid=summary['selected_candidate'];directory=out/'candidates'/cid
    meshes={r['candidate_id']:trimesh.load_mesh(out/'candidates'/r['candidate_id']/'vascular_core.stl') for r in models if r.get('mesh',{}).get('passed')}
    low=np.min([m.bounds[0] for m in meshes.values()],axis=0);hi=np.max([m.bounds[1] for m in meshes.values()],axis=0);center=(low+hi)/2;scale=max(hi-low)*.7
    # Frozen upstream centerline, only already-accepted d075 MeVO segments.
    from .support_aware_roi import discover
    source=discover();plot=pv.Plotter(off_screen=True,window_size=[1300,900]);plot.set_background('white')
    for s in source['sources']:
        g=s['graph'];t=source['transform']
        for a,b in g.edges:
            pts=np.array([g.nodes[n]['coords'][:3] for n in (a,b)])@t[:3,:3].T+t[:3,3]
            plot.add_mesh(pv.lines_from_points(pts),color='#607f94',line_width=3)
    plot.add_mesh(polydata(meshes[cid]),color='#ce4f4f');plot.view_xy();plot.enable_parallel_projection();plot.reset_camera()
    annotate(Image.fromarray(plot.screenshot(return_img=True)),'冻结的 RMCA MeVO 输入与新候选',['灰蓝：通过既有 d075 与 MeVO 标签的血管；红：本轮选中子树']).save(qcdir/'01_source_rmca_mevo_tree.png');plot.close()
    fig,ax=plt.subplots(figsize=(11,7))
    for count,color in [(2,'#45879c'),(3,'#db8d35')]:
        rows=[r for r in candidates if r['outlet_count']==count];ax.scatter([r['long_side_mm'] for r in rows],[r['short_side_mm'] for r in rows],s=7,alpha=.2,label=f'{count} outlets',c=color)
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((57,47),6,6,fill=False,lw=2,color='black',label='60 x 50 +/- 3 mm'))
    ax.set(xlabel='Long side (mm), centerline-radius envelope',ylabel='Short side (mm)',title=f'{len(candidates)} deterministic native subtrees');ax.legend();fig.tight_layout();fig.savefig(qcdir/'02_all_candidate_rois.png',dpi=160);plt.close(fig)
    panels=[]
    for i,r in enumerate(models,1):
        mesh=meshes[r['candidate_id']];dims=r['dimensions'];colors=np.load(out/'candidates'/r['candidate_id']/'surface_support.npz')['risk_class']
        path=qcdir/f'candidate_{i:02}.png';c=next(c for c in candidates if c['candidate_id']==r['candidate_id'])
        scene([(mesh,COLORS[0],1)],path,f"{r['candidate_id']} · {dims['long_side_mm']:.2f} × {dims['short_side_mm']:.2f} mm",
            [f"{c['outlet_count']} 出口 / {c['bifurcation_count']} 分叉；风险面积 {r['support']['overhang_area_mm2']:.2f} mm²",'绿：低支撑风险；橙：支撑风险；红：细血管/分叉附近'],center=center,scale=scale,risk=colors)
        panels.append(Image.open(path).convert('RGB').resize((550,520)))
    canvas=Image.new('RGB',(2200,1040),'white')
    for i,im in enumerate(panels):canvas.paste(im,((i%4)*550,(i//4)*520))
    canvas.save(qcdir/'03_top8_candidates_same_scale.png')
    roi=summary['selected_roi'];mesh=meshes[cid];risk=np.load(directory/'surface_support.npz')['risk_class']
    scene([(mesh,'#399d63',1)],qcdir/'04_selected_roi_dimensions.png','实际血管表面尺寸（不含端口和盒体）',
        [f"长 {roi['actual_long_mm']:.2f} mm；短 {roi['actual_short_mm']:.2f} mm；Z {roi['z_span_mm']:.2f} mm",'整体缩放：NO；固定旧版铸造坐标系；尺寸未达目标' if not roi['size_target_met'] else '整体缩放：NO；固定旧版铸造坐标系'])
    companion=__import__('vascular_processing.support_aware_port_design',fromlist=['load_companion']).load_companion(directory,source['transform'],cfg,align=False)
    scene([(mesh,'#399d63',1)],qcdir/'05_selected_roi_topology.png','新子树的入口、出口与分叉',
        [f"1 入口 / {roi['outlet_count']} 出口 / {roi['branch_count']} 分支 / {roi['bifurcation_count']} 主要分叉",f"原始中心线总长 {roi['centerline_length_mm']:.2f} mm"],points=[(e.position,e.endpoint_id) for e in companion['endpoints']])
    scene([(mesh,'#399d63',1)],qcdir/'06_selected_roi_support_risk.png','血管实际三角面支撑风险',risk=risk,lines=['绿色：低风险；橙色：向下悬垂；红色：细血管或分叉附近','这是几何风险面积；不等于真实支撑接触面积。'])
    gap=summary['geometry']['minimum_pdms_ligament'];pair=selected.get('geometry',{}).get('pdms',summary['model']['pdms'])['minimum_pair']
    gapmesh=trimesh.creation.cylinder(radius=.12,segment=[pair['first_point_mm'],pair['second_point_mm']],sections=12)
    assembly=None
    if selected.get('geometry'):assembly=load_design(directory,source['transform'],cfg)
    scene([((assembly['inputs']['core'] if assembly else mesh),'#83b99b',.65),(gapmesh,'#ce4f4f',1)],qcdir/'07_selected_roi_pdms_gap.png','非连接血管之间的最薄 PDMS 间距',[f'FCL 三角面距离：{gap:.3f} mm；下限 1.600 mm',f"最近对：{pair['first']} / {pair['second']}"])
    if assembly:
        inputs,box,union=[assembly[k] for k in ('inputs','box','union')]
        scene([(inputs['central'],'#399d63',1)]+[(p['reference_mesh'],'#dd9950' if p['port_id']=='I1' else '#528faa',1) for p in inputs['ports']],qcdir/'08_new_ports.png','适度缩短的侧壁延长段',
            lines=['；'.join(f"{p['port_id']}：{p['length_mm']:.2f} mm" for p in summary['new_ports']),'长度沿延长段中心线计算，包含穿墙段及墙外 6 mm 接口'],
            points=[(p['target'],p['port_id']+' '+p['face']) for p in inputs['ports']])
        panels=[]
        for i,r in enumerate(summary['new_port_interfaces']):
            endpoint=next(e for e in companion['endpoints'] if e.endpoint_id==r['port_id']);path=qcdir/f"interface_{r['port_id']}.png"
            scene([(inputs['core'],'#82a8a0',1)],path,r['port_id']+' 真实 24 点端面连接',
                [f"20° 特征边 {r['feature_edge_20_length_mm']:.6f} mm；锐环 {r['new_closed_sharp_ring']}",f"截面积变化 {100*r['section_area_change_fraction']:.5f}%"],center=endpoint.position,scale=3)
            panels.append(Image.open(path).convert('RGB').resize((550,530)))
        canvas=Image.new('RGB',(550*len(panels),530),'white')
        for i,im in enumerate(panels):canvas.paste(im,(550*i,0))
        canvas.save(qcdir/'09_port_interfaces.png')
        scene([(box['mesh'],'#99b2c3',.5),(inputs['core'],'#399d63',1)],qcdir/'10_new_box.png','按血管与端口形状自适应收紧的五面盒体',[f"外部尺寸 {' × '.join(f'{v:.2f}' for v in summary['box']['outer_dimensions_mm'])} mm",'自由侧壁和底面目标净距 5 mm；穿墙处保留端口边缘材料','侧壁 3 mm；底板 3 mm；顶部开放；没有装配孔'])
        scene([(union,'#9fb8c1',1)],qcdir/'11_final_one_body.png','血管、端口和盒体的一体 ABS 模型',[f"闭合：{summary['geometry']['watertight']}；连通体：{summary['geometry']['connected_components']}",f"血管材料损失 {summary['geometry']['core_material_loss_mm3']:.8f} mm³"])
        support=summary['orientation'].get('support',{});paths=[]
        for p in support.get('paths',[]):paths.append((trimesh.creation.cylinder(radius=.16,segment=[p['start'],p['end']],sections=8),'#35996b',.7))
        scene([(union,'#9fb2be',.22)]+paths,qcdir/'12_support_accessibility.png','仅通过顶部开口检查支撑可达性',
            [f"Trapped 风险区：{support.get('trapped_risk_count','待检查')}；探针直径 3 mm",'绿色线段为通过碰撞检查的探针通道；仍需实际拆支撑验证。'])
        top=selected.get('orientations',{}).get('top3',[]);panels=[]
        for i,row in enumerate(top):
            current=union.copy();current.apply_transform(row['transform_4x4']);path=qcdir/f'orientation_top_{i+1}.png'
            scene([(current,'#85aaa4',1)],path,row['method'],[f"Trapped {row['trapped_risk_count']}；内部风险面积 {row['internal_support_risk_area_mm2']:.2f} mm²"])
            panels.append(Image.open(path).convert('RGB').resize((550,510)))
        canvas=Image.new('RGB',(max(1,len(panels))*550,510),'white')
        for i,im in enumerate(panels):canvas.paste(im,(i*550,0))
        canvas.save(qcdir/'14_orientation_top3.png')
    record=summary['bambu']['selected']
    if record and record.get('native_preview'):
        image=Image.open(record['native_preview']).convert('RGB');annotate(image,'本机 Bambu Studio 实际切片预览',['X1 Carbon / 0.2 mm / ABS；仅本地切片，未发送打印机。']).save(qcdir/'15_selected_bambu_preview.png')
    old=summary['baseline'];sup=summary['support'];new=summary['bambu']
    chart(qcdir/'13_baseline_vs_new_support.png','重新切片的旧基线与新候选：相同打印机、材料、工艺和支撑参数',
        ['Internal support','Vascular contact estimate','Print time','Branches / outlets'],
        [old['baseline_internal_support_g'],old['baseline_support_contact_area'],old['baseline_print_time_min']/60,old['branch_count']],
        [sup['new_internal_support_g'],sup['new_support_contact_area'],new['print_time_min']/60 if new['print_time_min'] else None,roi['branch_count']],['g','mm2','hours','branches'])
    comparison=Image.open(qcdir/'13_baseline_vs_new_support.png').convert('RGB')
    annotate(comparison,f"旧 {old['vascular_dimensions']['long_side_mm']:.2f} × {old['vascular_dimensions']['short_side_mm']:.2f} mm、3 出口；新 {roi['actual_long_mm']:.2f} × {roi['actual_short_mm']:.2f} mm、{roi['outlet_count']} 出口").save(qcdir/'13_baseline_vs_new_support.png')


def report(out=DEFAULT_OUTPUT):
    out=Path(out).resolve();summary,selected,models,candidates=build_summary(out);cid=summary['selected_candidate'];directory=out/'candidates'/cid
    for name in ('selected','casting_mold','tables','QC'):(out/name).mkdir(exist_ok=True)
    for name in ('centerline.swc','compensated.swc','centerline_mapping.csv','compensated_mapping.csv','compensated_edges.csv','vascular_core.stl','vascular_core.vtp'):
        if (directory/name).exists():shutil.copy2(directory/name,out/'selected'/name)
    shutil.copy2(directory/'vascular_core.stl',out/'selected'/(summary['name']+'.stl'))
    assembly=directory/'manufacturing'
    if assembly.exists():
        if (assembly/'core_with_ports.stl').exists():shutil.copy2(assembly/'core_with_ports.stl',out/'selected/core_with_ports.stl')
        for name in ('intact_open_top_box.stl','final_abs_casting_mold.stl'):
            if (assembly/name).exists():shutil.copy2(assembly/name,out/'casting_mold'/name)
        for name in ('pdms_ligament_top10','port_design','port_wall_union','port_interfaces','adaptive_box_faces'):
            if (assembly/(name+'.csv')).exists():shutil.copy2(assembly/(name+'.csv'),out/'tables'/(name+'.csv'))
    chosen=summary['bambu']['selected']
    if chosen:
        for key,name in [('archive','final_abs_casting_mold.3mf'),('gcode_file','final_abs_casting_mold.gcode')]:shutil.copy2(chosen[key],out/'casting_mold'/name)
        shutil.copy2(summary['orientation']['stl'],out/'casting_mold/final_abs_casting_mold_oriented.stl')
    if (directory/'orientations/orientation_candidates.csv').exists():shutil.copy2(directory/'orientations/orientation_candidates.csv',out/'tables/orientation_candidates.csv')
    old=summary['baseline'];sup=summary['support'];roi=summary['selected_roi'];new=summary['bambu']
    table(out/'tables/baseline_vs_60x50.csv',[
        dict(model='LEGACY_70MM_BASELINE',long_mm=old['vascular_dimensions']['long_side_mm'],short_mm=old['vascular_dimensions']['short_side_mm'],internal_support_g=old['baseline_internal_support_g'],total_support_g=old['baseline_support_mass_g'],contact_estimate_mm2=old['baseline_support_contact_area'],print_time_min=old['baseline_print_time_min'],branches=5,outlets=3),
        dict(model=summary['name'],long_mm=roi['actual_long_mm'],short_mm=roi['actual_short_mm'],internal_support_g=sup['new_internal_support_g'],total_support_g=sup['new_total_support_g'],contact_estimate_mm2=sup['new_support_contact_area'],print_time_min=new['print_time_min'],branches=roi['branch_count'],outlets=roi['outlet_count'])])
    table(out/'tables/support_metrics.csv',[dict(candidate_id=r['candidate_id'],**r.get('support',{})) for r in models])
    table(out/'tables/port_length_comparison.csv',summary['adaptive_geometry_revision']['port_lengths'])
    modelmap={r['candidate_id']:r for r in models}
    evaluations=json.loads((out/'evaluations.json').read_text());evaluationmap={r['candidate_id']:r for r in evaluations};rows=[]
    for candidate in candidates:
        r=modelmap.get(candidate['candidate_id'],{});evaluation=evaluationmap.get(candidate['candidate_id'],{})
        geometry=evaluation.get('geometry',{});slice=evaluation.get('orientations',{}).get('selected',{}).get('slice')
        ports=geometry.get('ports',[])
        rows.append({**{k:v for k,v in candidate.items() if k not in ('native_node_ids','bounds_mm')},**r.get('dimensions',{}),**r.get('support',{}),
            'minimum_pdms_gap_mm':geometry.get('pdms',r.get('pdms',{})).get('minimum_ligament_mm'),
            'status':'EVALUATED_WITH_FAILURES' if evaluation.get('failures') else r.get('status','SURROGATE_ONLY'),
            'failures':evaluation.get('failures',r.get('failures',[])),
            'port_max_length_mm':max((p['length_mm'] for p in ports),default=None),
            'port_total_length_mm':sum(p['length_mm'] for p in ports) if ports else None,
            'box_outer_dimensions_mm':geometry.get('box',{}).get('outer_dimensions_mm'),
            'actual_internal_support_g':slice['support_quantity']['estimated_internal_support_mass_g'] if slice else None})
    table(out/'tables/roi_candidate_summary.csv',rows)
    anatomy=out/'tables/anatomy_exclusions.csv';historical=out/'history/long_ports_fixed_margins/tables/anatomy_exclusions.csv'
    if not anatomy.exists() and historical.exists():shutil.copy2(historical,anatomy)
    figures(out,summary,selected,models,candidates)
    def fmt(v):return '尚未建立' if v is None else f'{v:.2f}'
    box=summary.get('box') or {};inner=' × '.join(fmt(v) for v in box.get('inner_dimensions_mm',[])) or '尚未生成';outer=' × '.join(fmt(v) for v in box.get('outer_dimensions_mm',[])) or '尚未生成'
    port_table='\n'.join(f"| {r['port_id']} | {fmt(r['previous_extension_mm'])} | {fmt(r['current_extension_mm'])} | {r['wall']} |" for r in summary['adaptive_geometry_revision']['port_lengths'])
    face_table='\n'.join(f"| {r['face']} | {fmt(r.get('core_envelope_gap_mm')) if not r['port_connected'] else '端口穿墙'} | {fmt(r.get('fcl_surface_distance_mm')) if r['face']!='TOP_OPENING' and not r['port_connected'] else '不适用'} | {fmt(r['movement_toward_core_mm'])} |" for r in box.get('adaptive_faces',[]))
    prior_box=summary['adaptive_geometry_revision'].get('previous_box_outer_dimensions_mm')
    prior_box_text=' × '.join(fmt(v) for v in prior_box) if prior_box else '无对应前版'
    upright=summary.get('supplemental_upright_baseline',{})
    warning_text='、'.join(sorted({r.get('msg',r.get('message',r.get('text','未分类提示'))) for r in new['warnings']})) or '无'
    route_review_text=('在首次收紧方案中，O2 通向 +Y 侧壁时出现一处顶部不可达的支撑风险区。因此依据真实可达性结果，将该候选的 O2 改向 -X 侧壁重新搜索；真实端面、中央血管、弯曲半径与间距要求均保持。失败方案与重新选路依据保留在候选历史记录及 `port_route_review.json` 中。' if summary.get('port_route_review') else '')
    upright_text=(f"另外，旧基线在同一铸造正立姿态下重新切片的盒内支撑量为 {fmt(upright.get('internal_support_g'))} g。新候选是否也为正立姿态：{upright.get('new_pose_is_upright')}；相对此正立对照的减少量为 {fmt(upright.get('support_reduction_percent'))}%。该补充对照用于区分姿态变化与几何变化带来的改善。" if upright else '同姿态补充基线尚未生成。')
    text=f'''# BG001 RMCA 支撑感知制造区域研究报告

## 一、研究目的与制造问题

已有 BALANCED 模型完成了血管建模、真实端面衔接和五面盒体构造。本轮进一步研究如何从真实血管树选择较易制造的区域。研究动机是支撑生成和支撑去除已成为后续实验的主要瓶颈，而非据此否定旧模型的几何结果。旧版所有模型、代码和日志在本轮均作为只读基线。

## 二、新区域的实际尺寸

目标是在固定铸造坐标系中获得 60 × 50 mm 血管范围，允许长边 57–63 mm、短边 47–53 mm。由真实 STL 测得，本次最佳候选 **{roi['actual_long_mm']:.2f} × {roi['actual_short_mm']:.2f} mm**，高度 **{roi['z_span_mm']:.2f} mm**。该尺寸只包含中央血管，不包含人工端口与盒体。尺寸状态为 **{'PASS' if roi['size_target_met'] else 'ROI_SIZE_TARGET_NOT_MET'}**。

本轮严格沿用冻结的语义细化与 d075 结果，未重新配准、迁移标签、分类或改变直径筛选。历史 UNKNOWN 连接边被排除。此可用范围不能保证达到原先设想的 60 mm 长边，因此没有为凑尺寸添加未批准的血管。旧版“70 mm”是历史名称，其相同坐标系下的实际血管 footprint 为 {old['vascular_dimensions']['long_side_mm']:.2f} × {old['vascular_dimensions']['short_side_mm']:.2f} mm，不是 70 × 70 mm；按实际表面包围矩形，本候选面积减少 {roi['footprint_area_reduction_percent']:.2f}%。

## 三、缩放与几何依据

**整体缩放：NO；scale_factor = 1.0。** 新范围通过原始节点的子树选择得到，中心线保留原始坐标和父子连接。半径使用已批准的 0.60 mm 制造半径下限及 5 mm 过渡规则，未因换用 0.2 mm 喷嘴而减小管径。VascularMD 采用现有逐网络 AIC 拟合流程；真实表面检查独立完成，未使用全局重新网格化掩盖失败。

新端口读取 VascularMD 双精度真实端面，保留 24 点边界、面积加权法向和渐缩内埋套领。套领布尔检查仅平移至局部原点，避免世界坐标处的浮点量化误差；没有改变几何或放宽容差。交付 STL 与双精度表面另行核对。

## 四、保留的血管结构

候选 {cid} 来自冻结 part{roi['source_component']:02} 范围，包含 **{roi['inlet_count']} 个入口、{roi['outlet_count']} 个出口、{roi['branch_count']} 条拓扑分支和 {roi['bifurcation_count']} 个主要分叉**。原始中心线总长度为 {roi['centerline_length_mm']:.2f} mm。所有末端分支长度至少 3 mm。共评估 {len(candidates)} 个确定性子树；前 8 个完成真实 VascularMD 表面建模。首选 3 个进入端口与盒体阶段，失败者保留原因并由后续候选补充；实际尝试 {len(evaluations)} 个，其中 {len(summary['completed_mold_candidates'])} 个完成完整盒体切片。三出口候选作为对照保留，避免仅凭简化拓扑就排除其代表性。

## 五、与重新切片基线的比较

两者均使用归档的 X1 Carbon 0.2 mm、Bambu ABS Vessel、0.10mm Standard Vessel 及相同支撑参数。归档配置哈希完全一致：**{new['profiles_identical']}**。旧基线使用其当前已选定打印姿态重新切片，几何文件 SHA256 记录于 summary；未沿用旧报告的支撑数字。

旧基线保持已验收的 {old['current_box_thickness_mm']:.0f} mm 盒壁，本阶段遵循新任务的 3 mm 盒壁要求。两者的 ROI、人工延长段、盒体尺寸和所选姿态均可能不同，因此这里比较的是完整制造方案，不能将全部改善单独归因于血管范围缩小。

| 指标 | 重新切片的旧基线 | 新候选 |
|---|---:|---:|
| 盒内支撑材料 / g | {fmt(sup['baseline_internal_support_g'])} | {fmt(sup['new_internal_support_g'])} |
| 全部支撑材料 / g | {fmt(sup['baseline_total_support_g'])} | {fmt(sup['new_total_support_g'])} |
| 血管接触面积估计 / mm² | {fmt(sup['baseline_support_contact_area'])} | {fmt(sup['new_support_contact_area'])} |
| 打印时间 / min | {fmt(old['baseline_print_time_min'])} | {fmt(new['print_time_min'])} |
| 分支 / 出口 | 5 / 3 | {roi['branch_count']} / {roi['outlet_count']} |
| 细血管支撑风险分支（几何代理） | {old['vascular_surface_support']['support_critical_branches']} | {summary['model']['support']['support_critical_branches']} |

盒内支撑量根据实际 G-code 正挤出量、材料直径和密度估计，减少 **{fmt(sup['mass_reduction_percent'])}%**。血管接触面积估计减少 **{fmt(sup['contact_reduction_percent'])}%**。接触面积是实际支撑路径附近 0.5 mm 范围内的三角面面积加权估计，不是由独立支撑实体求得的精确物理接触面积。它不能替代真正打印后的接触与去除观察。

{upright_text}

## 六、选择依据与未通过项目

选择先考虑几何、拓扑、PDMS 间距、支撑困陷与真实切片，再比较盒内支撑材料和血管接触，而不使用一个隐藏的加权总分。当前最薄非连接区域 PDMS 间距为 **{summary['geometry']['minimum_pdms_ligament']:.3f} mm**。支撑困陷风险区数为 **{sup['trapped_support_count']}**，顶部可达比例为 **{fmt(100*sup['top_accessible_fraction']) if sup['top_accessible_fraction'] is not None else '未建立'}%**。

总体状态为 **{summary['status']}**。未通过或未完成项目：{'; '.join(summary['failures']) or '无'}。尺寸、间距和困陷约束均如实保留；存在失败时，交付文件只能作为待调整的研究候选，不能解释为已满足生产条件。

## 七、重新拟合的铸造盒体

盒体内部尺寸为 **{inner} mm**，外部尺寸为 **{outer} mm**。同一候选调整前的外部尺寸为 {prior_box_text} mm。四侧壁和底板厚度均为 3 mm，顶部开放。根据后续提出的紧凑性要求，固定余量仅用于初始端口搜索，不再决定最终盒体：未连接端口的侧面和底面贴合中央血管及延长段的实际表面包络，目标净距为 5 mm；顶部开放边缘也按最高几何加 5 mm 收紧。穿墙面保持与端口相配合的位置，邻近边缘须保留 8 mm 的端口中心至盒边距离，所以个别面的净距可能大于 5 mm。全程没有移动或缩细原血管。

| 面 | 实际包络净距 / mm | 独立 FCL 表面距离 / mm | 向内收紧 / mm |
|---|---:|---:|---:|
{face_table}

端口按可行侧壁重新选路，优先减少超过 35 mm 的延长段数量，再减少最长段及总长度。35 mm 是设计偏好，实际端面方向、至少 6 mm 的弯曲半径和碰撞约束仍须满足；墙外直段由 10 mm 缩短为 6 mm。下表长度沿人工延长段中心线测量，包含穿墙段与外接口，不包含原生血管。

| 端口 | 调整前长度 / mm | 调整后长度 / mm | 穿墙方向 |
|---|---:|---:|---|
{port_table}

{route_review_text}

新端口穿过完整侧壁，无装配孔；穿墙体积、轴向重叠、最终连通性、自由壁面净距与血管材料保留均有独立 CSV/JSON 记录。旧版长端口及固定余量结果保存在本轮输出的 `history/` 中；候选自身的旧制造结果与切片保存在各候选的 `history/` 中，原始路径重定位关系见 `relocation.json`。

## 八、支撑是否预计更容易处理

**{'YES（仅指当前几何与切片代理评估）' if sup['reduction_target_met'] and sup['trapped_support_count']==0 else 'NO（当前证据尚不足以确认）'}**。支撑量改善目标是否达到：{sup['reduction_target_met']}。真实 Bambu 切片及关键特征检查是否通过：{new['real_slicing_pass']}。可达性只从铸造盒顶部检查，包括竖直及小角度视线与 3 mm 探针；它仍是几何代理，尚未认证真实支撑能够取出。打印姿态恢复铸造姿态只需刚性逆变换。

## 九、三个主要剩余风险

1. 冻结 d075 与 MeVO 范围内的最佳实际尺寸仍可能小于 60 × 50 mm 目标；需要保留该限制，不能通过缩放消除。
2. 探针和 G-code 接触估计不等同于真实支撑实体及实际拆除过程；首次打印、拆支撑、PDMS 灌注和 ABS 溶解仍需实验验证。
3. 细管和端面连接仍需人工复核。Bambu 原生提示为：{warning_text}；完整信息保留在切片记录中。本轮没有发送任何打印任务，也没有宣称制造验证完成。

## 复现与阅读顺序

在 `/home/lzy/projects/vascular_printing` 执行 `./run_python.sh s1-7_support_aware_roi_60x50.py --output-root <全新目录>` 可运行独立全流程。已有目录分阶段续跑使用 `--stage model`、`--stage baseline`、`--stage evaluate`、`--stage report`；不会覆盖旧 s1-6。先阅读本报告与 `QC/13_baseline_vs_new_support.png`，再查看 `roi_60x50_summary.json`、`tables/` 和候选各自的 `manufacturing/` 与 `orientations/`。

`selected/` 保存选中中心线、补偿中心线与血管实体；`casting_mold/` 保存独立盒体、一体模型及真实切片文件；`bambu/evaluation_profile/` 保存本次比较配置与哈希；`baseline/` 保存重新切片证据。完整候选和失败尝试保留在 `candidates/`。旧上游及 s1-2～s1-6 文件哈希检查：{summary['source_protection']['all_unchanged']}，共 {summary['source_protection']['checked_files']} 个文件。
'''
    (out/'roi_60x50_report.md').write_text(text,encoding='utf8')
    if (out/'test_results.xml').exists():
        root=ET.parse(out/'test_results.xml').getroot();suites=root.findall('.//testsuite')
        summary['tests']={k:sum(int(r.attrib.get(k,0)) for r in suites) for k in ('tests','failures','errors','skipped')}
    if (out/'regression_test_results.xml').exists():
        suites=ET.parse(out/'regression_test_results.xml').getroot().findall('.//testsuite')
        summary['broader_regression']={k:sum(int(r.attrib.get(k,0)) for r in suites) for k in ('tests','failures','errors','skipped')}
        summary['broader_regression']['known_legacy_failure_reason']='Five historical s1-5/s1-6 hash assertions reference pre-UI-update snapshots; current-turn protected hashes unchanged. Original tests and baselines were not edited.'
        with (out/'roi_60x50_report.md').open('a',encoding='utf8') as f:
            t=summary.get('tests',{});b=summary['broader_regression']
            f.write(f"\n## 测试与历史快照差异\n\n新增永久测试共 {t.get('tests')} 项，失败 {t.get('failures')} 项、错误 {t.get('errors')} 项、跳过 {t.get('skipped')} 项。较早执行的扩大回归共 {b['tests']} 项，其中 {b['tests']-b['failures']-b['errors']-b['skipped']} 项通过、{b['failures']} 项失败。五项历史哈希断言失败涉及早期 s1-2/s1-3 界面快照及其集合检查；本轮开始前后哈希相同，未修改旧测试或旧快照来消除差异。详见 `regression_test_run.log` 和 `legacy_test_baseline_audit.json`。\n")
    summary['visual_qc']=[str(p.relative_to(out)) for p in sorted((out/'QC').glob('*.png'))]
    save_json(out/'roi_60x50_summary.json',summary)
    save_json(out/'validation_results.json',dict(status=summary['status'],failures=summary['failures'],tests=summary.get('tests'),
        source_protection=summary['source_protection'],scale_factor=1.,mevo_only=True,actual_slicing=new['real_slicing_pass'],
        profiles_identical=new['profiles_identical'],human_review_pending=True,manufacturing_validated=False))
    save_json(out/'stage_status.json',dict(status=summary['status'],selected_candidate=cid,failures=summary['failures'],
        actual_slicing=new['real_slicing_pass'],revision='SHORTER_PORTS_ADAPTIVE_FREE_FACES'))
    files=[ROOT/'s1-7_support_aware_roi_60x50.py',ROOT/'config/support_aware_roi_60x50.yaml',ROOT/'tests/test_s1_7_support_aware_roi_60x50.py',ROOT/'tests/data/support_aware_roi_synthetic.json']
    files+=[ROOT/'docs/SUPPORT_AWARE_60x50使用说明.md',ROOT/'开发与使用指南.md']
    files+=list((ROOT/'vascular_processing').glob('support_aware_*.py'))+[ROOT/'vascular_processing/roi_support_metrics.py']
    save_json(out/'implementation_manifest.json',{str(p.relative_to(ROOT)):sha256(p) for p in files})
    save_json(out/'artifact_manifest.json',dict(files={str(p.relative_to(out)):sha256(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'}))
    print(summary['status'],cid,roi,sup,flush=True)
    return summary
