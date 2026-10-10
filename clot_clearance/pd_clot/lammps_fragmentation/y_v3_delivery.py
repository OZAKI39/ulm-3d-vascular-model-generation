"""Independent saved-coordinate audit, Chinese review and lightweight evidence package."""
import json,csv,time,resource,shutil,zipfile,subprocess,collections
from pathlib import Path
import numpy as np
import pyvista as pv
import imageio.v2 as imageio
from PIL import Image
from .config import ROOT,dump
from .y_clot_morphology_v3 import OUT,RESULTS,load_geometry
from .y_bifurcation_sources import digest
from .native_y_wall_query import NativeYWallQuery
from .fragment_analysis import read_dump
from .failure_map import csv_out

def audit():
    started=time.perf_counter();g=load_geometry();n=g['clot_count'];q=NativeYWallQuery();cases=[];allrows=[];ports=[]
    for name in ['C0_ZERO_LOAD_DYNAMIC','C1_GAUSSIAN_A200']:
        p=OUT/'cases'/name;stats=json.loads((p/'WALL_CONTACT_DIAGNOSTICS.json').read_text());err=0.;worst=0.;outside=0;forceerr=0.;sampled=0
        for step,cols,a in read_dump(p/'trajectory.lammpstrj'):
            X=a[:n,[cols.index(k) for k in ['x','y','z']]];query=q.query(X);clear=query['clearance_m'];stored=a[:n,cols.index('d_wclear')];err=max(err,float(abs(clear-stored).max()));worst=min(worst,float(clear.min()));outside+=int((clear<0).sum());F=a[:n,[cols.index(f'f_nativewall[{i}]') for i in [1,2,3]]];normal=query['inward_normal'];mag=np.einsum('ij,ij->i',F,normal);forceerr=max(forceerr,float(np.linalg.norm(F-mag[:,None]*normal,axis=1).max()));assert mag.min()>-1e-12
            assert abs(a[n:,[cols.index(f'f_nativewall[{i}]') for i in [1,2,3]]]).max()==0
            assert np.linalg.norm(F[clear>20e-6+1e-12],axis=1).max(initial=0)<1e-12;sampled+=len(X)
        assert err<1e-12 and outside==0 and forceerr<1e-12
        cases.append(dict(case=name,wall_contact=stats,independent_saved_coordinate_checks=dict(points=sampled,outside_nodes=outside,clearance_agreement_max_m=err,inward_force_alignment_error_N=forceerr,anchors_wall_force_zero=True),performance=json.loads((p/'PERFORMANCE.json').read_text())))
        allrows.extend(csv.DictReader((p/'WALL_CONTACT_DIAGNOSTICS.csv').open()));ports.extend(stats['open_port_exits'])
    csv_out(OUT/'WALL_CONTACT_DIAGNOSTICS.csv',allrows);dump(OUT/'WALL_CONTACT_DIAGNOSTICS.json',dict(status='PASS',cases=cases,center_penetration_tolerance_m=2e-6,actual_maximum_penetration_m=0.,legal_exit_method='Actual native cap triangle crossing; no artificial cap force'))
    hand=json.loads((OUT/'FRAGMENT_HANDOFF.json').read_text());h={r['component_id']:r for r in hand['components']};tr=np.load(OUT/'cases/C1_GAUSSIAN_A200/accepted_trajectory.npz');X0=tr['position_m'][-1];events=json.loads((OUT/'FRAGMENT_TRANSPORT_EVENTS.json').read_text());outside_count=0;legal_count=0;total=0;max_out=0.;batches={}
    for r in csv.DictReader((OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv').open()):batches.setdefault(int(r['component_id']),[]).append([float(r[a+'_m']) for a in 'xyz'])
    for cid,positions in batches.items():
        rec=h[cid];offset=X0[np.array(rec['member_atom_ids'])-1]-rec['COM_xyz_m'];X=(np.asarray(positions)[:,None,:]+offset[None,:,:]).reshape(-1,3);z=q.query(X);total+=len(X)
        for i in np.flatnonzero(~z['inside_lumen']):
            legal=False
            for name,p in q.ports.items():
                distance=np.dot(X[i]-p['center'],p['outward'])
                if distance>=-2e-9 and q.port_crossing(X[i]-p['outward']*(max(distance,0)+1e-6),X[i]+p['outward']*2e-9):legal=True;break
            if legal:legal_count+=1
            else:outside_count+=1;max_out=max(max_out,-z['clearance_m'][i])
    assert outside_count==0,(outside_count,max_out)
    for e in events['events']:
        if e['event'].startswith('EXIT_'):ports.append(dict(stage='TRANSPORT',component_id=e['component_id'],time_s=e['time_s'],port=e['event'][5:],xyz_m=e['xyz_m'],mass_kg=e['mass_kg']))
    with (OUT/'OPEN_PORT_EVENTS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['stage','component_id','time_s','port','xyz_m','mass_kg']);w.writeheader();w.writerows(ports)
    dump(OUT/'INDEPENDENT_WALL_AND_TRANSPORT_AUDIT.json',dict(status='PASS',LAMMPS_cases=cases,transport_all_accepted_step_member_positions_checked=total,transport_illegal_outside_points=outside_count,transport_legal_port_exterior_points=legal_count,transport_maximum_illegal_penetration_m=max_out,wall_seconds=time.perf_counter()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
    historical=json.loads((OUT/'HISTORICAL_STAT_GUARD.json').read_text());changed=[]
    for a in historical:
        p=Path(a['path'])
        if not p.exists() or p.stat().st_size!=a['bytes'] or p.stat().st_mtime_ns!=a['mtime_ns']:changed.append(str(p))
    assert not changed,changed
    dump(OUT/'HISTORICAL_ASSETS_UNCHANGED.json',dict(status='PASS',files=len(historical),changed_paths=changed,method='Size and mtime_ns identity against start, read-only access; prior SHA registers reused. No historical bulk rehash.'))
    print('INDEPENDENT_AUDIT_PASS',total,legal_count,flush=True)

def deliver():
    g=load_geometry();n=g['clot_count'];hand=json.loads((OUT/'FRAGMENT_HANDOFF.json').read_text());trans=json.loads((OUT/'FRAGMENT_TRANSPORT_EVENTS.json').read_text());audit=json.loads((OUT/'INDEPENDENT_WALL_AND_TRANSPORT_AUDIT.json').read_text());visual=json.loads((OUT/'VISUALIZATION_PROVENANCE_AND_QA.json').read_text());native=json.loads((OUT/'NATIVE_EVENT_TRANSIENT_SIZE_AUDIT.json').read_text());geom=json.loads((OUT/'GEOMETRY_SUMMARY.json').read_text());cfg=json.loads((OUT/'CASE_CONFIG_FROZEN.json').read_text());outcomes=trans['components'];c1=json.loads((OUT/'cases/C1_GAUSSIAN_A200/result.json').read_text());records=[]
    for f in sorted((RESULTS/'figures').glob('*.png')):
        im=np.asarray(Image.open(f).convert('RGB'));assert im.std()>10;records.append(dict(path=str(f),bytes=f.stat().st_size,sha256=digest(f),shape=list(im.shape)))
    assert len(records)==5 and len(list((RESULTS/'figures').glob('*.pdf')))==5
    for p in (RESULTS/'figures').glob('*.pdf'):assert p.read_bytes().startswith(b'%PDF')
    movie=RESULTS/'animations/y_clot_realistic_fracture_and_transport.mp4';reader=imageio.get_reader(movie);count=reader.count_frames();metadata=reader.get_meta_data();reader.close();manifest=json.loads((RESULTS/'movie_frame_manifest/frames.json').read_text());assert count==len(manifest)==visual['movie_frames'];gif=Image.open(RESULTS/'animations/y_clot_realistic_fracture_and_transport.gif')
    for k in [0,2320,4660]:
        blocks=pv.read(RESULTS/'paraview'/f'fracture_step_{k:04d}.vtm');a=blocks['all_material_and_anchor_nodes'];tr=np.load(OUT/'cases/C1_GAUSSIAN_A200/accepted_trajectory.npz');i=int(np.flatnonzero(tr['step']==k)[0]);assert np.allclose(a.points,tr['position_m'][i],rtol=0,atol=1e-17) and a.n_points==len(g['X'])
    visual.update(figures_checked=records,MP4_frame_count=count,MP4_fps=metadata['fps'],GIF_frame_count=gif.n_frames,paraview_roundtrip_all_nodes=True,preview_visually_inspected=True)
    dump(OUT/'VISUALIZATION_PROVENANCE_AND_QA.json',visual)
    categories=collections.Counter(r['final_status'] for r in outcomes);branch={}
    for b in ['A','B']:
        selected=[r for r in outcomes if 'CLASSIFICATION_'+b in r['classification_sections']];exited=[r for r in outcomes if r['final_status']=='EXIT_OUTLET_'+b];branch[b]=dict(classified_count=len(selected),classified_mass_kg=sum(r['mass_kg'] for r in selected),exited_count=len(exited),exited_mass_kg=sum(r['mass_kg'] for r in exited),resolved_component_ids=[r['component_id'] for r in selected if r['node_count']>=3],unresolved_packet_count=sum(r['node_count']<=2 for r in selected))
    risk=dict(initial_mass_kg=hand['initial_clot_mass_kg'],attached_mass_kg=hand['attached_mass_kg'],detached_mass_kg=hand['detached_mass_kg'],resolved_detached_mass_kg=sum(r['mass_kg'] for r in hand['components'] if r['node_count']>=3),unresolved_debris_mass_kg=sum(r['mass_kg'] for r in hand['components'] if r['node_count']<=2),detached_ge200_ever=native['ge200_ever'],detached_ge300_ever=native['ge300_ever'],ge200_minimum_nodes=34,ge300_minimum_nodes=114,maximum_transient_detached_nodes=native['maximum_detached_nodes'],largest_final_detached_nodes=max(r['node_count'] for r in hand['components']),branches=branch,transport_final_status_counts=dict(categories),clinical_size_or_risk_prediction='NOT_ESTABLISHED',reference_scale_caveat='200um literature reference,300um engineering large component; neither is a universal patient danger threshold')
    assert abs(risk['attached_mass_kg']+risk['detached_mass_kg']-risk['initial_mass_kg'])<1e-20
    dump(OUT/'FRAGMENT_MASS_AND_RISK_PROXIES.json',risk)
    perf=dict(cases=[dict(case=z['case'],**z['performance'],callback_seconds_rank_sum=z['wall_contact']['callback_seconds_rank_sum'],wall_query_seconds_rank_sum=z['wall_contact']['query_seconds_rank_sum']) for z in audit['LAMMPS_cases']],formal_LAMMPS_dynamic_solves=2,run0_solves=1,contact_parameter_corrections=0,transport_wall_seconds=trans['wall_seconds'],transport_peak_RSS_bytes=trans['peak_RSS_bytes'],transport_wall_query_seconds=trans['wall_query_seconds'],native_event_postprocess_wall_seconds=native['wall_seconds'],independent_audit_wall_seconds=audit['wall_seconds'],visualization=visual,no_benchmark_sweep=True)
    perf['static_figure_label_review_wall_seconds']=json.loads((OUT/'STATIC_FIGURE_LABEL_REVIEW.json').read_text())['wall_seconds']
    for row in perf['cases']:
        z=json.loads((OUT/'cases'/row['case']/'result.json').read_text());row['initial_postprocess_wall_seconds']=z['postprocess_wall_s'];row['lineage_postprocess_wall_seconds']=z['lineage_summary']['size_postprocess_wall_s']
    dump(OUT/'PERFORMANCE_REAL_CASES.json',perf)
    status=dict(task='Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1',task_state='COMPLETE_WITH_DECLARED_EXPLORATORY_LIMITATIONS',project_stage='ACTIVE',gates={'A':'PASS','B':'PASS','C':'PASS','D':'PASS'},V3_MORPHOLOGY='PASS',NATIVE_Y_WALL_FIT='PASS',LAMMPS_V3_RUN0='PASS',ZERO_LOAD_DYNAMICS='PASS',NATIVE_Y_WALL_CONTACT='PASS',NATIVE_Y_WALL_PENETRATION='WITHIN_DECLARED_TOLERANCE',OPEN_PORT_CLASSIFICATION='PASS',FRACTURE_C1='EXECUTED',DETACHED_COMPONENTS='OBSERVED',FRAGMENT_HANDOFF='PASS',ONE_WAY_BACKGROUND_FLOW_TRANSPORT='PASS',BRANCH_CLASSIFICATION={'A':'OBSERVED_A','B':'OBSERVED_B'},REALISTIC_MEDICAL_RENDER='PASS',SOLVER_VIEW_AND_MOVIE_TRACEABILITY='PASS',CLOT_PRESENT_CFD_SOLVED='NOT_RUN',BUBBLE_STREAMING_CONNECTED='NOT_RUN',BIOLOGICAL_QUANTITATIVE_VALIDATION='NOT_ESTABLISHED',PATIENT_SPECIFIC='NOT_ESTABLISHED',FRAGMENT_SIZE_CLINICAL_PREDICTION='NOT_ESTABLISHED',PRIMARY_FRAGMENTATION_SOLVER='LAMMPS_BPM_PERI_LPS',CONTACT_MODEL='MODIFIED_1P0DX / CONDITIONAL_EXPLORATORY_DEFAULT',WALL_CONTACT_MODEL='NATIVE_Y_WALL_REPULSION_V1',FRACTURE_SCOPE='EXPLORATORY_ONLY',B1_B2='TERMINATED_BY_PROJECT_POLICY_CHANGE',PD2_FRACTURE='STOPPED',PHASE_FIELD_FRACTURE='ARCHIVED_OPTIONAL_REFERENCE_ONLY',scope_completed=True,own_jobs_remaining=[],automatic_followup=False,commit_or_push=False)
    dump(OUT/'FINAL_STATUS.json',status)
    review=f'''# Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1 审查

本轮授权范围已完成，Gate A/B/C/D 均 PASS；这里的 PASS 是几何、程序和证据链通过，不是生物物理定量验证。仅执行 C0、C1 两次正式 LAMMPS 动力学和一次双 MPI 回调 run0。没有增加载荷、材料校准、形态扫描或 CFD 求解。

1. **V3 形态和真实壁面。** 在同一只读 Y WALL、x≈−6mm、+z 内壁，用确定性非对称足迹和三个306–330µm宽的高斯隆起构造连续贴壁实体。轮廓包括上游厚肩、渐薄尾部和非镜像侧缘；配置在动力学前冻结，没有用结果反调形。连续体积 {geom['continuous_mesh_volume_m3']*1e9:.9f}mm³，比V2小4.936%；1927个新50µm立方节点体积 {geom['clot_volume_m3']*1e9:.6f}mm³，比连续体小11.304%。1561个固定锚点分三层，重新生成typed bonds。27点体素支持检查、真实管腔membership和非法键排查通过。相机和实际比例相同的V2/V3图可肉眼区分三个宽隆起；力学仍有50µm离散分辨率。

2. **求解中的完整壁面接触。** 从原固定 MODIFIED_1P0DX 源码构建独立共享库，使用两rank OpenMPI的fix external逐步回调；先完成LPS力密度到N的转换，再加原生三角壁力。最近三角形法向由邻接原生四面体向内定向，内外判据来自原四面体，不把开放WALL冒充闭合体。查询覆盖完整母管、分叉和子管。中心碰撞代理半径20µm、刚度50N/m、仅接近时阻尼ζ=.5；dtω=.316。C0全程零力零速度无断键。C1有79节点实际受到壁力；最小真实间隙18.012069µm（0.36024dx），最大真实中心穿壁深度0、穿壁节点质量0、持续时间0、扫掠穿越0。没有穿壁发生位置。允许中心穿壁容差2µm；代理半径重叠约1.988µm不等于穿出血管。C1没有真实端口离域。每步rank记录和独立重读实际坐标核查保留，没有LAMMPS坐标夹回，没有接触参数修正。

3. **真实断裂。** C1人工Gaussian参考平均预算200Pa，冻结参考表面面积与方向，包络在当前空间位置计算；方向仍是初始露出表面投影，不是动态wet surface牵引。脉宽1.1605ms、dt=.5µs、观察2.33ms。参考局部峰值770.874Pa，保存帧实际局部峰值750.258Pa、实际平均预算峰值183.486Pa，200Pa并非均匀局部牵引。共12810次真实断键。bulk断键率{c1['bulk_broken_fraction']:.6%}，interface断键率{c1['interface_broken_fraction']:.6%}；附着{c1['attached_clot_volume_fraction']:.6%}，游离{1-c1['attached_clot_volume_fraction']:.6%}。真实材料身份保持Liu小应变LPS映射和porcine-skin peeling界面代理，并未变成患者血管粘附模型。2030个原生断键时刻都重建图，最大短暂游离组件13节点，末态最大6节点；≥200µm（34节点）及≥300µm（114节点）全程未出现。附着主体不计游离大块。更多断键/碎屑不等于更好治疗。

4. **单次handoff。** 预先指定 step4660、2.33ms末态，原生surviving bonds与缓存/lineage逐一核对。208个游离分量，共225节点，质量{hand['detached_mass_kg']:.12g}kg；其中一个6节点组件质量7.5e−10kg，其余207个1–2节点质量包共219节点、2.7375e−8kg。附着主体1702节点、{hand['attached_mass_kg']:.12g}kg保留在末态位置。FRAGMENT_HANDOFF保留全部atom IDs、父子成员身份、首次游离时刻、实际COM与速度。没有使用V2备用场景，没有把附着主体推入流线。

5. **真实背景CFD查询。** 从实际result_121.vtu读取Velocity/Pressure，原网格26712点、146230四面体，最大速度10.074017mm/s；每一步在当前COM作严格P1插值，不外推，原生四面体质心插值误差约1.7e−18m/s。原文件代表1.21s时保存的无血栓稳态水流。>=3节点采用等效水动力直径的Stokes惯性时间和Schiller–Naumann有限Re修正；唯一6节点组件初期Re代理最大67.09，因此没有假称全程线性Stokes适用。1–2节点质量包用声明的100µs数值松弛，不将网格直径当真实生物粒径。45个分量曾处于100µm近壁带，1个接近保留血栓节点。近壁梯度、保留血栓附近路径未经含血栓CFD认证，分别标记NEAR_WALL/RETAINED_CLOT_PROXIMITY_UNCERTAIN。没有流体反馈、碎片相撞、保留主体碰撞、升力、附加质量、旋转和变形。

6. **实际分流和出口。** A分类截面65分量，质量{branch['A']['classified_mass_kg']:.12g}kg；B分类截面76分量，质量{branch['B']['classified_mass_kg']:.12g}kg。A真实出口62分量、{branch['A']['exited_mass_kg']:.12g}kg；B真实出口71分量、{branch['B']['exited_mass_kg']:.12g}kg。分类截面在子管10mm，真实出口12mm；8个已分类但未出口，67个5s内未到分类截面（其中15个COM仍为x≤0，52个已为x>0，因此不能统称未进入分叉），受初始位置、近壁速度和有限观察时间影响，不能据此认定真实粘附或闭塞。唯一6节点组件1259在4.241322s到A出口，其余出口事件是未解析质量包。没有放大血速或投到中心线。统计是此简化模型的工程演示，不是临床分流/捕获率。

7. **输运阶段壁面。** 与LAMMPS共用native_y_wall_query。指数松弛每步dt≤5ms、位移≤20µm；在积分子步内，对组件所有真实成员偏置作单向无摩擦约束，保持刚体平移，去除向壁法向速度，中心安全代理5µm。1403次约束子步有显式日志；接受态最小间隙约5.00002µm，无非法中心穿壁，扫掠防护拒绝数{sum(r['swept_guard_rejections'] for r in outcomes)}。合法出口按实际三角端口判断；部分成员先越出口的状态会保留，质心越出口才结束整个组件跟踪。独立重读所有接受步的{audit['transport_all_accepted_step_member_positions_checked']}个成员坐标，确认非法越壁为0。出口后的画面保留终止位置标记和统计，不继续外推速度。

8. **医学展示与原始证据。** 每个真实组件独立从当前节点的立方体支持并集重建闭合外壳，25µm采样、有限角部平滑，顶点圆滑位移限制12.5µm并限制在本组件支持内；不跨组件合并，不凭空增添大块。初始圆滑壳体积约0.227076mm³，比节点积分体积小约5.73%，仅作显示、从不替代质量或粒径计算。亚25µm裂缝和微结构不具可视分辨率；实际断键仍在raw节点/全部surviving bonds/事件/锚点/载荷箭头和VTM中审查。小碎屑显示细点，质量全保留。输运使用handoff外壳仅平移，主体不动。视频244帧，先ms断裂后明确切换seconds one-way；全Y、局部、分叉三视角，实际尺度，无位移放大。五张PNG/PDF及同帧raw/medical对照、可重载ParaView VTM与逐帧来源清单齐全。

9. **实际成本与软件检查。** C0耗时{perf['cases'][0]['wall_seconds']:.3f}s、C1耗时{perf['cases'][1]['wall_seconds']:.3f}s；两rank壁面查询累计分别{perf['cases'][0]['wall_query_seconds_rank_sum']:.3f}s和{perf['cases'][1]['wall_query_seconds_rank_sum']:.3f}s（rank和，不能当成额外墙钟相加）。成功输运{trans['wall_seconds']:.3f}s，其中共享壁面查询{trans['wall_query_seconds']:.3f}s；原生事件尺寸审查{native['wall_seconds']:.3f}s；独立坐标审查{audit['wall_seconds']:.3f}s；图和VTM{visual['figures_paraview_wall_seconds']:.3f}s，视频{visual['movie_wall_seconds']:.3f}s。各过程峰值RSS见PERFORMANCE_REAL_CASES。安装仅新共享库及隔离mpi4py；初期插值检查中float32质心舍入问题改为先转float64，输运部分成员越端口的程序问题也已修复，两个失败启动日志保留。没有增加正式LAMMPS工况，没有重新标定/流场重算/压力测试/网格收敛/GPU或材料扫描。

10. **仍缺失的物理。** 真正血栓面traction需含血栓流域、可靠近壁分辨率与牵引映射；真实微泡streaming需声场/气泡动力学及与血栓的耦合；动态wet surfaces需随裂面更新的几何、流体边界和载荷；定量输运需验证碎片形态/尺度/密度、近壁水动力、接触/粘附、旋转变形、相互作用及阻塞引起的流场更新。还需要实验与不确定性验证。本轮不能称作FSI、患者特异模型、生物粒径预测或临床风险模型。

代码入口：y_clot_morphology_v3.py、native_y_wall_query.py、native_y_wall_contact.py、y_v3_cases.py、fragment_handoff.py、background_flow_sampler.py、fragment_transport.py、medical_fragment_renderer.py、y_v3_delivery.py。配置及源文件hash在证据包中，历史文件仅只读；正式结果留在原目录，大文件只在轻量ZIP中列路径/尺寸/SHA。

技术来源：[LAMMPS fix external](https://docs.lammps.org/fix_external.html)、[Python回调接口](https://docs.lammps.org/Python_module.html)、[MPI安装约定](https://docs.lammps.org/Python_install.html)、[OpenFOAM中Schiller–Naumann公式实现](https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html)。接口同时与固定本地源码核对。无新增LAMMPS solver或材料定律。

本轮到此停止，等待人工审查。B1/B2、PD2保持停止；CLOT_PRESENT_CFD_SOLVED=NOT_RUN，BUBBLE_STREAMING_CONNECTED=NOT_RUN。
'''
    name='REVIEW_Y_CLOT_REALISTIC_MORPHOLOGY_WALL_CONTACT_AND_FLOW_TRANSPORT_V1.md';(OUT/name).write_text(review);shutil.copy2(OUT/name,RESULTS/name);shutil.copy2(OUT/'FINAL_STATUS.json',RESULTS/'FINAL_STATUS.json')
    # Current source snapshots, minimal rerun recipe and all case identities remain reviewable.
    code=OUT/'source';code.mkdir(exist_ok=True)
    for f in (ROOT/'pd_clot/lammps_fragmentation').glob('*.py'):shutil.copy2(f,code/f.name)
    (OUT/'REPRODUCE.txt').write_text('Use existing scientific python with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1. New output directories only: frozen results refuse overwrite.\npython -m pd_clot.lammps_fragmentation.y_clot_morphology_v3 freeze / geometry\npython -m pd_clot.lammps_fragmentation.y_v3_cases freeze / run0 / run 0 / run 1\npython -m pd_clot.lammps_fragmentation.fragment_transport freeze (before dynamics)\npython -m pd_clot.lammps_fragmentation.fragment_handoff\npython -m pd_clot.lammps_fragmentation.fragment_transport run\npython -m pd_clot.lammps_fragmentation.y_v3_delivery audit\npython -m pd_clot.lammps_fragmentation.medical_fragment_renderer\npython -m pd_clot.lammps_fragmentation.medical_fragment_renderer captions\npython -m pd_clot.lammps_fragmentation.y_v3_delivery deliver\nExact MPI and library environment: cases/*/RUN_COMMAND.json. Shared library CMake args/logs and isolated mpi4py metadata preserved. Do not rerun frozen geometry or scientific cases in place.\n')
    selected=[];large=[]
    for base in [OUT,RESULTS]:
        for p in sorted(base.rglob('*')):
            if not p.is_file() or p.name in ['evidence.zip','EVIDENCE_PACKAGE.json','LARGE_ARRAY_AND_MOVIE_REGISTER.json'] or 'source_before' in p.parts:continue
            # Keep big native traces in place, with hashes. ZIP includes small event/lineage/key data.
            if p.stat().st_size>8_000_000 or p.suffix in ['.mp4','.gif','.restart']:large.append(dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)))
            else:selected.append(p)
    dump(OUT/'LARGE_ARRAY_AND_MOVIE_REGISTER.json',large);selected.append(OUT/'LARGE_ARRAY_AND_MOVIE_REGISTER.json')
    with zipfile.ZipFile(OUT/'evidence.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in selected:z.write(p,p.relative_to(ROOT))
    dump(OUT/'EVIDENCE_PACKAGE.json',dict(path=str(OUT/'evidence.zip'),bytes=(OUT/'evidence.zip').stat().st_size,sha256=digest(OUT/'evidence.zip'),included_files=len(selected),large_files_preserved_in_place=len(large)))
    print('DELIVERY_COMPLETE',status['gates'],dict(categories),flush=True)

if __name__=='__main__':
    import sys
    {'audit':audit,'deliver':deliver}[sys.argv[1]]()
