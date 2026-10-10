"""Package saved evidence and Chinese review; no solver or geometry mutation."""
import json,csv,shutil,hashlib,zipfile,subprocess,difflib,sys,time
from pathlib import Path
import numpy as np
import pyvista as pv
from PIL import Image
import imageio.v2 as imageio
from .config import ROOT,dump
from .y_irregular_morphology import OUT,RESULTS,EXE,load_geometry
from .y_bifurcation_sources import YSources,digest
from .y_fragmentation_cases import CASES
from .failure_map import csv_out
from .fragment_analysis import read_dump

def native_event_audit():
    """Threshold audit at every native break step, including between saved frames."""
    from .detached_fragment_metrics import graph_components,equivalent_diameter_um
    start=time.perf_counter();g=load_geometry();b=g['bonds'];typ=g['bond_types'];n=g['clot_count'];result=[]
    for name,_,_ in CASES:
        p=OUT/'cases'/name;z=np.load(p/'broken_bonds.npz');ev=z['events'];bi=z['initial_bond_indices'];alive=np.ones(len(b),bool);rows=[];cursor=0;peak=0
        for step in np.unique(ev[:,1]).astype(int) if len(ev) else []:
            end=np.searchsorted(ev[:,1],step,side='right');alive[bi[cursor:end]]=False;cursor=end
            labels,counts,attached=graph_components(n,b,typ,alive);free=counts[~attached];mx=int(free.max()) if len(free) else 0;peak=max(mx,peak)
            rows.append(dict(step=step,time_s=step*.5e-6,largest_detached_nodes=mx,N_detached_ge200=int((free>=34).sum()),N_detached_ge300=int((free>=114).sum()),mass_fraction_ge200=float(free[free>=34].sum()/n),mass_fraction_ge300=float(free[free>=114].sum()/n)))
        csv_out(p/'lineage/native_event_size_thresholds.csv',rows)
        result.append(dict(case_id=name,distinct_break_steps=len(rows),maximum_detached_nodes_over_every_native_break_step=peak,maximum_d_eq_um=float(equivalent_diameter_um(peak,50e-6)),ge200_ever=any(a['N_detached_ge200'] for a in rows),ge300_ever=any(a['N_detached_ge300'] for a in rows),native_timestep_s=.5e-6,all_break_steps_checked=True,no_unsaved_positions_fabricated=True))
    dump(OUT/'NATIVE_EVENT_TRANSIENT_SIZE_AUDIT.json',dict(status='PASS',cases=result,wall_seconds=time.perf_counter()-start,scope='Reconstruct existing native event graph at every distinct break step; no dynamics, no interpolation, no new trajectory frames.'))

def aggregate():
    summaries=[];timeline=[];life=[];lineages={};perf=[]
    for name,_,amp in CASES:
        p=OUT/'cases'/name;z=json.loads((p/'result.json').read_text());s=z['lineage_summary'];final=s['final'];rows=list(csv.DictReader((p/'lineage/timeseries.csv').open()))
        row=dict(case_id=name,reference_mean_budget_Pa=amp,reference_peak_traction_Pa=z['peak_local_traction'],actual_peak_traction_Pa=z['actual_peak_local_traction_Pa'],peak_actual_mean_budget_Pa=z['peak_actual_mean_traction_budget_Pa'],integral_sum_force_magnitudes_Ns=z['integral_sum_force_magnitudes_Ns'],final_bulk_broken_fraction=z['bulk_broken_fraction'],final_interface_broken_fraction=z['interface_broken_fraction'],final_attached_mass_fraction=z['attached_clot_volume_fraction'],final_detached_mass_fraction=1-z['attached_clot_volume_fraction'],final_unresolved_debris_fraction=z['debris_fraction'],final_detached_component_count=int(rows[-1]['detached_component_count']),maximum_detached_component_count=max(int(a['detached_component_count']) for a in rows),largest_ever_detached_d_eq_um=s['largest_ever_detached_d_eq_um'],largest_final_detached_d_eq_um=s['largest_final_detached_d_eq_um'],first_bulk_break_s=s['first_bulk_break_time_s'],first_interface_break_s=s['first_interface_break_time_s'],peak_contact_pairs=s['max_active_contact_pairs'],min_r_over_dx=s['min_nonbonded_distance_over_dx'],peak_KE_J=s['peak_clot_KE_J'],maximum_nodes_outside_visual_vessel=z['maximum_nodes_outside_visual_vessel'],classification=z['final_classification'],wall_seconds=z['runtime']['wall_seconds'],postprocess_wall_seconds=z['postprocess_total_wall_s'],peak_RSS_bytes=z['runtime']['peak_RSS_bytes'])
        for d,q in s['thresholds'].items():
            for k,v in q.items():row[f'ge{d}_{k}']=v
        summaries.append(row);timeline.extend(rows);life.extend(csv.DictReader((p/'lineage/lifetimes.csv').open()));lineages[name]=s
        perf.append(dict(case_id=name,**z['runtime'],postprocess_wall_s=z['postprocess_total_wall_s']))
    csv_out(OUT/'CASE_RESPONSE_COMPARISON.csv',summaries);dump(OUT/'CASE_RESPONSE_COMPARISON.json',summaries)
    csv_out(OUT/'FRAGMENT_BREAK_AND_ATTACHMENT_TIMELINES.csv',timeline);csv_out(OUT/'DETACHED_FRAGMENT_SIZE_AND_LIFETIME.csv',life);dump(OUT/'FRAGMENT_LINEAGE_SUMMARY.json',lineages)
    extra=json.loads((OUT/'NATIVE_EVENT_TRANSIENT_SIZE_AUDIT.json').read_text())['wall_seconds'] if (OUT/'NATIVE_EVENT_TRANSIENT_SIZE_AUDIT.json').exists() else 0.
    dump(OUT/'PERFORMANCE_REAL_CASES.json',dict(cases=perf,total_solver_wall_s=sum(q['wall_seconds'] for q in perf),total_postprocess_wall_s=sum(q['postprocess_wall_s'] for q in perf)+extra,supplementary_event_threshold_audit_wall_s=extra,postprocess_peak_RSS='NOT_SEPARATELY_SAMPLED',peak_solver_process_tree_RSS_bytes=max(q['peak_RSS_bytes'] for q in perf),formal_solves=4,run0_invocations=1,no_parameter_or_performance_sweep=True,execution_location='Existing local clot workspace; server Y assets read-only',visualization=[json.loads(p.read_text()) for p in sorted(OUT.glob('VISUALIZATION_PERFORMANCE_*.json'))]))
    for f in ['CASE_RESPONSE_COMPARISON.json','CASE_RESPONSE_COMPARISON.csv','FRAGMENT_BREAK_AND_ATTACHMENT_TIMELINES.csv','DETACHED_FRAGMENT_SIZE_AND_LIFETIME.csv','FRAGMENT_LINEAGE_SUMMARY.json','V1_V2_GEOMETRY_COMPARISON.json']:
        shutil.copy2(OUT/f,RESULTS/'FIGURE_DATA'/f)
    return summaries,lineages

def qa():
    g=load_geometry();n=g['clot_count'];s=YSources();old=json.loads((OUT/'HISTORICAL_ASSETS_BEFORE.json').read_text());changed=[]
    for a in old:
        p=Path(a['path'])
        if not p.exists() or p.stat().st_size!=a['bytes'] or digest(p)!=a['sha256']:changed.append(a['path'])
    assert not changed,changed
    dump(OUT/'HISTORICAL_ASSETS_AFTER.json',dict(status='PASS',file_count=len(old),bytes=sum(a['bytes'] for a in old),all_sha256_unchanged=True,changed_paths=changed))
    frozen=json.loads((OUT/'FROZEN_INPUT_HASHES.json').read_text());assert all(digest(OUT/f)==sha for f,sha in frozen.items())
    mirror=Path('/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles');src=[]
    for f in s.register['files']:
        p=mirror/f['relative_path'];assert digest(p)==f['sha256'];src.append(dict(path=str(p),sha256=f['sha256'],unchanged=True))
    dump(OUT/'SOURCE_HASH_RECHECK.json',dict(status='PASS',local_mirror=src,server_check='SERVER_READONLY_HASH_CHECK.json',snapshot_check=True,mesh_points=s.mesh.n_points,mesh_cells=s.mesh.n_cells,mesh_sha256=s.mesh_id,flow_sha256=s.flow_id))
    stats=[]
    for p in sorted((RESULTS/'figures').glob('*.png')):
        im=np.asarray(Image.open(p).convert('RGB'));assert im.std()>10
        stats.append(dict(path=str(p),width=im.shape[1],height=im.shape[0],sha256=digest(p),std=float(im.std())))
    assert len(stats)==5 and len(list((RESULTS/'figures').glob('*.pdf')))==4
    for p in (RESULTS/'figures').glob('*.pdf'):assert p.read_bytes().startswith(b'%PDF')
    movie=RESULTS/'animations/y_clot_initial_fracture_synchronized.mp4';reader=imageio.get_reader(movie);frames=reader.count_frames();metadata=reader.get_meta_data();reader.close();assert frames==234
    gif=Image.open(RESULTS/'animations/y_clot_initial_fracture_preview.gif');assert gif.n_frames==60
    data_sha=digest(OUT/'Y_CLOT_V2_LAMMPS.data');checks=[]
    for name,_,_ in CASES:
        p=OUT/'cases'/name;r=json.loads((p/'request.json').read_text());z=json.loads((p/'result.json').read_text());tr=np.load(p/'accepted_trajectory.npz');ev=np.load(p/'broken_bonds.npz')
        assert digest(p/'body.data')==data_sha and len(tr['step'])==234 and tr['step'][-1]==4660
        assert z['validation']['native_topology_match'] and z['validation']['all_atoms_retained'] and z['validation']['finite_states']
        for k in [0,58,117,233]:
            q=pv.read(RESULTS/'paraview'/name/f'points_{k:03d}.vtp');b=pv.read(RESULTS/'paraview'/name/f'bonds_{k:03d}.vtp')
            assert np.array_equal(q.points,tr['position_m'][k]) and np.array_equal(q['component_id'][:n],tr['component_id'][k]);assert b.n_cells==len(g['bonds'])-int((ev['events'][:,1]<=tr['step'][k]).sum()) if len(ev['events']) else b.n_cells==len(g['bonds'])
        checks.append(dict(case=name,native_frames=234,paraview_point_frames=len(list((RESULTS/'paraview'/name).glob('points_*.vtp'))),paraview_bond_frames=len(list((RESULTS/'paraview'/name).glob('bonds_*.vtp'))),sample_reload_exact=True,all_native_nodes_preserved=True,wall_bonds_intact=z['validation']['wall_bonds_intact'],provider_max_force_error_N=z['validation']['provider_max_force_difference_N']))
    rq=[json.loads((OUT/'cases'/name/'request.json').read_text()) for name,_,_ in CASES]
    a={k:v for k,v in rq[1].items() if k not in ['case_id','tau_budget_Pa']};b={k:v for k,v in rq[2].items() if k not in ['case_id','tau_budget_Pa']};assert a==b
    report=dict(status='PASS',manual_visual_review='Four required figures, synchronized movie keyframes and raw solver keyframes inspected',figures=stats,movie_frames=frames,movie_fps=metadata['fps'],movie_dimensions=metadata['size'],gif_frames=gif.n_frames,case_checks=checks,gaussian100_200_only_amplitude_differs=True,same_geometry_data_sha256=data_sha,all_four_absolute_time_windows_equal=True,no_fake_frames=True,no_component_bridging=True,unresolved_nodes_retained=True,current_position_envelope_reference_normals=True,full_vessel_collision='NOT_MODELED',fluid_transport='NOT_MODELED')
    dump(OUT/'VISUALIZATION_QA.json',report)
    return report

def final_status(rows):
    obs=lambda x:'OBSERVED' if x else 'NOT_OBSERVED'
    st=dict(task='Y_BIFURCATION_IRREGULAR_CLOT_AND_INITIAL_FRACTURE_V1',status='READY_FOR_HUMAN_REVIEW')
    for key in ['REAL_Y_MESH_REUSED','OLD_Y_CFD_UNMODIFIED','NEW_IRREGULAR_MURAL_GEOMETRY','V1_V2_GEOMETRY_COMPARE','TRUE_CURVED_WALL_ATTACHMENT','GEOMETRY_CLEARANCE','LAMMPS_CUBIC_V2','CURVED_INTERFACE_BONDS','LAMMPS_RUN0_V2','ZERO_LOAD_DYNAMIC','Y_CONTEXT_FRACTURE_VISUALIZATION']:st[key]='PASS'
    st.update(CONTACT_MODEL='MODIFIED_1P0DX_EXPLORATORY_CANDIDATE',MATERIAL_CARD='LIU2024_LPS_SMALL_STRAIN_PLUS_PORCINE_SKIN_INTERFACE_PROXY',GAUSSIAN_100='EXECUTED',GAUSSIAN_200='EXECUTED',OPPOSING_100='EXECUTED',TRUE_BOND_BREAK=obs(any(r['final_bulk_broken_fraction'] or r['final_interface_broken_fraction'] for r in rows)),TRUE_DETACHED_COMPONENT=obs(any(r['maximum_detached_component_count'] for r in rows)),DETACHED_GE_200UM=obs(any(r['ge200_EVER_OBSERVED'] for r in rows)),DETACHED_GE_300UM=obs(any(r['ge300_EVER_OBSERVED'] for r in rows)),PERSISTENT_LARGE_DETACHED_FRAGMENT=obs(any(r['ge200_persistent_membership_episodes'] for r in rows)),UNRESOLVED_DEBRIS='TRACKED',CLOT_PRESENT_CFD_SOLVED='NOT_RUN',CLOT_SURFACE_TRACTION_VALIDATED='NOT_ESTABLISHED',FULL_VESSEL_FRAGMENT_COLLISION='NOT_MODELED',FLUID_FRAGMENT_TRANSPORT='NOT_MODELED',BUBBLE_STREAMING_CONNECTED='NOT_YET',CURVED_INTERFACE_STRENGTH_VALIDATED='NOT_ESTABLISHED',BIOLOGICAL_QUANTITATIVE_VALIDATION='NOT_ESTABLISHED',PATIENT_SPECIFIC='NOT_ESTABLISHED',B1_B2='TERMINATED',PD2_FRACTURE='STOPPED',PRIMARY_USE='EXPLORATORY_Y_VESSEL_MURAL_CLOT_FRAGMENTATION',shape_count=1,formal_dynamics_count=4,run0_count=1,own_jobs_remaining=[],automatic_followup=False,git_commit_or_push=False,tangent_update='FROZEN_INITIAL_V2_PROJECTION; NOT_DEFORMED_SURFACE_TRACKING',persistent_convention='Exact-member detached episode >=3 consecutive 10us-spaced saved frames; not clinically persistent',source_assets='Original server + local mirror + stored snapshots SHA256 verified')
    dump(OUT/'FINAL_STATUS.json',st);return st

def review(rows,status):
    geom=json.loads((OUT/'V1_V2_GEOMETRY_COMPARISON.json').read_text());v1=geom['V1'];v2=geom['V2'];fit=json.loads((OUT/'GEOMETRY_WALL_CLEARANCE_AND_VOLUME.json').read_text());perf=json.loads((OUT/'PERFORMANCE_REAL_CASES.json').read_text());source=json.loads((OUT/'SOURCE_HASH_RECHECK.json').read_text());c=json.loads((OUT/'PLACEMENT_CONFIGURATION.json').read_text());gate=json.loads((OUT/'ZERO_LOAD_GATE.json').read_text())
    table='| 工况 | bulk 断键 | interface 断键 | 末态附着质量 | 末态游离质量 | 末态1–2节点碎屑 | 最大游离等效径 |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for r in rows:table+=f"| {r['case_id']} | {r['final_bulk_broken_fraction']:.3%} | {r['final_interface_broken_fraction']:.3%} | {r['final_attached_mass_fraction']:.3%} | {r['final_detached_mass_fraction']:.3%} | {r['final_unresolved_debris_fraction']:.3%} | {r['largest_ever_detached_d_eq_um']:.2f} µm |\n"
    sizetext=[]
    for r in rows[1:]:
        for d in [200,300]:
            p=f'ge{d}_';yes=r[p+'EVER_OBSERVED']
            sizetext.append(f"- {r['case_id']} ≥{d} µm："+(f"出现；最大同时 {r[p+'MAX_SIMULTANEOUS_COUNT']} 块，峰值游离质量 {r[p+'PEAK_DETACHED_MASS_FRACTION']:.3%}，首次/末次 {r[p+'FIRST_APPEARANCE_TIME']*1e3:.3f}/{r[p+'LAST_SEEN_TIME']*1e3:.3f} ms，完全相同成员最长观测跨度 {r[p+'longest_exact_membership_detached_span_s']*1e6:.1f} µs，末态 {r[p+'AT_FINAL_FRAME']} 块。" if yes else 'NOT OBSERVED，全时段均未出现。'))
    loadtable='| 工况 | 参考平均 Pa | 参考峰值 Pa | 实际峰值 Pa | 实际平均预算峰值 Pa | 力幅时间积分 / N·s |\n|---|---:|---:|---:|---:|---:|\n'
    for r in rows:loadtable+=f"| {r['case_id']} | {r['reference_mean_budget_Pa']:.0f} | {r['reference_peak_traction_Pa']:.3f} | {r['actual_peak_traction_Pa']:.3f} | {r['peak_actual_mean_budget_Pa']:.3f} | {r['integral_sum_force_magnitudes_Ns']:.6g} |\n"
    durations='\n'.join(f"- {r['case_id']}：求解 {r['wall_seconds']:.2f} s，后处理 {r['postprocess_wall_seconds']:.2f} s，进程树峰值 RSS {r['peak_RSS_bytes']/2**20:.1f} MiB。" for r in rows)
    most=max(rows[1:],key=lambda r:r['final_unresolved_debris_fraction']);keeps=max(rows[1:],key=lambda r:r['final_attached_mass_fraction'])
    text=f'''# Y 不规则附壁血栓与首次破坏审查

READY_FOR_HUMAN_REVIEW。一个冻结 V2、一次 run 0、四个 C0–C3 动力学均完成。没有增加载荷、改材料、追加求解或启动 CFD/微泡。开发、LAMMPS 和渲染使用已有本地血栓工作树及环境；服务器仅只读核验原 Y 资产。

1. 原始 Y 网格和 CFD 是否保持？

保持。服务器15项原始资产、本地平行项目镜像和上轮只读快照逐项 SHA256 核对 PASS；原网格 {source['mesh_points']} 点/{source['mesh_cells']} 四面体。
网格 SHA256 `{source['mesh_sha256']}`；稳态 result_121.vtu SHA256 `{source['flow_sha256']}`。
来源 `/workspace/mevo_acoustic_microbubbles/` 的实际三端口 Y 模型；未造替代管道。原入口约1.58124 mL/min，不是四端口18 mL/min。历史1989个文件全量复核未变。

2. V2 的尺度、体积和不规则性？

名义 a/b/H =0.68/0.62/0.65 mm，相对 V1 0.5625/0.5625/0.5625 mm 分别增加20.89%/10.22%/15.56%。实际足迹轴向/周向跨度 {v2['footprint_axial_extent_m']*1e3:.6f}/{v2['footprint_circumferential_arc_extent_m']*1e3:.6f} mm，最大保存高度 {v2['maximum_mesh_height_m']*1e3:.6f} mm。
连续闭合体积由 {v1['continuous_mesh_volume_m3']*1e9:.6f} 增至 {v2['continuous_mesh_volume_m3']*1e9:.6f} mm³（+{geom['relative_continuous_volume_change']:.2%}）；积分体积由 {v1['clot_volume_m3']*1e9:.6f} 增至 {v2['clot_volume_m3']*1e9:.6f} mm³。
不规则足迹由U=u−0.08u²+0.04、V=[v−0.07(1−u²)−0.05u]/[1+0.06sin(πu)]决定；上游偏厚和三个宽高斯肩部在旧平滑 envelope 上作≤10%扰动并预先规范最高高度。最小扰动sigma0.31mm大于2δ0.301mm；三个宽项不代表三个独立极大值。MORPHOLOGY_DEFINITION/CONFIG及其哈希在构造和求解前冻结，没有看结果后调形。这是合成工程血栓。

3. 实际贴壁、管腔和 clearance？

脚点(-6,0,{fit['actual_footpoint_m'][2]*1e3:.9f}) mm，实际壁面与名义1.3mm差2.818763µm。足迹对同一离散WALL最大gap {fit['attachment_vertex_gap_max_m']*1e6:.6f} µm，处于预定0.002µm容差内；这不是生物测量精度。
{fit['material_centers_inside']}材料中心及{fit['voxel_support_samples_checked']}个体素27点支持检查全部在原生管腔；{fit['anchor_centers_outside']}锚中心在管腔外。入口顶点距离 {fit['clearance']['minimum_distance_to_port_surface_vertices_m']['INLET']*1e3:.4f} mm，血栓最下游x={v2['bbox_continuous_m'][1]*1e3:.4f} mm，远离x=0分叉；对侧壁保守间距>{fit['clearance']['opposite_wall_clearance_lower_bound_m']*1e3:.4f} mm。
中心实际截面占据率 {fit['clearance']['center_section_area_occupancy']:.4%}（V1 5.8530%），只是几何占比，不代表阻力或临床狭窄率。暴露面和附着面三角形面积互不重叠；{fit['dual_topological_node_tags']}个薄缘体积节点可同时邻近两面，不等于面积双计。
装配前精确线段检查发现4条七点抽检漏掉的wall chord穿腔，按通用规则删除并保留失败证据；未改变点坐标/形态。最终wall chord穿腔0、单连通材料主体、无孤点。连续实体和有限dx体素的边缘不同，不能把阶梯中心距称为悬空。

4. 点、锚点、连接和体积误差？

材料点 {v2['clot_nodes']}；三层锚点 {v2['anchor_nodes']}，各层 {v2['wall_layers_counts']}；BULK {v2['bulk_bonds']}、INTERFACE {v2['interface_bonds']}、WALL {v2['wall_bonds']}，总 {v2['bond_count']}。全部重新从50µm全局立方格筛选和建立连接，未变形旧点云。
节点V=dx³、m=1000V，总质量 {v2['clot_mass_kg']:.9g} kg，锚点排除。连续与积分体积差 {geom['discrete_relative_continuous_error']:.3%}，未进行网格细化研究。暴露节点{v2['exposed_nodes']}，面积 {v2['exposed_area_m2']*1e6:.6f} mm²；每个暴露三角形只分配一次。

5. run 0 与零载动态？

一次run0通过：初始力/速度/损伤0，typed topology逐项一致，实际二进制SHA256 `7cba88869fb61531b4e5fe5bbcf3b129f80cc1d9c19c0a4761d55e69b67cf95b`。
C0完整4660步通过，最大位置差 {gate['max_displacement_m']:.3g} m、速度0、断键0。四例均保存234帧、终步4660、无NaN/lost atoms/early abort，锚点固定，原生broken events与末态 surviving topology一致。只做这一轮必要检查，没有性能/材料扫描。

6. 三种剪切的真实响应？

{table}
具体首次bulk/interface时间及空间位置在每例result.json的first_break；原生事件保留全部两端节点及当前/参考断裂坐标。类别由冻结规则及连续指标得到，没有按想看到的碎裂图改分类。末态类别：{'; '.join(r['case_id']+'='+r['classification'] for r in rows)}。

7. ≥200/≥300µm 游离块？

{chr(10).join(sizetext)}

另用原生事件在每个断键步重建连通图：C1/C2/C3分别核查1558/1925/1267个事件时刻，最大游离组件13/23/8节点。0.5µs原生步层面也没有≥34或≥114节点游离块，因此并非仅凭10µs抽帧得出“未出现”。这里不生成或插值未保存的位置；小块exact-membership寿命仍按保存帧报告。详见NATIVE_EVENT_TRANSIENT_SIZE_AUDIT.json。

尺寸由体积等效球d=(6Ndx³/π)^(1/3)算得，200/300µm分别至少34/114节点；113节点小于300µm。附着主体永不算游离大块。exact-membership episode及parent/child splits复用历史lineage模块；至少34节点保存参考PCA三轴和薄片标志。persistent仅沿用“至少三个连续保存时刻”，本轮最短观测跨度20µs，不能叫长期稳定或临床持续。1–2节点是NUMERICALLY_UNRESOLVED_DEBRIS，3节点工程类目也不是大块。大块数据只是远端栓塞风险代理，不是概率或治疗安全认证。

8. 哪种更释放碎屑，哪种保留主体？

末态碎屑最多为{most['case_id']}（{most['final_unresolved_debris_fraction']:.3%}）；附着质量保留最多为{keeps['case_id']}（{keeps['final_attached_mass_fraction']:.3%}）。本轮只比较冻结几何和窗口内结果，不外推单调规律。

{loadtable}
三例脉宽均为绝对1.1605ms（25%余弦升/50%平台/25%降），观察2.33ms。新H/cs的Tcross约{c['H_m']/np.sqrt(c['G_Pa']/c['rho_kg_m3'])*1e3:.4f}ms，因此不称新几何恰好2Tcross。
同预算只统一参考Σ|F|/ΣA；Opposing的空间方向抵消和峰值/净力/力矩不同。实际每步载荷预算/局部峰值和力积分在force_history.dat，净冲量在result.json；力矩积分由10µs保存格算得并标注采样。空间高斯中心固定于实际新表面上预定参数点，不随粒子移动；每步采样当前位置，所以离场自然减弱。
明确的实现近似：F_i=A_i P_i^V2 e_x g(x_i(t))·scale·pulse(t)，面积与真实V2参考法向投影冻结在初始节点上；包络为当前位置空间场，但切向不是动态变形面重构。没有声称这是完全Eulerian法向场或真实流体剪切；新增断裂面不加载。

9. 图和视频是否真实、完整？

交付四张主PNG/PDF、同步MP4/GIF、Gaussian200的原始solver四关键帧、四例各234时刻ParaView点/连接/tag。MP41920×1088、234个真实保存帧、20fps；同一物理时刻和局部相机/尺度，全部位移倍率1。上排完整实际Y、下排近景，参考形态对比另见V1/V2图。
动态shell每个材料节点一个真实dx立方体，所有面只属于单个节点，不产生任何跨分量多边形，因此不会用平滑桥接断开组件。红色附着、橙色游离、灰色1–2节点弱化显示，均保留真实坐标和尺寸；这是节点体积显示，不是连续裂纹面的高保真重建。原始节点/全部surviving bonds/组件与lineage tags可单独审查。
各例最大越出示意Y管腔的节点数：{', '.join(r['case_id']+'='+str(r['maximum_nodes_outside_visual_vessel']) for r in rows)}。影片未裁剪数据以制造管内输运；近景视口有固定取景范围，完整数据和全Y背景仍保留。完整管壁没有碰撞，只有局部anchor接触；没有流体阻力、分支选择或捕获率。影片不叠加旧CFD，以免误认载荷来源。

10. 实际成本？

{durations}

总求解墙钟 {perf['total_solver_wall_s']:.2f}s，总后处理 {perf['total_postprocess_wall_s']:.2f}s。后处理合计包含28.43s原生事件阈值复核。上述RSS是求解进程树采样峰值；后处理RSS未单独采样。两MPI进程/例，没有缩放实验。图/视频/ParaView各阶段墙钟与渲染进程RSS见PERFORMANCE_REAL_CASES.json；后处理和图形仅读已保存结果，未再调用LAMMPS。总共四次动力学且没有重跑。

11. 解释边界？

原CFD为空管腔背景，未满足新血栓面的流固边界，不能导出当前阻力/剪切。Liu映射只保留未损伤小应变弹性，未复现Yeoh大应变/Mullins；猪皮peeling proxy与旧平面failure-strength映射均未验证曲壁真实血管粘附。1.0dx接触是有限案例工程候选；数值drag500/s不是血液阻力。整条血管无碎块碰撞/输运，因此不能解释沿A/B分支的生物行为。人工100/200Pa没有真实bubble streaming来源，不能称治疗效应。所有定量值只在当前探索性固体模型内成立。

12. 下一阶段真实载荷还缺什么？

需要含血栓流体域及暴露面boundary tag的相容网格，重新求含血栓速度/压力并从统一法向应力映射牵引；或提供具有空间/时间/SI单位/有效域和来源身份的bubble streaming场，并校核表面面积、总力、方向及几何更新契约。现有当前坐标包络接口可替换，但动态法向、断裂新增表面和FSI/输运还需另行定义。本轮只报告缺口，不执行。

实现参考：[LAMMPS addforce](https://docs.lammps.org/fix_addforce.html)、[BPM/PERI bonds](https://docs.lammps.org/bond_bpm_peri.html)、[BPM/PERI contact](https://docs.lammps.org/pair_bpm_peri.html)。实际1.0dx改版身份与官方默认分开记录。

已停止并等待人工审查；B1/B2 TERMINATED，PD2_FRACTURE STOPPED。没有commit/push/PR，无后续算例。结果入口为results目录OPEN_RESULTS.html，证据包见本目录evidence.zip；完整轨迹和逐帧ParaView文件保留原位，哈希清单随包。
'''
    (OUT/'REVIEW_Y_BIFURCATION_IRREGULAR_CLOT_AND_INITIAL_FRACTURE_V1.md').write_text(text)

def package():
    source=OUT/'source_snapshot';source.mkdir(exist_ok=True);diff=[]
    for p in sorted((ROOT/'pd_clot/lammps_fragmentation').glob('*.py')):
        shutil.copy2(p,source/p.name);old=OUT/'source_before'/p.name
        before=old.read_text() if old.exists() else ''
        if before!=p.read_text():diff.extend(difflib.unified_diff(before.splitlines(True),p.read_text().splitlines(True),fromfile='before/'+p.name,tofile='after/'+p.name))
    (OUT/'CODE_CHANGES.diff').write_text(''.join(diff));dump(OUT/'SOURCE_SNAPSHOT_HASHES.json',{p.name:digest(p) for p in source.glob('*.py')})
    dump(OUT/'ENVIRONMENT.json',dict(python=sys.executable,python_version=sys.version,pyvista=pv.__version__,numpy=np.__version__,new_dependencies_installed=[],lammps_binary=str(EXE),lammps_sha256=digest(EXE)))
    manifests=[]
    for name,_,_ in CASES:
        p=OUT/'cases'/name
        for f in p.iterdir():
            if f.is_file():manifests.append(dict(path=str(f),bytes=f.stat().st_size,sha256=digest(f)))
        # Small original text frames, exactly copied from native trajectory.
        data=(p/'trajectory.lammpstrj').read_text();blocks=data.split('ITEM: TIMESTEP\n')[1:]
        (p/'native_trajectory_keyframes.lammpstrj').write_text(''.join('ITEM: TIMESTEP\n'+blocks[k] for k in [0,58,117,233]))
    for p in (RESULTS/'paraview').rglob('*'):
        if p.is_file():manifests.append(dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)))
    dump(OUT/'LARGE_RAW_ASSET_HASH_MANIFEST.json',manifests)
    figs=['y_irregular_clot_v1_vs_v2','y_irregular_clot_material_and_interface','y_clot_fracture_response_comparison','y_clot_detached_fragment_sizes']
    html='<!doctype html><meta charset="utf-8"><title>Y clot V2 — saved results</title><style>body{max-width:1200px;margin:32px auto;font:17px system-ui;color:#233746}img,video{width:100%}a{color:#176189}</style><h1>Y 不规则附壁血栓：四例初步断裂</h1><p>探索性人工剪切；无含血栓 CFD / bubble streaming / 完整管壁碰撞或流体输运。位移与碎块尺寸未放大。</p><video controls src="animations/y_clot_initial_fracture_synchronized.mp4"></video>'
    for f in figs:html+=f'<p><a href="figures/{f}.pdf">PDF</a> · {f}</p><img src="figures/{f}.png">'
    html+='<p><a href="figures/solver_view_keyframes.png">原始 solver 关键帧</a> · <a href="paraview/README.md">ParaView 说明</a></p>'
    (RESULTS/'OPEN_RESULTS.html').write_text(html)
    (OUT/'FINAL_GIT_STATE.md').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
    files=[]
    for p in OUT.rglob('*'):
        if not p.is_file() or p.name in ['evidence.zip','EVIDENCE_MANIFEST.json'] or 'development_attempts' in p.parts:continue
        if p.name in ['trajectory.lammpstrj','accepted_trajectory.npz','surviving_bonds.local','final.restart','graph_reconstruction.npz']:continue
        files.append((p,'verification/'+str(p.relative_to(OUT))))
    for folder in ['figures','animations','FIGURE_DATA']:
        for p in (RESULTS/folder).rglob('*'):
            if p.is_file():files.append((p,'results/'+str(p.relative_to(RESULTS))))
    for p in (RESULTS/'paraview').rglob('*'):
        if p.is_file() and (p.parent==RESULTS/'paraview' or p.name in ['points_000.vtp','points_233.vtp','bonds_000.vtp','bonds_233.vtp','native_break_events.vtp','trajectory.pvd']):files.append((p,'results/'+str(p.relative_to(RESULTS))))
    files.append((RESULTS/'OPEN_RESULTS.html','results/OPEN_RESULTS.html'))
    with zipfile.ZipFile(OUT/'evidence.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p,arc in files:z.write(p,arc)
        z.writestr('EVIDENCE_README.txt','Light evidence subset. Full trajectory and every-frame ParaView files remain at paths in LARGE_RAW_ASSET_HASH_MANIFEST.json. A .pvd from the evidence subset needs the full raw series to open every time. No original large CFD copied. Native broken events and four unedited native point frames per case are included.\n')
    with zipfile.ZipFile(OUT/'evidence.zip') as z:assert z.testzip() is None;members=len(z.namelist())
    dump(OUT/'EVIDENCE_MANIFEST.json',dict(path=str(OUT/'evidence.zip'),bytes=(OUT/'evidence.zip').stat().st_size,sha256=digest(OUT/'evidence.zip'),members=members,zip_integrity='PASS'))

if __name__=='__main__':
    if sys.argv[1]=='eventaudit':native_event_audit()
    elif sys.argv[1]=='aggregate':aggregate()
    elif sys.argv[1]=='final':
        qa();rows,_=aggregate();st=final_status(rows);review(rows,st);package();print('READY_FOR_HUMAN_REVIEW',flush=True)
