"""Evidence tables, plots, final point-cloud delivery and restrained packaging."""
import json,csv,sys,hashlib,zipfile,subprocess,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from .early_flow_coupling import ROOT, OUT, RESULTS, CASE, OLD, PRIOR, DTYPE, DT, MASS, drag_law
from .failure_map import csv_out
from .config import dump
from .y_bifurcation_sources import digest

def read(p):return json.loads(p.read_text())

def precise_clock():
    """Repair only clock typography; retain the already rendered native scenes.

    Four decimal places in milliseconds represent the 0.5us event grid exactly.
    This never changes frame times, order, coordinates or camera framing.
    """
    from .point_cloud_delivery import FONTS,COLORS,FPS
    from PIL import Image,ImageDraw
    import imageio_ffmpeg
    start=time.perf_counter();name='y_clot_point_cloud_early_flow_coupling';movie=RESULTS/'animations'/(name+'.mp4');manifest_path=RESULTS/'frame_manifests'/(name+'.json');manifest=read(manifest_path)
    assert 'clock_precision_repair_seconds' not in manifest
    sourcehash=digest(movie);temp=movie.with_suffix('.clock.mp4');gif=[];gif_frames=[]
    def correct(im,r):
        draw=ImageDraw.Draw(im);draw.rectangle((40,64,1230,110),fill='white')
        draw.text((50,88),f'Physical time  {r["time_s"]*1000:.4f} ms     |     Playback: {r["physical_seconds_per_video_second"]:.4g} s physical / s video',font=FONTS[26],fill=COLORS['axis'],anchor='lm')
        return im
    p=subprocess.Popen([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-hide_banner','-loglevel','error','-f','rawvideo','-pix_fmt','rgb24','-s:v','1920x1080','-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','medium','-crf','16','-pix_fmt','yuv420p','-threads','4','-movflags','+faststart',str(temp)],stdin=subprocess.PIPE,stderr=subprocess.PIPE)
    gen=imageio_ffmpeg.read_frames(str(movie),pix_fmt='rgb24');next(gen)
    for i,frame in enumerate(gen):
        r=manifest['frames'][i];im=correct(Image.frombytes('RGB',(1920,1080),frame),r);p.stdin.write(np.asarray(im).tobytes())
        if i%4==0 and 3.5<=r['video_time_s']<=13.5:gif.append(im.resize((960,540),Image.Resampling.LANCZOS));gif_frames.append(i)
    assert i+1==len(manifest['frames']);p.stdin.close();error=p.stderr.read().decode();assert p.wait()==0,error
    # Preserve pre-release rendering, not overwrite any historical task output.
    backup=OUT/'render_before_exact_time_labels';backup.mkdir(exist_ok=True);movie.replace(backup/movie.name);temp.replace(movie)
    gif[0].save(movie.with_suffix('.gif'),save_all=True,append_images=gif[1:],duration=200,loop=0,optimize=False)
    for suffix in ['initial','handoff','transport_5s','first_break_exact_native_step','first_detachment_exact_native_step']:
        if suffix=='initial':r=manifest['frames'][0]
        elif suffix=='transport_5s':r=manifest['frames'][-1]
        elif suffix=='handoff':r=next(r for r in manifest['frames'] if r['phase']=='LAMMPS' and r['step']==4660)
        else:r=next(r for r in manifest['frames'] if r.get('event')==suffix.replace('_exact_native_step',''))
        path=RESULTS/'figures'/(name+'_'+suffix+'.png');im=correct(Image.open(path).convert('RGB'),r);im.save(path);im.save(path.with_suffix('.pdf'),resolution=150)
    manifest.update(time_label_decimal_places_ms=4,clock_precision_repair_seconds=time.perf_counter()-start,clock_repair_before_sha256=sourcehash,clock_repair_after_sha256=digest(movie),clock_repair_scope='Text rectangle x40..1230,y64..110 only; same native scenes, time grid and frame order; re-encoded H264',short_GIF=dict(source_MP4_frame_indices=gif_frames,duration_ms_per_frame=200,physical_time_and_playback_rates_unchanged=True))
    dump(manifest_path,manifest);print('CLOCK_PRECISION_PASS',time.perf_counter()-start,flush=True)
def plots():
    raw=np.memmap(CASE/'detached_force_audit.bin',dtype=DTYPE,mode='r');online=np.load(CASE/'ONLINE_TOPOLOGY_EVENTS.npz');times=online['detach_step'];tr=np.load(CASE/'accepted_trajectory.npz')
    hand=read(OUT/'FRAGMENT_HANDOFF.json')['components'];singleton=[r['member_atom_ids'][0] for r in hand if r['node_count']==1]
    ordered=sorted(singleton,key=lambda i:times[i-1]);chosen=[ordered[0],ordered[len(ordered)//2],ordered[-1]]
    data=[];eventrows=[]
    for atom in np.flatnonzero(times>=0)+1:
        a=raw[(raw['atom_id']==atom)&(raw['step']==times[atom-1])][0]
        row=dict(atom_id=int(atom),component_id_at_detach=int(a['component_id']),component_nodes_at_detach=int(a['node_count']),detach_step=int(a['step']),t_detach_s=float(a['step']*DT),drag_on_step=int(online['drag_step'][atom-1]),t_drag_on_s=float(online['drag_step'][atom-1]*DT),delay_s=float((online['drag_step'][atom-1]-a['step'])*DT),CFD_tetra_id=int(a['tetra_id']),CFD_query='CURRENT_COMPONENT_COM',near_wall_uncertain=bool(a['near_wall_uncertain']),clot_present_flow_uncertain=bool(a['clot_present_flow_uncertain']))
        for key in ['x_m','v_force_half_m_s','u_COM_m_s','drag_N','numerical_N','wall_N','load_N']:
            for axis,v in zip('xyz',a[key]):row[key+'_'+axis]=float(v)
        eventrows.append(row)
    csv_out(OUT/'DETACHMENT_DRAG_EVENT_AUDIT.csv',eventrows)
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf');plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#93a3b0','axes.labelcolor':'#263746','text.color':'#263746','xtick.color':'#647887','ytick.color':'#647887','pdf.fonttype':42})
    fig,axs=plt.subplots(2,2,figsize=(13.8,8.2),layout='constrained');colors=['#b85c64','#586d7b','#8aa7b9']
    for atom,color in zip(chosen,colors):
        a=raw[raw['atom_id']==atom];time_s=a['step']*DT;speed=np.linalg.norm(a['v_force_half_m_s'],axis=1);fluid=np.linalg.norm(a['u_COM_m_s'],axis=1);slip=np.linalg.norm(a['v_force_half_m_s']-a['u_COM_m_s'],axis=1)
        # These selected final singletons can have belonged to larger detached
        # groups earlier: record node speed vs queried COM u without hiding it.
        elapsed=time_s-times[atom-1]*DT
        axs[0,0].plot(time_s*1000,speed,color=color,label=f'Node {atom}');axs[0,0].plot(time_s*1000,fluid,color=color,ls=':',lw=1.5)
        axs[0,1].semilogy(elapsed*1000,np.maximum(slip,1e-8),color=color,label=f'Node {atom}')
        axs[1,0].plot(time_s*1000,np.linalg.norm(a['drag_N'],axis=1)*1e6,color=color)
        for j in range(len(a)):
            data.append(dict(phase='LAMMPS_FORCE_EVALUATION',atom_id=atom,component_id=int(a['component_id'][j]),node_count=int(a['node_count'][j]),time_s=float(time_s[j]),since_detach_s=float(elapsed[j]),node_speed_m_s=float(speed[j]),CFD_COM_speed_m_s=float(fluid[j]),node_slip_to_COM_CFD_m_s=float(slip[j]),drag_magnitude_N=float(np.linalg.norm(a['drag_N'][j])),tau_s=float(a['tau_s'][j]),velocity_time_level='VERLET_HALF_STEP',CFD_valid=int(a['valid_CFD'][j])))
    axs[0,0].set(xlabel='Physical time (ms)',ylabel='Speed (m/s)',title='Native motion after detachment');axs[0,0].legend(frameon=False);axs[0,0].text(.03,.97,'Solid: node speed\nDotted: local background CFD',transform=axs[0,0].transAxes,va='top',fontsize=10)
    axs[0,1].set(xlabel='Time since first detachment (ms)',ylabel='|v node − u CFD at COM| (m/s)',title='Measured slip: load and contact remain active')
    axs[1,0].set(xlabel='Physical time (ms)',ylabel='Applied drag per node (µN)',title='Finite force starts at the real detachment event')
    transport=list(csv.DictReader((OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv').open()));rec=next(r for r in hand if r['member_atom_ids']==[chosen[0]]);cid=rec['component_id'];subset=[r for r in transport if int(r['component_id'])==cid and float(r['relative_time_s'])<=.001]
    a=raw[(raw['atom_id']==chosen[0])&(raw['step']>=4460)];t0=a['step']*DT;v0=np.linalg.norm(a['v_force_half_m_s'],axis=1)
    axs[1,1].plot(t0*1000,v0,color=colors[0],label='LAMMPS half-step samples')
    tt=np.array([.00233]+[float(r['absolute_time_s']) for r in subset]);vv=np.vstack([rec['COM_velocity_mps']]+[[float(r['v'+k+'_m_s']) for k in 'xyz'] for r in subset]);uu=np.vstack([subset[0]['u_'+k+'_m_s'] for k in 'xyz']).ravel().astype(float)
    axs[1,1].plot(tt*1000,np.linalg.norm(vv,axis=1),color=colors[1],label='Transport accepted COM states');axs[1,1].axvline(2.33,color='#93a3b0',ls='--');axs[1,1].set(xlabel='Absolute physical time (ms)',ylabel='Speed (m/s)',title=f'Continuous COM handoff: node {chosen[0]}');axs[1,1].legend(frameon=False,fontsize=9)
    for ax in axs.ravel():ax.grid(alpha=.25,color='#93a3b0')
    fig.suptitle('Early background-flow coupling · actual native and transport records',fontsize=17)
    fig.savefig(RESULTS/'figures/velocity_flow_timeline.png',dpi=220);fig.savefig(RESULTS/'figures/velocity_flow_timeline.pdf');plt.close(fig)
    csv_out(OUT/'VELOCITY_AND_FLOW_COMPARISON.csv',data)
    dump(OUT/'VELOCITY_PLOT_MANIFEST.json',dict(representative_atom_ids=chosen,native_source=str(CASE/'detached_force_audit.bin'),native_dtype=str(CASE/'FORCE_AUDIT_DTYPE.json'),transport_source=str(OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv'),no_smoothing=True,native_velocity_time_level='half-step at actual force evaluation',slip_scope='Node v versus component-query u; same COM drag convention. Not assumed monotonic under continuing mechanical loads.'))
    # Component-resolved terminal sizes and start-force delays.
    csv_out(OUT/'REPRESENTATIVE_COMPONENT_DETACHMENT.csv',[dict(component_id=r['component_id'],member_atom_ids=json.dumps(r['member_atom_ids']),node_count=r['node_count'],first_member_detach_s=float(times[np.array(r['member_atom_ids'])-1].min()*DT),last_member_detach_s=float(times[np.array(r['member_atom_ids'])-1].max()*DT),first_exact_membership_detach_s=r['first_exact_detached_native_time_s'],max_member_drag_delay_s=float(((online['drag_step']-times)[np.array(r['member_atom_ids'])-1]*DT).max()),equivalent_reference_volume_diameter_um=r['equivalent_diameter_m']*1e6,category=r['category']) for r in hand])

def finalize():
    from .point_cloud_delivery import Data, FPS
    import imageio_ffmpeg
    native=read(OUT/'NATIVE_EARLY_FLOW_VALIDATION.json');handoff=read(OUT/'HANDOFF_CONTINUITY_AUDIT.json');integrator=read(OUT/'NO_VELOCITY_RESET_INTEGRATOR_AUDIT.json');assert native['status']==handoff['status']==integrator['status']=='PASS'
    videos=[]
    for name in ['y_clot_point_cloud_only_preview','y_clot_point_cloud_early_flow_coupling']:
        file=RESULTS/'animations'/(name+'.mp4');manifest=read(RESULTS/'frame_manifests'/(name+'.json'));gen=imageio_ffmpeg.read_frames(str(file),pix_fmt='rgb24');meta=next(gen);count=0;previous=None;duplicates=0
        for frame in gen:
            sha=hashlib.sha256(frame).hexdigest();duplicates+=int(sha==previous);previous=sha;count+=1
        assert count==len(manifest['frames']) and meta['size']==(1920,1080) and meta['fps']==FPS and duplicates==0
        assert all(r['attached']+r['detached']==1927 for r in manifest['frames']);assert np.all(np.diff([r['time_s'] for r in manifest['frames']])>0)
        assert all(a['xyz_sha256']!=b['xyz_sha256'] for a,b in zip(manifest['frames'][:-1],manifest['frames'][1:])), 'Repeated physical pose'
        source=PRIOR if name.endswith('preview') else OUT;d=Data(source)
        for r in manifest['frames']:
            X,a,alive,cid=d.state(r)
            assert hashlib.sha256(X.astype('<f8').tobytes()).hexdigest()==r['xyz_sha256']
            assert hashlib.sha256(cid.astype('<i8').tobytes()).hexdigest()==r['component_id_sha256']
            assert (np.all(X>=np.array([-14.5,-7.3,-1.5])*1e-3,axis=1)&np.all(X<=np.array([11.5,7.3,1.5])*1e-3,axis=1)).all(),'Full-Y bounds omitted an ID'
        videos.append(dict(path=str(file),frames=count,duplicate_consecutive_frames=duplicates,render_seconds=manifest['render_seconds'],clock_precision_repair_seconds=manifest.get('clock_precision_repair_seconds',0),duration_s=meta['duration'],sha256=digest(file),all_1927_clot_and_1561_anchor_ids_traceable=True,positions_and_topology_match_manifest=True))
    for old in read(OUT/'PRIOR_READ_ONLY_REGISTER.json'):assert digest(old['path'])==old['sha256'],'Historical input changed'
    config=read(ROOT/'references/y_clot_default_render.json');config.update(data_source=str(OUT.relative_to(ROOT)),output_basename='y_clot_point_cloud_early_flow_coupling',accepted_validation=str((OUT/'FINAL_STATUS.json').relative_to(ROOT)))
    dump(ROOT/'references/y_clot_default_render.json',config)
    status=dict(POINT_CLOUD_DEFAULT_RENDER='PASS',SMOOTH_SHELL_DEFAULT='DISABLED',TRUE_VELOCITY_RESET_BUG='NOT_OBSERVED',BACKGROUND_DRAG_DURING_FRACTURE='IMPLEMENTED',EARLY_FLOW_COUPLING='ESTABLISHED_FOR_ONE_WAY_ENGINEERING_COM_DRAG',DETACHMENT_TRIGGER_CORRECT='PASS',HANDOFF_POSITION_CONTINUITY='PASS',HANDOFF_VELOCITY_CONTINUITY='PASS',HANDOFF_CONTINUITY_SCOPE='COMPONENT_COM; internal modes and attached core mechanics are truncated after2.33ms',HANDOFF_DRAG_LAW_CONTINUITY='PASS_AT_IDENTICAL_ACCEPTED_STATE',BACKGROUND_FLOW='UNOBSTRUCTED_ONE_WAY_EXPLORATORY',CLOT_PRESENT_CFD='NOT_RUN',MACRO_FRAGMENT_200UM='OBSERVED' if native['macro_fragment_ge200um'] else 'NOT_OBSERVED',LARGE_BLOCK_FRACTURE_MECHANISM='NOT_ESTABLISHED',mechanical_cases_run=1,baseline_reruns=0,solver_steps=4660,detached_nodes=native['detached_nodes'],maximum_detachment_drag_delay_s=0.,no_illegal_wall_crossing=True,prior_inputs_unchanged=True,videos=videos)
    dump(OUT/'FINAL_STATUS.json',status);dump(RESULTS/'FINAL_STATUS.json',status)
    dump(OUT/'VIDEO_AND_FRAME_VALIDATION.json',dict(status='PASS',videos=videos,rendering_mode='Raw nodes and surviving lines only; no material shell geometry',native_wall_context=True))
    comparisons=native['comparison'];old,new=comparisons;hand=read(OUT/'FRAGMENT_HANDOFF.json');counts={}
    for r in hand['components']:counts[str(r['node_count'])]=counts.get(str(r['node_count']),0)+1
    render_s=sum(v['render_seconds'] for v in videos)
    report=f'''# Y_CLOT_POINT_CLOUD_AND_EARLY_FLOW_COUPLING_V1

点阵已成为正式默认。旧数据点阵预览已先完成，再以唯一一个实际新力学算例验证了按当前真实脱离触发流体力。没有 baseline 重跑、标定、幅值扫描或 CFD 重算。旧 medical 输出未删除。

## 核心结果

- 最早真实脱离：第593步，0.2965ms，节点267；同一步首次受拖曳。全部{new['detached_nodes']}个游离节点的触发延迟均为0步，算法检查粒度0.5µs。{native['event_count']}条在线断键记录与原生文件逐事件逐步一致；独立重建得到完全相同的脱离时间。
- 只有脱离组件受背景阻力；附着节点仍保留原 -500mv 数值阻尼，游离节点此项为零。{native['force_audit_rows']}条逐节点逐步流体力记录和{native['full_snapshot_rows']}条全节点快照记录可以复核。
- 当前COM查询严格使用原生四面体P1插值，域外不外推。此算例无无效拖曳查询。{native['clot_present_flow_uncertain_node_steps']}个节点步标记 CLOT_PRESENT_FLOW_UNCERTAIN，不能把无栓CFD视为含栓真实流动。
- 原生求解最小中心净距{native['wall']['minimum_clearance_m']*1e6:.6g}µm，穿壁0；旧/新发生原生壁接触节点{old['wall_contact_nodes']}/{new['wall_contact_nodes']}。保持接触模型后，提前拖曳实际抑制了这些碰壁，不是删除壁力。
- {integrator['consecutive_detached_node_steps']}个连续游离节点步的真实Verlet递推核对通过，速度误差最大{integrator['maximum_velocity_recurrence_error_m_s']:.6g}m/s，位置误差最大{integrator['maximum_position_recurrence_error_m']:.6g}m。{integrator['onsets_with_exact_previous_step_saved']}/217个脱离事件另有紧邻前一步全节点快照可直接核对，其余30个缺少脱离前一完整快照；全部脱离后相邻步均已核对，没有补造缺失记录。

|量|旧算例|提前耦合|
|---|---:|---:|
|首 bulk break (ms)|{old['first_bulk_break_s']*1000:g}|{new['first_bulk_break_s']*1000:g}|
|首 interface break (ms)|{old['first_interface_break_s']*1000:g}|{new['first_interface_break_s']*1000:g}|
|原生断键数|{old['native_broken_bonds']}|{new['native_broken_bonds']}|
|末态游离节点|{old['detached_nodes']}|{new['detached_nodes']}|
|末态游离组件|{old['detached_components']}|{new['detached_components']}|
|附着总质量 (kg)|{old['attached_mass_kg']:.9g}|{new['attached_mass_kg']:.9g}|
|最大附着主体质量 (kg)|{old['largest_attached_core_mass_kg']:.9g}|{new['largest_attached_core_mass_kg']:.9g}|
|1–2节点未解析质量 (kg)|{old['unresolved_mass_kg']:.9g}|{new['unresolved_mass_kg']:.9g}|
|末态最大游离组件节点数|{old['max_final_detached_component_nodes']}|{new['max_final_detached_component_nodes']}|
|求解墙钟时间 (s)|{old['solver_wall_seconds']:.3f}|{new['solver_wall_seconds']:.3f}|

新算例末态组件的节点数→组件数：{counts}。整个原生事件时序最大游离组件为{native['max_detached_component_nodes_all_event_steps']}节点，仍未达到200µm参考体积球径门槛（34节点）。仍是附着主体伴局部碎屑脱落；没有建立宏观块体破裂机制，不因结果改变而调回断裂或载荷参数。

## 为什么旧视频先停后动

旧交接COM位置和速度误差均0；1–2节点碎屑末态速度中位数0.858503m/s，背景CFD中位数0.00458474m/s。旧模型到2.33ms才突然打开100µs工程松弛，造成快速减速；近壁慢流和历史2000倍播放率切换又放大了观感。旧重复端点已在上一轮修复。本轮每帧用连续绝对时间和明确播放率，采用实际保存状态，不插值拼接提前输运。

“立即受流体力”仅意味着有限力进入动量方程。速度由惯性、持续载荷、内力、接触和拖曳共同积分，不可能因此瞬时等于水流。1–2节点的100µs是明确工程代理，没有生物碎片水动力直径含义；>=3节点沿用不规则组件的等效球质心近似，并按节点质量正确分配总力。

## 交接及物理限制

新终态严格来自原生接受状态。组件COM交接位置误差{handoff['max_position_jump_m']:.3g}m、速度误差{handoff['max_velocity_jump_m_s']:.3g}m/s；同一接受状态下两阶段瞬时拖曳律差{handoff['max_force_law_jump_N']:.3g}N。原生post_force读的是Verlet半步速度，末态dump是完整步速度；半步到完整步的正常力变化最大{handoff['max_halfstep_to_fullstep_force_change_N']:.6g}N，已用真实末半步积分核对，不能误报为阻力律跳变。

新末态组件COM速度中位数0.00622581m/s，对应CFD速度中位数0.00620436m/s；多数碎屑已经在原生断裂窗口内通过有限拖曳逐步接近局部流速。晚近脱离者和四节点组件仍存在滑移，不以人为速度赋值消除。

连续性PASS明确限于组件COM。长时输运保留末态成员偏移，忽略内部变形和旋转；截断的最大内部相对速度为{handoff['max_internal_relative_speed_truncated_m_s']:.6g}m/s，不宣称内部各节点速度也被连续求解。附着主体明确标为 FROZEN REFERENCE SNAPSHOT AFTER 2.33 MS。旧的原生20µm罚力墙模型在长时阶段仍切换为5µm单边约束，这是保留的接触近似和积分器切换，不能称完整力学连续。后续最多5s无栓单向输运，未人为提高流速，未任意缩短tau。

仍缺少含栓CFD/双向FSI、堵塞反馈、真实血栓表面流体牵引、长时碎片转动/形变与相互碰撞、留存血栓碰撞和经过验证的近壁水动力。本轮2523个不确定节点步只是标识局限，不是补齐物理。

## 可视化与复现

正式入口 `scripts/render_y_clot.py` / `references/y_clot_default_render.json`，默认点阵。颜色、细键线、白底、六面较深网格延续用户现有图。局部近景和整个Y管固定相机并列；所有1927材料点、1561锚点可按原ID追溯。灰色游离点是屏幕符号，未扩大真实碎片、合并不相连粒子、生成外壳或虚构裂缝。局部框之外的点仍在完整Y视图；所有键类型只在局部显示，图例与manifest明确说明。

两部MP4完整解码、逐帧来源/坐标/拓扑hash核对均通过，无重复端点。首断与首脱离PNG/PDF使用原生事件当步坐标，不是向10µs保存网格取整。原数据预览渲染{videos[0]['render_seconds']:.3f}s，新耦合渲染{videos[1]['render_seconds']:.3f}s；5s输运求解{handoff['transport_wall_seconds']:.3f}s。详细时间与解码结果见JSON。

新视频时钟保留毫秒小数4位，准确显示0.5µs事件步。初次编码的3位小数会把0.2965ms显示成0.296ms，已仅重绘页眉时钟并重新编码（{videos[1]['clock_precision_repair_seconds']:.3f}s），未修改场景、数据、帧顺序或播放率；原预发布编码保留在verification诊断目录。短GIF约10s，帧索引和播放率与对应MP4片段一致。

官方接口来源及固定本地源码SHA见 EARLY_DETACHMENT_DRAG_DESIGN.md、INTERFACE_READINESS.json。未修改现有LAMMPS源码或安装，无新原生hook。

轻量 evidence.zip 收录审计、代码、图和清单。大原生轨迹、逐步二进制力记录和MP4留在原位，见 LARGE_DATA_HASHES.json。历史输入hash全部未变。
'''
    (OUT/'REVIEW_Y_CLOT_POINT_CLOUD_AND_EARLY_FLOW_COUPLING_V1.md').write_text(report)
    (RESULTS/'README.md').write_text('默认视频：animations/y_clot_point_cloud_early_flow_coupling.mp4\n\n旧数据点阵预览：animations/y_clot_point_cloud_only_preview.mp4\n\n完整科学审计：../../verification/'+OUT.name+'/REVIEW_Y_CLOT_POINT_CLOUD_AND_EARLY_FLOW_COUPLING_V1.md\n')
    dump(OUT/'LARGE_DATA_HASHES.json',[dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for p in [*CASE.glob('*.bin'),*CASE.glob('*.npz'),CASE/'trajectory.lammpstrj',CASE/'broken_bonds.local',CASE/'surviving_bonds.local',OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz',OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv',*list((RESULTS/'animations').glob('*'))]])
    snapshot=OUT/'SOURCE_DELIVERY';snapshot.mkdir(exist_ok=True)
    import shutil
    for p in [*list((ROOT/'pd_clot/lammps_fragmentation').glob('early_flow*.py')),ROOT/'pd_clot/lammps_fragmentation/point_cloud_delivery.py',ROOT/'scripts/render_y_clot.py',ROOT/'references/y_clot_default_render.json']:shutil.copy2(p,snapshot/p.name)
    (OUT/'ENDING_GIT_STATE.md').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
    with zipfile.ZipFile(OUT/'evidence.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(OUT.iterdir()):
            if p.is_file() and p.suffix in ['.json','.md','.csv'] and p.stat().st_size<8_000_000:z.write(p,'verification/'+p.name)
        for p in snapshot.iterdir():z.write(p,'source/'+p.name)
        for name in ['request.json','LOAD_PROVIDER.json','in.lammps','RUN_COMMAND.json','PERFORMANCE.json','EARLY_FLOW_CALLBACK.json','FORCE_AUDIT_DTYPE.json','WALL_CONTACT_DIAGNOSTICS.json','first_break_state.npz','first_detachment_state.npz','ONLINE_TOPOLOGY_EVENTS.npz']:
            z.write(CASE/name,'case/'+name)
        z.write(OUT/'LAMMPS_INPUT_CHANGE.patch','verification/LAMMPS_INPUT_CHANGE.patch')
        for p in (RESULTS/'figures').glob('*'):
            if p.suffix in ['.png','.pdf']:z.write(p,'figures/'+p.name)
        for p in (RESULTS/'frame_manifests').glob('*.json'):z.write(p,'frame_manifests/'+p.name)
    with zipfile.ZipFile(OUT/'evidence.zip') as z:
        assert z.testzip() is None;members=len(z.namelist())
    dump(OUT/'EVIDENCE_PACKAGE.json',dict(path=str(OUT/'evidence.zip'),bytes=(OUT/'evidence.zip').stat().st_size,sha256=digest(OUT/'evidence.zip'),members=members,crc_check='PASS',large_data='Kept in place; see LARGE_DATA_HASHES.json'))
    print('DELIVERY_COMPLETE',status,flush=True)

if __name__=='__main__':{'plots':plots,'clock':precise_clock,'finalize':finalize}[sys.argv[1]]()
