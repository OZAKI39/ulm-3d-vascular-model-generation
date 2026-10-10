"""Offline evidence, figures and delivery for the continuous-handoff audit."""
import csv,json,hashlib,subprocess,time,zipfile
from pathlib import Path
import numpy as np
import pyvista as pv
import imageio.v2 as imageio
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .continuous_handoff_audit import OUT,RESULTS,PRIOR,PRIOR_RESULTS,CASE,quant
from .config import ROOT,dump
from .failure_map import csv_out
from .y_bifurcation_sources import YSources,digest
from .medical_fragment_renderer import component_shell

def read(name):return json.loads((OUT/name).read_text())
def rows(name):return list(csv.DictReader((OUT/name).open()))
def col(data,key):return np.array([float(r[key]) for r in data])
def savefig(fig,name):
    for suffix in ['png','pdf']:fig.savefig(RESULTS/'figures'/(name+'.'+suffix),dpi=180,bbox_inches='tight')
    plt.close(fig)

def complete_handoff_table():
    tr=np.load(CASE/'accepted_trajectory.npz');trans=np.load(PRIOR/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz');hand=json.loads((PRIOR/'FRAGMENT_HANDOFF.json').read_text())['components'];table=rows('handoff_velocity_audit.csv');by={int(r['component_id']):r for r in table}
    for j,h in enumerate(hand):
        cid=h['component_id'];ids=np.asarray(h['member_atom_ids'])-1;x=tr['position_m'][-1,ids].mean(0);v=tr['velocity_m_s'][-1,ids].mean(0);r=by[cid]
        for prefix,a in [('native_handoff_COM_m',x),('native_handoff_velocity_m_s',v),('transport_t0_COM_m',trans['position_m'][0,j]),('transport_t0_velocity_m_s',trans['velocity_m_s'][0,j])]:
            for axis,value in zip('xyz',a):r[prefix+'_'+axis]=float(value)
        r['handoff_JSON_position_error_m']=float(np.linalg.norm(np.asarray(h['COM_xyz_m'])-x))
        assert r['handoff_JSON_position_error_m']==0 and np.array_equal(v,trans['velocity_m_s'][0,j])
        assert abs(h['mass_kg']-len(ids)*1.25e-10)<1e-22
    csv_out(OUT/'handoff_velocity_audit.csv',table)
    a=read('HANDOFF_CONTINUITY_AUDIT.json');a['sampling_after_500us_elapsed_s']=quant(col(table,'sample_after_500us_s'));a['speed_ratio_sampling_caveat']='First real accepted state at or after500us, not a common-time500us reconstruction; inspect each elapsed time in CSV.';a['max_handoff_JSON_position_error_m']=max(float(r['handoff_JSON_position_error_m']) for r in table);dump(OUT/'HANDOFF_CONTINUITY_AUDIT.json',a)

def plots():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'#fffdf9'})
    summary=rows('handoff_velocity_audit.csv');early=rows('early_transport_native_steps_with_CFD_and_drag.csv');hist=rows('last_five_LAMMPS_COM_frames.csv');m=read('RENDER_MANIFEST.json');fig,axs=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    for r in summary:
        cid=r['component_id'];sel=[e for e in early if e['component_id']==cid];old=[e for e in hist if e['component_id']==cid];six=int(r['node_count'])>=3;color='#8c1e35' if six else '#819b9c';alpha=1 if six else .13
        t=np.r_[(col(old,'physical_time_s')-.00233)*1000,0,col(sel,'transport_elapsed_s')*1000];v=np.r_[np.linalg.norm(np.array([[float(e['COM_v'+a+'_m_s']) for a in 'xyz'] for e in old]),axis=1),float(r['initial_speed_m_s']),col(sel,'speed_m_s')]
        axs[0,0].plot(t,v,color=color,alpha=alpha,lw=1.8 if six else .6)
        if six:
            axs[0,1].plot(col(sel,'transport_elapsed_s')*1000,col(sel,'speed_m_s'),'o-',ms=3,color=color,label='6-node component: accepted states')
            u=np.linalg.norm(np.array([[float(e['u_CFD_start_m_s_'+a]) for a in 'xyz'] for e in sel]),axis=1);axs[0,1].plot(col(sel,'transport_elapsed_s')*1000,u,'o--',ms=3,color='#347d81',label='CFD sampled at actual pre-step COM')
    axs[0,0].set(xscale='symlog',yscale='log',xlabel='Time relative to handoff (ms; symlog)',ylabel='COM speed (m/s)',title='Native final cohorts and real accepted transport states');axs[0,0].set_xscale('symlog',linthresh=.05);axs[0,0].set_xticks([-.1,0,.1,1,10],['−0.1','0','0.1','1','10']);axs[0,0].axvline(0,color='black',lw=.7);axs[0,0].text(.02,.03,'Gray: 207 mass packets; red: 6-node component\nLines connect recorded samples; no new integration.',transform=axs[0,0].transAxes,fontsize=9)
    axs[0,1].set(xscale='log',yscale='log',xlabel='Transport elapsed time (ms)',ylabel='Speed (m/s)',title='Initial speed is inherited, then drag decelerates it');axs[0,1].legend(fontsize=8)
    axs[1,0].scatter(col(summary,'initial_CFD_speed_m_s'),col(summary,'initial_speed_m_s'),s=12,color='#586f73',alpha=.6);six=[r for r in summary if int(r['node_count'])>=3];axs[1,0].scatter(col(six,'initial_CFD_speed_m_s'),col(six,'initial_speed_m_s'),s=60,color='#8c1e35',label='6-node component');axs[1,0].set(xscale='log',yscale='log',xlabel='Local background CFD speed (m/s)',ylabel='Handoff COM speed (m/s)',title='Background CFD is much slower than native debris');axs[1,0].legend(fontsize=8)
    rr=m['frames'];ix=np.arange(len(rr));axs[1,1].plot((ix-233)/m['fps'],[r['physical_seconds_per_playback_second'] for r in rr],color='#8c1e35',label='New: actual manifest');axs[1,1].plot([-5.85,-.05,0,6.25],[.0004,.0004,.8,.8],color='#347d81',ls='--',label='Prior: typical rates around handoff');axs[1,1].axvline(0,color='black',lw=.7);axs[1,1].set(yscale='log',xlabel='Playback time relative to handoff (s)',ylabel='Physical seconds / playback second',title='Playback rate is explicit; physical clock never resets');axs[1,1].legend(fontsize=8)
    fig.suptitle('Handoff audit: zero position jump and zero velocity jump',fontsize=15);savefig(fig,'handoff_speed_and_playback_audit')
    event=rows('native_break_spatial_distribution.csv');inode=rows('interface_attachment_node_map.csv');series=rows('native_component_timeline.csv');diag=read('FRAGMENTATION_MECHANISM_DIAGNOSIS.json');fig=plt.figure(figsize=(14,9),layout='constrained');axes=[]
    for j,btype in enumerate([1,2]):
        ax=fig.add_subplot(2,3,j+1,projection='3d');sel=[r for r in event if int(r['bond_type'])==btype];sc=ax.scatter(col(sel,'reference_x_m')*1000,col(sel,'reference_y_m')*1000,col(sel,'reference_z_m')*1000,c=col(sel,'physical_time_s')*1000,s=2,alpha=.55,cmap='magma',rasterized=True);ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=('Bulk' if btype==1 else 'Interface')+' breaks at reference bond midpoints');fig.colorbar(sc,ax=ax,shrink=.65,orientation='horizontal',pad=.1,label='Native break time (ms)');ax.view_init(22,-65)
    ax=fig.add_subplot(2,3,3);st=diag['depth_strata'];ax.bar(np.arange(5),[r['broken_fraction_in_stratum']*100 for r in st],color='#8c1e35');ax.set_xticks(np.arange(5),['0–50','50–100','100–150','150–250','>250'],rotation=25);ax.set(yscale='log',xlabel='Initial exposed-surface distance (µm)',ylabel='Broken / initial bulk bonds in stratum (%)',title='99.28% of broken bulk bonds lie within 100 µm')
    ax=fig.add_subplot(2,3,4);sc=ax.scatter(col(inode,'x_m')*1000,col(inode,'y_m')*1000,c=col(inode,'broken_fraction'),s=10,cmap='magma',vmin=0,vmax=1,rasterized=True);ax.set(xlabel='Reference x (mm)',ylabel='Reference y (mm)',title='Interface node map (x–y projection)');fig.colorbar(sc,ax=ax,label='Fraction of direct interface bonds lost')
    ax=fig.add_subplot(2,3,5);ax.step(col(series,'time_s')*1000,col(series,'largest_detached_nodes'),where='post',color='#8c1e35');ax.axhline(34,ls='--',color='#387c81',label='200 µm reference volume:34 nodes');ax.axhline(114,ls=':',color='#777777',label='300 µm reference volume:114 nodes');ax.set(xlabel='Native fracture time (ms)',ylabel='Largest detached component (nodes)',title='All 2030 native break times: peak 13 nodes');ax.legend(fontsize=8)
    ax=fig.add_subplot(2,3,6);ax.axis('off');ax.text(0,.98,'Retained connectivity\n\n1702 attached nodes:\n1701 core +1 anchored singleton\n70261 surviving core bulk bonds\n21220 surviving core interface bonds\n1282 distinct anchors connected to core\n\n45 central axis sections tested\nMinimum 184 surviving crossing bulk bonds\n\n291 exact-member detached lifetimes\n7 have no saved COM frame: explicitly unavailable\n\nNo established macroscopic peel front or cut surface.\nGraph connectivity is evidence, not unique causation.',va='top',fontsize=10)
    fig.suptitle('Native failure mechanism: local debris shedding with attached core',fontsize=15);savefig(fig,'failure_spatial_and_connectivity_diagnosis')

def validate():
    start=time.perf_counter();reg=read('SOURCE_AND_PRIOR_RESULT_REGISTER.json');bad=[]
    for r in read('PRIOR_STAT_GUARD.json'):
        p=Path(r['path'])
        if not p.exists() or p.stat().st_size!=r['bytes'] or p.stat().st_mtime_ns!=r['mtime_ns']:bad.append(str(p))
    checks=[]
    for r in reg['files']:
        p=Path(r['path']);checks.append(dict(path=str(p),sha256=digest(p),matches_registered_sha=digest(p)==r['sha256']))
    s=YSources();assert s.mesh_id==reg['native_mesh_sha256'] and s.flow_id==reg['native_flow_sha256'];assert not bad and all(r['matches_registered_sha'] for r in checks)
    m=read('RENDER_MANIFEST.json');rr=m['frames'];times=col(rr,'physical_time_s');assert len(rr)==484 and np.all(np.diff(times)>0)
    assert rr[233]['phase']=='LAMMPS_FRACTURE' and rr[234]['phase']=='ONE_WAY_BACKGROUND_FLOW';assert abs(times[233]-.00233)<1e-15 and abs(times[234]-.00234)<1e-15
    assert all(r['medical_surface_member_nodes']+r['medical_marker_member_nodes']==1927 for r in rr)
    assert all(r['fixed_cameras_verified'] for r in rr) and all(r['display_validation']['illegal_outside_nodes']==0 for r in rr[234:])
    bracket=rows('RENDER_TRANSPORT_SOURCE_BRACKETS.csv');assert len(bracket)==250*208
    for r in bracket:
        lo,hi=float(r['source_time_lo_s']),float(r['source_time_hi_s']);t=float(r['elapsed_s']);terminal=r['terminal_track']=='True'
        assert r['same_exact_members']=='True' and 0<=float(r['alpha'])<=1
        assert (terminal and lo==hi and t>=hi) or lo-1e-12<=t<=hi+1e-12
    video=RESULTS/'animations/y_clot_fracture_to_flow_seamless_v2.mp4';reader=imageio.get_reader(video);meta=reader.get_meta_data();count=reader.count_frames();assert count==484 and meta['fps']==25 and meta['size']==(1600,928)
    decoded={k:reader.get_data(k)[84:-44] for k in [0,116,233,234,300,483]};reader.close();pixelchecks=[]
    for a,b in [(0,116),(233,234),(234,300),(300,483)]:
        diff=float(np.mean(np.any(np.abs(decoded[a].astype(int)-decoded[b].astype(int))>8,axis=2)));assert diff>0;pixelchecks.append(dict(frame_a=a,frame_b=b,scene_pixel_fraction_difference_over8=diff))
    im=Image.open(RESULTS/'animations/y_clot_fracture_to_flow_seamless_v2.gif');assert im.n_frames==121
    tr=np.load(CASE/'accepted_trajectory.npz');meshchecks=[]
    for k in [0,233]:
        blocks=pv.read(RESULTS/'paraview'/f'raw_medical_step_{int(tr["step"][k]):04d}.vtm');raw=blocks['raw_all_nodes'];assert np.max(np.abs(raw.points-tr['position_m'][k]))<1e-17 and raw.n_points==3488
        represented=0
        for cid in np.unique(tr['component_id'][k]):
            ids=np.flatnonzero(tr['component_id'][k]==cid);represented+=len(ids)
            if len(ids)<3:continue
            prior=component_shell(tr['position_m'][k,ids]);new=pv.read(OUT/'shell_cache'/f'frame_{k:03d}_component_{cid}.vtp');same=np.array_equal(new.points,prior.points) and np.array_equal(new.faces,prior.faces);assert same and new.volume==prior.volume
            meshchecks.append(dict(frame=k,component_id=int(cid),members=len(ids),points_faces_exactly_equal_prior_reconstruction=same,volume_exactly_equal_prior_reconstruction=True,reference_d_eq_um=float((6*len(ids)*(50e-6)**3/np.pi)**(1/3)*1e6),display_extent_um=(np.ptp(new.points,axis=0)*1e6).tolist()))
        assert represented==1927
        if k==233:assert blocks['medical_unresolved_mass_markers'].n_points==220 and len(np.unique(blocks['medical_unresolved_mass_markers']['atom_id']))==220
    final=dict(status='PASS',history_stat_files_checked=len(read('PRIOR_STAT_GUARD.json')),history_stat_mismatches=bad,registered_hashes=checks,native_CFD_mesh_and_flow_identity_unchanged=True,new_movie_frame_count=count,new_movie_fps=meta['fps'],new_movie_size=meta['size'],new_gif_frames=im.n_frames,new_physical_clock_strictly_increasing=True,duplicate_physical_timestamps=0,scene_capture_unique_hashes=len(set(r['scene_pixels_sha256'] for r in rr)),decoded_scene_difference_checks=pixelchecks,camera_changes=0,render_interpolation_brackets_checked=len(bracket),display_interpolated_illegal_outside_nodes=0,all_1927_clot_nodes_represented_each_frame=True,final_medical_marker_nodes=220,final_free_unresolved_marker_nodes=219,final_anchored_singleton_marker_nodes=1,raw_includes_all_3488_native_nodes=True,mesh_geometry_checks=meshchecks,wall_physics_scope='PRIOR_SCOPE_ONLY; interpolation display center check does not extend prior physical wall validation',new_solver_runs=0,wall_seconds=time.perf_counter()-start)
    dump(OUT/'FINAL_INTEGRITY_CHECK.json',final)
    a=read('handoff_timeline_audit.json');a['audit']=read('HANDOFF_CONTINUITY_AUDIT.json');a['new']=dict(fps=25,frame_count=484,duration_s=19.36,absolute_physical_times_s=times,duplicate_timestamps=0,last_fracture_frame=233,last_fracture_physical_time_s=float(times[233]),first_transport_frame=234,first_transport_physical_time_s=float(times[234]),first_transport_elapsed_s=rr[234]['transport_elapsed_s'],first_transport_slowdown=1/rr[234]['physical_seconds_per_playback_second'],camera_changes=0,artificial_hold=False,playback_design=read('PLAYBACK_CLOCK_DESIGN.json'),scene_pixel_capture_verified=True,source_state_interpolation='Adjacent actual accepted transport states, fixed members only; no LAMMPS interpolation');dump(OUT/'handoff_timeline_audit.json',a)
    return final

def reports():
    a=read('HANDOFF_CONTINUITY_AUDIT.json');diag=read('FRAGMENTATION_MECHANISM_DIAGNOSIS.json');perf=read('PERFORMANCE_RENDER_ONLY.json');rejected=read('development_rejected_static_capture/PERFORMANCE_RENDER_ONLY.json');perf['rejected_development_render_wall_seconds']=rejected['total_render_and_static_wall_seconds'];perf['rejected_attempt_reason']='Off-screen stale capture buffer; preserved separately, excluded from deliverable; explicit render before capture fixed and decoded scene differences checked.';perf['handoff_readonly_audit_wall_seconds']=a['wall_seconds'];perf['mechanism_readonly_audit_wall_seconds']=diag['wall_seconds'];perf['accounting_scope']='Individual recorded audit/render stages; preview, plotting, packaging and interactive development are not a full task wall-clock benchmark.';dump(OUT/'PERFORMANCE_RENDER_ONLY.json',perf)
    status=dict(HANDOFF_POSITION_CONTINUITY='PASS',HANDOFF_VELOCITY_CONTINUITY='PASS',ARTIFICIAL_VIDEO_HOLD='FOUND_AND_FIXED',BACKGROUND_DRAG_DURING_FRACTURE='NOT_INCLUDED',UNRESOLVED_DEBRIS_RENDER_HONEST='PASS',CLOT_SURFACE_VISUAL_QA='PASS',WALL_NO_PENETRATION='PRIOR_SCOPE_ONLY',MACRO_FRAGMENT_200UM='NOT_OBSERVED',MACRO_FRAGMENT_300UM='NOT_OBSERVED',LARGE_BLOCK_FRACTURE_MECHANISM='NOT_ESTABLISHED',CLOT_PRESENT_CFD='NOT_RUN',ONE_WAY_TRANSPORT='EXPLORATORY_ONLY',BIOLOGICAL_FRAGMENT_SIZE_VALIDATION='NOT_ESTABLISHED',CURRENT_FAILURE_MODE='LOCAL_DEBRIS_SHEDDING_WITH_ATTACHED_CORE',BACKGROUND_FLOW_NOT_APPLIED_DURING_FRACTURE=True,TRUE_VELOCITY_RESET_BUG='NOT_OBSERVED',GATE_A='COMPLETED',GATE_B='COMPLETED',GATE_C='COMPLETED',GATE_D='NOT_TRIGGERED_NO_CONFIRMED_BUG',GATE_E='COMPLETED',new_LAMMPS_solves=0,new_transport_solves=0,new_CFD_solves=0,completion='STOPPED_FOR_HUMAN_REVIEW',limitations=['Attached body remains native final snapshot despite nonzero native velocity; no continuation of attached-body mechanics.','Three-dimensional walls and coordinates are original; display interpolation center checks do not expand earlier physical wall validation.','Medical visual QA is relative to prior renderer; no biological validation or real subsurface-scattering physics.'])
    dump(OUT/'FINAL_STATUS.json',status)
    (OUT/'DISPLAY_VS_PHYSICS_CONTRACT.md').write_text('''# 显示与物理契约

科学来源为原 C1、原 native surviving/broken bonds、原 handoff、原 accepted transport steps，均只读。本轮没有新的 LAMMPS、输运或 CFD 求解。

| 项目 | 本轮显示处理 | 保留的物理含义 / 限制 |
|---|---|---|
| LAMMPS 阶段 | 全部234个真实保存帧，0–2.330 ms；不插值 | 每10 µs显示一次；0.5 µs原生断键时间保留在事件表，帧间事件不伪造中间姿态 |
| 输运阶段 | 同一组件固定成员、真实相邻 accepted steps 间线性插值 COM，仅用于画面 | 每组件插值端点、alpha逐帧记录；并非新积分状态，差分显示速度不代替原始速度 |
| 播放时钟 | 25 fps；早期接续10 µs增量，再按已保存轨迹的质量加权路径长度逐渐调整 | 每帧标注真实绝对时间及慢放/压缩倍数；只改变选帧时刻，不改物理轨迹 |
| 附着主体及≥3节点组件 | 复用上轮各组件独立外壳；仅平滑法向和深红材质，柔和光照 | 顶点、面、外壳体积逐项相等；不跨组件bridge、不改裂口、不做几何膨胀；外壳体积仍是显示代理，不等于连续血栓体积 |
| 1–2节点组件 | 每真实节点一个2 px、opacity0.32的细小质量标记，无实心球和拖尾 | 保留全部质量、身份与坐标；像素尺寸不是生物碎片粒径；无浓度云或扩散模型 |
| 末态质量覆盖 | 1701节点主体+6节点游离壳+220标记=1927节点 | 220标记=219游离未解析节点+1附着单节点；225游离节点仍属208个组件 |
| 附着主体 Phase2 | 保持2.330 ms真实末态 | 1702附着节点的末态COM速度约0.054716 m/s，后续弹性运动未求解；静止画面不代表物理已停止 |
| 管壁和镜头 | 原生WALL；固定局部、全Y、分叉视角，坐标统一1x；局部剖切仅显示 | full/raw/ParaView保留完整WALL，不移动管壁，不修正物理碰撞 |
| 箭头 | Phase2才显示原冻结未含血栓CFD的速度方向 | 只指示该阶段模型启用，不暗示Phase1含流动或流固耦合 |
| 出口后的标记 | 保留真实终止轨迹位置 | 已退出的轨迹停止记录，不能解读为出口沉积；片尾字幕明确说明 |

材质使用局部法向平均、缓慢变化的RGB和柔和光照模拟视觉柔和度；没有声称实现经过验证的组织散射、真实粗糙度或医学影像。原50 µm支撑仍限制轮廓，不能靠shader消除离散模型的几何分辨率。

同时提供两时刻的raw/medical PNG、PDF和可切换ParaView VTM。raw包含3488原始节点（1927 clot+1561 anchors）、存活typed bonds、组件身份；医学组件各自独立，未解析标记独立block，SI坐标。具体开关见results/paraview/README.md。

输运插值检查所有225游离节点显示位置处于native lumen或真实端口外延。它仅约束渲染插值合法性；不能替代上一轮的积分内墙约束与扫掠审查，也不消除近壁流速、主体穿行风险、无流体反馈、无转动/变形等原限制。
''')
    (OUT/'no_artificial_pause_report.md').write_text('''# 交接处停顿审查与修正

真实数据没有速度重置：208游离组件的handoff位置跳变与速度跳变均为0；native末态、handoff JSON、transport t=0逐组件一致。相同最后成员集合在此前尚未独立时只记作cohort，不能倒推为已游离。

旧244帧、20 fps视频在零基帧117/118重复同一个2.330 ms物理端点，多出0.05 s，端点总显示0.10 s。没有证据支持更长的额外多帧hold。标签从2.330 ms切到0.00 s相对输运时间，属于显示时钟切换。通常播放比率同时从2500倍慢放变成1.25倍慢放，跨阶段相差2000倍；首个非零输运显示时间已到40 ms，100 µs松弛几乎完全被跳过。

新版484帧、25 fps、19.36 s只保留一次handoff端点；帧233=2.330 ms，帧234=2.340 ms（transport_elapsed=10 µs），首个跨阶段步仍为4000倍慢放。随后按已保存质量加权轨迹路径长度逐渐调整取帧时间，逐帧显示倍数，绝对物理时钟严格递增，固定镜头，无重复端点或人为hold。完整时钟与逐组件插值出处见handoff_timeline_audit.json、RENDER_MANIFEST.json及RENDER_TRANSPORT_SOURCE_BRACKETS.csv。

保留真实模型现象：背景drag仅在handoff后加入，1–2节点采用100 µs工程松弛，早期速度迅速下降并可能转向；局部CFD较慢且100 µs输运可低于像素分辨率。附着主体的后续力学没有求解，Phase2只显示其末态。已退出组件显示最后记录位置。上述情况不靠假轨迹掩盖，字幕明确交代。

开发中的首次新版成片因缓存VTK画面未刷新而被拒收，已保留于development_rejected_static_capture。修正为每帧capture前显式render，并比较不含字幕的原始及解码画面。最终交付只使用通过动态检查的版本；该开发问题与历史handoff物理数据无关。
''')
    (OUT/'FUTURE_BACKGROUND_DRAG_COUPLING_CONTRACT.md').write_text('''# Gate D：本轮不触发求解；未来契约尚未实现

位置/速度连续且无单位或初始化错误，因此没有输运重算。本轮也没有经验证的“组件脱离即施加CFD drag”实现，不执行可选LAMMPS对照。

未来若实现，需先完成以下共同契约：在每个原生断键时刻基于MPI全局一致存活图判定真实脱离、维护稳定成员/质量与拆分继承关系；只对已脱离组件启用已标记的单向流力，时间戳和COM速度连续。原500/s数值阻尼与CFD drag的关系必须明确，不能双计后悄然改称物理流阻。1–2节点的100 µs松弛仍标为工程代理，不由dx宣称粒径验证；≥3节点球形阻力修正也须保留形状近似说明。使用同一native P1 CFD、SI单位、原真实WALL和端口；单独验证分布到节点的总力/力矩、积分稳定性、墙约束及耗散，不在renderer中补偿力学。附着主体的流体加载和后续动力学范围必须独立说明。

此契约不是现有功能，也不是全程FSI。没有体内实测或微泡场可为新增力背书。本轮物理参数、载荷、几何和背景流身份全部冻结。
''')
    rate=read('RENDER_MANIFEST.json')['frames'];range500=a['sampling_after_500us_elapsed_s']
    (OUT/'REVIEW_Y_CLOT_CONTINUOUS_HANDOFF_AND_HONEST_RENDER_V1.md').write_text(f'''# Y_CLOT_CONTINUOUS_HANDOFF_AND_HONEST_RENDER_V1 审查

本轮已完成原生交接审查、原生断裂机制诊断、离线重制MP4/GIF及raw/medical对照。没有新增LAMMPS、输运或CFD求解；历史244帧视频、上轮数据与代码保持只读。结论仍是 **LOCAL_DEBRIS_SHEDDING_WITH_ATTACHED_CORE**，未实现“大块血栓破裂”。

## 直接查看

- [新版MP4](../../results/y_clot_continuous_handoff_and_honest_render_v1/animations/y_clot_fracture_to_flow_seamless_v2.mp4) / [GIF](../../results/y_clot_continuous_handoff_and_honest_render_v1/animations/y_clot_fracture_to_flow_seamless_v2.gif)
- [初态raw/medical](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/raw_vs_honest_medical_initial.png) / [交接态raw/medical](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/raw_vs_honest_medical_handoff.png)，均另有PDF。
- [速度与播放时钟诊断](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/handoff_speed_and_playback_audit.png) / [破坏机制诊断](../../results/y_clot_continuous_handoff_and_honest_render_v1/figures/failure_spatial_and_connectivity_diagnosis.png)
- [显示契约](DISPLAY_VS_PHYSICS_CONTRACT.md)、[最终状态](FINAL_STATUS.json)、[逐组件交接表](handoff_velocity_audit.csv)、[完整帧清单](RENDER_MANIFEST.json)。

## Gate A：没有真实速度重置

直接重读native LAMMPS最后5个保存帧（2.290–2.330 ms），逐组件核对handoff和transport t=0，并审查2471个早期真实accepted steps（每组件至少5步，含x、v、P1 CFD速度、阻力、墙约束标记）。208组件的位置最大跳变 **0 m**，速度最大跳变 **0 m/s**，包括6节点分量和207个1–2节点分量。其质量总计2.8125e−8 kg，成员225节点。

| 类别 | 数量 | handoff速度 min / median / max (m/s) | 所处CFD速度 min / median / max (m/s) |
|---|---:|---|---|
| 1–2节点质量包 | 207 | 0.095294 / 0.858503 / 1.568498 | 0.000320 / 0.004585 / 0.009877 |
| 6节点分量 | 1 | 0.599820 / 0.599820 / 0.599820 | 0.004657 / 0.004657 / 0.004657 |

原阻力在当前实际pre-step状态重建：F_start=m(u−v)/tau，平均力来自该冻结局部流速下指数松弛增量；它不是新增求解。未发生墙约束的早期步速度复现误差上界2.255e−16 m/s，3个带墙约束的步另标记，不能把约束冲量误算为drag。

在“各组件首个不早于500 µs的真实接受步”上，1–2节点速度/初速度比中位数0.006822，方向变化中位数18.727°；6节点比值0.112257、转向0.151°。**这些并非同一500 µs时刻**：采样时刻范围{range500['min']*1000:.3f}–{range500['max']*1000:.3f} ms，中位数{range500['median']*1000:.3f} ms；逐行精确时刻在CSV。它们支持快速松弛，不能包装成更高时间分辨率的实际轨迹。

原因分开记录：TRUE_VELOCITY_RESET_BUG=NOT_OBSERVED；背景drag在handoff开启的力模型切换、100 µs工程代理松弛、旧视频2000倍播放率切换均存在。19组件处局部CFD低于1 mm/s，全部208组件单按当地CFD估算的100 µs位移低于旧局部视图1像素；这只是背景平流的像素界限，并不否认handoff初始惯性位移。

另外，1702附着节点的末态COM速度仍为0.054716 m/s，而Phase2未继续求解其弹性运动，只保留末态。主体“停住”有这个明确的模型截断原因，不能称为输运速度清零bug，也不能声称其已物理静止。

## Gate B/C：视频与医学显示的修正

旧视频只查到一处重复端点：零基帧117/118多出50 ms显示，总端点显示100 ms，没有长时间人为hold证据。新版删除重复端点，234个原生保存帧直接接250个输运显示时刻，25 fps共19.36 s；帧233=2.330 ms、帧234=2.340 ms，首个跨阶段步10 µs，与前段4000倍慢放一致。之后仅改变选帧时刻，画面持续给出绝对物理时间与每帧慢放/压缩倍数，不重置为0。不能把自然减速全部“消除”。

输运画面只对相邻真实accepted steps之间同一成员组件的COM线性插值；52000条组件/帧出处可追溯，全部225游离成员显示中心均经native域/端口合法性检查。没有LAMMPS位置或断键拓扑插值。所有镜头固定，流场箭头仅在Phase2出现，字幕明确未含血栓的单向背景流、附着主体末态快照及出口终止标记。

医学显示复用原各组件独立外壳，仅改深红材质、表面法向与柔和光照；关键网格逐顶点/逐面/体积与原重建完全相同。细小组件不拼接、不膨胀。橙色实心球改为真实节点上的2 px半透明质量标记；末态220标记包含219游离未解析节点和1附着单节点。加上1701节点主体和6节点壳恰好覆盖1927节点，质量未删。6节点参考体积等效直径约112.725 µm；该变形组件显示包围盒为266.872×141.507×192.018 µm，与原外壳逐点相同。最长方向大于200 µm不等于体积等效直径达到200 µm，也不是渲染新增膨胀，不能冒充几百微米体积的大块或经生物验证的组织碎片。

首次开发渲染的画面缓存未刷新，已拒收保留；修复后检查不含时间字幕的实际画面变化及编码后解码帧，避免“只有时间在走”。材料视觉QA只确认比上轮柔和且忠于数据；50 µm离散支撑仍限制轮廓，并非真实组织光学验证。

## Gate D：未触发额外求解

没有速度/单位/初始化bug，故不重算transport。也没有经过验证的全程detachment-triggered drag实现，未运行可选LAMMPS+transport。未来最小契约见FUTURE_BACKGROUND_DRAG_COUPLING_CONTRACT.md，包含原生脱离图、数值阻尼与流阻关系、质量/力/单位、墙约束和附着主体范围。本轮没有声称FSI。

## Gate E：为什么还没有成块破裂

11991/82270条bulk键、819/22040条interface键断裂。按参考键中点到真实初始露出面的距离，**99.2828% bulk断键位于100 µm内**；0–50 µm层断键率53.016%，50–100 µm层19.552%，100–150 µm层0.449%。因此浅表剥落不只是渲染错觉。该深度是初始表面距离代理，不是连续裂面厚度。

960节点初有interface连接，16节点失去全部直接interface键，未建立宏观脱粘前缘。总附着1702节点需精确区分为 **1701节点最大核心+1附着单节点**；上轮汇总的88.323819%是总附着质量。最大核心仍有70261条bulk键、21220条interface键并连接1282个anchor。三轴中央80%范围的45个截面至少仍有184条核心bulk键穿过，说明保留通路，不能据此声称任意方向不存在局部裂纹。

重建全部2030个原生断键时刻，得到291个精确成员集合生命周期；峰值游离块13节点、参考体积等效直径145.865 µm。末态208游离分量仍只有一个6节点，其余207个1–2节点。没有≥200 µm（至少34节点）或≥300 µm（至少114节点）的游离大块。7个极短生命周期没有匹配保存帧COM，明确置空；native触发断键位置、参考COM、成员、质量和生灭时刻保留，不伪造首次出现的COM。

目前Gaussian参考200 Pa、Liu-LPS critical stretch、曲面锚定组合没有形成隔离宏观体积的切割/剥离路径。连通图说明残存路径，不证明唯一因果；提高力幅也可能加重节点级粉化，不能保证大块。下一阶段只列一个独立候选：同一真实Y、固定材料上的“宽尺度梯度切向剥离对照”，先定义载荷身份及连贯前缘指标，本轮不运行。

## 冻结范围、开销和交付

K=33333.333333333336 Pa、G=940 Pa、s00_bulk=0.15、s00_interface=0.140624859375、contact_distance_factor=1.0dx、原三层anchor及native WALL均未改。V3连续体积0.271574388 mm³，离散1927节点/0.240875 mm³保持。B1/B2 TERMINATED，PD2 fracture STOPPED，FSI与bubble streaming未接入；含血栓CFD NOT_RUN。

原有近壁不确定、保留主体附近路径不确定、无流体反馈、无片间碰撞和无转动/变形等限制仍在。WALL_NO_PENETRATION=PRIOR_SCOPE_ONLY；本轮显示插值合法性检查不扩张上一轮碰撞验证范围。生物碎片尺寸仍NOT_ESTABLISHED，ONE_WAY_TRANSPORT=EXPLORATORY_ONLY。

记录的只读交接审查约{a['wall_seconds']:.1f} s，机制审查{diag['wall_seconds']:.1f} s；通过版本渲染/静态图{perf['total_render_and_static_wall_seconds']:.1f} s，峰值renderer RSS {perf['peak_renderer_RSS_bytes']/1024**3:.2f} GiB（不含encoder子进程）。另保留被拒首次渲染{rejected['total_render_and_static_wall_seconds']:.1f} s，不将开发成本隐去。均不是性能扫描；详见PERFORMANCE_RENDER_ONLY.json。

568个历史文件stat保护及关键数据SHA核对结果见FINAL_INTEGRITY_CHECK.json；源码、CSV、PNG/PDF、ParaView、报告与视频分别保留。evidence.zip不复制大型原轨迹、全部中间shell缓存或被拒视频；原路径与hash完整记录，MP4/GIF作为独立交付并列hash。本轮结束，等待人工审查。
''')
    (RESULTS/'paraview/README.md').write_text('''# 原生/医学同状态对照

用ParaView打开raw_medical_step_0000.vtm或raw_medical_step_4660.vtm，坐标单位m。两个文件均为真实LAMMPS保存时刻，没有输运插值。

raw_all_nodes包含全部3488节点的atom_id、component_id和attached；raw_surviving_typed_bonds是当时真实存活键。medical_component_*各自为一个真实连通组件；medical_unresolved_mass_markers保留每个1–2节点组件的真实节点、ID、质量及附着标记。末态220标记=219游离节点+1附着节点。

在Multiblock Inspector中分别显示raw节点/键或medical壳/标记。标记选Points、Size2、Opacity0.32，不使用球体Glyph。壳按medical_rgb着色；native_WALL_uncut按需透明。显示外壳为上轮重建，当前仅法向/材质改动；不将其体积当机械体积。单独隐藏一个block只是审阅开关，正式视频始终包含全部质量表示。
''')

def package():
    src=ROOT/'pd_clot/lammps_fragmentation';names=['continuous_handoff_audit.py','clot_failure_mechanism_audit.py','honest_clot_renderer.py','continuous_handoff_delivery.py'];snapshot=OUT/'SOURCE_DELIVERY';snapshot.mkdir(exist_ok=True)
    for name in names:(snapshot/name).write_bytes((src/name).read_bytes())
    register=[]
    for p in [*snapshot.glob('*.py'),*RESULTS.rglob('*')]:
        if p.is_file():register.append(dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)))
    dump(OUT/'DELIVERY_ARTIFACT_REGISTER.json',dict(files=register))
    git=[]
    for cmd in [['git','branch','--show-current'],['git','rev-parse','HEAD'],['git','status','--short']]:git+=['$ '+' '.join(cmd),subprocess.check_output(cmd,cwd=ROOT,text=True)]
    (OUT/'ENDING_GIT_STATE.md').write_text('\n'.join(git))
    archive=OUT/'evidence.zip';files=[]
    for p in OUT.rglob('*'):
        if p.is_file() and not set(p.relative_to(OUT).parts)&{'shell_cache','development_rejected_static_capture'} and p.name not in ['evidence.zip','EVIDENCE_PACKAGE.json','DELIVERY_RUN.log']:files.append(p)
    files.extend(p for p in RESULTS.rglob('*') if p.is_file() and p.suffix not in ['.mp4','.gif'])
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(files)):z.write(p,p.relative_to(ROOT))
        assert z.testzip() is None
    dump(OUT/'EVIDENCE_PACKAGE.json',dict(path=str(archive),bytes=archive.stat().st_size,sha256=digest(archive),files=len(set(files)),excluded='Large original trajectories, intermediate shell cache, rejected movie, final MP4/GIF (delivered separately with hash).'))
    print('DELIVERED',str(archive),archive.stat().st_size,flush=True)

def main():
    complete_handoff_table();plots();validate();reports();package()

if __name__=='__main__':main()
