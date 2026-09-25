#!/usr/bin/env python3
"""Source-bound Chinese review, validation record, local gallery and inventory."""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import argparse,csv,hashlib,html,subprocess,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle8_replay import REPO,read,digest
from particle_3d.particle81_simulation import OUTPUT,dump,lock_upstream
from particle_3d.particle81_replay import Scene,TITLE,INDEPENDENT,FOOTERS,OUTLETS
from particle_3d.particle81_visuals import ANIMATIONS
from particle_3d.particle81_figures import FIGURES
from verify_particle81 import source_files,artifact_files


def require(condition,message):
    if not condition:raise RuntimeError(message)


def finalize(commit):
    root=OUTPUT;scene=Scene(root);cat=scene.catalog;births=read(root/'data/birth_ledger.json')
    tests=read(root/'data/test_run.json');physics=read(root/'data/physics_audit.json');audit=read(root/'data/trajectory_audit.json')
    media=read(root/'data/media_audit.json');visual=read(root/'data/visual_inspection.json');repeat=read(root/'data/render_determinism.json')
    require(tests['all_pass'] and physics['all_pass'] and audit['all_pass'],'Automated checks failed')
    require(cat['completed']>=300,'Fewer than 300 real completed inlet-to-outlet trajectories')
    require(tests['source_sha256']==source_files(),'Source changed since permanent tests')
    require(tests['artifact_sha256']==artifact_files(),'Evidence changed since tests')
    for suite in tests['suites']:
        for key in ['junit','log']:require(digest(root/suite[key])==suite[key+'_sha256'],'Test receipt changed')
    for path,sha in source_files().items():
        saved=subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=REPO)
        require(hashlib.sha256(saved).hexdigest()==sha,'Source commit mismatch: '+path)
    require(media['all_pass'] and len(media['animations'])==4,'Incomplete media decoding')
    require(visual['result']=='PASS_WITH_LIMITATIONS','Missing actual visual inspection')
    require(set(visual['figures'])==set(FIGURES) and set(visual['animations'])==set(ANIMATIONS),'Incomplete figure/video inspection')
    for path,sha in visual['reviewed_file_sha256'].items():require(digest(root/path)==sha,'Visual review stale: '+path)
    sensitivity=read(root/'data/timestep_sensitivity.json')
    quarter={r['particle_id']:r for r in sensitivity['rows'] if r['dt_factor']==4}
    preserved=lock_upstream();completed=[scene.entries[i] for i in scene.completed]
    samples=np.concatenate([scene.arrays[i] for i in scene.completed]);length=np.array([e['path_length_m'] for e in completed]);res=np.array([e['residence_time_s'] for e in completed])
    failures=Counter(str(e['failure_detail']) for e in scene.entries.values() if e['failure_detail'])
    limits=[*FOOTERS,
        'Completed means first center crossing of an official open outlet cap; not full finite-body clearance.',
        'Numerical safety stops and finite-size admission guards censor the cohort; outlet fractions are not physiological probabilities.',
        'Unobserved branches are not reconstructed. Displaying the full mesh does not establish full network coverage.',
        'Long constant nominal concentration is not a post-bolus pharmacokinetic prediction.',
        'No acoustic image formation, localization noise, detection probability or measured ULM reconstruction validation.',
        'Handoff event-time convergence and production timestep are not established.']
    total_tests=sum(s['tests'] for s in tests['suites']);covered=sum(bool(n) for n in cat['outlet_counts'].values())
    provenance=read(root/'data/frozen_provenance.json')['provenance']
    outlet_flow={o:provenance['boundary_manifest']['boundaries'][o]['signed_flux_m3_s'] for o in OUTLETS}
    status=dict(AUTOMATED_CHECKS='PASS',MINIMUM_COMPLETED_SCALE='PASS',PREFERRED_1000_SCALE='PASS' if cat['completed']>=1000 else 'NOT_REACHED',
        FULL_NETWORK_BRANCH_COVERAGE='ALL_OUTLETS_OBSERVED' if covered==3 else 'INCOMPLETE',
        SAVED_SCENE_REPLAY='PASS',VISUALIZATION='PASS_WITH_LIMITATIONS',PHYSIOLOGICAL_COMPLETION='NOT_ESTABLISHED')
    record=dict(stage='Particle-8.1 Full-vessel ULM microbubble trajectory simulation',title=TITLE,created_utc=datetime.now(timezone.utc).isoformat(),
        statuses=status,source_commit=commit,source_commit_role='EXACT_SOURCE_MATCH_BOUND_TO_PERMANENT_TESTS',source_sha256=source_files(),
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=REPO,text=True).strip(),user_acceptance='PENDING_USER_REVIEW',
        independent_trajectory_superposition=True,independence_statement=INDEPENDENT,scientific_limitations=limits,
        scheduled=cat['scheduled'],admitted=cat['admitted'],completed=cat['completed'],outlet_counts=cat['outlet_counts'],end_reasons=cat['end_reasons'],
        failure_detail_counts=dict(failures),unclassified_successful_exits=0,outlet_frozen_volume_flux_m3_s=outlet_flow,
        outlet_frozen_volume_flux_fraction={o:q/births['Q_in_m3_s'] for o,q in outlet_flow.items()},
        acquisition_birth_window_s=cat['acquisition_birth_window_s'],acquisition_replay_end_s=cat['replay_end_time_s'],
        physical_sample_count=cat['physical_samples'],completed_path_length_m=dict(min=float(length.min()),median=float(np.median(length)),max=float(length.max())),
        completed_residence_time_s=dict(min=float(res.min()),median=float(np.median(res)),max=float(res.max())),
        completed_saved_speed_m_s=dict(min=float(np.linalg.norm(samples[:,4:7],axis=1).min()),max=float(np.linalg.norm(samples[:,4:7],axis=1).max())),
        clock=dict(Q_in_m3_s=births['Q_in_m3_s'],C_MB_m3=births['C_MB_m3'],H_D_feed=births['H_D_feed'],
            interarrival_s=1/(births['Q_in_m3_s']*births['C_MB_m3']),mode=births['scheduler'],orientation='ISOTROPIC_V0_MODEL_ASSUMPTION'),
        integration_config=completed[0]['integration_config'],full_geometry_wall_triangles=45221,source_scene_sha256=scene.sha256,
        automatic_tests=tests['suites'],total_tests_run_this_stage=total_tests,prior_all_stage_tests_rerun=False,
        trajectory_audit='data/trajectory_audit.json',physics_audit='data/physics_audit.json',physics_audit_sha256=digest(root/'data/physics_audit.json'),
        timestep_sensitivity=sensitivity,software_environment=read(root/'data/software_environment.json'),repeat_render=repeat,media_audit=media,visual_inspection=visual,display_contract=cat['display_only'],
        preserved_upstream_files=preserved,upstream_baseline_commit=read(root/'data/upstream_lock.json')['baseline_commit'],
        evidence_sha256=tests['artifact_sha256'])
    dump(root/'PARTICLE8_1_VALIDATION.json',record)
    outlet_text='、'.join(f'{o}：{cat["outlet_counts"][o]}' for o in OUTLETS)
    lines=['# Particle-8.1：完整血管微泡轨迹可视化层——中文审核报告','',TITLE,'',
        f'已实际生成 **{cat["completed"]} 条真实 INLET→OUTLET 完成轨迹**，保存全量接受子步、原始属性和失败记录，输出十张静态图、四组 MP4/GIF、逐帧账本与 CSV/JSON。达到 frozen-flow full-vessel ULM visualization 的可展示标准；不是 RBC-coupled physiological completion。',
        f'**覆盖范围：{covered}/3 个出口分支有完成轨迹。完整 Frozen 网格均已显示，未观测分支不算轨迹重建完成。**','',
        '## 九项审核问题','',
        '1. **是否真正入口到出口？** 是。每条完成轨迹均用原 P6.5 积分器在官方 Frozen FEM 中推进，并独立复核末接受线段与真实 OUTLET 三角面相交。出口语义是球心首次穿过开放端盖，不宣称整个有限球体完全离开。无人工直线补齐、短轨迹延长或 synthetic 通道。',
        f'2. **完成数量？** scheduled={cat["scheduled"]}，admitted={cat["admitted"]}，completed={cat["completed"]}。最低 300 条要求已达到；优选 1000 条状态：{status["PREFERRED_1000_SCALE"]}。所有未完成 ID 保留原始尺寸和明确终止原因。',
        f'3. **采集时间？** 从 t=0 到末次计划出生为 {cat["acquisition_birth_window_s"]:.9f} s（{cat["acquisition_birth_window_s"]/3600:.6f} h）；最终记账到 {cat["replay_end_time_s"]:.9f} s。相邻出生间隔 {record["clock"]["interarrival_s"]:.9f} s。长时保持名义浓度是此数据集的假设，不是一次注射后的有效 PK 持续时间预测。',
        f'4. **出口数量？** {outlet_text}。未分类的成功退出=0。终止/入场条件对样本有筛选，不能把此出口比例当作未截尾人群的生理分流概率。',
        '5. **独立叠加？** 是。'+INDEPENDENT+' 每个 MB 使用自身稳定 ID/出生时间/原始尺寸/朝向/入口随机流独立积分；不含 MB–MB、RBC–MB 同时水动力耦合。',
        '6. **旋转与累积？** 已输出四组完整血管旋转、ULM 式累积、三条代表性完整 journey、出口着色动画；均由同一 scene 导出。完整网格包含全部 45,221 个 WALL 三角形及四个开放端盖。累积只纳入当时已经完成的轨迹；未发生的未来轨迹不提前出现。',
        '7. **补上了什么？** 此前真实 P8 只有入口附近 1 ms 单 MB smoke replay。本阶段新增长采集出生账本、全血管独立轨迹积分、每 ID 终止审计、数百/上千条完整 passage 累积与完整单 MB journey。没有把此前 synthetic lifecycle 作为真实通行证据。',
        '8. **还缺什么？** 真实 RBC passage、RBC-resolved suspension、完整 PK、新 CFD、Particle-7.5、新 RBC deformation 均未建立。P6.5 仍限 sphere-normal regularization，非球形 lubrication、多体 mobility、production dt/neighbors 与 handoff event-time 收敛未冻结。没有声学成像或真实定位噪声模型；未覆盖出口分支仍缺轨迹。',
        '9. **是否可展示？** 是，作为真实 Frozen 三维几何中可追溯的独立微泡完整路径可视化层；应连同覆盖缺口、计算截尾与上述模型边界展示。完整生理物理完成、整棵血管所有分支重建不能由此自动推出。','',
        '## 数据与数值证据','',
        f'- 共 {cat["physical_samples"]} 个原始保存状态；完成轨迹长度 min/median/max={length.min()*1e6:.3f}/{np.median(length)*1e6:.3f}/{length.max()*1e6:.3f} µm；residence={res.min()*1e3:.3f}/{np.median(res)*1e3:.3f}/{res.max()*1e3:.3f} ms。',
        f'- Q_IN={births["Q_in_m3_s"]:.12g} m³/s，C_MB=8.5e12 m⁻³，H_D=0.45。P7 deterministic cumulative flux clock 原样复用；未改成 Poisson。H_D 是 feed volume fraction，本任务不模拟 RBC，也不强设 tube Hct。',
        '- 原 SonoVue inverse CDF 对每 ID 只抽尺寸一次；P7 isotropic orientation V0 明确是模型假设。位置采用原 P7 局部法向流量权重 proposal；FiniteSizeAdmission 最多 512 次重试位置，不重抽尺寸，不增加入口偏移。通过者分布已被有限尺寸条件化；Figure 02 分开显示初始 proposals 与 admitted 位置。',
        '- 每 ID NPZ 保留所有接受子步，含 elapsed t、位置、上一接受区间使用的速度/角速度、q、h_geom、g_nf、h_lower 和 nearfield state。首行速度是入场自由流诊断；不能当作已执行位移。CSV.gz 另含物理绝对时间、速度范数和 stable ID。绝对时钟用于 acquisition，elapsed t 保留长采集下的子步精度。',
        '- 原 P6.5 接受/拒绝、连续 wall 证书与二分时间细化保持原样；只用精确输入查询缓存减少重复计算。批量开发中扩展了相同 BVH 请求缓存，较早记录使用仅 wall-gap 缓存；两版均不改变物理算术，完整接受数组、终止原因与 provider 次数的对照见 cache_v2_parity.json。永久测试比较 native 与缓存适配器，原 full native ID 1/4 的完整原始数组也作为逐位一致对照；安全停止 ID 13 的全部保留状态与原生 1.5 s 长计算的前缀逐位一致，证明该保护只截断而未改变已接受轨迹。可选全局运动界优化已测但正式队列未开启。',
        '- dt=0.25 ms 沿用 P7 验证配置，只是本数据集数值步长，不是 production freeze；每 MB elapsed horizon=1.5 s。确定性计算保护包括 16,000 provider calls、512 次无接受推进、32 个名义步完全不动、64 个接受子步仅有 roundoff 量级位移或总时间推进低于名义 dt 的 1%。保护只停止计算并保存此前状态，不改物理力、位置、尺寸或接受规则。',
        '- 近壁计算停止不能解释为真实捕获或排除通行：原规则有 2 nm handoff 附近固定状态和纳秒级细分。诊断归档保留原完整 1.5 s 停滞例与开发参考；这些归档不混入正式队列。终止比例依赖计算保护，尚未建立 handoff 事件时间/轨迹完成率收敛。','',
        '| 终止类型 | 数量 |','|---|---:|']
    lines += [f'| {k} | {v} |' for k,v in cat['end_reasons'].items()]
    lines += ['', 'Frozen 的三个出口均有非零体积流量；某出口没有完成 MB 轨迹，不表示该分支没有流动。有限尺寸入场条件化和计算截尾使该队列不能直接替代无偏体积通量示踪。', '',
        '| 出口 | Frozen 流量 (m³/s) | Q_OUT / Q_IN | 完成 MB |', '|---|---:|---:|---:|']
    lines += [f'| {o} | {outlet_flow[o]:.9g} | {outlet_flow[o]/births["Q_in_m3_s"]:.4%} | {cat["outlet_counts"][o]} |' for o in OUTLETS]
    lines += ['', '另对预先选定 ID 1、4（基准完成）及 7、13（基准计算停止）分别重算 dt/2 与 dt/4，结果见 [步长敏感性](data/timestep_sensitivity.json)。这些重算不计入正式队列。比较只在共同 elapsed 时间内线性插值计算误差，不外推。四例有限对照不能建立所有近壁事件时间或完成率收敛。',
        f'ID 1/4 在 dt/4 的驻留时间分别为 {quarter[1]["refined_elapsed_s"]*1000:.6f}/{quarter[4]["refined_elapsed_s"]*1000:.6f} ms；相对基准变化 {quarter[1]["residence_difference_s"]*1000:.6f}/{quarter[4]["residence_difference_s"]*1000:.6f} ms；共同时间轨迹最大偏差为 {quarter[1]["common_time_position_difference_max_m"]*1e6:.6f}/{quarter[4]["common_time_position_difference_max_m"]*1e6:.6f} µm。两例仍从 OUTLET_02 退出，但不可把约亚微米偏差表述成零误差。ID 7/13 仍计算停止；其中 32 个名义步的静止保护随 dt 改变等待时长，因此停止时刻不是物理捕获时刻。',
        '', '## 显示变换与回放','',
        '- 所有可视化只读同一 saved scene：线性位置插值、P8 SLERP 朝向插值只作用于显示，不外推、不回写物理状态。轨迹使用保存折线，无样条、无路径抽稀、无轴向拉伸。仅平移显示原点并把米换算为 µm，当前球保留原半径。MB ID 采用屏幕空间引线标注以避开端盖文字；引线是注释，不计入物理轨迹。',
        '- 相机方位 30°→120°、俯仰 28°，正交等比例显示。1/4 动画慢放代表性真实出生到退出窗口并明确 CUT 空闲时间，最后切至最终记账；2 将完整采集窗口压缩为 240 帧；3 慢放并顺序显示三条独立的完整 birth→exit，并标记跨 MB 空闲时段切换；每次到达出口后保持已计算端点状态 12 个显示帧（0.8 视频秒），用于观察退出，不推进物理时间。单 MB 动画侧栏明确只统计当前选定 ID，其余三类侧栏统计完整 acquisition cohort。真实出生稀疏，未通过复制活跃 MB 制造粒子云。',
        '- MP4 1600×1000、15 fps；GIF 960×600，每五帧及首/中/末关键帧取样，333 ms 预览节奏。物理时钟以逐帧 JSON 为准。视频 active 只统计处于有效保存计算区间的轨迹。停止后不再外推；这些微泡后续状态未知，不等同于物理死亡、排出或真实悬浮液中消失。历史线透明度仅是显示淡化，不改变历史样本；accumulation 只画已完成历史，不把未完成尾迹当作成功路径。',
        '- XY/XZ/YZ 累积图按原接受 dt 加权，避免近壁大量细分造成虚假的计数密度。它是计算轨迹占据图，不是超声测量重建。','',
        '## 验证与审查范围','',
        f'本轮实际运行 {total_tests} 项永久测试，全部通过；包括 P8.1 与 P1/P6.5/P7/P8/此前 full3d 相关回归，并未声称重新执行全部历史阶段测试。锁定的 {preserved} 个上游文件 SHA256 保持不变。',
        '原始尺寸、isotropic q、每 ID 出生时刻、首次入口 proposal 由保存种子重新生成逐项精确比较；每条完成轨迹末段重新分类到官方出口。另对每 admitted MB 最多九个保存的非出口点复查原 FEM 内部性，这项是抽样，非全步额外点定位。',
        '每个 MP4/GIF 已全帧解码；保存帧账本与 scene 状态逐帧比较，验证 IDs、位置、未来历史排除、完成/失败记账和时钟。重复建立两个独立 renderer 导出同一帧，RGB 与物理状态精确一致。这项像素一致性覆盖当前 software_environment.json 记录的运行环境。',
        '代理实际查看十张图及四组已导出视频的首/中/末帧/关键帧，审查说明与文件散列见 visual_inspection.json；这不是用户人工验收，也未声称人眼看完每一视频帧。用户验收状态为 PENDING_USER_REVIEW。','',
        f'代码提交：`{commit}`；所有代码与测试及产物绑定。开发中停止的试算保留在 data/reference_*，不计入正式数据量，最终以根 trajectories、trajectory_catalog 与本报告为准。','',
        '## 可直接打开的产物','',
        '- [本地可视化总览](index.html) · [验证 JSON](PARTICLE8_1_VALIDATION.json)',
        '- [轨迹目录 JSON](data/particle8_1_trajectory_catalog.json) · [逐 MB CSV](data/trajectory_catalog.csv) · [全部样本 CSV.gz](data/trajectory_samples.csv.gz)',
        '- [事件 JSON](data/particle8_1_events.json) · [事件 CSV](data/events.csv) · [物理审计](data/physics_audit.json)',
        '- [测试记录](data/test_run.json) · [媒体解码](data/media_audit.json) · [视觉检查](data/visual_inspection.json)',
        '- [完整文件列表](OUTPUT_FILES.json) · [SHA256 清单](ARTIFACT_SHA256.json)','']
    for kind in ANIMATIONS:
        name='particle8_1_anim_'+kind;lines.append(f'- {kind}：[MP4](animations/{name}.mp4) · [GIF](animations/{name}.gif) · [帧记录](frames/{name}_frames.json)')
    (root/'PARTICLE8_1_REVIEW.md').write_text('\n'.join(lines),encoding='utf-8')
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Particle-8.1 完整血管 ULM 轨迹</title><style>body{font:16px/1.7 system-ui;color:#183449;background:#f5f8fa;margin:0}main{max-width:1500px;margin:auto;padding:32px}article{background:white;padding:22px;margin:24px 0;border:1px solid #ccd9df;border-radius:8px}video,img{width:100%;height:auto}a{color:#208ca6}aside{background:#fff7f0;padding:16px;border-left:4px solid #b56743}.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px}@media(max-width:800px){.grid{grid-template-columns:1fr}main{padding:12px}}</style><main><h1>Particle-8.1 / 完整血管微泡轨迹可视化层</h1>'''
    page+=f'<p>{TITLE}</p><p><b>{cat["completed"]} completed / {cat["scheduled"]} scheduled · {cat["acquisition_birth_window_s"]/3600:.3f} h acquisition · {total_tests} tests PASS</b></p><aside>{outlet_text}。有轨迹覆盖 {covered}/3 个出口分支；完整网格显示不等于全网络重建。独立 Frozen-flow 轨迹，不是 RBC-coupled physiology。计算截尾不代表生理捕获。</aside><p><a href="PARTICLE8_1_REVIEW.md">中文审核报告</a> · <a href="PARTICLE8_1_VALIDATION.json">验证 JSON</a> · <a href="data/trajectory_catalog.csv">轨迹 CSV</a> · <a href="OUTPUT_FILES.json">文件清单</a></p>'
    for kind in ANIMATIONS:
        name='particle8_1_anim_'+kind;m=read(root/'frames'/(name+'_frames.json'));middle=m['keyframe_indices'][1]
        page+=f'<article><h2>{html.escape(kind)}</h2><video controls preload="metadata" poster="inspection/{name}_decoded_{middle:04d}.png"><source src="animations/{name}.mp4" type="video/mp4"></video><p><a href="animations/{name}.mp4">MP4</a> · <a href="animations/{name}.gif">GIF</a> · <a href="frames/{name}_frames.json">逐帧证据</a></p></article>'
    page+='<div class="grid">'+''.join(f'<article><a href="figures/particle8_1_fig_{name}.png"><img loading="lazy" src="figures/particle8_1_fig_{name}.png" alt="{name}"></a></article>' for name in FIGURES)+'</div><p>'+html.escape(INDEPENDENT)+'</p></main></html>'
    (root/'index.html').write_text(page,encoding='utf-8')
    files=[p for p in sorted(root.rglob('*')) if p.is_file() and p.name not in ['OUTPUT_FILES.json','ARTIFACT_SHA256.json'] and '__pycache__' not in str(p)]
    dump(root/'OUTPUT_FILES.json',dict(root=str(root),self_exclusion=['OUTPUT_FILES.json','ARTIFACT_SHA256.json'],
        formal_cohort='trajectories/ only; data/reference_* and data/*benchmark are development diagnostics, not additional cohort MB',
        files=[dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=digest(p)) for p in files]))
    files.append(root/'OUTPUT_FILES.json')
    dump(root/'ARTIFACT_SHA256.json',dict(self_exclusion='ARTIFACT_SHA256.json',sha256={str(p.relative_to(root)):digest(p) for p in sorted(files)}))
    return status


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-commit',required=True);args=p.parse_args()
    print(finalize(args.source_commit))
