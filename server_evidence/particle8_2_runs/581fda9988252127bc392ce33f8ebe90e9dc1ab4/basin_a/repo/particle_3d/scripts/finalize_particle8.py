#!/usr/bin/env python3
"""Audit replay/media and bind the Chinese review to source and artifact hashes."""
from pathlib import Path
import argparse,json,sys,subprocess,xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import imageio_ffmpeg
from particle_3d.particle8_replay import REPO,P7,DEFAULT_OUTPUT,read,write,digest,ASSUMPTIONS,LIMITATIONS,snapshot
from particle_3d.particle8_visuals import FIGURES,ANIMATIONS


def audit_media(root):
    rows=[]
    for name in ANIMATIONS:
        path=root/'animations'/(name+'.mp4');reader=imageio_ffmpeg.read_frames(str(path),pix_fmt='rgb24')
        meta=next(reader);n=0;previous=None;changed=0
        for pixels in reader:
            if previous is not None and pixels!=previous:changed+=1
            previous=pixels;n+=1
        m=read(root/'data'/(name+'_frames.json'))
        ok=n==m['frame_count'] and meta['size']==(1440,900) and meta['fps']==15 and changed>30
        rows.append(dict(path=str(path.relative_to(root)),sha256=digest(path),decoded_frames=n,
            changed_frame_pairs=changed,duration_s=n/15,width=meta['size'][0],height=meta['size'][1],fps=meta['fps'],
            codec=meta['codec'],full_decode_pass=ok))
    write(root/'data/media_decode_audit.json',dict(encoder_version=imageio_ffmpeg.get_ffmpeg_version(),
        imageio_ffmpeg_version=imageio_ffmpeg.__version__,animations=rows,all_pass=all(x['full_decode_pass'] for x in rows)))
    return rows


def finalize(root,final_gate=False,decode=False):
    root=Path(root);data=root/'data';logs=root/'logs';p7=read(P7/'PARTICLE7_VALIDATION.json')
    media=audit_media(root) if decode or not (data/'media_decode_audit.json').exists() else read(data/'media_decode_audit.json')['animations']
    for m in media: assert digest(root/m['path'])==m['sha256'],'Media changed since decode'
    s=read(data/'synthetic_scene.json');real=read(data/'real_mixed_scene.json');single=read(data/'real_single_mb_scene.json')
    counts=snapshot(s,s['time_range_s'][1])['counts'];realcounts=snapshot(real,real['time_range_s'][1])['counts']
    queue=read(data/'pending_identity.json');restart=read(data/'restart_exact_fields.json')
    sched=np.genfromtxt(data/'03_04_scheduler.csv',delimiter=',',names=True)
    sched_pass=bool(np.all(sched['rbc_residual_m3']>=-1e-28) and np.all(sched['rbc_residual_m3']<sched['next_rbc_volume_m3']))
    mb_pass=bool(np.all(sched['mb_error']>-1) and np.all(sched['mb_error']<=1e-12))
    time_test=read(data/'08_timestep_parity.json')
    tests=[]
    for path in sorted(logs.glob('final_*.xml')):
        suites=list(ET.parse(path).getroot().iter('testsuite'))
        tests.append(dict(suite=path.stem.removeprefix('final_'),tests=sum(int(x.attrib['tests']) for x in suites),
            failures=sum(int(x.attrib.get('failures',0))+int(x.attrib.get('errors',0)) for x in suites),
            skipped=sum(int(x.attrib.get('skipped',0)) for x in suites),junit=str(path.relative_to(root))))
    commands=read(logs/'test_commands.json') if (logs/'test_commands.json').exists() else []
    inspection=read(data/'visual_inspection.json') if (data/'visual_inspection.json').exists() else dict(status='PENDING',reviewer='CODEX_AGENT')
    required_visuals={str(Path('figures')/(n+'.png')) for n in FIGURES}
    required_visuals.update(str(Path('animations')/(n+'.mp4')) for n in ANIMATIONS)
    required_visuals.update(str(Path('inspection')/(n+'_contact_sheet.png')) for n in ANIMATIONS)
    reviewed=inspection.get('reviewed_artifact_sha256',{})
    inspection_valid=(inspection['status']=='PASS' and required_visuals<=set(reviewed)
        and all((root/p).exists() and digest(root/p)==h for p,h in reviewed.items()))
    passed=(len(tests)==11 and all(not r['failures'] and not r['skipped'] for r in tests)
        and len(commands)==13 and all(c['returncode']==0 for c in commands)
        and sched_pass and mb_pass and queue['identity_exact'] and all(restart['exact_match_fields'].values())
        and all(m['full_decode_pass'] for m in media) and inspection_valid)
    if final_gate and not passed: raise ValueError('Final gate not satisfied: tests, media, or visual receipt incomplete')
    source_paths=[]
    for folder in ['src/particle_3d','scripts','tests/particle8']:
        source_paths.extend(p for p in (REPO/'particle_3d'/folder).rglob('*.py') if 'particle8' in p.name or 'tests/particle8' in str(p))
    source_paths.extend([REPO/'particle_3d/PARTICLE8_README.md',REPO/'particle_3d/requirements-particle8.txt'])
    sources={str(p.relative_to(REPO)):digest(p) for p in sorted(set(source_paths))}
    git=lambda *args:subprocess.check_output(['git',*args],cwd=REPO,text=True).strip()
    outputs=dict(figures=[str(Path('figures')/(n+'.png')) for n in FIGURES],
        animations_mp4=[str(Path('animations')/(n+'.mp4')) for n in ANIMATIONS],
        animations_gif=[str(Path('animations')/(n+'.gif')) for n in ANIMATIONS])
    final=dict(stage_name='Particle-8 / Continuous mixed-flow visualization and replay',
        git_commit=git('rev-parse','HEAD'),git_commit_role='TESTED_SOURCE_COMMIT_WHEN_SOURCE_HASHES_MATCH_COMMIT',git_branch=git('branch','--show-current'),
        upstream_dependencies=read(data/'upstream_provenance.json'),assumptions=ASSUMPTIONS,
        real_geometry=dict(inlet_flux_m3_s=p7['inlet_Q_m3_s'],flux_balance_residual=p7['mass_balance'],
            real_single_mb_entry_pass=p7['isolated_real_mb_transport']['status']=='PASS',
            real_single_mb_displacement_m=p7['isolated_real_mb_transport']['displacement_m'],
            real_single_mb_time_span_s=single['time_range_s'],real_mixed_inlet_smoke_counts=realcounts,
            real_rbc_admission_established=False,real_rbc_passage_established=False,
            admitted_rbc_interpretation='ONE_CAPSULE_SURROGATE_IS_NOT_ESTABLISHED_CONTINUOUS_ADMISSION_OR_PASSAGE',
            mixed_smoke_motion='STATIC_ADMITTED_OBSTACLE_NO_TRANSPORT'),
        synthetic_control=dict(lifecycle_demo_pass=all(c['scheduled']==c['pending']+c['active']+c['deleted'] for c in counts.values()),
            insert_delete_pass=p7['outlet_bookkeeping_status']=='PASS' and p7['restart']['status']=='PASS',
            no_duplicate_id_pass=len(s['records'])==len(set(r['particle_id'] for r in s['records'])),
            no_particle_loss_pass=all(c['scheduled']==c['pending']+c['active']+c['deleted'] for c in counts.values()),
            counts=counts,control_geometry=s['control_geometry'],role='P7_ACTUAL_LAMMPS_COMMON_PLUG_OPEN_SECTION_REPLAY'),
        scheduler=dict(rbc_volume_scheduler_pass=sched_pass,mb_scheduler_pass=mb_pass,
            timestep_independence_pass=all(x['exact_equal'] for x in time_test['partitions']),
            physical_times_unchanged=True,production_timestep_selection=False),
        pending_queue=dict(identity_preserved_pass=queue['identity_exact'],geometry_preserved_pass=queue['identity_exact'],orientation_preserved_pass=queue['identity_exact'],pending_lumen_positions_invented=False),
        restart=dict(restart_parity_pass=all(restart['exact_match_fields'].values()),exact_match_fields=restart['exact_match_fields'],
            checkpoint_time_s=.125,actual_destroy_and_binary_restore=restart['actual_destroy_and_binary_restore'],
            full_rng_states_exact=restart['full_rng_states_exact']),visualization_outputs=outputs,
        visualization_policy=dict(positions='SAVED_STATES_LINEAR_DISPLAY_INTERPOLATION_NO_EXTRAPOLATION',
            synthetic_time='THREE_PHYSICAL_WINDOWS_AND_EXPLICIT_CUT_TO_FINAL_AUDIT',
            synthetic_shape_display='STRETCHED_Z_CENTER_PATH_GLYPHS_PLUS_TRUE_SCALE_XY_SHAPE_PROJECTION',
            flux_time='DISPLAY_CLOCK_ONLY_NO_PHYSICAL_TIMESTAMPS_FOR_DIAGNOSTIC_SAMPLES',
            real_mixed='NO_INVENTED_RBC_MOTION_NO_PENDING_PARTICLES_PLACED_IN_LUMEN'),
        carried_forward_limitations=LIMITATIONS,automated_checks_pass=bool(passed),
        manual_visual_review='PASS' if inspection_valid else 'PENDING',
        manual_visual_review_role='AGENT_VISUAL_INSPECTION_OF_ALL_FIGURES_AND_DECODED_VIDEO_CONTACT_SHEETS_NOT_USER_ACCEPTANCE',
        user_acceptance='PENDING_USER_REVIEW',stage_result='PASS_WITH_CARRY_FORWARD_LIMITATIONS' if passed else 'PENDING_FINAL_GATE',
        tests=tests,total_tests=sum(x['tests'] for x in tests),media_decode_audit=media,source_sha256=sources,
        data_sha256={str(p.relative_to(root)):digest(p) for p in sorted(data.rglob('*')) if p.is_file()},
        figure_sha256={p:digest(root/p) for p in outputs['figures']},
        animation_sha256={p:digest(root/p) for p in outputs['animations_mp4']+outputs['animations_gif']})
    write(root/'PARTICLE8_VALIDATION.json',final)
    table='\n'.join(f'| {r["suite"]} | {r["tests"]} | {"PASS" if not r["failures"] else "FAIL"} |' for r in tests)
    animation_notes=[
      '真实 Frozen 入口；单 MB 在原入口面出生，无人工偏移，按原 P6.5 保存状态在 1 ms 内移动 0.440609 µm。左右分别为局部 3D 和等物理尺度投影。开始前 0.1 ms 留空。中间帧只是显示插值，没有新增动力学证据。',
      '真实 Frozen 混合入口接纳 smoke；时钟推进、预约队列增长，唯一已接纳的 P3 胶囊替代 RBC 保持静止。最终 RBC 1110/1/1109、MB 1/0/1（预约/接纳/排队）。队列在图表中显示，不伪造其空间位置。',
      '合成控制区，宽 200 µm、长 10 nm，公共速度 25 µm/s，中心驻留时间 0.4 ms。三段真实物理时间窗口展示启动、第一只及第二只 MB；窗口之间明示剪切，末帧切到 0.25 s 记账。左图拉伸主流向并用 ID 图标表示中心，右图保留原尺寸、原姿态的 xy 投影。有限形状跨越开口端面，不能称为完整粒子通过血管。',
      '独立的连续与重启保存数据分屏回放，跨越 0.125 s checkpoint。原 P7 确实销毁 LAMMPS 并读取二进制；本阶段读已有证据，不重新演出或生成随机历史。逐帧位置、ID、队列与计数一致，完整 scheduler/RNG 等字段也逐项比较。',
      '真实入口几何上的采样诊断。100000 个原 P7 候选点逐步累积，按三角片面积归一化为密度。10 s 是播放时钟，不是物理出生时钟；这些候选点不是已接纳粒子。']
    anim_text='\n\n'.join(f'### 动画 {i+1:02d}\n\n[MP4]({outputs["animations_mp4"][i]}) · [GIF 预览]({outputs["animations_gif"][i]})\n\n{note}' for i,note in enumerate(animation_notes))
    figure_notes=[
      '两条证据线和共用回放层分开；真实几何结果没有被合成动画代替。',
      '真实入口速度、三角片通量与面积归一化采样密度并列；质量差保留原值。',
      '累计目标与每只原 RBC 的体积之和一致，残差小于下一只原样本体积。',
      'MB 在累计数量跨整数时预约，误差严格大于 −1；不是 Poisson。',
      '真实入口排队达到 1109 个 RBC 和 1 个 MB；静态已接纳形状与排队计数分开。',
      '合成生命周期闭合，最终 4 个 RBC 活动；退出与删除同步，ID 不复用。',
      '等待前后 ID、几何、体积、四元数、来源完全一致，没有重抽。',
      '六类重启字段全部相等；0.125 s 后未来事件保持一致。',
      '从实际动画场景抽取关键帧；不是另外绘制的示意结果。',
      '真实 RBC passage、生产参数、完整悬浮液、PK 等边界醒目标明。']
    fig_text='\n\n'.join(f'### 图 {i:02d}\n\n![{FIGURES[i]}]({outputs["figures"][i]})\n\n审核观察：{note}' for i,note in enumerate(figure_notes))
    total=final['total_tests'];auto='PASS' if passed else 'PENDING_FINAL_GATE';manual=final['manual_visual_review']
    report=f'''# Particle-8 连续混合入口可视化与回放审核报告

## 1. 阶段定义

Particle-8 是可视化与回放层：把原始事件、稳定 ID、形状、姿态、物理时刻、pending、接纳、活动、出口、删除和重启记录变成可审查的图和动画。本轮新增 replay/scene/frame 导出、编码、自动验证及报告，未修改 P0–P7 科学代码，没有启动 Particle-7.5 或新的 RBC 变形模型。

## 2. 继承依赖

用户指定的 `/mnt/data/PARTICLE6_5_*` 和 `/mnt/data/PARTICLE7_*` 在本 WSL 中没有副本，采用仓库中相应的已提交文件。上游文件路径和 SHA256 均在 `data/upstream_provenance.json`。P7 的 285 个源码文件和 53 个数据/来源文件逐项核对；不改写上游历史报告中的审核快照。当前 Particle-8 指令明确把 P6.5、P7 作为已接受基础。

沿用 H_D=0.45（feed volume fraction）、C_MB=8.5e12 m⁻³（名义浓度估算）、确定性累计通量调度、暂定各向同性 SO(3) 姿态。沿用原 P2 RBC 与 SonoVue 样本、原 P6.5 sphere-normal near-field、P7 FIFO 和 LAMMPS 二进制重启。

## 3. 关键结果

真实入口 Q={p7['inlet_Q_m3_s']:.12g} m³/s；有符号入口减出口差为 {p7['mass_balance']['signed_in_minus_out_m3_s']:.6g} m³/s，相对差 {p7['mass_balance']['relative_signed_residual']:.6g}。保留原 Frozen 数值，未重新运行 CFD。

RBC 目标体积为 0.45×累计全血体积，逐只使用原采样体积；MB 按累计期望数量跨整数生成预约。调度物理时间与视频帧率无关。pending 前后几何、身份、姿态和来源不变，排队不产生虚假的管腔内位置。

实际 LAMMPS 合成回放记录有 2358 个预约/接纳事件：2356 RBC、2 MB。退出并删除 2352 RBC 和 2 MB，最终 4 RBC 活动；没有重复 ID 或粒子丢失。checkpoint 在 0.125 s、1180 个事件时保存；销毁恢复后的未来事件、完整 RNG、下一只 RBC、scheduler 与活动状态均与连续运行完全一致。

独立真实单 MB 保留 5 个接受状态，1 ms 位移 0.440609 µm。真实混合 smoke 接纳 1 个 RBC 胶囊替代形状，1109 RBC 和 1 MB 排队，未发生出口事件；这不能证明真实 RBC 连续 admission 或 passage。

已输出 10 张 PNG、5 个 H.264 MP4 和 5 个 GIF 预览。MP4 为 1440×900、15 fps，逐帧完整解码验证；GIF 是 3 fps 缩小预览。视频时间不是生产 timestep。每部动画有带物理/显示时钟角色、场景状态、ID、计数和状态 SHA256 的逐帧 manifest。

## 4. 两条证据线应该怎样理解

证据线 A 使用真实 Frozen 边界和已保存的数值轨迹，负责回答真实几何上目前观察到了什么。它诚实显示拥堵，而不是把 RBC 画成通过了血管。

证据线 B 使用 P7 合成开口控制区，负责展示持续预约、接纳、移动、出界、删除和重启这些机制。10 nm 控制区比粒子短，有限形状跨越开口；中心穿过控制区不是整个形状穿过生理血管。图中主流向拉伸只影响显示，真实位置始终按 SI 保存；原尺度的横截面投影并列保留。

真实混合 smoke 中，既有形状作为静态障碍保留；原始尺寸的后续粒子无法通过有限尺寸接纳检查就必须等待。保留 pending 是系统正确处理限制的表现，不能用缩小 RBC、降低 H_D 或重抽小粒子掩盖。它说明入口基础设施能诚实报告限制，不说明真实混合输运已完成。

Tube Hct 是独立诊断。原真实全腔中心归属快照约 3.273%，并没有强制为 45%；跨开口粒子按中心归属计完整体积，也不等于几何交集体积分数。

## 5. 动画审核

{anim_text}

## 6. 静态图审核

{fig_text}

## 7. 风险与限制

- 真实 RBC 连续接纳与穿管仍未建立；一个 P3 胶囊替代形状接纳不能外推成功穿管。
- Synthetic control 不是真实生理管腔证明；真实 smoke 没有混合输运。
- 非球形润滑、生产 timestep、neighbor cutoff/skin 均未冻结。
- P6.5 handoff 事件时刻与事件拓扑收敛未建立；事件计数仍是求解器诊断。
- 当前只有球形法向近场修正，没有完整多体远场流体动力学、切向/旋转润滑耦合或 RBC 润滑。
- 没有糖萼、黏附、表面粗糙度、真实分子接触或壳层力学模型；未建立 full suspension。
- LAMMPS 沿用单 MPI rank 的存储、邻居和重启用途，不承担粒子力积分。
- MB 浓度和 RBC 各向同性取向是模型假设，非该小鼠直接测量；没有完整 PK、清除或破坏模型。
- 没有 CFD、没有 Particle-7.5，也没有以可视化为由修改物理。
- 单 MB 中间视频帧为原接受状态之间的线性显示插值，不是新增求解状态或独立的连续壁间隙证明。
- 合成动画有明确时间剪切，流向坐标拉伸，左图粒子符号不按物理尺寸；不能据此衡量真实血管速度或粒径。

## 8. 自动测试与复现

| 测试组 | 数量 | 结果 |
|---|---:|---|
{table}

合计 {total} 项。完整命令和 JUnit 在 `logs/test_commands.json` 和各 `final_*.xml` 中。P8 永久测试验证缓存重导出一致、原几何/姿态逐项保留、精确出口删除、禁止外推和虚假 pending 位置、原轨迹回放、重启全状态一致、每帧再生、PNG 来源映射、MP4 完整解码、连续运动和 GIF 输出。开发中修复了相同状态插值引入浮点微扰的问题，并保留原失败日志及永久回归。

复现入口见 [PARTICLE8_README.md](../../PARTICLE8_README.md)。只重新作图无需重跑 P7 模拟；全部输出路径在 [OUTPUT_FILES.json](OUTPUT_FILES.json)，浏览入口为 [index.html](index.html)。

## 9. 审核结论

AUTOMATED_CHECKS = {auto}

MANUAL_VISUAL_REVIEW = {manual}

STAGE_RESULT = {final['stage_result']}

视觉检查角色：代理逐张查看 PNG，并查看从每个实际编码 MP4 解码抽取的联系表，另以全部帧解码和状态比较检查连续性。此处 PASS 仅指代理视觉检查完成，不代表用户已作人工验收；USER_ACCEPTANCE = PENDING_USER_REVIEW。

通过依据是回放层完成、全部动画可解码、调度/队列/生命周期/重启证据保持完整、分类及限制明确；不是 full scientific pass。

分支 `{final['git_branch']}`；代码证据绑定 `{final['git_commit']}`，逐文件源码 SHA256 为准。用户未要求 push 或合并，本阶段只作本地交付。
'''
    (root/'PARTICLE8_REVIEW.md').write_text(report)
    # Local browser index makes every primary deliverable accessible without a service.
    html=['<!doctype html><html lang="zh"><meta charset="utf-8"><title>Particle-8 回放审核</title>',
      '<style>body{font:17px system-ui;background:#f6f8fa;color:#183449;max-width:1100px;margin:40px auto;padding:0 24px}video,img{width:100%;background:white}article{padding:22px;background:white;margin:24px 0;border-radius:10px}h1{font-size:32px}a{color:#147e9d}.limit{border-left:5px solid #d85b47;padding:18px;background:#fff}</style>',
      '<h1>Particle-8：连续入口可视化与回放</h1><p class="limit">真实 RBC 连续通行仍未建立。合成控制区动画验证生命周期；真实混合 smoke 展示静态接纳障碍和排队。H_D=0.45 为入口体积流率比例，C_MB 为名义浓度假设。</p>',
      '<p><a href="PARTICLE8_REVIEW.md">中文审核报告</a> · <a href="PARTICLE8_VALIDATION.json">机器验证</a> · <a href="OUTPUT_FILES.json">全部路径</a></p>']
    for i,n in enumerate(ANIMATIONS):
        html.append(f'<article><h2>动画 {i+1:02d}</h2><p>{animation_notes[i]}</p><video controls preload="metadata" src="animations/{n}.mp4"></video><p><a href="animations/{n}.mp4">MP4</a> · <a href="animations/{n}.gif">GIF</a> · <a href="data/{n}_frames.json">逐帧记录</a></p></article>')
    for i,n in enumerate(FIGURES):html.append(f'<article><h2>图 {i:02d}</h2><p>{figure_notes[i]}</p><a href="figures/{n}.png"><img loading="lazy" src="figures/{n}.png" alt="{n}"></a></article>')
    html.append('</html>');(root/'index.html').write_text('\n'.join(html))
    files=[dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(root.rglob('*'))
           if p.is_file() and p.name not in ['OUTPUT_FILES.json','ARTIFACT_SHA256.json'] and 'logs' not in p.relative_to(root).parts]
    write(root/'OUTPUT_FILES.json',dict(root=str(root),outputs=files,implementation_source_files=sources,
        log_files=[str(p.relative_to(root)) for p in sorted(logs.rglob('*')) if p.is_file()],
        inventory_and_manifest=['OUTPUT_FILES.json','ARTIFACT_SHA256.json'],
        logs_directory=str(logs),reproduction_readme=str(REPO/'particle_3d/PARTICLE8_README.md')))
    write(root/'ARTIFACT_SHA256.json',{f['path']:f['sha256'] for f in files})
    print(final['stage_result'],total,'tests;',len(files),'artifacts',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT);p.add_argument('--final-gate',action='store_true');p.add_argument('--decode',action='store_true');a=p.parse_args()
    finalize(a.output,a.final_gate,a.decode)
