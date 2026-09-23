#!/usr/bin/env python3
"""Finalize reviewable evidence only when tests, media and inspection agree."""
from pathlib import Path
import argparse,csv,hashlib,html,subprocess,sys
from datetime import datetime,timezone
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle8_replay import REPO,read,write,digest,canonical_hash
from particle_3d.particle8_full3d_data import OUTPUT,CASES,LABELS,verify_preservation,scene_from_name,replay_frame
from verify_particle8_full3d import source_files,artifact_files


DETAILS={
 'A_real_single_mb_orbit':('真实 Frozen / 单 MB；保存状态来自 P6.5，经 P7/P8 回放。',
  '仅入口附近 1 ms、5 个保存状态，位移约 0.440609 µm。主视图完整血管，副图保持原尺寸放大观察。',
  '不能解释为贯穿完整血管；与 B 是独立 fixture，不能拼成真实混合输运。出生前无粒子。'),
 'A_real_single_mb_fixed':('同旋转 A，物理帧逐项一致，固定相机对照。',
  '同一 1 ms 单 MB 保存轨迹；固定视角便于区分短距离位移与相机引起的投影变化。',
  '固定与旋转版都只回放已有状态，不证明新时间步收敛或全血管通行。'),
 'B_real_mixed_smoke_orbit':('真实 Frozen / mixed inlet smoke；仅绘制实际 admitted 且 active 的粒子。',
  '0–42.986308099 s。结尾 RBC scheduled=1110、admitted=active=1、pending=1109；MB scheduled=pending=1、admitted=0。',
  '唯一 admitted RBC 是保存的静止 capsule surrogate。尾迹长度为零，不根据诊断速度制造运动。pending 从不画入血管。'),
 'C_synthetic_lifecycle_orbit':('SYNTHETIC CONTROL / NOT REAL FROZEN LUMEN；P8 synthetic scene。',
  '宽 200 µm、长 10 nm 的开放 section。三个连续观察窗口与最终记账帧之间有明确 CUT 标签。',
  '最终 2356 RBC 与 2 MB admitted；2352 RBC 与 2 MB exited/deleted；4 RBC active。pending=0 是该控制案例的实际结果，未虚构排队。有限形状跨越开放端盖，不是生理通过证据。'),
 'D_restart_continuity_orbit':('SYNTHETIC CONTROL / RESTART AUDIT；continuous 与 restarted 两份 P8 scene。',
  '0.124–0.126 s，checkpoint=0.125 s，并切到 0.25 s 最终记账。两视图采用相同拉伸和相机。',
  '所有显示帧状态散列相同、mismatch=0。使用已接受的 P7 binary restart 证据缓存，本轮没有重跑物理重启或产生新解。'),
}


def require(condition,message):
    if not condition:raise RuntimeError(message)


def finalize(commit):
    root=OUTPUT;tests=read(root/'data/test_run.json');media=read(root/'data/media_audit.json');visual=read(root/'data/visual_inspection.json')
    require(tests['all_pass'],'Permanent tests failed')
    require(tests['source_sha256']==source_files(),'Source changed since tests')
    require(tests['artifact_sha256']==artifact_files(),'Rendered evidence changed since tests')
    for suite in tests['suites']:
        for key in ['junit','log']:require(digest(root/suite[key])==suite[key+'_sha256'],'Test receipt changed')
    for path,sha in source_files().items():
        saved=subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=REPO)
        require(hashlib.sha256(saved).hexdigest()==sha,'Source does not match specified code commit: '+path)
    require(media['all_pass'] and len(media['animations'])==5,'Incomplete media decode')
    require(visual['result']=='PASS' and set(visual['cases'])==set(CASES),'Missing actual visual inspection')
    for group in ['reviewed_file_sha256','reviewed_media_sha256']:
        require(bool(visual[group]),'Empty visual inspection receipt')
        for path,sha in visual[group].items():require(digest(root/path)==sha,'Visual inspection is stale: '+path)
    preserved=verify_preservation(root);manifests={};frame_rows=[];particle_rows=[]
    for case,spec in CASES.items():
        m=read(root/'frames'/('particle8_full3d_'+case+'_manifest.json'));scene=scene_from_name(root,spec['scene']);manifests[case]=m
        for f in m['frames']:
            expected=replay_frame(scene,f['physical_time_s'],spec['tail_window_s'])
            require(f['state']==expected and f['state_sha256']==canonical_hash(expected),'Replay fidelity failed')
            require(f['classification']==spec['classification'] and f['labels']==LABELS,'Scientific labels changed')
            ids=[p['particle_id'] for p in expected['particles']]
            require(f['rendered_ids_main']==f['rendered_ids_secondary']==ids,'Unexpected rendered particle')
            require(not set(ids)&set(expected['pending_ids']) and not expected['pending_drawn_ids'],'Pending drawn as active')
            if case.startswith('D_'):require(f['restart_mismatches']==0 and f['restarted_state_sha256']==f['state_sha256'],'Restart mismatch')
            row=[case,f['frame_index'],f['video_time_s'],f['physical_time_s'],f['window_label'],f['camera_main']['azimuth_deg'],
                 f['camera_secondary']['azimuth_deg'],28,f['synthetic_z_display_factor'],f['state_sha256']]
            row.extend(expected['counts'][species][key] for species in ['RBC','MB'] for key in ['scheduled','admitted','pending','active','exited','deleted'])
            frame_rows.append(row)
            for p in expected['particles']:
                particle_rows.append([case,f['frame_index'],f['physical_time_s'],p['particle_id'],p['species'],
                    *p['position_m'],*p['q'],p['geometry_sha256'],p['admitted_shape_sha256'],len(p['tail'])])
    for name,header,rows in [
        ('frame_audit.csv',['case','frame_index','video_time_s','physical_time_s','window_label','main_azimuth_deg','inset_azimuth_deg','elevation_deg','display_z_factor','state_sha256']+
         [f'{s}_{k}' for s in ['RBC','MB'] for k in ['scheduled','admitted','pending','active','exited','deleted']],frame_rows),
        ('active_particles.csv',['case','frame_index','physical_time_s','particle_id','species','x_m','y_m','z_m','q_w','q_x','q_y','q_z','original_geometry_sha256','admitted_shape_sha256','tail_knots'],particle_rows)]:
        with (root/'data'/name).open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(header);writer.writerows(rows)
    total_tests=sum(s['tests'] for s in tests['suites'])
    statuses=dict(AUTOMATED_CHECKS='PASS',VISUAL_OUTPUTS='PASS',SCIENTIFIC_LABEL='PASS',REPLAY_FIDELITY='PASS',STAGE_RESULT='PASS_WITH_LIMITATIONS')
    old=read(REPO/'particle_3d/reports/particle8/PARTICLE8_VALIDATION.json')
    record=dict(stage='Particle-8 / Full 3D visualization enhancement only',created_utc=datetime.now(timezone.utc).isoformat(),
        statuses=statuses,source_commit=commit,source_sha256=source_files(),branch=subprocess.check_output(['git','branch','--show-current'],cwd=REPO,text=True).strip(),
        source_commit_role='EXACT_SOURCE_MATCH_AND_BOUND_TO_TEST_RUN',user_acceptance='PENDING_USER_REVIEW',
        automated_tests=tests['suites'],total_tests_run_this_stage=total_tests,
        prior_particle8_accepted_tests=old['total_tests'],prior_641_tests_rerun_in_full=False,
        upstream_preserved_files=preserved,upstream_baseline_commit=read(root/'data/upstream_lock.json')['baseline_commit'],
        upstream_sources=read(root/'data/discovered_inputs.json'),full_geometry=read(root/'data/full_geometry_manifest.json'),
        display_contract=read(root/'data/display_contract.json'),media_audit=media,
        rendered_frames=len(frame_rows),gif_decoded_frames=sum(x['decoded_gif_frames'] for x in media['animations']),
        replay_frames_checked=len(frame_rows),visual_inspection=visual,
        carried_forward_limitations=old['carried_forward_limitations'],
        scientific_labels=LABELS,evidence_sha256=tests['artifact_sha256'],
        csv_sha256={name:digest(root/'data'/name) for name in ['frame_audit.csv','active_particles.csv']})
    write(root/'PARTICLE8_FULL3D_VALIDATION.json',record)
    lines=['# Particle-8 完整三维可视化增强：中文审核报告','',
        '本阶段已完成显示层代码、永久测试、五组 MP4/GIF、关键帧、storyboard 和机器记录。审核结论仅适用于已有证据的回放展示，真实 RBC 连续通行仍未建立。','',
        '```']+[f'{k} = {v}' for k,v in statuses.items()]+['```','',
        f'代码提交：`{commit}`。新输出：`{root}`。既有 {preserved} 个受锁文件 SHA256 保持一致。','',
        f'本轮实际运行 {total_tests} 项测试（新增 full3d '+str(tests['suites'][0]['tests'])+' 项，原 P8 36 项），全部通过；895 个 MP4 帧和 179 个 GIF 帧均成功解码。先前 P0–P8 的 641 项是继承证据，本轮未重跑整套物理测试。','',
        '初次检查发现离散胶囊包围盒中心不是物理中心，测试已改为解析形状表面约束；同时改进 MB 显示网格的双精度坐标构造。抽查视频发现 A 出生前误用合成 section 说明，已修正分支、补回归测试并重新导出。早期 `logs/core_tests.*` 保留为开发记录，最终以 `logs/full3d_tests.*`、`logs/upstream_particle8.*` 和 `data/test_run.json` 为准。','',
        '自动发现并复用原 P8 四个 scene、事件/生命周期账本、稳定 ID、原始几何、姿态、物理时钟和重启缓存；复用 P8 snapshot / interpolate 与 figure 配色。Frozen 模型从官方 boundary manifest 读取，主视图包含全部 45,221 个 wall 三角形以及 INLET、OUTLET_01/02/03，无抽稀和局部截取。','',
        '渲染：PyVista/VTK 离屏、正交投影；MP4 H.264/yuv420p、1600×1000、15 fps，GIF 960×600、每五帧采样、约 3 fps。GIF 时间量化用于预览，以画面 physical clock 和 MP4 帧记录为准。','',
        '相机：主视图方位角按视频帧线性从 35° 增到 115°，俯仰恒定 28°；速度约 6.2–7.5°/视频秒。真实入口副图比主视图多 90°。固定 A 主视图为 35°、副图 125°，两版本物理状态散列完全相同。相机、原点平移和 µm 单位换算均仅用于显示，不旋转物理坐标。','',
        '插值：复用 P8 线性位置与 SLERP 姿态，仅用于显示；不外推额外物理解、不改变 birth/admit/delete 时序。尾迹为保存折线加近期窗口和当前时刻的裁切插值端点，无样条平滑。默认真实 0.35 ms、合成 0.25 ms，只有 active 粒子有尾迹，删除时移除；静止 RBC 不制造非零轨迹。支持 `--tail-mode full` 另行导出，但本交付审核针对 recent。','',
        '真实粒子保持原始尺寸与姿态；capsule 使用保存的 capsule_axis。A 放大图以原半径淡色轮廓、中心标记显示 MB，让不足一个半径的运动可见。合成主图 z×20,000 且采用中心 ID glyph，均显式标记 display-only；C 副图显示未拉伸原形，D 副图同尺度显示 restarted 对照。合成中心轨迹不代表有限形状完全穿过 section。','',
        '## 逐动画审查','']
    cards=[]
    for case,spec in CASES.items():
        name='particle8_full3d_'+case;m=manifests[case];desc,scope,limit=DETAILS[case];keys=m['keyframe_indices']
        links=' · '.join(f'[{title}](keyframes/{name}_decoded_{i:04d}.png)' for title,i in zip(['首帧','中帧','末帧'],keys))
        lines += [f'### {case}','',f'- 数据源：`{m["source_scene"]["path"]}`；{desc}',f'- 物理范围：{scope}',
            '- 显示插值：是，仅保存状态之间。',f'- 相机：{spec["camera"]}；尾迹：recent，{spec["tail_window_s"]*1000:g} ms。',
            f'- 限制：{limit}',f'- 输出：[MP4](mp4/{name}.mp4) · [GIF](gif/{name}.gif) · [storyboard](storyboard/{name}_storyboard.png)。',
            f'- 关键帧：{links}。','']
        thumbs=''.join(f'<a href="keyframes/{name}_decoded_{i:04d}.png"><img loading="lazy" src="keyframes/{name}_decoded_{i:04d}.png" alt="{case} frame {i}"></a>' for i in keys)
        cards.append(f'<article><h2>{html.escape(case)}</h2><p class="tag">{html.escape(spec["classification"])}</p><p>{desc} {scope}</p><p>{limit}</p>'
            f'<video controls preload="metadata" poster="keyframes/{name}_decoded_{keys[1]:04d}.png"><source src="mp4/{name}.mp4" type="video/mp4"></video>'
            f'<p><a href="mp4/{name}.mp4">MP4</a> · <a href="gif/{name}.gif">GIF</a> · <a href="storyboard/{name}_storyboard.png">三帧 storyboard</a> · <a href="frames/{name}_manifest.json">逐帧记录</a></p><div class="frames">{thumbs}</div></article>')
    lines+=['## 保留的科学边界','',
        'Real RBC passage remains unresolved。Synthetic control 仅用于 bookkeeping / lifecycle；真实 smoke 的静止 RBC 不能解释为连续输运。Production timestep 与 neighbor settings 未冻结。近场仍为 sphere-normal regularization，非球形 lubrication、切向/旋转耦合与完整多体 mobility 未冻结。H_D=0.45 是 feed volume fraction，不是瞬时 tube hematocrit；tube Hct 仅为独立诊断。没有 Particle-7.5、CFD、新 RBC 变形物理、full suspension 或 PK。原始 P8 全部 carry-forward limitations 原样保存在验证 JSON。','',
        '## 审核证据与范围','',
        '所有 MP4/GIF 已逐帧解码；逐帧数值审计核对位置、q、几何散列、stable IDs、pending 排除与 restart 一致性。代理视觉检查查看每组实际视频解码的首/中/末帧 contact sheet、中帧原尺寸图片、总览和科学范围图；未声称人工逐帧看完整视频，也不代替用户验收。用户验收状态为 PENDING_USER_REVIEW。查看文件和视频 SHA256、每组观察结论见 `data/visual_inspection.json`。','',
        '- [本地浏览与 storyboard 页面](index.html)',
        '- [总览 PNG](storyboard/full3d_storyboard.png) · [科学范围总结图](storyboard/scientific_scope_summary.png)',
        '- [机器验证](PARTICLE8_FULL3D_VALIDATION.json) · [测试执行记录](data/test_run.json) · [解码审计](data/media_audit.json)',
        '- [逐帧 CSV](data/frame_audit.csv) · [active 原始坐标/姿态 CSV](data/active_particles.csv)',
        '- [完整文件清单](OUTPUT_FILES.json) · [产物 SHA256](ARTIFACT_SHA256.json)','']
    (root/'PARTICLE8_FULL3D_REVIEW.md').write_text('\n'.join(lines),encoding='utf-8')
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Particle-8 完整三维回放审核</title><style>
    body{font:16px/1.65 system-ui,sans-serif;color:#243d50;background:#f7f9fb;margin:0}main{max-width:1400px;margin:auto;padding:28px}h1,h2{line-height:1.3}article{background:white;border:1px solid #dce4ea;border-radius:8px;padding:24px;margin:24px 0}.tag{font-weight:bold;color:#a74638}video{width:100%;max-height:850px;background:white}a{color:#14758c}.frames{display:grid;grid-template-columns:repeat(3,1fr);gap:8px}.frames img{width:100%}.summary{width:100%;max-width:1200px}aside{border-left:4px solid #b35647;padding:12px 18px;background:#fff8f6}@media(max-width:700px){.frames{grid-template-columns:1fr}main{padding:12px}article{padding:12px}}
    </style><main><h1>Particle-8 | 完整三维回放审核</h1><p>真实 Frozen 几何与合成控制分开显示；五组动画、首中末关键帧与可追溯状态记录。</p><aside>真实 RBC 连续通行尚未建立。真实 smoke 只有一个静止 admitted RBC；pending 仅计数。合成 10 nm section 仅为生命周期记账演示，z 拉伸与插值均为显示用途。</aside>'''
    page+=f'<p>{total_tests} tests PASS · {len(frame_rows)} MP4 frames decoded · PASS_WITH_LIMITATIONS · 用户验收待审</p>'
    page+='<p><a href="PARTICLE8_FULL3D_REVIEW.md">中文审核报告</a> · <a href="PARTICLE8_FULL3D_VALIDATION.json">机器记录</a> · <a href="OUTPUT_FILES.json">完整文件清单</a> · <a href="data/frame_audit.csv">逐帧 CSV</a></p>'
    page+=''.join(cards)+'<h2>总览与科学范围</h2><a href="storyboard/full3d_storyboard.png"><img class="summary" src="storyboard/full3d_storyboard.png" alt="五组动画首中末帧总览"></a><img class="summary" src="storyboard/scientific_scope_summary.png" alt="科学边界总结"><p>'+html.escape(' '.join(LABELS))+'</p></main></html>'
    (root/'index.html').write_text(page,encoding='utf-8')
    files=[p for p in sorted(root.rglob('*')) if p.is_file() and p.name not in ['OUTPUT_FILES.json','ARTIFACT_SHA256.json']]
    write(root/'OUTPUT_FILES.json',dict(root=str(root),self_exclusion=['OUTPUT_FILES.json','ARTIFACT_SHA256.json'],
        files=[dict(path=str(p.relative_to(root)),absolute_path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for p in files]))
    files.append(root/'OUTPUT_FILES.json')
    write(root/'ARTIFACT_SHA256.json',dict(self_exclusion='ARTIFACT_SHA256.json',sha256={str(p.relative_to(root)):digest(p) for p in sorted(files)}))
    return statuses


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-commit',required=True);args=p.parse_args()
    for key,value in finalize(args.source_commit).items():print(f'{key} = {value}')
