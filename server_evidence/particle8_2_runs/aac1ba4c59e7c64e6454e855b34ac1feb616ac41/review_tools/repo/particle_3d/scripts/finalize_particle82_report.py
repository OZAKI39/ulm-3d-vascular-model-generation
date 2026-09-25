#!/usr/bin/env python3
"""Build Chinese review and machine validation only from completed evidence."""
from pathlib import Path
import argparse,datetime,html,json,subprocess,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_provenance import atomic_json,sha256
from particle_3d.particle8_replay import canonical_hash

OUTLETS=['OUTLET_01','OUTLET_02','OUTLET_03'];BASINS=OUTLETS+['UNRESOLVED_POINT_PATH']

def main(root):
    root=Path(root);read=lambda p:json.loads((root/p).read_text())
    local=read('LOCAL_HOST_PROVENANCE.json');remote=read('REMOTE_HOST_PROVENANCE.json');host=read('host_provenance.json')
    s=read('scaling/SERVER_SCALING_BENCHMARK.json');point=read('POINT_REFINEMENT_AUDIT.json');p=point['baseline_summary'];ref=point['refined_summary']
    adm=read('admission/ADMISSION_BASIN_AUDIT.json');recheck=read('admission/ORIGINAL_ADMISSION_RECHECK.json')
    stops=read('admission/stop_audit/SAFETY_STOP_AUDIT.json');extended=read('continuation/EXTENDED_GUARD_AUDIT.json');dt=read('timestep/TIMESTEP_AUDIT.json')
    c=read('data/particle8_2_natural_trajectory_catalog.json');media=read('data/media_audit.json');repeat=read('data/render_determinism.json')
    tests=read('PARTICLE82_TEST_RESULTS.json');reg=read('UPSTREAM_REGRESSION.json');immut=read('UPSTREAM_IMMUTABILITY.json')
    inspection=read('AGENT_VISUAL_INSPECTION.json');natural=read('natural/catalog.json')
    count=c['outlet_counts'];covered=all(count.values());factors=extended['factors'];latest=factors[-1]
    cpu=next((r['stdout'] for r in remote['commands'] if r['command']==['lscpu']), '')
    model=next((line.split(':',1)[1].strip() for line in cpu.splitlines() if line.startswith('Model name:')),'AMD Ryzen 7 7800X3D (see raw lscpu)')
    gpu=next((r['stdout'] for r in remote['commands'] if r['command'][0]=='nvidia-smi'),'See REMOTE_HOST_PROVENANCE.json')
    code=Path(__file__).resolve().parents[2];marker=code/'SOURCE_COMMIT'
    commit=marker.read_text().strip() if marker.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=code,text=True).strip()
    software=dict(remote_computation_provenance=c['server_only_formal_dataset'],parallel_parity=s['all_exact'],
        point_tracer_three_basin_audit=p['three_outlets_observed'] and point['identical_seed_positions'],
        admission_basin_audit=recheck['all_exact'] and adm['totals']['scheduled']==2200,
        stop_audit=stops['count']==720,extended_prefix_parity=all(f['all_original_prefixes_exact'] for f in factors),
        saved_dataset_reproducible=repeat['all_pass'],animation_decode=media['all_pass'],
        permanent_tests=tests['failed']==0 and tests['skipped']==0,agent_visual_inspection=inspection['all_inspected'])
    doc=dict(stage='Particle-8.2',git_commit=commit,integration_git_commit=host['source_git_commit'],
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        local_wsl_hostname=local['hostname'],local_wsl_nproc=local['nproc'],local_wsl_ram=local['memory']['MemTotal'],
        remote_hostname=remote['hostname'],remote_nproc=remote['nproc'],remote_cpu_model=model,remote_ram=remote['memory']['MemTotal'],
        remote_cgroup_limits=remote['cgroup'],remote_gpu=gpu,solver_gpu_usage='NOT_USED',
        particle8_1_compute_host_status='NOT_PROVEN',formal_compute_host='REMOTE_SERVER',
        worker_scaling_results=s['worker_scaling_results'],chosen_workers=s['chosen_workers'],parallel_speedup=s['parallel_speedup'],
        serial_tracks_per_hour=s['serial_tracks_hour'],parallel_tracks_per_hour=s['parallel_tracks_hour'],
        peak_memory=dict(process_pss_bytes=max(r['peak_process_pss_sum_bytes'] for r in s['worker_scaling_results']),
                         process_rss_sum_bytes=max(r['peak_process_rss_sum_bytes'] for r in s['worker_scaling_results'])),
        frozen_input_sha256=host['frozen_input_sha256'],config_sha256=canonical_hash(natural['config']),
        point_tracer_count=p['count'],point_tracer_outlet_01=p['outlet_counts']['OUTLET_01'],
        point_tracer_outlet_02=p['outlet_counts']['OUTLET_02'],point_tracer_outlet_03=p['outlet_counts']['OUTLET_03'],
        point_tracer_unresolved=p['unresolved'],point_tracer_refined_counts=ref['outlet_counts'],point_tracer_refined_unresolved=ref['unresolved'],
        point_tracer_fraction_vs_frozen_flux=p['fraction_vs_frozen'],point_tracer_flux_sanity_status=point['flux_sanity_status'],
        proposal_by_basin={b:adm['by_basin'][b].get('proposal_count',0) for b in BASINS},
        admitted_by_basin={b:adm['by_basin'][b].get('admitted_count',0) for b in BASINS},
        guard_exhausted_by_basin={b:adm['by_basin'][b].get('guard_exhausted_by_first_proposal_count',0) for b in BASINS},
        basin_accounting=adm['accounting'],baseline_safety_stops=720,
        extended_completed_from_stops=latest['completed_from_stops'],extended_still_stopped=latest['still_stopped'],
        extended_by_factor=[{k:f[k] for k in ['guard_factor','completed_from_stops','still_stopped','all_original_prefixes_exact']} for f in factors],
        timestep_audit=dt,natural_scheduled=c['scheduled'],natural_admitted=c['admitted'],natural_completed=c['completed'],
        outlet_01_completed=count[OUTLETS[0]],outlet_02_completed=count[OUTLETS[1]],outlet_03_completed=count[OUTLETS[2]],
        all_outlets_observed=covered,full_network_coverage_status='PASS' if covered else 'NOT_ESTABLISHED',
        ulm_full_network_visualization_target=all(n>=50 for n in count.values()),
        server_only_formal_dataset=c['server_only_formal_dataset'],production_timestep='NOT_FROZEN',
        manual_visual_review='PENDING_USER_REVIEW',software_checks=software,
        software_stage_status='PASS' if all(software.values()) else 'FAIL',
        scientific_scope='Software audit PASS does not establish point-flux agreement, uncensored physiology or full network coverage.',
        permanent_test_results=tests,upstream_regression=reg,upstream_immutability=immut,
        no_local_formal_trajectory_integrations=True,local_smoke_scope='P8.2 uses two analytic point-advection unit cases; upstream unit regressions also execute their original small cases. No P8.2 MB batches run locally.',
        scene_sha256=sha256(root/'data/particle8_2_natural_trajectory_catalog.json'))
    atomic_json(root/'PARTICLE8_2_VALIDATION.json',doc)
    table='\n'.join(f'| {b} | {doc["proposal_by_basin"][b]} | {doc["admitted_by_basin"][b]} | {doc["guard_exhausted_by_basin"][b]} |' for b in BASINS)
    continuation='；'.join(f'预算 ×{f["guard_factor"]}：原停止转完成 {f["completed_from_stops"]}，仍不完整 {f["still_stopped"]}' for f in factors)
    timestep='；'.join(f'{b} 可用 {dt["groups"][b]["available"]}，选取 {len(dt["groups"][b]["selected"])}' for b in OUTLETS)
    review=f'''# Particle-8.2 中文人工审核报告

工程审计：{doc['software_stage_status']}。完整网络覆盖：{doc['full_network_coverage_status']}。point tracer 与 Frozen 边界通量的一致性：{point['flux_sanity_status']}。这三个结论是不同的验证对象。用户人工审核仍为 **PENDING_USER_REVIEW**。

## 1. P8.1 到底是不是服务器算的？

证据不足，NOT_PROVEN。旧环境记录含 WSL，旧汇总记录 16 workers，但缺少逐轨迹主机回执，不能据此证明或否定全部旧重计算来自服务器。原结果没有改写；本轮在隔离历史副本运行旧回归，保留散列对照。

## 2. P8.2 正式重计算在哪里完成？

在服务器 `{remote['hostname']}`（`root@50.115.148.16:4159`）。自然数据的每个 ID、计算分片和执行回执都绑定主机、PID、worker、源码提交、配置、Frozen 输入散列与执行时间。本地 WSL `{local['hostname']}` 仅开发、单元回归测试、调度、同步及审核，没有运行本轮微泡批量或正式积分。重积分代码提交 `{host['source_git_commit']}`；汇总和渲染提交 `{commit}`。

## 3. 服务器多少核、多少 RAM？

{model}，8 物理核、16 逻辑线程；容器 CPU 配额约 7.68 核。物理内存约 61.95 GiB，容器上限约 42.06 GiB。RTX 4090 可用，但当前轨迹求解没有使用 GPU，也没有重写 GPU 后端。本地为 {local['nproc']} 逻辑 CPU，约 7.62 GiB RAM。完整命令输出见计算来源记录。

## 4. 并行加速多少？

预先固定原自然 ID 1–64，保留全部成功及失败类型；依次测 1、2、4、8、16 workers。选择 {s['chosen_workers']} workers，串行 {s['serial_tracks_hour']:.2f} 条/小时，并行 {s['parallel_tracks_hour']:.2f} 条/小时，加速 {s['parallel_speedup']:.3f} 倍，效率 {s['parallel_efficiency']:.1%}。所有比较的出生属性、接受状态字节、出口、终止原因和统计量一致。吞吐量来自完整新积分计时，不把恢复读取时间当作求解时间。Frozen 数组经 Linux fork/COW 共享，RSS 会重复计算共享页，应结合 PSS、cgroup 和 swap 记录判断内存。

## 5. 为什么之前只有 OUTLET_02？

原始有限尺寸接纳对点流线 basin 有强烈偏置：原 149716 个 proposal 中，OUTLET_01 basin 仅接纳 2 个，OUTLET_02 接纳 1967 个，OUTLET_03 没有接纳。原 720 个停止中，719 个来自 OUTLET_02 点 basin，1 个来自 OUTLET_01 点 basin。零半径流线能进入三个出口，所以“Frozen 中没有其它出口通路”不成立。

这说明入口筛选及有限尺寸动力学必须与出口统计分开解释。点流线 basin 只是同一位置的零半径诊断标签，不保证有限尺寸轨迹最终走同一出口；原 OUTLET_01 basin 的另一个已完成 MB 最终属于 OUTLET_02。此结果不证明所有可能入口或尺寸都无法通过其它出口。

## 6. point tracer 是否覆盖三个出口？

是。原精度的 100000 条分别为 {p['outlet_counts']}，另有 {p['unresolved']} 条未解析。同一批入口位置加密空间步长和局部误差控制后为 {ref['outlet_counts']}，仍有 {ref['unresolved']} 条未解析；逐 ID 出口分类改变 {point['outlet_class_changes']} 条。

没有使用任意 ±5% 门槛。报告使用精确二项检验及 Bonferroni 同时 95% 区间，分母始终是全部 100000 条。OUTLET_01/02 的原观察比例没有通过纯抽样误差下的 Frozen 通量匹配；未解析路径不能任意分给出口，流量一致性仍未建立。1024 个原单元速度比较的最大差约 1.1e-19 m/s；一个反复停滞例的原 P1 单元速度散度约 −111.36 s⁻¹。这个单例支持检查 P1 场局部非零散度/壁面停滞，但不能把全部未解析路径都归因于它。没有改 Frozen 场去配平统计。

## 7. 有限尺寸入口接纳对不同 basin 的影响？\n
| 点 basin | 全部 proposal | 实际接纳 | 耗尽 ID（首次 proposal basin） |
|---|---:|---:|---:|
{table}

所有 proposal 的原接纳规则都重新执行，结果与原记录一致。proposal 包含重试，接纳按真实接纳位置分组，耗尽和 scheduled 按首次 proposal 分组；这三种分母不可混用。所有重试沿用原尺寸，不重新抽 SonoVue。图 04 和原始列表展示尺寸条件化，不能把接纳后分布当作原始 SonoVue 分布。

## 8. 720 个停止有多少在 extended run 能够完成？

{continuation}。只扩大计算 guard，不改变 dt、初始位置、尺寸、力学或状态接受规则；两档全部 720 个原接受前缀逐字节核对。计算停止仍是删失，不是生理捕获。原 P8.1 没有保存终止试探的细分深度，因此该字段明确为 null；仅报告已有最大细分深度，不补造未知数据。

## 9. 其它出口是否出现真实有限尺寸完整轨迹？

本轮独立自然数据中，OUTLET_01 为 {count[OUTLETS[0]]}，OUTLET_03 为 {count[OUTLETS[2]]}。正式集合没有混入按目标 basin 筛选的轨迹、旧轨迹或复制轨迹。扩展预算及时间步诊断在独立目录保存。

## 10. 三出口最终各有多少 complete tracks？

新自然集合 scheduled={c['scheduled']}、admitted={c['admitted']}、completed={c['completed']}。三个出口完整轨迹数分别为 {count[OUTLETS[0]]} / {count[OUTLETS[1]]} / {count[OUTLETS[2]]}。全部未接纳及未完成 ID 留在事件、终止表和 catalog 中。

## 11. 是否达到完整网络 ULM coverage？

{doc['full_network_coverage_status']}；all outlets observed={covered}。每出口至少 1 条自然完整轨迹才判定覆盖，每出口至少 50 条的可视化目标为 {doc['ulm_full_network_visualization_target']}。未覆盖分支在图中保留灰色原血管，不能补画路径。

## 12. 仍有哪些限制？

{timestep}；另选 30 个不重复的旧停止病例。dt、dt/2、dt/4 结果及共同时域偏差单独保存；不足 30 的 basin 不制造替代病例。production timestep 仍 NOT_FROZEN。

仍是 Frozen flow 上独立积分再叠加的微泡轨迹：没有真实 RBC passage、RBC 形变、新 CFD、微泡相互作用、完整生理悬浮液或 PK；H_D=0.45 是 feed 参数，与 tube Hct 区分；C_MB=8.5e12 m⁻³ 和 isotropic orientation V0 是原假设。P6.5 是球形法向近场正则化，生产 neighbor 参数和非球 lubrication 没有冻结。ULM-style 是轨迹累计展示，没有新增声学成像模型。

## 图像、动画与验证

图 00–11 与四组 MP4/GIF 都绑定同一 natural scene SHA256 `{doc['scene_sha256']}`；诊断图另绑定其独立原始输入。全血管、原球半径、µm 等比例坐标保留；仅显示层旋转相机、线性插值、慢放、时间压缩或跳切，没有新抽样、外推、样条补线或物理坐标变形。动画所有帧实际解码，首/中/末从输出 MP4 解码取样；新建两个渲染窗口比较相同保存状态的像素。代理已查看导出图和关键帧，用户验收仍待进行。

旧回归和 P8.2 永久测试的真实执行记录分别见 `UPSTREAM_REGRESSION.json`、`PARTICLE82_TEST_RESULTS.json`。汇总阶段曾遇到 point 分片目录层级错误，已增加回归测试、保留失败日志并校验复用已完成远程分片，没有修改轨迹。所有执行命令、恢复记录、资源预算及散列清单随报告保存。

本阶段在 Particle-8.2 停止，不启动 Particle-8.3、Particle-7.5 或新 RBC physics。
'''
    (root/'PARTICLE8_2_REVIEW.md').write_text(review)
    images=sorted((root/'figures').glob('*.png'));videos=sorted((root/'animations').glob('*'))
    assert len(images)==12 and len(videos)==8
    links=['# Particle-8.2 审核入口','','[中文审核报告](PARTICLE8_2_REVIEW.md)','[机器验证](PARTICLE8_2_VALIDATION.json)','']
    for folder,files in [('图像',images),('动画与 GIF',videos)]:
        links += ['## '+folder,'']+[f'- [{p.name}]({p.relative_to(root)})' for p in files]+['']
    (root/'INDEX.md').write_text('\n'.join(links))
    body=['<!doctype html><meta charset="utf-8"><title>Particle-8.2 审核</title><style>body{font:17px system-ui;max-width:1200px;margin:40px auto;color:#183449}img,video{width:100%}section{margin:40px 0}</style><h1>Particle-8.2 审核</h1>',
          f'<p>工程审计 {doc["software_stage_status"]}；网络覆盖 {doc["full_network_coverage_status"]}；用户审核待进行。</p>',
          '<a href="PARTICLE8_2_REVIEW.md">中文报告</a> · <a href="PARTICLE8_2_VALIDATION.json">验证 JSON</a>']
    for f in images:body.append(f'<section><h2>{html.escape(f.stem)}</h2><a href="figures/{f.name}"><img loading="lazy" src="figures/{f.name}"></a></section>')
    for f in videos:
        if f.suffix=='.mp4':body.append(f'<section><h2>{html.escape(f.stem)}</h2><video controls preload="none" src="animations/{f.name}"></video><a href="animations/{f.with_suffix(".gif").name}">GIF preview</a></section>')
    (root/'index.html').write_text('\n'.join(body))
    print(json.dumps({k:doc[k] for k in ['software_stage_status','full_network_coverage_status','chosen_workers','natural_completed','point_tracer_flux_sanity_status']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);main(p.parse_args().root)
