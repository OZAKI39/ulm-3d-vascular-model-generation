"""Render or open stored authorized repair results. Never starts a solver."""
import argparse
import html
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from py_scripts.fluid_physics.common import atomic_state,now,sha256_file
from .workflow import load_config,paths


def render():
    start=time.perf_counter();c=load_config();b,runs,out=paths(c)
    read=lambda p:json.loads(Path(p).read_text())
    auth=read(b/'authorization.json');assert auth['approved']
    attempts=read(runs/'gpu/budget_ledger.json')['attempts'];build_attempts=read(runs/'build/budget_ledger.json')['attempts']
    assert not any(x['status']=='RUNNING' for x in attempts+build_attempts)
    gpu=sum(x['charged_s'] for x in attempts);build_cost=sum(x['charged_s'] for x in build_attempts)
    failed=sum(x['charged_s'] for x in attempts if x['status']!='COMPLETED')
    records=[];analysis_costs=[]
    for a in attempts:
        task=a['task_id'];d=Path(a['directory']);execution=read(d/'execution.json');spec=read(d/'spec.json')
        folder=b/'runtime_reviews'/task
        review_file=folder/('review/review.json' if (folder/'review/review.json').exists() else 'review.json')
        review=read(review_file);assert review['execution']==execution
        from py_scripts.single_rbc_benchmark.analysis import read_csv
        metrics=read_csv(review_file.parent/'frame_metrics.csv')
        records.append(dict(task=task,directory=str(d),status=execution['status'],process_s=a['charged_s'],cached=False,
            dt=spec['dt'],planned_prep_steps=spec['prep_steps'],last_saved_step=review['last_saved_frame']['step'],
            last_saved_elapsed_time_star=review['last_saved_frame']['step']*spec['dt'],strain=0.,
            exact_failure_step_zero_based=6407 if task=='A5_continuous_preparation_30' else None,
            short_screen=review['short_screen'],sampled_area_max_relative=float(max(metrics['area_relative_drift'])),
            sampled_volume_max_relative=float(max(metrics['volume_relative_drift'])),review=str(review_file),
            implicit_reclassification_possible=task.startswith(('A0_','A1_','A2_','A3_')),
            completed_common_endpoint=False,completion=read(d/'completion.json') if (d/'completion.json').exists() else None))
        analysis_costs.append(dict(kind='saved_state_review',task=task,seconds=review['cpu_analysis_s']))
        if (folder/'conversion.json').exists():analysis_costs.append(dict(kind='HDF5_conversion',task=task,seconds=read(folder/'conversion.json')['cpu_conversion_s']))
    for name in ('A4_full_population_membership.json','A5_full_population_membership.json','A6_full_population_membership.json','A5_failure_step_geometry.json','A6_preparation_quality.json'):
        analysis_costs.append(dict(kind=name,seconds=read(b/name)['cpu_analysis_s']))
    member=read(b/'A6_full_population_membership.json');prep=read(b/'A6_preparation_quality.json')
    failure=read(b/'A5_failure_step_geometry.json');native=read(b/'native_trace_review.json')
    figures=read(b/'authorized_figures.json');analysis_costs.append(dict(kind='scientific_figures',seconds=figures['cpu_render_s']))
    costs=dict(authorized_s=dict(gpu=8500,cpu_solver=600,native_build=1800),charged_s=dict(gpu=gpu,cpu_solver=0,native_build=build_cost),
        remaining_s=dict(gpu=8500-gpu,cpu_solver=600,native_build=1800-build_cost),gpu_allocation_used_s=dict(A=gpu,B=0,C=0),
        gpu_allocation_remaining_s=dict(A=1600-gpu,B=1900,C=5000),failed_gpu_s=failed,successful_diagnostic_gpu_s=gpu-failed,
        gpu_allocations_s=auth['gpu_allocation_s'],no_automatic_retries=True,old_ledger_unchanged=True,
        known_separately_measured_cpu_analysis=analysis_costs,known_separately_measured_cpu_analysis_s=sum(x['seconds'] for x in analysis_costs),
        cost_semantics='GPU/build charges are monitored process wall seconds, including failed attempts. Diagnostic setup/output is included in that charge. Wall preparation is a setup subset and is not added twice. Listed CPU analyses are separate from solver execution. Untimed code-reading/native-table aggregation is not fabricated as zero. Browser validation has its own post-render receipt.')
    result=read(b/'comparison_results_before_authorized_execution.json')
    result.update(schema_version=3,recorded_at=now(),status='GATE_A_UNRESOLVED_HALF_DT_PREPARATION_COMPLETED_QUALITY_FAILED',
        root_cause_status='GEOMETRY_BEFORE_A5_OVERFLOW_CONFIRMED_INITIAL_TRIGGER_UNRESOLVED',runtime_fix_status='RUN_BOUNDARY_OBSERVATION_FIXED_SHORT_RANGE_STABILITY_ONLY',
        material_match='NOT_MATCHED',same_quality_same_endpoint_complete=False,qualified_speedup=None,
        new_solver_s=dict(gpu=gpu,cpu=0),new_compile_s=build_cost,budget=costs,new_raw_solver_records=records,
        new_cold_runs=dict(HemoCell=[],Mirheo=[]),new_end_to_end_s=dict(HemoCell=None,Mirheo=None),
        authorization=dict(record=str(b/'authorization.json'),approved=True,sha256=sha256_file(b/'authorization.json')),
        isolated_build=read(b/'isolated_build.json'),run_boundary_audit=str(b/'run_boundary_audit.json'),
        geometry_failure=failure,preparation_quality=prep,
        latest_complete_population_check=dict(record=str(b/'A6_full_population_membership.json'),stored_frames=len(member['frames']),
            total_confirmed_mismatched_point_frames=member['total_confirmed_mismatched_point_frames'],constant_species_ids=member['constant_species_ids'],
            uncertain_point_frames=sum(x['uncertain_count'] for x in member['frames']),scope=member['scope']),
        blockers=[
            'The ordered continuous dt=0.001 control still failed at zero-shear step 6407; its old pre-bounce membrane was already self-intersecting and beyond the WLC domain.',
            'The dt=0.0005 control completed only zero-shear t*=30, Gamma=0. Saved area drift reached 2.294581%, above the unchanged 2% gate; its final three residual windows failed the 0.002 criterion.',
            'Initial bending stress, independent membrane response, inertia, prepared-state comparability and the dt-dependent viscosity have no new matched acceptance evidence.',
            'No repaired configuration has passed the old shear-failure range or Gamma=4; no new matched cold repetitions or stricter shear window exists.'
        ],
        next_minimal_action='Resolve the remaining preparation failure before B/C: capture per-step membrane force and WLC extension before the blow-up and separate initial bending stress, thermal forcing and collision recoil in one-factor diagnostics. Then test at most two native material candidates against independent responses. Do not expand capacity to continue the already torn membrane, lower quality thresholds, or run the full timing queue.',
        browser_check='SEE_POST_RENDER_DELIVERY_RECEIPT',human_review='PENDING',git_commit_or_push=False)
    atomic_state(b/'comparison_results.json',result)
    for name in ('root_cause_evidence.json','failure_timeline.json','candidate_comparability.json','fix_log.md'):
        old=b/(Path(name).stem+'_before_authorized_controls'+Path(name).suffix)
        if not old.exists():shutil.copyfile(b/name,old)
    root=read(b/'root_cause_evidence_before_authorized_controls.json')
    root.update(recorded_at=now(),overall_status=result['root_cause_status'],runtime_fix_status=result['runtime_fix_status'],new_gpu_seconds=gpu)
    root['authorized_controls']=records
    root['hypotheses']+= [
        dict(id='repeated_run_observation_side_effect',status='CONFIRMED',finding_zh='每次 run() 都重新执行初始分类并清除对象力；correct_every=0 不阻止这种操作。A2 内侧数从 642 变为 643。连续观测修正已实际执行；旧 completion 中的 classifier_corrections=0 不是完整分类账目。',evidence=['run_boundary_audit.json']),
        dict(id='geometry_before_current_overflow',status='CONFIRMED',finding_zh='A5 第 6407 步进入前的旧膜已有 182 对独立确认的自交、最大 WLC 比 7.486；当前反弹前比值 88.425，单步位移 101.623。此例几何先于本次溢出损坏；最初触发与其准确起始步仍未捕获。',evidence=['A5_failure_step_geometry.json']),
        dict(id='half_dt_preparation_stability',status='SUPPORTED_NOT_CONFIRMED',finding_zh='顺序修正加连续调用的半步长对照完成零剪切 t*=30；这支持积分稳定性影响，不能单独证明所有旧故障根因或 Γ=4 稳定性，面积和残余准备质量仍失败。',evidence=['A6_preparation_quality.json','A6_full_population_membership.json'])]
    atomic_state(b/'root_cause_evidence.json',root)
    timeline=read(b/'failure_timeline_before_authorized_controls.json');timeline['authorized_control_runs']=records
    timeline['authorized_A5_exact_failure']=failure;timeline['scope']='Archived timeline retained, with separately identified authorized control evidence.'
    atomic_state(b/'failure_timeline.json',timeline)
    candidate=read(b/'candidate_comparability_before_authorized_controls.json')
    candidate.update(gate_A='UNRESOLVED_PREPARATION_QUALITY_FAILED_SHEAR_RANGE_NOT_TESTED',gate_B='NOT_MATCHED',gate_C='NOT_STARTED',
        runtime_variants=['ordered native + continuous dt=0.001: failed at preparation step 6407','ordered native + continuous dt=0.0005: t*=30 completed, quality failed'],
        new_material_candidates_fitted=0,preparation_quality=str(b/'A6_preparation_quality.json'),qualified_speedup=None)
    atomic_state(b/'candidate_comparability.json',candidate)
    residual=[x['residual_shape_rms_over_a_per_time'] for x in prep['transitions'][-3:]]
    table='\n'.join(f"| {x['task']} | {x['status']} | {x['last_saved_step']} | {x['last_saved_elapsed_time_star']:.4f} | {x['process_s']:.6f} |" for x in records)
    text=f'''# 单红细胞修复：实际运行与关口结果

**尚未完成同等质量比较，qualified_speedup=null。** 授权已使用；本次停止原因是物理质量关口失败，不是等待授权或显存不足。

原生库已隔离编译，未替换原安装。新增 local→halo 顺序、空 halo 早退和原生出错状态记录；另外发现并修正项目分段观测：Mirheo 每次 run() 都重新分类并清除对象力，因此改为每阶段一次连续调用、由原生插件保存中间状态。旧代码与错误的 classifier_corrections=0 原始字段保留，解释见 [run_boundary_audit.json](run_boundary_audit.json)。一个 CPU 行为回归先失败、修复后通过。

| 任务 | 原生执行状态 | 最后常规存帧步数 | 准备已计算 t*（存帧） | GPU 计费秒数 |
|---|---|---:|---:|---:|
{table}

以上全是 Γ=0 的诊断；常规存帧不是精确崩溃步。A0/A1 分别以 8659/7458 个候选超过 6400 失败。A2/A3 的短完成不能证明全程修复。A4 连续 5000 步完成，三个全粒子存帧未见可靠成员不符。

A5 在零基第 6407 步、t*=6.407 的 local outer 候选生成后报错，13,675 > 6,400；常规最后存帧为 6000 步。出错前旧位置已经有 **182 对自交、WLC 最大伸长比 7.486**；当前反弹前为 **68 对自交、最大伸长比 88.425**，顶点 189 单步移动 101.623。旧/新位置与 ID 同时保存，独立 CPU 法复核，未周期折回这些失效坐标。[出错一步的原始证据与几何](A5_failure_step_geometry.json)。因此本次不是已证实可安全扩容的正常碰撞负载；最初失稳发生在更早的未保存步骤，尚不能断言是弯曲应力、热噪声或反冲中的哪一个首先触发。

A6 把 dt 减半为 0.0005，同步保持物理准备时长，完成 **60,000 步、t*=30、Γ=0**，没有原生溢出。60 个完整液体存帧（t*=0 到 29.5）中，110,592 个粒子及内侧 642 个 ID 保持不变；确认成员不符 {member['total_confirmed_mismatched_point_frames']} 个点帧，近膜不确定 {sum(x['uncertain_count'] for x in member['frames'])} 个点帧。未保存的时间、最终 t*=30 的全液体状态不在这项全量结论内，终点另有 192 个探针。[全粒子核查](A6_full_population_membership.json)

但 A6 在 t*=28.5 的面积偏差 **2.294581% > 2%**；最后三组相邻准备窗口的形状残余率为 **{residual[0]:.6f}、{residual[1]:.6f}、{residual[2]:.6f} > 0.002**。末帧相对共享参考形状的中心去除 RMS/a 为 {prep['final_reference_shape_rms_over_a']:.6f}，这是参考形状诊断，不冒充双方新配对误差。延长至原定最大准备时间仍没有通过质量要求。[准备判据及逐窗口结果](A6_preparation_quality.json)

两边细胞尚未证明足够相同。本轮未开始新的材料拟合；原始伸长响应差 17.82%、常参考角造成的初始应力差、约 0.989 的膜/排液质量比和半步黏度差 6.47% 仍待处理。半步稳定性不能替代新的黏度和膜响应匹配。材料训练/留出标准仍见 [material_matching.json](material_matching.json)，没有改阈值，也没有声称已穷尽两组候选。

两边本轮到 Γ=4 的合格冷启动次数均为 0，端到端耗时均为空。历史 HemoCell 的 58.171937、60.782516 秒仅保留为旧配置结果，不能与本轮准备诊断相除。A6 的 {records[-1]['process_s']:.6f} 秒仅是这项零剪切诊断的完整进程费用。

本轮实际 GPU **{gpu:.6f}/8500 秒**（其中失败 {failed:.6f} 秒），CPU 求解 **0/600 秒**，隔离编译 **{build_cost:.6f}/1800 秒**。剩余 GPU {8500-gpu:.6f} 秒、CPU 求解 600 秒、编译 {1800-build_cost:.6f} 秒；A/B/C 的剩余额度分别为 {1600-gpu:.6f}/1900/5000 秒，未转移额度、未重置历史账本、未自动重试。单独 CPU 分析、图表与浏览器费用另列，编译控制器时间和墙粒子准备子项没有重复相加。

下一项最小行动是定位 WLC 越界之前的每步膜力与双侧反冲，分离初始弯曲应力、热扰动和积分响应，再决定最多两组原生材料候选。当前 A 的质量关口未通过，B 标定和 C 的 Γ=4 计时按原停止规则不启动。禁止靠容量扩展、删粒子、强制重分类或放宽阈值取得速度比。

[当前结果 JSON](comparison_results.json) · [根因状态](root_cause_evidence.json) · [候选准入](candidate_comparability.json) · [实际 GPU 账本](../../../runs/single_rbc_repair/{c['campaign_id']}/gpu/budget_ledger.json) · [离线中文 HTML](../../../test_code/outputs/single_rbc_repair/{c['campaign_id']}/comparison_review.html)

查看现有页面（不启动求解，不重新渲染）：

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.single_rbc_repair.authorized_result_report --open
```

CPU 预检命令仍是 `.venv/bin/python -B -m py_scripts.repair_single_rbc_benchmark --config py_scripts/single_rbc_benchmark_repaired.yaml --preflight-only`。原冻结 --execute 队列已有失败记录，不能作为自动重试入口。新页面渲染命令为 `.venv/bin/python -B -m py_scripts.single_rbc_repair.authorized_result_report --render`，仅读取已经保存的结果。浏览器验收绑定当前 HTML 哈希，单独见 delivery_receipt.json；人工验收 PENDING。本轮未提交或推送。
'''
    (b/'report_zh.md').write_text(text)
    (b/'fix_log.md').write_text('# 授权后的修复与复验\n\n'+text.split('| 任务')[0].split('\n',2)[-1]+'\n\n实际控制、首次复发、半步结果与质量失败详见 [report_zh.md](report_zh.md)。原生补丁为 [local_before_halo_v2.patch](native/local_before_halo_v2.patch)，未改容量或反弹公式。新 Python 连续观测代码位于 py_scripts/single_rbc_repair/continuous_worker.py 与 continuous_protocol.py，各次实际源码、spec 和执行记录保存在对应 runs 目录。\n\n[授权前修复记录原文](fix_log_before_authorized_controls.md) 保留，不将短测通过称为完整修复。\n')
    def link(name,label):return '<a href="'+html.escape(os.path.relpath(b/name,out))+'">'+html.escape(label)+'</a>'
    links=' · '.join(link(n,l) for n,l in [('report_zh.md','完整中文报告'),('comparison_results.json','当前结果与费用'),('root_cause_evidence.json','根因状态'),('A5_failure_step_geometry.json','出错一步原始几何'),('A6_full_population_membership.json','全粒子成员核查'),('A6_preparation_quality.json','准备窗口判据'),('isolated_build.json','隔离库与哈希'),('delivery_receipt.json','独立浏览器验收')])
    tr=''.join('<tr><td>'+html.escape(x['task'])+'</td><td>'+('完成准备诊断' if x['status']=='COMPLETED' else '原生失败')+f"</td><td>{x['last_saved_step']}</td><td>{x['process_s']:.3f}</td></tr>" for x in records)
    svg=[]
    for filename in ('stored_preparation.svg','failure_geometry.svg'):
        source=(out/'authorized_figures'/filename).read_text();svg.append(source[source.index('<svg'):])
    panel=f'''<section id="repair-current" style="padding:24px;background:#f3f7f8;border:3px solid #28627a;margin:16px">
<style>#repair-current table{{border-collapse:collapse;width:100%}}#repair-current td,#repair-current th{{border:1px solid #bccbd1;padding:7px;text-align:left}}#repair-current .evidence-figure>svg{{max-width:100%;height:auto}}</style>
<h1>单红细胞修复：连续半步准备已完成，质量关口未通过</h1>
<p><strong>尚未完成同等质量比较；qualified_speedup = null。</strong>隔离构建与 7 个 GPU 对照已实际执行。最远完成零剪切 60,000 步、t*=30、Γ=0；没有新的 Γ=4 冷启动结果。</p>
<p>已修正分段 run() 触发重复分类与对象力清除的观测方式。顺序与空 halo 补丁仍不足以消除原步长故障：A5 第 6407 步有 13,675 个候选超过 6,400，进入该步的旧膜已经自交 182 对、WLC 比 7.486；本次几何先于候选溢出损坏，最初触发仍待确定。</p>
<p>A6 半步长完成 t*=30，但面积最大偏差 <strong>2.2946% &gt; 2%</strong>；最后三组形状残余率为 {residual[0]:.5f}、{residual[1]:.5f}、{residual[2]:.5f}，均高于 0.002。<strong>关口 A 质量失败，B 材料未匹配，C 未启动。</strong>没有改变通过门槛。</p>
<p>60 个全液体存帧未见可靠成员不符，17 个近膜点帧保留为不确定；不把存帧结果当作连续时间的不穿透证明。双方初态、膜本构、惯性和新 dt 下黏度仍缺合格匹配。</p>
<table><thead><tr><th>实际任务</th><th>执行</th><th>最后常规存帧步数</th><th>GPU 秒</th></tr></thead><tbody>{tr}</tbody></table>
<p>GPU 已用 {gpu:.6f}/8500 秒（失败 {failed:.6f} 秒）；CPU 求解 0/600 秒；隔离编译 {build_cost:.6f}/1800 秒。剩余 {8500-gpu:.6f}、600、{1800-build_cost:.6f} 秒。费用含失败，未重试或挪用旧余额。停止原因是科学判据失败。</p>
<p>双方本轮完整端到端耗时为空。历史 HemoCell 的 58.17/60.78 秒属于旧配置；下方历史动画不代表新库完成剪切。原环境与 20,978 个旧文件哈希核对通过。人工验收 PENDING。{links}</p>
<h2>本轮实际准备存帧</h2><p>所有横轴都是零剪切准备时间，Γ 始终为 0。曲线只连实际存帧；失败后不延长。</p><div class="evidence-figure">{svg[0]}</div>
<h2>A5 原生出错一步：旧位置与当前反弹前位置</h2><p>按真实顶点 ID 和三角连接绘制，未折回已经损坏的膜；左右坐标范围不同并明确标出。</p><div class="evidence-figure">{svg[1]}</div>
</section>'''
    page=(out/'comparison_review_before_authorized_execution.html').read_text();begin=page.index('<section id="repair-current"');end=page.index('<section id="repair-legacy"',begin)
    page=page[:begin]+panel+page[end:];encoded=json.dumps(result,ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    page=re.sub(r'(<script id="repair-results" type="application/json">).*?(</script>)',lambda m:m[1]+encoded+m[2],page,count=1,flags=re.S)
    page=page.replace('修复 CPU 核查与历史记录 · ','修复实际运行核查 · ',1);(out/'comparison_review.html').write_text(page)
    record=dict(recorded_at=now(),html=str(out/'comparison_review.html'),html_sha256=sha256_file(out/'comparison_review.html'),
        comparison_results_sha256=sha256_file(b/'comparison_results.json'),renderer_sha256=sha256_file(__file__),render_s=time.perf_counter()-start,
        solver_started=False,compiler_started=False,human_review='PENDING',browser_check='SEE_POST_RENDER_DELIVERY_RECEIPT')
    atomic_state(b/'review_latest.json',record);atomic_state(b/'authorized_result_render.json',record)
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__);m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--render',action='store_true',help='Rebuild current report from completed stored analyses only')
    m.add_argument('--open',action='store_true',help='Open the existing HTML; no rendering or solver')
    m.add_argument('--verify',action='store_true',help='Check stored HTML and result hashes; no rendering or solver')
    a=p.parse_args();c=load_config();b,_,out=paths(c)
    if a.render:print(json.dumps(render(),ensure_ascii=False,indent=2));return
    r=json.loads((b/'review_latest.json').read_text());assert sha256_file(r['html'])==r['html_sha256']
    assert sha256_file(b/'comparison_results.json')==r['comparison_results_sha256']
    if a.open:
        target=subprocess.check_output(['/usr/bin/wslpath','-w',r['html']],text=True).strip()
        subprocess.Popen(['/mnt/c/Windows/explorer.exe',target],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    print(json.dumps(dict(status='STORED_HASHES_PASS',html=r['html'],solver_started=False,compiler_started=False),ensure_ascii=False))


if __name__=='__main__':main()
