"""Render authorized build/resource status from stored evidence; no solver call."""
import html
import json
import os
from pathlib import Path
import re
import shutil
import time

from py_scripts.fluid_physics.common import atomic_state,now,sha256_file,write_json
from .workflow import load_config,paths


def render():
    start=time.perf_counter();c=load_config();b,runs,out=paths(c)
    read=lambda p:json.loads(Path(p).read_text())
    auth=read(b/'authorization.json');assert auth['approved'] is True
    build=read(b/'isolated_build.json');assert build['status']=='BUILT'
    assert not list((runs/'gpu').rglob('execution.json')) and not list((runs/'cpu').rglob('execution.json'))
    resource=read(b/'gpu_resource_after_build.json');assert resource['status']=='GPU_RESOURCE_BLOCKED_NOT_LAUNCHED'
    charged=sum(a['charged_s'] for a in read(runs/'build/budget_ledger.json')['attempts'])
    budget=dict(authorized_s=dict(gpu=8500,cpu_solver=600,native_build=1800),charged_s=dict(gpu=0,cpu_solver=0,native_build=charged),remaining_s=dict(gpu=8500,cpu_solver=600,native_build=1800-charged),gpu_allocations_s=auth['gpu_allocation_s'],old_ledger_unchanged=True,no_automatic_retries=True)
    result=read(b/'comparison_results_before_authorized_execution.json')
    result.update(recorded_at=now(),status='AUTHORIZED_NATIVE_BUILT_GPU_RESOURCE_BLOCKED',new_compile_s=charged,
                  authorization=dict(approved=True,record=str(b/'authorization.json'),sha256=sha256_file(b/'authorization.json')),
                  budget=budget,isolated_build=build,gpu_resource=resource,
                  cpu_runtime_reader_tests=read(b/'cpu_runtime_reader_tests.json'),
                  blockers=['GPU free memory below frozen 3336 MiB launch guard; no GPU attempt was launched or charged',
                            'Native scheduling candidate compiled but causal short controls and Gamma=4 range have not run',
                            'Material/preparation/viscosity matching and complete repeated comparisons remain unverified'],
                  browser_check='NOT_RUN_ON_THIS_RENDER',git_commit_or_push=False)
    atomic_state(b/'comparison_results.json',result)
    used=resource['resources']['gpus'][0]['free_MiB'];need=resource['minimum_additional_free_MiB']
    status=f'''# 授权后的构建与资源状态

GPU 8500 秒、CPU 求解 600 秒、隔离编译 1800 秒已由用户明确批准，记录在 [authorization.json](authorization.json)。申请原件中的 approved=false 是申请时状态，原件未改写。

**隔离 Mirheo 库已编译成功，但尚未完成同等质量比较，qualified_speedup=null。** 新 GPU 和 CPU 求解均未启动；首批四个各最多 90 秒的对照仍待执行。当前显存采样为 {used} MiB 可用，冻结门槛为 3336 MiB（1536 MiB 保留 + 1800 MiB 估计任务用量），至少还需 {need} MiB。没有降低门槛或结束用户的其他程序。

编译实际计费 {charged:.9f} 秒，剩余 {1800-charged:.9f} 秒；GPU 剩余 8500 秒，CPU 求解剩余 600 秒。控制器完整耗时另见 [构建控制记录](isolated_build_controller_result.json)，不与已计费子任务相加。两次构建子任务均完成，没有重试。新增 4 项 CPU 存帧读取回归通过；不计为原生求解或 GPU 验证。

新库 SHA-256：`{build['library_sha256']}`。新库位置、准确命令和费用见 [isolated_build.json](isolated_build.json)、[构建账本](../../../runs/single_rbc_repair/{c['campaign_id']}/build/budget_ledger.json)、[编译原始日志](../../../runs/single_rbc_repair/{c['campaign_id']}/build/compile_isolated/console.log)。旧安装库 SHA-256 仍为 `{build['installed_library_sha256']}`，20,978 个受保护旧文件核对通过。

关口 A 的运行因果尚未验证，B 材料仍未匹配，C 没有启动。本轮没有新的共同终点、冷启动耗时、膜轨迹或速度比。旧四次 Mirheo 候选溢出与 CPU 几何证据保持原结论；编译成功不能证明调度假设解释了这些失败。

[实际资源阻塞记录](gpu_resource_after_build.json) · [当前结果 JSON](comparison_results.json) · [存帧读取 CPU 测试](cpu_runtime_reader_tests.json) · [离线页面](../../../test_code/outputs/single_rbc_repair/{c['campaign_id']}/comparison_review.html)

显存满足条件后，已有授权可继续首批冻结队列，无需新增预算：

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.repair_single_rbc_benchmark --config py_scripts/single_rbc_benchmark_repaired.yaml --execute
```

人工验收保持 PENDING；新的浏览器结果另存于当前输出目录的 browser_after_build。以下是授权前 CPU 报告原文，其中“等待授权”“未编译”等措辞仅描述当时状态。

---

'''
    atomic_report=b/'report_zh.md';atomic_report.write_text(status+(b/'report_zh_before_authorized_execution.md').read_text())
    source=out/'comparison_review_before_authorized_execution.html'
    if not source.exists():shutil.copyfile(out/'comparison_review.html',source)
    page=source.read_text()
    def link(name,label):return '<a href="'+html.escape(os.path.relpath(b/name,out))+'">'+html.escape(label)+'</a>'
    links=' · '.join(link(n,t) for n,t in [('authorization.json','用户授权'),('isolated_build.json','新库与构建用量'),('gpu_resource_after_build.json','显存门槛与实测'),('comparison_results.json','当前结果'),('failure_timeline.json','历史失败时间线'),('root_cause_evidence.json','根因证据'),('material_matching.json','材料未匹配项'),('cpu_runtime_reader_tests.json','CPU 读取回归'),('report_zh.md','中文报告')])
    panel=f'''<section id="repair-current" style="padding:24px;background:#eef6fa;border:3px solid #28627a;margin:16px;max-width:none">
<h1>单红细胞剪切流修复：隔离库已编译，GPU 等待显存</h1>
<p><strong>尚未完成同等质量比较；qualified_speedup = null。</strong>GPU、CPU 求解和隔离编译授权均已记录。新的隔离 Mirheo 库编译成功，原安装库没有替换；四个短对照尚未启动。</p>
<p>当前可用显存 {used} MiB，冻结启动门槛 3336 MiB，还需至少 {need} MiB。门槛包含 1536 MiB 保留和 1800 MiB 任务估计。没有 GPU 求解尝试或 GPU 计费，授权与资源就绪是两个独立条件。</p>
<table><thead><tr><th>预算类别</th><th>已授权秒数</th><th>本轮实际计费秒数</th><th>剩余秒数</th></tr></thead><tbody><tr><td>GPU 求解</td><td>8500</td><td>0</td><td>8500</td></tr><tr><td>CPU 求解</td><td>600</td><td>0</td><td>600</td></tr><tr><td>隔离编译</td><td>1800</td><td>{charged:.6f}</td><td>{1800-charged:.6f}</td></tr></tbody></table>
<p>两个构建子任务均已完成，无重试。旧预算和旧账本保持原样。新增 4 项 CPU 存帧读取回归通过，检查 ID、相位和实际完成步数；这些测试不代表 GPU 稳定性验证。</p>
<h2>物理结论保持证据边界</h2><p>local/halo 调度候选已编译，但实际故障因果仍为 SUPPORTED_NOT_CONFIRMED。旧成员异常、自交和弯曲初应力差继续保留；材料 NOT_MATCHED，关口 A 尚待运行验证，B 未通过，C 未启动。</p>
<p>两边新冷启动完成次数均为 0，新端到端耗时均为空。下方真实动画和计时全部属于旧 campaign，不代表新库已经运行到 Γ=4。</p>
<p>人工验收 PENDING。{links}</p></section>
'''
    begin=page.index('<section id="repair-current"');end=page.index('<section id="repair-legacy"',begin)
    page=page[:begin]+panel+page[end:]
    encoded=json.dumps(result,ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    page=re.sub(r'(<script id="repair-results" type="application/json">).*?(</script>)',lambda m:m[1]+encoded+m[2],page,count=1,flags=re.S)
    page=page.replace('修复 CPU 核查与历史记录 · ','修复构建与资源核查 · ',1)
    (out/'comparison_review.html').write_text(page)
    # Retain the historical link contained in the reused legacy page.
    companion=out/'browser_delivery_verified/browser_check.json';companion.parent.mkdir(exist_ok=True)
    old=b.parents[2]/'test_code/outputs/single_rbc_benchmark/rbc_shear_20260910/browser_delivery_verified/browser_check.json'
    if not companion.exists():shutil.copyfile(old,companion)
    assert sha256_file(companion)==sha256_file(old)
    record=dict(recorded_at=now(),html=str(out/'comparison_review.html'),html_sha256=sha256_file(out/'comparison_review.html'),previous_html=str(source),previous_html_sha256=sha256_file(source),render_s=time.perf_counter()-start,solver_started=False,human_review='PENDING',browser_check='NOT_RUN_ON_THIS_RENDER',historical_animation=True,authorization_approved=True)
    atomic_state(b/'review_latest.json',record)
    write_json(b/'authorized_build_status_render.json',record)
    return record


if __name__=='__main__':print(json.dumps(render(),ensure_ascii=False,indent=2))
