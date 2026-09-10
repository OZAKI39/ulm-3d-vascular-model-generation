"""Reuse the archived offline Plotly page, explicitly labelled as history."""
import html
import json
import os
import re
from pathlib import Path
import time
from py_scripts.fluid_physics.common import atomic_state,sha256_file
from .workflow import paths
from .forensics import ROOT,read


def export(c):
    started=time.perf_counter();b,_,out=paths(c);out.mkdir(parents=True,exist_ok=True)
    result=read(b/'comparison_results.json');timeline=read(b/'failure_timeline.json');causes=read(b/'root_cause_evidence.json')
    def link(name):
        return '<a href="'+html.escape(os.path.relpath(b/name,out))+'">'+html.escape(name)+'</a>'
    rows=[]
    for r in timeline['runs']:
        p,x=r['first_confirmed_membership_discrepancy'],r['first_confirmed_nonadjacent_intersection']
        a='未观察到' if p is None else ('准备阶段第 '+str(p['global_step'])+' 步' if p['phase']=='relaxation' else 'Γ='+str(p['strain']))
        z='存帧未检出' if x is None else 'Γ='+str(x['strain'])+'，'+str(x['confirmed_count'])+' 对'
        rows.append('<tr>'+''.join('<td>'+html.escape(str(s))+'</td>' for s in [r['task'],a,z,r['last_saved_frame']['strain'],str(r['error']['failing_step_interval'])])+'</tr>')
    cause_rows=''.join('<tr><td>'+html.escape(x['id'])+'</td><td>'+html.escape(x['status'])+'</td><td>'+html.escape(x['finding_zh'])+'</td></tr>' for x in causes['hypotheses'])
    panel='''<section id="repair-current" style="padding:24px;background:#eef6fa;border:3px solid #28627a;margin:16px;max-width:none">
<h1>单红细胞剪切流修复：本轮 CPU 核查</h1>
<p><strong>尚未完成同等质量比较；qualified_speedup = null。</strong>本轮新增求解用量为 0，隔离补丁尚未编译。旧授权没有被当作新授权。</p>
<p>已确认默认配置在零剪切准备阶段发生远离膜面的成员不符；受限队列诊断 Γ=2.2 出现 21 对可靠自交和 WLC 超伸长。local/halo 缓冲复用存在缺少执行依赖的源码证据，尚待实际短对照验证因果。</p>
<p>初始节点力 RMS：HemoCell 约 3.9×10⁻¹⁰，Mirheo 约 1.77。独立 Kantor 能量梯度与 Mirheo 原生初始总力的相对差约 0.14%，支持弯曲参考应力造成初态失配；这不是运行崩溃原因的证明。</p>
<table><thead><tr><th>历史任务</th><th>首次可靠成员不符</th><th>首次可靠自交</th><th>最后存帧 Γ</th><th>失败步区间（全局，1 起算）</th></tr></thead><tbody>'''+''.join(rows)+'''</tbody></table>
<p>区间来自保存下界和下一输出边界；精确失败步没有旧记录。近膜带为 10⁻⁴，自交容差为 10⁻⁵，均为公共长度单位。探针统计是点帧次数；不作全部粒子泄漏率。</p>
<table><thead><tr><th>原因</th><th>证据等级</th><th>当前结论</th></tr></thead><tbody>'''+cause_rows+'''</tbody></table>
<h2>修复与可比性关口</h2><p>已准备 local→halo 顺序依赖、空 halo 早退出、共享/独立 bouncer 对照、受控原生状态记录和独立材料探针。原库、旧案例、旧数据均保留。A 未通过，B 未通过，C 未启动；两边新冷启动完成次数均为 0。</p>
<p>申请上限：GPU 求解 8500 秒（A 1600，B 1900，C 5000），CPU 求解 600 秒，隔离编译 1800 秒。失败计费，单 GPU、无自动重试；前关失败即停止后续昂贵实验。正式参数尚未冻结；只有规则与第一批短对照已冻结。</p>
<p>旧 HemoCell 冷任务为 58.17、60.78 秒；旧 Mirheo 仅有失败成本。本轮新完整耗时均为空。CPU 重分析、编译、诊断、匹配和渲染成本分别记录。</p>
<p>人工验收：PENDING。'''+ ' · '.join(link(n) for n in ['failure_timeline.json','root_cause_evidence.json','initial_bending_force_audit.json','material_matching.json','candidate_comparability.json','frozen_benchmark_plan.json','authorization_request.json','fix_log.md','report_zh.md'])+'''</p>
</section><section id="repair-legacy" style="border-top:8px solid #b78b27;padding-top:20px">
<div style="padding:20px;background:#fff1cb"><h2>以下整页为旧 campaign 的原始离线结果</h2><p>下方动画、D、倾角、A/V、计时与“本轮”等原文均属于 rbc_shear_20260910 历史记录，不代表隔离修复后的新结果。复用其真实存帧与既有交互代码，缺失 Γ=4 仍保持空白。</p></div>
'''
    legacy=ROOT/'test_code/outputs/single_rbc_benchmark/rbc_shear_20260910/single_rbc_review.html'
    page=legacy.read_text()
    body=re.search(r'<body\b[^>]*>',page)
    main=re.search(r'<main\b[^>]*>',page)
    if body is None and main is None:raise ValueError('LEGACY_TEMPLATE_BODY_NOT_FOUND')
    encoded=json.dumps(result,ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    insertion=body.end() if body else main.start()
    page=page[:insertion]+panel+page[insertion:]
    ending='</body>' if '</body>' in page else '</html>'
    page=page.replace(ending,'</section><script id="repair-results" type="application/json">'+encoded+'</script>'+ending)
    page=page.replace('<title>','<title>修复 CPU 核查与历史记录 · ',1)
    target=out/'comparison_review.html';target.write_text(page,encoding='utf8')
    record=dict(html=str(target),html_sha256=sha256_file(target),legacy_html_sha256=sha256_file(legacy),
                render_s=time.perf_counter()-started,solver_started=False,human_review='PENDING',
                browser_check='NOT_RUN_ON_THIS_RENDER',historical_animation=True)
    atomic_state(b/'review_latest.json',record)
    return record
