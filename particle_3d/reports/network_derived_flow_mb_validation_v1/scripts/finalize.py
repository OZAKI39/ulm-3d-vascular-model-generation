"""Assemble Chinese review and immutable-input evidence from measured artifacts."""
from pathlib import Path
import datetime,hashlib,json,html,xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
R=Path(__file__).resolve().parents[1]
read=lambda p:json.loads(Path(p).read_text())
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def table(head,rows):return '| '+' | '.join(head)+' |\n| '+' | '.join(['---']*len(head))+' |\n'+''.join('| '+' | '.join(map(str,row))+' |\n' for row in rows)
def ms(x):return '—' if x is None else f'{x*1e3:.3f}'

def main():
    s=read(R/'data/analysis_summary.json');ctx=read(R/'run_context.json');cohort=read(R/'data/cohort_provenance.json')
    protection=read(R/'data/particle_science_protection_final.json');remote=read(R/'logs/remote_final_snapshot_verification.json')
    det=read(R/'data/worker_count_determinism_subset.json');trans=read(R/'data/paired_outcome_transition.json')
    old=read(R/'data/old_paired_metrics.json');new=read(R/'data/new_paired_metrics.json');point=read(R/'data/point_vs_mb.json');audit=read(R/'data/saved_state_safety_audit.json')
    suites=list(ET.parse(R/'logs/tests.xml').getroot().iter('testsuite'))
    tests={k:sum(int(t.get(k,'0')) for t in suites) for k in ['tests','failures','errors','skipped']};tests['elapsed_s']=sum(float(t.get('time','0')) for t in suites)
    assert tests['tests']>=12 and tests['failures']==tests['errors']==tests['skipped']==0
    assert det['status']=='PASS' and protection['science_files_changed']==remote['source_files_changed']==0
    artifacts=[]
    for p in sorted((R/'figures').glob('*.png')):
        with Image.open(p) as im:
            im.verify()
        with Image.open(p) as im:width,height=im.size
        pdf=p.with_suffix('.pdf');assert pdf.read_bytes().startswith(b'%PDF') and pdf.stat().st_size>1000
        artifacts.append(dict(png=str(p),pdf=str(pdf),width=width,height=height,png_sha256=sha(p),pdf_sha256=sha(pdf)))
    assert len(artifacts)==7
    receipt=[]
    for label in ['OLD','NEW']:
        for p in sorted((R/'outputs'/label/'trajectories').glob('*.json')):
            m=read(p);npz=p.with_suffix('.npz');assert sha(npz)==m['samples_sha256']
            receipt.append(dict(label=label,bubble_id=m['particle_id'],meta_sha256=sha(p),samples_sha256=sha(npz)))
    assert len(receipt)==60
    for p in (R/'outputs/NEW_POINT').glob('point_*.json'):
        assert sha(p.with_suffix('.npz'))==read(p)['path_sha256']
    dump(R/'data/trajectory_artifact_verification.json',dict(status='PASS',paired_trajectory_count=60,point_trajectory_count=30,rows=receipt))
    dump(R/'data/figure_artifact_checks.json',dict(status='PASS',figures=artifacts,visual_review='Figures 01, 02, 04 directly inspected; 04 label spacing corrected without data changes',dpi=320))
    runtimes={name:read(R/'outputs'/name/'completed.json')['wall_seconds'] for name in ['determinism_w1','determinism_w3','NEW','NEW_POINT']}
    history=R.parents[1]/'outputs/particle9a1_2mmps/provenance/smoke30_completed.json'
    runtimes['historical_old_smoke_batch_s']=read(history)['wall_seconds']
    now=datetime.datetime.now(datetime.timezone.utc);began=datetime.datetime.strptime(ctx['stamp'],'%Y%m%dT%H%M%SZ').replace(tzinfo=datetime.timezone.utc)
    s.update(science_files_changed=0,protected_science_files_checked=protection['files_checked'],remote_snapshot_files_checked=remote['source_files_checked'],
        tests=tests,worker_count_determinism=det,worker_6_parity='PASS: permanent test compares all samples and full decompressed logs for IDs 3/15/22 in 1, 3 and 6 workers',
        runtime_seconds=runtimes,task_elapsed_seconds_to_report=(now-began).total_seconds(),completed_utc=now.isoformat(),
        remote_workspace=ctx['remote'],report_directory=str(R),figures=artifacts,
        figure02_path=str(R/'figures/02_paired_old_new_trajectories.png'),figure03_path=str(R/'figures/03_paired_outcome_transition.png'),
        main_report_path=str(R/'NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md'),
        nan_inf_count_point=0,point_count=30,point_outlet_counts=dict(__import__('collections').Counter(r['point_outcome'] for r in point)),
        final_safety_gate_pass=all(s[k]==0 for k in ['penetration_count_new','handoff_violation_count_new','continuous_certificate_violations_new','solver_fail_count_new','nan_inf_count_new','accepted_diagnostic_nan_inf_count_new','unclassified_solver_corruption_new']),
        production_500_launched=False,old_results_deleted_or_overwritten=False,taylor_hood_active_processes_at_stop=0)
    assert s['final_safety_gate_pass']
    dump(R/'data/final_summary.json',s)
    counts=table(['N=30 固定配对组','completed','stationary','solver fail','O1','O2','O3'],[[label,s[f'{k}_completed'],s[f'{k}_stationary'],s[f'{k}_fail'],*[s[f'{k}_outlet_counts'][o] for o in ['O1','O2','O3']]] for label,k in [('OLD','old'),('NEW','new')]])
    transitions=table(['OLD → NEW','颗数'],[[k,v] for k,v in s['transition_matrix'].items()])
    pointtable=table(['NEW point tracer → NEW MB','颗数'],[[k,v] for k,v in s['point_vs_MB_transition'].items()])
    stopped='Taylor-Hood 已安全停止并保留全部已有进度。原 ILU2 超过 25 分钟且主存约 36 GiB，RCM 后 ILU2 出现 GPU OOM；停止前最后一次单步 ILU1 smoke 已自行正常退出，但仍没有完整稳态 P2 解。因此按资源工程成本停止，不构成 Taylor-Hood 科学失败。详见 [停止说明](TAYLOR_HOOD_STOP_NOTE_ZH.md)。'
    main=f'''# 一句话结论

仅替换为 network-derived 背景流后，这 30 个固定 Bubble 得到 **29 completed、1 stationary、0 solver fail**；与 OLD 配对相比，**19 颗从 O2 改走 O3**。穿透、未解决 handoff 违规、NaN/Inf、未分类求解器损坏均为 0。最终状态：`{s['final_status']}`。

{counts}
# 为什么停止 Taylor-Hood

{stopped}

状态固定为 `TAYLOR_HOOD_VALIDATION_STOPPED_RESOURCE_COST`。本机 37,689 个、服务器 2,685 个先前受保护文件检查均为 0 变更；没有删除失败日志，也没有继续调试或运行完整 P2 求解。

# 新旧流场有什么区别

OLD 使用原 equal-pressure 2.0 mm/s 场，NEW 使用网络推导的出口压力场。NEW authoritative VTU 由 case 的 `reports/physics_validation_H0.json` 中 PASS export 记录定位并核对 SHA，未按名字猜测。

重新从实际 VTU 对原始边界三角形进行 P1 精确积分：OLD 的 O1/O2/O3 = **4.257917/85.205026/10.537057%**；NEW = **6.540825/44.820122/48.639054%**。两者 Qin = 1.5513591604402322e-14 m³/s，原始双精度求和相对总质量残差约 2.03e-16，均通过。独立已存审计使用补偿求和得到另一舍入量级的残差；本轮不将残差强制置零。

OLD SHA256：`{s['old_flow_sha']}`。NEW SHA256：`{s['new_flow_sha']}`。

两者均为 70,363 点、371,402 个 TET4。坐标、连接、WALL/INLET/O1/O2/O3 及 Global IDs 的 numeric equality 和 bitwise equality 均通过。[几何契约](data/old_new_flow_contract.json)及[重新积分结果](data/reintegrated_flow_flux.json)保存全部路径和哈希。

# 这一轮哪些东西没有改变

同一 ID 的直径、物理半径、初始位置、四元数、保存的初始速度与角速度、seed/stream、初始 tetra 全部相同。原 P9-A.1/P6.5、nearfield、lubrication、壁面流体动力学、translation–rotation、contact、handoff、active-set、积分、终止与出口分类均保持原实现。dt=0.00025 s、物理 horizon=1.5 s、所有原有安全终止门槛不变。SonoVue、浓度、半径和 clearance/floor 未改。

本轮新增的 adapter 仅负责选择已冻结输入和装入已保存的初始状态，调用原 `integrate_admitted`、`Particle9AStepper` 与原可选诊断观察器。NEW 后续每一步读取 NEW 的速度/梯度。没有速度平滑、divergence 修正、重新求 FEM、入口改动或任何动力学调参。

本机保护清单涵盖 206 个原有包内文件（含全部源码及既有缓存），最终变更与新增均为 0；服务器独立快照 114 个源码/源码资产文件也为 0 变更。旧服务器目录少了若干后来增加的诊断模块，故本轮同步并验证了完整本地快照，而未使用旧服务器目录静默替代。

# 30 个 paired Bubble 是什么

完全恢复历史 **P9-A.1 smoke30** 的 30 个稳定 ID，master seed=202609249。沿用历史 Method B 外部 inlet 的已接受事件；没有重新生成 population，没有使用 Method C 或 P9-A.3 interior injection。历史随机分位数没有单独落盘，但精确 anchor、position seed、完整事件元数据和已保存的状态都可恢复，所以无需重抽样。

唯一永久 cohort：[paired_cohort_manifest.csv](data/paired_cohort_manifest.csv)，SHA `{s['cohort_manifest_sha']}`。其配套 [paired_events.json](data/paired_events.json) 和 [paired_initial_states.json](data/paired_initial_states.json)保存全部状态。OLD/NEW 的 finite sphere admission、fluid inclusion、wall clearance、open-cap 规则与 initial tetra 对所有 30 颗都相同且可接受。

历史 OLD 的 flow SHA、当时全部 {cohort['historical_science_files_verified']} 个 Python 源码哈希、各事件哈希、初始数组及轨迹文件哈希均匹配，故直接复用其完整轨迹与诊断日志，原输出不改写。配对表并未把两批不同抽样的泡混在一起。

# 同一颗 Bubble 的路线发生了什么

{transitions}
改变出口的 ID：{', '.join(map(str,s['outlet_changed_ids']))}。没有 O3→O2，也没有 stationary 与 completed 之间的转变。[逐泡配对表](data/paired_outcome_transition.csv)保留所有 30 行。

![配对轨迹](figures/02_paired_old_new_trajectories.png)

Figure 02 的八颗示例按 outcome 类别、稳定 ID 和直径跨度选取，规则在 [representative_selection.json](data/representative_selection.json)，未按轨迹外观挑选。蓝色虚线为 OLD，橙色实线为 NEW；重合段可能互相覆盖，不表示丢失 OLD 轨迹。

![迁移矩阵](figures/03_paired_outcome_transition.png)

# stationary 怎么变了

OLD 的唯一 stationary 为 **ID 15，D=2.489958 µm**，NEW 中仍 stationary，恢复数 0，新增加数 0。末态位置约 (92.8863, 46.6312, 114.7719) µm，新旧差距处于几何舍入预算内；wall gap 约 2 nm，NEW 局部流速约 6.0812 mm/s，并非当地流体停滞。

两组末段连续 32 个已接受求解均有 3 个独立接触约束，rank=3、非负乘子，平移速度处于原求解器 velocity budget 内。证据支持当前刚性有限尺寸模型下的接触静止结果，不是软件失败，也不等于体内生理捕获的实验证据。没有做虚拟半径或变形实验。

# transit time 怎么变了

仅比较两组均 completed 的 29 颗：平均 **{s['mean_transit_time_old']*1e3:.3f}→{s['mean_transit_time_new']*1e3:.3f} ms**；配对差值中位数 +{s['paired_transit_time_statistics']['delta_s']['median']*1e3:.3f} ms，NEW/OLD 比值中位数 {s['paired_transit_time_statistics']['new_over_old']['median']:.3f}。28 颗变慢，1 颗变快。ID 18 虽仍走 O3，时间由 367.163 ms 降至 86.910 ms，所以均值和中位数描述不同。stationary 的 33 ms 等已模拟时间不当作 transit time。

# near-wall exposure 怎么变了

30 颗各自最小 wall gap 的中位数在两组都约 2 nm；均值 **3.618→2.320 nm**。这与原 molecular/handoff floor 保持不变相容，不能据此宣称发生穿透。最小 h/a 中位数约 0.003344→0.003284。基于已保存接受步的 h/a≤0.1 暴露时间，均值约 **10.713→17.543 ms**；这只是离散路径诊断，不是新积分或新壁面模型。

# point tracer 和 Bubble 有什么区别

没有找到与本 cohort 30 个中心逐一匹配的既存 NEW tracer 记录，因此只从这些相同中心运行原有 double-precision tetra-adjacency P1 RK23 诊断器，step=0.2 µm、error=1e-11、path horizon=2 mm。没有重新审核全场。全部 30 tracer 到达原始出口三角形，并由原分类器复核。

{pointtable}
**3 颗两者均完成但出口不同**：ID 26、27（point O2，MB O3），ID 46（point O1，MB O2）。另有 ID 15：point O3，MB stationary。因此总 outcome difference 为 4，不能把 stationary 也计作“换出口”。这表明同一初始位置的有限尺寸与壁面动力学能进一步改变背景流路由。

# 有没有 penetration / solver failure

NEW 的 penetration=0、采样态 handoff violation=0、全部连续证书（含分段 union 的子证书）违规=0、NaN/Inf=0、未分类 solver corruption=0、solver fail=0、inlet escape=0。ID 15 的 stationary 单独按接触证据审核，未伪装成 completed。

1 worker 与 3 workers 在 ID 3/15/22 上的完整 samples、终点、状态、出口、时间和解压事件日志逐位一致，正式 6 workers 的同三颗也一致。永久测试 **{tests['tests']} passed，{tests['elapsed_s']:.2f} s**。首次新测试只识别单段证书，遇到 union 证书产生 KeyError；已保留初次日志并补充逐段审核，未改动力学、数值参数或轨迹。

NEW 30 泡服务器轨迹阶段 **{runtimes['NEW']:.2f} s（6 workers）**；同中心 point tracer 阶段 {runtimes['NEW_POINT']:.3f} s。1-worker/3-worker 预检分别 {runtimes['determinism_w1']:.2f}/{runtimes['determinism_w3']:.2f} s，并行启动；这些阶段计时不含输入加载、传输、制图和报告时间。OLD 本轮复用，不重复积分。

# P1 局部 conservation limitation

NEW 总流量与分流验证通过，但 **P1/P1+VMS 局部守恒仍不完美**。实际已有审计的 599 个合法 root 截面（612 个候选中 13 个按几何排除）：max flux deviation **3.454212%**，RMS **1.789924%**。引用的原始 JSON/CSV 已按 SHA 复制为 [known flow limitation](data/known_flow_limitation_sources.json)。本轮接受并显式保留这个模型限制，没有修复或掩盖它；也未运行 P2 对照。

# 这一轮能说明什么

在原有 Particle 科学实现、几何、固定初始状态和容差下，新出口边界条件得到的冻结 3D 流场如何改变这一组微泡的出口、通过时间及近壁经历；并支持当前 30 泡的数值安全性与 worker 确定性。

# 这一轮不能说明什么

N=30 是历史固定 paired validation cohort，不是新 entering population；其 O1/O2/O3 数量不能当作真实微泡分流，更不能与 fluid 6.54/44.82/48.64% 做统计等价判断。网络推导边界条件没有因本轮通过而成为 physiological ground truth。stationary 也不能单独证明体内捕获或排除真实可变形效应。

# 是否建议进入正式大样本

建议 **`{s['large_sample_recommendation']}`**：数值安全门槛、确定性与失败审计通过，仍需用户人工审核本报告和 ID 15 的模型含义后决定下一阶段。**本轮未启动 500 MB production。**

全部机器可读结果：[final_summary.json](data/final_summary.json)。详细[配对审核](PAIRED_OLD_NEW_TRAJECTORY_REVIEW_ZH.md)及[安全审核](NEW_FLOW_MB_SAFETY_AUDIT_ZH.md)。服务器独立工作区：`{ctx['remote']}`。
'''
    (R/'NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md').write_text(main)
    pairrows=[[t['bubble_id'],f"{t['diameter_um']:.4f}",t['old_outcome'],t['new_outcome'],ms(t['old_transit_time_s']),ms(t['new_transit_time_s']),f"{t['old_min_gap_m']*1e9:.4f}",f"{t['new_min_gap_m']*1e9:.4f}"] for t in trans]
    detail='# 逐泡配对审核\n\n同一颗泡只替换背景流；N=30，非 population statistics。所有旧结果 SHA 和完整事件元数据匹配才复用。\n\n'+table(['ID','直径 µm','OLD','NEW','OLD transit ms','NEW transit ms','OLD min gap nm','NEW min gap nm'],pairrows)
    detail+='\n## 指标语义\n\nCOMPLETED 只用于原出口三角形首次中心穿越；stationary 保留为 STATIONARY，transit time 留空，物理已模拟时间单独记录。contact_events 统计已接受诊断序列中“无接触→有接触”的连续接触段起点；contact_accepted_steps 与最大同时约束数另外保存。handoff_events 只计已接受事件。solver iteration statistics 指 contact active-set 迭代，直接线性求解无迭代数时不编造该指标。planar wall correction 是已存壁面权重乘以目标切向速度与 bulk 切向速度之差的模，另存全部 unconstrained hydro−free 平移修正。所有这些指标只读取原诊断。\n\n'
    detail+='原始 20 列轨迹（含速度、角速度、姿态、gap、g_nf、h_lower、状态、dt）和完整 accepted/rejected 诊断保存在 `outputs/OLD`、`outputs/NEW`。所有连续 union 证书的分段范围必须无缝覆盖 [0,1]，逐个叶证书审核。\n\n'
    detail+='## point tracer 差异\n\n'+table(['ID','point','MB','出口改变（两者都完成）'],[[r['bubble_id'],r['point_outcome'],r['mb_outcome'],r['both_exited_route_differs']] for r in point if r['outcome_differs']])
    detail+='\n## 图片\n\n'+''.join(f'[{Path(a["png"]).stem} PNG](figures/{Path(a["png"]).name}) · [PDF](figures/{Path(a["pdf"]).name})\n\n' for a in artifacts)
    (R/'PAIRED_OLD_NEW_TRAJECTORY_REVIEW_ZH.md').write_text(detail)
    cert_count=sum(a['continuous_certificate_count'] for a in audit if a['flow']=='NEW');top_count=sum(a['continuous_top_level_certificate_count'] for a in audit if a['flow']=='NEW')
    safety=f'''# NEW 30 泡安全审核

最终安全门槛 PASS。穿透、采样态 handoff 违规、连续证书违规、NaN/Inf、未分类 solver corruption、solver fail 均为 0；无 inlet escape。所有 30 个结果都有终态记录，没有删除不完成轨迹。

共检查 {top_count} 个顶层连续证书，递归得到 {cert_count} 个原始 support-plane-minus-h_lower 叶证书。每个 union 的子区间无重叠空隙地覆盖 [0,1]，且声明 held-velocity path 未改；每个叶证书的 minimum_g_nf_bound 均满足原 roundoff budget。采样态 gap≥−roundoff 与 g_nf≥−roundoff 分别复核；roundoff=2.009411463123824e-17 m。

## 唯一 stationary：ID 15

OLD/NEW 中均为 D=2.489958478 µm。最终 32 个已接受求解具有 rank=3 的三个独立接触约束，乘子非负，平移速度小于原 KKT velocity budget。末态坐标约 (92.886303,46.631184,114.771900) µm，gap≈2 nm；NEW 局部流速 6.081224 mm/s。故支持当前刚性有限尺寸接触静止，不能声称流体停滞、软件失败、体内捕获证据或真实微泡必然无法通过。无半径替代实验。

## 确定性与科学保护

1/3 workers 的 ID 3、15、22 全 samples 与解压事件日志逐位一致；永久测试又与 6 workers 正式结果逐位比较。12 项测试全部通过。首次测试的证书 schema KeyError 是新观察代码对 union 证书支持不足，原始失败日志已保留于 `logs/initial_tests.log/xml`，修复限于本轮新测试与本轮新审计器。

本地包内原有 206 个文件、服务器源码快照 114 个文件均 SHA 未改变。NEW/OLD formal VTU、复用 OLD 轨迹及原报告均未覆盖，新增报告、adapter 和测试位于独立新目录。未执行任何科学源码变更、FEM 修正或新增 500 泡任务。

## 证据

[所有末态及原求解器诊断](data/saved_state_safety_audit.json)、[NEW 指标](data/new_paired_metrics.json)、[12 tests](logs/tests.log)、[science protection](data/particle_science_protection_final.json)、[remote snapshot](logs/remote_final_snapshot_verification.json)、[worker determinism](data/worker_count_determinism_subset.json)。

NEW P1 局部截面误差 max=3.454212%、RMS=1.789924% 仍是本实验已知模型限制。建议 `{s['large_sample_recommendation']}`，需人工审核后再决定大样本。
'''
    (R/'NEW_FLOW_MB_SAFETY_AUDIT_ZH.md').write_text(safety)
    preview='<!doctype html><html lang="zh"><meta charset="utf-8"><title>Network flow paired 30 review</title><style>body{font:17px/1.7 sans-serif;max-width:1200px;margin:30px auto;background:white;color:#222}img{max-width:100%}section{margin:35px 0}code{overflow-wrap:anywhere}</style><h1>Network-derived flow：30 泡配对验证</h1>'
    preview+='<p>OLD / NEW 均 29 completed、1 stationary、0 solver fail；19 颗 O2→O3。N=30 不能当作 population split。</p><p><a href="NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md">中文主报告</a> · <a href="data/final_summary.json">机器可读汇总</a></p>'
    for a in artifacts:
        n=Path(a['png']).name;stem=Path(n).stem;preview+=f'<section><h2>{html.escape(stem)}</h2><a href="figures/{stem}.pdf">PDF</a><br><img src="figures/{n}" alt="{html.escape(stem)}"></section>'
    (R/'OPEN_RESULTS.html').write_text(preview+'</html>')
    readme=f'''# 可复现入口

本地原科学包：`{R.parents[1]/'src/particle_3d'}`（只读）。新增 adapter/分析/绘图位于本目录 `scripts/`；永久测试位于 `particle_3d/tests/network_flow_mb_validation_v1/`。

服务器工作区：`{ctx['remote']}`。输入包 SHA：`{remote['input_bundle_sha']}`。该工作区每次运行前逐一验证 `data/code_snapshot_manifest.json`、flow SHA、共享 cohort。所有 30 泡事件固定，不执行 sampler。

已执行的命令形态：
```
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1 TMPDIR={ctx['remote']}/tmp /root/particle8_2_runs/env/bin/python {ctx['remote']}/scripts/runner.py --root {ctx['remote']} --label NEW --workers 6 --ids all --output NEW
```
既有 output 名存在时 runner 主动拒绝覆写。复跑必须使用新的、独立的 workspace/output；不要重跑历史路径。determinism 用 `--ids 3,15,22 --workers 1/3`；point 用 `--mode point --workers 1`。本轮 OLD 通过源码/flow/event/样本 SHA 核验复用，记录在 cohort_provenance.json。

本地只读永久检查：
```
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python -m pytest -q -p no:cacheprovider {R.parents[1]/'tests/network_flow_mb_validation_v1'}
```
本目录 `archive/` 保存原始传输包、Taylor-Hood 最后已完成单步 smoke，以及修正图标签之前的原图。`server_bundle/` 保留科学源码快照与 frozen_reference。所有正式输出均由未经修改的 Particle 科学模块生成。
'''
    (R/'REPRODUCE.md').write_text(readme)
    print(json.dumps(dict(status=s['final_status'],tests=tests,figures=len(artifacts),source_changes=0,report=str(R)),indent=2))
if __name__=='__main__':main()
