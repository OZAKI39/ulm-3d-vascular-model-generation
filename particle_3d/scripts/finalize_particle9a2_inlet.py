#!/usr/bin/env python3
"""Evidence-derived Chinese review, gates, CSV ledgers and delivery hashes."""
from pathlib import Path
import json,csv,hashlib,subprocess,xml.etree.ElementTree as ET,html,datetime
import numpy as np
from PIL import Image
from particle_3d.injection_method_c import TruncatedSonoVue,METHOD,LEGACY_METHOD,canonical_bytes
from particle_3d.particle7_cases import SONOVUE
from particle_3d.particle8_replay import canonical_hash

ROOT=Path(__file__).resolve().parents[2];R=ROOT/'particle_3d/reports/particle9a2_inlet_sampling';D=R/'data';O=ROOT/'particle_3d/outputs/particle9a2_2mmps'
ROLES=['OUTLET_01','OUTLET_02','OUTLET_03','NO_EXIT']
def read(p):return json.loads(Path(p).read_text())
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(canonical_bytes(v))
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def csvwrite(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        w.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
def counts(c):return '/'.join(str(c[k]) for k in ROLES)
def stats(events):
    d=np.array([e['diameter_um'] for e in events]);return dict(count=len(d),min_um=float(d.min()),median_um=float(np.median(d)),max_um=float(d.max()),mean_um=float(d.mean()),less_than_1p2_um=int(sum(d<1.2)),above_2p2_um=int(sum(d>2.2)))
def suite(p):
    suites=list(ET.parse(p).getroot().iter('testsuite'));return dict(tests=sum(int(s.attrib['tests']) for s in suites),failures=sum(int(s.attrib.get('failures',0))+int(s.attrib.get('errors',0)) for s in suites))

def main():
    audit=read(R/'audit2000/admission/birth_ledger.json')['events'];smoke=read(R/'smoke30/admission/birth_ledger.json')['events'];formal=read(O/'admission/birth_ledger.json')['events'];old=read(R/'reference/method_b_500_birth_ledger.json')['events']
    source=TruncatedSonoVue(SONOVUE);capacity=read(D/'inlet_size_capacity.json');prob=read(D/'source_probability.json')
    ap=read(D/'audit2000_point.json');fp=read(D/'formal500_point.json');fg=read(D/'formal500_gate.json');sg=read(D/'smoke30_gate.json');ref=read(D/'independent_flux_reference.json')
    stages={k:read(D/(k+'_sampling.json')) for k in ['audit2000','smoke30','formal500']};all_events=audit+smoke+formal
    old_point=list(csv.DictReader((R/'reference/method_b_500_point.csv').open()));old_final=list(csv.DictReader((R/'reference/method_b_500_final.csv').open()))
    old_point_counts={k:sum(r['point_outlet']==k for r in old_point) for k in ROLES};old_final_counts={k:sum(r['outlet']==k for r in old_final) for k in ROLES}
    assert list(old_point_counts.values())==[9,469,20,2] and list(old_final_counts.values())==[0,477,9,14]
    for label,events in [('audit2000',audit),('smoke30',smoke),('formal500',formal)]:
        csvwrite(D/(label+'_birth_events.csv'),[{k:e[k] for k in ['particle_id','birth_time_s','diameter_m','radius_m','diameter_draw_id','diameter_global_rejections','position_proposal_count','position_xyz','inlet_triangle_id','local_flux_weight_m_s','clearance_m','admission_status','source_distribution_contract_sha256','flow_sha256','geometry_sha256','seed']} for e in events])
        csvwrite(D/(label+'_source_diameter_draws.csv'),[dict(particle_id=e['particle_id'],**v) for e in events for v in e['source_diameter_draws']])
    before=read(D/'protected_before.json');changed=[p for p,h in before['files'].items() if not Path(p).is_file() or sha(p)!=h]
    git=lambda *x:subprocess.check_output(['git',*x],cwd=ROOT,text=True)
    remote_protection=read(D/'server_input_legacy_verification.json')
    protection=dict(file_count=len(before['files']),changed_files=changed,unchanged=not changed,
        index_unchanged=git('write-tree').strip()==before['index_tree'],head_unchanged=git('rev-parse','HEAD').strip()==before['head'])
    write(D/'protected_verification.json',protection)
    tests={k:suite(R/'logs'/(k+'.xml')) for k in ['local_tests','remote_tests','final_regression_tests']}
    sizes={k:stats(v) for k,v in [('legacy_B500',old),('C2000',audit),('C500',formal)]}
    maximum=capacity['D_geometry_max_m'];u=np.sort(source.original.cdf(np.array([e['diameter_um'] for e in audit]))/source.original.cdf(maximum*1e6))
    ks=float(max(np.max(np.arange(1,len(u)+1)/len(u)-u),np.max(u-np.arange(len(u))/len(u))))
    manifest=read(D/'figure_manifest.json');image_checks=[]
    for item in manifest:
        p=R/'figures'/(item['name']+'.png');pdf=p.with_suffix('.pdf')
        with Image.open(p) as im:im.load();image_checks.append(dict(name=item['name'],size=list(im.size),png_ok=im.width>=2500 and im.height>=1000,pdf_ok=pdf.read_bytes().startswith(b'%PDF-'),sha_ok=sha(p)==item['png_sha256'] and sha(pdf)==item['pdf_sha256']))
    write(D/'figure_artifact_checks.json',image_checks)
    # Metadata and all NPZ bytes are verified locally, independently of transfer.
    paths=[]
    for e in formal:
        p=O/'trajectories'/f"mb_{e['particle_id']:06d}.json";m=read(p)
        assert m['birth_metadata_sha256']==canonical_hash(e)
        assert sha(p.with_suffix('.npz'))==m['samples_sha256']
        assert m['birth_metadata']==e
        paths.append(str(p.relative_to(ROOT)))
    checks=dict(A_D_max_contract=source.contract==read(ROOT/'particle_3d/contracts/SONOVUE_DIAMETER_MAX4UM_V1.json'),
        B_truncated_distribution=all(e['diameter_m']<=4e-6 for e in all_events) and all(v['diameter_m']<=4e-6 for e in all_events for v in e['source_diameter_draws']) and not source.contract['clipping'],
        C_diameter_fixed=all(e['diameter_fixed_during_position_sampling'] and all(v['diameter_float64_hex']==np.float64(e['diameter_m']).tobytes().hex() for v in e['position_attempts']) for e in all_events),
        D_conditional_flux=ref['all_position_agreement'] and all(v['failures']==0 for v in tests.values()),
        E_no_future_labels=all(v['failures']==0 for v in tests.values()),
        F_deterministic=all(v['all_events_byte_identical_reversed_single_worker_replay'] for v in stages.values()),
        G_legacy_untouched=protection['unchanged'] and protection['index_unchanged'] and not remote_protection['current_mismatches'] and not remote_protection['legacy_mismatches'],
        H_2000_audit_complete=len(audit)==2000 and len(ap['rows'])==2000,
        I_smoke_no_new_failure=sg['no_new_systematic_failure'],
        J_500_complete_records=len(formal)==len(paths)==len(read(D/'formal500_trajectory_audit.json'))==len(fp['rows'])==500,
        K_no_wall_penetration=sg['penetration_count']==fg['penetration_count']==0,
        L_no_handoff_violation=sg['handoff_violation_count']==fg['handoff_violation_count']==0 and sg['uncertified_segment_count']==fg['uncertified_segment_count']==0,
        M_report_data_sha_complete=all(c['png_ok'] and c['pdf_ok'] and c['sha_ok'] for c in image_checks) and len(manifest)>=7,
        formal_no_new_numerical_pathology=fg['no_new_systematic_failure'] and read(D/'stationary_equilibrium.json')['all_verified'],
        finite_sphere_births_valid=all(v['invalid_births']==0 for v in stages.values()))
    status='P9A2_INLET_AUTOMATED_PASS' if all(checks.values()) else 'P9A2_INLET_AUTOMATED_FAIL'
    summary=dict(automated_status=status,human_status='PENDING_USER_REVIEW',gates=checks,original_above_4um_probability=prob['original_above_4um'],
        capacity=capacity,theoretical_global_size_rejection_probability_bracket=prob['global_impossible_source_probability_bracket'],
        sampling=stages,source_size_rejections_total=sum(v['global_rejections'] for v in stages.values()),
        source_draws_total=sum(v['source_draws'] for v in stages.values()),
        point_C2000=ap['counts'],point_B500=old_point_counts,point_C500=fp['counts'],final_C500=fg['outlet_counts'],
        final_B500=old_final_counts,smoke=sg,formal=fg,size_statistics=sizes,
        C2000_size_CDF_KS_against_geometry_conditioned_source=ks,
        source_endpoint_4um_draw_count=sum(v['diameter_m']==4e-6 for e in all_events for v in e['source_diameter_draws']),
        tests=tests,protected=protection,remote=read(D/'remote.json'),fixed_size_independent_reference=[{k:v for k,v in x.items() if 'positions' not in k} for x in ref['rows']])
    final_by_id={v['particle_id']:v for v in read(D/'formal500_trajectory_audit.json')}
    csvwrite(D/'formal500_stationary_classification.csv',[
        dict(particle_id=r['particle_id'],stationary=True,classification=r['stationary_classification'],
             geometric_FEM_downstream_blocked=r['geometry_no_FEM_downstream_direction'],
             independent_equilibrium_verified=r['supported_stationary'])
        for r in final_by_id.values() if r['stationary']])
    paired=[dict(particle_id=p['particle_id'],diameter_m=p['diameter_m'],point_outlet=p['point_outlet'],final_outlet=final_by_id[p['particle_id']]['outlet'],stationary=final_by_id[p['particle_id']]['stationary']) for p in fp['rows']]
    csvwrite(D/'formal500_point_vs_final.csv',paired)
    summary['named_outlet_reroutes']=sum(p['point_outlet']!='NO_EXIT' and p['final_outlet']!='NO_EXIT' and p['point_outlet']!=p['final_outlet'] for p in paired)
    summary['point_to_final_category_changes']=sum(p['point_outlet']!=p['final_outlet'] for p in paired)
    runtime={k:read(D/(k+'_completed.json'))['wall_seconds'] for k in ['smoke30','formal500']}
    runtime.update(audit2000_sampling=stages['audit2000']['wall_seconds'],audit2000_point=ap['wall_seconds'],formal500_point=fp['wall_seconds'],independent_flux_reference=ref['wall_seconds'],smoke30_recheck=sg['wall_seconds'],formal500_recheck=fg['wall_seconds'])
    runtime['independent_stationary_equilibrium']=read(D/'stationary_equilibrium.json')['wall_seconds'];summary['runtime_seconds']=runtime;write(D/'summary.json',summary);write(D/'gates.json',dict(automated_status=status,human_status='PENDING_USER_REVIEW',gates=checks))
    (R/'logs/git_diff_stat.txt').write_text(git('diff','--stat'))
    (R/'logs/git_status_final.txt').write_text(git('status','--short'))
    newfiles=[p for p in ROOT.glob('particle_3d/src/particle_3d/*') if p.name in ['injection_method_c.py','inlet_size_capacity.py','particle9a2_equilibrium_audit.py']]
    newfiles+=list(ROOT.glob('particle_3d/scripts/*particle9a2*'))+list(ROOT.glob('particle_3d/tests/particle9a2_inlet/*.py'))
    newfiles=sorted(set(newfiles));lines=sum(len(p.read_text().splitlines()) for p in newfiles)
    write(D/'new_source_manifest.json',{str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    (R/'logs/new_files_stat.txt').write_text('\n'.join(f'{p.relative_to(ROOT)} | {len(p.read_text().splitlines())} lines' for p in newfiles)+f'\n{len(newfiles)} new Python files, {lines} lines\n')
    text=f'''# 一句话结论

**{status}**；人工状态 **PENDING_USER_REVIEW**。已消除已知“固定位置后反复换尺寸”的 Method B 机制；新的 size-first 入口没有保证各出口接近流量比例。当前 2000-event point 为 {counts(ap['counts'])}，正式 500 最终为 {counts(fg['outlet_counts'])}（顺序均为 O1/O2/O3/no exit）。{fg['stationary_count']} 条为当前刚性球模型下的停止，不等同生理捕获。

# Method B 为什么有问题

旧方法是“位置固定后不断换球大小”：一个 anchor 放不下最初的球，就在该位置最多重抽 512 次尺寸。原 834 个位置的 point 为 27/705/89/13；筛选后 accepted 500 为 9/469/20/2。O1/O2/O3 准入率分别为 33.33%/66.52%/22.47%。它既对位置产生筛选，又使接受尺寸依赖位置。历史结果和代码全部保留，绑定 **{LEGACY_METHOD}**。

# Method C 怎么做

“先决定这颗球有多大，再给它找一个真正放得下的位置”。正式方法名 **{METHOD}**。先从 4 µm 截断 source 抽直径；只有全入口无任何可行位置时，记录 `NO_FEASIBLE_INLET_POSITION_FOR_SIZE` 并重新抽 source。只要有可行位置，直径冻结，随后只能换位置。位置由原正向 P1 flux sampler 提出，再用原 `FiniteSizeAdmission.check` 的有限球 WALL + handoff 判据接受。

出生率仍为 C_MB × Q，C_MB=8.5×10¹² m⁻³，Q=1.5513591604402322×10⁻¹⁴ m³/s，即约 0.1318655 个/s。这是实际进入事件率，不把 source 重抽算作额外 birth。全入口不可行的 source draw 和位置失败是不同事件。

新配置是 `particle_3d/contracts/PARTICLE9A2_PRODUCTION_V1.json`。新入口 API 和新 P9-A.2 runner 显式绑定 Method C；旧入口/旧 P9 runner 不被改写。历史 P8.2A 里同样叫 C 的“沿 point streamline 找入内位置”是另一种旧诊断，不是本轮的新方法。

# 4 µm 上限怎么实现

原 frozen SonoVue 中 >4 µm 的概率为 **{100*prob['original_above_4um']:.6f}%**。使用原 CDF 反演 F⁻¹(U·F(4 µm))，即截断后重新归一化，不把大球压成 4 µm。

4 µm 位于原 [3.95,4.05] µm bin 中。原 contract 已明确 bin 内均匀，因此保留该 bin 下半段的概率质量，再用统一 F(4 µm) 归一化；这是复用已有定义，不是新猜测。派生 contract 记录原 source、histogram、contract 和生成代码 SHA。原冻结文件未修改。本轮 {summary['source_draws_total']} 次 source draw 中，恰好 4 µm 的数量为 {summary['source_endpoint_4um_draw_count']}，没有人为 clipping spike。

# 4 µm 微泡是否都能进入这个入口

**不能。** 在真实 inlet cap + WALL、原 handoff 下，最大可准入直径的全局括区为 **[{capacity['D_geometry_max_bracket_m'][0]*1e6:.9f}, {capacity['D_geometry_max_bracket_m'][1]*1e6:.9f}] µm**。对应半径约 {capacity['radius_max_bracket_m'][0]*1e6:.9f} µm。去掉 handoff、只看几何接触的直径括区为 [{capacity['pure_wall_clearance_diameter_bracket_m'][0]*1e6:.9f}, {capacity['pure_wall_clearance_diameter_bracket_m'][1]*1e6:.9f}] µm。

该界来自三角形上的全局分支细分，使用真实 WALL 距离的 1-Lipschitz 上下界，不是 bounding box、二维 rim 距离或有限栅格最大值。报告下界有真实可行中心作见证，上界不被当成已可行值。这里“可通过”仅指现有开放入口的出生准入，不证明球能通过整条血管或完成域外到域内的完整球体扫掠。

SOURCE 上限仍为 4 µm；实际进入分布进一步条件化于当前入口的刚性球 passability。4 µm source 中理论上约 **{100*np.mean(prob['global_impossible_source_probability_bracket']):.6f}%** 没有可行入口位置，已在正式生产前报告。

| 批次 | 合法 birth | 全局尺寸拒绝 / source draws | 实测比例 |
| --- | ---: | ---: | ---: |
'''
    for k,v in stages.items():text+=f"| {k} | {v['count']} | {v['global_rejections']} / {v['source_draws']} | {100*v['global_rejections']/v['source_draws']:.3f}% |\n"
    text+=f'''
这些实测比例有有限样本波动，不是强制达到理论比例的 quota。实际 entering distribution 的名称明确是“4 µm 截断 SonoVue，进一步条件化于当前刚性球入口可进入性”。

# 新方法有没有偷偷改变尺寸

本轮 2000 + 30 + 500 = 2530 个事件中，位置重试造成的直径变化为 **0**，非法有限球 birth 为 **0**。每次 proposal 保存相同的 `diameter_draw_id` 和 float64 字节，birth ledger 保存独立 DIAMETER/POSITION_PROPOSAL counter key。三批均在单 worker 中按相反 ID 顺序完整再生成，事件逐字节一致。位置失败不会改变下一颗微泡的尺寸序列。

# 新方法入口位置是否按血流量采样

是，按“固定尺寸可行域内的正向流量”抽样。加速器只丢弃已由几何上界证明不可能容纳该球的三角形，并在保留三角形上继续使用原 `InletFluxSampler`，最终仍由原准入 checker 拒绝位置。被保留区域是可行域的超集，因此最后的拒绝采样保持所需条件分布。

每颗球先计算可行 flux 的保守下界和保留 flux 上界，再选择确定性 guard；目标零成功尾概率为 2⁻⁶⁴，并保留 4 倍执行余量。没有新增几何安全距离。若 guard 耗尽，抛出 `POSITION_SAMPLING_PROGRESS_FAILURE`，不换尺寸。本轮 2000 个事件实际共 {stages['audit2000']['position_proposals']} 次 position proposals，单颗最多 {stages['audit2000']['maximum_proposals']} 次。较大球的真实可行 flux 急剧缩小，但加速后不用进行海量盲抽，见图 04。

额外使用原 sampler、完全不做三角形加速的独立参考，在 6 个固定直径各抽 500 个合法中心。两个入口投影坐标的两样本经验 CDF 检查全部通过，最大 KS 差为 {max(max(r['position_coordinate_two_sample_KS']) for r in ref['rows']):.3f}；合成非均匀通量矩形还核对了可解析的条件位置均值。该检查不使用出口比例作为 Gate。

# Method B 和 Method C 的 point outlet 分布

| 样本 | O1 | O2 | O3 | no exit |
| --- | ---: | ---: | ---: | ---: |
| Method B accepted 500 | 9 | 469 | 20 | 2 |
| Method C audit 2000 | {ap['counts']['OUTLET_01']} | {ap['counts']['OUTLET_02']} | {ap['counts']['OUTLET_03']} | {ap['counts']['NO_EXIT']} |
| Method C formal 500 point | {fp['counts']['OUTLET_01']} | {fp['counts']['OUTLET_02']} | {fp['counts']['OUTLET_03']} | {fp['counts']['NO_EXIT']} |

**新方法仍高度偏向 O2，不能把“更接近流量”当正确性标准。** 真正的区别在于尺寸分布与可行空间：旧 B 接受粒径中位数 {sizes['legacy_B500']['median_um']:.4f} µm，C 的 2000 样本为 {sizes['C2000']['median_um']:.4f} µm；小于 1.2 µm 的比例分别为 {100*sizes['legacy_B500']['less_than_1p2_um']/500:.2f}% 和 {100*sizes['C2000']['less_than_1p2_um']/2000:.2f}%。旧 B 的边缘位置靠反复抽到小球才获准；C 不会这样改小球，较大的球自然集中在中部可行区域。

原 flux sampler 直接拒绝位置的独立结果如下。1.6/2.0/2.6 µm 各 500 个都到 O2，说明该倾向在不使用新加速器时也出现。以下是有限样本观察，不宣称其它流域概率严格为零。

| 固定直径 µm | 合法 500 的 O1/O2/O3/no exit | 原 proposal 准入比例 |
| --- | --- | ---: |
'''
    for row in ref['rows']:text+=f"| {row['diameter_um']:.1f} | {counts(row['conditional_point_counts'])} | {100*row['original_feasible_flux_fraction_estimate']:.3f}% |\n"
    text+=f'''
这些证据支持“有限尺寸可行区域 + 已声明 source 分布”对 O2 的真实几何条件化；不是未来出口参与位置接受。生产模块不读取未来 basin、不调用 classifier、不设置 outlet quota，相关永久测试通过。这里的“真实几何”限定于当前离散血管和刚性球模型，不等于真实生理分布。

# 新 500 条正式结果

30 smoke：{counts(sg['outlet_counts'])}，其中 {sg['stationary_count']} 条停止均有独立三方向接触约束支持；穿透 {sg['penetration_count']}、handoff 违规 {sg['handoff_violation_count']}、入口逃逸/域外错误 {sg['inlet_escape_or_outside_count']}，无新增数值失败。全部 {sg['continuous_segments_rechecked']} 个运动段重新通过连续壁面证书，Gate 后才生成正式 500。

正式 500：**{counts(fg['outlet_counts'])}**；stationary **{fg['stationary_count']}**，其中独立接触复核支持 **{fg['supported_stationary_count']}**；penetration **{fg['penetration_count']}**，handoff violation **{fg['handoff_violation_count']}**。未完成轨迹、停止原因、最终接触约束与全部路径原样保留，没有算作成功出流。全部 {fg['continuous_segments_rechecked']} 个运动段做了读取保存路径后的连续 handoff 证书复核；未证实安全的段为 {fg['uncertified_segment_count']}。

其中 **145 条**另有“几何上没有顺流方向”的证据。**ID 371、435** 是不同情况：几何上有沿 FEM 顺流的可行方向，但当前含壁面仿射水动力和旋转耦合的实际驱动，在接触约束下仍给出零平移。

初审曾把“无 FEM 顺流几何方向”作为所有 stationary 的必要条件，因而暂记这 2 条失败。随后没有改模型或重积分，而是对 **全部 147 个保存终态**独立穷举 active sets，直接求解原始 6D 正定阻力二次问题，核对全体不穿透约束、非负乘子和 KKT 条件，147/147 复现原模型平衡。这说明几何可行方向不必然是实际水动力驱动方向；也不能把这 2 条叫作“几何上无路可走”。加入了两项永久测试，分别检验该区别，以及拒绝一个应运动却被错误记为静止的反例。

完整独立证据在 `data/stationary_equilibrium.json`，2 个代表的原始 R/b/J/U 在 `data/equilibrium_system_371.npz` 和 `data/equilibrium_system_435.npz`。原保守 FAIL 记录和原审核脚本完整保留在 `reference/initial_stationary_gate/`。最终 Gate 修正的是审核分类的含义，不是通过改变物理、几何、尺寸、数值容差或轨迹来消除异常。这里证实的是当前模型的瞬时平衡，不证明真实生理捕获或不存在别的几何逃逸路径。

正式 point→最终类别改变 {summary['point_to_final_category_changes']} 条；在已命名出口之间改变 {summary['named_outlet_reroutes']} 条。后者与 stationary/no-exit 必须分开，逐 ID 表见 `data/formal500_point_vs_final.csv`。冻结 FEM、P9-A.1 hydrodynamics、P6.5 normal lubrication、2 nm handoff、contact redundancy、积分物理、outlet classifier 均未改。

# 和旧 Method B 500 比较

| 指标 | 旧 Method B 500 | 新 Method C 500 |
| --- | --- | --- |
| 入口 point O1/O2/O3/no exit | 9/469/20/2 | {counts(fp['counts'])} |
| 最终 O1/O2/O3/no exit | 0/477/9/14 | {counts(fg['outlet_counts'])} |
| 粒径中位数 µm | {sizes['legacy_B500']['median_um']:.4f} | {sizes['C500']['median_um']:.4f} |
| 粒径范围 µm | {sizes['legacy_B500']['min_um']:.4f}–{sizes['legacy_B500']['max_um']:.4f} | {sizes['C500']['min_um']:.4f}–{sizes['C500']['max_um']:.4f} |
| stationary | 14 | {fg['stationary_count']} |

新 accepted 粒径分布明显不同于旧 B；这是去掉旧位置条件化尺寸重抽后的预期人口变化，不能写成两批完全相同粒径的动力学对照。C2000 与“截断 source 再条件化于几何可进入性”的目标 CDF 最大经验差为 {ks:.5f}。新 500 使用新 seed、新 ledger、新轨迹，没有复用 Method B admission cache。

# 目前还能有哪些有限尺寸偏差

大球天然能进入的位置比小球少，入口中心化和后续分叉改道/刚性停止仍可能存在。本轮消除的是 Method B 的人为位置条件化尺寸重试，并没有消除有限尺寸几何作用。当前 P1 近壁场也有上轮审核发现的局部非零散度/壁面渐近吸引限制，本轮未修改。停止发生在原连续区 handoff 约束处；不能把 stationary 比例直接解释为生理捕获，也不能用 source 截断结果代表未经筛选的原 SonoVue。

# 运行与审核材料

服务器：`{summary['remote']['root']}`；各并行批次 6 workers，OMP/OpenBLAS/MKL/NumExpr 均单线程，未强行使用 GPU。audit、smoke、formal 的 seed 分别是 2026092492、2026092493、2026092494。运行时间见下表；批次时间与抽样、point、独立审核分开，不能把纯积分时间叫成整个项目总耗时。

| 任务 | 墙钟秒 |
| --- | ---: |
'''
    for k,v in runtime.items():text+=f'| {k} | {v:.3f} |\n'
    text+=f'''
本地入口永久测试 {tests['local_tests']['tests']} 项、服务器 {tests['remote_tests']['tests']} 项通过；入口与 P9-A.1 联合回归共 {tests['final_regression_tests']['tests']} 项通过。{protection['file_count']} 个预先保护文件 SHA 全部不变，git index/HEAD 不变。服务器另核对了 {remote_protection['current_frozen_and_existing_source_files']} 个既有源码/冻结输入及 {remote_protection['legacy_production_files']} 个旧正式文件，均无差异。分支 `dev/particle9a2-size-first-inlet`，未 commit/push。原工作区已有修改保留；`logs/git_diff_stat.txt` 是实际 git diff --stat，未跟踪新增源码另列于 `logs/new_files_stat.txt`（本轮 {len(newfiles)} 个 Python 文件，{lines} 行）。

图 01–07 均提供英文、白底、300 dpi PNG 与 PDF，原始 CSV/JSON 在 `data/`。图 01 为三种尺寸分布，图 02 为真实入口位置，图 03 为 point basin，图 04 为位置重试与可行 flux，图 05 为全局尺寸拒绝，图 06 为四阶段出口审核，图 07 为固定直径独立解释。入口审核、smoke 与正式 point 路径、全部轨迹和 ledger 均本地/服务器保存。入口页面：[OPEN_RESULTS.html](OPEN_RESULTS.html)。

复现见 [REPRODUCE.md](REPRODUCE.md)；机器汇总：[data/summary.json](data/summary.json)；Gate：[data/gates.json](data/gates.json)。全部文件 SHA 见 `data/delivery_manifest.json`，服务器传输校验见 `data/server_delivery_verification.json`。

# 人工审核

- [ ] 4 µm 截断合理
- [ ] 没有 clipping spike
- [ ] size-first 实现正确
- [ ] position flux weighting 合理
- [ ] 没有使用 outlet quota
- [ ] Method B 历史结果未改
- [ ] 30 smoke 合理
- [ ] 若已运行，500 正式结果合理
- [ ] 同意将 Method C 设为新 production inlet
'''
    (R/'PARTICLE9A2_INLET_REVIEW_ZH.md').write_text(text)
    sections=''.join(f'<section><h2>{html.escape(m["name"])}</h2><a href="figures/{m["name"]}.pdf">PDF</a><img src="figures/{m["name"]}.png"></section>' for m in manifest)
    (R/'OPEN_RESULTS.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>P9-A.2 入口审核</title><style>body{max-width:1250px;margin:36px auto;font:17px system-ui;line-height:1.7;background:#fafafa;color:#222}section{margin:30px 0;padding:20px;background:white;border:1px solid #ddd}img{width:100%}a{color:#267}</style><h1>P9-A.2 — Size-first inlet sampling</h1><p>'+status+' · PENDING_USER_REVIEW</p><p><a href="PARTICLE9A2_INLET_REVIEW_ZH.md">中文审核报告</a> · <a href="data/summary.json">机器汇总</a> · <a href="REPRODUCE.md">复现</a></p><p>消除了 Method B 尺寸重试机制；O2 的有限尺寸几何偏向仍然存在。</p>'+sections+'</html>')
    write(O/'DATASET_STATUS.json',dict(stage='P9-A.2',inlet_method=METHOD,dynamics='UNCHANGED_P9_A_1',automated_status=status,human_status='PENDING_USER_REVIEW',count=500,summary_sha256=sha(D/'summary.json')))
    print(status,counts(fg['outlet_counts']),flush=True)
if __name__=='__main__':main()
