"""Build the two figures and report from executed measurements, without CFD/dynamics."""
from pathlib import Path
import csv
import hashlib
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
IDS=[1189,1329,642,822]


def read(name):return json.loads((HERE/name).read_text())
def save(name,data):(HERE/name).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def arrays(name):
    with np.load(HERE/name) as z:return {k:z[k] for k in z.files}
def norm(v):return float(np.linalg.norm(v))


def key_interface(event):
    options=[]
    for interface in event['interfaces']['measurements']:
        probes=[p for p in interface['approaches'] if p['status']=='MEASURED']
        if interface['interface']['kind']=='TETRA' and len(probes)>=2 and probes[-1]['from_state']['wall_weight']>0:
            options.append((probes[-1]['delta_RHS_m_s'],interface,probes))
    return max(options,key=lambda v:v[0])[1:]


def main():
    events={pid:read(f'result_{pid}.json') for pid in IDS}
    replay=read('REPLAY_RESULTS.json');display=arrays('display_windows.npz')
    plan=read('FROZEN_PLAN.json');metadata=read('SERVER_EVIDENCE.json')
    stages=[]
    for p in sorted((HERE/'monitoring').glob('*/EXECUTION.json')):
        d=json.loads(p.read_text());log=p.parent/'job.log'
        stages.append(dict(stage=p.parent.name,elapsed_s=d['elapsed_s'],exit_code=d['exit_code'],
            scientific_stage_status='BRANCHES_STOPPED_SEE_LOG' if 'STOP_BRANCH' in log.read_text() else 'PROCESS_COMPLETED' if d['exit_code']==0 else 'FAILED_ATTEMPT_SEE_LOG',
            stop_reason=d['stop_reason'],peak_process_tree_GiB=d['peak_tree_RSS_GiB'],
            completed_UTC=d['completed_UTC'],log=log.read_text(),
            receipt_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    stages.sort(key=lambda row:row['completed_UTC'])
    receipt=dict(status='COMPLETED_WITH_EXPLICIT_SCIENTIFIC_LIMITS',stages=stages,
        supervised_compute_total_s=sum(d['elapsed_s'] for d in stages),budget_s=3600,
        retained_closed_replay_count=3,actual_engine_run_count=4,
        engineering_retries='One tuple-unpack adapter correction; one original CPU-fallback adapter correction; one half-dt rerun after an inconclusive global displacement isolation bound. Exact continuous neighbor checks subsequently passed. Failed attempts are included in compute time.',
        original_science_modified=False,full_population_started=False,CFD_started=False,
        tests=read('TEST_RESULTS.json') if (HERE/'TEST_RESULTS.json').exists() else None)
    save('EXECUTION_RECEIPT.json',receipt)
    claims={
        1189:dict(driver='近壁缓慢运动经有限壁面特征离开，进入更快背景流；G/curl 经壁面耦合另叠加小阶跃',
                  confidence='主过程为已复现模型机制；导数界面小阶跃有数值证据；闭合 dt 对照 NOT_EVALUABLE',
                  remedy='补齐可恢复的安全监视器 checkpoint，再仅复跑本窗；单独标记导数界面'),
        1329:dict(driver='另一壁面位置的离壁与背景流增速；局部 G 界面导致小幅减速折点',
                  confidence='主过程为已复现模型机制；局部导数阶跃已定位；闭合 dt 对照 NOT_EVALUABLE',
                  remedy='同窗 checkpoint/dt 对照；保留原场，避免把局部下降折点归为主加速'),
        642:dict(driver='离壁和背景流增速，叠加小 G 界面折点；显示掩码提前再显现强化观感',
                 confidence='模型机制、导数阶跃和显示强化均有证据；闭合 dt 对照 NOT_EVALUABLE',
                 remedy='补 checkpoint；显示说明注明 ±0.05 s 平均与隐藏状态，不能修改轨迹'),
        822:dict(driver='连续 P1 背景 u 的有限空间变化；显式 held 推进及记录粒度影响幅值/路径',
                 confidence='空间增速机制及 dt 敏感性已复现；不是壁面释放；未证明时间步收敛或局部生理精度',
                 remedy='评估将真实过面事件纳入步长控制，并保留此窗作回归；本轮不实施')}
    rows=[]
    for pid,e in events.items():
        before,after=(e['accounting'][0][k] for k in ['from_state','to_state'])
        rows.append(dict(particle_id=pid,time_before_s=before['time_s'],time_after_s=after['time_s'],
            speed_before_mm_s=before['speed_mm_s'],speed_after_mm_s=after['speed_mm_s'],
            gap_before_um=before['gap_um'],gap_after_um=after['gap_um'],
            u_before_mm_s=before['u_speed_mm_s'],u_after_mm_s=after['u_speed_mm_s'],
            direction_change_deg=e['accounting'][0]['direction_change_deg'],
            rhs_reproduction_max_error_m_s=e['maximum_rhs_vs_held_error_m_s'],
            alignment_max_error_m=e['maximum_alignment_error_m'],
            excluded='泡间近场、旧/新协议混淆、动画制造原始运动'+('、壁面阻力/接触' if pid==822 else '、单独硬接触释放解释整体加速'),
            **claims[pid]))
    with (HERE/'event_metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    compact_events={}
    for pid,event in events.items():
        compact_events[str(pid)]={k:v for k,v in event.items() if k!='interfaces'}
        compact_events[str(pid)]['interfaces']={k:v for k,v in event['interfaces'].items() if k not in ['all_segment_partitions','all_tetra_crossings']}
        compact_events[str(pid)]['full_segment_and_interface_evidence_file']=f'result_{pid}.json'
    compact_replay={k:v for k,v in replay.items() if k!='runs'}
    compact_replay['runs']=[{k:v for k,v in run.items() if k!='audit'} for run in replay['runs']]
    compact_replay['all_accepted_and_rejected_audit_rows_file']='REPLAY_RESULTS.json'
    cause=dict(status='DIAGNOSTIC_ONLY',date='2026-10-05',run='realization_0000_contactfix',
        flow_sha256='73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc',
        source_historical_commit='98d705f13babdfa0cba6c12331623292a588cb77',
        original_field_FAIL_preserved=True,new_fixed_window_population_run=False,
        event_conclusions=rows,events=compact_events,replay=compact_replay,
        original_ledgers={str(k):read(f'ledger_{k}.json') for k in IDS},
        display=read('DISPLAY_CONTRACT_CHECK.json'),input_integrity=read('FINAL_INPUT_INTEGRITY.json'),
        read_historical_evidence=plan['existing_evidence'],
        two_factor_mixed_packet_test=dict(status='NOT_RUN',reason='Assembly decomposition and strict spatial controls identify the entering terms; no mixed physical packet required'),
        unexecuted=['Closed re-integration of1189/1329/642: missing historical safety-monitor state',
                    'Full population or CFD', 'Physiological validation', 'Formal timestep convergence'],
        execution_receipt='EXECUTION_RECEIPT.json')
    save('CAUSE_RESULTS.json',cause)

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(4,2,figsize=(14,12),layout='constrained')
    for row,pid in enumerate(IDS):
        z=arrays(f'measured_{pid}.npz');t=z['time_s'];a,b=axes[row]
        a.step(t,np.linalg.norm(z['held_m_s'],axis=1)*1e3,where='post',color='#1675ba',lw=1.2,label='Original held bubble speed')
        a.plot(t,np.linalg.norm(z['flow'][:,:3],axis=1)*1e3,color='#d58116',ls='--',lw=1.1,label='FEM u at accepted positions')
        if pid!=822:a.set_yscale('log')
        a.set(title=f'ID {pid}: original accepted trajectory',ylabel='Speed (mm/s)',xlabel='Physical time (s)')
        if pid==642:
            ft=display['time_642'];means=display['rolling_speed_642']*1e3;shown=display['displayed_642']
            a.plot(ft,means,'.',color='#8441a7',ms=3,label='Display rolling mean (±0.05 s)')
            a.scatter(ft[~shown],means[~shown],marker='x',s=12,color='#ae2145',label='Hidden display frames')
            a.axvline(1.025287356321839,color='#8441a7',ls=':',lw=1)
        if pid!=822:
            interface,probes=key_interface(events[pid]);a.axvline(interface['located_time_s'],color='#b22934',ls=':',lw=1,label='Measured derivative interface')
        else:a.axvline(.651,color='#b22934',ls=':',lw=1,label='Largest original increment')
        a.legend(fontsize=7.5,loc='best');a.grid(alpha=.2)
        b.plot(t,z['wall'][:,0]*1e6,color='#1675ba',lw=1.2,label='Actual sphere-surface gap')
        b.axhline(events[pid]['radius_m']*1e6,color='#444444',ls='--',label='One bubble radius')
        if pid!=822:b.axhline(.002,color='#999999',ls=':',label='2 nm handoff');b.set_yscale('log')
        for interface in events[pid]['interfaces']['measurements']:
            if interface['interface']['kind']!='TETRA':b.axvline(interface['located_time_s'],color='#9b4f80',alpha=.35,lw=.7)
        b.set(title='Gap; purple marks: sampled wall/contact switches',xlabel='Physical time (s)',ylabel='Surface gap (µm)')
        b.grid(alpha=.2);b.legend(fontsize=8)
    fig.suptitle('Frozen BraVa field / historical population — measurements, not physical validation\nNo smoothing. Held-speed steps use actual accepted times; FEM lines join recorded samples.',fontsize=13)
    fig.savefig(HERE/'original_events.png',dpi=300);plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    colors=['#1675ba','#d58116','#29915b']
    for result,color in zip(replay['runs'],colors):
        dt=result['dt_upper_s'];z=arrays(f'replay_822_dt{dt*1e6:.0f}us.npz');t=z['time_s']
        mask=(t>=.6503)&(t<=.6525)
        axes[0,0].plot(t[mask],np.linalg.norm(z['RHS_m_s'][mask],axis=1)*1e3,'.-',ms=4,color=color,label=f'{dt*1e3:g} ms cap')
        if dt<.0005:
            common=z['common_indices'];axes[0,1].plot(t[common],z['path_delta_at_common_m']*1e6,color=color,label=f'{dt*1e3:g} ms vs original')
    axes[0,0].axvspan(.6505,.651,color='gray',alpha=.12)
    axes[0,0].set(title='ID822: RHS at actual accepted positions',xlabel='Physical time (s)',ylabel='Speed (mm/s)')
    axes[0,1].set(title='ID822: path differences at common saved times',xlabel='Physical time (s)',ylabel='Path difference (µm)')
    for pid,color in zip(IDS[:3],colors):
        interface,probes=key_interface(events[pid])
        axes[1,0].loglog([p['separation_m']*1e9 for p in probes],[p['delta_RHS_m_s']*1e3 for p in probes],'o-',color=color,label=f'ID {pid}, t={interface["located_time_s"]:.6f} s')
    axes[1,0].set(title='Near-wall P1 derivative interfaces: residual velocity difference',xlabel='Distance between valid probes (nm)',ylabel='Norm of translational RHS difference (mm/s)')
    projected=[]
    for pid in IDS[:3]:
        a=events[pid]['pathwise_increment_accounting'];parts=np.array(a['summed_translation_mm_s']);dv=np.array(a['actual_translation_delta_mm_s']);unit=dv/np.linalg.norm(dv)
        groups=np.array([parts[0],parts[1]+parts[2],parts[3]+parts[4],parts[5]])
        projected.append(groups@unit)
    matrix=np.array(projected);limit=np.max(abs(matrix));ax=axes[1,1]
    ax.imshow(matrix,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto')
    ax.set_xticks(range(4),['Background u','Curl + shear','Wall resistance','Constraint'])
    ax.set_yticks(range(3),[str(p) for p in IDS[:3]])
    for i in range(3):
        for j in range(4):ax.text(j,i,f'{matrix[i,j]:+.3f}',ha='center',va='center',color='white' if abs(matrix[i,j])>.6*limit else 'black')
    ax.set(title='Pathwise contributions projected on net Δv (mm/s)\nSigned algebra, not causal percentages or speed differences',ylabel='Particle ID')
    for ax in [axes[0,0],axes[0,1],axes[1,0]]:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Controls: original GPU engine, frozen external event clock, exact equation accounting\nLines connect measured states; they do not represent additional integration.',fontsize=13)
    fig.savefig(HERE/'cause_controls.png',dpi=300);plt.close(fig)

    report=['# 原始速度急变的原因定位\n\n2026-10-05；原冻结 BraVa 场与原正式历史群体；**DIAGNOSTIC_ONLY**。',
        '## 事件结论（先给答案）\n',
        '下表速度是在所列真实接受位置计算、且与下一行保存的 held velocity 对齐的瞬时速率，不是视频窗口平均。gap 是球面到有限三角壁面的间隙。',
        '| ID | 物理时刻 s | 微泡速度 mm/s | gap µm | 背景 u mm/s | 直接驱动及关键对照 | 排除项 | 可信度／最小处理建议 |',
        '|---|---|---|---|---|---|---|---|']
    for row in rows:
        pid=row['particle_id']
        control='原 GPU 原/半/四分之一步长；逐面双侧查询' if pid==822 else '原方程逐点重现＋分项增量＋严格界面双侧查询'
        report.append(f"| {pid} | {row['time_before_s']:.7g} → {row['time_after_s']:.7g} | {row['speed_before_mm_s']:.6g} → {row['speed_after_mm_s']:.6g} | {row['gap_before_um']:.6g} → {row['gap_after_um']:.6g} | {row['u_before_mm_s']:.6g} → {row['u_after_mm_s']:.6g} | {row['driver']}；{control} | {row['excluded']} | {row['confidence']}。{row['remedy']} |")
    report.extend(['',
        '**主结论：**1189、1329、642 的大幅慢→快过程，是当前模型在当前原始输入下的离壁及进入更快背景流过程；另有小幅、可定位到 P1 导数界面的离散速度阶跃。822 属于高速内部流动中的空间增速，单步幅度和整窗路径对 dt 敏感。动画会强化突然感，不能解释原始轨迹为何变化。以上均不等于真实生理运动已验证。',
        '## 身份、已有证据与本轮新增量',
        '服务器原始 `SUMMARY.json`、`largest_step_increments.json`、`largest_increment_context.json`、`display_review.json` 已实际读取，未重复扫描1565条或运行旧审计。原运行 `realization_0000_contactfix`：1565条，出口458/527/437、141观察截止、2数值失败；不是固定灌注1.240815660 s＋观察10 s的新运行。',
        '场 SHA256：`73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`。历史代码提交 `98d705f13babdfa0cba6c12331623292a588cb77`。142份原科学源码按该提交核对，仅 `population_v1/run.py` 是当前已变更的入口；本轮不导入它，也不调用新 fixed-window 入口。实际输入共195项通过结束后复核，见 `INPUT_SHA256.json`。原 field FAIL 和正式状态均保留。',
        '| 已有结论 | 实际证据文件 | 本轮补充／仍缺对照 |','|---|---|---|'])
    for e in plan['existing_evidence']:report.append(f"| {e['conclusion']} | `../{e['evidence']}` | {e['missing_control']} |")
    report.extend(['',
        '旧201候选是前后各约0.06908 s的视频窗口筛选，不是瞬时跳变定义；“增量>5 mm/s且速率翻倍”零命中不等于没有急变。636/1189/393 属同一空间区域，未重复它们的12个旧固定状态作为新成果。440/1505 的全零单元停滞结论保持，未作为本轮增速案例。',
        '## 时间语义、邻居与真实线段',
        '原 `gpu_engine._trial` 先在 `(x[i],t[i])` 解速度，再用实际 dt 更新位置，写到第 i+1 行；出生行存初始自由速度。本轮逐段检验 `x[i+1]−x[i]=dt[i]·velocity[i+1]`，将出生、事件切步及非名义步长单列。四窗最大位移恒等式误差为 '+f"{max(e['maximum_alignment_error_m'] for e in events.values()):.3g} m，最大 RHS 重现误差为 {max(e['maximum_rhs_vs_held_error_m_s'] for e in events.values()):.3g} m/s。",
        '窗口均冻结为中心±0.12 s，没有扩大。822在0.650308305496727 s出生，窗口中此前无粒子，未伪造出生前状态。原日志及每个共同接受线段上的连续相对运动检查均无泡间近场；四窗最小近场外余量分别约22.071、523.491、74.043、85.417 µm。全部真实接受线段均枚举了穿过的单元，包括中间单元和最后一个区间。',
        '| ID | 实际 dt 最小 / 中位 / 最大 µs | 全部线段过面数 | 原拒绝试步数（其中本粒子负责） |','|---|---|---|'])
    for pid,e in events.items():
        dt=e['actual_dt_s'];report.append(f"| {pid} | {dt['min']*1e6:.5g} / {dt['median']*1e6:.5g} / {dt['max']*1e6:.5g} | {len(e['interfaces']['all_tetra_crossings'])} | {e['original_rejected_trials_in_window']}（{e['original_rejects_responsible_target']}） |")
    report.extend(['',
        '原窗口没有目标粒子的 CPU 组件回退；界面诊断的个别探针触发原有 CPU 接触解算分支，记录器直接复用该权威函数。未删除实际邻居、投影位置、更改 gap、滤波速度或放宽原可接受性条件。',
        '## 哪一项造成变化：方程增量与界面控制',
        '`R U = b + q + r`，`q=Jᵀλ` 是运动学乘子修正，`r=RU−b−q`。原解算器符号已核对。逐步核算 `R₁ΔU=Δb−ΔR·U₀+Δq+Δr`，用原线性解算函数求响应，没有构造显式逆。',
        '`b` 分为自由背景平动、背景旋转（同一G的curl/2）以及壁面剪切仿射项；阻力分为法向近场和切向平动/旋转壁面块。完整平动 m/s、旋转 s⁻¹、力 N、力矩 N·m 分开保存。逐步分量以“各相邻状态的后状态R”为参考，存在抵消，不能当作唯一因果占比。',
        '下表是前后代表状态之间，逐步分量在净速度向量变化方向上的投影，单位 mm/s。它们相加对应速度**向量差**的模，不是标量速率差；舍入后可能略有差异。',
        '| ID | 背景平动u | curl＋壁面剪切 | 壁面阻力变化 | 实际约束修正 |','|---|---|---|---|'])
    for pid,vals in zip(IDS[:3],matrix):report.append('| '+str(pid)+' | '+' | '.join(f'{v:+.6f}' for v in vals)+' |')
    report.extend(['',
        '大幅后段增速主要由沿路径采样到更快的背景 u 直接驱动，壁面阻力变化帮助粒子从受抑制运动转向随流；不能将整个过程叫作硬接触弹开。1189旧固定状态取消硬约束而保留壁面阻力仅0.03030→0.03035 mm/s，本轮逐步约束项也无法解释最终约49 mm/s速度。',
        '### 已定位的小幅导数界面阶跃',
        '以下控制均沿真实接受线段，从界面两侧以递减间距查询；每个探针由原包含判据及严格非负重心坐标双重验证，没有向邻单元外推。',
        '| ID / 界面时刻 s | 最细有效探针间距 nm | 微泡左右速率 mm/s | 速度向量差 mm/s | 背景u差 mm/s | G差 s⁻¹ |','|---|---|---|---|---|---|'])
    for pid in IDS[:3]:
        interface,probes=key_interface(events[pid]);p=probes[-1]
        report.append(f"| {pid} / {interface['located_time_s']:.12g} | {p['separation_m']*1e9:.6g} | {p['from_state']['speed_mm_s']:.8g} → {p['to_state']['speed_mm_s']:.8g} | {p['delta_RHS_m_s']*1e3:.7g} | {p['delta_u_m_s']*1e3:.5g} | {p['delta_G_s_inv']:.6g} |")
    report.extend(['',
        '这些界面上，u差随间距缩小，G差保持，壁面法向不变或趋近，微泡速度差保留。分解把保留项定位到 **curl引起的背景旋转载荷及壁面剪切载荷**；平动u、阻力变化和约束修正趋小。1189的小阶跃约0.03188 mm/s，发生在主离壁释放之前；1329所列阶跃是小幅下降，不能把它称为整体加速原因；642约0.00832 mm/s的向量阶跃同样远小于其后段数十mm/s增速。',
        '这构成“离散导数进入粒子RHS导致局部数值阶跃”的直接证据，强于仅报告G不连续或非零散度；仍不能由它推导整段轨迹的生理误差。所测试的高速体内界面上，u与粒子RHS差随距离递减，即使G差仍保留，也没有对应的有限平动速度跃阶。',
        'FACE→EDGE附近实测法向及RHS两侧差很小，部分可解析壁面切换差随间距缩小。极接近2 nm接触壳的一些切换间距已低于严格定位余量，机器结果保留 NOT_EVALUABLE；不声称验证了所有几何/约束界面的解析极限，也不把三角形ID变化本身当根因。',
        '## 822：原引擎闭合回放与时间步影响',
        '原0.6505→0.651 s区间位移168.6564 µm，穿越3个内部面、4个四面体；场u的空间变化解释其速度向量增量 `[11.344734,35.863133,3.435389]` mm/s。其gap约236 µm，壁面权重及约束均为零。三个过面位置的两侧u/RHS差随间距下降，没有测得有限平动速度阶跃。',
        '从实际出生的x、v、omega、q、半径、出生时间和源事件开始，调用原 `GPUPopulationEngine`，保留原事件处理及拒绝规则。其它粒子的原接受时钟和路径被冻结，因此这是孤立组件的闭合局部对照，不是全群体dt认证。原步长轨迹最大位置重现差4.66×10⁻¹⁷ m。',
        '| dt上限 ms | 实际接受区间数 | 0.6505 / 0.651 s速度 mm/s | 同物理区间增速 mm/s | 附近最大相邻增速 mm/s | 全窗对原路径最大差 µm |','|---|---|---|---|---|---|'])
    for run,event in zip(replay['runs'],replay['event_comparison']):
        s=event['common_speed_mm_s'];report.append(f"| {run['dt_upper_s']*1e3:g} | {run['actual_dt_s']['accepted_count']} | {s[0]:.6f} / {s[1]:.6f} | {event['speed_increment_over_original_physical_interval_mm_s']:.6f} | {event['maximum_adjacent_scalar_increment_mm_s']:.6f} | {run['maximum_path_delta_at_common_m']*1e6:.6g} |")
    report.extend(['',
        '半步与四分之一步的路径最大差约200.302 µm，均在实际共同接受时刻比较；相应真实接受步确已变小。原路径终点横截面（原事件末位置的y平面）到达时刻为0.651、0.6509920884、0.6509874338 s。四分之一步将固定速率门槛的10%–90%上升历时限定在0.125–0.375 ms，仍是接受点括区间，不插值伪造积分步。',
        '结果支持有限空间增速持续存在，也支持“单条0.5 ms记录使变化显得更集中”；原0.5 ms显式推进对幅度、时机及长一些的窗内路径存在明显数值偏差。没有据三步长宣布收敛。所有其它历史粒子的连续相对路径均复核：三个步长下近场外最小余量约85.42、83.76、82.42 µm，故没有漏掉会改变该闭合对照的泡间组件。',
        '### 其它三个事件为何没有闭合dt结论',
        '原始状态块保存了x/v/omega/q、半径和时钟，但未保存中途恢复所需的GPUSafetyMonitor provider调用计数、no-acceptance计数及名义步/历史环形缓冲。它们会影响是否提前终止与接受性。重置为零会改变原引擎状态。本轮不伪造它们，也不从出生重算数秒轨迹来扩大窗口。',
        '因此1189/1329/642的RHS、逐步驱动项和导数阶跃可评价，但“离壁时刻/路径对dt的敏感性”仍不能与“该有限时间过程已时间准确”区分。最小缺失对照是包含这些字段和原外部事件时钟的可恢复checkpoint，再进行同一固定窗口原/半步长对照；无需先换CFD或换框架。',
        '## 显示只强化观感',
        '旧576帧动画的物理帧间距约0.0230275 s（约46个名义步）；当前0–1.2 s动画262帧，物理间隔约0.0045977 s。后者按t±0.05 s真实3D速率的时间加权均值≥2 mm/s筛选；出生/终止只取实际寿命交集。此均值有0.05 s前视，不能当成当前RHS。',
        '642在1.0252873563 s再显现，而其主要FACE→EDGE释放约1.0463330998 s，显示门槛提前约21 ms被未来快段抬高。原轨迹此前一直存在；没有粒子重新出生或动力学跳跃。1189和1329的所选事件不在当前0–1.2 s播放范围。新读取的显示数据与原始速度积分、阈值掩码逐点复核通过。',
        '## 最小建议与保留边界',
        '1. 输出与诊断明确保存RHS求值时刻、区间起止及held速度；保留出生速度独立语义，避免错位分析。\n2. 为原引擎增加完整可恢复checkpoint；之后只补三个已有窗口的dt对照。\n3. 对822评估把真实过面事件和局部位移误差纳入时间步控制，用本窗三步长结果作回归。\n4. 对近壁G/curl/剪切项保留界面回归；后续若修复派生场，必须保持(u,G,curl)一致并单独验证，不以速度平滑或动画美观作为通过标准。\n5. 显示说明区分低速隐藏、出生、出流与观察截止；不改原轨迹和账本。',
        '这些是建议，本轮未实施任何生产物理修复。没有混合场包的积分，没有新CFD、网格、求解器、全域审计或新population。',
        '## 文件、测试与执行记录',
        '[原轨迹综合图](original_events.png)；[界面/步长/增量控制图](cause_controls.png)；[机器结果](CAUSE_RESULTS.json)；[事件指标CSV](event_metrics.csv)；[执行receipt](EXECUTION_RECEIPT.json)；[复现方法](REPRODUCE.md)。',
        f"监督计算累计 {receipt['supervised_compute_total_s']:.2f} s，包含未完成的工程尝试，低于60 min；每个窗口及各个新增作业均低于10 min。保留3条不同dt的822回放，连同首次半步运行共实际执行4条，未超过10条上限。没有后台科学任务继续运行。",
        '测试仅针对本轮新增的时间对齐、增量重构和严格接口采样，以及已保存结果的闭合性；测试通过不代表生产物理资格通过。`TEST_RESULTS.json` 保存实际命令与结果。'])
    document=re.sub(r'\|\n\n\|','|\n|','\n\n'.join(report))
    document=re.sub(r'\n{3,}','\n\n',document)
    (HERE/'CAUSE_REPORT_ZH.md').write_text(document+'\n')
    print('BUILT_REPORT_AND_TWO_PNGS',receipt['supervised_compute_total_s'],flush=True)


if __name__=='__main__':main()
