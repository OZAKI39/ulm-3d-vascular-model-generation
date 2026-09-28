"""Export accepted-step data, regular 0.5-ms position samples, CSV and report."""
from pathlib import Path
import csv, hashlib, json
import numpy as np

HERE=Path(__file__).resolve().parents[1]


def main():
    complete=json.loads((HERE/'data/production_complete.json').read_text())
    metrics=json.loads((HERE/'data/metrics.json').read_text())
    cohort=json.loads((HERE/'data/cohort.json').read_text())
    prep=json.loads((HERE/'data/preparation.json').read_text())
    gpu=json.loads((HERE/'data/gpu_tracks_validation.json').read_text())
    parity=json.loads((HERE/'data/native_kernel_parity.json').read_text())
    exact=json.loads((HERE/'data/native_exact_bytes_verification.json').read_text())
    assert complete['integration_finished'] and len(metrics)==1500
    assert complete['numerical_gate'] and gpu['PASS'] and exact['PASS_FOR_ACTUAL_HANDOFF_DEPENDENCY']
    events={e['particle_id']:e for e in cohort['events']}
    columns=['particle_id','source_event_id','diameter_um','birth_time_s','status','outlet',
             'trajectory_age_s','sample_count','accepted_steps','rejected_trials',
             'path_length_m','maximum_speed_m_s','minimum_wall_gap_m','minimum_g_nf_m',
             'stationary_supported','penetration_count','handoff_violation_count',
             'inlet_escape_count','nan_inf_count','unclassified_corruption_count']
    with (HERE/'data/trajectory_catalog.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        writer.writerows({k:r[k] for k in columns} for r in metrics)
    counts=[]
    for outlet in ['O1','O2','O3']:
        subset=[r for r in metrics if r['outlet']==outlet]
        counts.append(dict(outlet=outlet,trajectory_count=len(subset),fraction_of_1500=len(subset)/1500,
            fraction_of_exited=len(subset)/max(sum(r['outlet'] is not None for r in metrics),1),
            fluid_fraction=prep['flux']['OUTLET_0'+outlet[1]]['signed_Q_m3_s']/prep['flux']['INLET']['signed_Q_m3_s']))
    with (HERE/'data/outlet_summary.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(counts[0]));w.writeheader();w.writerows(counts)
    samples=[];offsets=[0];cpu=0.;maxrss=0
    for row in metrics:
        folder=HERE/'tracks'/f"mb_{row['particle_id']:06d}"
        s=np.load(folder/'trajectory.npz')['samples'];event=events[row['particle_id']]
        # Exact linear accepted-segment position interpolation; terminal point
        # retained at its actual crossing time even when not on the 0.5 ms grid.
        times=np.arange(int(np.floor(s[-1,0]/.0005))+1)*.0005
        if abs(times[-1]-s[-1,0])>1e-12:times=np.r_[times,s[-1,0]]
        else:times[-1]=s[-1,0]
        xyz=np.column_stack([np.interp(times,s[:,0],s[:,k]) for k in range(1,4)])
        samples.append(np.column_stack((times,times+event['birth_time_s'],xyz)))
        offsets.append(offsets[-1]+len(times))
        receipt=json.loads((folder/'receipt.json').read_text())
        cpu+=receipt['cpu_seconds'];maxrss=max(maxrss,receipt['peak_rss_kib'])
    np.savez_compressed(HERE/'data/trajectories_dt0p5ms.npz',samples=np.vstack(samples),
        offsets=np.asarray(offsets),particle_ids=np.array([r['particle_id'] for r in metrics]),
        radius_m=np.array([r['radius_m'] for r in metrics]),
        columns=np.array(['age_s','physical_time_s','x_m','y_m','z_m']))
    summary=dict(complete,CPU_seconds=cpu,CPU_core_equivalents=cpu/complete['elapsed_s'],
        peak_worker_rss_kib=maxrss,outlet_summary=counts,GPU_diagnostics_pass=gpu['PASS'],
        all_outlets_observed=all(r['trajectory_count']>0 for r in counts),native_geometry_parity=parity,
        regular_export_samples=offsets[-1],birth_last_s=cohort['events'][-1]['birth_time_s'])
    outlet_table='\n'.join(f"| {r['outlet']} | {r['trajectory_count']} | {100*r['fraction_of_1500']:.3f}% | {100*r['fraction_of_exited']:.3f}% | {100*r['fluid_fraction']:.3f}% |" for r in counts)
    (HERE/'data/final_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    text=f"""# BraVa 四端口新流场微泡轨迹

使用同一经过审核的 18 mL/min 稳态 FEM 流场。原始有限尺寸微泡模块保持不变。

- 总数：1500；名义 dt=0.5 ms，近壁事件可能采用内部子步。
- 出口计数：{complete['outlets']}；终态：{complete['statuses']}。
- 三个出口实际均有微泡：{summary['all_outlets_observed']}。
- 安全检查计数：{complete['safety_totals']}；来源、完成标记和文件哈希：{complete['all_completion_hashes_verified']}。
- CPU 进程并行积分；CUDA FP64 检查和 NVIDIA EGL 绘图；不将 GPU 后处理称为 GPU 微泡动力学。
- 出口覆盖取自然模拟结果；不按出口筛选或改变 SonoVue ≤4 μm 原始条件粒径分布。
- 单向冻结流场、独立微泡，无微泡–微泡及 RBC 耦合。
- 1000/s 是有限数值样本的到达率设定，不代表已测得的人体剂量。动画按每条轨迹的年龄对齐，不代表同步注射。
- candidate 0、15 共用轨迹，只对显示位置作已核实的刚性变换。

| 出口 | 微泡数 | 全部样本占比 | 完成者占比 | 流体分流占比 |
|---|---:|---:|---:|---:|
{outlet_table}

轨迹计数含有限采样、有限尺寸及离散积分效应，不能要求等于流量占比。
本轮没有证明时间步或网格无关性。`SUPPORTED_STATIONARY` 与停留时间截断分别报告，均不自动等同于生理滞留。
内部截面与独立体积散度积分已确认，此冻结P1场存在约百分之几的局部速度通量缺陷；边界总流量守恒和轨迹安全检查不能消除这一限制。近壁停留、路径分配及局部WSS不应解读为已充分验证的物理预测。详细数据见 `../reports/internal_sections/`。
"""
    (HERE/'MICROBUBBLE_RESULTS_ZH.md').write_text(text)
    print(json.dumps({k:summary[k] for k in ['outlets','statuses','numerical_gate','CPU_core_equivalents','all_outlets_observed']},indent=2))
if __name__=='__main__':main()
