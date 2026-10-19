from pathlib import Path
import json,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/workspace/y_clot/project');TASK='clot_irregular_lps_discretization_v4';RUN='20261019T221000Z'
O=R/'verification'/TASK/RUN;A=R/'analysis'/TASK;F=R/'results'/TASK/RUN
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
    (F/'figures').mkdir(parents=True,exist_ok=True)
    baseline={};comparisons={}
    for window in ['intact','damaged']:
        before=json.loads((O/'p0'/f'PERFORMANCE_{window.upper()}.json').read_text());after=json.loads((O/'p1'/f'PERFORMANCE_{window.upper()}.json').read_text());baseline[window]=before
        b=O/'p0/windows'/window/'biphasic';a=O/'p1/windows'/window/'biphasic';checks={}
        for name in ['initial.npz','final.npz']:
            x=dict(np.load(b/'keys'/name));y=dict(np.load(a/'keys'/name))
            for key in ['x','velocity','force','material_force','kappa','d','atom_ids','external_force','surface_force','surface_active']:
                checks[name+':'+key]=bool(np.array_equal(x[key],y[key]));assert checks[name+':'+key]
            assert np.array_equal(x['scalars'][:6],y['scalars'][:6])
        assert (b/'history.csv').read_bytes()==(a/'history.csv').read_bytes()
        assert before['summary']['secondary_events']==after['summary']['secondary_events']
        comparisons[window]={'before_wall_s':before['wall_s'],'after_wall_s':after['wall_s'],'speedup':before['wall_s']/after['wall_s'],'before_cpu_s':before['cpu_s'],'after_cpu_s':after['cpu_s'],'before_peak_RSS_bytes':before['peak_RSS_bytes'],'after_peak_RSS_bytes':after['peak_RSS_bytes'],'accepted_before_after':[before['summary']['accepted_steps'],after['summary']['accepted_steps']],'prediction_rejects_before_after':[before['summary']['predictive_rejections'],after['summary']['predictive_rejections']],'state_bitwise_checks':checks,'all_history_csv_bytes_identical':True,'lineage_parent_children_identical':True,'exclusive_stage_savings_s':{k:v['wall_s']-after['python_exclusive'][k]['wall_s'] for k,v in before['python_exclusive'].items()}}
    dump(O/'PERFORMANCE_BASELINE_V3.json',{'windows':baseline,'prior_full_run_cost':json.loads((R/'verification/clot_fragment_surface_coupling_v3/20261018T170400Z/COMPUTE_COST.json').read_text()),'short_windows_not_full_run_extrapolation':True,'rollback_note':'V3 rejects predicted drift before native stepping. No committed-state rollback occurs in these windows; immutable input and geometry are not repeatedly copied by a rollback operation.','native_timings_are_subdivisions_not_added_to_python_totals':True})
    dump(O/'PERFORMANCE_OPTIMIZATION.json',{'status':'PASS','comparisons':comparisons,'changes':['Direct old-member to current-label lineage lookup instead of all component-pair intersections','Exact reference query cache','Exact current-tetra containment cache keyed by coordinates and solid support'],'timing_scope':'One serial paired eight-step sample per state. Includes initialization and dense first-event outputs. No claimed full-protocol speedup.','physics_changed':False})
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(15,6),layout='constrained')
    for ax,(window,q) in zip(axes,baseline.items()):
        rows=sorted(q['python_exclusive'].items(),key=lambda z:z[1]['wall_s'],reverse=True)
        labels=[k.replace('_',' ') for k,_ in rows];values=[v['wall_s'] for _,v in rows]
        ax.barh(labels,values,color='#9c2338');ax.invert_yaxis();ax.set_title(f'{window.title()}: 8 accepted steps');ax.set_xlabel('Exclusive wall time (s)')
    fig.savefig(F/'figures/FIG_01_NATIVE_PERFORMANCE_BREAKDOWN.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,(window,q) in zip(axes,comparisons.items()):
        vals=[q['before_wall_s'],q['after_wall_s']];ax.bar(['V3 original','Exact optimization'],vals,color=['#8499a6','#9c2338']);ax.set_title(f'{window.title()} window: {q["speedup"]:.2f}×');ax.set_ylabel('Total window wall time (s)')
        for i,v in enumerate(vals):ax.text(i,v+.4,f'{v:.2f}',ha='center')
    fig.savefig(F/'figures/FIG_02_PERFORMANCE_OPTIMIZATION.png',dpi=180);plt.close(fig)
    shutil.copy2(F/'figures/FIG_01_NATIVE_PERFORMANCE_BREAKDOWN.png',F/'figures/FIG_01_NATIVE_TIME_BREAKDOWN.png');shutil.copy2(F/'figures/FIG_02_PERFORMANCE_OPTIMIZATION.png',F/'figures/FIG_02_BEFORE_AFTER_TIMING.png')
    q=comparisons['damaged']
    (O/'PERFORMANCE_BOTTLENECKS.md').write_text(f'''# V3 实测性能与保持物理语义的优化\n\n读取正式 2.650125 ms 损伤态并继续八个原生接受步，未重跑历史完整工况。原窗口 {q['before_wall_s']:.6f} s，优化窗口 {q['after_wall_s']:.6f} s，约 {q['speedup']:.3f} 倍。原始谱系采用父分量与全部当前分量交叉匹配；现在按原始成员 ID 直接查询当前标签，保留每一次真实分裂与完整成员定义。初末态物理数组、完整标量 CSV 与父子事件完全一致。\n\n完整态窗口未测得加速：{comparisons['intact']['before_wall_s']:.6f}→{comparisons['intact']['after_wall_s']:.6f} s。各窗口包含初始化、首次事件集中保存及末态 restart；不能把短窗口加速比当作原 11.7 小时工况的已实测整体加速。计时为单次配对样本，没有统计误差区间。逐项 wall/CPU/次数/均值及原生子项见 JSON。\n\n原生力/能量计算并非此损伤窗口的最大热点。表面几何始终按当前坐标和真实失效拓扑更新；缓存严格比较输入，没有每 N 步延迟受力更新。预测拒绝发生在原生提交之前，因此本窗口没有已提交状态的 rollback。物理短窗只做两组前后配对，不新增参数工况。\n''')
    print('PERFORMANCE_EQUIVALENCE_PASS',comparisons['damaged']['speedup'],flush=True)
if __name__=='__main__':main()
