"""Source-audited residual plots: global evolution and within-step reduction.

Uses the first KSP history norm associated with each logged NS iteration,
exactly as petsc_impl.cpp supplies FSILS.RI.iNorm to Integrator::corrector.
No field solve, trajectory integration or video rendering is performed.
"""
from pathlib import Path
import json
import re
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from compute_field_diagnostics import ROOT, CASE, OUT, dump, sha, write_csv, parse_solver_log

ORANGE='#b97422'
BLUE='#19496b'
GRAY='#6c757e'


def nonlinear_stop(iteration, global_ratio, step_ratio, tolerance, minimum, maximum):
    """The exact OR rule in Integrator.cpp, including iteration limits."""
    return (iteration >= maximum or
            ((global_ratio <= tolerance or step_ratio <= tolerance) and iteration >= minimum))


def orders_reduced(start, end):
    if start <= 0 or end <= 0:
        raise ValueError('Residual norms must be positive')
    return float(np.log10(start/end))


def audit_nonlinear(history):
    tree=ET.parse(CASE/'run/solver.xml')
    eq=tree.find(".//Add_equation[@type='fluid']")
    tol=float(eq.find('Tolerance').text)
    minimum=int(eq.find('Min_iterations').text)
    maximum=int(eq.find('Max_iterations').text)
    assert tree.find('.//Continue_previous_simulation').text.lower()=='false'
    log=(CASE/'run/solver.log').read_text()
    events={}
    for line_number,line in enumerate(log.splitlines(),1):
        if line.startswith('SV13Q_BEGIN '):
            fields=dict(re.findall(r'(\w+)=([^\s]+)',line))
            assert fields['attempt']=='0', 'Recovery attempts require a separate mapping'
            fields['log_line']=line_number
            events[int(fields['logical'])]=fields
    solves=history['linear_solves']
    reference=solves[0]['petsc_monitor'][0]['residual_norm']
    assert solves[0]['step']==1 and solves[0]['nonlinear_iteration']==1
    groups={}
    iteration_rows=[]
    max_print_error=0.
    for solve in solves:
        first_monitor=solve['petsc_monitor'][0]
        assert first_monitor['iteration']==0
        norm=first_monitor['residual_norm']
        assert np.isclose(norm,first_monitor['true_residual_norm'],rtol=2e-12)
        groups.setdefault(solve['step'],[]).append((solve,norm))
    rows=[]
    for step,group in groups.items():
        initial=group[0][1]
        assert [s['nonlinear_iteration'] for s,n in group]==list(range(1,len(group)+1))
        for i,(solve,norm) in enumerate(group):
            global_ratio=norm/reference
            step_ratio=norm/initial
            for reconstructed,printed in [(global_ratio,solve['nonlinear_Ri_over_R0']),
                                          (step_ratio,solve['nonlinear_Ri_over_R1'])]:
                # NS text prints four significant figures; monitor prints more.
                assert np.isclose(reconstructed,printed,rtol=5.1e-4,atol=0)
                max_print_error=max(max_print_error,abs(reconstructed/printed-1))
            stop=nonlinear_stop(i+1,global_ratio,step_ratio,tol,minimum,maximum)
            assert stop==(i==len(group)-1), f'Stopping-rule mismatch at step {step}, iteration {i+1}'
            event=events[solve['linear_solve_index']]
            assert int(event['step'])==step and int(event['log_line'])<solve['log_line']
            assert solve['linear_converged'] and not solve['petsc_reason']['diverged']
            iteration_rows.append(dict(step=step,iteration=i+1,log_line=solve['log_line'],
                norm=norm,global_ratio=global_ratio,step_ratio=step_ratio,reported_global_ratio=solve['nonlinear_Ri_over_R0'],
                reported_step_ratio=solve['nonlinear_Ri_over_R1'],minimum_iterations_met=(i+1>=minimum),
                global_tolerance_met=(global_ratio<=tol),step_relative_tolerance_met=(step_ratio<=tol),
                maximum_iterations_reached=(i+1>=maximum),stop_rule_satisfied=stop,
                linear_reason=solve['petsc_reason']['reason'],preconditioner_reused=event['reuse']=='1',
                preconditioner_rebuild_reason=event['rebuild_reason'],preconditioner_log_line=event['log_line']))
        last=group[-1][1]
        final=iteration_rows[-1]
        row=dict(step=step,time_s=group[-1][0]['time_s'],iterations=len(group),start_norm=initial,end_norm=last,
            start_global=initial/reference,end_global=last/reference,end_over_start=last/initial,
            reduction_orders=orders_reduced(initial,last),
            logged_start_global=group[0][0]['nonlinear_Ri_over_R0'],
            logged_end_global=group[-1][0]['nonlinear_Ri_over_R0'],
            global_tolerance_met=final['global_tolerance_met'],step_relative_tolerance_met=final['step_relative_tolerance_met'],
            minimum_iterations_met=final['minimum_iterations_met'],maximum_iterations_reached=final['maximum_iterations_reached'],
            first_log_line=group[0][0]['log_line'],last_log_line=group[-1][0]['log_line'])
        assert row['minimum_iterations_met'] and not row['maximum_iterations_reached']
        rows.append(row)
    transitions=[]
    for previous,current in zip(rows,rows[1:]):
        if previous['iterations']==current['iterations']:continue
        step=current['step']
        current_iterations=[r for r in iteration_rows if r['step']==step]
        comparison=next(r for r in iteration_rows if r['step']==step-1 and r['iteration']==current['iterations'])
        final=current_iterations[-1]
        transitions.append(dict(step=step,previous_step=step-1,previous_iterations=previous['iterations'],
            current_iterations=current['iterations'],previous_final_global=previous['end_global'],
            current_final_global=current['end_global'],current_final_step_relative=current['end_over_start'],
            previous_same_iteration_global=comparison['global_ratio'],
            previous_same_iteration_step_relative=comparison['step_ratio'],
            previous_same_iteration_stop_rule=False,current_final_stop_rule=final['stop_rule_satisfied'],
            global_tolerance_met=final['global_tolerance_met'],step_relative_tolerance_met=final['step_relative_tolerance_met'],
            no_preconditioner_rebuild_in_current_step=all(r['preconditioner_reused'] for r in current_iterations),
            linear_reasons=[r['linear_reason'] for r in current_iterations],
            evidence_log_lines=[r['log_line'] for r in current_iterations],
            preconditioner_log_lines=[r['preconditioner_log_line'] for r in current_iterations]))
        assert not comparison['stop_rule_satisfied']
    assert [t['step'] for t in transitions]==[4,23]
    assert all(t['no_preconditioner_rebuild_in_current_step'] for t in transitions)
    source_names=['vendor/svMultiPhysics_stage_q/Code/Source/solver/'+n for n in
        ['petsc_impl.cpp','Integrator.cpp','output.cpp','initialize.cpp','utils.cpp']]
    late=[r for r in rows if r['step']>=55]
    audit=dict(all_pass=True,reference_norm=reference,reference_definition='R_ref = R^(1,1), fixed for the entire fresh run; k is 1-based',
        residual_measure='Initial KSP history norm for each nonlinear correction, in the solver scaled system; this is the NS Ri measure',
        record_timing='Last = last logged initial norm of a linear correction, not an extra residual evaluation after the final correction',
        plotted_values='Ratios reconstructed from high-precision initial KSP norms; match all printed NS Ri/R0 and Ri/R1 within print rounding',
        max_relative_difference_from_rounded_NS_rows=max_print_error,
        nonlinear_tolerance=tol,minimum_iterations=minimum,maximum_iterations=maximum,
        stopping_rule='k >= max_iterations OR (k >= min_iterations AND (R/R_ref <= tol OR R/R_step_start <= tol))',
        all_71_steps_match_first_eligible_stopping_record=True,maximum_iteration_cap_never_reached=True,
        nonlinear_records=len(iteration_rows),time_steps=len(rows),
        last_record_below_global_tolerance_steps=sum(r['global_tolerance_met'] for r in rows),
        last_record_below_step_relative_tolerance_steps=sum(r['step_relative_tolerance_met'] for r in rows),
        start_already_below_global_tolerance_steps=[r['step'] for r in rows if r['start_global']<=tol],
        observed_late_end_global_range=[min(r['end_global'] for r in late),max(r['end_global'] for r in late)],
        numerical_floor_established=False,residual_is_not_solution_error=True,
        near_zero_substitution_threshold_approx=10*np.finfo(float).eps**2,
        minimum_initial_KSP_norm=min(r['norm'] for r in iteration_rows),
        near_zero_substitution_triggered=False,transitions=transitions,
        log_sha256=sha(CASE/'run/solver.log'),solver_xml_sha256=sha(CASE/'run/solver.xml'),
        source_code_sha256={n:sha(ROOT/n) for n in source_names})
    assert audit['minimum_initial_KSP_norm']>audit['near_zero_substitution_threshold_approx']
    write_csv(OUT/'data/nonlinear_step_audit.csv',rows)
    write_csv(OUT/'data/nonlinear_iteration_audit.csv',iteration_rows)
    write_csv(OUT/'data/nonlinear_transition_audit.csv',transitions)
    excerpt_lines=set()
    for row in iteration_rows:
        if row['step'] in [3,4,22,23]:
            excerpt_lines.update([row['log_line'],row['preconditioner_log_line']])
    log_lines=log.splitlines()
    (OUT/'data/nonlinear_transition_log_excerpt.txt').write_text(
        'Source: ../run/solver.log; original 1-based line numbers are retained.\n'+
        '\n'.join(f'{n}: {log_lines[n-1]}' for n in sorted(excerpt_lines))+'\n')
    dump(OUT/'NONLINEAR_RESIDUAL_AUDIT.json',audit)
    return rows,audit


def style_axis(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',which='major',color='#dfe3e6',lw=.6)
    ax.tick_params(labelsize=9)


def nonlinear_panels(axes,rows,audit):
    x=np.array([r['step'] for r in rows])
    start=np.array([r['start_global'] for r in rows])
    reduction=np.array([r['reduction_orders'] for r in rows])
    ax=axes[0]
    ax.semilogy(x,start,color=ORANGE,lw=1.6,marker='o',ms=2.5,markevery=3,label='First logged residual in each step')
    ax.axhline(audit['nonlinear_tolerance'],color=GRAY,ls='--',lw=1.2,label=r'Global nonlinear tolerance: $10^{-10}$')
    ax.set(title='(A) Global transient evolution',xlabel='Simulation time step',
        ylabel=r'Step-start residual, $R_{\mathrm{start}}^{(n)}/R_{\mathrm{ref}}$',xlim=(.2,72),ylim=(1e-15,20))
    ax.legend(frameon=False,fontsize=8.5,loc='upper right')
    ax.annotate('Observed late-step plateau\nNumerical floor not established',xy=(64,start[63]),xytext=(31,5e-6),
        fontsize=8.5,color=GRAY,arrowprops=dict(arrowstyle='->',color=GRAY,lw=.9),va='center')
    ax=axes[1]
    ax.bar(x,reduction,width=.78,color=[ORANGE if n in [4,23] else BLUE for n in x],zorder=3)
    ax.set(title='(B) Nonlinear reduction within each step',xlabel='Simulation time step',
        ylabel=r'Orders reduced: $\log_{10}(R_{\mathrm{start}}^{(n)}/R_{\mathrm{last}}^{(n)})$',xlim=(.2,72),ylim=(0,15.4))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True,nbins=6))
    for step,tx,ty,label in [(4,10,14.2,'Step 4: 4 → 3 iterations'),(23,29,11.4,'Step 23: 3 → 2 iterations')]:
        ax.annotate(label,xy=(step,reduction[step-1]+.2),xytext=(tx,ty),fontsize=8.5,color=ORANGE,
            arrowprops=dict(arrowstyle='->',color=ORANGE,lw=.9),va='center')
    ax.text(.50,.51,'10 orders = a ten-billion-fold\nreduction within one step',transform=ax.transAxes,fontsize=8.5,color=GRAY)
    ax.annotate('Late steps start below\nthe global tolerance',xy=(63,reduction[62]+.2),xytext=(43,4.2),
        fontsize=8.5,color=GRAY,arrowprops=dict(arrowstyle='->',color=GRAY,lw=.9),va='center')
    for ax in axes:style_axis(ax)


def linear_panels(axes,history,summary,rtol):
    ax=axes[0]
    for reason,color,marker,label in [('CONVERGED_RTOL','#176b91','o','Relative-tolerance stop'),
                                     ('CONVERGED_ATOL','#c55c26','^','Absolute-tolerance stop')]:
        rs=[r for r in summary if r['reason']==reason]
        ax.semilogy([r['linear_solve_index'] for r in rs],[r['final_true_relative_residual'] for r in rs],
            marker,ms=3,color=color,label=label)
    ax.axhline(rtol,color=GRAY,lw=1,ls='--',label=r'Linear relative tolerance: $10^{-10}$')
    ax.set(xlabel='Linear-solve index',ylabel=r'Final true residual, $\|b-Ax\|/\|b\|$',
        title='(C) Linear solve termination',xlim=(0,len(summary)+1),ylim=(7e-11,5e-8))
    ax.legend(frameon=False,fontsize=8.5,loc='upper left')
    ax=axes[1]
    for k,color,label in [(0,BLUE,'Step 1, iteration 1'),(83,'#379778','Step 30, iteration 1'),(166,'#c55c26','Step 71, iteration 2')]:
        monitors=history['linear_solves'][k]['petsc_monitor']
        ax.semilogy([m['iteration'] for m in monitors],[m['true_relative_residual'] for m in monitors],
            color=color,lw=1.3,label=label)
    ax.axhline(rtol,color=GRAY,lw=1,ls='--')
    ax.set(xlabel='GMRES iteration',ylabel=r'True residual, $\|b-Ax\|/\|b\|$',title='(D) Selected GMRES histories')
    ax.legend(frameon=False,fontsize=8.5,loc='upper right')
    for ax in axes:style_axis(ax)


def save(fig,name):
    for suffix in ['png','pdf','svg']:
        fig.savefig(OUT/'figures'/f'{name}.{suffix}',dpi=320,facecolor='white')


def render_residuals(history,summary,steps,rtol,atol):
    rows,audit=audit_nonlinear(history)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,'axes.titlesize':12,
        'axes.linewidth':.8,'xtick.direction':'out','ytick.direction':'out','figure.facecolor':'white',
        'axes.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    # Primary figure: distinct quantities on distinct panels, with no connected
    # last-record envelope that could be mistaken for one iterative history.
    fig,axes=plt.subplots(1,2,figsize=(13.0,6.1))
    nonlinear_panels(axes,rows,audit)
    fig.subplots_adjust(left=.072,right=.98,bottom=.28,top=.78,wspace=.30)
    fig.suptitle('Nonlinear Residual Evolution and Within-Step Reduction',y=.975,fontsize=15)
    fig.text(.5,.887,r'$R_{\mathrm{ref}}=R^{(1,1)}$ is fixed for the entire run; '
        r'$R_{\mathrm{norm}}^{(n,k)}=R^{(n,k)}/R_{\mathrm{ref}}$ ($k$ is 1-based).',ha='center',fontsize=11)
    fig.text(.5,.835,'Start / last = first / last logged records in each step. R is the nonlinear solver\'s scaled residual measure.',ha='center',fontsize=8.7,color=GRAY)
    fig.text(.072,.181,'Stopping rule: at least 2 iterations, then global OR step-relative normalized residual ≤ 10⁻¹⁰.',fontsize=9)
    fig.text(.072,.128,'Steps 4 and 23 use one fewer iteration than the preceding step; tolerance is met, with no preconditioner rebuild.',fontsize=9)
    fig.text(.072,.075,'A lower nonlinear residual indicates better satisfaction of the discretized nonlinear equations within that time step.',fontsize=8.7)
    fig.text(.072,.029,'Residual magnitude is not a CFD solution-error estimate. The late plateau is observed, not a verified floating-point floor.',fontsize=8.7,color=GRAY)
    save(fig,'residual_nonlinear');plt.close(fig)
    # Supplement: each vertical connector represents only ONE time step.
    fig,ax=plt.subplots(figsize=(12.2,5.9))
    x=np.array([r['step'] for r in rows]);starts=np.array([r['start_global'] for r in rows]);ends=np.array([r['end_global'] for r in rows])
    ax.set_yscale('log');ax.vlines(x,ends,starts,color='#b6bdc4',lw=.9,zorder=1)
    ax.scatter(x,starts,s=16,facecolors='white',edgecolors=ORANGE,lw=1,zorder=3,label='First logged record')
    ax.scatter(x,ends,s=13,color=BLUE,zorder=3,label='Last logged record')
    ax.axhline(audit['nonlinear_tolerance'],color=GRAY,lw=1.1,ls='--',label=r'Global nonlinear tolerance: $10^{-10}$')
    ax.set(xlim=(.2,72),ylim=(1e-15,20),xlabel='Simulation time step',ylabel=r'Logged residual / fixed $R_{\mathrm{ref}}$')
    ax.legend(frameon=False,fontsize=9,loc='upper right')
    for step,tx,ty,text in [(4,12,2,'Step 4: 3 iterations\nGlobal + step-relative tolerance met'),
                           (23,31,1e-5,'Step 23: 2 iterations\nGlobal tolerance met')]:
        ax.scatter([step],[ends[step-1]],s=62,facecolors='none',edgecolors='#bc582d',lw=1.4,zorder=5)
        ax.annotate(text,xy=(step,ends[step-1]*1.5),xytext=(tx,ty),fontsize=9,color='#a55029',
            bbox=dict(facecolor='white',edgecolor='none',alpha=.94,pad=3),
            arrowprops=dict(arrowstyle='->',color='#bc582d',lw=1.))
    style_axis(ax)
    fig.subplots_adjust(left=.09,right=.98,bottom=.24,top=.79)
    fig.suptitle('Nonlinear Residual Reduction Within Each Time Step',y=.97,fontsize=15)
    fig.text(.5,.89,r'Each vertical connector pairs the first and last records of ONE step. $R_{\mathrm{ref}}=R^{(1,1)}$ is fixed.',ha='center',fontsize=10)
    fig.text(.09,.13,'No lines connect last records across steps. At steps 4 and 23, fewer iterations explain the higher last-record values.',fontsize=9)
    fig.text(.09,.075,'Neither event rebuilds the preconditioner. All steps satisfy the actual nonlinear stopping rule; the iteration cap is never reached.',fontsize=8.7)
    fig.text(.09,.025,'The low residual plateau is not an established numerical floor or a physical-accuracy estimate.',fontsize=9,color=GRAY)
    save(fig,'residual_nonlinear_pairs');plt.close(fig)
    # All residual diagnostics, retained for the existing presentation and index.
    fig,axes=plt.subplots(2,2,figsize=(13.0,9.4))
    nonlinear_panels(axes[0],rows,audit);linear_panels(axes[1],history,summary,rtol)
    fig.subplots_adjust(left=.075,right=.98,bottom=.13,top=.86,wspace=.30,hspace=.51)
    fig.suptitle('Solver Residual Diagnostics: Nonlinear and Linear Convergence',y=.975,fontsize=15)
    fig.text(.5,.918,r'Nonlinear reference: $R_{\mathrm{ref}}=R^{(1,1)}$, fixed for all 71 steps. Within-step reduction is evaluated separately.',ha='center',fontsize=10)
    fig.text(.075,.065,'Nonlinear stop: at least 2 iterations and global OR step-relative ratio ≤ 10⁻¹⁰. Linear stop: relative 10⁻¹⁰ OR absolute 10⁻²⁴.',fontsize=8.7)
    fig.text(.075,.024,'Residuals assess the discretized algebraic solve; they do not establish CFD solution accuracy or a floating-point floor.',fontsize=9,color=GRAY)
    save(fig,'residual_convergence')
    fig.canvas.draw()
    for name,ax in zip(['linear_termination','gmres_histories'],axes[1]):
        bounds=ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).expanded(1.05,1.08)
        for suffix in ['png','pdf','svg']:
            fig.savefig(OUT/'figures'/f'residual_{name}.{suffix}',dpi=320,bbox_inches=bounds,facecolor='white')
    plt.close(fig)
    dump(OUT/'RESIDUAL_FIGURE_VALIDATION.json',dict(all_pass=True,white_background=True,
        nonlinear_main_panels=['Global step-start evolution','Within-step log10 reduction'],
        paired_connector_supplement=True,global_reference_explicit=True,actual_nonlinear_tolerance_shown=True,
        numerical_floor_not_claimed=True,residual_not_interpreted_as_solution_error=True,steps_4_23_source_and_log_audited=True,
        plotted_norms_match_printed_NS_rows=True,high_precision_source='Initial KSP history norms, matching solver code',
        raw_log_sha256=audit['log_sha256'],audit_sha256=sha(OUT/'NONLINEAR_RESIDUAL_AUDIT.json'),
        plotter_sha256=sha(__file__),figures={str(p.relative_to(OUT)):sha(p) for p in (OUT/'figures').glob('residual_*')}))
    write_review(rows,audit)
    return audit


def write_review(rows,audit):
    by_step={r['step']:r for r in rows}
    examples='\n'.join(f"| {n} | {by_step[n]['iterations']} | {by_step[n]['start_global']:.4e} | {by_step[n]['end_global']:.4e} | {by_step[n]['reduction_orders']:.3f} |" for n in [3,4,10,20,22,23,50,71])
    text=f'''# 残差图科学表达修订与原始日志核验

本次只读取现有求解日志并更新残差图、说明和 PPT。原速度、压力、WSS、网格及三段固定轴旋转视频保持不变。

## 主图与辅助图

- [非线性双面板主图 PNG](figures/residual_nonlinear.png) / [PDF](figures/residual_nonlinear.pdf) / [SVG](figures/residual_nonlinear.svg)：A 仅展示每步初始残差随全程演化；B 展示每步内部残差下降的数量级 `log10(start/last)`，柱高 10 即该步下降 10 个数量级。
- [逐步首末配对图 PNG](figures/residual_nonlinear_pairs.png) / [PDF](figures/residual_nonlinear_pairs.pdf) / [SVG](figures/residual_nonlinear_pairs.svg)：每条竖线只连接同一步的第一条与最后一条记录，不连接不同步的末值。
- [四面板诊断总图](figures/residual_convergence.png)：A/B 为上述两个不同过程，C/D 为原始线性求解结束残差与代表性 GMRES 历史。

图中统一使用白底和英文。不存在两条首末残差包络同时连续连接的旧画法。

## R 的含义与固定参考值

此处 R 是本次求解器用于非线性收敛判断的**缩放线性系统初始残差范数**，对应每次非线性修正关联的 KSP 历史第一项。`petsc_impl.cpp` 将它传入 `FSILS.RI.iNorm`；`Integrator.cpp` 使用它判断非线性停止，`output.cpp` 将相对量打印成 NS 行的 `Ri/R0` 和 `Ri/R1`。它不是分别报告的动量残差或质量残差，也不是以 Pa 表示的压力误差。

当前是 fresh run（`Continue_previous_simulation=false`）。全程参考量只在第一次确定：

`R_ref = R^(1,1) = {audit['reference_norm']:.12e}`，其中非线性记录编号 k 从 1 开始。

所有时间步均使用同一个 `R_ref`。因此 Panel A 的每步第一条记录无需等于 1；如果除以该步自己的第一条记录，才应等于 1。Panel B 使用同一步的首末比值，公共参考量会抵消。

新图采用 KSP 监测输出的高精度初始范数重构比值。全部 167 条结果均与原 NS 文本的 `Ri/R0`、`Ri/R1` 在四位有效数字的打印误差内一致；最大相对差约 {audit['max_relative_difference_from_rounded_NS_rows']:.3e}。没有改变原日志，也没有以平滑或插值替代迭代数据。

**最后记录的时点：**记录在本次线性修正执行后打印，但其中用于非线性判断的 R 来自本次修正开始时。它是最后一条已有非线性记录，不是最后一次更新后额外重新装配并评估的残差。图中使用 `last`，避免把它表述为新增评估得到的精确终态误差。

## 实际停止条件与容差线

配置为最少 2 次、最多 12 次非线性迭代，非线性容差为 **1e-10**。源码停止规则为：达到最大迭代次数，或者在满足最少次数后，以下任一条件成立：

1. `R/R_ref ≤ 1e-10`：相对全程固定参考的条件。
2. `R/R_step_start ≤ 1e-10`：相对本步初始参考的条件。

图中的虚线是第一项的真实非线性容差，**不是浮点下限**，也不是 PETSc 线性绝对容差。全 71 步都在第一次满足上述收敛条件的记录处结束；没有任何一步触及 12 次上限。所有末记录均满足全程参考条件，其中 {audit['last_record_below_step_relative_tolerance_steps']} 步同时满足步内相对条件。

后期每步开始时已经低于全程容差，但仍必须执行最少 2 次迭代。因此 Panel B 后期下降量很小不表示该步未收敛，也不能要求每一步都再下降 10 个数量级。

## 第 4、23 步末值跳高的核验

| 时间步 | 非线性记录数 | 起始 R/R_ref | 最后 R/R_ref | 步内下降数量级 |
|---:|---:|---:|---:|---:|
{examples}

**第 3 → 4 步：**第 3 步在第 3 条记录时，全程相对残差约为 1.412e-10，步内相对残差约为 1.120e-10，均未达到阈值，必须继续第 4 次。第 4 步在第 3 条记录时，这两项已分别约为 7.756e-11 和 7.096e-11，因此在 3 次后停止。末记录高于前一步的末记录，是少进行一次迭代后的正常停止结果。

**第 22 → 23 步：**第 22 步的第 2 条记录仍为 1.594e-10（全程参考）和 1.824e-6（步内参考），均未达到阈值，因此继续第 3 次。第 23 步的第 2 条记录为 8.586e-11（全程参考）和 1.835e-6（步内参考），全程条件已经满足，且达到最少 2 次，因此正常停止。不能将此误写成“第 23 步满足了 1e-10 的步内相对下降”。

这些步骤的 `SV13Q_BEGIN` 均记录 `reuse=1`、`rebuild_reason=REUSE`，没有发生预条件器重建；第 4、23 步的线性求解全部为 `CONVERGED_RTOL`，没有失败或恢复重试。残差记录定义在这些点也没有改变。

证据：[带原始行号的日志摘录](data/nonlinear_transition_log_excerpt.txt)、[逐非线性记录 CSV](data/nonlinear_iteration_audit.csv)、[逐步 CSV](data/nonlinear_step_audit.csv)、[跳变核验 CSV](data/nonlinear_transition_audit.csv)、[完整机器核验 JSON](NONLINEAR_RESIDUAL_AUDIT.json)。

## 低残差平台的正确解释

约 1e-14 的平台只标注为**观察到的平台**。当前没有通过高精度计算、缩放对比或其他独立实验确定其具体来源，因此不将它命名为已证实的 round-off floor、归一化下限或双精度极限。源码的近零替换逻辑也没有被本日志的初始 KSP 范数触发。

图注使用：**A lower nonlinear residual indicates better satisfaction of the discretized nonlinear equations within that time step.**

残差不是 CFD 解相对真实解的误差估计。很小的代数残差不能证明网格足够细、时间步足够小、边界条件符合真实生理，也不能声称解达到 1e-14 的物理精度。本次没有开展新的网格或时间步收敛研究。

## 源码依据与复现

- [Integrator.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/Integrator.cpp)：约 947–989 行，固定参考、本步参考与停止条件；约 151–162 行，求解、修正与日志输出的先后顺序。
- [petsc_impl.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/petsc_impl.cpp)：约 318–321 行，KSP 历史第一项赋给 `initNorm`。
- [output.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/output.cpp)：约 124–126 行，`Ri/R0`、`Ri/R1` 的打印定义。
- [utils.cpp](../../../vendor/svMultiPhysics_stage_q/Code/Source/solver/utils.cpp)：约 141–159 行，近零判定。

源码与日志 SHA-256 均保存在审计 JSON。新增绘图入口为 `scripts/flow_2mmps/residual_figures.py`，读取原日志后即可重绘，不需要 CFD、WSS 或轨迹重算。永久测试位于 `tests/flow_2mmps/test_residual_semantics.py`。
'''
    (OUT/'RESIDUAL_REVIEW_ZH.md').write_text(text)


def main():
    import csv
    dt=float(ET.parse(CASE/'run/solver.xml').find('.//Time_step_size').text)
    history=parse_solver_log((CASE/'run/solver.log').read_text(),dt)
    summary=[]
    with (OUT/'data/residual_linear_solves.csv').open() as f:
        for row in csv.DictReader(f):
            row['linear_solve_index']=int(row['linear_solve_index'])
            row['final_true_relative_residual']=float(row['final_true_relative_residual'])
            summary.append(row)
    options=(CASE/'run/PETSC_OPTIONS.txt').read_text().split()
    rtol,atol=[float(options[options.index(k)+1]) for k in ['-ksp_rtol','-ksp_atol']]
    audit=render_residuals(history,summary,None,rtol,atol)
    compute=json.loads((OUT/'COMPUTE_VALIDATION.json').read_text())
    generated=list((OUT/'data').iterdir())+list((OUT/'figures').glob('residual_*'))
    compute['outputs_sha256']={str(p.relative_to(OUT)):sha(p) for p in generated if p.is_file()}
    compute['residual_scientific_revision']=dict(audit_file='NONLINEAR_RESIDUAL_AUDIT.json',
        audit_sha256=sha(OUT/'NONLINEAR_RESIDUAL_AUDIT.json'),plotter_sha256=sha(__file__),no_field_recomputation=True)
    dump(OUT/'COMPUTE_VALIDATION.json',compute)
    media=json.loads((OUT/'MEDIA_VALIDATION.json').read_text())
    for video in media['videos']:assert sha(OUT/video['file'])==video['sha256']
    media['figures']=[dict(file=str(p.relative_to(OUT)),sha256=sha(p)) for p in sorted((OUT/'figures').glob('*'))]
    media['residual_figures_refreshed_without_video_rerender']=True
    media['residual_figure_validation_sha256']=sha(OUT/'RESIDUAL_FIGURE_VALIDATION.json')
    dump(OUT/'MEDIA_VALIDATION.json',media)
    print(json.dumps(dict(all_pass=True,reference=audit['reference_norm'],transitions=audit['transitions'],
        floor_established=False,step_relative_tolerance_passed=audit['last_record_below_step_relative_tolerance_steps']),indent=2))


if __name__=='__main__':main()
