#!/usr/bin/env python3
"""New CPU-rendered diagnostics from complete logs and measured saved fields."""
import json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
os.environ['MPLCONFIGDIR']=str(ROOT/'outputs/sv1_1/plot_cache');os.environ['LIBGL_ALWAYS_SOFTWARE']='1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from sv_validation.sv11 import load,REPORT
from sv_validation.provenance import write_json,sha256
font=Path('/mnt/c/Windows/Fonts/msyh.ttc')
if font.exists():fm.fontManager.addfont(str(font));plt.rcParams['font.family']=fm.FontProperties(fname=str(font)).get_name()
plt.rcParams.update({'axes.unicode_minus':False,'font.size':11,'figure.dpi':130,'savefig.dpi':170})
base=load('baseline_solver_diagnosis');old=base['linear_solves']
run=load('petsc_short_execution') if (REPORT/'petsc_short_execution.json').exists() else None
new=run['history']['linear_solves'] if run else []
qc=load('petsc_short_qc') if (REPORT/'petsc_short_qc.json').exists() else {'boundary_history':[],'states':[],'steady_intervals':[]}
files=[]
def save(fig,name,scope='diagnostic comparison'):
    fig.tight_layout();path=REPORT/name;fig.savefig(path);plt.close(fig)
    files.append({'path':str(path.relative_to(ROOT)),'sha256':sha256(path),'scope':scope})
def step_summary(rows,key,fn):
    steps=sorted({r['step'] for r in rows});return steps,[fn([r[key] for r in rows if r['step']==s]) for s in steps]
colors=['#b45309','#006d9c']
fig,axes=plt.subplots(3,1,figsize=(11,10),sharex=True)
fig.suptitle('换用 PETSc 后，每一步的线性方程还能正常解出来吗？',fontsize=16)
for rows,label,color in zip([old,new],['SV1 FSILS 历史','SV1.1 PETSc'],colors):
    if not rows:continue
    steps,failed=step_summary(rows,'linear_converged',lambda v:sum(not x for x in v))
    axes[0].plot(steps,failed,'o-',label=label,color=color)
    steps,it=step_summary(rows,'linear_iterations',np.mean);axes[1].plot(steps,it,'o-',label=label,color=color)
    # Ill-conditioned FSILS residuals reset to zero are not precision evidence.
    good=[r for r in rows if r['residual_trustworthy'] and r['linear_residual_relative'] is not None and r['linear_residual_relative']>0]
    axes[2].scatter([r['step']+(r['nonlinear_iteration']-2)*.1 for r in good],[r['linear_residual_relative'] for r in good],label=label,color=color,alpha=.7)
axes[0].set_ylabel('未收敛次数');axes[1].set_ylabel('每次求解平均迭代数');axes[2].set_ylabel('记录的终止残差比');axes[2].set_yscale('log');axes[2].set_xlabel('时间步')
axes[2].axhline(1e-10,color='gray',ls='--',label='相对容差 1e-10（另有绝对容差）')
for ax in axes:ax.grid(alpha=.25);ax.legend(fontsize=9)
fig.text(.5,.003,'FSILS 外层 NS 次数与 PETSc GMRES 次数不可等同；病态警告触发的 FSILS 零残差不作为精度证据。',ha='center',fontsize=9)
save(fig,'linear_solver_before_after.png')

fig,ax=plt.subplots(figsize=(11,5));ax.set_title('每一步需要多少次线性迭代？')
for rows,label,color in zip([old,new],['FSILS NS 外层（历史）','PETSc GMRES'],colors):
    ax.plot([r['step']+(r['nonlinear_iteration']-1)/5 for r in rows],[r['linear_iterations'] for r in rows],'o-',ms=3,label=label,color=color)
ax.set(xlabel='时间步（同一步按非线性迭代错开）',ylabel='线性迭代次数');ax.legend();ax.grid(alpha=.25)
save(fig,'linear_iterations_over_time.png')

fig,ax=plt.subplots(figsize=(11,5));ax.set_title('质量守恒误差有没有持续下降？')
old_boundary=base.get('boundary_history',[])
if not old_boundary:
    old_boundary=list({r['step']:r['boundary_at_end_of_step'] for r in old if r.get('boundary_at_end_of_step')}.values())
for rows,label,color in zip([old_boundary,qc['boundary_history']],['FSILS 原生积分（历史）','PETSc 原生积分'],colors):
    if rows:ax.semilogy([r['step'] for r in rows],[max(r['epsilon_mass'],1e-17) for r in rows],
                       marker='s' if '历史' in label else 'o',linestyle='-' if '历史' in label else '--',
                       markerfacecolor='none',markersize=8 if '历史' in label else 4,label=label,color=color)
if qc['states']:ax.scatter([s['step'] for s in qc['states']],[s['epsilon_mass'] for s in qc['states']],s=85,marker='x',color='#d62828',label='PETSc VTU 独立积分',zorder=5)
ax.axhline(1e-6,color='#d62828',ls='--',label='最终验收线 1e-6');ax.set(xlabel='时间步',ylabel='|Qout − Qin| / Qtarget');ax.legend();ax.grid(alpha=.25)
ax.text(.98,.15,'重合的轨迹表示两组误差接近，并非缺少曲线',ha='right',transform=ax.transAxes,fontsize=9)
save(fig,'mass_balance_over_time.png')

fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True);fig.suptitle('外层流体迭代是否正常下降？')
for ax,rows,title in zip(axes,[old,new],['FSILS 历史','PETSc']):
    for step in sorted({r['step'] for r in rows}):
        selected=[r for r in rows if r['step']==step]
        ax.semilogy([r['nonlinear_iteration'] for r in selected],[r['nonlinear_Ri_over_R0'] for r in selected],'.-',label=f'步 {step}')
    ax.set(title=title,xlabel='非线性迭代',ylabel='Ri / R0');ax.axhline(1e-10,color='gray',ls='--');ax.grid(alpha=.25)
    if rows:ax.legend(ncol=2,fontsize=8)
save(fig,'nonlinear_convergence.png')

fig,ax=plt.subplots(figsize=(10,4));ax.set_title('最终是否真正达到稳态？')
intervals=qc['steady_intervals']
if intervals:
    ax.semilogy([r['step'] for r in intervals],[r['E_u'] for r in intervals],'o-',label='E_u')
    ax.semilogy([r['step'] for r in intervals],[r['E_Q'] for r in intervals],'o-',label='E_Q');ax.legend()
else:
    ax.axis('off');ax.text(.5,.55,f'未达到稳态验收\n只有 {len(qc["states"])} 个 PETSc 已保存状态\n需要至少 6 个状态形成连续 5 个区间',ha='center',va='center',fontsize=16,transform=ax.transAxes)
save(fig,'steady_convergence.png','availability / steady gate evidence')

if qc['states']:
    q=qc['states'][-1];fig,ax=plt.subplots(figsize=(10,5));ax.set_title(f'第 {q["step"]} 步各端口流量（瞬态诊断，未经稳态验收）')
    values=[q['Q_in_m3_s'],*[q['outlet_flows_m3_s'][n] for n in ['OUTLET_01','OUTLET_02','OUTLET_03']],q['Q_out_total_m3_s']]
    ax.bar(['Qin','OUTLET_01','OUTLET_02','OUTLET_03','Qout total'],np.array(values)/q['Q_target_m3_s'],color=['#006d9c','#609a63','#609a63','#609a63','#609a63'])
    ax.axhline(1,color='gray',ls='--');ax.set_ylabel('实际积分流量 / Qtarget');save(fig,'flux_balance.png','unaccepted transient diagnostic')

fig,axes=plt.subplots(1,2,figsize=(11,5));fig.suptitle('线性求解器的时间和内存消耗（描述性比较）')
prior=load('flow_execution','sv1')['runs'][0]
labels=['FSILS 历史','PETSc 短程'];times=[prior['elapsed_s'],run['elapsed_s'] if run else np.nan];rss=[prior['peak_rss_kib']/1024,run['peak_rss_kib']/1024 if run else np.nan]
axes[0].bar(labels,times,color=colors);axes[0].set_ylabel('Python 单调时钟 / 秒');axes[1].bar(labels,rss,color=colors);axes[1].set_ylabel('最大单进程 RSS / MiB')
fig.text(.5,.006,'均为 4 MPI ranks、OMP=1；RSS 不是 MPI 总和。历史 GNU wall 与单调时钟不一致，不据此推断加速比。',ha='center',fontsize=9)
save(fig,'solver_resource_usage.png')
write_json(REPORT/'visuals.json',{'files':files,'accepted_solution_available':(REPORT/'accepted_solution.json').exists(),
           'withheld_until_accepted':['velocity_global.png','velocity_slices.png','pressure_global.png','pressure_sections.png','outlet_flow_split.png'],
           'renderer':'Matplotlib Agg CPU; no solver/GPU rendering','chinese_font':str(font)})
print('Generated',len(files),'new evidence figures')
