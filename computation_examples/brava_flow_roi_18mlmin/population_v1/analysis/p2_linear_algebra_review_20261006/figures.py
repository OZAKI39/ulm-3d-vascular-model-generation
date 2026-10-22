from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
g=read(HERE/'P2_PETSC_GRAPH_MEASUREMENTS.json');a=read(HERE/'ASM_LU_MEMORY_AUDIT.json');p=read(HERE/'P1_LINEAR_SOLVER_PARITY.json');f=read(HERE/'P2_LINEAR_SYSTEM_RESOURCE_RESULT.json');res=read(HERE/'RESOURCES.json');limit=int(res['memory.max'])/2**30;reserve=.2*limit
fig,axes=plt.subplots(1,2,figsize=(13,5.6),layout='constrained')
labels=['Native init','Mat + RHS','ASM overlap 0','ASM overlap 1','ASM overlap 2'];values=[g['RSS_before_MatSetUp_GiB'],g['RSS_after_MatSetUp_GiB']]+[x['measured_total_rank_RSS_after_overlap_GiB'] for x in a['overlap']]
axes[0].bar(labels,values,color='#2477a8');axes[0].tick_params(axis='x',rotation=25);axes[0].set(ylabel='Sum of 8 rank RSS [GiB]',title='Measured full P2 graph stages (no LU)');axes[0].set_ylim(0,max(values)*1.27)
for i,v in enumerate(values):axes[0].text(i,v+.12,f'{v:.2f}',ha='center')
payload=[g['CSR_value_column_row_payload_bytes']/2**30,(g['owned_rhs_bytes']+g['ghost_rhs_bytes'])/2**30,a['overlap'][2]['submatrix_value_index_payload_GiB']]
axes[1].bar(['Matrix CSR','RHS + ghosts','ASM2 submatrices'],payload,color='white',edgecolor='#a66421',hatch='///');axes[1].set(ylabel='Payload only [GiB]',title='Calculated bytes from measured counts');axes[1].set_ylim(0,max(payload)*1.55)
for i,v in enumerate(payload):axes[1].text(i,v+.05,f'{v:.3f}',ha='center')
axes[1].text(.03,.93,'LU symbolic structures: NOT_RUN\nNumeric LU factors: NOT_RUN\nMatGetInfo memory=0 means unavailable',transform=axes[1].transAxes,va='top',fontsize=10)
fig.suptitle(f'P2 memory accounting | cgroup {limit:.2f} GiB | required reserve {reserve:.2f} GiB\nRSS snapshots and calculated payloads are different quantities; do not add the two panels.',fontsize=12)
fig.savefig(HERE/'01_P2_MEMORY_BREAKDOWN.png',dpi=300);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(13,5.6),layout='constrained');x=np.arange(3)
axes[0].bar(x-.18,[v['row_ratio'] for v in a['overlap']],.36,label='Subdomain rows / global rows',color='#2477a8');axes[0].bar(x+.18,[v['submatrix_nnz_global_ratio'] for v in a['overlap']],.36,label='Submatrix NNZ / global NNZ',color='#57a58d');axes[0].set(xticks=x,xticklabels=['Overlap 0','Overlap 1','Overlap 2'],ylabel='Measured ratio',title='Native ownership + PETSc overlap');axes[0].axhline(1,color='black',lw=.8);axes[0].legend(fontsize=9)
est=[v['estimated_LU_value_index_GiB'] for v in a['overlap']];axes[1].bar(x,est,facecolor='white',edgecolor='#a66421',hatch='///');axes[1].set(xticks=x,xticklabels=['Overlap 0','Overlap 1','Overlap 2'],ylabel='LU value/index estimate [GiB]',title='Estimated factors only: NOT measured P2 fill');axes[1].axhline(limit,color='#9b2525',label=f'Cgroup total {limit:.2f} GiB');axes[1].axhline(limit-reserve,color='#9b2525',ls='--',label=f'Total less reserve {limit-reserve:.2f} GiB');axes[1].legend(fontsize=9,loc='upper left');axes[1].set_ylim(0,max(est)*1.3)
for i,v in enumerate(est):axes[1].text(i,v+.5,f'{v:.2f}',ha='center')
fig.suptitle('Actual ASM graph growth + inherited P1 LU fill sensitivity\nLU symbolic allocates numeric arrays in this PETSc build; unsafe factorization was not attempted.',fontsize=12)
fig.savefig(HERE/'02_ASM_OVERLAP_AND_FILL.png',dpi=300);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(16,5.8),layout='constrained');names=list(p['results']);rr=[p['results'][n]['peak_tree_RSS_GiB'] for n in names];axes[0].bar(names,rr,color=['#63686f','#2477a8','#57a58d'][:len(names)]);axes[0].set(title='Measured P1 parity jobs',ylabel='Peak process-tree RSS [GiB]');axes[0].set_ylim(0,max(rr)*1.35)
for i,n in enumerate(names):axes[0].text(i,rr[i]+.12,f'{rr[i]:.2f}\n'+('PASS' if p['results'][n]['PASS'] else 'FAIL'),ha='center',fontsize=10)
if f.get('actual_P2_job_peak_tree_RSS_GiB') is not None:
 vals=[f['actual_P2_job_peak_tree_RSS_GiB'],f['peak_nonreclaimable_upper_GiB'],f['minimum_conservative_headroom_GiB']];axes[1].bar(['P2 job RSS','Cgroup nonrec.','Headroom'],vals,color=['#2477a8','#63686f','#57a58d']);axes[1].axhline(limit,color='#9b2525',label='Cgroup total');axes[1].axhline(reserve,color='#9b2525',ls='--',label='Required reserve');axes[1].legend(fontsize=9);axes[1].set_ylim(0,limit*1.12)
 if f.get('CUDA_allocation_failure'):axes[1].text(.03,.80,'Host samples only. GPU allocation failed.\nPost-setup reserve was NOT certified.',transform=axes[1].transAxes,va='top',fontsize=9,color='#9b2525')
 for i,v in enumerate(vals):axes[1].text(i,v+.3,f'{v:.2f}',ha='center')
 ranks=f.get('rank_peak_RSS_bytes',{});axes[2].bar([str(i) for i in range(8)],[ranks.get(str(i),0)/2**30 for i in range(8)],color='#2477a8');axes[2].set(xlabel='MPI rank',ylabel='Process high-water RSS [GiB]',title='P2 observed rank high-water\nBefore failed PC setup; lower bounds')
else:
 axes[1].text(.5,.5,'P2 PC setup / solve\nNOT_RUN\nNo measured peak available',ha='center',va='center',transform=axes[1].transAxes);axes[1].set(ylim=(0,limit*1.12));axes[1].axhline(limit,color='#9b2525',label='Cgroup total');axes[1].axhline(reserve,color='#9b2525',ls='--',label='Required reserve');axes[1].legend(fontsize=9);axes[2].axis('off');axes[2].text(.05,.5,f.get('status','NOT_RUN'),wrap=True,fontsize=10,transform=axes[2].transAxes)
axes[1].set(title='P2 resource gate: FAILED setup' if f.get('CUDA_allocation_failure') else 'P2 single-system resource gate',ylabel='GiB (distinct quantities, not additive)');axes[1].tick_params(axis='x',rotation=20)
fig.suptitle(f'Solver resource comparison | cgroup {limit:.2f} GiB | reserve {reserve:.2f} GiB\nP1 parity does not establish P2 field accuracy. P2 setup failed; polling missed the transient GPU allocation peak.',fontsize=12)
fig.savefig(HERE/'03_P2_SOLVER_RESOURCE_COMPARISON.png',dpi=300);plt.close(fig)
put('FIGURE_DATA.json',dict(graph_stage_RSS_GiB=dict(zip(labels,values)),CSR_payload_GiB=payload,overlap=a['overlap'],P1_results=p['results'],P2_system=f,cgroup_limit_GiB=limit,required_reserve_GiB=reserve,dpi=300))
(HERE/'OPEN_RESULTS.html').write_text('''<!doctype html><html lang="zh"><meta charset="utf-8"><title>P2/P1 线性代数与资源审计</title><style>body{max-width:1500px;margin:32px auto;font:18px/1.6 sans-serif;padding:0 20px;color:#182d39}img{width:100%;height:auto}a{color:#126a9b}figure{margin:35px 0}</style><h1>P2/P1 线性代数与资源审计</h1><p><a href="P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md">中文主报告</a> · <a href="P2_LINEAR_ALGEBRA_DECISION_GATE.json">机器判定</a> · <a href="FIGURE_DATA.json">图表原始数据</a></p><p>蓝/绿实色表示测量；斜线表示由计数计算或继承填充率估计。未运行的阶段明确标为 NOT_RUN。内存变小不等于流场精度提高。</p>'''+''.join(f'<figure><img src="{n}.png"><figcaption>{caption}</figcaption></figure>' for n,caption in [('01_P2_MEMORY_BREAKDOWN','图1：完整曲面 P2 原生图结构的实测 RSS 与理论有效载荷，两个面板不能相加。'),('02_ASM_OVERLAP_AND_FILL','图2：实测 ASM 行/结构项重复与原 P1 LU 填充率敏感性。未进行不安全的 LU symbolic/numeric 分配。'),('03_P2_SOLVER_RESOURCE_COMPARISON','图3：P1 回归与 P2 单系统资源门；不代表实际 P2 稳态流场认证。')])+'</html>')
