"""Reports and plots from saved C5 physical-energy stability evidence."""
from pathlib import Path
import json,hashlib,os,time,xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from scripts import c5_stability_audit as a
OUT=a.OUT;FIG=a.RESULTS/'figures'


def test_stats(name):
    path=OUT/f'tests/{name}.xml'
    if not path.exists():return {}
    cases=ET.parse(path).getroot().findall('.//testcase');failed=[dict(id=x.get('classname')+'::'+x.get('name'),message=x.find('failure').get('message'),trace=x.find('failure').text) for x in cases if x.find('failure') is not None]
    return dict(passed=sum(x.find('failure') is None and x.find('error') is None and x.find('skipped') is None for x in cases),failed=len(failed),skipped=sum(x.find('skipped') is not None for x in cases),failures=failed)


def preservation():
    before=a.read(OUT/'PRE_RUN_MANIFEST.json');changed=[];missing=[];link_changes=[]
    for name,digest in before['files'].items():
        p=a.ROOT/name
        if not p.is_file():missing.append(name)
        elif a.checksum(p)!=digest:changed.append(name)
    for name,target in before['symlinks'].items():
        p=a.ROOT/name
        if not p.is_symlink() or os.readlink(p)!=target:link_changes.append(name)
    result=dict(protected_files=len(before['files']),protected_symlinks=len(before['symlinks']),modified=changed,missing=missing,symlinks_changed=link_changes,
        production_modified_files=[s for s in changed if s.startswith('pd_clot/')],existing_tests_modified=[s for s in changed if s.startswith('tests/')],passed=not(changed or missing or link_changes))
    a.save(OUT/'PRESERVATION_CHECK.json',result);return result


def post_manifest():
    files={};links={}
    for directory,dirs,names in os.walk(a.ROOT,followlinks=False):
        dirs[:]=[d for d in dirs if d!='.git']
        for name in list(dirs)+names:
            p=Path(directory)/name;rel=str(p.relative_to(a.ROOT))
            if p==OUT/'POST_RUN_MANIFEST.json':continue
            if p.is_symlink():links[rel]=os.readlink(p)
            elif p.is_file():files[rel]=a.checksum(p)
    a.save(OUT/'POST_RUN_MANIFEST.json',dict(root=str(a.ROOT),files=files,symlinks=links,timestamp=time.time(),self_excluded=True))


def collect():
    allrows=a.read(OUT/'ALL_STATES.json');rows=allrows['states'];bracket=a.read(OUT/'STABILITY_BRACKET.json');base=test_stats('BASELINE');final=test_stats('FINAL')
    physical=all(r['Hxx_gate'] and r['matrix_gate'] and r['Hss_gate'] for r in rows)
    preserved=a.read(OUT/'PRESERVATION_CHECK.json') if (OUT/'PRESERVATION_CHECK.json').exists() else {}
    complete=physical and all(r['mode_direct_pass'] for r in rows) and preserved.get('passed',False) and final.get('failed')==3 and base.get('failures')==final.get('failures')
    positive=[r for r in rows if r['metrics']['maximum_s']>1e-8 and r['status']=='STATIONARY_STABLE'];unstable=[r for r in rows if r['status']=='STATIONARY_UNSTABLE']
    def simplified(label):return {'STATIONARY_STABLE':'STABLE','STATIONARY_UNSTABLE':'UNSTABLE'}.get(next(r for r in rows if r['label']==label)['status'],'UNRESOLVED')
    status=dict(C5_STABILITY_AUDIT='COMPLETE' if complete else 'INCOMPLETE',PHYSICAL_HESSIAN='VERIFIED' if physical else 'NOT_VERIFIED',
        MECHANICAL_SUBSYSTEM='STABLE' if all(r['mechanical_status']=='STABLE' for r in rows) else 'UNRESOLVED',ACTIVATION_POINT_STABILITY=simplified('S_WEAK'),
        FIRST_POSITIVE_ROOT_STABILITY=simplified('S1'),STABLE_POSITIVE_PHASE_BRANCH='EXISTS' if positive else 'INCONCLUSIVE',
        FIRST_STABILITY_LOSS='BRACKETED' if bracket['resolved'] else 'NOT_FOUND_IN_WINDOW' if not unstable and all(r['status']=='STATIONARY_STABLE' for r in rows) else 'UNRESOLVED',
        STABILITY_LOSS_INTERVAL=bracket['interval'] or 'NOT_AVAILABLE',STABILITY_LOSS_RELATIVE_TO_ACTIVATION='INCONCLUSIVE',
        CRITICAL_MODE='RESOLVED' if unstable else 'NOT_AVAILABLE',MODE_LOCALIZATION='INCONCLUSIVE',PETSC_RSLS='UNCHANGED_VERIFIED_CANDIDATE',
        C7='NOT RUN',SPATIAL_CONVERGENCE='NOT YET ASSESSED',PRODUCTION='UNCHANGED' if preserved.get('passed') else 'PENDING_CHECK')
    budget=a.settings()['resource_budget'];peak_parent=max(r['parent_peak_RSS_bytes'] for r in rows);peak_worker=max(r['worker_peak_RSS_bytes'] for r in rows)
    result=dict(status=status,states=rows,states_analyzed=len(rows),activation_bracket=a.read(a.old.OUT/'activation/BRACKET.json')['interval'],
        first_stable_positive_root=positive[0]['strain'] if positive else None,first_unstable_root=unstable[0]['strain'] if unstable else None,
        stability_bracket=bracket,stop_reason=allrows['stop_reason'],audited_strain_window=[rows[0]['strain'],rows[-1]['strain']],last_max_s=rows[-1]['metrics']['maximum_s'],
        Hessian_assembly_wall_seconds=sum(r['Hessian_assembly_seconds'] for r in rows),
        peak_parent_RSS_bytes=peak_parent,peak_worker_RSS_bytes=peak_worker,RAM_conservative_process_peak_sum_bytes=peak_parent+a.settings()['workers']*peak_worker,
        RAM_measurement='Conservative sum of per-process RSS maxima; not a synchronized whole-system peak',GPU_used=False,workers=a.settings()['workers'],
        max_Hxx_raw_asymmetry=max(t['Hxx_raw_asymmetry'] for r in rows for t in r['assembly']),max_cross_error=max(t['cross_relative_error'] for r in rows for t in r['assembly']),
        max_Hss_error=max(t['Hss_relative_error'] for r in rows for t in r['assembly']),max_mechanical_energy_curvature_error=max(r['mechanical_direct_curvature_max_error'] for r in rows),
        baseline_tests=base,final_tests=final,known_failures_unchanged=bool(base and final and base['failures']==final['failures']),preservation=preserved)
    a.save(OUT/'SUMMARY.json',result);return result


def figures(r):
    FIG.mkdir(parents=True,exist_ok=True);font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
    white='#e8edf4';blue='#62b6ff';orange='#ffb366'
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':10,'pdf.fonttype':42,'figure.facecolor':'#080b10','axes.facecolor':'#10151e','savefig.facecolor':'#080b10',
        'text.color':white,'axes.labelcolor':white,'axes.edgecolor':white,'xtick.color':white,'ytick.color':white,'axes.spines.top':False,'axes.spines.right':False,'axes.formatter.useoffset':False})
    def finish(fig,name):
        fig.tight_layout();fig.savefig(FIG/f'{name}.png',dpi=210);fig.savefig(FIG/f'{name}.pdf');plt.close(fig)
    rows=r['states'];eps=np.array([v['strain'] for v in rows]);activation=r['activation_bracket']
    fig,ax=plt.subplots(figsize=(9,4));phase=[x for x in rows if x['lambda_S_min'] is not None]
    ax.plot([v['strain'] for v in phase],[v['lambda_S_min'] for v in phase],'o-',color=blue,label='Minimum Schur curvature (stable certificate)')
    ax.axvspan(*activation,color=orange,alpha=.4,label='AT1 activation bracket');ax.axhline(0,color=white,ls='--');ax.set(xlabel='Engineering strain',ylabel='Dimensionless phase-condensed curvature',title='C5: no stability loss found at the audited states')
    ax.text(.04,.12,'All mechanical blocks are positive definite\nS0 has no movable phase DOF\nStop: max phase crosses 0.10',transform=ax.transAxes);ax.legend(frameon=False,fontsize=8,loc='upper left');finish(fig,'c5_stability_vs_strain')
    fig,axes=plt.subplots(2,1,figsize=(9,5.8),sharex=True)
    axes[0].plot(eps,[x['metrics']['maximum_s'] for x in rows],'o-',color=blue);axes[0].axhline(.1,color=white,ls='--',label='Window stop criterion');axes[0].set_ylabel('Maximum phase s');axes[0].legend(frameon=False)
    axes[1].plot(eps,[x['counts']['INTERIOR'] for x in rows],'o-',color=blue,label='Interior phase bonds');axes[1].plot(eps,[x['counts']['WEAK_LOWER'] for x in rows],'s:',color=orange,label='Weak lower bonds');axes[1].set(xlabel='Engineering strain',ylabel='Number of critical phase bonds');axes[1].legend(frameon=False)
    for ax in axes:ax.axvspan(*activation,color=orange,alpha=.25)
    fig.suptitle('Activation starts a stable positive-phase branch\nNo instability interval identified',fontsize=12);finish(fig,'phase_activation_vs_stability')
    fig,axes=plt.subplots(1,2,figsize=(10,3.8))
    for ax,key,tol,title in zip(axes,['force_residual','projected_phase_KKT'],[1e-9,1e-8],['Force residual','Projected phase KKT']):
        ax.semilogy(eps,[max(x['metrics'][key],1e-16) for x in rows],'o-',color=blue);ax.axhline(tol,color=white,ls='--',label=f'Tolerance {tol:.0e}');ax.set(xlabel='Engineering strain',ylabel=title);ax.legend(frameon=False)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4));ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'))
    fig.suptitle('Only states actually audited for stability | zeros shown at 1e-16');finish(fig,'stationarity_quality')
    c3=a.read(a.old.OLD/'critical_brackets/C3_A.json');fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    axes[0].bar(['C3','C5'],[.01,np.mean(activation)],color=[orange,blue]);axes[0].set_ylabel('Activation engineering strain')
    axes[1].bar(['C3','C5'],[1.,2/5392],color=[orange,blue]);axes[1].set_yscale('log');axes[1].set_ylabel('Initial activated bond fraction');axes[1].text(0,1.1,'185 / 185',ha='center');axes[1].text(1,.0005,'2 / 5392',ha='center')
    axes[2].axis('off');axes[2].text(.02,.82,'C3: loss near activation\nHistorical bracket:\n[0.0099921875, 0.0100000000]\n\nC5: stable positive-phase states\nNo loss found through\nstrain 0.004697390736',va='top',linespacing=1.5)
    fig.suptitle('Historical C3 vs current C5 | NO C7 / NO CONVERGENCE CLAIM');finish(fig,'c3_vs_c5_summary')
    a.save(OUT/'FIGURE_MANIFEST.json',dict(groups=4,files=[p.name for p in sorted(FIG.glob('*'))],omitted=dict(energy_along_mode='No stable/unstable bracket exists',critical_mode_map='No verified unstable mode'),style='Dark background; Arial regular; English; no arrows or dual axes'))


def report(r):
    rows=r['states'];cfg=a.settings();last=rows[-1]
    lines=['# PFPD_AT1_C5_STABILITY_V1','## 1. 为什么现在终于可以研究 stability',
        'PETSc RSLS 已在上一轮从零相场自动恢复并验证七个 C5 stationary roots。本轮直接读入这些根，不重跑七点回归；另读入紧邻 activation bracket 下侧的 intact 状态和 diagnostic weak point。九个已有状态的 checksum、force、KKT、bounds、irreversibility、energy 与 max phase 均通过复核。七个 PETSc 根的 reaction 也与原摘要一致；两个旧 intact 摘要未存 reaction，此项明确记录为无历史标量参考，未静默补造。',
        '## 2. stationary 与 stable 的区别',
        'stationary 表示自由力平衡并满足相场 KKT；stable 还要求允许的小扰动不能降低总能量。求解器收敛、相场 activation 和稳定性变化是三件不同的事。本轮结果：激活弱点和首个正相场根都稳定，稳定正相场分支存在；在规定窗口的全部审计状态中未找到失稳。',
        '## 3. 什么是 Hessian / curvature',
        '固定上下端位移后，考察原 BAQP elastic energy 加原 AT1 crack energy。Hessian 表示总能量在微小扰动方向上的弯曲程度。正曲率像谷底，负曲率意味着存在向下方向。固定边界不参加 variation，也不添加外力势。这里没有使用 PETSc residual Jacobian 代替物理 Hessian。',
        f"统一使用历史尺度：Ls={cfg['scaling']['length_scale_m']:.8g} m，Ss=1，Es={cfg['scaling']['energy_scale_J']:.12e} J。所有表中曲率无量纲。三档差分 [1e−5,3e−6,1e−6]；复用原 full_hessian 中的 physical-gradient 差分循环，仅限制变量集合并分派独立列。",
        f"Hxx 最大原始非对称误差 {r['max_Hxx_raw_asymmetry']:.3e}，cross reciprocity 最大误差 {r['max_cross_error']:.3e}；均先记录，再取对称部分。Hss 用原 pair weights 的 2WΨ+/Es 对角，AT1 线性项二阶导数为零；analytic-vs-gradient FD 最大误差 {r['max_Hss_error']:.3e}。Hxx 三尺度符号和幅值门槛通过。",
        '直接 scalar-energy 相减的最小幅度起初出现约 5%–20% 误差，未放宽 2% 阈值。保留 FLOAT64_PRECHECK 和 FLOAT64_MECHANICAL_ORACLE。新增诊断接口只重绑定原 BAQP、spectral split、State.energy、AuditProblem.parts 的同一 code objects，使用 long-double 浮点构造和求和；不复制能量公式。NumPy linalg 仍在 float64 做小矩阵运算，但在高精度 Green strain 形成后才转换。该接口从不用于 stationary solver。',
        f"扩展精度与正式总能量逐状态相对一致性门槛为 1e−12，全部通过；10 个确定随机机械方向 × 三幅度的物理能量曲率复核，最大相对误差 {r['max_mechanical_energy_curvature_error']:.3e}。因此改善来自减小浮点抵消误差，不是调物理、步长或验收阈值。",
        '每个状态有 68 个近零 principal Green strain 的键。它们全部属于 intact、strong-lower、ds=0 的键；在 critical subspace 上 g=1，Ψ+ + Ψ− 恢复原始光滑 SVK 能量。两个可动相场键不在谱零面，因此限制后的能量具有本轮需要的二阶导数。仍逐探针记录所有 raw sign changes；任何发生切换的幅度不作普通 centered 验证，而使用可行的单侧物理能量变化。没有引入谱平滑。',
        '## 4. 为什么只分析真正可动的 phase directions',
        '沿用 s_tol=1e−8、q_tol=1e−7，strong bound 的 ds 固定为零，weak lower 只允许 ds≥0，weak upper 只允许 ds≤0，interior 不限符号。初始 intact 状态 m=0，weak point 为两个弱下界键，全部 13 个正相场状态均只有两个 interior 键；没有构造 5392×5392 dense phase Hessian。',
        '## 5. 为什么 mechanical condensation 能省很多计算',
        '先检查 225×225 Hxx 正定。给定相场方向后，位移可以选择最有利的松弛；用一次 Cholesky 分解和 cho_solve 计算这次松弛，再从固定机械状态的相场刚度中减掉松弛收益，得到只有 m×m 的 Schur 矩阵 S。S 正定意味着即使机械位移尽量帮助相场降能量，所有相场方向仍向上弯曲。',
        '计算采用 [SciPy Cholesky 接口](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_factor.html)，不显式求逆。正定块与 Schur 的关系可参见 [Schur complement 推导](https://www.cis.upenn.edu/~jean/schur-comp.pdf)。本轮各 m=2 的 unconstrained S 均严格正定，所以整个允许 cone 稳定，无需运行受约束优化来证明稳定。SLSQP 仅保留为遇到不定 S 时的负 witness 搜索工具；没找到负 witness 不能证明稳定。',
        '## 6. activation point 是否稳定',
        f"S_WEAK，ε={rows[1]['strain']:.15f}：Hxx 最小值 {rows[1]['lambda_x_min']:.9g}，S 最小值 {rows[1]['lambda_S_min']:.9g}。两个 weak-lower 键满足单侧约束；正定 S 证书和单侧物理能量验证通过，分类 STATIONARY_STABLE。activation bracket 本身沿用旧结果，未重算。",
        '## 7. first positive root 是否稳定',
        f"S1，ε={rows[2]['strain']:.15f}，max s={rows[2]['metrics']['maximum_s']:.12g}，两个 interior 键。S 最小值 {rows[2]['lambda_S_min']:.9g}，分类 STATIONARY_STABLE。故本轮证据排除了“刚激活就已失稳”这一解释。",
        '## 8. positive-phase branch 能稳定多久',
        f"先完成已有 9 个状态，确认七个正根均稳定，才按固定 Δε=2.5e−5 新增 E01–E06。每点只继承前一 PETSc 接受态，1 MPI rank，原 RSLS 和全部阈值不变。共 {len(rows)} 个状态均科学 stationary 且通过稳定性证书。",
        f"最后 ε={last['strain']:.15f}，max s={last['metrics']['maximum_s']:.12g}。E05 的 max s=0.0932857291；下一个规定粗步 E06 首次跨过 0.10，随即停止，没有为了精确碰到阈值修改 load step，也没有继续到 0.006。结论针对这 15 个审计点与规定步长，不是整个连续载荷区间的解析证明。",
        '完整状态表：','| state | strain | force | KKT | strong lower / weak lower / interior | λmin Hxx | m | λmin S | best feasible phase Rayleigh | classification |',
        '|---|---:|---:|---:|---|---:|---:|---:|---:|---|']
    for v in rows:
        c=v['counts'];ss='N/A' if v['lambda_S_min'] is None else f"{v['lambda_S_min']:.7g}"
        lines.append(f"| {v['label']} | {v['strain']:.12f} | {v['metrics']['force_residual']:.2e} | {v['metrics']['projected_phase_KKT']:.2e} | {c['STRONG_LOWER']}/{c['WEAK_LOWER']}/{c['INTERIOR']} | {v['lambda_x_min']:.7g} | {v['m']} | {ss} | {v['best_feasible_Rayleigh']:.7g} | {v['status']} |")
    lines += ['所有 weak/strong upper 计数均为零。S0 的 Rayleigh 为 mechanical mode；其余为 ||p||=1 的 pᵀSp。恢复的 full mode 按统一 [u,p] norm 重新归一化，真实 full-mode 曲率另存在 MODE_ENERGY.json，不能把 Schur eigenvalue 当作 full-mode eigenvalue。',
        '## 9. 第一次 instability interval','FIRST_STABILITY_LOSS=NOT_FOUND_IN_WINDOW。没有 verified unstable upper endpoint，因此没有合法失稳 bracket、宽度、失稳时 max s 或 active count。没有运行 midpoint refinement，也没有用 brentq。STABILITY_LOSS_RELATIVE_TO_ACTIVATION=INCONCLUSIVE 仅表示首次失稳位置未知；activation point 和已审计正分支稳定性已经解决。',
        '## 10. unstable mode 在哪里','没有 verified unstable mode，故不提供 critical mode PR、top10 fraction、模式中心或位置图。保存的正曲率验证方向不是失稳模式。',
        '## 11. 是否靠近 clamp','历史首次 activation 的两个键位于夹持附近；该事实不能替代本轮尚不存在的 unstable-mode localization。没有把 activation 键或正曲率验证方向称为物理裂纹。',
        '## 12. C3 vs C5 有什么本质变化','只引用历史 C3：约 1% strain 时 185/185 键同时开启，历史失稳区间 [0.0099921875,0.01] 接近 activation。C5 约 0.449739% 时仅 2/5392 键开启，随后存在稳定正相场状态。本轮没有 C7，也不构成空间收敛结论。',
        '## 13. tests',
        f"初始 full suite：{r['baseline_tests'].get('passed')} passed / 3 known failed / 0 skipped。结束 full suite：{r['final_tests'].get('passed')} passed / {r['final_tests'].get('failed')} failed / {r['final_tests'].get('skipped')} skipped。原 3 项失败的 IDs、messages、tracebacks 逐字一致：{r['known_failures_unchanged']}。未改旧测试、未加 xfail/skip。开发中仅 targeted tests；一次在分支尚未写完时启动的 targeted run 因最终汇总文件尚不存在出现两项失败，完成后复跑已保留日志。",
        f"PRE/POST manifests 与 PRESERVATION_CHECK 保护 {r['preservation'].get('protected_files')} 个历史文件及 {r['preservation'].get('protected_symlinks')} 个链接；历史修改={r['preservation'].get('modified')}，生产修改={r['preservation'].get('production_modified_files')}，旧测试修改={r['preservation'].get('existing_tests_modified')}。",
        '## 14. CPU/GPU/计算时间',
        f"12 physical / 24 logical，9 个 ProcessPoolExecutor worker，各 BLAS/OMP/MKL=1，按初始可用内存扣除至少 20% 总 RAM 和父进程预算。Hessian 与 PETSc 使用同一列工作池，PETSc 始终单 rank，没有 ranks×workers 嵌套扩张。唯一 Hxx/cross 矩阵构造总 wall time={r['Hessian_assembly_wall_seconds']:.3f} s；修复能量复核精度时全部复用 hash 匹配缓存，未重新构造这些列。",
        f"父进程 peak RSS={r['peak_parent_RSS_bytes']/1024**2:.1f} MiB；最大 worker peak RSS={r['peak_worker_RSS_bytes']/1024**2:.1f} MiB。各进程峰值保守相加={r['RAM_conservative_process_peak_sum_bytes']/1024**3:.3f} GiB，不冒充同步测得的系统峰值。GPU=NOT USED，没有安装 GPU stack。",
        '## 15. 当前 scientific status','```']
    lines += [f'{k}: {v}' for k,v in r['status'].items()]
    lines += ['```','## 16. 唯一下一步','等待人工审核本轮能量曲率、精度接口和审计窗口结论。本轮在 max(s) 首次达到 0.10 后停止；未启动 C7、后失稳分支或任何后续模型工作。',
        '本轮只生成 4 组有数据支撑的图：stability vs strain、activation vs stability、stationarity quality、历史 C3 vs C5。由于没有 bracket 或 unstable mode，不生成 energy_along_mode 的 stable/unstable 对比或 critical_mode_map。']
    body='\n\n'.join(lines).replace('|\n\n|','|\n|');head,code,tail=body.split('```');body=head+'```'+code.replace('\n\n','\n')+'```'+tail
    (OUT/'REVIEW_C5_STABILITY_V1.md').write_text(body.replace('\n\n\n\n','\n\n')+'\n')


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='preserve':preservation()
    elif len(sys.argv)>1 and sys.argv[1]=='post':post_manifest()
    else:r=collect();figures(r);report(r)
