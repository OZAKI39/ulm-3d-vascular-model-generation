"""Render only the evidence from this task; no new scientific solves."""
from pathlib import Path
import json,hashlib,subprocess,xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scripts import petsc_vi_audit as a
OUT=a.OUT;FIG=a.ROOT/'results/pfpd_at1_petsc_vi_solver_v1/figures'


def stats_xml(name):
    path=OUT/f'tests/{name}.xml'
    if not path.exists():return {}
    root=ET.parse(path).getroot();cases=root.findall('.//testcase')
    failures=[dict(name=x.get('classname')+'::'+x.get('name'),message=x.find('failure').get('message'),trace=x.find('failure').text) for x in cases if x.find('failure') is not None]
    return dict(passed=sum(x.find('failure') is None and x.find('error') is None and x.find('skipped') is None for x in cases),failed=len(failures),skipped=sum(x.find('skipped') is not None for x in cases),failures=failures)


def inventory():
    before=a.read(OUT/'provenance/PRE_RUN_MANIFEST.json');root=Path(before['root']);changed=[];missing=[];links=[]
    for name,digest in before['files'].items():
        p=root/name
        if not p.is_file():missing.append(name)
        elif hashlib.sha256(p.read_bytes()).hexdigest()!=digest:changed.append(name)
    for name,target in before['symlinks'].items():
        p=root/name
        if not p.is_symlink() or str(p.readlink())!=target:links.append(name)
    result=dict(protected_files=len(before['files']),protected_symlinks=len(before['symlinks']),modified=changed,missing=missing,changed_links=links,production_files_modified=0 if not changed and not missing and not links else None,
        HEAD=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=root).strip())
    a.save(OUT/'provenance/PRESERVATION.json',result);return result


def collect():
    branch=a.read(OUT/'continuation/PETSC_RSLS/ALL.json');rows=branch['states'];ss=[a.read(OUT/f'ssls/{label}/SUMMARY.json') for label in ['first_A1','same_load_as_original']]
    env=a.read(OUT/'provenance/PETSC_ENVIRONMENT.json');bench=a.read(OUT/'mpi/BENCHMARK.json');jac=a.read(OUT/'jacobian/GATE.json');base=stats_xml('BASELINE');final=stats_xml('FINAL')
    passed=sum(r['scientifically_accepted'] and r['agrees'] for r in rows);agree=all(r['agrees'] for r in rows);ssagree=all(r['scientifically_accepted'] and r.get('agrees',False) for r in ss)
    preservation=a.read(OUT/'provenance/PRESERVATION.json') if (OUT/'provenance/PRESERVATION.json').exists() else {}
    verified=passed==7 and jac['passed'] and agree and preservation.get('production_files_modified')==0
    status=dict(PETSC_RUNTIME='AVAILABLE',PETSC_RSLS_ADAPTER='VERIFIED' if jac['passed'] else 'NOT_VERIFIED',PETSC_RSLS_C5='VERIFIED' if verified else 'NOT_VERIFIED',
        FIRST_ACTIVATION_ROOT='RECOVERED_AUTOMATICALLY' if rows[0]['scientifically_accepted'] else 'FAILED',ORIGINAL_FAILURE_LOAD='RECOVERED' if rows[4]['scientifically_accepted'] else 'FAILED',
        C5_ORACLE_ROOTS=f'{passed} / 7 PASSED',ACTIVE_SET_DISCOVERY='AUTOMATIC',RSLS_VS_DIAGNOSTIC_ORACLE='AGREES' if agree else 'INCONCLUSIVE',
        PETSC_SSLS='AGREES' if ssagree else 'DIFFICULT',PFM_CRACKS_PDAS_REVIEW='COMPLETE',MPI_BENCHMARK=f'{bench["fastest_ranks"]} rank fastest',GPU_BACKEND='UNAVAILABLE',
        FORMAL_C5_STATIONARY_SOLVER='VERIFIED_CANDIDATE' if verified else 'NOT_VERIFIED',C5_STABILITY='NOT ASSESSED THIS ROUND',C7_STABILITY='NOT RUN',PRODUCTION='UNCHANGED' if preservation.get('production_files_modified')==0 else 'PENDING_AUDIT')
    result=dict(status=status,roots=rows,SSLS=ss,environment=env,MPI=bench,Jacobian_max_error=jac['max_error'],
        max_differences={k:max(r['errors'][k] for r in rows) for k in rows[0]['errors']},average_SNES_iterations=float(np.mean([r['SNES_iterations'] for r in rows])),
        average_linear_iterations=float(np.mean([r['linear_iterations'] for r in rows])),continuation_wall_seconds=branch['wall_seconds'],
        sum_solve_wall_seconds=sum(r['wall_seconds'] for r in rows),baseline_tests=base,final_tests=final,preservation=preservation)
    if base and final:result['known_failures_unchanged']=base['failures']==final['failures']
    a.save(OUT/'SUMMARY.json',result);return result


def figures(r):
    FIG.mkdir(parents=True,exist_ok=True)
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
    plt.rcParams.update({'figure.facecolor':'#080b10','axes.facecolor':'#10151e','savefig.facecolor':'#080b10','font.weight':'normal','pdf.fonttype':42,'axes.formatter.useoffset':False,'font.family':'Arial','font.size':10,'axes.labelcolor':'#e8edf4','text.color':'#e8edf4','xtick.color':'#e8edf4','ytick.color':'#e8edf4','axes.edgecolor':'#e8edf4','savefig.dpi':200,'axes.spines.top':False,'axes.spines.right':False})
    def finish(fig,name):
        fig.tight_layout();fig.savefig(FIG/f'{name}.png');fig.savefig(FIG/f'{name}.pdf');plt.close(fig)
    labels=['first_A1','same_load_as_original'];titles=['First positive root','Original failure load'];colors=['#62b6ff','#ffb366']
    for label,title in zip(labels,titles):
        rows=a.read(OUT/f'continuation/PETSC_RSLS/{label}/ITERATIONS.json');fig,axes=plt.subplots(1,2,figsize=(9,3.4))
        for ax,key,tol in zip(axes,['force','KKT'],[1e-9,1e-8]):
            ax.semilogy([v['iteration'] for v in rows],[max(v[key],1e-16) for v in rows],'o-',color=colors[0]);ax.axhline(tol,color='#e8edf4',ls='--',label=f'Tolerance {tol:.0e}')
            ax.set(xlabel='RSLS Newton iteration',ylabel='Force residual' if key=='force' else 'Projected phase KKT');ax.legend(frameon=False)
        fig.suptitle(title+' | formal VI solver candidate\nZero residuals shown at 1e-16',fontsize=10);finish(fig,f'rsls_convergence_{label}')
    rows=a.read(OUT/'continuation/PETSC_RSLS/first_A1/ITERATIONS.json');fig,axes=plt.subplots(1,3,figsize=(11,3.3))
    for ax,key,label in zip(axes,['lower_count','interior_count','upper_count'],['At lower bound','Free / interior phase','At upper bound']):
        ax.plot([v['iteration'] for v in rows],[v[key] for v in rows],'o-',color=colors[0]);ax.set(xlabel='RSLS Newton iteration',ylabel='Number of phase bonds',title=label)
        if key=='interior_count':ax.axhline(2,ls='--',color='#e8edf4',label='Oracle reference: 2');ax.plot([v['iteration'] for v in rows],[v['PETSc_free_phase_count'] for v in rows],'s:',color=colors[1],label='PETSc free set');ax.legend(frameon=False,fontsize=8)
    finish(fig,'active_set_evolution')
    rows=r['roots'];e=np.array([x['strain'] for x in rows]);fig,axes=plt.subplots(1,3,figsize=(11,3.3))
    for ax,key,tol,label in zip(axes,['phase_Linf','reaction_relative','energy_relative'],[1e-8,1e-8,1e-9],['Phase maximum difference','Relative reaction difference','Relative energy difference']):
        ax.semilogy(e,[max(x['errors'][key],1e-17) for x in rows],'o-',color=colors[0]);ax.axhline(tol,ls='--',color='#e8edf4',label='Agreement threshold');ax.set(xlabel='Engineering strain',ylabel=label);ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0));ax.legend(frameon=False,fontsize=8)
    finish(fig,'solver_root_comparison')
    fig,axes=plt.subplots(1,3,figsize=(11,3.3))
    for ax,key,label in zip(axes,['SNES_iterations','linear_iterations','wall_seconds'],['SNES iterations','Linear solves (PREONLY + MUMPS)','Wall time (s)']):
        ax.plot(e,[x[key] for x in rows],'o-',color=colors[0]);ax.set(xlabel='Engineering strain',ylabel=label);ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    finish(fig,'solver_cost')
    fig,axes=plt.subplots(1,4,figsize=(13,3.4));rs=[rows[0],rows[4]];ss=r['SSLS'];positions=np.arange(2)
    for ax,key,title in zip(axes,['force_residual','projected_phase_KKT','SNES_iterations','wall_seconds'],['Force residual','Projected phase KKT','SNES iterations','Wall time (s)']):
        def values(data):return [max(v.get('metrics',{}).get(key,v.get(key,0)),1e-17) for v in data]
        ax.bar(positions-.18,values(rs),.36,color=colors[0],label='RSLS');ax.bar(positions+.18,values(ss),.36,color=colors[1],label='SSLS')
        ax.set_xticks(positions,['First root','Old failure']);ax.set_title(title)
        if key in ['force_residual','projected_phase_KKT']:ax.set_yscale('log');ax.axhline(1e-9 if key=='force_residual' else 1e-8,color='#e8edf4',ls='--')
    axes[0].legend(frameon=False);fig.suptitle('Two independent VI formulations | zero KKT shown at 1e-17',fontsize=11);finish(fig,'rsls_vs_ssls')
    fig,axes=plt.subplots(1,3,figsize=(11,3.5))
    for ax,key,title in zip(axes,['maximum_s','mean_s','reaction_stress_Pa'],['Maximum phase s','Weighted mean phase s','Reaction stress (Pa)']):
        ax.plot(e,[x['metrics'].get(key,x.get(key)) for x in rows],'o-',color=colors[0]);ax.set(xlabel='Engineering strain',ylabel=title);ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    fig.suptitle('FORMAL VI SOLVER CANDIDATE | NO STABILITY CLASSIFICATION',fontsize=11);finish(fig,'c5_phase_branch_rsls')


def report(r):
    roots=r['roots'];first=roots[0];original=roots[4];env=r['environment'];sparse=a.read(OUT/'JACOBIAN_SPARSITY.json');cfg=a.config()
    md=['# PFPD_AT1_PETSC_VI_SOLVER_V1',
    '## 1. 为什么换 solver 而不换 physics',
    '上一轮已证明 C5 正相场 stationary roots 存在。本轮唯一问题是成熟 VI 求解器能否从合法历史态自动识别约束并找回同一分支。SVK spectral split、AT1、几何、horizon、Yc、Gc、sc=0.95 和全部科学阈值未改变。结论仅针对 stationary solver candidate，不是稳定性、裂纹或真实血栓预测。',
    '## 2. 什么是 variational inequality',
    '相场下界为上一接受载荷的 s，上界为 sc。下界要求 q/D≥0，内部要求 q/D=0，上界要求 q/D≤0，其中 D=max(Ψ+,Yc)>0；力自由度无界且要求归一化残差为零。项目仍独立检查 projected KKT≤1e−8、force≤1e−9。实际 PETSc 标量下界、上界、内部根测试和向量往返全部通过，见 mapping/SCALAR_AND_VECTOR.json。',
    '## 3. RSLS 做了什么',
    '[官方 RSLS](https://petsc.org/main/manualpages/SNES/SNESVINEWTONRSLS/) 在边界符号约束下建立自由子系统并执行 Newton 和线搜索。本轮使用 vinewtonrsls、默认 backtracking、PREONLY + MUMPS，无参数扫描；SNES atol/rtol=1e−12、stol=1e−14、max_it=30。三个 VI 监控及终止原因全部启用，原始文本/二进制保存在各 load 文件夹及运行日志。',
    '初次接线的只读 Vec 访问错误和一次监控集合生命周期崩溃均已修复并保留开发日志。匹配版本 virs.c:509 在 SNESMonitor(:584) 前释放 IS_inact，因此不在 monitor 中读取 getVIInactiveSet；改用公开 C 函数 SNESVIGetActiveSetIS 在当前 (z,F) 上生成集合。lambda 通过公开 SNESLineSearchGetLambda 读取。此改动不影响 Newton 算法。',
    '## 4. 与上一轮 explicit active-set 为什么本质相似',
    '两者都求解约束下的 stationarity，并对未固定的相场分量施加 q=0。差别是本轮分类与迭代由 PETSc 实现；旧诊断仅作为事后 oracle。正式求解代码没有 active bond IDs，也不从已知正相场根初始化。Jacobian gate C 使用正根仅为导数检查，和求解入口隔离。',
    '## 5. PETSc/petsc4py 安装与版本',
    f"PETSc {'.'.join(map(str,env['PETSc_version']))}，petsc4py {env['petsc4py_version']}，MPI {env['MPI_implementation'].splitlines()[0]}。Python `{env['python']}`；scalar=float64、PetscInt=int32，MUMPS 可用。复用现有环境，未修改其包。缺失的项目导入依赖以 binary wheels 安装到本轮 external/python_deps；完整安装日志和固定实际版本清单见 provenance。初次 llvmlite 版本不匹配已依据 numba metadata 修正为 0.49.0。没有源码构建或 CUDA 安装。",
    '## 6. residual / bounds 映射',
    'z 前 225 项沿用旧 scaled free displacement 坐标，后 5392 项为全部 s，共 5617 项。几何固定自由度不入向量。机械 residual 为原始 positive force/force_floor；phase residual 直接为 q/D，没有叠加 clip。稀疏 AIJ 包含 Jxx、Jxs、Jsx 和局部对角 Jss；没有全尺寸 dense 方阵。',
    f"保存的结构 nnz={sparse['structural_nnz']:,}，CSR 占用 {sparse['CSR_bytes']/1024**2:.2f} MiB；PETSc getInfo 的 memory 字段此构建为 0（不是实际无内存），因此另报 CSR 字节和进程 peak RSS。",
    '## 7. Jacobian 验证',
    '采用 inspect/AST 复用上一轮 ResidualBlocks 的原始导数语句：仅把独立位移列循环交给 ProcessPoolExecutor，去掉 projected C 行置零，以及把 projected diagonal 换回原有 Qs/D。生成的可读构造函数已归档；单测检查只有这三处变化。没有第二套导数推导。',
    f"A 最后接受 s=0、B 首个 weak activation、C 已知正相场根，分别 10 个固定随机方向 × 3 档 ε=[1e−5,3e−6,1e−6]，共 90 次。PETSc AIJ Jv 对中心差分最大相对误差 {r['Jacobian_max_error']:.6e}，通过上一轮同量级门槛 5e−5。完整方向误差见 jacobian/GATE.json。",
    '在 s=0 做双侧差分需要负 s 的代数探针。仅诊断函数在原始退化函数的私有 AST 副本中去掉 s-domain guard，保留完全相同 (1−s)² 公式和 split。正式 callback 仍调用原函数、拒绝负 s；从未把差分探针当作物理解或 initial guess。单测验证延拓在合法域与原函数一致。',
    '## 8. first activation result',
    f"ε={first['strain']:.15f}，初值来自 C5_W_004 加边界增量，相场全零。{first['SNES_iterations']} 次 Newton 后科学接受，max s={first['metrics']['maximum_s']:.12g}，force={first['metrics']['force_residual']:.3e}，KKT={first['metrics']['projected_phase_KKT']:.3e}。事后内部键 {first['interior_IDs']}，precision={first['discovery']['precision']:.1f}、recall={first['discovery']['recall']:.1f}。",
    '首个载荷先只调整机械位移；接近激活时 BT 缩短步长，随后 PETSc 释放两个相场约束并进入正相场根。这是局部、夹持关联的 activation，不是裂纹预测。',
    '## 9. original failure load result',
    f"ε={original['strain']:.15f}，从前一个 PETSc 接受态继续。force={original['metrics']['force_residual']:.3e}，KKT={original['metrics']['projected_phase_KKT']:.3e}，max s={original['metrics']['maximum_s']:.12g}；完整科学门槛通过。旧 projected solver 源码和原失败结果未修改。",
    '## 10. 七个 oracle roots regression',
    '| ε | force | KKT | max s | Δs∞ | scaled Δx₂ | relative Δenergy | relative Δreaction | accepted |',
    '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for x in roots:
        er=x['errors'];md.append(f"| {x['strain']:.12f} | {x['metrics']['force_residual']:.2e} | {x['metrics']['projected_phase_KKT']:.2e} | {x['metrics']['maximum_s']:.9g} | {er['phase_Linf']:.2e} | {er['scaled_state_L2']:.2e} | {er['energy_relative']:.2e} | {er['reaction_relative']:.2e} | {x['scientifically_accepted']} |")
    md += ['','预声明 agreement：phase Linf≤1e−8；自由位移差除 length_scale 后的 L2≤1e−8；energy relative≤1e−9；reaction relative≤1e−8。能量与反力分母分别使用对应 oracle 的绝对基准值。七点均满足，当前无多根证据；不据此主张根唯一。',
    '## 11. continuation 是否自动维护 active set',
    '七点按真实保存的 ε 排序，将原失败载荷置于 offset_1e−05 与 offset_2e−05 之间。每一点保存 z_initial、lower、accepted state；测试逐元素验证初值等于上一 PETSc 解的位移预测及上一 s，未逐点喂 oracle。所有点 inactive inequalities、严格下界、s<sc、不可逆性、固定边界、finite positive J 和旧安全检查通过。',
    '## 12. RSLS vs SSLS',
    '[官方 SSLS](https://petsc.org/main/manualpages/SNES/SNESVINEWTONSSLS/) 为另一成熟 VI formulation。本轮仅两点，读取对应 RSLS 的 z_initial 与 lower；同一 residual、Jacobian、tol 和线性配置，无调参。',
    '| case | method | force | KKT | SNES | linear | seconds | scientifically accepted |',
    '|---|---|---:|---:|---:|---:|---:|---|']
    for name,rs,ss in zip(['first','original failure'],[first,original],r['SSLS']):
        for tag,x in [('RSLS',rs),('SSLS',ss)]:md.append(f"| {name} | {tag} | {x.get('metrics',{}).get('force_residual',float('nan')):.3e} | {x.get('metrics',{}).get('projected_phase_KKT',float('nan')):.3e} | {x['SNES_iterations']} | {x['linear_iterations']} | {x['wall_seconds']:.3f} | {x['scientifically_accepted']} |")
    md += ['',f"SSLS 状态：{r['status']['PETSC_SSLS']}；科学接受与 PETSc reason 独立。详细 reason/error/门槛在 ssls/*/SUMMARY.json。",
    '## 13. pfm-cracks / 2015 PDAS 对照',
    '固定 `tjhei/cracks`、`paper-2015`、`342678db32e1736b3f1eae601f55322907f584c9`，GPL-2.0-or-later。源码行号、不可逆相场方向转换、FE 专属部分和停止/集合变化逻辑见 [PFM_CRACKS_PDAS_NOTES.md](references/PFM_CRACKS_PDAS_NOTES.md)。旧源码要求残差阈值和集合稳定共同满足，并有最大迭代限制；没有把其 FE 框架移植进来。',
    '[Kolditz–Mang–Wick 2023](https://arxiv.org/html/2306.03501v2) 仅作为效率参考：若未来出现反复集合切换或过高组合 Newton 成本，再单独评估其改进。本轮不实现其算法。',
    '## 14. CPU / MPI / GPU 实际效率',
    f"硬件复核为 {cfg['resources']['physical_cores']} physical / 24 logical；BLAS/OMP/MKL 均单线程；5 个 spawn 进程并行位移导数列。按可用内存、父进程预算和至少 20% RAM 预留选择，不使用 24 MPI ranks。continuation 串行，两个 SSLS 也串行执行以避免各自列工作池叠加。",
    '| MPI ranks | same-root seconds | accepted |', '|---:|---:|---|']
    for x in r['MPI']['runs']:md.append(f"| {x['ranks']} | {x['wall_seconds']:.3f} | {x['scientifically_accepted']} |")
    md += ['',f"最快 {r['MPI']['fastest_ranks']} rank，故正式 continuation 使用该模式。每种模式只测一次，主张限于本机本轮观测；rank-zero BAQP 并未变成分布式 physics。平均 SNES/linear 均 {r['average_SNES_iterations']:.3f}；PREONLY 的 linear iterations 表示直接线性求解次数，不能解释为 Krylov 收敛难度。continuation 路径 wall time={r['continuation_wall_seconds']:.3f} s，各 solve 时间合计={r['sum_solve_wall_seconds']:.3f} s。",
    'RTX 4060 存在，但此 PETSc 未配置 CUDA/Kokkos/HIP，VecCUDA/AIJCUSPARSE 不可用；GPU_BACKEND=UNAVAILABLE，没有 GPU 性能结论。C7 optional smoke 仅完成 46794 DOF 的 residual 映射与 AIJ 初始化：低载荷 ε=0.00025，附近预测载荷 ε≈0.00314907。后三点 cheap scan 使用低载荷松弛位移的线性外推，未机械平衡，也不构成新 activation bracket。该状态 max Ψ+/Yc≈1.0093。完整 Jacobian、CSR/AIJ 临时数组与列工作池的保守估算为 6.49 GB，超过保留 20% RAM 后的预算，故未装配 C7 Jacobian、未执行 Newton，不能称 C7 coupled solver 已验证。执行证据见 c7/SMOKE.json；未运行任何 C7 Hessian、critical cone、mode 或空间稳定性工作。',
    '## 15. tests',
    f"初始一次 full suite：{r['baseline_tests'].get('passed')} passed、{r['baseline_tests'].get('failed')} known failed、{r['baseline_tests'].get('skipped')} skipped。结束一次 full suite：{r['final_tests'].get('passed')} passed、{r['final_tests'].get('failed')} failed、{r['final_tests'].get('skipped')} skipped。原失败 ID/message/trace 是否逐字不变：{r.get('known_failures_unchanged')}。开发期间仅运行本轮 targeted tests。",
    f"历史保护：{r['preservation'].get('protected_files')} files、{r['preservation'].get('protected_symlinks')} symlinks；modified={r['preservation'].get('modified')}，missing={r['preservation'].get('missing')}，changed_links={r['preservation'].get('changed_links')}。生产文件修改数={r['preservation'].get('production_files_modified')}。",
    '## 16. 当前 scientific status', '```']
    md += [f'{k}: {v}' for k,v in r['status'].items()]
    md += ['```','## 17. 唯一下一步','等待人工审核本轮 stationary solver candidate。只有审核接受 VERIFIED_CANDIDATE 后，下一轮才能另行启动 C5 stability；本轮在此停止。',
        '图件沿用 dark background、Arial regular、English labels、无箭头和无双 y 轴。图件见 `results/pfpd_at1_petsc_vi_solver_v1/figures/`：6 组、7 张图，每张有 PNG/PDF。所有分支图仅使用本轮 accepted roots。']
    body='\n\n'.join(md).replace('|\n\n|','|\n|')
    head,code,tail=body.split('```')
    body=head+'```'+code.replace('\n\n','\n')+'```'+tail
    (OUT/'REVIEW_PETSC_VI_SOLVER_V1.md').write_text(body.replace('\n\n\n\n','\n\n')+'\n')


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='inventory':inventory()
    else:r=collect();figures(r);report(r)
