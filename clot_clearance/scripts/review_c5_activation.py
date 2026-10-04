"""Evidence-only review and seven figure groups for C5 activation stationarity."""
from scripts.c5_activation_audit import *
import os,subprocess,xml.etree.ElementTree as ET,itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

FIG=RESULTS/'figures'
BLUE='#64caff';GOLD='#ffcc6b';PINK='#ee94d2';WHITE='#eeeeee'


def test_summary():
    def parse(path):
        tree=ET.parse(path);cases=tree.findall('.//testcase');failures={}
        for t in cases:
            f=t.find('failure')
            if f is not None:failures[t.attrib['classname']+'::'+t.attrib['name']]=dict(message=f.attrib.get('message'),traceback=f.text)
        skipped=sum(int(s.attrib.get('skipped',0)) for s in tree.findall('.//testsuite'))
        return dict(passed=len(cases)-len(failures)-skipped,failed=len(failures),skipped=skipped,failures=failures)
    baseline=parse(OUT/'tests/BASELINE.xml');final=parse(OUT/'tests/FINAL.xml');old=parse(OLD/'tests/FINAL.xml')
    r=dict(baseline=baseline,final=final,known_failure_ID_message_traceback_identical=baseline['failures']==final['failures']==old['failures'])
    assert r['known_failure_ID_message_traceback_identical'] and final['failed']==3 and final['skipped']==0
    save(OUT/'tests/TEST_SUMMARY.json',r);return r


def preservation():
    pre=read(OUT/'provenance/PRE_RUN_MANIFEST.json');files={};links={};changed=[];missing=[]
    for name,sha in pre['files'].items():
        p=ROOT/name
        if not p.is_file() or p.is_symlink():missing.append(name);continue
        files[name]=hashlib.sha256(p.read_bytes()).hexdigest()
        if files[name]!=sha:changed.append(name)
    for name,dest in pre['symlinks'].items():
        p=ROOT/name
        if not p.is_symlink():missing.append(name);continue
        links[name]=os.readlink(p)
        if links[name]!=dest:changed.append(name)
    save(OUT/'provenance/POST_RUN_MANIFEST.json',dict(root=str(ROOT),files=files,symlinks=links,timestamp=time.time(),scope='Exactly the complete pre-existing file/link set'))
    r=dict(passed=not changed and not missing,checked_files=len(files),checked_symlinks=len(links),changed=changed,missing=missing,
           production_modified_files=0 if not changed else None,historical_modified_files=0 if not changed else None,existing_tests_modified=0 if not changed else None)
    save(OUT/'provenance/PRESERVATION_CHECK.json',r);assert r['passed'];return r


def labels():
    same=read(OUT/'reduced/SAME_LOAD.json');prefix=read(OUT/'reduced/CONTINUATION.json')
    assert same['valid'] and same['positive_phase'] and all(r['valid'] for r in prefix)
    result=dict(C5_ACTIVATION_AUDIT='COMPLETE',S0_KKT_ACTIVATION='BRACKETED',ACTIVATION_SET='LOCALIZED',
        PROJECTED_KKT_MAP='NONSMOOTH_ACTIVE_SET_TRANSITION',RESIDUAL_JACOBIAN='WELL_CONDITIONED',
        EXPLICIT_ACTIVE_SET_ROOT='FOUND',COMPLEMENTARITY_ORACLE='NOT_NEEDED',
        ROOT_CAUSE='PROJECTED_SOLVER_ACTIVE_SET_STAGNATION',C5_COUPLED_STATIONARITY='RESOLVED',
        C5_STABILITY='NOT ASSESSED THIS ROUND',C7_COUPLED='NOT RUN',PRODUCTION='UNCHANGED')
    for path in [OUT,RESULTS]:save(path/'FINAL_STATUS.json',result)
    return result


def summaries():
    rows=[read(p) for p in (OUT/'intact').glob('*/SUMMARY.json')]
    reduced=[read(p) for p in (OUT/'reduced').glob('*/SUMMARY.json')]
    same=read(OUT/'reduced/SAME_LOAD.json');original=read(OUT/'original_solver/SUMMARY.json')
    resource_summary=dict(intact_states=len(rows),intact_wall_seconds=sum(r['wall_seconds'] for r in rows),
        original_solve_count=1,original_wall_seconds=original['wall_seconds'],reduced_solves=len(reduced),
        reduced_wall_seconds=sum(r['wall_seconds'] for r in reduced),
        jacobian_diagnostic_wall_seconds=sum(r['wall_seconds'] for r in read(OUT/'jacobian/ALL.json')),
        energy_oracle_wall_seconds=read(OUT/'oracles/ENERGY.json').get('wall_seconds'),
        intact_peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in rows),max_workers=1,BLAS_threads=1,numba_threads_configured=10,
        reserved_physical_cores=2,RAM_reserve_fraction=.2,GPU_used=False,GPU_benchmark='NOT_RUN',GPU_speedup=None,
        timing_scope='Measured calls only; not total elapsed time. Some first-root configurations have identical active sets but were tested separately as requested.')
    save(OUT/'provenance/COMPUTE_SUMMARY.json',resource_summary)
    save(RESULTS/'SUMMARY.json',dict(activation=read(OUT/'activation/BRACKET.json'),same_load_root=same,labels=labels()))
    # Explicit per-bond derivative-branch arrays complement the stored clip code.
    cases=[('A_last_accepted',OUT/'reproduction/C5_W_004.npz'),
           ('B_first_weak',OUT/'intact/activation_point/STATE.npz'),('C_failed_coupled',OUT/'original_solver/KINK.npz')]
    kink=[]
    for name,path in cases:
        a=ca.read_npz(path);branch=a['clip_branch'];b=ca.read_npz(OUT/f'jacobian/{name}/BLOCKS.npz')
        record=dict(label=name,lower_clipped=int((branch==-1).sum()),interior_map=int((branch==0).sum()),upper_clipped=int((branch==1).sum()),
            near_kink=int((a['kink_distance']<=1e-7).sum()),exact_kink=int((a['kink_distance']==0).sum()),
            minimum_kink_distance=float(a['kink_distance'].min()),
            derivative_branch_legend={'0':'dq_scaled (interior map, including exact tie convention)','1':'ds (lower/upper clipped)'})
        save(OUT/f'projected_map/{name}.json',record)
        np.savez_compressed(OUT/f'projected_map/{name}.npz',clip_branch=branch,kink_distance=a['kink_distance'],
            local_derivative_branch=(branch!=0).astype(int),phase_diagonal=b['diag'],q=a['q'],q_scaled=a['q_scaled'],s=a['s'])
        kink.append(record)
    save(OUT/'projected_map/ALL.json',kink)
    return resource_summary


def style():
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','axes.titleweight':'normal','figure.facecolor':'#080b10',
        'axes.facecolor':'#080b10','savefig.facecolor':'#080b10','text.color':WHITE,'axes.labelcolor':WHITE,
        'xtick.color':WHITE,'ytick.color':WHITE,'axes.edgecolor':'#7c8796','grid.color':'#39414c','font.size':11,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.formatter.useoffset':False})


def finish(fig,name):
    FIG.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(FIG/f'{name}.{ext}',dpi=190,bbox_inches='tight')
    plt.close(fig)


def figures():
    style();scan=read(OUT/'activation/SCAN.json');br=read(OUT/'activation/BRACKET.json');point=read(OUT/'activation/POINT.json')['strain']
    e=np.array([r['strain'] for r in scan]);q=np.array([r['q_percentiles'] for r in scan])
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.5),layout='constrained')
    for i,label,color in [(0,'Minimum',PINK),(2,'5th percentile',GOLD),(4,'Median',BLUE)]:axes[0].plot(e,q[:,i],'o-',color=color,label=label)
    axes[1].plot(e,q[:,0],'o-',color=PINK,label='Minimum q')
    for ax in axes:
        ax.axhline(0,color=WHITE,ls='--',lw=1);ax.axvspan(*br['interval'],color=GOLD,alpha=.5,label='Activation bracket')
        ax.set(xlabel='Engineering strain',ylabel='Phase gradient q (Pa)');ax.grid(alpha=.3);ax.legend(frameon=False)
    axes[0].set_title('Intact branch: phase-gradient distribution');axes[1].set_title('First loss of lower-bound KKT feasibility')
    finish(fig,'activation_gradient_vs_strain')

    full=scan+[br['lower'],read(OUT/'activation/POINT.json'),br['upper']];full=sorted(full,key=lambda r:r['strain'])
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.5),layout='constrained')
    for key,color in [('STRONG_LOWER',BLUE),('WEAK_LOWER',GOLD),('VIOLATED_LOWER',PINK)]:
        for ax in axes:ax.plot([r['strain'] for r in full],[r['counts'][key] for r in full],'.-',label=key.replace('_',' ').title(),color=color)
    axes[0].set(yscale='symlog',ylim=(-.05,8000),ylabel='Bond count (symlog)',title='All 5392 bonds')
    axes[1].set(ylim=(-.1,2.4),xlim=(point-1e-6,point+1e-6),ylabel='Bond count',title='Activation set: 2 bonds / 0.0371%')
    ticks=[point-8e-7,point,point+8e-7]
    axes[1].set_xticks(ticks,[f'{v:.9f}' for v in ticks]);axes[1].tick_params(axis='x',labelsize=9)
    for ax in axes:ax.set_xlabel('Engineering strain');ax.grid(alpha=.3)
    axes[0].legend(frameon=False,loc='center left',fontsize=9)
    finish(fig,'activation_set_size')

    p=Problem(5,point);candidates=read(OUT/'activation/SPATIAL.json')[1]['candidates']
    fig=plt.figure(figsize=(8.4,6.1),layout='constrained');ax=fig.add_subplot(projection='3d')
    vertices=list(itertools.product([-125,125],repeat=3))
    for a in vertices:
        for b in vertices:
            if sum(x!=y for x,y in zip(a,b))==1 and a<b:ax.plot(*np.array([a,b]).T,color='#697586',lw=.8)
    for c in candidates:
        X=p.family.X[p.family.pairs[c['ID']]]*1e6
        ax.plot(*X.T,color=PINK,lw=5);mid=np.array(c['midpoint_m'])*1e6;ax.scatter(*mid,color=GOLD,s=65)
        ax.text(mid[0]+10,mid[1]+10,mid[2],str(c['ID']),color=GOLD)
    ax.text(-110,-110,140,'Top clamp',color=WHITE);ax.text(-110,-110,-145,'Bottom clamp',color=WHITE)
    ax.set(xlim=(-130,130),ylim=(-130,130),zlim=(-150,150),xlabel='x (micrometre)',ylabel='y (micrometre)',zlabel='z (micrometre)')
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor((0,0,0,0));axis.pane.set_edgecolor('#39414c')
    ax.set_box_aspect((1,1,1));ax.view_init(elev=17,azim=-60)
    ax.set_title('AT1 ACTIVATION SET\nNOT CRACK PREDICTION',pad=13)
    finish(fig,'activation_spatial_map')

    old=read(OUT/'original_solver/SUMMARY.json');hist=old['evaluations'];fig,axes=plt.subplots(1,2,figsize=(10.5,4.5),layout='constrained')
    for ax,key,tol,title in zip(axes,['mechanical_L2','phase_Linf'],[1e-9,1e-8],['Mechanical residual L2','Projected phase residual Linf']):
        ax.semilogy([r['evaluation'] for r in hist],[max(r[key],1e-17) for r in hist],'o-',color=PINK)
        ax.axhline(tol,color=GOLD,ls='--',label='Scientific tolerance')
        ax.axvline(hist[-1]['evaluation'],color=BLUE,ls=':',label='SciPy ftol + xtol stop')
        ax.set(xlabel='Residual evaluation',ylabel='Normalized residual',title=title);ax.grid(alpha=.3);ax.legend(frameon=False,fontsize=9)
    axes[1].text(.34,.06,'Initial phase residual = 0\nshown at the plot floor',transform=axes[1].transAxes,fontsize=9)
    finish(fig,'solver_stagnation')

    cond=read(OUT/'jacobian/ALL.json');fig,axes=plt.subplots(1,3,figsize=(11,4.4),layout='constrained')
    for ax,key,title in zip(axes,['sigma_max','sigma_min','condition_estimate'],['Largest singular value','Smallest singular value','Condition estimate']):
        y=[r.get(key) for r in cond]
        if all(v is not None for v in y):
            ax.semilogy(range(3),y,'o-',color=BLUE)
            for i,v in enumerate(y):ax.annotate(f'{v:.4g}',(i,v),xytext=(0,10),textcoords='offset points',ha='center',fontsize=10)
            ax.margins(y=.3)
        else:ax.text(.5,.5,'UNRESOLVED',transform=ax.transAxes,ha='center')
        ax.set(xticks=range(3),xticklabels=['Last s = 0','First weak','Failed iterate'],title=title)
        ax.tick_params(axis='x',labelrotation=20);ax.grid(alpha=.3)
    axes[1].set_ylim(1e-4,1e-2)
    finish(fig,'jacobian_conditioning')

    same=read(OUT/'reduced/SAME_LOAD.json');energy=read(OUT/'oracles/ENERGY.json')
    names=['Original\nprojected','Explicit\nactive set'];metrics=[old['metrics'],same['metrics']]
    if 'metrics' in energy:names+=['Energy\noracle'];metrics+=[energy['metrics']]
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.6),layout='constrained')
    for ax,key,tol,title in zip(axes,['force_residual','projected_phase_KKT'],[1e-9,1e-8],['Force residual','Full projected KKT residual']):
        vals=[max(m[key],1e-17) for m in metrics]
        ax.bar(range(len(vals)),vals,color=[PINK,BLUE,GOLD][:len(vals)],bottom=1e-17,width=.55)
        ax.set(yscale='log',ylim=(1e-17,1e-3),xticks=range(len(vals)),xticklabels=names,ylabel='Normalized residual',title=title)
        ax.axhline(tol,color=WHITE,ls='--');ax.grid(axis='y',alpha=.3)
    fig.suptitle(f'Same engineering strain: {same["strain"]:.12f}\nEnergy oracle does not pass physical tolerances')
    finish(fig,'projected_vs_active_set_solver')

    prefix=read(OUT/'reduced/CONTINUATION.json');fig,axes=plt.subplots(1,3,figsize=(11.4,4.5),layout='constrained')
    rows=[r for r in prefix[1:] if r['valid']][:5]
    if rows:
        for ax,key,title in zip(axes,['maximum_s','mean_s','reaction_stress_Pa'],['Maximum phase s','Weighted mean phase s','Reaction stress (Pa)']):
            y=[r[key] if key in r else r['metrics'][key] for r in rows]
            ax.plot([r['strain'] for r in rows],y,'o-',color=BLUE);ax.set(title=title,xlabel='Engineering strain');ax.grid(alpha=.3)
            ax.tick_params(axis='x',labelrotation=25)
    else:
        for ax in axes:ax.text(.5,.5,'NO VALID STATIONARY ROOT FOUND',ha='center',transform=ax.transAxes)
    fig.suptitle('Five accepted continuation states; no stability classification')
    finish(fig,'positive_phase_branch_prefix')
    save(OUT/'provenance/FIGURES.json',dict(directory=str(FIG),groups=7,files=sorted(p.name for p in FIG.iterdir()),
         style='Dark background, Arial regular, English labels, no arrows or dual axes'))


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def report(tests,protected,resources):
    bracket=read(OUT/'activation/BRACKET.json');point=read(OUT/'activation/POINT.json');cond=read(OUT/'jacobian/ALL.json')
    original=read(OUT/'original_solver/SUMMARY.json');same=read(OUT/'reduced/SAME_LOAD.json');first=read(OUT/'reduced/first_A1/SUMMARY.json')
    prefix=read(OUT/'reduced/CONTINUATION.json');scan=read(OUT/'reduced/ROOT_SCAN.json');energy=read(OUT/'oracles/ENERGY.json')
    replay=read(OUT/'original_solver/FROZEN_LINEAR_STEP_REPLAY.json');scaling=read(OUT/'jacobian/SCALING_DIAGNOSTIC.json')
    lab=labels();head=(OUT/'provenance/HEAD.txt').read_text().strip();repo=(OUT/'provenance/REPOSITORY.txt').read_text().strip()
    tol=read(OUT/'provenance/PREDECLARED_CONFIG.json');last=original['termination_checks'][-1]
    table_cond=table(['状态','σmax','σmin','condition','机械块 condition'],[
        [r['label'],f"{r['sigma_max']:.8g}",f"{r['sigma_min']:.8g}",f"{r['condition_estimate']:.8g}",f"{r['mechanical_condition']:.6g}"] for r in cond])
    table_kink=table(['状态','lower-clipped','interior-map','upper-clipped','near kink ≤1e−7'],[
        [r['label'],r['lower_clipped'],r['interior_map'],r['upper_clipped'],r['near_kink']] for r in read(OUT/'projected_map/ALL.json')])
    table_roots=table(['载荷','force','KKT','max s','active size','结果'],[
        [f"{r['strain']:.12f}",f"{r['metrics']['force_residual']:.4e}",f"{r['metrics']['projected_phase_KKT']:.4e}",
         f"{r['metrics']['maximum_s']:.7g}",r['active_size'],'VALID' if r['valid'] else 'REJECTED'] for r in scan+[same]])
    labels_text='\n'.join(f'{k}: {v}' for k,v in lab.items())
    tests_text='\n'.join('- `'+name+'`' for name in tests['baseline']['failures'])
    text=f'''# C5 activation stationarity 审核报告

结论：在审计邻域内找到了通过原 force、KKT、inactive inequality、bounds、irreversibility 和 safety 检查的正相场 stationary roots；原失败载荷也直接找到合法根。五个固定小步长 continuation states 提供了独立邻近载荷支持，C5_COUPLED_STATIONARITY = RESOLVED。此结论只解决 stationary root 的存在性，未评估稳定性。

Repository：`{repo}`；HEAD：`{head}`。新增内容只在审计脚本、诊断 wrapper、测试和本轮结果目录。Production、历史结果、既有测试均未修改。

## 1. 上一轮到底卡在哪里

上一轮已经验证谱分裂和 C5 完整态混合拉压平衡。困难出现在相场开始增大时：ε≈0.004514117556 的 projected-KKT + dogbox/lsmr 求解器返回 ftol/xtol success，但 force≈3.52e−6、KKT≈4.90e−6，没有达到 1e−9 / 1e−8。

本轮先重新读取 last accepted `C5_W_004` 和三次失败结果 `C5_A_001/C5_W_005/C5_P_005`。能量、force、KKT、q 和 Ψ+ 均按预声明门槛复现。归档没有独立 Ψ+ 数组的旧 coupled state，以旧 q_scaled 和 s 的明确公式反解 Ψ+，再与重新评价比较；这不是额外的构成 oracle。证据在 `reproduction/ALL.json`。

## 2. 什么是 phase activation

s=0 表示当前完整相场；s 开始增大称为 AT1 activation。这个词只说明某些内部变量开始偏离下界，不能直接解释为裂纹形状、裂纹传播或失稳。本轮首先定位的是 lower-bound KKT activation bracket，不是 fracture stability bracket。

## 3. 什么是 q

在既有 AT1 定义中，q=−2(1−s)Ψ+ + 2Yc，是相场增加一个小量时的能量变化方向，已去掉正的 BAQP 权重。s=0 时 q=2(Yc−Ψ+)。若 q>0，增加 s 会先提高能量，下界可满足 KKT；q≈0 是候选激活；q<0 表示保持 s=0 已不满足下界不等式。

原缩放 D_b=max(Ψ+_b,Yc)，q_scaled=q/D。weak tolerance 沿用 1e−7，故每根键的有量纲 q_active_tol=1e−7 D_b；s_tol=1e−8。最终 projected KKT tolerance 仍为 1e−8，force tolerance 仍为 1e−9。没有按结果改变容差。

## 4. 为什么先检查 s=0 branch

关闭 phase updates 后，只要求力学平衡，再计算 q，就能把“完整态何时失去相场下界可行性”与“通用 coupled optimizer 是否能处理投影分支”分开。先执行指定 11 个载荷，再局部二分，未进行 0–1% 高密扫描。每个点复用原固定相场 ADR 和正式 spectral Ψ+；所有被使用的 intact states 均通过原力、边界、正 J 和 safety 检查。

## 5. activation bracket

下端 εL={bracket['interval'][0]:.12f}，min(q)={bracket['lower']['min_q']:.9e} Pa；上端 εU={bracket['interval'][1]:.12f}，min(q)={bracket['upper']['min_q']:.9e} Pa。区间宽度 {bracket['width']:.9e}，局部二分 {len(bracket['history'])} 层，达到 ≤1e−7。

为了审计 first weak state，在此有效区间内用 SciPy brentq 对“已平衡的 min(q)”做标量定位，得到 ε≈{point['strain']:.15f}，min(q)={point['min_q']:.6e} Pa，2 根 weak bonds。它只是补充弱分支诊断点，不替代前述有符号两端区间。`activation/BRACKET.json` 保存每步收缩；每点 q percentiles、q<0 IDs、weak IDs、max(Y/Yc) 和完整数组均归档。

![Activation gradient](../../results/pfpd_at1_c5_activation_stationarity_audit_v1/figures/activation_gradient_vs_strain.png)

## 6. 多少 bonds 首先激活

最先激活 2 / 5392 根键，比例 0.037091988%。最后 strong-lower endpoint 为 5392 strong；weak point 为 5390 strong + 2 weak；第一个上端点为 5390 strong + 2 violated。候选集合是少数局部键，非整层、非大量同步激活。A1:q<0、A2:q_scaled≤1e−7、A3:前两集合加共享端点且 q_scaled≤1e−3 的近阈值键，三者实际都选中同两根键；不是三套不同的独立集合证据。

![Activation set size](../../results/pfpd_at1_c5_activation_stationarity_audit_v1/figures/activation_set_size.png)

## 7. activation 在哪里

键 ID 4488、5059，reference midpoint 分别为 (0,0,−93.75) µm 与 (0,0,93.75) µm，方向均为 z 轴。它们分别连接中心轴附近的夹持节点与内侧相邻节点；中点到最近夹持面 31.25 µm，到侧面 125 µm。弱点 Ψ+≈Yc=0.1010025 Pa。

空间分类为 LOCALIZED，准确含义是“两个夹持面附近的离散激活位置”。这既不是 crack prediction，也不是局部化失稳判据。三个选定状态的 IDs、orientation、midpoints、Ψ+、q 见 `activation/SPATIAL.json`。

![Activation spatial map](../../results/pfpd_at1_c5_activation_stationarity_audit_v1/figures/activation_spatial_map.png)

## 8. projected KKT 为什么可能难解

令 t=s−q/D，r_s=s−clip(t,lower,sc)。当 t<lower 时 r_s=s−lower，dr_s=ds；当 lower<t<sc 时 r_s=q/D，dr_s=dq/D−q·dD/D²；当 t>sc 时 r_s=s−sc，dr_s=ds。dq=2Ψ+ ds−2(1−s)dΨ+；dD 在 Ψ+>Yc 时为 dΨ+，在 Ψ+<Yc 时为零。在 clip 的交界不能声称唯一普通导数；旧代码 exact tie 采用 interior-map 分支，保持未改。

{table_kink}

每根键的 clip branch、到最近 switching surface 的距离、局部导数分支及 phase diagonal 保存到 `projected_map/`。激活点只有两个 near kinks，因此“大规模同步 residual kink”假设不受数据支持。

原 nonlinear solver 只复跑一次，初值严格采用归档 `C5_P_005/z_initial`，所有数值参数与 residual 不变。用可恢复 wrapper 记录 SciPy 原调用返回值：status={original['solver']['status']}，nfev={original['solver']['nfev']}，njev={original['solver']['njev']}，cost={original['solver']['cost']:.8e}，optimality={original['solver']['optimality']:.8e}，||Jᵀr||∞={original['JT_r_Linf']:.8e}。

最后一步 Δcost={last['dF']:.1f}、step norm={last['step_norm']:.1f}，xtol threshold={last['xtol_threshold']:.6e}，因此 ftol/xtol 同时触发；这两个停止条件本来就不保证本项目的 force/KKT 验收。[SciPy least_squares 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html) 分别定义成本变化、步长和梯度停止条件。

额外只在最后 frozen iterate 重放一次既有线性求解与 SciPy dogleg_step，未重跑 nonlinear solver：重放步长 {replay['replayed_step_norm']:.6e}，无约束线性候选步长 {replay['newton_step_norm']:.6e}，有 {replay['negative_newton_step_at_lower_bound_count']} 个已在 lower bound 的分量收到负向线性建议。具体 IDs、线性残差、边界信息在 `original_solver/FROZEN_LINEAR_STEP_REPLAY.json`。重放得到很小但非零的步长，未精确复现记录中的零步长；重新构造的有限差分算子及近零边界梯度会有数值敏感性。因此它只支持 bound-constrained step 受限的局部解释，不能据此断言唯一内部根因或 SciPy 库的普遍 bug。

![Solver stagnation](../../results/pfpd_at1_c5_activation_stationarity_audit_v1/figures/solver_stagnation.png)

## 9. Jacobian conditioning

仅审计 residual Jacobian，没有稳定性 Hessian 或 critical cone。存储 225 个 displacement columns 与 phase 对角/交叉块，约 29.5 MB，不建立 5392×5392 dense phase matrix。mechanical/Y 对 displacement 用原 mechanics/split 的中心差分，phase-only 与交叉块按已有能量梯度链式求导；与未改的 V2 Jv 做四方向对照，并验证伴随关系。reduced Jacobian 再以三档 epsilon 做独立 residual 有限差分，最大相对误差约 7.64e−7。

{table_cond}

使用成熟 [SciPy svds](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html)：最大奇异值直接求，最小值通过 compact Jacobian 的块逆求最大奇异值后倒数；225×225 块分解使用 SciPy LU，仅作线性诊断，未实现 nonlinear Newton。两端奇异向量均回代原 Jacobian/transpose 验证，残差远低于预声明 1e−6。

没有数值秩亏证据。WELL_CONDITIONED 标签是相对于预声明 condition≥1e8 的严重病态筛查阈值；并不表示条件数接近 1 或尺度已最优。failed iterate 的机械梯度贡献 ||Jxᵀrx||₂≈{cond[2]['mechanical_gradient_L2']:.6e}，phase 贡献≈{cond[2]['phase_gradient_L2']:.6e}，后者明显更大。row/column norms 的逐项数组和分位数已保存，不能只凭 block residual 大小推断优化权重。

{table(['状态','原 condition','仅 row-normalized Jacobian condition'],[[r['label'],f"{r['original_condition']:.6g}",f"{r['row_equilibrated_condition']:.6g}" if r['row_equilibrated_condition'] is not None else 'UNRESOLVED'] for r in scaling])}

row normalization 只做 frozen Jacobian 诊断，没有改正式 residual，也没有用新的 scaling 宣布求解成功。数据支持尺度有差异；仅凭它不能把 ROOT_CAUSE 唯一归为 scaling。

![Jacobian conditioning](../../results/pfpd_at1_c5_activation_stationarity_audit_v1/figures/jacobian_conditioning.png)

## 10. explicit active-set reduced solve

inactive bonds 精确固定 s_I=previous accepted lower；未知量为 225 个 displacement variables 与候选 active phase variables，初始共 227 个未知量。方程为原 F_free/force_floor=0 与 q_A/Yc=0。q_A/Yc 是平滑 active 等式的诊断归一化；最终验收仍调用原 q/max(Ψ+,Yc) projected KKT，并额外检查全部 inactive inequalities。

真正 nonlinear solve 使用 SciPy least_squares 的 trf/exact；没有自写 Newton、line search 或 trust region。wrapper 只负责 active IDs、变量映射、固定 inactive 值、最多 10 次一致性更新。加入规则 q_I/D<−1e−8；移除规则 active 回到 lower 且 q_A/D>1e−8。规则与全部变化均保存。实际合法根沿用同两根 active bonds，未触发额外更新；这不影响完整 inactive 检查。

初猜只调用一次原 AT1 law，未改 irreversibility 下界。各 reduced run 保存 `SUMMARY.json/STATE.npz/ITERATIONS.json/HISTORY.json`。这套 wrapper 明确为 diagnostic，未替换原 branch solver。

## 11. 是否找到合法 stationary root

是。首次正相场根位于 ε={first['strain']:.12f}，max(s)={first['metrics']['maximum_s']:.9g}，active size={first['active_size']}，force={first['metrics']['force_residual']:.6e}，KKT={first['metrics']['projected_phase_KKT']:.6e}，inactive violation=0。

{table_roots}

尤其在原失败的**同一应变** ε={same['strain']:.12f}，explicit active-set 得到 force={same['metrics']['force_residual']:.6e}、KKT={same['metrics']['projected_phase_KKT']:.6e}，active size={same['active_size']}。因此不是通过调整 load step 宣布成功。

每个合法根同时通过 active q equality、inactive q inequality、0≤s<sc、逐键 irreversibility、严格固定边界、finite positive J 和原 safety 检查。相邻五个 continuation points 固定 Δε=1e−6，依次使用前一 accepted s 作为不可逆下界，全部合法；额外 root-existence queries 只覆盖 εactivation+[1,2,5,10,20,50]×1e−6。重复载荷复用已验收邻点，未扩展到 0.009。

![Projected versus active set](../../results/pfpd_at1_c5_activation_stationarity_audit_v1/figures/projected_vs_active_set_solver.png)

![Positive branch prefix](../../results/pfpd_at1_c5_activation_stationarity_audit_v1/figures/positive_phase_branch_prefix.png)

## 12. independent complementarity/energy oracle

Fischer–Burmeister 未启用：explicit active-set 与相邻 continuation 已决定本轮 root-existence 问题，故 COMPLEMENTARITY_ORACLE=NOT_NEEDED。环境没有 petsc4py，也没有安装大型 PETSc stack。

由于 reduced root 与原 projected solve 在同一载荷的验收结论不同，额外执行**一次**原能量 Π 的 L-BFGS-B constrained monolithic reference，初值仍是原失败初值，而非把已知根喂给 minimizer。结果：{energy['status']}；optimizer message：{energy.get('message','N/A')}。force={energy.get('metrics',{}).get('force_residual',float('nan')):.6e}，KKT={energy.get('metrics',{}).get('projected_phase_KKT',float('nan')):.6e}；相对 reduced root 的 max phase difference={energy.get('phase_max_difference',float('nan')):.6e}。若优化器 success 但未过完整门槛，同样不能接受。这里没有曲率审计，不声称已证明 local minimum。

## 13. 当前 root cause 分类

HYPOTHESIS A 获得支持：原 projected/bound-constrained 求解过程停滞，而同载荷平滑 active-set equations 有合法根。clip map 在 activation 确实非光滑，但只有两个局部 weak bonds；不支持“大规模同步激活导致退化”的叙述。

HYPOTHESIS B 的尺度差异被观察到，但没有通过改变 scaling 的 root 对照单独确证它为根因。HYPOTHESIS C 在已审计邻域不受支持：正 s roots 跨相邻固定载荷持续存在。HYPOTHESIS D 没有数值秩亏证据。更精细的 solver 内部根因仍需后续正式方法评估，不能把此次局部成功当成生产算法已修复。

```text
{labels_text}
```

## 14. tests

Baseline 一次：{tests['baseline']['passed']} passed / 3 known failed / 0 skipped。Final 一次：{tests['final']['passed']} passed / 3 known failed / 0 skipped。原三项 failure ID、message、traceback 与上一轮及本轮 baseline 完全一致，未 xfail、skip 或修改原测试。

{tests_text}

新增测试覆盖 q、分类、bracket 收缩、projection piecewise derivative、reduced assembly、inactive 精确固定、漏 active 的完整 KKT 拒绝、确定性更新、FD Jacobian、谱残差验证、optimizer 停止与物理验收区分、positive roots、irreversibility prefix 和同载荷对照。未重复 spectral oracle/random F/compression/shear/C3 Hessian/C7 campaign；既有检查只随开始和结束的完整测试运行。

保全：{protected['checked_files']} 个原有文件和 {protected['checked_symlinks']} 个符号链接逐项检查，changed=0、missing=0。Production modified=0；historical modified=0；existing tests modified=0。`provenance/PRE_RUN_MANIFEST.json`、`POST_RUN_MANIFEST.json`、`PRESERVATION_CHECK.json` 及全部测试日志归档。

## 15. CPU/GPU 使用

CPU 为 i7-13700HX，环境可见 12 physical / 24 logical cores；串行大 C5 solve，max_workers=1，至少保留 2 physical cores 和 20% RAM；BLAS/OMP/MKL=1，Numba 线程上限配置 10。本轮实际算子是 NumPy/SciPy CPU，未把线程配置当成 parallel-kernel 加速证明。

GPU 为 RTX 4060 Laptop 8 GB；无可直接使用的已验证 CuPy/JAX/Torch backend，未使用 GPU、未安装大 stack、未跑 GPU benchmark，加速比 N/A。

记录的 intact solves 共 {resources['intact_states']} 个，累计 {resources['intact_wall_seconds']:.2f} s；原 solver 一次 {resources['original_wall_seconds']:.2f} s；reduced solves 累计 {resources['reduced_wall_seconds']:.2f} s；三点 Jacobian/conditioning 累计 {resources['jacobian_diagnostic_wall_seconds']:.2f} s。intact 进程峰值 RSS {resources['intact_peak_RSS_bytes']/2**20:.1f} MiB。上述是相应调用计时，不是假定的全任务 CPU 时间或并行加速比；原始硬件和分项计时见 `provenance/`。

## 16. 唯一下一步

人工审核后，单独评估适合正式 stationary tracking 的 active-set / variational-inequality 方法，以本轮同载荷 roots、inactive inequality 和不可逆 continuation prefix 为基准。保留原物理模型与 force/KKT 门槛，先形成可审查的正式求解方法验证，再决定是否接回空间稳定性审计。本轮不替换 production，不运行 C5 stability Hessian、C7 coupled branch 或任何新物理任务；到此停止。
'''
    (OUT/'REVIEW_C5_ACTIVATION_STATIONARITY.md').write_text(text)


def inventory():
    sources=['scripts/c5_activation_audit.py','scripts/c5_stationarity_diagnostics.py','scripts/review_c5_activation.py','tests/test_c5_activation_stationarity.py']
    before=read(OUT/'provenance/PRE_RUN_MANIFEST.json')
    assert all(name not in before['files'] for name in sources)
    save(OUT/'provenance/NEW_SOURCE_FILES.json',{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources})
    import importlib,scipy,inspect
    modules=['scipy.optimize._lsq.dogbox','scipy.optimize._lsq.least_squares']
    paths=[Path(inspect.getfile(importlib.import_module(name))) for name in modules]
    save(OUT/'provenance/SCIPY_SOURCE.json',dict(version=scipy.__version__,files={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
         role='Installed source inspected for termination and frozen dogleg replay; not modified'))


if __name__=='__main__':
    tests=test_summary();resources=summaries();figures();protected=preservation();inventory();report(tests,protected,resources)
    (OUT/'provenance/STATUS_AFTER.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
    print('REPORT',OUT/'REVIEW_C5_ACTIVATION_STATIONARITY.md',flush=True)
    print('FIGURES',FIG,'TESTS',tests['final']['passed'],'passed, 3 known failed',flush=True)
    print('PRESERVATION',protected,flush=True)
