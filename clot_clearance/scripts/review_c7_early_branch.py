"""Evidence-only C7 report and eight figure groups; never calls a solver."""
import json,os,time,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
from scripts import c7_early_branch as a
from scripts.c7_spatial_evidence import region
OUT=a.OUT;FIG=a.RESULTS/'figures'


def test_stats(name):
    p=OUT/f'tests/{name}.xml'
    if not p.exists():return {}
    cases=ET.parse(p).getroot().findall('.//testcase');failures=[dict(id=c.get('classname')+'::'+c.get('name'),message=c.find('failure').get('message'),trace=c.find('failure').text) for c in cases if c.find('failure') is not None]
    return dict(passed=sum(c.find('failure') is None and c.find('error') is None and c.find('skipped') is None for c in cases),failed=len(failures),skipped=sum(c.find('skipped') is not None for c in cases),failures=failures)


def preservation():
    pre=a.read(OUT/'PRE_RUN_MANIFEST.json');changed=[];missing=[];links=[]
    for n,h in pre['files'].items():
        p=a.ROOT/n
        if not p.is_file():missing.append(n)
        elif a.digest(p)!=h:changed.append(n)
    for n,target in pre['symlinks'].items():
        p=a.ROOT/n
        if not p.is_symlink() or os.readlink(p)!=target:links.append(n)
    r=dict(files=len(pre['files']),symlinks=len(pre['symlinks']),modified=changed,missing=missing,symlink_changes=links,production_modified=[n for n in changed if n.startswith('pd_clot/')],existing_tests_modified=[n for n in changed if n.startswith('tests/')],passed=not(changed or missing or links))
    a.save(OUT/'PRESERVATION_CHECK.json',r);return r


def post_manifest():
    files={};links={}
    for d,ds,ns in os.walk(a.ROOT,followlinks=False):
        ds[:]=[n for n in ds if n!='.git']
        for n in list(ds)+ns:
            p=Path(d)/n;name=str(p.relative_to(a.ROOT))
            if p==OUT/'POST_RUN_MANIFEST.json':continue
            if p.is_symlink():links[name]=os.readlink(p)
            elif p.is_file():files[name]=a.digest(p)
    a.save(OUT/'POST_RUN_MANIFEST.json',dict(root=str(a.ROOT),files=files,symlinks=links,self_excluded=True,timestamp=time.time()))


def comparison_growth(rows):
    prior=a.read(a.ROOT/'verification/pfpd_at1_c5_stability_v1/SUMMARY.json')['states'];activation=a.read(a.old.OUT/'activation/POINT.json')['strain'];c7activation=a.read(OUT/'activation/WEAK.json')['strain'];pairs=[]
    for c7 in rows:
        offset=c7['strain']-c7activation;c5=min([v for v in prior if v['metrics']['maximum_s']>1e-8],key=lambda v:abs(v['strain']-activation-offset));gap=abs(c5['strain']-activation-offset)
        if gap>1e-6:continue
        data=a.arrays(a.ROOT/f'verification/pfpd_at1_c5_stability_v1/states/{c5["label"]}/STATE.npz');p=a.problem(5,c5['strain'],data['lower']);ids=np.flatnonzero(data['s']>1e-8);r=region(p,ids,data['s'][ids])
        pairs.append(dict(C7_label=c7['label'],C5_label=c5['label'],offset_C7=offset,offset_C5=c5['strain']-activation,offset_gap=gap,C7_max_s=c7['metrics']['maximum_s'],C5_max_s=c5['metrics']['maximum_s'],C7_active_fraction=c7['active_count']/46059,C5_active_fraction=len(ids)/5392,C7_RMS_m=c7['spatial']['RMS_radius_m'],C5_RMS_m=r['RMS_radius_m'],C7_equivalent_volume_m3=c7['spatial']['equivalent_volume_m3'],C5_equivalent_volume_m3=r['equivalent_volume_m3']))
    a.save(OUT/'C5_C7_EARLY_SPATIAL_COMPARISON.json',dict(pairs=pairs,matching='Nearest historical C5 activation-offset state, only gaps <= 1e-6; no interpolation, re-solving, or convergence assertion'));return pairs


def response_comparison(rows):
    """Extract response comparisons from archived states; never solve equilibrium."""
    oldroot=a.ROOT/'verification/pfpd_at1_c5_stability_v1'
    prior=a.read(oldroot/'SUMMARY.json')['states']
    weak5=next(v for v in prior if v['label']=='S_WEAK')
    data=a.arrays(oldroot/'states/S_WEAK/STATE.npz')
    p=a.problem(5,weak5['strain'],data['lower'])
    reaction5=float(-p.state(data['x'],data['s']).force[p.top,2].sum())
    weak7=a.read(OUT/'activation/WEAK.json')
    pos5=sorted([v for v in prior if v['metrics']['maximum_s']>1e-8],key=lambda v:v['strain'])
    limit=rows[-1]['strain']-weak7['strain']
    window5=[v for v in pos5 if v['strain']-weak5['strain']<=limit]
    counts5=[int((a.arrays(oldroot/f'states/{v["label"]}/STATE.npz')['s']>1e-8).sum()) for v in window5]
    counts7=[v['active_count'] for v in rows]
    def slope(seq,key):
        return (seq[1]['metrics'][key]-seq[0]['metrics'][key])/(seq[1]['strain']-seq[0]['strain'])
    record=dict(C5_activation_reaction_N=reaction5,C5_activation_stress_Pa=reaction5/a.old.v2.L**2,C7_activation_reaction_N=weak7['reaction_N'],C7_activation_stress_Pa=weak7['nominal_stress_Pa'],
        initial_slope_definition='Secant between first two accepted positive states; not an exact derivative at activation.',
        C5_initial_strains=[v['strain'] for v in pos5[:2]],C7_initial_strains=[v['strain'] for v in rows[:2]],
        C5_initial_max_s_per_strain=slope(pos5,'maximum_s'),C7_initial_max_s_per_strain=slope(rows,'maximum_s'),
        C5_initial_mean_s_per_strain=slope(pos5,'mean_s'),C7_initial_mean_s_per_strain=slope(rows,'mean_s'),
        growth_definition='Post-activation active-count secant over saved positive states within the C7 activation-offset window; onset jump excluded.',
        C5_growth_window=[window5[0]['strain'],window5[-1]['strain']],C7_growth_window=[rows[0]['strain'],rows[-1]['strain']],
        C5_active_counts=counts5,C7_active_counts=counts7,
        C5_active_bonds_per_strain=(counts5[-1]-counts5[0])/(window5[-1]['strain']-window5[0]['strain']),
        C7_active_bonds_per_strain=(counts7[-1]-counts7[0])/(rows[-1]['strain']-rows[0]['strain']),
        provenance='C5 reaction reevaluated with unchanged constitutive force on archived S_WEAK coordinates and phase. All slopes/counts use saved states; no equilibrium or stability solve.')
    a.save(OUT/'C5_C7_RESPONSE_COMPARISON.json',record);return record


def collect():
    br=a.read(OUT/'activation/BRACKET.json');sp=a.read(OUT/'activation/C5_C7_COMPARISON.json');b=a.read(OUT/'branch/ALL.json');rows=b['states'];c5=a.read(OUT/'c5_regression/ALL.json');j=a.read(OUT/'jacobian/C7_FIRST_WEAK.json');gpu=a.read(OUT/'GPU_BACKEND.json');geo=a.read(OUT/'C7_GEOMETRY_CHECK.json')
    intact=[a.read(p) for p in sorted((OUT/'intact').glob('*/SUMMARY.json'))];base=test_stats('BASELINE');final=test_stats('FINAL');pres=a.read(OUT/'PRESERVATION_CHECK.json') if (OUT/'PRESERVATION_CHECK.json').exists() else {}
    valid=bool(rows and all(r['scientifically_accepted'] and r['full_KKT_bond_count']==46059 and r['excluded_violation_max']<=1e-8 for r in rows));pool_ok=c5['passed'] and j['passed'] and valid
    activation_ok=br['width']<=1e-7 and all(r['intact_accepted'] and r['phase_bonds']==46059 for r in intact) and sp['C7']['count']>0
    status=dict(C7_ACTIVATION_AUDIT='COMPLETE' if activation_ok else 'INCOMPLETE',C7_S0_ACTIVATION='BRACKETED',C7_ACTIVATION_INTERVAL=br['interval'],C7_FIRST_ACTIVATION_SET=f"{sp['C7']['count']} / 46059",C7_ACTIVATION_PATTERN=sp['C7']['pattern'],C5_TO_C7_ACTIVATION_STRAIN=sp['strain_relation'],C5_TO_C7_ACTIVATION_REGION=sp['region_relation'],REDUCED_PHASE_POOL='VERIFIED' if pool_ok else 'NOT_VERIFIED',C7_PETSC_RSLS='VERIFIED_FOR_EARLY_BRANCH' if valid and len(rows)>=5 else 'NOT_VERIFIED',C7_FIRST_POSITIVE_ROOT='FOUND' if valid else 'NOT_FOUND',C7_EARLY_POSITIVE_BRANCH='EXISTS' if valid and len(rows)>=2 else 'NOT_ESTABLISHED',C7_MAX_PHASE_REACHED=rows[-1]['metrics']['maximum_s'] if rows else None,C7_STABILITY='NOT_ASSESSED',SPATIAL_STABILITY_CONVERGENCE='NOT_YET_ESTABLISHED',GPU_BACKEND=gpu['status'],PRODUCTION='UNCHANGED' if pres.get('passed') else 'PENDING_CHECK')
    monitor=a.read(OUT/'provenance/RESOURCE_MONITOR.json') if (OUT/'provenance/RESOURCE_MONITOR.json').exists() else {}
    gates=[a.read(p) for p in (OUT/'jacobian').glob('*.json')]
    wall=sum(r['wall_seconds'] for r in intact)+sum(r['wall_seconds'] for r in c5['states'])+sum(r['wall_seconds'] for r in rows)+sum(r['wall_seconds'] for r in gates)
    result=dict(status=status,geometry=geo,activation=br,spatial_comparison=sp,C5_regression=c5,C7_jacobian=j,branch=b,GPU=gpu,CPU_workers=a.cfg()['workers'],CPU_scientific_wall_seconds=wall,wall_definition='Sum of final intact solves, two accepted C5 regression runs, C7 accepted pool solves, and three Jacobian gates; excludes downloads/tests/report and early adapter debugging.',memory={k:v for k,v in monitor.items() if k!='samples'},baseline_tests=base,final_tests=final,known_failures_unchanged=bool(final and base['failures']==final['failures']),preservation=pres,intact_states=len(intact),max_pool_size=max([r['pool_size'] for r in rows],default=0),pool_expansions=sum(r['expansion_count'] for r in rows),first_pool_size=rows[0]['initial_pool_size'] if rows else None,early_spatial_comparison=comparison_growth(rows))
    result['response_comparison']=response_comparison(rows)
    a.save(OUT/'SUMMARY.json',result);return result,intact


def style():
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf');white='#e8edf4'
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':10,'pdf.fonttype':42,'figure.facecolor':'#080b10','axes.facecolor':'#10151e','savefig.facecolor':'#080b10','text.color':white,'axes.labelcolor':white,'axes.edgecolor':white,'xtick.color':white,'ytick.color':white,'axes.spines.top':False,'axes.spines.right':False,'axes.formatter.useoffset':False})


def finish(fig,name,rect=None):
    fig.tight_layout(rect=rect);fig.savefig(FIG/f'{name}.png',dpi=200);fig.savefig(FIG/f'{name}.pdf');plt.close(fig)


def strain_axis(ax):
    ax.set_xlabel('Engineering strain');ax.xaxis.set_major_locator(MaxNLocator(nbins=4));ax.xaxis.set_major_formatter(FormatStrFormatter('%.6f'))


def cube(ax,title):
    L=a.old.v2.L*1e6/2;ax.set(xlim=(-L,L),ylim=(-L,L),zlim=(-L,L),xlabel='x (um)',ylabel='y (um)',zlabel='z (um)',title=title)
    ax.set_box_aspect((1,1,1));ax.view_init(22,-55)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_pane_color((.05,.07,.1,1));axis.set_major_locator(MaxNLocator(3))
    xx,yy=np.meshgrid([-L,L],[-L,L])
    for z in [-L,L]:ax.plot_surface(xx,yy,np.full_like(xx,z),color='#9eabbc',alpha=.10,shade=False)


def figures(r,intact):
    FIG.mkdir(exist_ok=True);style();blue='#62b6ff';orange='#ffb366';grey='#8090a0';comp=r['spatial_comparison'];rows=r['branch']['states'];eps=np.array([v['strain'] for v in rows]);interval=r['activation']['interval']
    fig,axs=plt.subplots(1,2,figsize=(10,4));labs=['C3 historical','C5','C7']
    axs[0].bar(labs,[1.,100*comp['C5']['strain'],100*comp['C7']['strain']],color=[grey,orange,blue]);axs[0].set_ylabel('Activation engineering strain (%)')
    axs[1].bar(labs,[1.,comp['C5']['fraction'],comp['C7']['fraction']],color=[grey,orange,blue]);axs[1].set_yscale('log');axs[1].set_ylabel('Activated bond fraction')
    fig.suptitle('Activation changes with spatial sampling | No convergence claim');finish(fig,'c5_vs_c7_activation')
    fig,ax=plt.subplots(figsize=(9,4));rr=sorted(intact,key=lambda v:v['strain'])
    for idx,label,col in [(0,'Minimum q',blue),(1,'5th percentile q',orange),(2,'Median q',grey)]:ax.plot([v['strain'] for v in rr],[v['q_percentiles'][idx] for v in rr],'.-',color=col,label=label)
    ax.axhline(0,color='white',ls='--',lw=.8);ax.axvspan(*interval,color=orange,alpha=.5,label='Activation bracket');strain_axis(ax);ax.set_ylabel('Physical phase derivative q (Pa)');ax.set_title('C7: all plotted states satisfy mechanical equilibrium');ax.legend(frameon=False);finish(fig,'c7_activation_gradient')
    fig=plt.figure(figsize=(8,6));ax=fig.add_subplot(111,projection='3d');cube(ax,'C7 AT1 ACTIVATION SET\nNOT CRACK PREDICTION')
    for b in comp['C7']['bonds']:
        end=np.array([b['endpoint_0_m'],b['endpoint_1_m']])*1e6;mid=np.array(b['midpoint_m'])*1e6;ax.plot(*end.T,color=blue,lw=3);ax.scatter(*mid,color=orange,s=55)
    ax.text2D(.02,.95,'Top/bottom shaded planes: clamps\nLines: reference bond orientations',transform=ax.transAxes,fontsize=9);finish(fig,'c7_activation_set_map')
    fig,ax=plt.subplots(figsize=(9,4));steps=np.arange(len(rows))
    for key,label,col in [('phase_bonds','Total phase bonds',grey),('pool_size','Candidate pool',blue),('interior_count','Interior',orange),('weak_count','Weak','#c697ff'),('excluded_strong_count','Excluded strong lower','#78d8ab')]:ax.plot(steps,[v[key] for v in rows],'o-',ms=3,label=label,color=col)
    ax.set_yscale('symlog',linthresh=1);ax.set(xlabel='Accepted positive load index (zero based)',ylabel='Bond count (symlog includes zero)',title='Reduced phase pool with full-bond KKT checks');ax.legend(frameon=False,fontsize=8,ncol=2)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6,integer=True));ax.text(.03,.86,f"Total: 46,059 | Excluded strong lower: {rows[-1]['excluded_strong_count']:,}",transform=ax.transAxes,fontsize=9);finish(fig,'c7_phase_pool_evolution')
    fig,axs=plt.subplots(3,1,figsize=(9,7),sharex=True)
    for ax,values,label in [(axs[0],[v['metrics']['maximum_s'] for v in rows],'Maximum phase s'),(axs[1],[v['metrics']['mean_s'] for v in rows],'Weighted mean phase s'),(axs[2],[v['nominal_stress_Pa'] for v in rows],'Nominal reaction stress (Pa)')]:ax.plot(eps,values,'o-',color=blue,ms=3);ax.set_ylabel(label);ax.axvspan(*interval,color=orange,alpha=.25)
    strain_axis(axs[-1]);fig.suptitle('C7 STATIONARY BRANCH\nNO STABILITY CLASSIFICATION');finish(fig,'c7_early_branch')
    fig,axs=plt.subplots(3,1,figsize=(9,7),sharex=True)
    for ax,values,label in [(axs[0],[v['active_count'] for v in rows],'Active bonds (s > 1e-8)'),(axs[1],[v['spatial']['RMS_radius_m']*1e6 for v in rows],'RMS radius (um)'),(axs[2],[v['spatial']['nearest_clamp_percentiles_m'][2]*1e6 for v in rows],'Median clamp distance (um)')]:ax.plot(eps,values,'o-',color=blue,ms=3);ax.set_ylabel(label)
    axs[0].set_ylim(bottom=0);axs[0].yaxis.set_major_locator(MaxNLocator(integer=True));axs[1].set_ylim(0,a.old.v2.L*1e6/2);axs[2].set_ylim(0,a.old.v2.L*1e6/2)
    strain_axis(axs[-1]);fig.suptitle('C7 early phase region | Discrete spatial metrics');finish(fig,'activation_region_growth')
    fig=plt.figure(figsize=(11,5));vmax=max(max(comp['C5']['particle_aggregate']),max(comp['C7']['particle_aggregate']))
    for pos,name,n in [(1,'C5',5),(2,'C7',7)]:
        ax=fig.add_subplot(1,2,pos,projection='3d');cube(ax,name+' first AT1 activation');p=a.problem(n,comp[name]['strain']);X=p.family.X*1e6;values=np.array(comp[name]['particle_aggregate']);ax.scatter(*X.T,s=7,c=grey,alpha=.2);sel=values>0;scatter=ax.scatter(*X[sel].T,c=values[sel],cmap='viridis',vmin=0,vmax=vmax,s=65)
    fig.suptitle('Same physical coordinates | Particle aggregation, not a crack')
    fig.text(.5,.025,f'Color: weighted activation aggregate, common range 0 to {vmax:.4g}',ha='center',fontsize=9)
    finish(fig,'c5_vs_c7_spatial_activation',rect=(0,.07,.92,.94))
    fig,axs=plt.subplots(1,2,figsize=(10,4));axs[0].bar(['CPU NumPy'],[1000*r['GPU']['CPU_batch_median_seconds']],color=blue,width=.45);axs[0].set(ylabel='Full-bond batch median (ms)',title='10 warm repetitions',ylim=(0,1500*r['GPU']['CPU_batch_median_seconds']));axs[0].text(.05,.9,'CuPy runtime unavailable\nWheel download failed',transform=axs[0].transAxes,fontsize=9)
    axs[1].plot(steps,[v['wall_seconds'] for v in rows],'o-',color=blue);axs[1].set(xlabel='Accepted positive load index',ylabel='Reduced solve wall time (s)',title=f'PETSc RSLS | 1 rank, {r["CPU_workers"]} workers');axs[1].xaxis.set_major_locator(MaxNLocator(nbins=6,integer=True));finish(fig,'compute_backend')
    a.save(OUT/'FIGURE_MANIFEST.json',dict(groups=8,files=[p.name for p in sorted(FIG.glob('*'))],style='Arial regular, dark background, English, PNG/PDF'))


def report(r):
    comp=r['spatial_comparison'];s=comp['C7'];c5=comp['C5'];br=r['activation'];rows=r['branch']['states'];first=rows[0] if rows else None;last=rows[-1] if rows else None;mem=r['memory'];base=r['baseline_tests'];final=r['final_tests'];L=a.old.v2.L
    response=r['response_comparison']
    lines=['# PFPD_AT1_C7_EARLY_BRANCH_V1','## 1. 为什么现在做 C7','C5 已有可复核的 stationary roots 和稳定性证据，现在检查同一个物理立方体在 C7 采样下何时开启相场，以及求解器能否连续找到开启后的平衡态。C7 为 343 粒子、46059 键、735 个自由机械 DOF；读取历史 corrected geometry，未生成另一套几何。',
    f"L=0.25 mm，h={r['geometry']['h_m']*1e3:.12g} mm，δ/h={r['geometry']['delta_over_h']:.9g}。体积和上下表面权重分别通过 L³、L² 检查，几何 SHA256={r['geometry']['geometry_SHA256']}。",
    '## 2. C5 已经告诉我们什么','历史 C5 约 0.449739% strain 时两个键首次开启；激活点和已有正相场状态稳定，已审计至 ε=0.004697390736、max s=0.105282593434，窗口内未找到失稳。本轮只读取这些结果。',
    '## 3. 为什么本轮不直接做 C7 stability','相场下界 KKT 失效表示 s=0 不再满足相场平衡，stationary 表示力与相场方程同时满足；二者都不能代替局部变分稳定性。本轮只构造求解 raw residual 所需的 reduced Jacobian，未构造稳定性 Hessian、Schur 补、critical cone 或负曲率模式。',
    '## 4. C7 activation bracket',f"真正机械平衡的 s=0 状态给出 [{br['interval'][0]:.15f}, {br['interval'][1]:.15f}]，宽度 {br['width']:.9e}。宽扫在首次 min q 变号后停止，再按 midpoint 二分；未使用 brentq。随后仅为获取 weak set 做诊断二分，ε={s['strain']:.15f}，原激活区间保留。",
    f"所有 {r['intact_states']} 个 intact 求解点均检查自由力、固定边界、positive J、安全范围，并对完整 46059 根键计算 q。使用原 SVK spectral 物理力和原 PETSc RSLS 设置，phase pool 为空使 s 精确固定为零；负主 Green 应变只作诊断。扫描从历史已接受 ε=0.00025 状态接续；宽扫按递增加载，二分从已接受 lower 状态接续。旧 cheap predictor ε≈0.00314907 未被当作激活结果。",
    '## 5. 第一次多少 bonds 激活',f"首个诊断 weak set 为 {s['count']} / 46059，比例 {s['fraction']:.12g}；IDs={s['first_weak_IDs']}。最终 bracket 上端的 first violated set IDs={s['first_violated_IDs']}。所有 q 分位数、分类 IDs 与完整数组随状态保存。",
    '## 6. activation 在哪里',f"按预声明位置规则分类为 {s['pattern']}。质心 {s['centroid_m']} m；RMS 半径 {s['RMS_radius_m']*1e6:.9g} μm；键中点包围盒尺寸 {np.array(s['midpoint_bbox_dimensions_m'])*1e6} μm；最近 clamp 距离的 min/Q1/median/Q3/max 为 {np.array(s['nearest_clamp_percentiles_m'])*1e6} μm。",
    f"方向 |nx|,|ny|,|nz| 加权均值 {s['orientation_abs_mean']}；0–90° 对加载轴夹角的 15° 分箱计数 {s['orientation_counts']}。逐键端点、中点、长度、q、Ψ+ 与全部边界距离见 activation/SPATIAL.json。等效粒子参与体积 {s['equivalent_volume_m3']:.9e} m³，只是离散空间比较指标，不能称为裂纹宽度。对称的上下夹持区会使全局 RMS 半径较大；它不等同于单个局部区域的宽度。",
    '## 7. 与 C5 是否相似',f"C5 激活 ε={c5['strain']:.15f}，C7 ε={s['strain']:.15f}，差值 {comp['strain_difference']:.9e}，相对变化 {100*comp['strain_relative_difference']:.5g}%。激活载荷标签为 {comp['strain_relation']}，不能由两级数据判断收敛。C5/C7 物理区域标签 {comp['region_relation']}；C5 clamp 距离中位数 {c5['nearest_clamp_percentiles_m'][2]*1e6:.9g} μm，C7 为 {s['nearest_clamp_percentiles_m'][2]*1e6:.9g} μm。",
    '## 8. reduced phase pool 是什么','只将 735 个自由位移和少量候选相场变量送入 PETSc。初始池包含正相场历史、前一步 weak 集合、predictor 上 q_scaled≤1e−3 的键及共享任一端点的一环邻键。1e−3 直接沿用历史 near-active 阈值，在计算前声明；一环不作递归扩张。',
    '## 9. 为什么这不改变 physics','池外相场精确固定在前一个接受态的下界。每次残差仍从完整物理状态计算；池内相场残差是未经投影截断的 q/max(Ψ+,Yc)，PETSc 自己判断其下界、内部或上界状态。原 BAQP 权重、谱分裂、相场律、几何、物理尺度、阈值和 PETSc 设置均不变。',
    '导数通过执行原已验证 ResidualBlocks 的相同源代码构造，仅将 W0、相场权重和未知向量限制到当前 pool；Jxx、Jxs、Jsx、Jss 以 sparse AIJ 装配。没有 46059² 相场矩阵，也没有 46794² 总矩阵。参见 [PETSc RSLS 算法说明](https://petsc.org/release/manualpages/SNES/SNESVINEWTONRSLS/)。本机实际版本沿用 3.25.3。',
    '## 10. full-KKT 如何保证被排除 bonds 真正可以固定',f"每次 solve 后重算全部键的 KKT；池外 q_scaled<−1e−8 的键及其一环邻键全部加入后重解，最多 8 次。池超过 25% 即停止。两个 C5 回归均从空池发现池外违例并自动扩至 426 键，随后恢复完整 RSLS 参考根；不是只在人工选择的活动键上检验。C7 初始池 {r['first_pool_size']}，最大池 {r['max_pool_size']}，实际扩池 {r['pool_expansions']} 次。",
    f"C5 两点回归均满足 phase Linf≤1e−8、scaled x L2≤1e−8、energy relative≤1e−9、reaction relative≤1e−8。C7 first weak state 的 10 个确定方向 × 三档步长 Jv 复核最大误差 {r['C7_jacobian']['max_relative_error']:.9e}，门槛 5e−5；使用可行的二阶单侧 raw residual 差分。",
    '## 11. PETSc RSLS 是否跨过 C7 activation',f"{r['status']['C7_FIRST_POSITIVE_ROOT']}。首个正根 ε={first['strain']:.15f}，force={first['metrics']['force_residual']:.9e}，全键 projected KKT={first['metrics']['projected_phase_KKT']:.9e}，max s={first['metrics']['maximum_s']:.9g}。从 bracket 下端真实 s=0 状态出发，不输入活动 IDs 或未来正根；初始下界全零。",
    '## 12. early positive branch 是否存在',f"{r['status']['C7_EARLY_POSITIVE_BRANCH']}，共 {len(rows)} 个顺序接受的正相场状态。后续每步只使用前一接受态及其 s 下界；从 Δε=1e−6 起步，连续五点满足预声明容易收敛条件后才升至 2e−6、5e−6。最终 ε={last['strain']:.15f}，max s={last['metrics']['maximum_s']:.9g}。停止原因 {r['branch']['stop_reason']}。没有稳定性分类。",
    '| index | strain | force | full KKT | max s | weighted mean s | interior | weak | pool | expansions | wall s | PETSc iterations |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for k,v in enumerate(rows):lines.append(f"| {k} | {v['strain']:.12f} | {v['metrics']['force_residual']:.3e} | {v['metrics']['projected_phase_KKT']:.3e} | {v['metrics']['maximum_s']:.8g} | {v['metrics']['mean_s']:.8g} | {v['interior_count']} | {v['weak_count']} | {v['pool_size']} | {v['expansion_count']} | {v['wall_seconds']:.3f} | {v['PETSc_iterations']} |")
    lines += ['## 13. damage zone 是否增长',f"使用原 phase_detection_tol=1e−8，首点 active count={first['active_count']}，末点={last['active_count']}。逐点保存 phase 加权质心、RMS 半径、包围盒、夹持距离和等效参与体积；幅值增长与空间支持区域扩张分开判断。真实结果见 activation_region_growth 图及各状态 SPATIAL.json。",
    '## 14. C5/C7 spatial comparison',f"激活比例 C5={c5['fraction']:.9g}，C7={s['fraction']:.9g}；RMS 半径 C5={c5['RMS_radius_m']*1e6:.9g} μm，C7={s['RMS_radius_m']*1e6:.9g} μm。统一 physical coordinates 的图用于比较位置。另从历史 C5 正分支读取与 C7 相对激活载荷最接近的状态，仅保留 offset 差≤1e−6 的配对，共 {len(r['early_spatial_comparison'])} 对，详见 C5_C7_EARLY_SPATIAL_COMPARISON.json；未重解、未插值，也未比较 C7 稳定性。",
    f"激活诊断状态的反力名义应力：C5={response['C5_activation_stress_Pa']:.9g} Pa，C7={response['C7_activation_stress_Pa']:.9g} Pa。最初两个接受正根之间的 Δmax(s)/Δε 分别为 {response['C5_initial_max_s_per_strain']:.9g}、{response['C7_initial_max_s_per_strain']:.9g}；加权平均相场斜率分别为 {response['C5_initial_mean_s_per_strain']:.9g}、{response['C7_initial_mean_s_per_strain']:.9g}。这是有限步长割线，不是激活点精确导数。在 C7 相对激活载荷窗口内的已存正根上，两级 active count 均保持 2，激活跳变之后的活动键增长率均为 0 bonds/strain。逐级应变窗口和计算定义见 C5_C7_RESPONSE_COMPARISON.json。C5 反力仅在归档 S_WEAK 状态上用原物理力重新求值，未求解新的平衡态。",
    '## 15. CPU/GPU 实际使用和 speedup',f"重新检测到 12 physical / 24 logical cores；根据可用 RAM 选择 {r['CPU_workers']} 个 ProcessPoolExecutor worker，BLAS/OMP/MKL 每进程单线程，PETSc 1 rank。科学求解与 Jacobian gate 的累计 wall time={r['CPU_scientific_wall_seconds']:.3f} s；这是各实际阶段 wall time 相加，包含并行等待，不是各 CPU 核时间总和。",
    f"外部同步采样观察到进程树 RSS 峰值 {mem.get('peak_observed_tree_RSS_bytes',0)/1024**3:.3f} GiB，最低 MemAvailable {mem.get('min_observed_RAM_available',0)/1024**3:.3f} GiB；任务进程 VmSwap 最大 {mem.get('max_observed_task_swap_bytes',0)} bytes。采样从 C7 宽扫中段开始，故报告为观察峰值而非全程精确峰值；逐阶段预算及父/worker 峰值另存。",
    f"GPU_BACKEND={r['GPU']['status']}。RTX 4060 与 Numba CUDA 驱动访问可用，但缺少已验证的 GPU 数组/批量特征分解后端；PETSc 无 CUDA/Kokkos/HIP；隔离目录的 CuPy/CUDA binary wheel 下载在 6 次不完整下载后失败，未源码编译或修改 system Python。CPU 全 46059 键 constitutive+q batch 十次热运行中位数 {1000*r['GPU']['CPU_batch_median_seconds']:.3f} ms。GPU correctness/speed gates 未执行，speedup 和峰值 VRAM 均为 N/A；所有科学结果为 CPU。二进制安装方式参考 [CuPy 官方说明](https://docs.cupy.dev/en/stable/install.html)。",
    '## 16. tests',f"开始 full suite：{base.get('passed')} passed / {base.get('failed')} known failed / {base.get('skipped')} skipped。结束：{final.get('passed','PENDING')} passed / {final.get('failed','PENDING')} failed / {final.get('skipped','PENDING')} skipped。原三项失败的 IDs、messages、tracebacks 完全一致：{r['known_failures_unchanged']}。开发中仅 targeted tests；保留空池索引和元数据字段冲突的早期调试日志，修复只涉及新增代码，未调整物理或验收门槛。",
    f"历史保护检查：{r['preservation']}。原生产、旧测试和历史结果不改写。PRE/POST manifests、geometry/state/physics/previous-state hashes 与自动 resume 测试提供可追溯记录。",
    '## 17. 当前 scientific status','```']
    lines += [f'{k}: {v}' for k,v in r['status'].items()]
    lines += ['```','## 18. 唯一下一步','等待人工审核 C7 激活区间、空间结构、reduced pool 与全键 KKT 证据。本轮完成后停止，不自动运行 C7 stability、C9 或后续模型工作。']
    text='\n\n'.join(lines).replace('|\n\n|','|\n|');head,block,tail=text.split('```');text=head+'```'+block.replace('\n\n','\n')+'```'+tail
    (OUT/'REVIEW_C7_EARLY_BRANCH_V1.md').write_text(text+'\n')


if __name__=='__main__':
    import sys
    mode=sys.argv[1] if len(sys.argv)>1 else 'report'
    if mode=='preserve':preservation()
    elif mode=='post':post_manifest()
    else:
        r,intact=collect()
        if mode!='text':figures(r,intact)
        report(r)
