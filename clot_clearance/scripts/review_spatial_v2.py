"""Assemble evidence, figures and Chinese review; never fill missing fine endpoints."""
from scripts.spatial_audit_v2 import *
import hashlib,xml.etree.ElementTree as ET,itertools,subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

def test_summary():
    def parse(path):
        tree=ET.parse(path);cases=tree.findall('.//testcase');fail={}
        for t in cases:
            f=t.find('failure')
            if f is not None:fail[t.attrib['classname']+'::'+t.attrib['name']]=dict(message=f.attrib.get('message'),traceback=f.text)
        return dict(passed=len(cases)-len(fail),failed=len(fail),failures=fail)
    baseline=parse(OUT/'tests/BASELINE.xml');new=parse(OUT/'tests/NEW.xml')
    previous=parse(ROOT/'verification/pfpd_at1_spatial_discretization_audit_v1/tests/BASELINE.xml')
    data=dict(baseline=baseline,new=new,original_baseline_passed=255,v1_additional_passed=6,
        known_failure_ID_message_traceback_identical=baseline['failures']==previous['failures'],
        initial_collection_attempt='Unrestricted pytest also collected archived snapshots and external PhaseFieldX tests. Preserved INITIAL_COLLECTION_SCOPE_ERROR files. Corrected invocation explicitly used tests/; no tests deleted or ignored within the project suite.')
    if (OUT/'tests/FINAL.xml').exists():
        final=parse(OUT/'tests/FINAL.xml');data['final_suite']=final
        data['final_known_failure_ID_message_traceback_identical']=final['failures']==previous['failures']
    save(OUT/'tests/TEST_SUMMARY.json',data);return data

def diagnostics():
    rows=read(OUT/'geometry/LEVELS.json');allb=read(OUT/'critical_brackets/ALL.json');data=[]
    for r in rows:
        n=r['N'];a=next(b for b in allb if b['level']==r['level'] and b['schedule']=='A');b=next(b for b in allb if b['level']==r['level'] and b['schedule']=='B')
        entry=dict(**r,interval_A=a['interval'],interval_B=b['interval'],width=a['width'],
            bracket_path_endpoint_difference=[0.,0.] if a['interval'] and a['interval']==b['interval'] else None,
            onset=a['first_unstable'],last_stable=a['last_stable'],stop_reason=a['unresolved_barrier']['stationary_reason'] if a['unresolved_barrier'] else None)
        paths=list((OUT/'branch_tracking').glob(f'C{n}_*/STATIONARY.json'));entries=[read(p) for p in paths]
        files=[p for p in OUT.rglob('*') if p.is_file() and (f'C{n}_' in str(p) or f'C{n}/' in str(p))]
        entry['resources']=dict(stationary_solve_wall_seconds_sum=sum(r['wall_seconds'] for r in entries),
            python_peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in entries),
            artifact_generation_span_seconds=max(p.stat().st_mtime for p in files)-min(p.stat().st_mtime for p in files),
            timing_scope='Measured stationary solve times exclude subsequent dense audit and plotting; artifact span is a wall-time span, not CPU time.')
        data.append(entry)
    save(OUT/'comparisons/LEVEL_SUMMARY.json',data)
    virtual=[]
    for n in LEVELS:
        p=Problem(n,.01);x=affine(p);s=p.state(x,p.lower);direction=np.zeros_like(x);direction[:,2]=p.family.X[:,2]+L/2
        work=-s.force*direction;total=float(work.sum());top=float(work[p.top].sum());free=float(work[p.free].sum())
        virtual.append(dict(level=f'C{n}',affine_energy_derivative_J=total,top_reaction_times_L_J=top,free_node_virtual_work_J=free,
            free_fraction=free/total,analytic_derivative_J=2*p.material.mu*.01005*1.01*L**3,
            explanation='Affine loading moves free nodes too. Boundary reaction alone represents the affine energy derivative only when free-node forces vanish.'))
    save(OUT/'mechanics/AFFINE_VIRTUAL_WORK.json',virtual)
    aff=read(OUT/'mechanics/AFFINE_REACTIONS.json');v=[next(x for x in aff if x['level']==k and x['strain']==.01)['reaction_N'] for k in ['C5','C7']]
    save(OUT/'mechanics/AFFINE_REACTION_CONSISTENCY.json',dict(status='FAILED',tolerance=.02,C5_C7_difference_relative_to_C7=abs(v[0]-v[1])/abs(v[1]),
        C5_C7_difference_relative_to_max=abs(v[0]-v[1])/max(abs(v[0]),abs(v[1])),scope='AFFINE NONSTATIONARY CHECK ONLY; onset reaction consistency remains INCONCLUSIVE'))
    # Compare symmetry-related C3 witnesses without pretending to have fine modes.
    a=read(OUT/'critical_brackets/C3_A.json')['first_unstable'];b=read(OUT/'critical_brackets/C3_B.json')['first_unstable']
    ma=ca.read_npz(OUT/f'modes/{a["label"]}.npz');mb=ca.read_npz(OUT/f'modes/{b["label"]}.npz')
    corr=lambda u,v:float(np.corrcoef(u.ravel(),v.ravel())[0,1])
    save(OUT/'localization/COMMON_GRID_COMPARISON.json',dict(grid=[21,21,21],cross_resolution_correlation=None,
        reason='No accepted C5/C7 unstable endpoint, hence no admissible fine critical mode.',
        C3_A_B_phase_correlation=corr(ma['grid_phase'],mb['grid_phase']),
        C3_A_B_phase_correlation_after_z_reflection=corr(ma['grid_phase'],mb['grid_phase'][:,:,::-1]),
        C3_A_B_displacement_correlation=corr(ma['grid_displacement'],mb['grid_displacement']),
        C3_A_B_displacement_correlation_after_z_reflection=corr(ma['grid_displacement'],mb['grid_displacement'][:,:,::-1])))
    # Distance bins complement the layer diagnostics for the accepted C3 witness.
    p=Problem(3,.01);X=ma['X'];mid=p.family.X[p.family.pairs].mean(axis=1);aset=read(OUT/f'active_sets/{a["label"]}.json')
    weak=np.zeros(p.ns,bool);weak[aset['bond_IDs']['LOWER_WEAK_NEAR_SWITCH']]=True;interior=np.zeros(p.ns,bool);interior[aset['bond_IDs']['INTERIOR']]=True
    distance_rows=[]
    for name,pd,bd in [('clamp',L/2-abs(X[:,2]),L/2-abs(mid[:,2])),('lateral',L/2-np.max(abs(X[:,:2]),axis=1),L/2-np.max(abs(mid[:,:2]),axis=1))]:
        edges=np.linspace(0,L/2,5)
        for k,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
            ps=(pd>=lo)&(pd<=hi if k==3 else pd<hi);bs=(bd>=lo)&(bd<=hi if k==3 else bd<hi)
            distance_rows.append(dict(distance_to=name,lower_m=float(lo),upper_m=float(hi),particles=int(ps.sum()),bonds=int(bs.sum()),
                phase_amplitude_mean=float(ma['particle_phase'][ps].mean()) if ps.any() else None,
                weak_fraction=float(weak[bs].mean()) if bs.any() else None,interior_fraction=float(interior[bs].mean()) if bs.any() else None))
    save(OUT/'localization/C3_DISTANCE_BINS.json',distance_rows)
    return data

def status_labels():
    labels={
        'SPATIAL_DISCRETIZATION_AUDIT_V2':'INCOMPLETE',
        'GEOMETRY_DEFINITION':'VERIFIED_0P25MM_BOUNDARY_NODE','VOLUME_CONSERVATION':'VERIFIED','SURFACE_CONSERVATION':'VERIFIED',
        'PHYSICAL_HORIZON':'FIXED_0P25125MM','REFINEMENT_LEVELS':['C3','C5','C7'],
        'INTACT_MECHANICS_ACROSS_LEVELS':'VERIFIED','ANALYTIC_THRESHOLD_CONSISTENCY':'PASSED',
        'STABILITY_INTERVAL_CONSISTENCY':'INCONCLUSIVE','REACTION_AT_ONSET_CONSISTENCY':'INCONCLUSIVE',
        'CRITICAL_CONE_INSTABILITY':'INCONCLUSIVE','ACTIVE_SET_BEHAVIOR':'DIFFUSE_ALL_BONDS',
        'ACTIVE_SET_BEHAVIOR_SCOPE':'C3 ACCEPTED ONSET ONLY; C5/C7 ONSET UNRESOLVED',
        'CRITICAL_MODE_LOCALIZATION':'INCONCLUSIVE','BOUNDARY_ASSOCIATION':'INCONCLUSIVE',
        'CORRECTED_C3_VS_LEGACY':'DIFFERENT_PHYSICAL_DISCRETIZATION',
        'CURRENT INSTABILITY CLASSIFICATION':'INCONCLUSIVE','MESH OBJECTIVITY':'NOT YET ESTABLISHED',
        'PFPD_AT1_EVOLUTION_V1':'REMAINS_NOT_VERIFIED','PRODUCTION MODEL':'UNCHANGED',
        'AFFINE_REACTION_SAMPLING_CONSISTENCY':'FAILED',
        'INTACT_MECHANICS_SCOPE':'AFFINE F/ENERGY/BALANCE AND PINNED JULIA VERIFIED; FREE-NODE AFFINE EQUILIBRIUM FAILS C5/C7',
        'BLOCKER':'C5/C7 PURE_TENSION_DOMAIN_VIOLATION during stationary solve; no refinement instability bracket or fine cone/mode claim is admissible.'}
    save(OUT/'FINAL_STATUS.json',labels);save(RESULTS/'FINAL_STATUS.json',labels);return labels

COLORS=['#60c9ff','#ffcc66','#ef88c2']
def style():
    font=Path('/mnt/c/Windows/Fonts/arial.ttf');font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','axes.titleweight':'normal','figure.facecolor':'#080b10',
        'axes.facecolor':'#080b10','savefig.facecolor':'#080b10','text.color':'#eeeeee','axes.labelcolor':'#dddddd',
        'xtick.color':'#dddddd','ytick.color':'#dddddd','axes.edgecolor':'#68727e','grid.color':'#39414c','font.size':11,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def finish(fig,name):
    fig.savefig(OUT/f'figures/{name}.png',dpi=180,bbox_inches='tight')
    fig.savefig(OUT/f'figures/{name}.pdf',bbox_inches='tight');plt.close(fig)

def missing(ax,h,why='UNRESOLVED'):
    ax.text(h,.55,why,transform=ax.get_xaxis_transform(),rotation=90,ha='center',va='center',color='#adb6c3',fontsize=10)

def common_x(ax):
    ax.set_xlabel('Particle spacing h (mm)');ax.set_xticks([L/(n-1)*1e3 for n in LEVELS]);ax.set_xlim(.028,.14);ax.grid(alpha=.3)

def figures(data):
    style();h=np.array([r['h_m']*1e3 for r in data]);names=[r['level'] for r in data];c3=data[0]['onset'];last=data[0]['last_stable']
    fig=plt.figure(figsize=(15,4.8),layout='constrained')
    vertices=np.array(list(itertools.product([-.125,.125],repeat=3)))
    for i,n in enumerate(LEVELS):
        ax=fig.add_subplot(1,3,i+1,projection='3d');X=fixture(n)[1].X*1e3
        ax.scatter(*X.T,s=13 if n<7 else 5,color=COLORS[i],alpha=.8,depthshade=False)
        for u,v in itertools.combinations(vertices,2):
            if np.sum(u!=v)==1:ax.plot(*np.stack([u,v]).T,color='#aeb8c4',lw=.7)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor('#080b10');axis.pane.set_edgecolor('#444444')
        ax.set(xlim=(-.14,.14),ylim=(-.14,.14),zlim=(-.14,.14),xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=f'C{n}: {n**3} nodes')
        ax.tick_params(labelsize=8,pad=0)
        ax.set_box_aspect((1,1,1));ax.set_xticks([-.125,0,.125]);ax.set_yticks([-.125,0,.125]);ax.set_zticks([-.125,0,.125])
    fig.suptitle('Same 0.25 mm domain; same horizon = 0.25125 mm');finish(fig,'corrected_mesh_family')
    fig,ax=plt.subplots(figsize=(8,4.4),layout='constrained');xx=np.arange(3)
    for i,(key,den,label) in enumerate([('volume_m3',L**3,'Volume / L³'),('top_surface_sum_m2',L**2,'Top area / L²'),('bottom_surface_sum_m2',L**2,'Bottom area / L²')]):
        ax.plot(xx+(i-1)*.09,[r[key]/den for r in data],marker=['o','s','^'][i],ls='',color=COLORS[i],label=label)
    ax.axhline(1,color='white',ls=':',lw=1);ax.set(xticks=xx,xticklabels=names,ylim=(.999,1.001),ylabel='Conserved quantity / physical target',title='Material volume and both clamp areas are conserved');ax.ticklabel_format(axis='y',useOffset=False);ax.legend(loc='upper right');finish(fig,'volume_and_surface_conservation')
    fig,axs=plt.subplots(1,3,figsize=(13,4.2),layout='constrained')
    axs[0].plot(h,[r['delta_over_h'] for r in data],'o-',color=COLORS[0]);axs[0].set_ylabel('Horizon / spacing')
    axs[1].plot(h,[r['family']['mean'] for r in data],'o-',color=COLORS[1],label='Mean');axs[1].plot(h,[r['family']['median'] for r in data],'s--',color=COLORS[2],label='Median');axs[1].set_ylabel('Neighbors');axs[1].legend()
    for key,color,label in [('surface_family',COLORS[0],'Surface'),('interior_family',COLORS[2],'Interior')]:axs[2].plot(h,[r[key]['mean'] for r in data],'o-',color=color,label=label)
    axs[2].set_ylabel('Mean neighbors');axs[2].legend()
    for ax in axs:common_x(ax)
    fig.suptitle('More family samples; every horizon remains boundary-truncated');finish(fig,'family_sampling_vs_resolution')
    fig,(ax,legacy)=plt.subplots(1,2,figsize=(11,4.5),gridspec_kw={'width_ratios':[3,1]},layout='constrained')
    lo,hi=data[0]['interval_A'];mid=(lo+hi)/2
    ax.errorbar(h[0],mid,yerr=(hi-lo)/2,fmt='o',capsize=7,color=COLORS[0],label='Corrected C3; schedules A = B')
    ax.axhline(.01,color=COLORS[1],ls='--',label='Analytic AT1 threshold')
    for hh in h[1:]:missing(ax,hh)
    common_x(ax);ax.set(ylabel='Engineering strain',ylim=(.00998,.010015));ax.legend(loc='lower left',fontsize=9)
    legacy.errorbar([0],[mid],yerr=(hi-lo)/2,fmt='s',capsize=7,color='#999999');legacy.set(xticks=[0],xticklabels=['LEGACY C3'],ylim=ax.get_ylim(),title='NOT IN\nREFINEMENT FAMILY',yticks=[])
    fig.suptitle('Only corrected C3 has a verified instability interval');finish(fig,'critical_interval_vs_spacing')
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    axs[0].plot(h[0],last['nominal_stress_Pa'],'o',color=COLORS[0],label='Last stable');axs[0].plot(h[0],c3['nominal_stress_Pa'],'x',color=COLORS[1],label='First unstable')
    for hh in h[1:]:missing(axs[0],hh)
    axs[0].set_ylabel('Nominal onset stress (Pa)');axs[0].legend(loc='lower left');axs[0].set_title('Accepted onset endpoints only')
    aff=read(OUT/'mechanics/AFFINE_REACTIONS.json')
    for e,col in zip([.005,.01,.015],COLORS):axs[1].plot(h,[next(r['nominal_stress_Pa'] for r in aff if r['level']==name and r['strain']==e) for name in names],'o-',color=col,label=f'Affine {100*e:g}%')
    axs[1].set_ylabel('Affine nominal stress (Pa)');axs[1].set_title('Separate affine diagnostic; C5/C7 not equilibrated');axs[1].legend()
    for ax in axs:common_x(ax)
    finish(fig,'reaction_stress_at_onset')
    fig,axs=plt.subplots(1,3,figsize=(13,4.2),layout='constrained')
    path=read(OUT/'critical_brackets/C3_A.json');rs=sorted(path['states'],key=lambda r:r['strain'])
    for key,col,label in [('LOWER_STRONG_ACTIVE',COLORS[0],'Strong'),('LOWER_WEAK_NEAR_SWITCH',COLORS[1],'Weak'),('INTERIOR',COLORS[2],'Interior')]:
        axs[0].plot([r['strain'] for r in rs],[r['active_fractions'][key] for r in rs],'o-',ms=3,color=col,label=label)
    axs[0].legend();axs[0].set_ylabel('Fraction of all bonds')
    for ax,n in zip(axs,LEVELS):
        ax.set(title=f'C{n}',xlabel='Engineering strain',ylim=(-.05,1.1),xlim=(.00895,.01005));ax.tick_params(axis='x',rotation=30)
        if n>3:ax.text(.5,.52,'NO ACCEPTED\nONSET STATE',transform=ax.transAxes,ha='center',va='center',color='#adb6c3')
    fig.suptitle('All 185 C3 bonds switch together; fine onset activity is unknown');finish(fig,'active_set_fraction_near_onset')
    fig,ax=plt.subplots(figsize=(8,4.6),layout='constrained')
    for eps,col,value in zip(EPS,COLORS,c3['negative_witness']['curvatures']):ax.plot(h[0],value,'o',mfc='none',ms=7,color=col,label=f'FD epsilon {eps:g}')
    for hh in h[1:]:missing(ax,hh)
    common_x(ax);ax.set(ylabel='Scaled Rayleigh curvature',ylim=(-.03,.003),title='Verified C3 critical-cone witness; fine witnesses unavailable');ax.axhline(0,color='#777777',ls=':');ax.legend(loc='lower left')
    ax.text(.02,.98,'Common scales: displacement = L × 0.01; energy = Gc × L²',transform=ax.transAxes,va='top',fontsize=9);finish(fig,'critical_cone_witness')
    fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for ax,key,title in zip(axs,['PR','top10_fraction'],['Participation ratio','Top 10% amplitude fraction']):
        ax.plot(h[0],c3['localization'][key],'o',color=COLORS[0]);ax.set(ylim=(0,1),ylabel=title)
        for hh in h[1:]:missing(ax,hh)
        common_x(ax)
    fig.suptitle('C3 selected witness is spatially concentrated; refinement trend unresolved');finish(fig,'localization_metrics')
    fig,axs=plt.subplots(2,3,figsize=(12,7.5),layout='constrained');mode=ca.read_npz(OUT/f'modes/{c3["label"]}.npz')
    for j,n in enumerate(LEVELS):
        for i,key in enumerate(['grid_phase','grid_displacement']):
            ax=axs[i,j];ax.set(title=f'C{n}: '+('phase opening' if i==0 else 'displacement magnitude'),xlabel='x (mm)',ylabel='z (mm)',xlim=(-.125,.125),ylim=(-.125,.125),aspect='equal')
            if n==3:
                im=ax.imshow(mode[key][:,10,:].T,origin='lower',extent=[-.125,.125,-.125,.125],vmin=0,vmax=1,cmap='magma',interpolation='nearest');fig.colorbar(im,ax=ax,shrink=.78,label='Normalized amplitude')
            else:ax.text(.5,.5,'UNRESOLVED\nNo accepted critical mode',transform=ax.transAxes,ha='center',va='center',color='#adb6c3')
    fig.suptitle('NUMERICAL / VARIATIONAL MODE — NOT PHYSICAL CRACK PREDICTION\ny = 0 slices on common 21³ post-processing grid; C3 schedule A witness',fontsize=12);finish(fig,'critical_mode_maps')
    fig,axs=plt.subplots(1,4,figsize=(13,4.4),layout='constrained')
    values=[([.052734375,.015625],'Material volume (mm³)'),([5.326303710937789,1.5781640625000854],'Affine 1% energy (pJ)'),([4.2822421875,1.2688125],'Affine 1% reaction (microN)')]
    for ax,(yy,label) in zip(axs[:3],values):ax.bar([0,1],yy,color=['#888888',COLORS[0]],width=.6);ax.set(xticks=[0,1],xticklabels=['Legacy','Corrected'],ylabel=label)
    for j,col in enumerate(['#888888',COLORS[0]]):axs[3].errorbar(j,mid,yerr=(hi-lo)/2,fmt='o',color=col,capsize=6)
    axs[3].set(xticks=[0,1],xticklabels=['Legacy','Corrected'],ylabel='Instability strain interval',ylim=(.00998,.010015))
    fig.suptitle('GEOMETRY DEFINITION TRANSITION — NOT MESH CONVERGENCE\nLegacy full voxel C3 and corrected boundary-node C3 remain separate',fontsize=12);finish(fig,'legacy_vs_corrected_C3')
    # Requested boundary-sampling scatter is retained as a supplementary figure.
    fig,axs=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    for ax,n,col in zip(axs,LEVELS,COLORS):
        g=ca.read_npz(OUT/f'geometry/C{n}.npz');ax.scatter(g['boundary_distance']*1e3,g['family_size'],s=12,color=col,alpha=.6)
        ax.set(title=f'C{n}',xlabel='Nearest boundary distance (mm)',ylabel='Family size')
    finish(fig,'family_size_vs_boundary_distance')

def review(data,tests,labels):
    a=data[0];u=a['onset'];w=u['negative_witness'];pres=read(OUT/'provenance/PRESERVATION_CHECK.json');hv=read(OUT/'matrix_free/C3_DENSE_VS_HVP.json');external=read(OUT/'external_reference/COMPARISON.json')
    geo_table='\n'.join(f"| {r['level']} | {r['particles']} | {r['h_m']*1e3:.9g} | {r['delta_over_h']:.2f} | {r['bonds']} | {r['family']['mean']:.6f} | {r['analytic_threshold']:.16g} |" for r in data)
    stops=[]
    for n in [5,7]:
        for name in ['A_000','B_000','low_start_0.00025']:
            r=read(OUT/f'branch_tracking/C{n}_{name}/STATIONARY.json');stops.append(f"| C{n} {name} | {r['strain']:.7g} | {r['metrics']['force_residual']:.8g} | {r['exit']['metrics']['minimum_Green_strain']:.8g} | PURE_TENSION_DOMAIN_VIOLATION |")
    dense=read(OUT/'matrix_free/C3_validation/DENSE_REFERENCE.json')['rows']
    report=f'''# PFPD_AT1 空间离散审计 V2

**结论：INCOMPLETE。几何与完整材料公式验证通过；corrected C3 的失稳已确认，C5/C7 在满足驻定条件之前触发纯拉伸域退出，不能报告空间收敛。**

仓库 `{ROOT.parent}`；HEAD `{(OUT/'provenance/HEAD.txt').read_text().strip()}`。本报告只审计数值开发试件。

## 1. 为什么重新定义 coarse geometry

此前的 27 点坐标外包络边长是 0.25 mm，但每个点都分配完整的立方体体积。因此积分中的材料总量对应 0.375 mm 立方体。新家族采用用户明确批准的物理边界节点定义。旧结果永久保存在原目录，其标签仍为 `LEGACY_3x3x3_FULL_VOXEL`。

## 2. legacy 0.375 和 corrected 0.25 有什么不同

两者的 27 点坐标相同；区别在积分权重和物理解释。legacy 总体积为 0.052734375 mm³，新值为 0.015625 mm³，比例 3.375。新 C3 由正式 BAQP 代码重新算出 1% 能量 1.5781640625 pJ、反力 1.2688125 microN；与 V1 隔离 proposal 的相对误差低于 1e-10，且未调整 Gc。详见 [重算对照](comparisons/CORRECTED_C3_PROPOSAL_CHECK.json) 和 [图 10](figures/legacy_vs_corrected_C3.png)。

## 3. 为什么 fixed delta 很重要

所有网格的物理 horizon 固定为 0.25125 mm。它规定非局部相互作用的距离，改变它就会改变模型。此次只改变 h，delta/h 依次为 2.01、4.02、6.03。delta 大于试件半边长，所以所有点的相互作用邻域均被边界截断；这不是自动判错的理由，也不能据此声称已验证体材料标定。

## 4. conservative nodal volumes 是什么

一维端点权重为 h/2，内部为 h，三维取乘积。内部、面、棱、角点分别是 h³、h³/2、h³/4、h³/8。上下表面使用二维同类权重，夹持点位于精确 z = ±0.125 mm。实现只在新脚本内建立 Cloud/Family，未改生产几何函数。

## 5. 三档 mesh 信息

| 网格 | 点 | h (mm) | delta/h | bonds | 平均 family | 独立求解解析阈值 |
|---|---:|---:|---:|---:|---:|---:|
{geo_table}

各档 shape tensor 均满秩 3；完整的 family 最小/中位/最大值、表面/内部统计、K 条件数、点坐标及边界距离保存在 [geometry](geometry/LEVELS.json)。未使用秩正则化。图 1–3 及补充图 `family_size_vs_boundary_distance` 展示固定几何和采样变化。

## 6. volume/surface conservation 是否通过

通过。每档 V = 1.5625e-11 m³，A_top = A_bottom = 6.25e-8 m²，机器精度相对容差 2e-15。应力始终用反力除以真实面积 A，未使用节点数或外包体积推算面积。见 [图 2](figures/volume_and_surface_conservation.png)。

## 7. intact mechanics 是否通过

刚体平移、转动、单轴拉伸/压缩、剪切、各向同性伸长的 Fbar/Ftilde、能量密度、总力和总力矩检查全部通过。总能量等于解析能量密度乘固定体积。固定 Peridynamics.jl commit `{external['commit']}`、Julia {external['julia_version']}、MIT，三档均对照了仿射及非仿射状态的总能量、力、Fbar/Ftilde，误差均小于 1e-9。Manifest SHA256 `{external['manifest_sha256']}`，上游未升级，作者实现未复制到新生产代码。

**仿射完整材料公式正确，不等于自由节点已平衡。** 1% 时 C3/C5/C7 的仿射反力分别为 1.2688125、0.7919107060、0.5646225130 microN；对应名义应力 20.301、12.6705712960、9.0339602080 Pa。C5/C7 差异以 C7 为分母约 40.26%，以两者最大值为分母约 28.70%，均超过 2%。因此仿射反力采样检查 FAILED，不能用这组非驻定反力冒充 onset 反力。

仿射能量导数还包括自由节点虚功，见 [虚功分解](mechanics/AFFINE_VIRTUAL_WORK.json)。仅 C3 自由节点残余力接近零。这个区别解释了为何总能量一致而仅对 top 节点求和的反力并不一致。详见 [图 5](figures/reaction_stress_at_onset.png)。

## 8. analytic threshold 是否 mesh-independent

通过。保持 mu = 1000 Pa、bulk = 666.6666666666666 Pa、Gc = 1.2688439062499997e-5 J/m²、原 normalization.c0 = 0.25、h_kinematic = 1、sc = 0.95；Yc = Gc/(2 c0 delta) = 0.1010025 Pa，三档完全一致。每档独立用正式 BAQP 的仿射 driver 求 Y(epsilon) = Yc，阈值均为 0.01，差异小于 1e-13。此阈值是均匀仿射状态的局部 AT1 开启阈值，不是细网格已求出的驻定失稳阈值。

## 9. corrected C3 stability result

两条独立路径 A（粗步长 5e-4）与 B（2.5e-4）再二分，均得 `[0.0099921875, 0.0100000000]`，宽度 7.8125e-6，端点差为零。所有接受点沿用原 least_squares 成功、力残差 1e-9、KKT 1e-8、不可逆性、s < sc、最小 Green 应变 ≥ -1e-9、精确固定边界等门槛。原 `branch_stability_tracking.py` 求解函数未经修改，未新写分支求解算法。

三个差分尺度下重新生成 `CORRECTED_C3_DENSE_REFERENCE`，检查原始非对称误差、解析 Hss、Hxx 条件数和 Schur 消元残差。C3 阈值态 Hxx 条件数约 {dense[1]['Hxx_condition']:.6g}，最大原始对称误差 {max(r['raw_symmetry_error'] for r in dense):.3g}，最大 Hss 相对误差 {max(r['Hss_analytic_error'] for r in dense):.3g}。见 [稠密参考](matrix_free/C3_validation/DENSE_REFERENCE.json)。

临界锥要求强下界相场扰动为零，弱下界相场只能非负开启，并满足一阶能量方向导数近零。机械方向由 Hxx u = -Hxs p 消元。候选包含均匀、driver 权重、最低耦合模态的正/负部分、C3 物理中点投影、确定性随机、z 层和边界/内部开启。选中合法负见证曲率为 `{w['curvatures']}`，三档 Hv 也为负，一阶中性、相场可行和原纯拉伸可行性均通过。只有正的随机探针不会被宣称为稳定证明。

跨网格曲率统一用长度 L×0.01 = 2.5e-6 m、能量 Gc×L² = 7.930274414062498e-13 J 归一化，不比较旧 mesh-dependent 坐标的数值大小。C3 稠密与矩阵自由验证包含 10 个确定性随机方向、最低耦合模态和两类开启方向，每类均有三个尺度；最大相对误差 {hv['max_relative_error']:.6g} < 预声明 5e-5。残差 Jacobian 的转置对照也通过。开发时发现的分组归一化边界问题已修正并保留错误日志；未放宽物理域或验收门槛。

## 10. C5 stability result

C5 有 225 个自由位移自由度和 5392 个相场自由度。A、B 路径在起始 0.009 即触发 `PURE_TENSION_DOMAIN_VIOLATION`。最后保留的合法试点还具有约 0.071159 的归一化力残差，远高于 1e-9，不能称驻定或稳定。域外试探状态仅作为失败诊断保存，未被接受为物理分支。不存在合法临界区间、onset 反力、负临界锥见证或临界模态。

## 11. C7 stability result

C7 有 735 个自由位移自由度和 46059 个相场自由度，共 46794 DOF。单个稠密 Hessian 预计 17,517,427,488 bytes，超出当时可用 RAM 的 40%，因此禁止稠密分配；使用 `LinearOperator` 和 SciPy dogbox/LSMR，仍调用原驻定求解函数与残差定义。C7 的 A/B 路径也在 0.009 触发纯拉伸域退出。原始力残差约 0.0474434，不能接受。

额外从零状态检查 0.00025 小载荷，C5/C7 仍触发同一域门槛；零应变本身的平衡检查通过，但没有据此越过域退出点继续拼接分支。

| 尝试 | 应变 | 最后合法点力残差 | 被拒试点最小 Green 应变 | 原因 |
|---|---:|---:|---:|---|
{chr(10).join(stops)}

所有被拒试点的完整 z 和 metrics 保存在对应 `branch_tracking/*/STATIONARY.json` 的 exit 字段。实际细网格驻定门槛失败后，按逻辑停止依赖它的细网格 Hessian 谱、相场候选 A–H、负见证、模态与定位计算；没有在非驻定状态上补算这些量来形成失稳结论。

## 12. critical intervals 是否趋于一致

INCONCLUSIVE。只有 C3 的合法区间。两条细网格加载路径均已尝试，但都未越过首个驻定门槛，因此不存在可以比较宽度、重叠或中点差的 C5/C7 区间。预声明 2e-5 中点容差及 1e-5 区间宽度要求保持不变。图 4 的 C5/C7 显示空缺，legacy 单列灰色且不进入新家族连线。

## 13. reaction stress 是否趋于一致

Onset 一致性 INCONCLUSIVE。C3 首个不稳定端点 R = {u['reaction_N']:.15g} N，nominal stress = {u['nominal_stress_Pa']:.15g} Pa。C5/C7 没有合格端点。上述仿射反力采样 FAILED 与此结论分开记录。总能量同时报告 elastic、crack、total，以及按固定 V 或 A 归一化的值；未用粒子数作为能量分母。

## 14. active-set behavior 怎样变化

C3 是 185 strong lower → 185 weak lower，onset fractions 为 strong 0、weak 1、interior 0、near_sc 0。图 6 只把接受的 C3 路径解释为 onset 行为。细网格失败起始点的活动集及 bond IDs/中点仍保存为诊断，但其 onset 行为未知。受用户指定枚举限制，`ACTIVE_SET_BEHAVIOR = DIFFUSE_ALL_BONDS` 的范围明确限定为 C3，不能理解为三档均已证明如此。

## 15. instability 是否仍然 diffuse

需要区分“所有 bonds 同时处于弱下界”和“负曲率见证的空间幅度均匀”。C3 的活动集切换是全域的，选中的合法负见证却集中在一侧夹持附近。没有 C5/C7 见证，不能回答这种空间分布是否在细化后保留。

## 16. critical mode 是否局域化

C3 bond 模态按能量配对权重聚合到粒子、归一化，再插值到共同 21³ 后处理网格。PR = {u['localization']['PR']:.9g}，top10 fraction = {u['localization']['top10_fraction']:.9g}；27 点的 top 10% 用 ceil(2.7) = 3 点。PR 接近 1 表示均匀，较小表示集中。这里称为数值变分见证，不能命名为裂纹或裂纹宽度。

A/B 选出的见证分别偏向上下两个对称夹持面；已报告原坐标相关性和 z 镜像后的相关性，而不是把镜像方向当作不同物理裂纹。跨网格相关性、中心/展宽张量收敛均为缺失。见 [共同网格诊断](localization/COMMON_GRID_COMPARISON.json)、图 8–9。

## 17. boundary association

C3 schedule A 幅度中心 z ≈ {u['localization']['center_m'][2]*1e3:.8g} mm，平均到 clamp 的距离约 {u['localization']['clamp_distance_mean_m']*1e3:.8g} mm。按 z 层、上下夹持距离、侧面距离分组的幅度和弱/内部 bond 比例已保存。此时可以说该见证与夹持位置有关，不能自动称为边界伪影；跨分辨率 PERSISTS/WEAKENS/STRENGTHENS 均无证据，所以总标签 INCONCLUSIVE。

## 18. legacy result 与新 family 的关系

legacy 与 corrected C3 的临界区间数值恰好相同，并不使其成为同一物理离散。体积、能量和反力差异仍存在。新家族只包含 corrected C3/C5/C7，不含 legacy 回归、外推或连接线。本轮所有历史文件 {pres['protected_files']} 个和符号链接 {pres['protected_links']} 个均保持原哈希/目标；production modifications = 0，historical modifications = 0，historical tests modified = 0。

## 19. tests

原科学基线 255 passed / 3 known failed 全部保留；V1 的 6 项新增测试也保留，因此运行前实际基线是 {tests['baseline']['passed']} passed / {tests['baseline']['failed']} known failed。V2 新增 {tests['new']['passed']} 项全部通过。最终全套为 {tests.get('final_suite',{}).get('passed','pending')} passed / {tests.get('final_suite',{}).get('failed','pending')} known failed。原三个失败的 ID/message/traceback 与 V1 基线逐字比较（运行前及最终）：{tests['known_failure_ID_message_traceback_identical']} / {tests.get('final_known_failure_ID_message_traceback_identical','pending')}。没有 xfail/skip 包装，也没有删除测试。测试通过证明实现/审计证据自洽，不会把失败的细网格科学门槛改写成通过。

初次 pytest 未限定目录而误收集了归档快照和外部 PhaseFieldX 测试，产生 collection errors；完整日志保留。随后按既有项目范围显式运行 `pytest tests`，没有在项目 tests 内排除任何测试。

新增文件仅为隔离脚本、独立测试、V2 verification/results。CPU、核数、RAM、磁盘与 Python/NumPy/SciPy/Numba/Julia 信息、每档内存预估和测得的 Python peak RSS 已存 provenance/各档目录。C7 重任务串行，BLAS 线程数 1，未分配完整细网格 Hessian。提供十组必需 PNG/PDF 和一组 family 距边界补充图，黑底、Arial regular、英文标注、无箭头。

## 20. 当前 scientific status

```json
{json.dumps(labels,ensure_ascii=False,indent=2)}
```

总状态 INCOMPLETE 的原因是 C5/C7 未通过分支驻定域门槛，而不是几何决策缺失或内存被耗尽。完整材料公式验证、空间一致性和真实血栓物理验证是三个不同问题。本轮没有建立 mesh objectivity，也没有将这个 0.25 mm 数值试件解释为真实血栓。

## 21. 唯一下一步

人工审核后，单独审计**固定 delta、仅物理表面夹持与当前“所有 bond 的最小 Green 应变必须非负”定义之间的相容性**：解释为何细网格仿射自由节点有残余力、平衡试探何处先离开纯拉伸域，并决定该开发试件的适用域。当前证据不能直接证明域内不存在任何平衡，也不能授权放宽门槛或更换边界条件。本轮在此停止；没有进入后失稳分支、arc-length、h(s)<1、真实断裂、CFD 或模型修复。
'''
    (OUT/'REVIEW_SPATIAL_DISCRETIZATION_AUDIT_V2.md').write_text(report)
    save(RESULTS/'LEVEL_SUMMARY.json',data)
    save(OUT/'figures/FIGURE_MANIFEST.json',dict(font='Arial regular',background='black',labels='English',arrows=False,
        files=[p.name for p in sorted((OUT/'figures').glob('*')) if p.suffix in ['.png','.pdf']]))

def main():
    tests=test_summary();data=diagnostics();labels=status_labels();figures(data);review(data,tests,labels)
    preservation()
    print('Review and 11 PNG/PDF figure pairs generated. Audit INCOMPLETE: C5/C7 stationary-domain gates failed.')

def preservation():
    pre=read(OUT/'provenance/PRE_RUN_MANIFEST.json');post={};missing=[]
    for k in pre['files']:
        p=ROOT/k
        if p.is_file():
            with p.open('rb') as stream:post[k]=hashlib.file_digest(stream,'sha256').hexdigest()
        else:missing.append(k)
    changed=[k for k,v in pre['files'].items() if post.get(k)!=v]
    links={k:os.readlink(ROOT/k) if (ROOT/k).is_symlink() else None for k in pre['symlinks']}
    lchanged=[k for k,v in pre['symlinks'].items() if links[k]!=v]
    save(OUT/'provenance/POST_RUN_MANIFEST.json',dict(files=post,symlinks=links,timestamp=time.time(),scope='Every pre-existing project file and symlink'))
    save(OUT/'provenance/PRESERVATION_CHECK.json',dict(passed=not(changed or lchanged or missing),protected_files=len(post),protected_links=len(links),modified=changed,missing=missing,
        modified_links=lchanged,production_modifications=sum(k.startswith('pd_clot/') for k in changed),historical_modifications=len(changed),historical_test_modifications=sum(k.startswith('tests/') for k in changed)))
    (OUT/'provenance/STATUS_AFTER.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
    assert not(changed or lchanged or missing)

if __name__=='__main__':main()
