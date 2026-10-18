"""Physical-unit field generator and offline discretization audit; no dynamics."""
from common import *
from scipy.spatial import cKDTree
from scipy.ndimage import label,generate_binary_structure
from scipy.special import ndtri
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def correlated_phase(X,correlation_length_m,seed=20261018,strong_fraction=1/3,modes=256):
    """Stationary random Fourier approximation to a Gaussian covariance field.

    ell is in metres; this API has no desired fragment size/count parameters.
    The same coefficients define all grid samplings; no seed/outcome selection.
    """
    ell=np.broadcast_to(np.asarray(correlation_length_m,float),(3,)).copy()
    assert np.all(ell>0) and 0<strong_fraction<1
    rng=np.random.default_rng(seed);wave=rng.normal(size=(modes,3))/ell;phase=rng.uniform(0,2*np.pi,modes)
    values=np.empty(len(X))
    for i in range(0,len(X),4096):values[i:i+4096]=np.sqrt(2/modes)*np.cos(X[i:i+4096]@wave.T+phase).sum(1)
    return values>ndtri(1-strong_fraction),values,wave,phase

def audit(X,phase,domain,dx,horizon):
    nd=int(domain.max());records=[];b=cKDTree(X).query_pairs(horizon,output_type='ndarray');same=phase[b[:,0]]==phase[b[:,1]]
    distinct=phase[b[:,0]]&phase[b[:,1]]&(domain[b[:,0]]!=domain[b[:,1]])
    for c in range(1,nd+1):
        ids=np.flatnonzero(domain==c)
        if not len(ids):continue
        y=X[ids];other=X[phase&(domain!=c)];minimum=float(cKDTree(other).query(y)[0].min()) if len(other) else None
        inside=np.isin(b[:,0],ids)&np.isin(b[:,1],ids)
        records.append(dict(domain_id=c,nodes=len(ids),volume_m3=len(ids)*dx**3,diameter_volume_m=(6*len(ids)*dx**3/np.pi)**(1/3),
            center_span_m=np.ptp(y,axis=0),cell_envelope_span_m=np.ptp(y,axis=0)+dx,center_m=y.mean(0),radius_gyration_m=float(np.sqrt(np.mean(np.sum((y-y.mean(0))**2,axis=1)))),
            minimum_other_domain_node_distance_m=minimum,minimum_nominal_weak_gap_m=max(0.,minimum-dx) if minimum is not None else None,
            internal_original_bonds=int(inside.sum())))
    return dict(grid_pitch_m=dx,horizon_m=horizon,material_nodes=len(X),strong_nodes=int(phase.sum()),strong_fraction=float(phase.mean()),
        spatial_domains=nd,domains=records,bulk_original_bonds=len(b),cross_phase_bonds=int((~same).sum()),direct_interdomain_strong_bonds=int(distinct.sum()),
        gap_definition='Nearest reference cell-envelope gap, estimated as center distance minus dx; not a measured histological septum or geodesic matrix thickness')

def main():
    ref=load_npz(V2_O/'provenance/reference.npz');f=load_npz(V2_O/'BIPHASIC_FIELD.npz');n=int(ref['clot_count']);X=ref['X'][:n]
    baseline=audit(X,f['material_phase'][:n]==1,f['strong_domain_id'][:n],DX,HORIZON)
    b=ref['bonds'];p=f['material_phase'];dom=f['strong_domain_id'];bulk=ref['bond_types']==1
    baseline['same_original_graph_counts']={key:int(mask.sum()) for key,mask in {
        'bulk':bulk,'cross_phase':bulk&(p[b[:,0]]!=p[b[:,1]]),'live_cross_phase':bulk&(p[b[:,0]]!=p[b[:,1]])&~ref['initial_cut'],
        'direct_interdomain_SS':bulk&(p[b[:,0]]==1)&(p[b[:,1]]==1)&(dom[b[:,0]]!=dom[b[:,1]]),
        'live_direct_interdomain_SS':bulk&(p[b[:,0]]==1)&(p[b[:,1]]==1)&(dom[b[:,0]]!=dom[b[:,1]])&~ref['initial_cut']}.items()}
    candidates=[];folder=O/'offline_material_candidates';folder.mkdir(parents=True,exist_ok=True)
    # Illustrative correlation hypotheses bracket broad sub-millimetre region
    # organization, not a desired 20–200-node output. No candidate is admitted.
    for ell in [100e-6,200e-6,400e-6]:
        for dx in [25e-6,50e-6,100e-6]:
            axes=[np.linspace(-span/2,span/2,round(span/dx)+1) for span in [1.2e-3,1.2e-3,.6e-3]]
            xx=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3);strong,values,wave,phase=correlated_phase(xx,ell)
            domain,nc=label(strong.reshape(tuple(map(len,axes))),generate_binary_structure(3,1));domain=domain.ravel()
            # Full neighbor lists on the largest grid are unnecessary for a
            # sampling sensitivity diagnostic; evaluate bonds only on dx=50um.
            if dx==DX:
                q=audit(xx,strong,domain,dx,HORIZON)
                np.savez_compressed(folder/f'hypothesis_ell_{ell*1e6:.0f}um_dx50um.npz',X=xx,strong=strong,domain=domain,field=values,wave_vectors_per_m=wave,phases=phase)
            else:
                sizes=np.bincount(domain)[1:];q=dict(grid_pitch_m=dx,material_nodes=len(xx),strong_fraction=float(strong.mean()),spatial_domains=int(nc),
                    domain_volume_m3=sizes*dx**3,bond_audit='NOT_COMPUTED_ON_REFINED_SAMPLING; no refined mechanical run')
            q.update(correlation_length_m=ell,seed=20261018,status='HYPOTHETICAL_OFFLINE_ONLY',formal_admission=False,ell_over_dx=ell/dx,ell_over_horizon=ell/HORIZON)
            candidates.append(q)
    provenance=dict(status='MICROSTRUCTURE_SCALE_UNCALIBRATED',formal_material='EXACT_FROZEN_V2_10_DOMAIN_BASELINE',baseline_field_sha256=sha(V2_O/'BIPHASIC_FIELD.npz'),
        evidence=[dict(source='Liu et al., Fracture mechanics of blood clots: Measurements of toughness and critical length scales',local_file='/home/lzy/projects/Fracture mechanics of blood clots - Measurements of toughness and critical length scales.pdf',section='3.4',scale_m=[.1e-3,1e-3],object='Fracture dissipative/process-zone and nonlinear elastic scales',not_equivalent_to='RBC/fibrin phase correlation length'),
        dict(source='Non contrast enhanced volumetric histology of blood clots through high resolution propagation-based X-ray microtomography (2022)',doi='10.1038/s41598-022-06623-8',url='https://www.nature.com/articles/s41598-022-06623-8',archive='literature/PMC8854637.xml',object='Compact RBC-rich aggregates in dried retrieved stroke-clot sample 1',scale_description='hundreds of micrometres; porous boundaries tens of micrometres',numerical_correlation_length_m=None,limits='Qualitative regional-size observation in a processed specimen; not a calibrated strong fibrin/platelet-domain covariance or universal distribution'),
        dict(source='Same 2022 microCT study',object='Individual red blood cells',scale_m=[4e-6,7e-6],imaging_voxel_m=1.23e-6,not_equivalent_to='Mechanical material-domain correlation length'),
        dict(source='Karmakar et al., Towards a Model of Thrombus Embolization (2025)',doi='10.1002/cnm.70118',local_file='/home/lzy/projects/Towards a Model of Thrombus Embolization  Structural Response and Failure of.pdf',sections=['4.2','4.3'],finding='Discusses horizon interpretation and lack of systematic structural/fracture scale data; material model in that paper is homogeneous',numerical_correlation_length_m=None)],
        conclusion='Evidence supports heterogeneous regions and several distinct scales, but does not identify the phase covariance/threshold mechanics for this clot. No physical domain diameter is fitted from fragment targets.',
        candidate_rule='Gaussian random Fourier field in physical metres; fixed seed; fixed theoretical strong fraction threshold; no prescribed domain count or node count; the 100/200/400um ell values are uncalibrated hypotheses, not measured diameters.',baseline=baseline,candidates=candidates)
    dump(O/'DOMAIN_SCALE_PROVENANCE.json',provenance)
    text='''# 物理材料尺度核查\n\n状态：MICROSTRUCTURE_SCALE_UNCALIBRATED。正式 V3 使用原样冻结的 V2 10 强区，不把目标 20–200 节点反推成材料团尺寸。\n\nLiu 的 0.1–1 mm 是断裂过程区／非线性变形尺度，不是组织相区相关长度。Karmakar 2025 的 horizon 讨论也没有给出当前双相材料的实测相关长度。\n\n2022 年 microCT/SEM 原始研究在一个处理后的取出血栓中观察到数百微米级 RBC 聚集区和数十微米级多孔边界；单个 RBC 约 4–7 µm。前者是 RBC 区域观察，并非本模型强纤维蛋白区的已标定直径或三维协方差。该研究也展示分层和交织结构，不能把所有血栓表示为相同尺寸孤岛。来源：https://www.nature.com/articles/s41598-022-06623-8 。开放全文 XML 已归档，测量对象与局限见 DOMAIN_SCALE_PROVENANCE.json。\n\n当前 dx=50 µm 无法解析单个 RBC，也不足以稳定解析数十微米边界；horizon=150.5 µm 是当前离散非局部参数，不能当成组织学证据。原 10 区的体积等效直径约 401 µm、中心跨度 250–350 µm，最小参考单元间隙约 50 µm；存在跨相及跨域直接非局部连接。它们不等同于最终游离块。完整节点体积 V=dx³ 保持，因此端点取样体积不等于中心包围盒体积。\n\n新增生成器以米制相关长度 ell、各向异性、强相概率和固定 seed 为输入，不接受目标碎块节点数／目标块数。离线展示 ell=100/200/400 µm 三种假设以及同一连续随机场的 dx=25/50/100 µm 取样差异；ell 是协方差尺度，不能直接叫作团直径。所有候选仅用于结构与分辨率分析，不参与正式动力学，也不按碎裂结果选择 seed。\n\n图中材料相区只按参考六邻接识别，原生断裂始终采用全部存活原始材料键。区间弱基质厚度报告为参考单元包围盒最近间隙估计，而非组织学测量或真实最薄曲折通路。\n'''
    (O/'PHYSICAL_DOMAIN_SCALE_ASSESSMENT.md').write_text(text)
    fig=plt.figure(figsize=(16,10),layout='constrained');ax=fig.add_subplot(221,projection='3d');ids=np.flatnonzero(f['material_phase'][:n]==1)
    ax.scatter(*(X[ids]*1e6).T,c=f['strong_domain_id'][ids],cmap='tab10',s=6,depthshade=False);ax.view_init(17,-64);ax.set_box_aspect((2,2,1));ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',title='Formal baseline retained: 10 frozen regions')
    ax=fig.add_subplot(222);ds=baseline['domains'];ax.bar(np.arange(1,11),[q['diameter_volume_m']*1e6 for q in ds],color='#39869c');ax.axhline(HORIZON*1e6,color='#b64a36',label='numerical horizon');ax.set(xlabel='Frozen region ID',ylabel='Volume-equivalent diameter (µm)',title='Region size is not fragment size');ax.legend()
    ax=fig.add_subplot(223)
    for ell in [100e-6,200e-6,400e-6]:
        qs=[q for q in candidates if q['correlation_length_m']==ell];ax.plot([q['grid_pitch_m']*1e6 for q in qs],[q['spatial_domains'] for q in qs],'-o',label=f'hypothesis ell={ell*1e6:.0f} µm')
    ax.set(xlabel='Sampling pitch (µm)',ylabel='Six-neighbor strong-region count',title='Same physical random field, different sampling');ax.legend()
    ax=fig.add_subplot(224);ax.axis('off');ax.text(0,1,'Independent evidence and resolution\n\n• RBC: 4–7 µm; unresolved at dx=50 µm\n• RBC aggregates: broadly hundreds of µm\n• Porous boundaries: broadly tens of µm\n• Fracture process / nonlinear zone: 0.1–1 mm\n\nThese describe different physical objects.\nNo measured phase covariance is available for this clot.\n\nFormal material unchanged.\nCandidates: OFFLINE HYPOTHESES ONLY.\nMICROSTRUCTURE_SCALE_UNCALIBRATED',va='top',fontsize=13,linespacing=1.5)
    fig.suptitle('Material length scales: independent evidence, frozen baseline, and sampling limits',fontsize=19)
    (F/'figures').mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(F/'figures'/f'FIG_02_MATERIAL_LENGTH_SCALE_ANALYSIS.{ext}',dpi=180)
    plt.close(fig);print('PHYSICAL_SCALE_AUDIT_COMPLETE',baseline['same_original_graph_counts'],flush=True)

if __name__=='__main__':main()
