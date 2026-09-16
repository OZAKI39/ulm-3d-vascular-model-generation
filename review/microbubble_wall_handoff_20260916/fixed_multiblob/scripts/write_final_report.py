from pathlib import Path
import json,csv,hashlib,numpy as np
R=Path(__file__).resolve().parents[1];S=json.loads((R/'validation/INDEPENDENT_FINALIZER.json').read_text());C=json.loads((R/'contracts/FIXED_MULTIBLOB_WALL_FEASIBILITY_CONTRACT.json').read_text())
with (R/'LEAKAGE_SCREENING.csv').open() as f:L=list(csv.DictReader(f))
with (R/'PATCH_CONVERGENCE.csv').open() as f:P=list(csv.DictReader(f))
perf=[r for r in S['performance'] if r['query_count']>0];p100=next(p for p in perf if p['N_wall']==100);pair=json.loads((R/'remote_raw/PAIR_ACTIVE_BASELINE_BENCHMARK.json').read_text());pairtime=float(np.median(pair['seconds_per_query']));ratio=p100['median_query_seconds']/pairtime
best=min(S['configurations'],key=lambda x:x['max_action']);pct=lambda x:f'{100*x:.3f}%'
patchcore=[r for r in P if float(r['patch_radius'])==12 and .01<=float(r['epsilon'])<=.2]
patchpass=sum(r['pass']=='True' for r in patchcore)
peak=max((p.get('max_rss_bytes') or 0) for p in S['performance'])/2**30
replay=json.loads((R/'validation/REMOTE_TO_WSL_INTEGRITY.json').read_text());immut=json.loads((R/'validation/IMMUTABILITY_AFTER.json').read_text());assert replay['status']=='PASS' and immut['status']=='PASS'
figs=list(R.glob('VIS_*.png'));vtps=list(R.glob('*.vtp'));assert len(figs)==12 and len(vtps)==4
def property_status(key):
 good=[p for p in S['configurations'] if p[key]]
 if not good:return 'FAIL_TESTED_CONFIGURATIONS'
 if len(good)==len(S['configurations']):return 'PASS_TESTED_CONFIGURATIONS'
 return 'PARTIAL_PASS; accepted subtests at beta '+','.join(str(b) for b in sorted({p['beta'] for p in good}))
final={
'FIXED_MULTIBLOB_WALL_FEASIBILITY_AUDIT':S['status'],
'TOOL_CAPABILITY':'PASS_PRESCRIBED_MOTION; large-matrix failure recorded separately',
'CORE_NEAR_WALL_ACCURACY':'FAIL',
'PRIMARY_IMPLEMENTATION':'Pecnut stokesian-dynamics + thin prescribed-motion/Schur adapter',
'PRIMARY_REPO':'https://github.com/Pecnut/stokesian-dynamics',
'PRIMARY_COMMIT':C['primary_commit'],'PRIMARY_LICENSE':'MIT',
'SECONDARY_IMPLEMENTATION':'NONE_RUNTIME; RMBW and RigidBodyIB capability audit only',
'WALL_REPRESENTATION':'FIXED_HARD_SPHERES','RMBW_REFERENCE_SHA256':C['table_sha256'],
'FREE_SPACE_TRANSLATION':'PASS','FREE_SPACE_ROTATION':'PASS',
'PLANAR_CONTINUOUS_WALL_APPROXIMATION':'FAIL_WITHIN_TESTED_RESOLUTIONS',
'SELECTED_LAYOUT':'NONE_ACCEPTED','SELECTED_BETA':'NONE_ACCEPTED','SELECTED_SPACING_RATIO':'NONE_ACCEPTED (tested 1.01, 1.05, 1.1)','SELECTED_LAYER_COUNT':'NONE_ACCEPTED',
'PATCH_CONVERGENCE':f'PARTIAL_COVERAGE; {patchpass}/{len(patchcore)} completed core comparisons pass',
'CORE_NEAR_WALL_RANGE':'0.01–0.2','CORE_NEAR_MAX_NORMAL_ERROR':S['core_max_normal'],'CORE_NEAR_MAX_TANGENTIAL_ERROR':S['core_max_tangent'],'CORE_NEAR_MAX_MATRIX_ACTION_ERROR':S['core_max_action'],'CORE_NEAR_MAX_RR_ERROR':S['core_max_RR'],'CORE_NEAR_MAX_TR_ERROR':S['core_max_TR'],
'ULTRA_NEAR_RANGE':'0.001–0.01','ULTRA_NEAR_WALL_ACCURACY':'FAIL; NOT_VALIDATED_FOR_NEAR_CONTACT_ADHESION',
'MAX_LATERAL_PHASE_SPREAD':S['max_core_TT_phase'],'LATERAL_PHASE_INVARIANCE':property_status('phase_pass'),'MAX_TANGENTIAL_ANISOTROPY':S['max_core_HEX_TT_anisotropy'],'LATTICE_ISOTROPY':property_status('anisotropy_pass'),
'GEOMETRIC_ROUGHNESS_OVER_A':S['max_surface_roughness_RMS'],'LEAKAGE_METRIC':f"ABS(normalized lower-probe Vz) range [{S['min_leakage']}, {S['max_leakage']}]",'LEAKAGE_CONTROL':('PASS_TESTED_CONFIGURATIONS' if S['leakage_layer_gate_pass'] else 'FAIL_LAYER_TREND; ABSOLUTE_PROXY_GATE_PASS') if S['leakage_absolute_gate_pass'] else 'FAIL_ABSOLUTE_PROXY_GATE',
'MAX_RECIPROCITY_ERROR':S['max_reciprocity'],'MIN_SCALED_EIGENVALUE':S['min_eigenvalue'],'MAX_CONDITION_NUMBER':S['max_condition'],'DISSIPATION':'PASS_COMPLETED_VALID_MATRICES',
'CYLINDER_TEST':'NOT_RUN_PLANAR_GATE_FAILED','CYLINDER_R_OVER_A_TESTED':'NONE','CURVED_WALL_SELF_CONVERGENCE':'UNVERIFIED','REAL_STL_PATCHES_TESTED':0,'REAL_STL_STATIC_FEASIBILITY':'UNVERIFIED_PLANAR_GATE_FAILED',
'MAX_WALL_OBJECT_COUNT_TESTED':2000,'MAX_WALL_OBJECT_COUNT_WITH_VALID_QUERY':S['max_wall_count_measured'],
'FIRST_QUERY_TIME_S':p100['first_including_precompute_seconds'],'REUSED_QUERY_TIME_S':p100['median_query_seconds'],'MEMORY_GB':peak*2**30/1e9,
'FIXED_WALL_FIRST_TO_REUSED_SPEEDUP':p100['first_including_precompute_seconds']/p100['median_query_seconds'],
'CURRENT_PAIR_STATIC_QUERY_TIME_S':pairtime,'FIXED_WALL_TO_CURRENT_PAIR_COST_RATIO':ratio,
'ONLINE_TIMESTEP_FEASIBILITY':'UNVERIFIED','OFFLINE_REFERENCE_FEASIBILITY':'UNVERIFIED_CURVED_NOT_RUN; current planar representation not qualified',
'GPU_BACKEND_TESTED':'NO','INDEPENDENT_FINALIZER':'PASS_WITH_EXECUTION_FAILURES_SEPARATELY_RETAINED','VISUALIZATION_GENERATION':'PASS; 12 PNG and 4 VTP, blocked stages explicitly labelled','HUMAN_VISUAL_REVIEW':'PENDING',
'PALABOS_BASELINE_MODIFIED':'NO','LAMMPS_CORE_MODIFIED':'NO','CURRENT_WALL_V0_MODIFIED':'NO','RMBW_SOURCE_MODIFIED':'NO','ADHESION':'OFF','RBC':'OFF','BUOYANCY':'OFF','LIFT':'OFF','REMOTE_TO_WSL_INTEGRITY':'PASS',
'REPORT':str(R/'FIXED_MULTIBLOB_WALL_FEASIBILITY_AUDIT_REPORT.md'),'LOCAL_RESULT_DIR':str(R),'NEXT_STAGE':'RECOMMENDED_ONLY — BEM_CURVED_WALL_FEASIBILITY_AUDIT'}
(R/'FINAL_TERMINAL_SUMMARY.txt').write_text('\n'.join(f'{k} = {v}' for k,v in final.items())+'\n');(R/'validation/FINAL_CLASSIFICATION.json').write_text(json.dumps(final,indent=2,ensure_ascii=False)+'\n')
res='\n'.join(f"|{r['beta']}|{r['core_queries']}|{pct(r['worst_action'])}|{pct(r['mean_action'])}|{pct(r['worst_TT_phase'])}|" for r in S['resolutions'])
def perf_line(p):
 status='VALID_QUERIES' if p['query_count'] else ('QUERY_FAILED' if p.get('wall_assembly_seconds') is not None else 'PREFLIGHT_BLOCKED')
 num=lambda v:'UNVERIFIED' if v is None else f'{v:.6g}'
 return f"|{p['N_wall']}|{status}|{num(p.get('wall_assembly_seconds'))}|{num(p.get('wall_factor_seconds'))}|{p['query_count']}|{num(p['median_query_seconds'])}|{num(p['max_rss_bytes']/2**30 if p.get('max_rss_bytes') else None)}|"
pt='\n'.join(perf_line(p) for p in sorted(S['performance'],key=lambda p:(p['N_wall'],bool(p.get('wall_assembly_seconds')))))

text=f'''这个方法不是把弯曲血管壁强行当成平面。它用很多固定不动的数值小球/水动力 blob 沿真实墙面排列。如果这些固定点足够密，它们可能共同表现得像连续的 no-slip wall。本轮首先用平面墙与已经验证的 RMBW 对照，检查这种离散墙是否真的趋近连续墙；只有平墙通过后，才测试曲壁和真实 vascular STL。

# 固定小球墙可行性核查 V0

**结论：{S['status']}。范围限定为本轮 Pecnut 固定硬球实现、所测分辨率和资源界限；这不是对全部 regularized-blob 方法或未测试极限的否定。** 主工具能力、自由球标度、有效矩阵的数值性质通过；连续平墙物理精度未通过。没有选择可投入生产的布局。

1. **是否趋近 RMBW？** 细化有改善，但所完成配置均未同时满足核心精度、相位、各向同性、截断和泄漏门槛。最小“最坏作用误差”的诊断配置为 β={best['beta']}、{best['layout']}、{best['layers']} 层，仍为 {pct(best['max_action'])}；它不是合格配置。
2. **需要什么 β / spacing / layers？** 尚未得到合格组合。正式间距比 1.01，诊断 1.05/1.1；HEX/SQUARE、1–3 层按内存可行性逐项记录。β=.0625 及可选 .03125 的固定 4a 墙面已超过本轮内存界限，未缩小墙面来冒充细化。
3. **相位敏感性多大？** 所测核心区 TT 最大相位跨度 {pct(S['max_core_TT_phase'])}；下表给出随细化的变化。每个晶胞包括 bead/bridge/pore 和 16 个规则位置，没有只挑球顶位置。
4. **人工各向异性多大？** HEX 核心区方向 0°/45°/90° 的 TT 最大跨度 {pct(S['max_core_HEX_TT_anisotropy'])}，RR 最大跨度 {pct(S['max_core_HEX_RR_anisotropy'])}。完整逐位置数据保留。
5. **泄漏是否可控制？** 冻结双探针指标的幅值实测范围 {S['min_leakage']:.6g}–{S['max_leakage']:.6g}，门槛 .05；当前状态 {final['LEAKAGE_CONTROL']}。这是有限墙面跨墙传递代理，包含绕边效应，不冒充无限墙渗透率。
6. **超近壁 .001–.01 是否可靠？** 不可靠；最大作用误差 {pct(S['ultra_max_action'])}，不能用于近接触黏附。
7. **圆柱曲壁是否自收敛？** 未运行，平墙门槛未通过。
8. **真实 STL patch 能否稳定构造？** 本轮未验证。原 STL 保留且 SHA 相同，未构造伪造试验点或改造几何。此前几何有效性指标与真实水动力误差的对应关系仍未回答。
9. **在线每步是否可承受？** 没有合格物理配置，且未冻结正式粒子数量/时步吞吐目标，故在线可承受性 UNVERIFIED。一个目标、100 个墙球、CPU 4 线程：首次含预计算 {p100['first_including_precompute_seconds']:.6g} s，复用查询中位数 {p100['median_query_seconds']:.6g} s；约为当前一个活跃粒子对静态算子的 {ratio:.6g} 倍成本。
10. **推荐什么？** 当前固定硬球实现不迁入正式墙模型，也不作为已验证曲壁离线参考。建议下一独立审查为 BEM_CURVED_WALL_FEASIBILITY_AUDIT；仅建议，未启动。真正 regularized-blob 后端也未被本轮测试排除。

## 工具与科学定义

主工具 [Pecnut/stokesian-dynamics](https://github.com/Pecnut/stokesian-dynamics)，提交 `{C['primary_commit']}`，MIT。完整能力证据见 `OPEN_SOURCE_TOOL_CAPABILITY_AUDIT.md`。使用原生自由空间 M∞、两球润滑 excess、11 DOF/球的应力偶耦合。固定墙 U/Ω/E 严格为零，目标球 E=0；六个力/力矩响应来自完整算子的规定运动列。标准 Schur 消元和 Cholesky 缓存不改变水动力数学。原生完整矩阵与缓存的对照原始数据在 H5 和 NPZ 内。墙球接触诊断的完整阻力矩阵条件数约 4.37×10⁵，主间距 1.01 时约 245；目标 6×6 约 4.27，原生/缓存误差仍在 4×10⁻¹⁶ 内。`WALL_WALL_LUBRICATION_DIAGNOSTIC.csv` 由独立 finalizer 从完整原始矩阵重算，不把高刚度壁球接触等同于目标响应不稳定。该消元结论仅适用于严格固定墙。

这是一层或多层 **FIXED_HARD_SPHERES**，不是 regularized blobs、BEM、LBM 网格或直接 no-slip 三角面。第一层球心 z=−a_w，名义墙 z=0，目标球心 z=a+h。没有有效墙高拟合，没有逐间隙 offset，没有补平孔隙。上游 M∞ 内部 s'<2.001 的近接触保护与 scalar 表最低 s'=2.00001 夹取均记录；主间距 1.01、正式目标间隙 ≥.001a 不借此进行墙面校准。接触墙球只用于原生/缓存一致性诊断。

μ=.001 Pa·s；使用冻结 SonoVue d10/d50/d90 的**半径**。9 个原生自由球 SI 矩阵覆盖三种半径和 μ×.5/1/2。无拟合 bulk 校正；TT、RR、TR/RT 由独立 finalizer 重新归一化核查。各静态无量纲响应还保存了三种尺寸的 SI 相似标度转换，转换数据不冒充三次独立大墙求解。

RMBW 表 SHA `{C['table_sha256']}`。核心比较范围 .01–.2；参考的有限间隙 RR/TR 资格限制、远区模型分支限制继续保留。参考表不是独立曲壁真值。Native Stokesian Dynamics 的两球 scalar 插值和有限多极截断也构成方法误差，因此“误差”代表当前实现与冻结平墙参考之差，不能全部归结为纯几何粗糙度。

## 细化、相位与有限墙面

|β|完成核心查询数|最坏矩阵作用误差|平均矩阵作用误差|最大 TT 相位跨度|
|---|---:|---:|---:|---:|
{res}

β=.5→.25 的平均误差改善为 {S['mean_action_improvement']:.6g}（无量纲绝对差值，即百分点除以 100）；冻结的三条件止损规则结果为 `{S['coarse_to_fine_stop']}`。没有事后改变 .05 阈值。β=.125 的表中平均值只包含可执行的层数，不能与全六配置平均值直接当作严格配对比较；`MATCHED_REFINEMENT.csv` 另给相同布局、层数、间隙和相位的配对比较。继续到能通过资源预检的细化配置；未完成的大模型在 `GEOMETRY_EXECUTION.csv`、`FINER_RESOLUTION_PREFLIGHT.json`、执行 receipt 中保留，不计为已测收敛。

平墙面半径按 4/6/8/12a 构造，补丁中心不随目标位置移动。完成的核心 8a→12a 比较中 {patchpass}/{len(patchcore)} 通过各分量门槛；未完成配置不参与此分母，也不因此授予最终 patch PASS。`PATCH_CONVERGENCE.csv` 保存每个间隙/相位的法向、切向、TT、RR、TR 变化。有限墙面变化小不等于连续墙误差小。

粗糙度由真实球面上包络计算，同时保存未覆盖竖直柱比例；没有对孔隙做插值填充。最大已计算表面 RMS/a={S['max_surface_roughness_RMS']:.6g}。第一接触高度按目标球与每个墙球的实际非重叠条件求解，允许名义 h/a 为负；这是粗糙墙几何结果，不是重新定义冻结间隙。

## 泄漏与矩阵性质

上方单位法向技术力驱动球、下方零力零力矩探针，墙全部固定；用同位置无墙双球响应归一化。上下探针半径均为 a，下探针随墙底面位置设置，归一化显式保留距离影响。原始跨墙响应允许方向反转，CSV 保留符号；筛漏强度和层数改善门槛使用幅值，不能借负号把恶化误判为通过。绝对传递很小与随层数单调改善是两个不同检查。`LEAKAGE_SCREENING.csv` 给出原始上下响应、层数比值和 spacing 扫描；不假装 native 工具给出了未计算的 Eulerian 速度场。

有效查询的最大互易误差 {S['max_reciprocity']:.6g}、最小缩放特征值 {S['min_eigenvalue']:.6g}、最大条件数 {S['max_condition']:.6g}，最小测试耗散 {S['minimum_dissipation']:.6g}，最大线性残差 {S['max_linear_residual']:.6g}。每个有效 6×6 矩阵用冻结 256 个混合动作及 10,000 个耗散向量独立核查；H5 保存每向量响应/误差，包含总阻力和 excess。数值矩阵正定不意味着连续墙近似正确。

## 计算成本与执行限制

AMD Ryzen 7 7800X3D，BLAS 4 线程，NumPy 1.26.4 / SciPy 1.11.4 / Numba 0.59.1。RTX 4090 存在但未测试 GPU 后端。墙 M∞ 可一次分解并复用；100 个目标位置使用冻结随机种子，完整位置和时间保留。

|墙球数|状态|装配 s|分解 s|有效查询数|查询中位数 s|进程峰值 RSS GiB|
|---:|---|---:|---:|---:|---:|---:|
{pt}

原生逐墙球交叉块调用、目标应力偶消元、标量润滑与矩阵求解均包含在 fixed-wall 查询时间中。当前粒子对基线复制了原 `rigid_math.cpp`，在同一 CPU 单独编译，活跃近场 gap=.05a、两个 d50 球，检查确有一个近场事件，5×10000 次算子装配中位数 {pairtime:.9g} s。它提取固定邻球的目标 6×6 块；不是两个引擎的完整 timestep 等价基准。最初 gap=.1a 落在原生 cutoff 边界的无活跃对计时作为控制数据保留，**不用于主成本比**。

2,000 球最初被保守内存预检拒绝；根据 1,000 球实测内存重新估算后，在不改变 20 GiB 上限的前提下实际执行一次。装配/分解完成，但三个查询报告特征值不收敛，进程以 -6 退出并打印 `free(): invalid size`。失败 H5、日志和 receipt 原样保留；未计入有效响应或成功查询速度。小规模本地/远端对照仍一致，现阶段大矩阵失败根因 UNVERIFIED，不能归因为墙面物理，也没有隐藏重跑。5,000/10,000 球按原资源门槛未执行。

## 核查、交付与边界

独立脚本 `scripts/finalize_fixed_multiblob_wall_audit.py` 从原始矩阵重算自由空间标度、RMBW 比较、作用误差、互易性、特征值、耗散、相位、各向异性、patch 差异、粗糙度和泄漏；不读取求解器 PASS 标签代替计算。CPU 大矩阵失败另列，不能被平均统计隐藏。

- 主数据：`FIXED_MULTIBLOB_WALL_AUDIT.h5`；不可变远端原始证据：`remote_raw/`。
- 必需 8 项 CSV 均已生成，另有 patch、接触、粗糙度、失败与执行范围表。
- 12 张 PNG 与 4 个 VTP 已生成。曲壁/真实 STL 图为明确标注未运行的占位图，对应 3 个 VTP 是带状态字段的空数据；没有伪造 patch 或轨迹。平墙 VTP 为失败算例示意，坐标按冻结 d50 半径转换为米。
- `VISUALIZATION_PROVENANCE.json` 记录脚本、数据、PNG、VTP SHA。HUMAN_VISUAL_REVIEW=PENDING。
- 远端生成清单，下载后实际 `sha256sum -c`：PASS。现有 wall V0、原参考表、SonoVue、pair baseline 和上游源码 SHA 核查 PASS。Palabos、LAMMPS、RBC、几何和正式源码均未修改。
- 没有动态 smoke、adhesion、buoyancy、lift、RBC、two-way Palabos coupling 或后续生产迁移。

重现入口见 `README_REPRODUCE.md`。完整终端摘要见 `FINAL_TERMINAL_SUMMARY.txt`。所有推荐仅供人工审阅，未自动进入下一阶段。
'''
(R/'FIXED_MULTIBLOB_WALL_FEASIBILITY_AUDIT_REPORT.md').write_text(text);print((R/'FINAL_TERMINAL_SUMMARY.txt').read_text())
