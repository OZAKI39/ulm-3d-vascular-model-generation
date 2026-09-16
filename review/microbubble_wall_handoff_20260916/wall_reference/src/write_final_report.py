"""Generate conclusions from frozen audit evidence, without package evaluations."""
from pathlib import Path
import csv,json,math,hashlib,subprocess
R=Path(__file__).resolve().parents[1]
W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['work'])
def read(name):return json.loads((R/name).read_text())
def save(name,obj):(R/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def csvout(name,rows):
    with (R/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
S=read('validation/MATRIX_AUDIT_SUMMARY.json');SC=read('WALL_REFERENCE_SCALING_AUDIT.json')
C=read('WALL_REFERENCE_CONVENTION.json');A=read('validation/INDEPENDENT_FINALIZER.json')
assert A['status']=='PASS'
rows=list(csv.DictReader((R/'WALL_REFERENCE_RAW_RESULTS.csv').open()))
main={label:{float(x['epsilon']):x for x in rows if x['implementation']==label and x['size']=='d50' and float(x['viscosity_Pa_s'])==.001} for label in S}
cmp=list(csv.DictReader((R/'WALL_REFERENCE_CROSS_COMPARISON.csv').open()))
cross={mode:max(float(r['symmetric_relative_difference']) for r in cmp if r['reference_A']=='RMBW_LUBRICATION' and r['mode']==mode) for mode in ['normal_TT','tangential_TT','RR_parallel','RR_normal','TR_x_Ty','RT_y_Fx','TR_y_Tx','RT_x_Fy']}
oneill=read('raw/RMBW_LUBRICATION_ONEILL_TABLE.json');on=[x for x in oneill if x['epsilon']<=.2]
on_errors={field:max(x[field] for x in on) for field in ['Fstar_relative_difference','Gstar_relative_difference']}
q=main['RMBW_LUBRICATION'];q0=q[.001];q1=q[.002]
log_slope=lambda col:(float(q0[col])-float(q1[col]))/math.log(2)
support={
    'normal':{'status':'PASS','reference':'Brenner series, independently evaluated from accessed Ascoli Eq.(29)',
        'relative_resistance_error_all_14_gap_max':S['RMBW_LUBRICATION']['max_normal_error_vs_Brenner'],
        'relative_resistance_error_epsilon_0p001':S['RMBW_LUBRICATION']['normal_error_at_smallest_gap'],
        'R_over_bulk_times_epsilon_at_0p001':S['RMBW_LUBRICATION']['normal_R_epsilon_at_smallest_gap'],
        'mobility_monotonic_with_gap':S['RMBW_LUBRICATION']['normal_monotonic_with_gap']},
    'tangential':{'status':'SUPPORTED_NEAR_FIELD_WITH_TABLE_COMPARISON','source':'Sprinkle Table I tracing GCB; independent ONeill1964 translation table',
        'log_slope_R_over_bulk_epsilon_0p001_to_0p002':log_slope('Rhat_tangential_TT'),'leading_log_coefficient':8/15,
        'ONeill_near_table_max_relative_Fstar_difference':on_errors['Fstar_relative_difference']},
    'RR_parallel':{'status':'PARTIAL_THEORY_SUPPORT','log_slope_R_over_bulk_epsilon_0p001_to_0p002':log_slope('Rhat_RR_parallel'),
        'leading_log_coefficient':.4,'source':'Sprinkle Table I and Liu-Prosperetti Eq.(5.1)',
        'constant_discrepancy':{'Sprinkle_RMBW':.3817,'LiuProsperetti':.3709,'status':'UNRESOLVED_CONSTANT_DISCREPANCY'},
        'independent_finite_gap_accuracy_bound':'UNVERIFIED'},
    'RR_normal':{'status':'CONTACT_LIMIT_SUPPORTED','reference':'Liu-Prosperetti Eq.(4.1)',
        'expected_contact_limit_R_over_bulk':1.202056903159594,
        'R_over_bulk_epsilon_0p001':float(q0['Rhat_RR_normal']),
        'independent_finite_gap_accuracy_bound':'UNVERIFIED'},
    'TR_RT':{'status':'NEAR_FIELD_SOURCE_AND_INDEPENDENT_TRANSLATION_TABLE_SUPPORTED',
        'primary_reciprocity':'PASS','ONeill_near_table_max_relative_Gstar_difference':on_errors['Gstar_relative_difference'],
        'original_author_linear_fit_retained':True,'new_fit_performed':False,
        'independent_continuous_gap_accuracy_bound':'UNVERIFIED','far_epsilon_above_9p018296':'TR_RT_ZERO_FALLBACK_NOT_HIGH_ORDER_REFERENCE'},
    'comparison_definition':'abs(A-B)/max(abs(A),abs(B)); signed coefficients; zeros/zeros=0',
    'cross_package_max_differences':cross,
    'accuracy_gates':'No universal 1% or agreement-based pass; reported differences and source scope retained'}
save('validation/THEORY_SUPPORT_AUDIT.json',support)
unverified=['RR_parallel finite-gap absolute accuracy; literature constant 0.3817 vs 0.3709 unresolved',
    'RR_normal finite-gap absolute accuracy outside independently supported contact asymptotic',
    'TR/RT continuous-gap accuracy beyond accessed ONeill table points and near asymptotic',
    'High-order TR/RT/RR at h/a > 9.018296 (upstream zero-coupling/bulk-rotation fallback)']
selection={'status':'PASS_WITH_LIMITATIONS','primary':'RigidMultiblobsWall Lubrication single-sphere total resistance reference',
    'primary_API':'Lubrication_Class.ResistCOO_wall(Sup_if_true=True) plus original subtracted bulk diagonal, six impulse solves',
    'secondary':'PyStokes wallBounded diagonal far-field diagnostic ONLY; complete matrix rejected for nonreciprocity',
    'additional_scalar_reference':'Brenner independently evaluated series; ONeill1964 translation table',
    'sphere_semianalytical_disposition':'REJECT_FULL_MATRIX; retain normal scalar and all failed outputs as diagnostics',
    'tested_h_over_a':C['gaps_h_over_a'],'recommended_near_field_h_over_a':[.001,.2],
    'qualified_reference_tested_points_h_over_a':[e for e in C['gaps_h_over_a'] if e<=5],
    'range_interpretation':'0.001–5 is a qualified reference interval on tested points, not a blanket independent error guarantee for every mode or every intervening gap',
    'far_field_exclusion':'At h/a 10 and 20 retain outputs, but do not use zero TR or bulk RR as a high-order reference',
    'radii_verified_m':C['radii_m'],'unverified_modes':unverified,'GPL_role':'REFERENCE_ONLY',
    'criteria':{'source_theory_traceability':'Classical sources, 2562-blob table and author fit explicitly separated',
        'near_wall_coverage':'All modes available in original lubrication path', 'normal':'Brenner check PASS',
        'tangential':'Log behavior and ONeill table support','rotation':'Two distinct axes; partial independent accuracy support',
        'coupling':'Source and ONeill translation torque support; reciprocal',
        'reciprocity':'PASS','positivity':'PASS','stability':'126/126 finite; max scaled condition <836',
        'assumptions':'Rigid no-slip sphere, single no-slip plane, Stokes flow'},
    'selection_based_on_speed_or_installation_ease':False,'production_integration':'NOT_STARTED'}
save('PRIMARY_REFERENCE_SELECTION.json',selection)
save('validation/STAGE_GATES.json',{
    'status':'PASS_WITH_LIMITATIONS', 'gate_scope':'Final selected primary is RMBW_LUBRICATION; package-specific failures are not erased',
    'RMBW_SOURCE_PROVENANCE':'PASS','RMBW_INSTALL':'PASS','RMBW_REFERENCE_RUN':'PASS',
    'PYSTOKES_SOURCE_PROVENANCE':'PASS','PYSTOKES_INSTALL':'PASS','PYSTOKES_SHORT_TEST':'PASS','PYSTOKES_REFERENCE_RUN':'PASS',
    'UNIT_CONVENTION':'PASS','6DOF_MATRIX_EXTRACTION':'PASS','RECIPROCITY_AUDIT':'PASS_PRIMARY_ONLY',
    'POSITIVITY_AUDIT':'PASS_PRIMARY_ONLY','RADIUS_SCALING_AUDIT':'PASS_PRIMARY_ONLY','VISCOSITY_SCALING_AUDIT':'PASS',
    'NEAR_WALL_BEHAVIOR_AUDIT':'PASS_PRIMARY_WITH_MODE_LIMITATIONS','PRIMARY_REFERENCE_SELECTED':'YES',
    'VISUALIZATION_GENERATION':'PASS','INDEPENDENT_ALL_MODE_FINITE_GAP_ACCURACY':'UNVERIFIED',
    'retained_failures':{'RMBW_SPHERE':['negative modes at five near gaps','radius scaling','planar symmetry','unguarded table extrapolation'],
                        'PYSTOKES':['TR/RT reciprocity','normal and tangential lubrication limits']},
    'unverified_modes':unverified})

# Machine-readable mode map separates API availability, matrix validity and accuracy.
validity=[]
for label in S:
 for e,row in main[label].items():
  for mode in ['normal_TT','tangential_TT','TR_RT','RR_parallel','RR_normal']:
   if label=='RMBW_LUBRICATION':
    scope='QUALIFIED_NEAR_FIELD_REFERENCE' if e<=.2 else ('QUALIFIED_TABULATED_REFERENCE' if e<=5 else 'LOW_ORDER_FALLBACK')
    if mode.startswith('RR') or mode=='TR_RT':scope+='; INDEPENDENT_FINITE_GAP_ERROR_BOUND_UNVERIFIED'
   elif label=='PYSTOKES':scope='LOW_ORDER_DIAGONAL_FAR_FIELD_DIAGNOSTIC' if mode!='TR_RT' else 'FAIL_RECIPROCITY'
   else:scope='INVALID_FULL_MATRIX; NORMAL_SCALAR_DIAGNOSTIC_ONLY' if mode=='normal_TT' else 'INVALID_FULL_MATRIX'
   validity.append({'implementation':label,'epsilon':e,'mode':mode,'API_SUPPORT':'SUPPORTED','finite':row['finite'],
     'full_matrix_valid':row['physical_matrix_valid'],'scope':scope,'unavailable':'NO'})
csvout('validation/MODE_VALIDITY_MAP.csv',validity)

initial_attempt={'status':'FINALIZER_CHECK_CORRECTED_WITHOUT_NUMERICAL_RERUN',
    'initial_assertion':'all rows have no warnings',
    'observation':'First RMBW_SPHERE import records deprecated imp and unclosed table-file warnings; 1/126 sphere rows, two warnings',
    'correction':'Preserve and report warnings. The user requires warning recording, not suppression.',
    'raw_output_changed':False,'package_evaluations_repeated':0}
save('provenance/FINALIZER_INITIAL_CHECK_RECORD.json',initial_attempt)
downloads=read('theory/ADDITIONAL_DOWNLOAD_PROVENANCE.json')
for d in downloads:
 if d['name']=='LiuProsperetti2010':
  d.update(status='WEB_FULL_TEXT_VERIFIED_LOCAL_DOWNLOAD_FAILED',local_download_status='HTTP403',
    web_access_status='PASS',verified_equations=['4.1','5.1','5.2','5.3','5.4'],pdf_pages_one_based=[10,13],
    constant_discrepancy_record='Eq5.1 0.3709 vs Sprinkle TableI 0.3817 retained, unresolved')
save('theory/ADDITIONAL_DOWNLOAD_PROVENANCE.json',downloads)
save('provenance/FINAL_BUILD_ENVIRONMENT_NOTES.json',{
    'actual_Eigen_include':'/usr/include/eigen3','actual_Eigen_debian_package':'libeigen3-dev 3.4.0-4build0.1',
    'environment_Eigen_not_used_for_this_compile':'conda Eigen 5; -I/usr/include/eigen3 was explicit',
    'RMBW_Boost':'1.82.0 Python38/Numpy38','CPU_threads':1,'GPU_USED':False,
    'Vast_access_used':False,'final_python_manifests':['provenance/wallref_rmbw_pip_freeze_FINAL.txt','provenance/wallref_pystokes_pip_freeze_FINAL.txt'],
    'PyMuPDF_added_for_QA_only':'1.24.10; after first environment freeze; not used in package numeric runs'})

summary={
 'WALL_HYDRODYNAMICS_REFERENCE_AUDIT':'PASS_WITH_LIMITATIONS',
 'RMBW_REPOSITORY':'https://github.com/stochasticHydroTools/RigidMultiblobsWall',
 'RMBW_COMMIT':read('RIGIDMULTIBLOBSWALL_PROVENANCE.json')['commit'], 'RMBW_LICENSE':'GPL-3.0',
 'RMBW_INSTALL':'PASS','RMBW_REFERENCE_IMPLEMENTATION':'LUBRICATION',
 'RMBW_NORMAL_TRANSLATION':'SUPPORTED; INDEPENDENT_BRENNER_CHECK_PASS',
 'RMBW_TANGENTIAL_TRANSLATION':'SUPPORTED; NEAR_FIELD_THEORY_AND_TABLE_SUPPORTED',
 'RMBW_ROTATION_PARALLEL':'SUPPORTED_WITH_LIMITATIONS', 'RMBW_ROTATION_NORMAL':'SUPPORTED_WITH_LIMITATIONS',
 'RMBW_TRANSLATION_ROTATION_COUPLING':'SUPPORTED_WITH_LIMITATIONS; RECIPROCITY_PASS',
 'PYSTOKES_REPOSITORY':'https://github.com/rajeshrinet/pystokes',
 'PYSTOKES_COMMIT':read('PYSTOKES_PROVENANCE.json')['commit'], 'PYSTOKES_VERSION':'2.3.3','PYSTOKES_LICENSE':'MIT',
 'PYSTOKES_INSTALL':'PASS','PYSTOKES_SHORT_TEST':'PASS (13 tests)',
 'PYSTOKES_NORMAL_TRANSLATION':'API_SUPPORTED; NEAR_WALL_ASYMPTOTIC_FAIL',
 'PYSTOKES_TANGENTIAL_TRANSLATION':'API_SUPPORTED; NEAR_WALL_ASYMPTOTIC_FAIL',
 'PYSTOKES_ROTATION_PARALLEL':'API_SUPPORTED; LOW_ORDER_APPROXIMATION',
 'PYSTOKES_ROTATION_NORMAL':'API_SUPPORTED; LOW_ORDER_APPROXIMATION',
 'PYSTOKES_TRANSLATION_ROTATION_COUPLING':'API_SUPPORTED; RECIPROCITY_FAIL',
 'PRIMARY_WALL_REFERENCE':selection['primary'], 'SECONDARY_WALL_REFERENCE':selection['secondary'],
 'REFERENCE_ASSUMPTION':C['assumption'], 'GAP_RANGE_TESTED_H_OVER_A':'0.001 ... 20 (14 frozen points)',
 'PRIMARY_REFERENCE_VALID_H_OVER_A_RANGE':'0.001–5 QUALIFIED; near-field focus 0.001–0.2; finite-gap RR/TR error bounds remain unverified',
 'MAX_PRIMARY_RECIPROCITY_ERROR':S['RMBW_LUBRICATION']['reciprocity_max_relative'],
 'MAX_PRIMARY_RECIPROCITY_ERROR_DEFINITION':'relative; absolute SI cross-block maximum = 0.001953125 (roundoff on ~1e13 entries)',
 'MIN_PRIMARY_SCALED_MOBILITY_EIGENVALUE':S['RMBW_LUBRICATION']['min_scaled_eigenvalue'],
 'RADIUS_SCALING_AUDIT':'PASS (primary); RMBW_SPHERE=FAIL', 'VISCOSITY_SCALING_AUDIT':'PASS (all three implementations)',
 'NORMAL_NEAR_WALL_ASYMPTOTIC':'PASS (primary)', 'TANGENTIAL_NEAR_WALL_REFERENCE':'SUPPORTED_BY_THEORY_AND_ONEILL_TABLE',
 'ROTATIONAL_REFERENCE':'PARTIAL_INDEPENDENT_SUPPORT; FINITE_GAP_ACCURACY_UNVERIFIED',
 'TRANSLATION_ROTATION_REFERENCE':'NEAR_FIELD_SUPPORTED; FAR_FIELD_LIMITED',
 'CROSS_PACKAGE_NORMAL_MAX_DIFFERENCE':cross['normal_TT'],
 'CROSS_PACKAGE_TANGENTIAL_MAX_DIFFERENCE':cross['tangential_TT'],
 'CROSS_PACKAGE_ROTATION_MAX_DIFFERENCE':max(cross['RR_parallel'],cross['RR_normal']),
 'CROSS_PACKAGE_TR_MAX_DIFFERENCE':max(cross[m] for m in cross if m.startswith('TR') or m.startswith('RT')),
 'CROSS_PACKAGE_DIFFERENCE_DEFINITION':support['comparison_definition'], 'UNVERIFIED_MODES':unverified}
vis=read('VISUALIZATION_PROVENANCE.json')
for item in vis['artifacts']:summary[Path(item['file']).stem]=str(R/item['file'])
summary.update(VISUALIZATION_GENERATION='PASS',HUMAN_VISUAL_REVIEW='PENDING',
    PALABOS_BASELINE_MODIFIED='NO',LAMMPS_BASELINE_MODIFIED='NO',BUBBLE_BUBBLE_BASELINE_MODIFIED='NO',SONOVUE_SAMPLER_MODIFIED='NO',
    PRODUCTION_WALL_MODEL_IMPLEMENTED='NO',MICROBUBBLE_WALL_HYDRODYNAMICS='PENDING',MICROBUBBLE_ADHESION='PENDING',
    MICROBUBBLE_WALL_EXCLUSION='PENDING',WALL_LOCAL_PLANE_APPROXIMATION='PENDING',CURVED_STL_WALL_MODEL='PENDING',
    RBC='OFF',GPU_PERFORMANCE_READY='NO',SOURCE_INTEGRITY='PASS',BASELINE_INTEGRITY='PASS (1392 local manifest files)',
    REPORT_DIR=str(R),REPORT=str(R/'WALL_HYDRODYNAMICS_REFERENCE_AUDIT_REPORT.md'),
    NEXT_STAGE='USER_REVIEW_BEFORE_MICROBUBBLE_WALL_HYDRODYNAMICS_V0')
save('FINAL_SUMMARY.json',summary)
terminal='\n'.join(f'{k} = '+(json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else str(v)) for k,v in summary.items())+'\n'
(R/'FINAL_TERMINAL_SUMMARY.txt').write_text(terminal)

table='| 实现 | finite | 互易性最大相对误差 | 最小归一化特征值 | 半径缩放 | 完整矩阵 |\n|---|---:|---:|---:|---|---|\n'
for label,s in S.items():table+=f"| {label} | {s['finite_rows']}/126 | {s['reciprocity_max_relative']:.6g} | {s['min_scaled_eigenvalue']:.9g} | {s['radius_status']} | {'通过必要门' if s['all_physical_matrix_valid'] else '拒绝'} |\n"
near='| h/a | RMBW LUB M⊥/MT0 | PyStokes M⊥/MT0 | RMBW LUB M∥/MT0 | PyStokes M∥/MT0 |\n|---|---:|---:|---:|---:|\n'
for e in [.001,.005,.01,.2,5,20]:
 r=q[e];p=main['PYSTOKES'][e]
 near+=f"| {e:g} | {float(r['Mhat_normal_TT']):.9g} | {float(p['Mhat_normal_TT']):.9g} | {float(r['Mhat_tangential_TT']):.9g} | {float(p['Mhat_tangential_TT']):.9g} |\n"
ctable='| 迁移率模式 | 最大相对差 | 说明 |\n|---|---:|---|\n'
for mode in ['normal_TT','tangential_TT','RR_parallel','RR_normal','TR_x_Ty','RT_y_Fx']:
 ctable+=f"| {mode} | {100*cross[mode]:.6f}% | {'保留符号；TR/RT 不互易' if mode.startswith(('TR','RT')) else '描述差异，不作为准确性门'} |\n"
radius='\n'.join(f"| {size} | {2*a*1e6:.15g} | {a:.15g} |" for size,a in C['radii_m'].items())
figures='\n\n'.join(f"[{Path(x['file']).name}]({x['file']})\n\n![{Path(x['file']).stem}]({x['file']})" for x in vis['artifacts'])
report=f'''# WALL_HYDRODYNAMICS_REFERENCE_AUDIT

**结论：PASS_WITH_LIMITATIONS。** 主参考选择 **RigidMultiblobsWall 的 Lubrication 单球壁面阻力实现**；PyStokes 仅保留为远场对角迁移率诊断。两个工具在极近壁面并不一致，且本次冻结的 PyStokes 完整矩阵不满足互易性。没有把软件一致性当成物理正确性，也没有修改上游源码或输出符号。

主参考建议用于已测试的 h/a=0.001–5 点，重点为 0.001–0.2 的近壁范围；这是有条件的参考资格，**不是该连续区间所有模式的独立精度认证**。h/a=10、20 的输出仍归档，但不能把其中零 TR 和 bulk RR 用作高阶参考。后续阶段等待人工审阅。

## 冻结问题与可复现输入

本地 WSL、CPU 单线程；未使用 Vast/GPU，也没有模拟 timestep。对象为不可变形 rigid no-slip sphere，壁面 z=0、流体 z>0，球心 (0,0,a+h)，Stokes regime。主黏度 0.001 Pa·s，另用 0.0005 和 0.002 Pa·s 核查缩放。没有引入密度、LBM dt、流场或真实 STL。

尺寸从现有 SonoVue contract 的 target quantiles 读取：

| 尺寸 | diameter (µm) | radius (m) |
|---|---:|---:|
{radius}

原 contract SHA256：`82ca08f51b3efef290c53beda53707b969b0e45e49ee6bb3bc636f3fb2f792d9`。三半径维度曲线缩放仅在这些点核查，不建立第二套手填尺寸。完整 gap、单位和门限见 [冻结约定](WALL_REFERENCE_CONVENTION.json)。每个实现 3 半径 × 14 gap × 3 黏度 = 126 个矩阵，共 **378 个首次正式矩阵**；另在 O’Neill 原表 8 个 alpha 点上对两个实现补充共 16 个明确标注的历史表格诊断，没有替换主 sweep。

## 软件身份、能力与实际实现

RMBW 官方仓库 [{summary['RMBW_COMMIT']}]({summary['RMBW_REPOSITORY']}/tree/{summary['RMBW_COMMIT']})，GPL-3.0，参考用途。独立 Python 3.8 / NumPy 1.24 / SciPy 1.10 环境；原始 CPU C++ binding 使用 GCC 13.3、Boost 1.82、系统 Eigen 3.4 构建。两个未初始化 submodule 的 gitlink SHA 已冻结；此单球路径不需要它们。

`sphere_mobility` 默认是低阶 blob 自项，并不是“best”单球函数。优先核查的 `sphere_best_mobility_known` 将 Huang 法向近似、Goldman/Faucheux 切向关系和 **162-blob 历史表的 RR/TR 样条**混用，不能称所有模式都是解析精确解。实测其完整矩阵失败：Goldman cross prefactor 在源码写成 `6*pi*a*2`，直接 SI 小半径时量纲缩放不成立；同时四个平面交叉项全同号，违反绕 z 旋转 90° 的几何协变。未修补。h/a≤0.02 的五个 gap 产生共 45 个负特征值结果；其 RR/TR 表在 h/a>9 无保护外推，h/a=20 的 RR_parallel 甚至超过 bulk。该入口不能作为完整主参考，法向标量结果可单独保留。

选择的 `Lubrication_Class.ResistCOO_wall(Sup_if_true=True)` 提供原始单球 wall **excess resistance**。wrapper 加回源码原先扣除的 isolated bulk 对角，恢复 total resistance，再解六个单位负载；这不是自行补齐缺失模式，也不是运行多体 RPY+润滑组合求解器。它在近壁使用公开渐近关系，中间区间使用仓库已有 **2562-blob 表格**；本轮没有生成新的 multiblob 球体、改变分辨率或拟合。较小 gap 的 TR 和 RR_parallel 包含原论文作者拟合的线性项，归档明确标注。

PyStokes 官方仓库 [{summary['PYSTOKES_COMMIT']}]({summary['PYSTOKES_REPOSITORY']}/tree/{summary['PYSTOKES_COMMIT']})，2.3.3、MIT；独立 Python 3.12.3 / NumPy 1.26.4 / SciPy 1.12.0。官方 short tests **13/13 PASS**，但这些测试不覆盖本轮发现的 TR/RT 互易性问题。wallBounded 的 TT/TR/RT/RR API 都存在，采用有限阶镜像/多极自项近似；N=1 坐标序列为 xyz，API 使用累加语义，因此每个单位负载均从零输出数组开始。

## 完整矩阵和单位核查

统一 q=(Vx,Vy,Vz,Ωx,Ωy,Ωz)，g=(Fx,Fy,Fz,Tx,Ty,Tz)。转矩绕球心，力、转矩为外加负载，采用同一个右手坐标系。MT0=1/(6πμa)，MR0=1/(8πμa³)，D=diag(√MT0,√MT0,√MT0,√MR0,√MR0,√MR0)。使用 B=D⁻¹MD⁻¹；q̂=D⁻¹q、ĝ=Dg 保持功率，R=D⁻¹ solve(B,I) D⁻¹。没有通过强制对称化改变矩阵。

RMBW lubrication 使用固定 native L*=1 µm、η*=1 mPa·s，实际半径和黏度仍随 case 改变。原因是原 COO API 会丢弃绝对值小于 1e-12 的 native entries；直接 SI 会丢失旋转项。独立反向核对表明 SI resistance 转换最大相对残差 **{A['SI_conversion_independent_checks']['rmbw_native_to_SI_resistance']:.3g}**。debye_cut=1e-5 在所有测试 gap 都不激活。转换不修改物理模型。

{table}
合法非零结构为六个对角元，以及 (0,4),(4,0),(1,3),(3,1)。x/y 两个切向块应满足一对交叉项与另一对相反号。RMBW lubrication 的 unexpected entries 和平面协变误差均为零；最高缩放条件数 **835.209**，所有 solve 残差 ≤{A['SI_conversion_independent_checks']['scaled_resistance_identity']:.3g}。原始 SI 交叉项约 10¹³，最大互易性绝对误差 0.001953125 对应相对误差 2.455e-16；不能单看不同量纲块混合下的大数字。

PyStokes 自项令 s=a/(a+h)，实测 M_TR[Vx,Ty]=−s⁴/(64πμa²)，M_RT[Ωy,Fx]=+s⁴/(64πμa²)，所以互易性相对误差为 **2**。这不是两个软件坐标差异；同一个包的两种互易实验已经异号。功率共轭坐标变换不能把非对称矩阵变成对称矩阵。对称部分特征值为正并不能抵消该失败。

全部 378 结果 finite、无异常、无 NaN。旧 sphere 第一次导入有两条警告（deprecated imp、table file 未关闭），保存在原始第一行和 finalizer；其余结果没有捕获警告。没有 retry 到正常，也没有删除极近 gap。所有完整矩阵和阻力都保存；对不合法 M 的代数逆只标记 `ALGEBRAIC_DIAGNOSTIC_ONLY_INVALID_M`，不得当作物理阻力。主参考半径缩放误差 {SC['RMBW_LUBRICATION']['max_radius_relative_error']:.3g}，黏度缩放误差 {SC['RMBW_LUBRICATION']['max_viscosity_relative_error']:.3g}；PyStokes 缩放也通过，旧 sphere 的半径缩放失败约 {SC['RMBW_SPHERE']['max_radius_relative_error']:.6g}。

## 经典理论、近壁结果及其限制

独立法向参考采用 [Ascoli 学位论文 Eq.(29)](https://thesis.caltech.edu/4428/3/Ascoli_EP_1988.pdf) 复现的 Brenner 级数，60/80 位精度分别求和；不是调用软件包公式。原 Brenner / GCB 期刊全文访问未确认，已单列 `UNVERIFIED_SOURCE_ACCESS`，未凭记忆编造完整解。

在 h/a=0.001，主参考 R⊥/(6πμa)=**1002.352831055907**，独立 Brenner=**1002.3533268805824**，相对差 **4.9466e-7**；全 14 gap 最大相对差 **0.829481%**。R⊥/Rbulk × epsilon 趋向 1，法向 mobility 随 gap 减小而下降。PyStokes 在该点仍有 M⊥/MT0=0.250250，有限阶公式的接触极限为 1/4，不满足 lubrication suppression；切向接触极限同样趋于 1/2，不能代替正确的对数阻力行为。

{near}
平行平移的主导对数系数在最小两 gap 实测为 **{support['tangential']['log_slope_R_over_bulk_epsilon_0p001_to_0p002']:.9g}**，与 8/15 一致。[O’Neill 1964 原表](https://doi.org/10.1112/S0025579300003508) 使用不旋转球的 F*、G*；本轮按 cosh(alpha) 恢复精确高度，并把流体转矩符号转换为外加负载。h/a≤0.2 的历史表点，F* 最大相对差 **{100*on_errors['Fstar_relative_difference']:.6f}%**，G* 最大相对差 **{100*on_errors['Gstar_relative_difference']:.6f}%**。表中远场小转矩有效位数有限，全部差异仍保留；GCB 曾修正的 Dean–O’Neill 旧旋转列没有用作裁判。

旋转轴平行于壁面的阻力具有对数发散；轴垂直于壁面则趋于有限值 ζ(3)≈1.2020569。两者分别记录，单球绕壁法向自转不等同于上一阶段的 bubble–bubble center-line twist。[Sprinkle 2020 Table I](https://arxiv.org/pdf/2005.06002) 的可读公式支持主参考的来源追踪；该论文自己的拟合/表格不能充当其实现的独立精度证据。[Liu–Prosperetti 2010](https://pages.jh.edu/aprospe1/publications/PapersPublished/Rotating/LiuJFM_2010.pdf) 支持法向轴接触极限和共同的平行轴主导对数项，但 Eq.(5.1) 的常数 0.3709 与 RMBW/Sprinkle 的 0.3817 不同，差异未裁定。没有改参考常数。

有限 gap 的 RR_parallel、RR_normal 和连续 gap 的交叉项尚无完整独立高精度误差界。h/a>9.018296，原 lubrication 分支使用 TT 低阶远场项、TR=0、RR=bulk；这解释 h/a=10、20 的差异，不能声称高阶耦合已验证。主参考在 h/a=5、10、20 的最大对角偏离 bulk 依次为 0.185201、0.102273、0.053571，趋势正确但不意味着 h/a=5 应已等于 bulk。

## 软件差异与最终选择

相对差定义为 abs(A−B)/max(abs(A),abs(B))，零/零记 0；先保留系数符号再作差。下表跨全部主网格/三半径，A=RMBW lubrication、B=PyStokes，μ=0.001。<2%、2–10%、>10% 仅为 agreement map 标签。

{ctable}
近壁 normal、tangential、RR_parallel 和耦合差异显著；RR_normal 相对接近也不能证明双方在有限 gap 都精确。远场对角项趋于一致，但 PyStokes 耦合仍存在互易性问题。主参考选择依据完整模式、来源、近壁理论、互易性、正性、稳定性和刚性球假设，未按速度或安装难度排名。

**未独立核定的模式/精度范围：**

'''+'\n'.join('- '+x for x in unverified)+f'''

因此本阶段采用 **PASS_WITH_LIMITATIONS**。建议后续生产 wall module 从公开理论独立实现，normal 直接对照 Brenner；近壁平移/耦合对照 O’Neill/GCB 来源与 RMBW lubrication；旋转先保留这里的常数与有限 gap 精度限制。RMBW 为 **REFERENCE_ONLY**，归档许可证保留，不将 GPL 源码直接并入声称独立许可的生产插件；PyStokes MIT 也没有集成进生产系统。

## 归档、复现与边界

两源码 checkout 和 PyStokes 构建副本的全部 tracked files、运行 binding 均核验。CSV/HDF5 的 378×36 个条目与首次输出逐项一致，独立单位检查通过。六组原有本地 baseline 共 **1392 个 manifest 文件**，运行前后 hashes 一致；没有访问或修改服务器。本次范围内 Palabos、LAMMPS、SonoVue、passive transport、bubble–bubble 数据未改动。详细记录见 [finalizer](validation/INDEPENDENT_FINALIZER.json)、[源码终审](provenance/SOURCE_INTEGRITY_FINAL.json) 和 [baseline 终审](provenance/BASELINE_IDENTITY_AFTER.json)。环境初始与最终 manifest 均保留，PyMuPDF 只用于读图。

所有输出、脚本、原始源码 tar、文献/图像和日志纳入最终 `SHA256SUMS`；复现命令见 [README](README_REPRODUCE.md)。没有生成生产 wall force、hard wall、adhesion 或真实 STL case。MICROBUBBLE_WALL_HYDRODYNAMICS、WALL_EXCLUSION、ADHESION、LOCAL_PLANE、CURVED_STL 均为 PENDING；RBC=OFF、GPU_PERFORMANCE_READY=NO，已有 bubble–bubble twist 仍 PENDING。

## 八张 QA 图

曲线由独立脚本直接读取冻结 CSV/HDF5。没有手工改曲线。`HUMAN_VISUAL_REVIEW=PENDING`；以下结果可供用户审阅。

{figures}

完整机器可读结论：[FINAL_SUMMARY.json](FINAL_SUMMARY.json)。完整终端摘要：[FINAL_TERMINAL_SUMMARY.txt](FINAL_TERMINAL_SUMMARY.txt)。
'''
(R/'WALL_HYDRODYNAMICS_REFERENCE_AUDIT_REPORT.md').write_text(report)
print('REPORT_GENERATED',str(R/'WALL_HYDRODYNAMICS_REFERENCE_AUDIT_REPORT.md'))
