"""Assemble the causal audit, design obligations and provenance from saved data."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/baqp_near_boundary_kinematic_energy_audit_v1'
RESULT=ROOT/'results/baqp_near_boundary_kinematic_energy_audit_v1'
REPORT='REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md'
def read(name):return json.loads((OUT/name).read_text())
def write(name,value):(OUT/name).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def fmt(x):return f'{x:.8e}' if isinstance(x,(float,np.floating)) else str(x)
def table(headers,rows):
    def line(r):return '| '+' | '.join(fmt(x).replace('|',r'\|') for x in r)+' |\n'
    return line(headers)+line(['---']*len(headers))+''.join(line(r) for r in rows)
def norm(a):return np.sqrt(np.sum(a*a,axis=tuple(range(1,a.ndim))))


def main():
    cfg=read('PREDECLARED_CONFIG.json');full=read('CONTINUOUS_FULL_BALL_REFERENCE.json');lat=read('FULL_LATTICE_REFERENCE.json')
    source=read('SOURCE_DECOMPOSITION.json');energy=read('ENERGY_SOURCE_DECOMPOSITION.json')
    closures={v['field']:v for v in read('CAUSAL_IDENTITY_SUMMARY.json')};cross=read('AXIAL_SHEAR_CROSSCHECK_SUMMARY.json')
    resource=read('RESOURCE_USAGE.json');git=read('STARTING_GIT_STATE.json');m=np.load(OUT/'MOMENT_TENSORS.npz')
    assert all(v['PASS'] for v in closures.values()) and cross['PASS'] and lat['PASS']
    assert '10 passed' in (OUT/'TARGETED_TESTS.log').read_text()
    write('TARGETED_TESTS.json',{'command':'PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B -m pytest tests/test_near_boundary_kinematic_energy_audit_v1.py -q -p no:cacheprovider',
        'passed':10,'failed':0,'seconds':.67,'scope':'Only the ten requested focused tests; full pytest not run','log':'TARGETED_TESTS.log'})
    pdf=Path('/home/lzy/projects/temp_storage/nonlocal_boundary_literature_20261005/shojaei.pdf')
    digest=hashlib.sha256(pdf.read_bytes()).hexdigest()
    assert digest=='f9109effa980245a9a7f9f927cc8a18a27cd928a1c41e1507bde643ec4439958'
    literature={'title':'Surface and Boundary Corrections in Peridynamics Using Optimized Nodal Influence Weights',
        'authors':'Shojaei, Ongaro, Hermann, Silling, Trovalusci, Cyron','year':2026,'doi':'10.1002/nme.70385',
        'primary_url':'https://tore.tuhh.de/bitstreams/8fda7c78-cb40-401e-b8ce-206cd986fafb/download',
        'local_pdf':str(pdf),'sha256':digest,'size':pdf.stat().st_size,
        'access':'Read previously retrieved primary publisher PDF and local extracted text; web open this round exceeded 10 MiB tool limit. No new download or repository copy.',
        'uploaded_copy_identity':'Not byte-verified against a user-uploaded PDF; prior primary-source provenance retained.',
        'relevant_locations':['PDF page 3: PD-to-PD rationale','PDF page 9, Eq. 33: full symmetric odd moments','PDF page 10, Eq. 40: half-horizon odd moments','Section 4, Eqs. 44-56: formulation-specific optimized influence weights'],
        'scope':'Motivation only. The BAQP M3/K and Ftilde projection identities are derived here from unchanged local BAQP code. BB/OSB weights are not BAQP targets.'}
    write('LITERATURE_PROVENANCE.json',literature)
    levels=np.rint(m['distance_face_m']/cfg['h_m']).astype(int);factor=cfg['beta_s']/cfg['L_shear_m']/cfg['gradient_scale']
    summaries=[]
    for group in ['bulk-like','face-near','edge-near','corner-near']:
        q=m['geometry_class']==group
        summaries.append([group,int(q.sum()),int(m['family_size'][q].min()),int(m['family_size'][q].max()),float(m['condition_K'][q].max()),float(norm(m['M3'][q]).max()),float(norm(m['G_physical'][q]).max()*factor)])
    with (OUT/'GEOMETRY_CLASS_SUMMARY.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['class','nodes','family_min','family_max','condition_K_max','M3_norm_max_m6','predicted_bias_max_scaled']);w.writerows(summaries)
    # Signed representative tensors: select exact physical lattice points, no averaging across opposing faces.
    X=m['X'];B=9*cfg['h_m'];examples=[]
    for name,coord in [('bottom_z_face',[0,0,-B]),('top_z_face',[0,0,B]),('positive_x_face',[B,0,0]),('positive_xy_edge',[B,B,0]),('positive_corner',[B,B,B])]:
        k=int(np.argmin(norm(X-np.array(coord))));assert np.linalg.norm(X[k]-coord)<1e-15
        examples.append({'location':name,'node':k,'X_m':X[k].tolist(),'M3_zz_m6':m['M3'][k,2,2,:].tolist(),'G_m':m['G_physical'][k].tolist(),'predicted_shear_row_scaled':(factor*m['G_physical'][k]).tolist()})
    q=levels==4
    geometry_extra={'signed_examples':examples,'depth_4h':{'nodes':int(q.sum()),'family_size_unique':np.unique(m['family_size'][q]).tolist(),
        'R2_max_scaled':float(norm(m['G_uniform'][q]).max()*factor),'R3_max_scaled':float(norm(m['G_physical'][q]).max()*factor)},
        'class_rule':cfg['node_geometry_rule'],'interpretation':'Complete retained neighbor count is not complete weighted symmetry when family intersects boundary half-volumes.'}
    write('GEOMETRY_CAUSAL_DETAILS.json',geometry_extra)
    statuses={
        'CONTINUOUS_FULL_BALL_M3':'ZERO','DISCRETE_FULL_LATTICE_M3':'NUMERICAL_ZERO','FULL_LATTICE_QUADRATIC_REPRODUCTION':'PASS',
        'FBAR_MOMENT_IDENTITY':'PASS','FBAR_BIAS_PRIMARY':'TRUNCATED_ODD_MOMENT','BOUNDARY_VOLUME_WEIGHTING_EFFECT':source['volume_contribution'],
        'K_ANISOTROPY_MODIFIER':{'continuous_reference_range':source['K_modifier_continuous_range'],'lattice_reference_range':source['K_modifier_lattice_range']},
        'FTILDE_PROJECTION_IDENTITY':'PASS','FTILDE_BIAS_PRIMARY':'PROJECTED_NODAL_RECONSTRUCTION_BIAS',
        'AXIAL_SHEAR_GEOMETRY_INVARIANCE':'PASS','PSI_PLUS_CAUSAL_RECONSTRUCTION':'PASS','PSI_PLUS_BIAS_PRIMARY':'KINEMATIC_RECONSTRUCTION',
        'GLOBAL_ENERGY_KINEMATIC_BIAS':{f:{c:energy[f][c]['kinematic_bias_J'] for c in ['intact','positive']} for f in energy},
        'GLOBAL_ENERGY_WEIGHT_QUADRATURE_BIAS':{f:{c:energy[f][c]['weight_quadrature_bias_J'] for c in ['intact','positive']} for f in energy},
        'SURFACE_OPERATOR_CORRECTION_READINESS':'READY_FOR_BAQP_SPECIFIC_DESIGN','BOUNDARY_DATA_CONVERSION':'UNCHANGED_VERIFIED',
        'C5_C7_ACTIVATION':'NOT_RUN','CLOT_GEOMETRY':'NOT_RUN','FRACTURE':'NOT_RUN'}
    write('FINAL_STATUS.json',statuses)
    sections=[]
    def add(title,text):sections.append(f'## {len(sections)+1}. {title}\n\n{text.strip()}\n')
    add('Why boundary-force conversion is no longer the main blocker',r'''
已有轴向和剪切非仿射制造基准的边界数据局部化验证保持有效，范围仍是这两个完整态基准。当前问题出现在**精确制造目标节点坐标本身**，无需新的平衡求解即可复现，因此本轮检查的是运动学算子。旧边界转换、载荷、残差和 C5/C7 接受态均未改。

读取来源：[剪切审计](../nonlocal_boundary_nonaffine_shear_localization_v1/REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md)、[轴向审计](../nonlocal_boundary_nonaffine_localization_v1/REVIEW_NONAFFINE_BOUNDARY_LOCALIZATION_V1.md)、[本轮冻结配置](PREDECLARED_CONFIG.json)、[目标态复现](SAVED_TARGET_REPRODUCTION.json)。两个保存的 X、V、目标坐标和解析 F 均逐字节匹配；用原 `block_family` 恢复相同拓扑，仅各调用一次固定态 `baqp.evaluate()`，不求平衡。

N19=6,859 个物理粒子、690,020 个无向键，L=1.125 mm，h=62.5 µm，δ=251.25 µm，δ/h=4.02；quadratic 核。SVK μ=1000 Pa、λ=0、bulk=666.6666666666666 Pa。PHASE=DISABLED；equilibrium、fracture、clot、C5/C7 均未运行。
''')
    observed=[]
    for f,v in closures.items():
        for channel in ['Fbar','Ftilde','Psi_Fbar','Psi_Ftilde']:
            z=v[channel];observed.append([f,channel,z['actual_max'],z['actual_max_scaled'],z['predicted_max_scaled']])
    add('Observed Fbar/Ftilde/Psi_plus bias',r'''
F 偏差使用 Frobenius 范数并除以固定 0.004；这不是除以局部 F 的相对误差。Ψ⁺ 列为最大绝对能量密度偏差，单位 Pa，并除以各场同一个解析节点最大值。剪切 Ψ_ref=0.004016032032 Pa，轴向 Ψ_ref=0.016064064 Pa；节点和键共用各自场的参考。
'''+'\n'+table(['Field','Channel','Actual max (F dimensionless / Ψ Pa)','Actual / scale','Predicted / scale'],observed)+r'''
剪切近边界 Fbar 约 2.976%、Ftilde 约 2.943%、Ψ⁺ 约 6.43%，描述的是目标态的重构偏差，不能与上一轮平衡求解的目标恢复精度混合。[全量因果摘要](CAUSAL_IDENTITY_SUMMARY.json)。
''')
    add('Quadratic BAQP reconstruction identity',r'''
对任意二次映射 y，Taylor 展开没有余项：

$$\eta_a=F^{exact}_{ac}(X_i)\xi_c+\tfrac12H_{acd}\xi_c\xi_d.$$

原 BAQP 定义 $A_{ab}=\sum_j\omega\eta_a\xi_bV_j$，故

$$A_{ab}=F^{exact}_{ac}K_{cb}+\tfrac12H_{acd}M_{3,cdb},\qquad
\boxed{\Delta\bar F_{an}=\tfrac12H_{acd}M_{3,cdb}(K^{-1})_{bn}}.$$

这是由当前 [baqp.py](../../pd_clot/pfpd_candidate/baqp.py) 的 $\bar F=AK^{-1}$ 推导的有限和恒等式，条件是 K 可逆；不是从 BB/OSB 论文移植的公式。实际实现的位移形式 $I+A_uK^{-1}$ 在精确算术下等价。

剪切 $u_x=\gamma Z+\beta_s Z^2/L$，$\gamma=0.003,\beta_s=0.001$；轴向 $u_z=\epsilon Z+\beta Z^2/L$，$\epsilon=0.003,\beta=0.001$。令 $a=x$ 或 $z$：

$$H_{azz}=2\beta/L,\qquad \Delta\bar F_{an}=(\beta/L)M_{3,zzb}(K^{-1})_{bn}.$$

其余输出行只有浮点误差。通用 Hessian 收缩与此专门形式分别计算并核对。
''')
    add('Second- and third-moment definitions',r'''
$$K_{cb}=\sum_{j\in H_i}\omega_{ij}\xi_c\xi_bV_j,\qquad
M_{3,cdb}=\sum_{j\in H_i}\omega_{ij}\xi_c\xi_d\xi_bV_j.$$

ω 无量纲、ξ 单位 m、V 单位 m³，所以 K、M3、G=M3_zz: K⁻¹ 单位分别为 m⁵、m⁶、m；H 单位 m⁻¹。独立累加器按标量分量与端点 scatter 求和，反向键的 K 项同号、M3 项异号；不借用 evaluate 的 F 输出生成预测。

[MOMENT_TENSORS.npz](MOMENT_TENSORS.npz) 保存全部物理节点的 K、invK、M3、均匀体积对照和 G，并包含最近面、棱、角距离，family size、ω0、condition(K) 和几何类。分类只数“距某轴最近面小于 δ”的轴数，0/1/2/3 分别为 bulk-like/face-near/edge-near/corner-near；它不宣称多跳力模板完全闭合。
'''+'\n'+table(['Geometry','Nodes','Family min','Family max','max cond(K)','max ||M3|| [m6]','max predicted / 0.004'],summaries))
    add('Continuous full-ball reference',r'''
对当前 $\omega(r)=(1-r/\delta)^2$，$r<\delta$：

$$K_{full}=\frac{4\pi}{3}I\int_0^\delta r^4(1-r/\delta)^2\,dr
=\frac{4\pi\delta^5}{315}I,\qquad \omega_{0,full}=\frac{2\pi\delta^3}{15}.$$

球域和径向核在 $\xi\mapsto-\xi$ 下不变，而每个三阶单项式变号，因此所有 M3 分量积分为零，包括 zzx、zzy、zzz。二次 Fbar 重构由第 3 节恒等式精确成立。
'''+f'\nK 标量系数 **{fmt(full["K_analytic"][0][0])} m⁵**；解析 ||M3||=0。独立 8×12×24 极坐标 Gauss/方位积分的 K 相对误差 {fmt(full["K_relative_error"])}，M3 数值范数 {fmt(full["M3_numerical_norm"])} m⁶，按 ||K||δ 缩放为 {fmt(full["M3_numerical_scaled"])}。\n\n[连续球参考](CONTINUOUS_FULL_BALL_REFERENCE.json)。')
    add('Discrete full-lattice reference',f'''
原有 inversion-symmetric offset 集含 {lat['neighbors']} 个邻居，所有体积为 h³，使用同一核和 δ。K 相对连续球参考误差为 **{fmt(lat['K_relative_error_vs_continuum'])}**；||M3||={fmt(lat['M3_norm'])} m⁶（按 ||K||δ 缩放 {fmt(lat['M3_scaled'])}）。

在全部 6,859 个物理观测位置外延同一多项式，以完整 offset 家族直接重构，shear/axial 最大 Fbar 误差/0.004 分别为 **{fmt(lat['full_reconstruction']['shear']['max_scaled'])} / {fmt(lat['full_reconstruction']['axial']['max_scaled'])}**，低于 1e-10。

K 的 0.3724% 积分误差不妨碍此二次场重构，因为同一个离散 K 在 A 与 K⁻¹ 中配对，奇矩才决定二次误差。H3/H4 作为本轮偏差来源不获支持；这不证明 BAQP 对所有非仿射场都精确。

[完整晶格参考](FULL_LATTICE_REFERENCE.json)、[全观测点直接重构](FULL_LATTICE_RECONSTRUCTION.npz)。
''')
    add('Truncated-support reference',r'''
R1 为完整对称晶格；R2 保留实际物理邻居坐标，诊断体积统一设 h³；R3 为实际梯形体积。R2/R3 不改变已有生产状态，只建立独立诊断有限和。
'''+'\n'+table(['Variant','max error / 0.004','V-weighted RMS / 0.004'],[[k,v['max_scaled'],v['weighted_RMS_scaled']] for k,v in source['variants'].items()])+r'''
R2 单独已产生约 2.85% 最大偏差，R3 约 2.98%；支持截断足以产生主体偏差。连续与晶格 full references 均为零或数值零，符合用户 CASE A 的判据。[逐节点三变体表](FAMILY_VARIANT_COMPARISON.csv)。
''')
    add('Effect of physical boundary volumes',r'''
$$A=E_{support}-E_{full},\quad B=E_{physical}-E_{support},\quad C=E_{physical}-E_{full},\quad A+B=C.$$

这些是逐节点向量贡献。各自的最大范数和 RMS **不具有标量加法关系**，不能用两列最大值相减定义体积贡献。
'''+'\n'+table(['Contribution','max / 0.004','V-weighted RMS / 0.004'],[[k,source[k]['max_scaled'],source[k]['weighted_RMS_scaled']] for k in ['support_contribution','volume_contribution','total_contribution']])+f'''
向量闭合残差/0.004={fmt(source['decomposition_residual_scaled'])}。体积贡献最大约 0.583%、RMS 约 0.275%，可测且次于支持贡献，不能忽略。距边界 4h 的 {geometry_extra['depth_4h']['nodes']} 个节点均有 256 邻居，R2 最大 {fmt(geometry_extra['depth_4h']['R2_max_scaled'])}，R3 仍为 {fmt(geometry_extra['depth_4h']['R3_max_scaled'])}（均 /0.004）：边界半体积可破坏完整坐标集合的加权对称性。

已检查相对两端面、侧面、棱、角的有符号 M3_zz: 和 G，见 [GEOMETRY_CAUSAL_DETAILS.json](GEOMETRY_CAUSAL_DETAILS.json)。Z 轴中心线上上下端偏差反号；侧面可产生输入 x/y 分量，不能把所有误差都看成 zz 分量或赋予统一符号。[向量分解](SOURCE_CONTRIBUTIONS.npz)、[摘要](SOURCE_DECOMPOSITION.json)。
''')
    add('Exact Fbar error reconstruction',table(['Field','max residual / 0.004','V-weighted RMS / 0.004'],[[f,v['Fbar']['residual_max_scaled'],v['Fbar']['residual_RMS_scaled']] for f,v in closures.items()])+r'''
全部节点通过预声明 max≤1e-10、weighted RMS≤1e-11；未改门槛。实际残差是完整张量之差的范数，不是两个误差范数之差。通用/专门预测之间的误差也接近舍入尺度。独立 K 与原 evaluate 的 K 对照及保存目标态 bitwise 检查通过。

[逐节点闭合](FBAR_MOMENT_CLOSURE.csv)、[深度统计](FBAR_DEPTH_PROFILES.csv)、[因果摘要](CAUSAL_IDENTITY_SUMMARY.json)。
''')
    add('Role of K anisotropy',r'''
保留物理 M3，分别使用物理 K、连续完整球 K 及完整晶格 K：

$$R_K=\frac{\|\tfrac12HM_3K_{physical}^{-1}\|}{\|\tfrac12HM_3K_{full}^{-1}\|}.$$
'''+f'''
连续球参考的 R_K 范围 **{source['K_modifier_continuous_range']}**；晶格参考 **{source['K_modifier_lattice_range']}**；有效节点 {source['modifier_valid_count']}。安全条件为分母/0.004>1e-10，否则存 NaN，不能将 0/0 当放大率。condition(K) 范围 {source['condition_K_range']}。

这个按用户定义命名的 “K anisotropy modifier” 同时包含 K 幅值损失与各向异性，不是纯各向异性量。它会放大或转向已有的奇矩误差；对任意可逆 K，若 M3=0，二次误差仍为零，因此 K 本身不是本轮根因。
''')
    add('Exact Ftilde projection identity',r'''
二次 y 的 secant 等于中点解析梯度乘 ξ，且解析 F 沿位置线性：

$$\eta=F^{exact}_{mid}\xi,\qquad F^{exact}_{mid}=\tfrac12(F^{exact}_i+F^{exact}_j).$$

代入原 $\widetilde F=F_{av}+(\eta-F_{av}\xi)\otimes\xi/|\xi|^2$，令 $D=(\Delta\bar F_i+\Delta\bar F_j)/2$、$n=\xi/|\xi|$，得到

$$\boxed{\Delta\widetilde F=D(I-n\otimes n)}.$$

先用实际节点 ΔF 独立核验投影恒等式；再用**独立矩预测**替换两个节点 ΔF，验证从 M3 一直到键的完整因果链。
'''+'\n'+table(['Field','Prediction','max residual / 0.004','Weighted RMS / 0.004'],[[f,k,v[k]['residual_max_scaled'],v[k]['residual_RMS_scaled']] for f,v in closures.items() for k in ['Ftilde','Ftilde_moment_chain']])+r'''
全部 690,020 个键满足相同 1e-10/1e-11 门槛；RMS 权重为原 BAQP 能量权重 $2\bar\omega_{ij}V_iV_j$。结论是节点重构偏差经横向投影保留，不是本轮证据所支持的额外键修正缺陷。
''')
    with (OUT/'BOND_ORIENTATION_PROFILES.csv').open() as f:o=[r for r in csv.DictReader(f) if r['field']=='shear']
    add('Bond-orientation dependence',r'''
逐键保存 |nx|、|ny|、|nz|、r/δ、中点最近面深度/δ、实际和预测 ΔFtilde 范数、D 范数及削减比。投影的 Frobenius 范数满足 $\|D(I-nn^T)\|_F\le\|D\|_F$；仅在 D 安全非零时报告比值。
'''+'\n'+table(['|nz| bin','Bonds','Actual max / 0.004','Weighted RMS / 0.004','Weighted mean reduction ratio'],[[f"[{v['lower']},{v['upper']}{']' if v['bin']=='4' else ')'}",int(v['bonds']),float(v['actual_F_max_scaled']),float(v['actual_F_RMS_scaled']),float(v['reduction_ratio_weighted_mean'])] for v in o])+r'''
接近横向的键组保留较大偏差；但削减比不是单纯的 $\sqrt{1-n_z^2}$，因为实际 G 在棱角含 x/y 分量且 D 随空间变化。比值均值也不必按 |nz| 单调。接近 z 的键仍可能保留横向输入分量。[键方向表](BOND_ORIENTATION_PROFILES.csv)、[逐键几何](BOND_GEOMETRY.npz)。
''')
    add('Axial/shear geometry cross-check',r'''
$$G_{i,n}=M_{3,i,zzb}(K_i^{-1})_{bn},\qquad
\Delta\bar F^{shear}_{xn}=(\beta_s/L_s)G_{i,n},\quad
\Delta\bar F^{axial}_{zn}=(\beta/L_a)G_{i,n}.$$
'''+f'''
映射输出行后的实际差最大/0.004={fmt(cross['max_mapped_difference_scaled'])}，体积加权 RMS={fmt(cross['weighted_RMS_difference_scaled'])}；均低于 1e-10。矩预测映射差仅 {fmt(cross['predicted_row_difference_scaled'])}。保留历史输入 L_axial=0.001125、L_shear=18h=0.0011250000000000001 的一 ulp 差异，未为求完全重合修改输入。

几何与核相同，G 相同；这个交叉结果将根因定位到运动学算子，不能归因于剪切本构或牵引专属现象。[逐节点映射](AXIAL_SHEAR_GEOMETRY_CROSSCHECK.csv)。
''')
    add('Psi_plus causal reconstruction',r'''
沿用 [svk_spectral_split.py](../../pd_clot/pfpd_candidate/svk_spectral_split.py)：

$$E=\tfrac12(F^TF-I),\quad E_+=\sum_\alpha\max(e_\alpha,0)q_\alpha\otimes q_\alpha,$$
$$\Psi^+(F)=\tfrac\lambda2\langle\operatorname{tr}E\rangle_+^2+\mu\,E_+:E_+.$$

当前候选使用 Green 应变/SVK 谱分裂；公式来源是已冻结实现及其原谱分裂审计，本轮不重新宣称它等同论文的完整 PFPD 模型。节点取 $F_{pred}=F_{exact}+\Delta\bar F_{M3}$，键取 $F_{pred}=F_{exact,mid}+\Delta\widetilde F_{M3}$，完整非线性调用同一谱分裂。预测不是对实际 ΔΨ 拟合，也不采用线性灵敏度作为证据。
'''+'\n'+table(['Field','Sample','Actual max bias [Pa]','Predicted max bias [Pa]','max causal residual / Psi_ref'],[[f,k,v[k]['actual_max'],v[k]['predicted_max'],v[k]['residual_max_scaled']] for f,v in closures.items() for k in ['Psi_Fbar','Psi_Ftilde']])+r'''
四组因果残差均远低于 1e-8；因而本轮约 6.4% 的剪切正能量偏差可由重构 F 解释，没有独立谱映射缺陷的证据。这是这些状态上的因果核验，不是对所有谱分裂分支或一般材料的证明。[深度数据](PSI_BIAS_RECONSTRUCTION.csv)。
''')
    erows=[]
    for f,e in energy.items():
        for c in ['intact','positive']:
            v=e[c];erows.append([f,c]+[v[k] for k in ['U0_J','U1_J','U2_J','kinematic_bias_J','weight_quadrature_bias_J','total_bias_J']])
    add('Global energy source decomposition',r'''
以无向键计数，保留原能量权重 $b_{ij}=2\bar\omega_{ij}V_iV_j$：

$$U_0=\sum_{i<j}b_{ij}W(\widetilde F_{ij}),\quad
U_1=\sum_{i<j}b_{ij}W(F^{exact}_{mid}),\quad
U_2=\int_\Omega W(F^{exact}(X))\,dV.$$

$$\Delta U_{kin}=U_0-U_1,\quad\Delta U_{weight/quad}=U_1-U_2,\quad
U_0-U_2=\Delta U_{kin}+\Delta U_{weight/quad}.$$

U1 只在本构输入中替换 F，未创建新生产模型。正能量类比将 W 换成 Ψ⁺，没有相场演化。下表全为 J。
'''+'\n'+table(['Field','Channel','U0 [J]','U1 [J]','U2 [J]','U0-U1 [J]','U1-U2 [J]','U0-U2 [J]'],erows)+f'''
剪切 intact 的运动学、权重/积分及总相对偏差分别为 {100*energy['shear']['intact']['kinematic_relative_to_U2']:.8f}%、{100*energy['shear']['intact']['weight_quadrature_relative_to_U2']:.8f}%、{100*energy['shear']['intact']['total_relative_to_U2']:.8f}%。全局负号只表示净能量差，并不意味着各节点正能量偏差都为负。

32 点 Gauss 对照 64 点的最大相对差 {fmt(max(max(e['Gauss32_vs64_relative']) for e in energy.values()))}；四组分解闭合均为 0 J，重复求和确定性检查通过。深度壳层按键中点半 h 分组，连续参照用嵌套中心立方体能量差，详见 [ENERGY_DEPTH_PROFILES.csv](ENERGY_DEPTH_PROFILES.csv)。不能把壳层积分或整体能量吻合等同于局部无偏。
''')
    add('Primary root cause',r'''
在这两个固定二次场上，完整连续球与完整晶格的奇矩为零；截断均匀体积家族产生显著奇矩；物理体积改变其幅度和分布；独立矩恒等式、投影恒等式、跨场行映射和非线性能量重建全部通过。因此：

**HORIZON TRUNCATION → NONZERO ODD MOMENTS → FBAR BIAS → PROJECTED FTILDE BIAS → PSI_PLUS BIAS**。

H1 是主体来源；H2 为可测的次级贡献。H3 的“全晶格仍有显著二次重构误差”及 H4 的“完整 BAQP 对此二次场天然不精确”不成立于本基准；H5/H6 的额外键或谱映射误差解释不获支持；H7 无需作为此次结论。这个结论不扩展为所有非仿射场、任意边界或非规则点云的普遍误差定理。

完整状态标签保存为 [FINAL_STATUS.json](FINAL_STATUS.json)：
'''+'\n'+table(['Status','Value'],[[k,json.dumps(v,ensure_ascii=False) if isinstance(v,dict) else v] for k,v in statuses.items()]))
    add('What is and is not justified about surface correction',r'''
证据足以进入 BAQP 专属算子修正的设计，但尚未选取或验证任何修正。

Shojaei 等（2026）讨论截断非局部算子向完整 horizon PD 算子的匹配；其 Eq.33/40 分别展示完整/半域奇矩行为。这为本轮诊断提供动机。其 BB/OSB 算子与权重目标具有模型依赖性，不能直接作为 BAQP 权重。来源：[原始出版 PDF](https://tore.tuhh.de/bitstreams/8fda7c78-cb40-401e-b8ce-206cd986fafb/download)，[DOI](https://doi.org/10.1002/nme.70385)，[缓存哈希及阅读位置](LITERATURE_PROVENANCE.json)。本轮的三维 M3/K 及 Ftilde 恒等式由本项目 BAQP 定义直接推导。

公开原始全文的既有缓存经 SHA256 复核，本轮网页打开受工具 10 MiB 限制；未重新下载 50 MB 文件，也未声称与用户上传 PDF 逐字节核对。没有应用 perifit、BB/OSB 权重、ghost/mirror、修改 K/M3 或改变物理 horizon。
''')
    add('BAQP-specific future correction requirements',r'''
设计义务完整列在 [BAQP_OPERATOR_CORRECTION_TARGETS.md](BAQP_OPERATOR_CORRECTION_TARGETS.md)：仿射精确、K 可逆与条件性、适当的二次奇矩目标、互易性、共同能量一阶变分、保持物理 horizon、区分边界数据转换、非规则点云可用性、不针对 C5/C7 首激活调参。

还需正视一个直接代数约束：平直单侧边界上，若所有权重非负且法向 ξ 全同号，则法向三阶矩不可能在法向二阶矩为正时变为零。仅靠现有家族内正标量权重强制“所有奇矩为零”可能不可行。此处是 BAQP 设计约束推导，不能误引为论文结论；见设计说明中的条件和证明。未实现任何替代构造。
''')
    add('Implication for future clot-wall damage initiation',r'''
当前候选的损伤驱动力依赖 Ψ⁺，因此真实 clot-wall 附近的首激活可能同时受真实应变和截断重构偏差影响；这是基于本轮算子证据的后续风险推断。不能据此断言现存 C5/C7 热点已完全解释、首激活应变已修正、mesh objectivity 已恢复，或真实血栓将在哪里断裂。

本轮没有 clot 几何、相场、断裂、CFD 或 activation 计算。未来须在完成算子设计并单独验证后，再区分表面误差与物理驱动；旧 C5/C7 科学状态保持原样。
''')
    add('One next task only',f'''
唯一建议下一任务：**BAQP_SURFACE_OPERATOR_CORRECTION_DESIGN_V1**。先研究 BAQP 多项式一致性目标、单侧正权重可行性及统一能量—力约束，形成可审阅设计。本轮不执行该任务。

10 项指定定向测试全部通过（0.67 s）；未运行全量 pytest。主诊断使用最多 8 CPU workers，BLAS=1，墙钟 {resource['wall_seconds']:.3f} s；采样峰值进程树 RSS={resource['peak_process_tree_RSS_bytes']} B，主进程高水位={resource['peak_main_RSS_bytes']} B，较大者 {max(resource['peak_process_tree_RSS_bytes'],resource['peak_main_RSS_bytes'])/2**30:.3f} GiB。最小可用系统 RAM 比例 {100*resource['minimum_available_RAM_fraction']:.2f}%（要求 ≥20%）；GPU_NOT_USED_BY_DESIGN。资源值对应科学诊断阶段，绘图和报告不包含在这段墙钟中。

425 个冻结输入/旧源码/历史结果的 SHA256 在计算前后保持不变，最终检查另见 FINAL_ARTIFACT_QA.json；新模块仅作诊断，不接求解器。旧 manifest、旧 assertions、kernel、volumes、材料、谱分裂、能量、边界转换均未改。分支 `{git['branch']}`，HEAD `{git['HEAD']}`；不提交、不推送。

交付索引：[数据字典](DATA_DICTIONARY.md)、[六组图](../../results/baqp_near_boundary_kinematic_energy_audit_v1/README.md)、[定向测试](TARGETED_TESTS.json)、[历史完整性](HISTORICAL_INTEGRITY.json)、[资源记录](RESOURCE_USAGE.json)。所有图为英文 600 dpi PNG 与矢量 PDF；散点抽样规则见 FIGURE_SAMPLING.json，全部 gate 使用全量数据。

**停止并等待人工审核；不自动开始修正或后续科学计算。**
''')
    (OUT/REPORT).write_text('# BAQP near-boundary kinematic / energy audit V1\n\n2026-10-06 · 固定态诊断 · 原物理公式未改\n\n'+ '\n'.join(sections))
    design=r'''# BAQP operator correction targets — design only

状态：READY_FOR_BAQP_SPECIFIC_DESIGN。没有选择或实现算子修正。

本轮证据及适用范围见 [REVIEW](REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md)。以当前 quadratic kernel、固定 h/δ、两种 Z-Z 二次场为证据；任意 Hessian、非规则点云及曲面边界尚须另行设计验证。

## Target algebra and invariants

当前 $\Delta\bar F=\tfrac12HM_3K^{-1}$。保持当前线性重构形式时，对任意二次 H 消除误差需要合适的三阶矩目标；仅将 K 调成连续 K 不会消除 M3。二阶一致性与可逆性需同时满足，不能通过削弱法向支持“制造”零误差。

| Obligation | Future design acceptance evidence |
|---|---|
| Affine exactness | 对任意常数 F 验证重构精确，不限于一次轴向拉伸 |
| Shape tensor invertibility | 保持 K 满秩，分析条件数与接近退化点云的行为 |
| Quadratic / odd-moment consistency | 明确目标是完整有限 horizon 算子，推导一般 Hessian 的误差；区分逐行 G 与完整 M3 条件 |
| Pair reciprocity | 端点重编号及键反向不改变标量总能量与一致的成对力分配；不能假设节点权重各自拟合自然互易 |
| One energy and first variation | 修正后的依赖关系须进入同一总能量，完整变分得到所有力项；若权重依赖构形须包含其导数 |
| Physical horizon | 不通过扩大/改变 δ 或 kernel 来掩盖原问题；未来设计必须明确支持和物理长度 |
| Boundary-data conversion | 载荷转换与表面算子分别定义；若改变算子，重新验证转换与该能量的匹配，不直接沿用旧数值结果 |
| Irregular clouds | 验证体积、邻域各向异性、曲率和不规则布点的可解性；正则化不得未经说明改变物理 |
| No activation fitting | 不以 C5/C7 首激活位置/应变或旧热点为权重优化目标 |
| Intact regression and energy quadrature | 仿射、二次与运动学—权重能量分解分别验收；U1 不是替代模型 |

## Feasibility constraint for nonnegative scalar weights

设平直边界法向 n 的保留键都满足 $\xi_n\ge0$，有效标量权重 $a_j=\omega_j V_j c_j\ge0$。若 $K_{nn}=\sum a_j\xi_n^2>0$，则至少一项满足 $a_j>0,\xi_n>0$，因此

$$M_{3,nnn}=\sum_j a_j\xi_n^3>0.$$

反向单侧家族同理为负。若强制该三阶矩为零，则在这些条件下只能令所有非零法向键权重为零，导致法向 K 退化。因此“现有单侧家族 + 非负标量权重 + K 可逆 + 所有二次奇矩消失”不能同时满足。对于当前纯法向二阶导数，M3_nn: K⁻¹=0 等价于 M3_nn:=0（K 可逆），也受此约束。

这是本轮从有限和直接推出的必要性结果，未声称 Shojaei 论文给出此 BAQP 定理。未来可比较更丰富的重构空间或显式高阶项等数学途径的可行性，但不得由此默认采用有符号权重、ghost、扩大 horizon 或任何未审核方案。对称完整 lattice 外延仅是本轮参考，不是 ghost 实现。

## Relation to Shojaei PD-to-PD matching

Shojaei 等的完整 horizon 目标匹配思想提供设计方向，BB/OSB 权重与距离幂次服务于各自算子。BAQP 的目标需从 K、M3、Ftilde、能量及其一阶变分重新推导。原始来源和具体位置见 [LITERATURE_PROVENANCE.json](LITERATURE_PROVENANCE.json) 及 [原论文](https://doi.org/10.1002/nme.70385)。未调用 perifit，没有验证其权重适用于 BAQP。

## One next task

BAQP_SURFACE_OPERATOR_CORRECTION_DESIGN_V1：给出满足上述义务且数学可行的 BAQP 专属候选设计和独立验收计划，交人工审核。当前不执行，不启动 activation、fracture、clot-wall 或 CFD。
'''
    (OUT/'BAQP_OPERATOR_CORRECTION_TARGETS.md').write_text(design)
    dictionary=r'''# Data dictionary and reproduction

本目录是固定态诊断。科学入口 `scripts/near_boundary_kinematic_energy_audit.py`，独立诊断模块 `pd_clot/pfpd_candidate/nonlocal_boundary/kinematic_energy_diagnostics.py`；绘图和报告入口分别为 `scripts/plot_near_boundary_kinematic_energy.py` 与 `scripts/report_near_boundary_kinematic_energy.py`，后两者只读取已存科学数组。

运行环境：`/home/lzy/projects/computation_examples/_shared/python_environment/bin/python`；项目目录 `PYTHONPATH=.`，`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`，`python -B`。重跑科学入口会覆盖本轮输出，需先使用新输出目录保存独立复现；本轮没有重跑旧 equilibrium。

## Index and units

节点索引沿用保存的 X（6,859 行），键索引与 BOND_GEOMETRY.pairs（690,020 行）完全对齐。数组均未为画图重排。MOMENT_TENSORS 的 K/invK/M3 单位 m⁵/m⁻⁵/m⁶，G 单位 m，X/depth distances 单位 m，volume/omega0/bond_weight 单位 m³。F 及其偏差无量纲，Ψ 单位 Pa，逐键能量 U 单位 J；方向余弦、condition(K)、比值无量纲。

MOMENT_TENSORS 包含物理 R3、uniform R2、lattice R1 及 continuous 参考；节点 invK 直接保存，uniform invK 也保存。SOURCE_CONTRIBUTIONS 是单个非零输出行的三个输入方向分量，support+volume=total；用于剪切，轴向保持其历史一 ulp L 差异。geometry_class 数轴方向近面数量，不用于力学模板闭合判定。

NODE_DIAGNOSTICS 保存完整实际/预测/残差 3×3 张量；标量 F_error 为 Frobenius 范数。BOND_DIAGNOSTICS 的 predicted_F_error 来自实际节点 D 投影；moment_predicted_F_error 来自独立 M3 预测，是 Ψ 因果重建的输入。projection_residual / moment_chain_residual 由完整张量之差计算再取范数。为节省空间未保存重复的全键 Ftilde 张量。

Psi_actual_bias=Psi_actual−Psi_exact；Psi_predicted_bias=Psi_predicted−Psi_exact，均保留符号。Psi_causal_residual=Psi_actual−Psi_predicted。CSV 的 actual_max/predicted_max/residual_max 取绝对值最大；节点/键 RMS 分别用 V_i / 2*wbar*Vi*Vj 权重。

所有 F 的 scaled 字段除以 0.004；Ψ scaled 字段除以各场 Psi_ref_Pa（解析物理节点最大，节点和键共享）。不使用局部趋零分母。K modifier 的安全条件为分母/0.004>1e-10；projection reduction ratio 的条件为 ||D||/0.004>1e-10。未定义值保留 NaN，并有有效计数/标记。

FBAR_MOMENT_CLOSURE.csv 每节点每场一行；FTILDE_PROJECTION_CLOSURE.csv 按深度聚合，完整逐键量在相应 NPZ。PSI_BIAS_RECONSTRUCTION.csv 按场、节点/键、深度聚合。ENERGY_DEPTH_PROFILES.csv 的每层连续 U2 来自内外两个中心立方体积分差；壳界=半 h 中点±h/4，裁剪到 [0,L/2]。全部壳层加总核对全局能量。最大深度半 h 壳可能没有键，需保留完整连续体积分（本次诊断已核对总和）。

## Provenance and figure policy

PREDECLARED_CONFIG 保存全部固定门槛、资源估算及几何/体积规则。PROTECTED_INPUT_HASHES 保护 425 个旧文件，HISTORICAL_INTEGRITY 和最终 QA 验证不变；旧 manifest 未改。SAVED_TARGET_REPRODUCTION 记录 bitwise 目标态复现和 Ftilde SHA。TASK_REQUEST 是本轮原请求。LITERATURE_PROVENANCE 指向原文缓存，不复制大型 PDF。

六幅图见 results 同名目录。FIGURE_SAMPLING 对每个散点面板记录固定 stride，不随机、不修改值；所有 gate 和 profiles 用全量数据。FIGURE_LAYOUT_QA 检查可见文字边界；最终 QA 检查 600 dpi PNG 与 PDF 无栅格 image 对象。资源统计仅覆盖科学诊断，非完整对话或绘图用时。

## Saved array inventory

下表是实际 NPZ schema；完整 dtype/shape 另存 ARRAY_SCHEMA.json。CSV 列名存 CSV_SCHEMA.json。

'''
    arrays={};csvs={};arrayrows=[]
    for path in sorted(OUT.glob('*.npz')):
        with np.load(path,allow_pickle=False) as a:
            arrays[path.name]={k:{'shape':list(a[k].shape),'dtype':str(a[k].dtype)} for k in a.files}
            arrayrows.append([path.name,', '.join(a.files)])
    for path in sorted(OUT.glob('*.csv')):
        with path.open() as f:
            r=csv.reader(f);head=next(r);count=sum(1 for _ in r)
            csvs[path.name]={'columns':head,'data_rows':count}
    write('ARRAY_SCHEMA.json',arrays);write('CSV_SCHEMA.json',csvs)
    (OUT/'DATA_DICTIONARY.md').write_text(dictionary+table(['File','Arrays'],arrayrows))
    titles=[('third_moment_fbar_causality','Third-moment Fbar causality'),('full_vs_truncated_moment_comparison','Full vs truncated moments'),
        ('ftilde_projection_bias','Ftilde transverse projection'),('psi_plus_bias_causal_reconstruction','Positive-energy causal reconstruction'),
        ('energy_error_source_decomposition','Global energy source decomposition'),('axial_shear_geometric_bias_invariance','Axial/shear geometric invariance')]
    index='# BAQP near-boundary kinematic / energy audit — figures\n\nEnglish, 600 dpi PNG and vector PDF. All plotted values come from saved audit data; scatter sampling is deterministic and recorded.\n\n'
    for n,(stem,title) in enumerate(titles,1):index+=f'{n}. **{title}**: [PNG](figures/{stem}.png) · [PDF](figures/{stem}.pdf)\n\n'
    index+='[Scientific report](../../verification/baqp_near_boundary_kinematic_energy_audit_v1/'+REPORT+') · [Data dictionary](../../verification/baqp_near_boundary_kinematic_energy_audit_v1/DATA_DICTIONARY.md)\n'
    (RESULT/'README.md').write_text(index)
    print('REPORT, DESIGN TARGETS, STATUS, DATA DICTIONARY AND FIGURE INDEX WRITTEN')


if __name__=='__main__':main()
