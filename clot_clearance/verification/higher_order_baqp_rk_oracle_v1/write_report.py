"""Render the review and data dictionary from saved results; no scientific runs."""
import csv
import json
import re
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
RESULTS = ROOT / 'results/higher_order_baqp_rk_oracle_v1'
FIGURES = '../../results/higher_order_baqp_rk_oracle_v1/figures/'

def read(name):
    return json.loads((OUT / name).read_text())

def dump(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')

def rows(name):
    with (OUT / name).open() as f:
        return list(csv.DictReader(f))

def table(headers, records):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, r)) + ' |' for r in records])

def fmt(x):
    return f'{float(x):.6e}'

def figure(name):
    return f'[600 dpi PNG]({FIGURES}{name}.png) · [Vector PDF]({FIGURES}{name}.pdf)'

def main():
    s = read('C1_PD2_SUMMARY.json')
    c = read('CONDITIONING_SUMMARY.json')['PD2']
    p = read('POLYNOMIAL_SUMMARY.json')
    affine = read('PD2_AFFINE_REGRESSION.json')
    julia = read('PERIDYNAMICSJL_CROSSCHECK.json')
    cpu, jcpu = read('RESOURCE_USAGE.json'), read('JULIA_RESOURCE_USAGE.json')
    git = read('STARTING_GIT_STATE.json')
    identity = read('INPUT_IDENTITY.json')
    regression = read('C1_BASELINE_REGRESSION.json')
    assert read('PYTHON_AUDIT_STATUS.json')['PASS']
    assert julia['status'] == 'PASS'
    log = (OUT / 'TARGETED_TESTS.log').read_text()
    testmatch = re.search(r'(\d+) passed in ([0-9.]+)s', log)
    assert testmatch and int(testmatch[1]) == 12
    tests = dict(passed=12, failed=0, skipped=0, wall_seconds=float(testmatch[2]),
                 suite='verification/higher_order_baqp_rk_oracle_v1/tests/test_pd2_oracle.py',
                 full_project_suite_run=False, log='TARGETED_TESTS.log')
    dump('TARGETED_TESTS.json', tests)
    statuses = dict(
        C1_BASELINE_REGRESSION='PASS', PD2_MOMENT_MATRIX='FULL_RANK',
        PD2_CONDITIONING='WELL_CONDITIONED', PD2_POLYNOMIAL_REPRODUCTION='VERIFIED',
        PD2_AFFINE_REGRESSION='PASS', PD2_AXIAL_QUADRATIC='PASS',
        PD2_SHEAR_QUADRATIC='PASS', PD2_CROSS_TERM='PASS', RK2_DIAG_CROSS_TERM='FAIL',
        FULL_3D_QUADRATIC_BASIS_REQUIRED='YES', PD2_FTILDE_DIAGNOSTIC='PASS',
        PD2_PSI_PLUS_DIAGNOSTIC='PASS', PERIDYNAMICSJL_CROSSCHECK='PASS', GMLS_REQUIRED='NO',
        HIGHER_ORDER_KINEMATIC_CANDIDATE='VERIFIED', PRODUCTION_FORCE='NOT_IMPLEMENTED',
        PRODUCTION_BAQP='UNCHANGED', BOUNDARY_DATA_CONVERSION='UNCHANGED_VERIFIED',
        C5_C7_ACTIVATION='NOT_RUN', PHASE='NOT_RUN', FRACTURE='NOT_RUN', CLOT_GEOMETRY='NOT_RUN')
    dump('FINAL_STATUS.json', statuses)
    decision = f'''# GMLS trigger decision

GMLS_REQUIRED = **NO**.

| Predeclared trigger | Evidence | Triggered |
| --- | --- | --- |
| A: rank-deficient physical support | 0 / 6,859 nodes; minimum rank 9 | NO |
| B: severe conditioning causing failed reproduction | Maximum condition {c['condition_max']:.9f}; all nine reproduction gates pass | NO |
| C: published RK implementation and independent package disagree | Installed Peridynamics.jl full PD2 passes all three representative-node comparisons | NO |
| D: explicit request for an irregular-cloud alternative | This task requests the saved N19 regular cloud | NO |

No GMLS solver, moment regularization, horizon change, or support repair was implemented.
RK2_DIAG fails the XZ field because its polynomial space omits XZ; this is the intended
negative control, not a PD2 failure or a GMLS trigger.

Evidence: [conditioning](CONDITIONING_SUMMARY.json), [reproduction](POLYNOMIAL_SUMMARY.json),
[Julia](PERIDYNAMICSJL_CROSSCHECK.json), [configuration declared before computation](PD2_CONFIGURATION.json).

Only proposed next task: **HIGHER_ORDER_BAQP_ENERGY_FORCE_CANDIDATE_V1**.
It has not been started. Stop for human review.
'''
    (OUT / 'GMLS_TRIGGER_DECISION.md').write_text(decision)

    sections = []
    def section(title, text):
        sections.append(f'## {len(sections)+1}. {title}\n\n{text.strip()}\n')

    section('Why this round uses published higher-order reconstruction instead of inventing a surface correction', r'''
本轮回答已经发表的完整二次 RK 梯度能否消除截断邻域的运动学偏差。答案是：**在已测试的 N19 规则网格及指定多项式场上，可以；HIGHER_ORDER_KINEMATIC_CANDIDATE = VERIFIED。**
固定坐标、邻居、核、视界、物理体积与旧目标位置，只改变重构空间的阶数。没有开发高阶力、改变表面权重或增加虚拟节点。

数学来源是 [Hillman, Pasetto, Zhou, Generalized reproducing kernel peridynamics, §§4.1–4.2, Eqs. (63)–(70), (77)–(81)](https://www.mcharleshillman.com/uploads/7/1/8/9/71894909/2019_generalized_reproducing_kernel_peridynamics.pdf)
中的差分非局部导数与矩条件，以及 [Behzadinasab, Trask, Bazilevs, Part I: core methods, §2.2, Eqs. (3), (5), Example 1](https://arxiv.org/pdf/2004.11477)
中的三维完整二次梯度离散。缓存 PDF 与页码、SHA256 见 [SOURCE_PROVENANCE.json](references/SOURCE_PROVENANCE.json)。本文公式是这些梯度定义在当前物理体积离散与无量纲局部坐标下的实现。

[Part II: wave propagation and enforcement of stress boundary conditions](https://arxiv.org/abs/2004.11478) 仅作为未来自然边界/FSI 文献记录，本轮没有采用其力或边界算法。
''')
    section('Current C1 root cause', r'''
前轮证据见 [原审计报告](../baqp_near_boundary_kinematic_energy_audit_v1/REVIEW_BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_V1.md)。对二次位移，C1 节点误差含二阶导数与加权三阶矩的收缩：

$$
K_{bd}^{(i)}=\sum_j\omega_{ij}\xi_b\xi_dV_j,\qquad
M_{3,bcd}^{(i)}=\sum_j\omega_{ij}\xi_b\xi_c\xi_dV_j,
$$
$$
(\bar F_{C1}-F_{exact})_{ak}
=\tfrac12 u_{a,bc}M_{3,bcd}^{(i)}(K_i^{-1})_{dk}.
$$

完整对称球或对称离散邻域的奇次矩消失；物理边界截断后一般不消失。该节点偏差通过当前键投影进入 $\widetilde F$，继而影响 $\Psi^+$。
本轮 C1 仅提取旧 `baqp.py` 的节点代数；没有调用会同时计算力的 `evaluate()`。两个旧场的目标位置、解析梯度和 C1 节点梯度均与旧存档 **逐位一致**。键误差指标的重算差异仅为浮点运算顺序差异，详见 [C1_BASELINE_REGRESSION.json](C1_BASELINE_REGRESSION.json)。
''')
    section('Published PD2/RK mathematical formulation', r'''
定义 $\xi_{ij}=X_j-X_i$，$q=\xi_{ij}/\delta$，列向量

$$Q(q)=[q_x,q_y,q_z,q_x^2,q_y^2,q_z^2,q_xq_y,q_xq_z,q_yq_z]^T.$$

原论文 Part I 的 Example 1 使用同一九项空间。常数由 $p_j-p_i$ 自动消去，因此没有常数基项。对每个物理节点：

$$M_i=\sum_{j\in H_i}\omega_{ij}Q_{ij}Q_{ij}^TV_j,\qquad
S=\frac1\delta\begin{bmatrix}I_3\\0_{6\times3}\end{bmatrix},\qquad M_i C_i=S,$$
$$g_{ij}=\omega_{ij}V_j C_i^TQ_{ij},\qquad
(\nabla_hp)_i=\sum_j(p_j-p_i)g_{ij},$$
$$\bar F_i=I+\sum_j(u_j-u_i)\otimes g_{ij}.$$

直接线性求解 $M_i C_i=S$；没有显式矩阵逆、伪逆或正则化。上式的物理梯度权重已经包含 $V_j$。

**离散再现论证。** 若局部支持满秩，$\sum_jQ_{ij}g_{ij}^T=M_iC_i=S$。任何二次多项式差分可写作 $p_j-p_i=c_i^TQ_{ij}$，因而 $\sum_j(p_j-p_i)g_{ij}^T=c_i^TS$，恰为节点处的物理梯度。此代数结论不要求邻域对称；数值实现仍需逐节点检查满秩与浮点残差。

这并没有把原始几何 $M_3$ 改成零。它通过完整二次重构让有效导数权重满足二次矩约束，与先前关于单侧邻域正标量权重的限制不矛盾。
''')
    section('Implementation conventions and coordinate scaling', r'''
主代码位于 [higher_order_rk](../../pd_clot/pfpd_candidate/higher_order_rk/)，入口是 [run_audit.py](run_audit.py)。矩阵行是位移输出分量、列是参考坐标导数：$F_{ab}=\delta_{ab}+\partial u_a/\partial X_b$。

| 量 | 单位/约定 |
| --- | --- |
| $X,\xi,\delta$ | m |
| $Q,\omega$ | 无量纲 |
| $V_j,M_i$ | m³ |
| $S$ | m⁻¹ |
| $C_i$ | m⁻⁴ |
| $g_{ij}$ | m⁻¹，已包含邻点体积 |
| $F$ | 无量纲 |

因此不能再次乘邻点体积，也不能遗漏选择矩阵的 $1/\delta$。直接物理多项式测试及独立单位缩放测试覆盖这一约定。

N19 表示每轴 19 点，共 6,859 点、690,020 条无向键、1,380,040 个有向邻接。$h=62.5\,\mu m$，$\delta=251.25\,\mu m$，$\delta/h=4.02$，$L=18h=1.125\,mm$。保留原核 $\omega(r)=(1-r/\delta)^2$ 及内部/面/棱/角的梯形体积 $h^3,h^3/2,h^3/4,h^3/8$。几何、拓扑、核和体积的原数组哈希全部相符，见 [INPUT_IDENTITY.json](INPUT_IDENTITY.json)。

轴向与剪切场直接复用旧类。历史轴向 $L=0.001125$ 与剪切/几何 $18h=0.0011250000000000001$ 的一 ULP 表示差异也予以保留，以确保目标数组逐位相等。

所有门槛、归一化与条件数标志在计算前写入 [PD2_CONFIGURATION.json](PD2_CONFIGURATION.json)。节点 RMS 用 $V_i$；键 RMS 用原能量权重 $2\bar w_{ij}V_iV_j$。各自定义为 $\sqrt{\sum w e^2/\sum w}$。节点深度取最近物理面距离，键深度取键中点的最近面距离。
参考缓存只依赖几何、邻域、体积、核、视界及基，SHA256 见 [CACHE_PROVENANCE.json](CACHE_PROVENANCE.json)。未保存每个场的完整键张量副本，只保存标量键误差及可复算节点状态。
''')
    groups = rows('PD2_CONDITIONING_GROUPS.csv')[:5]
    section('PD2 moment-matrix rank / conditioning',
        '全部节点满秩。SVD 判秩门槛为 $9\\epsilon_{mach}\\sigma_{max}$；没有剔除边界节点。\n\n' +
        table(['几何类', '节点数', '最小秩', '最小 σ_min [m³]', 'κ median', 'κ P95', 'κ P99', 'κ max'],
              [[r['group'], r['nodes'], r['rank_min'], fmt(r['sigma_min'])] +
               [f'{float(r[k]):.6f}' for k in ['condition_median','condition_P95','condition_P99','condition_max']] for r in groups]) +
        '\n\n类别只按几何确定：三个坐标方向各自到最近对应面的距离中，少于 δ 的方向数为 0/1/2/3，分别记为 interior-like/face-near/edge-near/corner-near。它不是力模板完全闭合的判据。\n\n'
        '预先规定的报告标志：κ≤10⁶ 为 WELL_CONDITIONED，10⁶<κ≤10¹⁰ 为 ELEVATED，κ>10¹⁰ 为 SEVERE；这是开发描述，不是普适稳定性定理。绝对奇异值很小与 M 的 m³ 单位有关，不能单独据此宣称病态。\n\n'
        '逐节点全部九个奇异值见 [PD2_NODE_CONDITIONING.csv](PD2_NODE_CONDITIONING.csv)，逐深度分组见 [PD2_CONDITIONING_GROUPS.csv](PD2_CONDITIONING_GROUPS.csv)。\n\n' + figure('pd2_conditioning'))
    section('Nine-term polynomial reproduction',
        '在全部 6,859 个物理节点直接计算九个**物理坐标多项式**的解析梯度。线性残差除以 1，二次残差除以 L（m）；不除以局部可能为零的解析梯度。门槛为 max≤10⁻¹⁰、体积加权 RMS≤10⁻¹¹。\n\n' +
        table(['p', 'max normalized residual', 'V-weighted RMS', '结果'],
              [[name, fmt(v['max_scaled']), fmt(v['weighted_RMS_scaled']), 'PASS'] for name,v in p.items()]) +
        '\n\n[PD2_POLYNOMIAL_REPRODUCTION.csv](PD2_POLYNOMIAL_REPRODUCTION.csv) 保存全部 61,731 个节点×多项式记录，含物理残差分量。热图每格取该几何类全部节点最大值。\n\n' + figure('pd2_polynomial_reproduction'))
    section('Affine regression',
        '先于二次场评估，使用旧基准的轴向与剪切仿射梯度 0.003。以下 F 误差均为 Frobenius 范数除以 0.004；最大值门槛 10⁻¹⁰。\n\n' +
        table(['场', 'PD2 max', 'PD2 V-weighted RMS', '结果'],
              [[name, fmt(v['max_scaled']), fmt(v['weighted_RMS_scaled']), 'PASS'] for name,v in affine.items()]) +
        '\n\n[PD2_AFFINE_REGRESSION.json](PD2_AFFINE_REGRESSION.json)。')

    definitions = {
        'axial': r'复用 $u_z=0.003Z+0.001Z^2/L$，其余位移为零；$F_{zz}=1+0.003+0.002Z/L$，其余为单位矩阵对应项。',
        'shear': r'复用 $u_x=0.003Z+0.001Z^2/L$，其余位移为零；$F_{xz}=0.003+0.002Z/L$，其余为单位矩阵对应项。',
        'cross': r'唯一新增场是 $u_x=0.004XZ/L$；$F_{xx}=1+0.004Z/L$，$F_{xz}=0.004X/L$，其余为单位矩阵对应项。没有调节幅值。'}
    for field,title in [('axial','Quadratic axial result'),('shear','Quadratic shear result'),('cross','XZ cross-term discriminator')]:
        text = definitions[field] + '\n\n' + table(['重构', 'max F error / 0.004', 'V-weighted RMS / 0.004'],
            [[model,fmt(s[field][model]['Fbar']['max_scaled']),fmt(s[field][model]['Fbar']['weighted_RMS_scaled'])] for model in ['C1','PD2','RK2_DIAG']])
        text += '\n\n在精确目标位置直接评估，没有平衡求解。PD2 的 max≤10⁻¹⁰、RMS≤10⁻¹¹ 均通过。'
        if field == 'cross':
            text += '\n\nRK2_DIAG 只是同一实现选择前六个基项，省略 XY/XZ/YZ。它通过旧轴向与剪切场，却在 XZ 场出现约 6.055% 的梯度尺度误差，构成预期负对照：**FULL_3D_QUADRATIC_BASIS_REQUIRED = YES**。该控制不等同于 Peridynamics.jl 名称为 `:RK2` 的含常数七项空间。'
        section(title,text)

    gr = rows('NODE_GEOMETRY_STATISTICS.csv')
    records = []
    for field in ['axial','shear','cross']:
        for cls in ['interior-like','face-near','edge-near','corner-near']:
            models = {r['model']:r for r in gr if r['field']==field and r['geometry_class']==cls}
            records.append([field,cls,fmt(models['C1']['max_scaled']),fmt(models['PD2']['max_scaled'])])
    section('C1 vs PD2 near-boundary comparison',
        '下表的每个最大值涵盖该类全部节点，F 误差归一化为 0.004。完整类内 RMS 见 [NODE_GEOMETRY_STATISTICS.csv](NODE_GEOMETRY_STATISTICS.csv)。\n\n' +
        table(['场','几何类','C1 max','PD2 max'],records) +
        '\n\n两个旧场的 C1 最坏节点误差保持约 2.976%，而 PD2 的三个场均达到浮点精度量级。XZ 也通过，排除了只靠轴向平方项偶然通过测试的解释。分层曲线使用全部节点且没有人为误差下限；若存在严格零值，log 图省略零标记，原始值仍保留在 CSV。\n\n' + figure('c1_vs_pd2_boundary_bias'))
    section('Diagnostic bond-associated Ftilde result', r'''
仅在三个 PD2 节点场全部通过后，使用原架构的诊断投影：

$$F_{avg}=\tfrac12(\bar F_i+\bar F_j),\quad
\widetilde F=F_{avg}+
\frac{(\eta_{ij}-F_{avg}\xi_{ij})\otimes\xi_{ij}}{|\xi_{ij}|^2},\quad
\eta_{ij}=x_j-x_i.$$

对二次位移，$F_{exact}$ 是坐标的一次函数，所以 $F_{exact}(X_{mid})=(F_{exact}(X_i)+F_{exact}(X_j))/2$。沿键积分也给出 $\eta_{ij}=F_{exact}(X_{mid})\xi_{ij}$，投影修正应为零。这是当前三个场的精确性论证，不推广到任意非线性场。

''' + table(['场','C1 max / 0.004','PD2 max / 0.004','PD2 energy-weighted RMS / 0.004'],
        [[field,fmt(s[field]['C1']['Ftilde']['max_scaled']),fmt(s[field]['PD2']['Ftilde']['max_scaled']),fmt(s[field]['PD2']['Ftilde']['weighted_RMS_scaled'])] for field in s]) +
        '\n\n全部 690,020 条键参与；max≤10⁻¹⁰、RMS≤10⁻¹¹ 均通过。[PD2_FTILDE_DIAGNOSTIC.csv](PD2_FTILDE_DIAGNOSTIC.csv) 为逐层结果。没有由此计算应力散度、Z、反力或任何高阶力。')
    section('Psi_plus result', r'''
仅在所有 PD2 键梯度通过后调用原有谱分裂程序 [svk_spectral_split.py](../../pd_clot/pfpd_candidate/svk_spectral_split.py)：

$$E=\tfrac12(F^TF-I),\qquad E^+=\sum_a\max(e_a,0)n_a\otimes n_a,$$
$$\Psi^+(F)=\frac\lambda2\langle\operatorname{tr}E\rangle_+^2+\mu E^+:E^+.$$

材料保持 $\lambda=0$、$\mu=1000\,Pa$。同一谱程序分别作用于解析 F、C1 F 与 PD2 F，没有改变应变度量、特征值处理或能量定义。每场统一取 $\Psi_{ref}=\max_i\Psi^+(F_{exact}(X_i))$，节点与键共用该物理域参考尺度，避免逐点除以近零能量。

''' + table(['场','Ψ_ref [Pa]'],[[field,fmt(s[field]['Psi_ref_Pa'])] for field in s]) + '\n\n' +
        table(['场','采样','C1 max relative bias','PD2 max relative bias','PD2 weighted RMS'],
            [[field,sample,fmt(s[field]['C1'][key]['max_scaled']),fmt(s[field]['PD2'][key]['max_scaled']),fmt(s[field]['PD2'][key]['weighted_RMS_scaled'])]
             for field in s for sample,key in [('node','Psi_Fbar'),('bond','Psi_Ftilde')]]) +
        '\n\n节点采用体积权重，键采用原能量权重。全部 max≤10⁻⁸、RMS≤10⁻⁹ 通过。旧剪切 C1 偏差约 6.43%，新增 XZ 场约 14.29%，PD2 三场均消除至数值精度量级。\n\n' + figure('pd2_ftilde_and_energy_bias'))
    section('Peridynamics.jl cross-check',
        f'独立对照调用现有未修改包：Julia {julia["Julia_version"]}，Peridynamics.jl {julia["package_version"]}，commit `{julia["commit"]}`。路径为 `verification/pfpd_candidate_v1/external/Peridynamics.jl`。包的 115 个源文件/项目文件前后哈希一致。\n\n'
        '从已安装源码确认 `RKCMaterial(monomial=:PD2)`，直接调用 `rkc_weights!` 与 `rkc_defgrad!`。包基项顺序为 `[x,y,z,x²,xy,xz,y²,yz,z²]`，与 Python 是同一空间的排列；两个选择矩阵都选前三个线性项。使用相同二次核、δ、物理 Vj。包实现将 1/δ 放在核系数中，Python 放在选择矩阵中，物理权重等价。\n\n'
        '包现有 `invreg` 采用 SVD 倒数路径；epsilon/lambda/beta 都设为 0，本轮满秩支持上没有截断、正则化或伪逆修复。主 Python 路径仍只做线性求解。构造的 `VelocityVerlet` 对象仅用于选择储存布局，未初始化、积分或调用力程序；没有升级或安装依赖。\n\n'
        '四个代表节点的 Python 零基编号为 3429（内部）、3420（面）、171（棱）、0（角）。Julia 输出文件使用一基编号，比较时明确转换。\n\n' +
        table(['场','max ‖F_Python−F_Julia‖ / 0.004','门槛','结果'],
            [[field,fmt(julia['case_results'][field]['max_F_difference_scaled']),'1e-10','PASS'] for field in s]) +
        f'\n\n还逐邻居核对 551 个物理有向权重，最大 $\\|g_P-g_J\\|\\delta$ 为 {fmt(julia["max_gradient_weight_difference_times_delta"])}。F 的微小差异与线性解法及 `Δu` 运算顺序有关。四点对照是独立软件检查；全节点科学依据仍是九项再现与解析场验证。\n\n'
        '[原始张量、配置与资源](PERIDYNAMICSJL_CROSSCHECK.json)；[原包适配脚本](julia_pd2_oracle.jl)。\n\n' + figure('python_vs_peridynamicsjl_pd2'))
    section('Whether GMLS is required',
        '**GMLS_REQUIRED = NO。** 没有秩亏、严重病态导致再现失败、RK/Julia 分歧或本轮不规则点云替代要求。未实现 GMLS。完整触发条件与证据见 [GMLS_TRIGGER_DECISION.md](GMLS_TRIGGER_DECISION.md)。')
    section('What this proves',
        '**已发表的完整二次 RK 梯度，在本轮固定规则网格和指定场上，消除了先前截断邻域造成的非仿射 F、诊断 Ftilde 和 Ψ_plus 偏差。** 这是 CASE A，允许将其标记为经过验证的高阶运动学候选。\n\n' +
        table(['状态','结果'],statuses.items()) +
        '\n\nRK2_DIAG 的 FAIL 是缺失交叉基项的预期负对照；BOUNDARY_DATA_CONVERSION 的 UNCHANGED_VERIFIED 表示沿用先前验证状态，本轮没有重新集成边界力。\n\n'
        f'12 项指定测试全部通过（{tests["wall_seconds"]:.2f} s），没有执行完整项目测试套件。测试日志见 [TARGETED_TESTS.log](TARGETED_TESTS.log)。Python 使用 {cpu["workers"]} 个线程工作单元、BLAS=1，墙钟 {cpu["wall_seconds"]:.3f} s；Julia 单工作单元 {jcpu["wall_seconds"]:.3f} s。两计算阶段合计 {cpu["wall_seconds"]+jcpu["wall_seconds"]:.3f} s，最大监测进程树 RSS {max(cpu["peak_process_tree_RSS_bytes"],jcpu["peak_process_tree_RSS_bytes"])/2**30:.3f} GiB；这是计算阶段时间，不含查阅、绘图与写报告。两个阶段最小可用内存占比均超过 61%。GPU_NOT_USED_BY_DESIGN。\n\n'
        '输出字段、单位及重绘入口见 [DATA_DICTIONARY.md](DATA_DICTIONARY.md)。历史受保护文件校验见 [HISTORICAL_INTEGRITY.json](HISTORICAL_INTEGRITY.json)，最终交付核验见 [FINAL_ARTIFACT_QA.json](FINAL_ARTIFACT_QA.json)。')
    section('What this does NOT prove',
        '本轮没有实现高阶生产力，不能宣称高阶 BAQP 力学、能量—力一致性、稳定性、平衡、网格客观性、断裂或血栓壁面模型已经验证。也没有进行 C5/C7 激活、相场演化、真实血栓几何、CFD、Ogden、黏附或微泡加载。\n\n'
        '完整二次再现本身不保证力密度算子的一致性或稳定性；这一限制也见 [Part I](https://arxiv.org/pdf/2004.11477) 的研究动机。本轮的 Ftilde 是当前 BAQP 架构的运动学诊断，不等同于移植了论文中的完整高阶动力学模型。现有生产源码、旧结果和边界转换没有改变。')
    section('Implication for future irregular clot clouds',
        '完整 P2 再现依赖局部支持的多项式可辨识性与数值条件。规则网格面/棱/角处通过，不自动意味着不规则血栓点云、薄结构或任意单侧邻域也通过。未来使用实际点云时仍应逐节点检查九项秩、条件数、物理体积与解析场再现，保留失败节点，并依据既定触发条件决定是否需要 GMLS；本轮没有进行该工作，也没有自动扩大视界。')
    section('One next task only',
        '**HIGHER_ORDER_BAQP_ENERGY_FORCE_CANDIDATE_V1**。未来从同一高阶运动学定义能量并由其一阶变分推导力，再检查客观性、全局力/力矩、能量—力有限差分、仿射与二次边界 patch tests，之后才重新集成已验证的边界转换。\n\n'
        '本轮在运动学 oracle 报告完成后停止，等待人工审核；以上下一任务没有启动。')
    assert len(sections)==19
    header = '# HIGHER_ORDER_BAQP_RK_ORACLE_V1\n\n'
    header += '范围：KINEMATIC / ENERGY-DRIVER ORACLE ONLY。\n\n'
    header += f'开始状态见 [STARTING_GIT_STATE.json](STARTING_GIT_STATE.json)。配置和输入身份在计算前冻结，原始任务见 [TASK_REQUEST.txt](TASK_REQUEST.txt)。\n\n'
    (OUT / 'REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md').write_text(header+'\n'.join(sections))

    schema = {}
    for file in sorted(OUT.rglob('*.npz')):
        with np.load(file, allow_pickle=False) as a:
            schema[str(file.relative_to(OUT))] = {k:dict(shape=list(a[k].shape),dtype=str(a[k].dtype)) for k in a.files}
    for file in sorted(OUT.glob('*.csv')):
        with file.open() as f:
            reader=csv.reader(f);columns=next(reader);count=sum(1 for _ in reader)
        schema[file.name]=dict(columns=columns,rows=count)
    dump('DATA_SCHEMA.json',schema)
    dictionary = '''# Data dictionary and reproduction

All data describe fixed exact target configurations. No equilibrium or time history was generated.
Python IDs are zero-based. Tensor axes are node, displacement component, physical derivative.
The complete array shapes/dtypes and CSV headers/row counts are in [DATA_SCHEMA.json](DATA_SCHEMA.json).

| Files | Meaning / units |
| --- | --- |
| REFERENCE_GEOMETRY.npz | X in m; volume in m³; undirected pairs; directed source/target/ptr; saved energy weights; node/bond depth divided by δ; geometric classes |
| PD2_REFERENCE_CACHE.npz | M in m³, nine descending singular values in m³, rank, condition, coefficients in m⁻⁴, directed g in m⁻¹ including Vj |
| RK2_DIAG_REFERENCE_CACHE.npz | Same conventions for the six-term negative control |
| PD2_NODE_CONDITIONING.csv | One record per node, all singular values, condition flag and family size |
| PD2_CONDITIONING_GROUPS.csv | All-node, four geometric-class, and nearest-face-depth statistics |
| POLYNOMIAL_REPRODUCTION.npz / PD2_POLYNOMIAL_REPRODUCTION.csv | Exact/computed physical gradients and residuals; linear scale 1, quadratic scale L in m; all 6,859 × 9 entries |
| PD2_AFFINE_REGRESSION.npz / .json | Two affine computed/exact nodal gradients and normalized errors |
| *_NODAL_STATE.npz | X, x in m; analytic F and C1/PD2/RK2_DIAG nodal F, dimensionless |
| PD2_AXIAL_QUADRATIC.csv / PD2_SHEAR_QUADRATIC.csv / PD2_CROSS_TERM.csv | All nodes for all three reconstruction spaces; F error divided by 0.004 |
| NODE_GEOMETRY_STATISTICS.csv | Maximum and physical-volume weighted RMS per field/model/geometric class |
| F_DEPTH_PROFILES.csv | All-node/all-bond layer reductions; Fbar uses nodal volume, Ftilde uses original energy weights |
| *_BOND_DIAGNOSTICS.npz | All 690,020 bonds: F error is nonnegative Frobenius norm; Psi error is signed computed-minus-exact density in Pa |
| *_NODE_PSI.npz | Analytic positive-energy density and signed C1/PD2 errors in Pa |
| PD2_FTILDE_DIAGNOSTIC.csv | Bond-layer F errors, absolute and divided by 0.004 |
| PD2_PSI_PLUS_DIAGNOSTIC.csv | Node/bond-layer Psi errors; max and weighted RMS are absolute magnitudes normalized by the single analytic reference per field |
| C1_PD2_SUMMARY.json | Global max and weighted RMS, absolute and normalized; per-field Psi_ref in Pa |
| julia_crosscheck/PYTHON_REPRESENTATIVES.npz | Only 4 representative nodes and their 551 directed weights, not a duplicate of the entire reference cache |
| julia_crosscheck/*_F.csv | Four × nine row-major tensor entries; Python/Julia pairs |
| julia_crosscheck/julia_gradient_weights.csv | No header: original i/j in one-based indexing, then physical gx, gy, gz |
| PERIDYNAMICSJL_CROSSCHECK.json | Four representative tensors per field, exact package provenance and weight comparison |
| PROTECTED_INPUT_HASHES.json | Pre-run hashes of 496 protected historical files; old manifests remain untouched |
| OUTPUT_MANIFEST.json | SHA256 and byte size of every file in the three new namespaces except the manifest itself |

All layer reductions use the full relevant population, without random sampling. Geometry classes
count coordinate axes with nearest opposite-face distance < δ; they do not assert force-stencil closure.
The P2 cache is reference-only and displacement-independent. Cache/config hashes are recorded in
[CACHE_PROVENANCE.json](CACHE_PROVENANCE.json); input-array identity is in [INPUT_IDENTITY.json](INPUT_IDENTITY.json).

Review figures use only saved outputs. To regenerate them from the project root:

```bash
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \\
  /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B \\
  verification/higher_order_baqp_rk_oracle_v1/plot_figures.py
```

The scientific computation entry point is `run_audit.py`; the independent installed-package comparison
is `run_julia_crosscheck.py`. These are retained for reproducibility and were not rerun during report
generation. The 12 focused tests are in `tests/test_pd2_oracle.py`; their completed result is recorded
in [TARGETED_TESTS.json](TARGETED_TESTS.json). No full project suite was run.

Reference PDFs are retained with source URLs and SHA256 in [SOURCE_PROVENANCE.json](references/SOURCE_PROVENANCE.json).
The current audit does not replace or edit historical scientific manifests, active context, or results.
'''
    (OUT / 'DATA_DICTIONARY.md').write_text(dictionary)
    names = ['c1_vs_pd2_boundary_bias','pd2_polynomial_reproduction','pd2_conditioning','pd2_ftilde_and_energy_bias','python_vs_peridynamicsjl_pd2']
    readme = '# Higher-order RK kinematic oracle figures\n\n'
    readme += 'Published full-quadratic RK removes the identified boundary kinematic and positive-energy bias for the tested N19 fields. Production force is not implemented.\n\n'
    readme += '[Scientific review](../../verification/higher_order_baqp_rk_oracle_v1/REVIEW_HIGHER_ORDER_BAQP_RK_ORACLE_V1.md) · [Data dictionary](../../verification/higher_order_baqp_rk_oracle_v1/DATA_DICTIONARY.md)\n\n'
    readme += table(['Figure','600 dpi PNG','Vector PDF'],[[name,f'[PNG](figures/{name}.png)',f'[PDF](figures/{name}.pdf)'] for name in names])
    (RESULTS / 'README.md').write_text(readme+'\n')
    print('Report: 19 sections; status, GMLS decision, data schemas and figure index written.')

if __name__ == '__main__':
    main()
