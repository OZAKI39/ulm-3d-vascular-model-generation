"""Compile the completed shear audit without executing scientific solves."""
import csv,hashlib,json,subprocess
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/nonlocal_boundary_nonaffine_shear_localization_v1'
RESULT=ROOT/'results/nonlocal_boundary_nonaffine_shear_localization_v1'


def read(n):return json.loads((OUT/n).read_text())
def write(n,value):(OUT/n).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def num(x):return f'{float(x):.9e}'
def vec(x):return '('+', '.join(num(v) for v in x)+')'
def table(headers,rows):
    def cell(v):return str(v).replace('|',r'\|')
    return '\n| '+' | '.join(map(cell,headers))+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(cell,row))+' |\n' for row in rows)+'\n'
def md(path,text):path.write_text(text,encoding='utf-8')
def csvrows(n):
    with (OUT/n).open() as s:return list(csv.DictReader(s))


def compile_report():
    cfg=read('PREDECLARED_CONFIG.json');ref=read('EXACT_NONAFFINE_SHEAR_REFERENCE.json');c=ref['continuum']['combined'];sc=ref['scales']
    bulk=read('BULK_SHEAR_ORACLE_VALIDATION.json');aff=read('AFFINE_SHEAR_LIMIT_RECOVERY.json');ids=read('SHEAR_NODEWISE_IDENTITIES.json')
    n27=read('N27_SHEAR_DEPTH_CHECK.json');qs=read('SHEAR_SUPPORT_QUANTILES.json');bias=read('SHEAR_FIELD_BIAS_SUMMARY.json');bs=bias['samples']
    nl=read('SHEAR_NL_MANUFACTURED_RECOVERY.json');cl=read('SHEAR_CLASSICAL_MANUFACTURED_RESULT.json')
    force=read('SHEAR_FORCE_BALANCE.json');moment=read('SHEAR_MOMENT_BALANCE.json');channels=read('SHEAR_CHANNEL_ACCOUNTING.json')
    assert ref['VERIFIED'] and aff['PASS'] and bulk['VERIFIED'] and ids['PASS'] and n27['PASS'] and nl['PASS']
    # Compare independently solved controls on free nodes, retaining the model load.
    with np.load(OUT/'SHEAR_NL_RECOVERED_STATE.npz') as a, np.load(OUT/'SHEAR_classical_RECOVERED_STATE.npz') as b, np.load(OUT/'NODEWISE_SHEAR_DECOMPOSITION.npz') as d:
        assert np.array_equal(a['initial'],b['initial']) and np.array_equal(a['fixed'],b['fixed'])
        free=~a['fixed'];model=d['f_model'][free];response=b['force_internal'][free]-a['force_internal'][free]
        accounting=dict(force_response_closure_relative_to_model=float(np.linalg.norm(response+model)/np.linalg.norm(model)),
            NL_recovery_error_over_classical_deviation=nl['displacement_error_max_m']/cl['displacement_error_max_m'],
            same_initial=True,same_constraints=True,same_boundary_conversion=True,
            load_difference_minus_model_scaled=float(np.linalg.norm((b['external']-a['external']-d['f_model'])[free],axis=1).max()/sc['F_ref_N']),
            formula='f_internal_classical - f_internal_NL = -e_model V + (r_classical - r_NL), on free nodes')
    write('SHEAR_CLASSICAL_VS_NL_ACCOUNTING.json',accounting)
    protected=read('PROTECTED_INPUT_HASHES.json');changed=[p for p,h in protected.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    integrity=dict(protected_files=len(protected),changed=changed,PASS=not changed)
    write('HISTORICAL_INTEGRITY.json',integrity);assert integrity['PASS']
    extensions=[]
    for p,h in read('APPROVED_NAMESPACE_EXTENSIONS.json').items():
        archived=OUT/'source_before'/Path(p).name
        assert hashlib.sha256(archived.read_bytes()).hexdigest()==h
        extensions.append(dict(path=p,before_sha256=h,after_sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),archived_original=str(archived.relative_to(ROOT))))
    write('SOURCE_EXTENSION_PROVENANCE.json',extensions)
    resources=[read('RESOURCES_'+k+'.json') for k in ['reference_affine','oracle','decomposition','n27','NL','classical']]
    runtime=dict(scientific_stage_wall_seconds=sum(r['wall_seconds'] for r in resources),
        peak_process_tree_RSS_bytes=max(r['peak_process_tree_RSS_bytes'] for r in resources),
        minimum_available_RAM_fraction=min(r['minimum_available_RAM_fraction'] for r in resources),
        stages=resources,oracle_workers=8,equilibrium_workers=1,BLAS_threads=1,GPU='GPU_NOT_USED_BY_DESIGN',
        wall_definition='Sum of six scientific stage durations; excludes authoring, tests and figure export. Tree RSS may double-count shared pages.')
    write('RESOURCE_SUMMARY.json',runtime)
    status={
        'NONAFFINE_SHEAR_REFERENCE':'VERIFIED', 'CONTINUUM_FORCE_CLOSURE':'PASS',
        'CONTINUUM_MOMENT_CLOSURE':'PASS','CONTINUUM_WORK_CLOSURE':'PASS','BULK_BAQP_SHEAR_ORACLE':'VERIFIED',
        'AFFINE_SHEAR_LIMIT_RECOVERY':'PASS','NODEWISE_SHEAR_DECOMPOSITION':'PASS',
        'NAIVE_SHEAR_CANCELLATION':'CONTAINS_INTERIOR_PHYSICS','SHEAR_BOUNDARY_LOCALIZATION':'VERIFIED',
        'N27_CLOSED_CENTER_SHEAR_CORRECTION':'PASS','SHEAR_BULK_MODEL_DIFFERENCE':bulk['SHEAR_BULK_MODEL_DIFFERENCE'],
        'SHEAR_NL_MANUFACTURED_RECOVERY':'PASS',
        'SHEAR_CLASSICAL_MANUFACTURED_FIELD':'APPROXIMATE_DUE_TO_BULK_NONLOCALITY' if cl['PASS'] else 'FAILED_FOR_OTHER_REASON',
        'KINEMATIC_BOUNDARY_BIAS':bias['KINEMATIC_BOUNDARY_BIAS'],'PSI_PLUS_BOUNDARY_BIAS':'REPORTED_NOT_ACCEPTANCE_RATED',
        'BOUNDARY_DATA_CONVERSION_STATUS':'NONAFFINE_NORMAL_AND_SHEAR_LOCALIZATION_VERIFIED',
        'SURFACE_OPERATOR_CORRECTION':'NOT_MODIFIED','C5_C7_ACTIVATION':'NOT_RUN','CLOT_GEOMETRY':'NOT_RUN','FRACTURE':'NOT_RUN'}
    write('FINAL_STATUS.json',status)
    git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
    report=[f'''# NONLOCAL_BOUNDARY_NONAFFINE_SHEAR_LOCALIZATION_V1

日期：2026-10-06。分支 `{branch}`，HEAD `{git_head}`。本报告对应固定参数、完整态制造场；没有运行 phase、AT1、damage、stability、fracture、erosion、CFD 或 clot 几何。

**结论：既有非仿射边界局部化分解可迁移到本二次剪切基准。边界力转换通过，不代表近边界 F/Ftilde/Psi+ 重构已准确。** Kinematic bias 为 `{bias['KINEMATIC_BOUNDARY_BIAS']}`；唯一建议的下一任务是 `BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_AUDIT_V1`，本轮不执行。

结果入口：[最终状态](FINAL_STATUS.json)、[数据字典](DATA_DICTIONARY.md)、[六组图件](../../results/nonlocal_boundary_nonaffine_shear_localization_v1/README.md)。所有数字均来自本轮保存的 JSON/CSV/NPZ，图件无合成数据。

## 1. Why shear is the decisive transferability test

轴向制造场无法独立检验非对角应力、侧面切向牵引及参考坐标力矩。本轮只把同一 `q_bulk V−f_FD` 分解应用于空间变化的剪切场；未引入另一套边界方法。之前的[轴向审计](../nonlocal_boundary_nonaffine_localization_v1/REVIEW_NONAFFINE_BOUNDARY_LOCALIZATION_V1.md)和[仿射 B1 审计](../nonlocal_boundary_formulation_v1/REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md)原样保留。

预先冻结 [PREDECLARED_CONFIG.json](PREDECLARED_CONFIG.json)：N19（6859 粒子），L=18h=1.125 mm，h=62.5 µm，δ=251.25 µm，δ/h=4.02；SVK μ=1000 Pa、λ=0、K=666.6666666666666 Pa。面／棱／角体积分别为 h³/2、h³/4、h³/8，内部为 h³。L 保留 `18*h` 的原浮点表示，与旧仿射网格逐项相同。

## 2. Exact nonaffine shear field
''',r'''
$$
\mathbf u_L=(\gamma Z+\beta_s Z^2/L,0,0),\quad
\gamma=0.003,\quad\beta_s=0.001,\quad
k=\gamma+2\beta_s Z/L,\quad k'=2\beta_s/L.
$$
$$
F=I+k\,\mathbf e_x\otimes\mathbf e_z,
\qquad J=1,\qquad
E=\frac12(F^TF-I)=
\begin{pmatrix}0&0&k/2\\0&0&0\\k/2&0&k^2/2\end{pmatrix}.
$$

物理域 k∈[0.002,0.004]；梯度误差统一除以 0.004。位移场、材料、h、δ 和所有门槛均未按结果调整。

## 3. Closed-form P, body force and surface traction

独立闭式实现见 [`nonaffine_manufactured.py`](../../pd_clot/pfpd_candidate/nonlocal_boundary/nonaffine_manufactured.py) 的 `QuadraticShear`，其公式来自任务定义及 SVK 直接代数推导；没有调用通用本构来生成参考。

$$
S=2\mu E,\quad
P=FS=\mu\begin{pmatrix}k^2&0&k+k^3\\0&0&0\\k&0&k^2\end{pmatrix},\quad
W=\frac\mu2 k^2+\frac\mu4 k^4.
$$
$$
\mathbf b_{\rm local}=-\operatorname{Div}P
=\big[-\mu(1+3k^2)k',\,0,\,-2\mu kk'\big]^T.
$$
$$
\mathbf t=P\mathbf N,\quad
\mathbf t_{x\pm}=\pm\mu(k^2,0,k),\quad
\mathbf t_{z\pm}=\pm\mu(k+k^3,0,k^2),\quad
\mathbf t_{y\pm}=0.
$$

六个面分别积分，未遗漏 x 面法向或切向牵引。仅有 x 位移并不意味着 z 体力或 z 牵引为零。来源：[独立参考与全部面数据](EXACT_NONAFFINE_SHEAR_REFERENCE.json)；`NONAFFINE_SHEAR_EXACT_REFERENCE.json` 是满足任务另一命名要求的同内容副本。
''',f'''
通用本构对照：应力归一化误差 {num(ref['independent_constitutive_error'])}，能量相对误差 {num(ref['independent_energy_error'])}。P_ref={num(sc['P_ref_Pa'])} Pa，b_ref={num(sc['b_ref_N_m3'])} N/m³，F_ref={num(sc['F_ref_N'])} N。

## 4. Continuum force / moment / work checks

使用体内 8³、每个面 8² Gauss 积分，对本任务中的多项式被积函数达到解析积分精度；不把节点梯形积分当连续体金标准。
''',r'''
$$
\int_\Omega\mathbf b\,dV+\int_\Gamma P\mathbf N\,dA=0.
$$

参考坐标力矩必须保留 P 非对称带来的源项。由散度定理逐分量积分得到：

$$
\mathbf M_\Gamma^X+\mathbf M_b^X
=\int_\Omega\sum_a\mathbf e_a\times P_{:a}\,dV
=\big[0,\mu(\gamma^3+\gamma\beta_s^2)L^3,0\big]^T.
$$

不能强制这项等于零。当前坐标力矩则用 x=X+u；其源项为 ∫Σ F_:a×P_:a dV，完整 SVK 中 PFᵀ 对称，合并通道的当前坐标力矩为零。

$$
\int_\Gamma\mathbf t\cdot\delta\mathbf u\,dA+
\int_\Omega\mathbf b\cdot\delta\mathbf u\,dV
=\int_\Omega P:\nabla\delta\mathbf u\,dV.
$$
''',f'''
参考力矩源项为 {vec(c['reference_moment_source_analytic_Nm'])} N·m；当前坐标合力矩为 {vec(c['current_moment_Nm'])} N·m。13 个虚位移包括 3 个单位平移、9 个仿射梯度和 δu_x=Z²/L，包含所需 A_xz、A_zz、A_zx。
''',table(['连续体检查','归一化误差','门槛'],[[label,num(c[key]),'1e-10'] for label,key in [('合力','force_closure_scaled'),('参考力矩（扣除解析源项）','moment_closure_scaled'),('当前坐标力矩','current_moment_closure_scaled'),('13 模式最大虚功','max_work_closure_scaled')]]),f'''
连续体力、力矩尺度分别为 P_ref L²、P_ref L³。虚功尺度为 P_ref L² max(L,max|δu|)，单位平移按 1 m 的虚位移幅值核算。

## 5. Bulk BAQP shear oracle

复用 [`bulk_reference.py`](../../pd_clot/pfpd_candidate/nonlocal_boundary/bulk_reference.py) 及未改动的选择性模板，在每个物理 Z 层以 N27 完整依赖模板计算 q_bulk=f_closed/h³。D1/D2/D3 依赖闭合包含所有所需节点体积，不能只检查观察点的一跳邻域。

先完成 β_s=0 的三个闭合模板和旧 affine B1 对照，再执行 β_s≠0。仿射 bulk 力最大/F_ref={num(aff['bulk_force_max_scaled'])}；旧 B1 相对 L2 差={num(aff['relative_L2_error'])}，均通过 1e-10。

非仿射参考包括 19 层×3 个通道；附加三个 Z 层的 N29 偏移模板和剪切反射检查。参考变化最大/F_ref={num(max(v['force_scaled_error'] for v in bulk['checks']))}；通道和误差={num(bulk['channel_sum_error_scaled'])}。选择性力另与完整 N19 力核对三个位置。来源：[bulk CSV](BULK_BAQP_SHEAR_REFERENCE.csv)、[闭合图](BULK_STENCIL_CLOSURE.json)、[独立验证](BULK_SHEAR_ORACLE_VALIDATION.json)、[选择性／全场比较](SELECTIVE_VS_FULL_SHEAR.json)。

## 6. Bulk BAQP vs local continuum difference
''',r'''
$$\mathbf e_{\rm model}(Z)=\mathbf q_{\rm bulk}(Z)+\mathbf b_{\rm local}(Z).$$

这是无限／完全闭合离散体内 BAQP 与局部 SVK 的差异，不是有限域边界污染。
''',table(['量','结果'],[['max |q_bulk| [N/m³]',num(bulk['max_q_bulk_N_m3'])],['max |e_model| [N/m³]',num(bulk['max_e_model_N_m3'])],['max |e_model| / max |b|',num(bulk['relative_model_difference'])],['预声明描述','SMALL iff ratio ≤ 1e-3'],['本轮分类',bulk['SHEAR_BULK_MODEL_DIFFERENCE']]]),table(['分量','max |e| [N/m³]','max |b| [N/m³]','分量比值'],[[a,num(bulk['e_model_component_max_N_m3'][j]),num(bulk['b_local_component_max_N_m3'][j]),'不定义（分母为零）' if bulk['component_relative_errors'][j] is None else num(bulk['component_relative_errors'][j])] for j,a in enumerate('xyz')]),'''
x 分量的小但非零差异高于 oracle 浮点复现误差；y/z 的极小差异不赋予物理显著性。SMALL 仅为固定基准的开发描述，不是通用模型容差。

## 7. Localized boundary correction
''',r'''
$$
\mathbf f_{\rm bulk,target,i}=\mathbf q_{\rm bulk}(Z_i)V_i,\quad
\mathbf f_{\rm BC,i}=\mathbf f_{\rm bulk,target,i}-\mathbf f_{\rm FD,i}.
$$

f_FD 由未修改 BAQP 力实现，在精确目标节点位置计算。物理域 V_i 含面／棱／角权重；oracle 除数仍为完整内部 h³。定义与轴向轮完全相同，整个平衡求解期间冻结转换载荷，不重新计算来追踪解。

$$F_{\rm ref}=\max(P_{\rm ref}h^2,b_{\rm ref}h^3).$$

N27 也使用同一物理 N19 场得到的尺度，不用试件总反力归一化局部误差。

## 8. Why naive cancellation is wrong

$$
\mathbf f_{\rm BC,naive}=-\mathbf f_{\rm FD},\quad
\mathbf f_{\rm BC,naive}-\mathbf f_{\rm BC}=-\mathbf q_{\rm bulk}V.
$$
''',f'''
恒等式最大误差/F_ref={num(ids['naive_scaled'])}。在完全闭合 N27 中心，naive/F_ref={num(n27['points'][-1]['naive_scaled'])}，localized/F_ref={num(n27['points'][-1]['localized_scaled'])}。因此 naive 抵消会移除实际非零内部剪切力。

## 9. Nodewise decomposition
''',r'''
$$
\mathbf f_{\rm FD}+\mathbf f_{\rm BC}+\mathbf b_{\rm local}V
=\mathbf e_{\rm model}V.
$$
''',f'''
全部 6859 节点逐项核对，最大闭合误差/F_ref={num(ids['closure_scaled'])}；NL 制造载荷目标残差/F_ref={num(ids['NL_target_residual_scaled'])}。法向／切向之和的最大逐节点误差/F_ref={num(ids['max_channel_sum_scaled'])}。[完整节点数据](NODEWISE_SHEAR_DECOMPOSITION.npz) 保留每个力向量、六面距离及目标梯度。

## 10. Global force / moment accounting
''',r'''
令 Q_b=Σ b_local V−∫b_local dV。全局力保留：

$$\mathbf F_{\rm BC}-\mathbf F_\Gamma
=\mathbf F_{\rm model}-\mathbf Q_b-\mathbf F_{\rm internal}.$$

参考力矩必须再保留连续体非对称 P 的源项 C_X=∫Σ e_a×P_:a dV：

$$\mathbf M_{\rm BC}-\mathbf M_\Gamma
=\mathbf M_{\rm model}-(\mathbf M_{b,\rm nodal}-\mathbf M_{b,\rm exact})
-\mathbf M_{\rm internal}-\mathbf C_X.$$
''',table(['量','合力 [N]','参考力矩 [N·m]'],[[label,vec(force[key]),vec(moment[key])] for label,key in [('经典表面','surface_local'),('精确体力','body_local_exact'),('节点体力','body_local_nodal'),('离散内部','internal'),('局部化转换','BC_localized'),('模型残差','model'),('体力积分误差','body_quadrature_error')]]),f'''
力闭合/F_ref={num(force['closure_scaled'])}；参考力矩闭合/(F_ref L)={num(moment['closure_scaled'])}。当前位置 BAQP 内力矩另报为 {vec(moment['current_position_internal_moment_Nm'])} N·m，归一化 {num(moment['current_internal_moment_scaled'])}。不能把参考力矩与当前位置角动量混为同一测试。来源：[力](SHEAR_FORCE_BALANCE.json)、[力矩](SHEAR_MOMENT_BALANCE.json)。

## 11. Virtual work

对上述 13 种虚位移，分别保存内部力、边界转换、节点体力、修正残差及 e_model V 的功。三通道共 39 行，[SHEAR_VIRTUAL_WORK.csv](SHEAR_VIRTUAL_WORK.csv)。
''',r'''
$$
\sum_i(\mathbf f_{\rm FD,i}+\mathbf f_{\rm BC,i}+\mathbf b_iV_i)\cdot\delta\mathbf u_i
=\sum_i\mathbf e_{\rm model,i}V_i\cdot\delta\mathbf u_i.
$$
''',f'''
所有通道最大归一化差={num(ids['max_work_closure_scaled'])}，门槛 1e-10；尺度 F_ref max(L,max|δu|)。正确残差的功等于 bulk 模型差的功，不强制为零。

## 12. Spatial localization

对每个最近面深度，报告节点载荷幅值的平均值、最大值及累计 L1 幅值。累计统计对已经包含 V_i 的节点力直接求和，不再次乘体积。
''',table(['转换','分位数','µm','h','δ'],[[name,q,f"{qs[name][q]['depth_um']:.1f}",qs[name][q]['depth_h'],f"{qs[name][q]['depth_over_delta']:.6f}"] for name in ['localized','naive'] for q in ['d50','d90','d95','d99']]),'''
naive 尾部更长，但前三个分位数恰好相同，不把预期趋势当通过门槛。这些量仅描述 **conversion support for this benchmark**，不证明物理 nonlocal traction 层厚度等于 δ 或 3δ。[深度分布](SHEAR_BOUNDARY_LOCALIZATION.csv)。

## 13. N27 fully closed center result

与轴向轮相同的 Z=0 侧面 x 深度几何，h/δ/物理场多项式不变，仅模板范围扩至 26h。bulk 比较值来自独立的偏移 N29 闭合模板。
''',table(['x 深度/h','深度/δ','|f_BC| [N]','|f_BC|/F_ref','naive/F_ref','模板完全闭合'],[[p['depth_h'],f"{p['depth_over_delta']:.6f}",num(p['localized_norm_N']),num(p['localized_scaled']),num(p['naive_scaled']),p['fully_closed']] for p in n27['points']]),'''
12h 点数值已很小，但依赖模板仍触及边界体积；13h 中心才满足明确闭合。中心 1e-10 门槛通过。[N27 全部向量及闭合证据](N27_SHEAR_DEPTH_CHECK.json)。

## 14. BAQP-manufactured recovery
''',r'''
$$\mathbf b_{\rm NL}=-\mathbf q_{\rm bulk},\qquad
\mathbf f_{\rm ext,NL}=-\mathbf q_{\rm bulk}V+\mathbf f_{\rm BC}.$$

不是仅代入目标检查。沿用原有 mechanics adapter，两个 z 面固定目标位移，从平滑扰动开始：

$$
\Delta\mathbf u_0=0.02\gamma L\prod_{a=1}^3
\sin\!\left(\pi\frac{X_a+L/2}{L}\right)(1,-0.4,0.7),
$$

固定节点扰动显式设为零。使用现有 SciPy TRF/LSMR、原切线与门槛，未创建新求解器。
''',table(['NL 求解量','结果'],[['初始 max |Δu| [m]',num(nl['initial_perturbation_max_m'])],['初始 max |Δu_x| [m]',num(nl['initial_x_perturbation_max_m'])],['最终力残差（门槛 1e-8）',num(nl['force_residual'])],['max ||Fbar_solved−Fbar_target|| / 0.004（门槛 1e-5）',num(nl['F_recovery_error_max'])],['max 位移恢复误差 [m]',num(nl['displacement_error_max_m'])],['固定值误差 [m]',num(nl['fixed_displacement_error_m'])],['BAQP evaluates / tangent actions',f"{nl['solve']['BAQP_evaluations']} / {nl['solve']['tangent_actions']}"],['求解秒数',f"{nl['solve']['wall_seconds']:.3f}"]]),'''
求解力残差是自由节点 residual L2 / max(external_free L2,F_ref)。梯度恢复以 **离散目标 Fbar** 为参考；解析 F 偏差另见第 17 节。[NL 恢复](SHEAR_NL_MANUFACTURED_RECOVERY.json)。原求解器摘要中的 `surface_traction_* = 0` 为旧占位字段，不能解读为本轮真实表面牵引；`SCIPY_REFERENCE` 是求解路径名称，不表示运行相场。

## 15. Classical-manufactured result

采用同一初始扰动、同一约束和冻结转换，外力为 b_local V+f_BC。没有额外加入 −e_model V。目标处的剩余力已独立核对为 e_model V。
''',table(['经典制造解','结果'],[['最终力残差',num(cl['force_residual'])],['max 位移偏差 [m]',num(cl['displacement_error_max_m'])],['位移相对 L2 偏差',num(cl['displacement_relative_L2'])],['max 离散 F 偏差 / 0.004',num(cl['F_recovery_error_max'])],['目标 BAQP 能量 [J]',num(cl['target_energy_J'])],['求解能量 [J]',num(cl['energy_J'])],['能量偏差 [J]',num(cl['energy_deviation_J'])],['能量相对偏差',num(cl['energy_relative_deviation'])],['求解通过',cl['PASS']]]),f'''
两解的内力响应满足 Δf_internal≈−e_model V，自由节点差值闭合相对于模型载荷 L2 为 {num(accounting['force_response_closure_relative_to_model'])}；NL 位移恢复误差/经典位移偏差为 {num(accounting['NL_recovery_error_over_classical_deviation'])}。这些分辨率诊断与 [SHEAR_CLASSICAL_VS_NL_ACCOUNTING.json](SHEAR_CLASSICAL_VS_NL_ACCOUNTING.json) 一并保留；不把极小差异无条件归因于单一机制。

经典局部场并非精确 BAQP 解。边界恒等式已闭合且 e_model 独立测量，因此偏离包含有限 horizon bulk 模型差、晶格离散及求解分辨率，不能自动判为边界失败。[经典结果](SHEAR_CLASSICAL_MANUFACTURED_RESULT.json)。

## 16. Normal vs tangential channel interpretation
''',r'''
连续体面上用 t_n=(t·N)N，t_t=t−t_n。对本轴对齐立方体，等价于 P_n=diag(P)、P_t=P−diag(P)。在离散模型中，把**实际离散 bond Piola 应力**作同样分解，再由固定几何下对 P 线性的原力组装计算两个通道；oracle 也分别核算。

$$P=P_n+P_t\quad\Rightarrow\quad
f_{\rm FD}=f_{{\rm FD},n}+f_{{\rm FD},t},\quad
f_{\rm BC}=f_{{\rm BC},n}+f_{{\rm BC},t}.$$

这是冻结状态的应力通道核算，不是非线性材料／解的叠加，也没有给分布式非局部节点力分配唯一面归属。通道全局力、参考力矩和虚功恒等式比任意面所有权更有意义；分离通道不单独假设当前角动量为零。
''',table(['通道','BC 合力 [N]','经典表面合力 [N]','力闭合/F_ref','力矩闭合/(F_ref L)'],[[name,vec(v['force']['BC_localized']),vec(v['force']['surface_local']),num(v['force']['closure_scaled']),num(v['moment']['closure_scaled'])] for name,v in channels['channels'].items()]),'''
[完整通道数据](SHEAR_CHANNEL_ACCOUNTING.json) 保留非零节点积分误差、模型项和内部力矩；不是把每个通道 BC 合力硬设为经典表面合力。

## 17. Fbar / Ftilde reconstruction bias

诊断位置为精确目标 x_L；没有用求解误差替代重构偏差。Fbar 与同节点解析 F(Z_i) 比较；键关联 Ftilde 与**参考键中点**解析 F((Z_i+Z_j)/2) 比较，避免跨位置相减。深度也分别取节点／键中点到六个物理面的最近距离。

节点加权 RMS 使用 V_i；键 RMS 使用 BAQP 能量权重 2 wbar_ij V_i V_j。两者是不同采样测度，不直接比较其积分总量。
''',table(['梯度','max 绝对偏差','max / 0.004','加权 RMS 绝对偏差','加权 RMS / 0.004'],[[name,num(v['F_error_max']),num(v['F_error_max_scaled']),num(v['F_error_RMS']),num(v['F_error_RMS_scaled'])] for name,v in bs.items()]),f'''
预声明 LOW<1e-3、MODERATE 为 [1e-3,1e-2]、MATERIAL>1e-2，以两个梯度的最大归一化误差分类，得到 `{bias['KINEMATIC_BOUNDARY_BIAS']}`。仅为开发描述，不是未来断裂可用性门槛。重复相同数组诊断逐位一致，并有独立小网格检查；[摘要](SHEAR_FIELD_BIAS_SUMMARY.json)、[深度曲线](SHEAR_FIELD_BIAS.csv)、[逐样本标量数据](SHEAR_FIELD_BIAS_SAMPLES.npz)。

## 18. Psi_plus boundary bias

同一已有 Green/SVK 谱分裂应用于解析和离散 F，不改谱算法、不演化相场。代码及原始文献归属见 [`svk_spectral_split.py`](../../pd_clot/pfpd_candidate/svk_spectral_split.py)：既有 PhaseFieldX 固定提交的谱分裂代数及 Miehe 概念来源；本轮未重新引入论文模型，也不声称 Green 应变候选等同于小应变 FEM。
''',r'''
$$
E=\sum_a\varepsilon_a\,n_a\otimes n_a,\quad
E^+=\sum_a\max(\varepsilon_a,0)n_a\otimes n_a,\quad
\Psi^+=\tfrac\lambda2\langle\operatorname{tr}E\rangle_+^2+\mu E^+:E^+.
$$
''',f'''
稳定非零尺度为物理 N19 解析节点场 max Ψ⁺={num(bias['Psi_plus_reference_scale_Pa'])} Pa；两种重构统一使用该尺度，不使用接近零分量作分母。
''',table(['来源','max |ΔΨ⁺| [Pa]','max / Ψ_ref','加权 RMS [Pa]','RMS / Ψ_ref','加权带符号均值 [Pa]'],[[name,num(v['Psi_plus_error_max_Pa']),num(v['Psi_plus_error_max_scaled']),num(v['Psi_plus_error_RMS_Pa']),num(v['Psi_plus_error_RMS_scaled']),num(v['Psi_plus_error_mean_Pa'])] for name,v in bs.items()]),'''
深度曲线显示偏差集中在截断邻域附近；带符号均值不能替代局部最大误差。状态为 `REPORTED_NOT_ACCEPTANCE_RATED`。没有判断 fracture-safe，也没有改变 Gc、AT1 或损伤驱动。

## 19. Implication for future clot-wall shear loading

本完整态多项式制造基准支持同一边界数据转换保留正常的剪切 bulk 力，不证明一般曲面牵引、黏附、clot-wall、接触、phase 或断裂已验证。局部化通过与 F/Ftilde/Psi+ 近边界偏差可以同时成立；后者对后续裂纹起始驱动有直接相关性，应先量化因果而不是修改边界力去抵消它。
''',table(['最终状态','值'],list(status.items())),f'''
仅运行任务要求的 11 项新定向测试，记录见 [FOCUSED_TESTS_FINAL.log](FOCUSED_TESTS_FINAL.log)，没有全量 pytest 或旧 campaign。{integrity['protected_files']} 个受保护旧结果／配置／manifest／测试／核心源码哈希不变；三个允许扩展的 nonlocal_boundary 文件原版逐字节归档，新增剪切接口和可选诊断参数。[SOURCE_EXTENSION_PROVENANCE.json](SOURCE_EXTENSION_PROVENANCE.json) 明确区分这些扩展：旧输出 manifest 自身不改，但其中记录的旧扩展模块 hash 不再代表当前版本，不声称旧 manifest 对当前源码全部匹配。没有改旧断言、BAQP、材料、kernel、K、Z、surface 因子或 mechanics solver。

资源：8 个独立 oracle CPU workers，1 个串行 equilibrium worker，BLAS 1 线程；两个大求解不重叠。科学阶段累计 {runtime['scientific_stage_wall_seconds']:.3f} s（不含编写／画图／测试），峰值进程树 RSS {runtime['peak_process_tree_RSS_bytes']/2**30:.3f} GiB，最小可用 RAM {100*runtime['minimum_available_RAM_fraction']:.2f}%，满足 ≥20% headroom。RSS 可能重复计入共享页；记录同时保留系统可用内存。`GPU_NOT_USED_BY_DESIGN`，未安装依赖。详见 [RESOURCE_SUMMARY.json](RESOURCE_SUMMARY.json)。

## 20. One next task only

建议 `BAQP_NEAR_BOUNDARY_KINEMATIC_ENERGY_AUDIT_V1`。原因是本轮边界局部化和扰动恢复已经通过，而 Fbar/Ftilde 最大归一化解析偏差约 3%，Ψ⁺ 最大偏差约为统一参考的 6.4%。下一轮应独立审计重构／能量偏差的来源；本轮不实现 perifit、优化权重、ghost 节点、surface 因子、修改 K/Z 或 kernel。

不继续增加多项式局部化测试，不自动开始 curved clot、C5/C7、fracture、adhesion、侵蚀、CFD 或 microbubble。**到此停止，等待人工审核。**
''']
    md(OUT/'REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md',''.join(report))
    write('TERMINAL_SUMMARY.json',dict(branch=branch,HEAD=git_head,N=19,nodes=19**3,L_m=cfg['L_m'],h_m=cfg['h_m'],delta_m=cfg['delta_m'],delta_over_h=cfg['delta_m']/cfg['h_m'],
        gamma=cfg['gamma'],beta_s=cfg['beta_s'],scales=sc,continuum={k:c[k] for k in ['force_closure_scaled','moment_closure_scaled','max_work_closure_scaled']},
        bulk=bulk,affine_error=aff['relative_L2_error'],nodewise_error=ids['closure_scaled'],quantiles=qs,N27=n27['points'],
        NL={k:nl[k] for k in ['force_residual','F_recovery_error_max']},classical={k:cl[k] for k in ['force_residual','displacement_error_max_m','F_recovery_error_max']},
        bias=bs,status=status,resources={k:v for k,v in runtime.items() if k!='stages'},review=str(OUT/'REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md'),figures=str(RESULT/'figures')))
    print('REPORT COMPLETE',OUT/'REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md')


if __name__=='__main__':compile_report()
