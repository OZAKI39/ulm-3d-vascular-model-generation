#!/usr/bin/env python3
"""Run final offline Particle-2 checks and write evidence; never grant manual PASS."""
from pathlib import Path
import argparse
import csv
from datetime import datetime, timezone
from importlib.metadata import version
import json
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET
import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
REPO = PACKAGE.parent
REPORT = PACKAGE / "reports/particle2"
DATA = REPORT / "data"
LOGS = REPORT / "logs"
sys.path.insert(0, str(PACKAGE / "src"))
from particle_3d.audit import read_frozen, sha256, FEM_BRANCH
from particle_3d.particle2_audit import check_dependencies, P0_COMMIT, P1_COMMIT, FROZEN_COMMIT, P2_BRANCH
from particle_3d.rbc_distribution import load_contract, DEFAULT_CONTRACT
from particle_3d.particle2_cases import NORM_ATOL, IDENTITY_RELATIVE_BUDGET


def git(*args, cwd=REPO):
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_data(name):
    return json.loads((DATA / name).read_text())


def run(name, args, cwd):
    print(f"Checking {name}", flush=True)
    result = subprocess.run(args, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    path = LOGS / f"{name}.log"
    path.write_text(result.stdout)
    return dict(name=name, command=args, cwd=str(cwd), returncode=result.returncode, log=str(path.relative_to(REPO)))


def junit(path):
    cases = list(ET.parse(path).iter("testcase"))
    failed = sum(c.find("failure") is not None or c.find("error") is not None for c in cases)
    skipped = sum(c.find("skipped") is not None for c in cases)
    return dict(total=len(cases), passed=len(cases)-failed-skipped, failed=failed, skipped=skipped)


def numerical_evidence(shear, wide, real):
    """Aggregate full-run synthetic maxima and independently recheck saved real rows."""
    pop = np.genfromtxt(DATA / "C57BL6_RBC_GEOMETRY_VALIDATION_100000.csv", delimiter=",", names=True)
    reconstructed = (4/3)*np.pi*pop["a_m"]*pop["b_m"]*pop["c_m"]
    volume_error = float(np.max(np.abs(reconstructed/pop["volume_m3"]-1)))
    geometry = dict(N=len(pop), all_numeric_finite=bool(all(np.isfinite(pop[n]).all() for n in pop.dtype.names)),
                    axes_positive=bool(np.all(pop["c_m"]>0)), a_equals_b=bool(np.array_equal(pop["a_m"],pop["b_m"])),
                    c_less_than_a=bool(np.all(pop["c_m"]<pop["a_m"])),
                    r_in_open_zero_one=bool(np.all((pop["r"]>0)&(pop["r"]<1))),
                    lambda_in_open_minus_one_zero=bool(np.all((pop["jeffery_lambda"]>-1)&(pop["jeffery_lambda"]<0))),
                    max_volume_reconstruction_relative_error=volume_error)
    qerr=max(m["max_quaternion_norm_error"] for m in shear+wide+real)
    perr=max(m["max_p_norm_error"] for m in shear+wide+real)
    identity=max(m["max_jeffery_identity_error"] for m in shear+wide)
    real_identity=0.; real_scaled=0.; real_omega_error=0.
    for path in [DATA/"04_static.csv",DATA/"05_rotation.csv"]+sorted(DATA.glob("08_real_g*_dt*.csv")):
        with path.open(newline="") as stream:
            rows=list(csv.DictReader(stream))
        p=np.array([[float(row[f"p_{a}"]) for a in "xyz"] for row in rows])
        q=np.array([[float(row[f"q_{a}"]) for a in "wxyz"] for row in rows])
        g=np.array([[[float(row[f"G_{i}{j}_s_inv"]) for j in range(3)] for i in range(3)] for row in rows])
        omega=np.array([[float(row[f"Omega_{a}_s_inv"]) for a in "xyz"] for row in rows])
        curl=np.array([[g0[2,1]-g0[1,2],g0[0,2]-g0[2,0],g0[1,0]-g0[0,1]] for g0 in g])
        lam=np.array([float(row["jeffery_lambda"]) for row in rows])[:,None]
        e=.5*(g+g.transpose(0,2,1)); w=.5*(g-g.transpose(0,2,1))
        ep=np.einsum("nij,nj->ni",e,p)
        rhs=np.einsum("nij,nj->ni",w,p)+lam*(ep-np.sum(p*ep,axis=1)[:,None]*p)
        errors=np.max(np.abs(np.cross(omega,p)-rhs),axis=1)
        local=float(np.max(errors)); identity=max(identity,local)
        qerr=max(qerr,float(np.max(np.abs(np.linalg.norm(q,axis=1)-1))))
        perr=max(perr,float(np.max(np.abs(np.linalg.norm(p,axis=1)-1))))
        if path.name.startswith("08_"):
            real_identity=max(real_identity,local)
            real_scaled=max(real_scaled,float(np.max(errors/np.maximum(1.,np.linalg.norm(g,axis=(1,2))))))
            real_omega_error=max(real_omega_error,float(np.max(np.abs(omega-(.5*curl+lam*np.cross(p,ep))))))
    return dict(accepted_geometry_checks=geometry,max_quaternion_norm_error=qerr,max_p_norm_error=perr,
                max_jeffery_identity_error=identity,jeffery_identity_error_unit="s^-1; maximum absolute component",
                real_max_jeffery_identity_error_s_inv=real_identity,
                real_max_identity_error_scaled_by_max_one_gradient_frobenius=real_scaled,
                real_max_angular_velocity_formula_error_s_inv=real_omega_error,
                max_jeffery_period_relative_error=max(m["relative_error"] for m in shear+wide),
                distribution_wide_max_period_relative_error=max(m["relative_error"] for m in wide),
                max_error_scope="All synthetic steps via saved maxima; every saved static/rigid/real step independently rechecked.")


def statistics_table(stats):
    lines=["| 模型最终接受量 | mean | SD | 1% | 5% | 25% | 50% / median | 75% | 95% | 99% |",
           "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for key,label in [("D_um","D (µm)"),("V_fL","V (fL)"),("a_um","a=b (µm)"),("c_um","c (µm)"),
                      ("full_thickness_um","2c (µm)"),("r","r"),("jeffery_lambda","λ")]:
        s=stats[key]
        numbers=[s["mean"],s["sd"]]+[s["quantiles"][p] for p in ["1","5","25","50","75","95","99"]]
        lines.append("| "+label+" | "+" | ".join(f"{x:.8f}" for x in numbers)+" |")
    return "\n".join(lines)


def review(v,stats,static,rotation):
    tests=v["tests"];status=v["automated_checks"];initial=v["real_initialization"]
    comparisons=v["validation_timestep_comparison"]
    rows=[]
    for g in range(5):
        coarse,middle,_=comparisons[3*g:3*g+3]
        values=[np.degrees(s[k]) for s in [coarse,middle] for k in ["max_axis_difference_rad","median_axis_difference_rad"]]
        rows.append(f"| {coarse['r']:.8f} | "+" | ".join(f"{x:.6f}" for x in values)+" |")
    comparison_table="\n".join(rows)
    figure_notes=[
        "应该看什么：蓝色是文献报告，橙色是 V0 假设，绿色是推导几何。实际看到了什么：从 D/V 抽样到 a、b、c 的来源可追溯，PREPRINT 和模型 CV 标记清楚。有没有异常：文献 RDW 不统一，14.8% 保留为用户指定模型选择，不能写成该直径样本的实测值。",
        "应该看什么：D/V 柱状图与 latent、guard-only 曲线的区别。实际看到了什么：最终接受分布在小直径尾部受形状拒绝影响，平均值和 SD 与 latent 参数稍有差异。有没有异常：这是拒绝采样的预期结果，保护范围不是实验极值，未强迫 accepted SD 等于 0.93。",
        "应该看什么：a=b、c、完整厚度 2c、r 和 λ 的分布。实际看到了什么：全部样本满足 a=b>c>0、0<r<1 和 −1<λ<0。有没有异常：c/r/λ 均为模型推导量，不能当成文献直接测量的厚度或形状分布。",
        "应该看什么：D 对 V、2c、r 的散点，以及被拒绝的红叉。实际看到了什么：c≥a 的组合被原样拒绝，较小 D 和较大 V 更容易过厚。有没有异常：最终样本不再被要求独立；散点只作确定性抽样显示，统计使用全部 100,000 个样本。",
        f"应该看什么：零梯度时姿态是否保持不变。实际看到了什么：五种样本的旋转矩阵误差为 {static['max_rotation_matrix_error']:.1e}，平移最大误差 {static['max_position_error_m']:.3e} m。有没有异常：没有数值异常；这里静止指流体不旋转且无梯度，另有恒定平移来同时验证 V=u。",
        f"应该看什么：p 三分量与解析整体旋转是否重合，以及 q 的单位长度。实际看到了什么：五种 r 的结果一致，最大轴分量误差 {rotation['max_axis_vector_error']:.3e}。有没有异常：误差在预先设定的舍入预算内，未漏掉 curl 的 1/2。",
        f"应该看什么：五种不同 r 的翻转速度，以及测量周期和解析周期。实际看到了什么：较小 r 周期更长，三种验证 dt 的最大周期相对误差为 {v['max_jeffery_period_relative_error']:.3e}。有没有异常：没有超过预设预算；图中是 p 与 −p 等价的形状周期，有向 p 的完整周期还要乘 2。",
        f"应该看什么：覆盖分布的 64 个样本是否沿解析曲线排列。实际看到了什么：全部通过，最大周期相对误差 {v['distribution_wide_max_period_relative_error']:.3e}。有没有异常：未发现仅中等形状可用而较薄/较厚样本失效的情况；这不是 64 个相互作用粒子的模拟。",
        "应该看什么：真实血管、INLET/三个 OUTLET、中心线和蓝色短轴箭头。实际看到了什么：图示中位 r 样本的最细验证重放，五种形状中心线重合，15 次均从 OUTLET_02 离开。有没有异常：中心检查通过；箭头长度为显示比例且不是速度或物理半径，整个椭球的壁面间隙仍未验证。",
        "应该看什么：短轴分量、q 长度误差、角速度大小及 tetra 切换。实际看到了什么：各状态有限，q 误差在浮点舍入量级；角速度随原始梯度切换而有细碎变化。有没有异常：没有平滑 G、vorticity 或 Ω；有限性和单位长度并不证明真实姿态时间精度。",
        f"应该看什么：粗/中步长相对最细重放的短轴夹角差。实际看到了什么：最大差从 {v['real_dt_max_difference_coarse_deg']:.3f}° 降到 {v['real_dt_max_difference_middle_deg']:.3f}°，最薄样本的中位差略增，其余中位差下降。有没有异常：局部差异仍明显，真实姿态对 dt 敏感，未达到生产收敛证明；最细步长与自身为零仅是参照定义。"
    ]
    figures="\n\n".join(f"![{Path(path).stem}](figures/{Path(path).name})\n\n{note}" for path,note in zip(v["figures"],figure_notes))
    content=f"""# Particle-2 人工审核报告

## 1. 这一阶段做了什么

第一次建立 C57BL/6 RBC 尺寸分布，再让分布中不同形状的单个刚性 RBC 平移和旋转。
100,000 个样本只验证几何分布；五个代表样本、64 个分层样本及真实 FEM 重放都是独立单 RBC 算例。
没有用一只“平均 RBC”替代分布。自动检查 **{status}**，本阶段人工审核仍待用户。

分支：`{v['git_branch']}`。被测试的实现 commit：`{v['git_commit']}`。
机器记录指向已提交源码；之后仅保存报告和日志的证据提交不改变该实现。
P1 用户审核先以 `{v['particle1_acceptance_evidence_commit']}` 单独提交：USER_CHAT_REVIEW、PASS、blocking_issue=NONE。
P0 依赖 `{P0_COMMIT}`；P1 依赖 `{P1_COMMIT}`；FEM base `{FROZEN_COMMIT}`。
P0/P1 旧源码、测试、图片及数据逐文件 SHA 核验一致；FEM 冻结数据未变，也没有运行 CFD。

## 2. 文献真正告诉了我们什么

Moss 等 2025 年的 [PREPRINT](https://pmc.ncbi.nlm.nih.gov/articles/PMC12132290/) Table 3 报告
WT C57BL/6Case 成熟 RBC 共 1,156,720 个：直径 mean=6.79、SD=0.93、mode=6.67、median=6.67、IQR=1.33 µm。
Table 4 的 MCV 是 47.9±2.6 fL，其中 ±2.6 是不同小鼠 MCV 的变化，不能当成单细胞体积 SD。
这份来源不是经过同行评审的最终分布结论。

[Rivera 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3656420/) Table 1 的 C57BL/6J 雌/雄 MCV 为 47.8/48.4 fL；
[De Franceschi 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1895196/) Table 1 的 C57BL6 control 为 MCV 49.3±0.9 fL、RDW 13.2±0.9%。
后两者为同行评审文献，只交叉核查量级；不同亚系和测量体系没有被当成完全相同。
原文复核还发现 Moss Table 4 的 RDW=17.3±0.5%，**并非本项目的 14.8%**。
14.8% 按用户明确批准的 V0 选择保留，角色是 MODEL CHOICE，不冒充直径数据集的单细胞体积测量。

直接报告量是上述文献统计；正态分布、CV=14.8%、D/V 独立和三个标准差范围是模型近似/假设；
c、2c、r、λ 及 accepted population 的统计是模型推导结果。
主要 D 来源是 C57BL/6Case，交叉核查涉及其它 C57BL/6 亚系，合同标记 MIXED_SOURCE_WITHIN_C57BL6_FAMILY。
完整文献信息及 source role 见 [sources.json](literature/sources.json) 和 [文献审核](00_distribution_literature_report.md)，未复制论文全文。

## 3. RBC distribution 怎么生成

先独立抽取 D~Normal(6.79,0.93²) µm 和 V~Normal(47.9,7.0892²) fL；每个候选同时具有 D/V。
固定 chunk=4096，按原顺序依次做 D guard、V guard、算 a=b=D/2 与 c=3V/(πD²)，再检查 c<a。
D guard=[4.00,9.58] µm，V guard=[26.6324,69.1676] fL，都是 3-SIGMA PARAMETRIC MODEL GUARD，**不是实验 min/max**。
过厚的 c≥a 组合整对丢弃，继续抽下一个；不 clip c，不改 D 或 V。
D/V 独立是缺少 matched single-cell joint data 时的临时假设，形状拒绝后可产生相关性。

合同 [C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json](../../contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json)
SHA256：`{v['distribution_contract_sha256']}`。
N={v['validation_population_N']}，seed={v['validation_population_seed']}，角色 RBC_GEOMETRY_VALIDATION_ONLY；不是 hematocrit 或生产粒子群。
实际检查 {v['candidate_count']} 个候选，D guard 拒绝 {v['diameter_guard_rejections']}，V guard 拒绝 {v['volume_guard_rejections']}，
合计 guard 拒绝 {v['guard_rejection_count']}，shape 拒绝 {v['shape_rejections']}，接受率 {v['acceptance_rate']:.10%}。
这些拒绝计数按先后顺序互斥；整块生成 {v['generated_candidate_count']} 个候选，末尾 {v['unused_generated_tail_count']} 个未检查，不计入 candidate_count。

进入状态前 axes×1e−6 转成 m，V×1e−18 转成 m³；1 fL=1 µm³。
所有 accepted geometry 均 finite、a=b>c>0，体积重建最大相对误差 {v['accepted_geometry_checks']['max_volume_reconstruction_relative_error']:.3e}。
同 N/seed/环境的两个独立进程输出 float64 数组、CSV 和候选台账完全一致：{v['reproducibility_status']}。
环境 Python {v['environment']['python_version']}、NumPy {v['environment']['packages']['numpy']}、PCG64；没有声称跨 NumPy 版本逐位一致。
人口 CSV SHA：`{v['population_sha256']}`，相邻 metadata 保存 D/V/axes 数组 SHA、来源、单位和版本。

## 4. 分布长什么样

下表来自全部最终接受样本，SD 使用 ddof=1，50% 等于 median；这些是 **FINAL-DERIVED statistics**。
保护范围和形状拒绝改变了最终分布，所以不要求最终 SD 精确等于 latent SD。
accepted D/V 的样本相关系数为 {stats['accepted_D_V_sample_correlation']:.8f}；这不是生物学独立性证据。

{statistics_table(stats)}

![01 直径和体积](figures/01_rbc_diameter_volume_distribution.png)

![02 推导几何](figures/02_rbc_derived_geometry_distribution.png)

![03 联合散点及拒绝候选](figures/03_rbc_joint_geometry_scatter.png)

图 01 保留 latent 与只做 guard 的目标曲线，图 02 区分半厚度 c 和完整厚度 2c。
图 03 的红叉保留原始拒绝值，绘图只确定性抽取部分点，统计不抽样。
latent 和 guard-only 在形状拒绝之前分别用独立 erf CDF 做 smoke 检查，预设 DKW α=1e−6，全部通过；
没有把最终 accepted 分布拿去强行拟合未拒绝的正态分布。

## 5. RBC 姿态怎样表示

RBCGeometry 保存 SI 轴长、体积、r、λ 与来源；RBCState 保存编号、位置、四元数、速度、角速度及 geometry。
状态不可变，数组使用独立只读内存，正式几何由 sampler provenance 构造；手填轴长只允许明确标记 SYNTHETIC_ONLY 的数学测试。
四元数固定 **q=(w,x,y,z)**，R(q) 把 body 坐标变到 world 坐标。
body 短轴 e3=(0,0,1)，世界短轴 **p=R(q)e3**，表示厚度方向。
q 与 −q 是同一旋转；轴对称椭球的 p 与 −p 是同一形状朝向，这是两种不同的等价关系。
每步用 world Ω 的增量从左乘：q_new=normalize(delta_q ⊗ q_old)，不使用欧拉角。
显示时选择与前一 q 点积非负的等价符号，不改变物理旋转；Ω=0 时姿态不变。
真实初始短轴为 (1,2,3)/sqrt(14)，仅为 VALIDATION_INITIAL_ORIENTATION_ONLY，没有冻结入口姿态分布。

## 6. RBC 为什么会旋转

当地流体整体转动会带动 RBC；当地的拉伸和剪切也会改变扁球短轴方向，影响大小由本个样本的 r 决定。
G 分为 E=(G+G.T)/2 的拉伸/剪切部分与 W=(G−G.T)/2 的整体旋转部分。
r=c/a<1，λ=(r²−1)/(r²+1)<0，世界角速度为 Ω=0.5*vorticity+λ cross(p,E@p)。
独立数学测试确认 cross(Ω,p)=W@p+λ(E@p−(p.T@E@p)p)，还以四元数增量的中心差分核查 p_dot。
位置继续使用旧位置速度的显式 Euler，V=u_inf；姿态用旧 Ω 的有限旋转增量，新位置只刷新 V/Ω。
没有质量/惯性、旧 drift、其它力或额外动力学。

简单剪切 γ=20 s⁻¹ 只是 SYNTHETIC_VALIDATION_PARAMETER。
形状轴 Jeffery 周期 T_axis=π/γ(r+1/r)，对应 p 到 −p 后形状重复；有向 p 的完整周期为 2*T_axis。
测量周期由数值未折叠角第一次经过 −π 的相邻时间点插值得到，未把解析值直接当成测量值。

## 7. 自动测试

| 检查 | PASS/FAIL | 最大误差 / 范围 | test file |
|---|---|---|---|
| P0 原回归 | {v['particle0_regression']} | {tests['particle0']['passed']} passed，原源码/测试/数据/图 SHA 不变 | tests/particle0/ |
| P1 原回归 | {v['particle1_regression']} | {tests['particle1']['passed']} passed，原源码/测试/数据/图 SHA 不变 | tests/particle1/ |
| P2 总计 | {status} | {tests['particle2']['passed']} passed / {tests['particle2']['failed']} failed / {tests['particle2']['skipped']} skipped | tests/particle2/ |
| Frozen handoff | {v['frozen_handoff']} | {tests['frozen']['passed']} passed；原冻结分支执行 | 原 frozen checkout tests/ |
| Frozen integrity | {v['frozen_integrity']} | validate_frozen + verify_manifest；4602 payload 文件、31 科学 SHA | frozen 脚本及日志 |
| 合同与来源 | {status} | PREPRINT、参数、guard 和 CV role 核验 | test_rbc_distribution_contract.py |
| RNG/guard/复现/拒绝 | {status} | 独立进程逐位一致；独立 CDF smoke | test_rbc_distribution_reproducibility.py、test_rbc_distribution_guards.py、test_rbc_shape_rejection.py |
| 几何、SI、分布来源 | {status} | 体积相对误差 {v['accepted_geometry_checks']['max_volume_reconstruction_relative_error']:.3e} | test_rbc_volume_reconstruction.py、test_rbc_si_units.py、test_rbc_geometry.py、test_no_average_rbc_shortcut.py |
| 四元数约定/不可变/等价 | {status} | body→world、左乘、方向、q/−q、Ω=0 | test_rbc_state.py、test_quaternion_convention.py、test_quaternion_sign_equivalence.py |
| q 和 p 单位长度 | {status} | max q {v['max_quaternion_norm_error']:.3e}；max p {v['max_p_norm_error']:.3e} | test_quaternion_norm.py、test_real_rbc_quaternion_norm.py |
| Jeffery 独立恒等式 | {status} | 全部场最大绝对分量误差 {v['max_jeffery_identity_error']:.3e} s⁻¹ | test_gradient_decomposition.py、test_jeffery_angular_velocity_identity.py |
| 无梯度静止姿态 | {v['static_flow_test']} | 旋转矩阵误差 {static['max_rotation_matrix_error']:.1e} | test_rbc_static_orientation.py、test_rbc_translation_follows_flow.py |
| 整体旋转 | {v['rigid_rotation_test']} | 轴分量误差 {rotation['max_axis_vector_error']:.3e}；Ω 误差 0 | test_rigid_rotation_orientation.py |
| 五几何/三 dt 剪切周期 | {v['simple_shear_test']} | 最大周期相对误差 {v['max_jeffery_period_relative_error']:.3e} | test_simple_shear_jeffery_period.py |
| 64 几何覆盖分布 | {v['distribution_wide_orientation_test']} | 最大周期相对误差 {v['distribution_wide_max_period_relative_error']:.3e} | test_distribution_wide_jeffery.py |
| 真实中心/V=u | {status} | {v['real_row_count']} 行 / {v['real_checked_segment_count']} 段；V−u=0 | test_real_rbc_center_finite.py、test_real_rbc_velocity_matches_fem.py |
| 真实不同形状/姿态 | {status} | 5×3 次；全部有限；每行公式与每步增量核查 | test_real_rbc_distribution_geometries.py、test_real_rbc_orientation_finite.py |
| 时间步对比与 11 图 | {status} | acos(abs(dot))，源数据 SHA，PNG 分辨率，重新作图 | test_particle2_timestep_comparison.py、test_particle2_review_artifacts.py、test_particle2_plot_generation.py |

以上 P2 文件均在 [tests/particle2](../../tests/particle2/)。完整命令见 [final_check_commands.json](logs/final_check_commands.json)，
各 suite 日志及 JUnit 见 [TEST_SUMMARY.md](logs/TEST_SUMMARY.md)。所有最终 suite 无失败、无跳过。
norm 预算在运行前固定为 512eps={NORM_ATOL:.6e}；合成恒等式预算为 1024eps×γ；
γdt=1/512、1/1024、1/2048 的周期相对误差预算为 4γdt，未根据结果放宽。
整体旋转的舍入界为 {rotation['axis_error_bound']:.6e}。
真实恒等式误差除以 max(1,||G||_F) 后最大 {v['real_max_identity_error_scaled_by_max_one_gradient_frobenius']:.3e}。

开发中修复了 Ω=0 时非单位输入 quaternion 的 API 归一化，并加永久回归测试；重新生成受影响静止流证据。
图 08 调整标题后自动 PNG 分辨率检查曾为 1 failed / 73 passed，已扩大画布并保留失败日志；
最终重新检查通过，没有改测试阈值或科学参数。还修复了字体缺字和图例重叠等显示问题。

真实起点沿用 P1 的入口相邻 tetra 选择规则，在 {initial['candidate_count']} 个候选中选 canonical tetra {initial['initial_tetra_id']}。
同一中心起点为 {initial['initial_position_m']} m；五个样本取 r 的 5/25/50/75/95% 附近原始样本。
64 样本按 r 分成等数量区间，seed=2026092064 每层选择一个。
真实三个 VALIDATION_ONLY dt (s)：{initial['validation_timesteps_s']}，与 P1 已验证值相同，未为改善结果试调。
每段检查首次表面交点；出口处按实际分段 dt 更新姿态，保存终点场并停止，不投影、不反射。
15 次均从 OUTLET_02 离开，同 dt 不同形状的中心线相同。三个出口时间为
{', '.join(f'{v["real_fem_cases"][i]["event"]["time_s"]:.12f}' for i in range(3))} s。

真实姿态与最细重放的比较如下，单位为度。采用共同时间、最细 p 的线性插值后归一化与 acos(abs(dot))；
其中还包含中心路径不同导致的采样梯度不同，最细重放不是解析真解。

| r | 粗/最细 max (°) | 粗/最细 median (°) | 中/最细 max (°) | 中/最细 median (°) |
|---|---:|---:|---:|---:|
{comparison_table}

五个最大差均下降；最薄样本 median 略增，其余下降，不能写成所有误差都单调。
**最大差仍可达 {v['real_dt_max_difference_middle_deg']:.3f}°，真实姿态时间步敏感，未证明生产收敛。**
线性 tetra 中 G 是分片常数，跨 tetra 可跳变；保留原 G/vorticity/Ω，不进行平滑，也没有修改 frozen FEM。
自动 PASS 表示满足本阶段规定的算法、有限性、单位长度、合成解析及只读要求，不表示已满足生产姿态精度。

## 8. 人工审核图

每张图的输入 CSV/JSON 及 SHA 在 data/00_figure_sources.json 至 data/10_figure_sources.json。
剪切时间序列为绘图稀疏保存，误差最大值计算覆盖每一步；真实轨迹保存每一步。

{figures}

## 9. 当前限制

- diameter 来源是 C57BL/6Case PREPRINT；不同 C57BL/6 亚系没有被当成相同群体。
- D 正态是 V0 approximation，volume SD 来自 MCV 与所选 RDW-CV 的模型组合，不是直接单细胞联合测量。
- D/V independence 是 provisional assumption；guard 是模型范围，accepted 统计是模型推导结果。
- RBC 是 rigid spheroid，非真实 biconcave membrane；没有 deformation、膜节点、弹簧或弯曲能。
- 没有 wall gap、wall force、contact、lubrication、reflection 或 projection。FINITE_SIZE_WALL_CLEARANCE=NOT_VALIDATED_PARTICLE3。
- 没有 RBC-RBC、RBC-MB、多体阻力、Brownian、lift、重力或旧 reduced-order RBC drift。
- 没有 LAMMPS，没有 hematocrit injection，没有 production particle population。
- 没有 production timestep，也没有 production inlet orientation distribution；真实 dt 敏感性仍明显。
- 毛细血管中真实 RBC 会变形，本 V0 不能描述；没有通过缩小 RBC 或修改 distribution 来隐藏此限制。

## 10. 人工审核状态

AUTOMATED_CHECKS = {status}

MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW

STAGE_RESULT = AUTOMATED_PASS_PENDING_USER_REVIEW

FINITE_SIZE_WALL_CLEARANCE = NOT_VALIDATED_PARTICLE3

PRODUCTION_PARTICLE_TIMESTEP_FROZEN = false

PRODUCTION_ORIENTATION_DISTRIBUTION_FROZEN = false

NO_CFD_EXECUTED = true

PARTICLE3_STARTED = false

Particle-2 在此停止。只做本地提交，没有 push 或 merge main。
全部变更见 [FILES_CHANGED.txt](FILES_CHANGED.txt)，复现方法见 [PARTICLE2_README.md](../../PARTICLE2_README.md)，
机器验证记录见 [PARTICLE2_VALIDATION.json](PARTICLE2_VALIDATION.json)。人工审核需由用户完成。
"""
    (REPORT/"PARTICLE2_REVIEW.md").write_text(content,encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handoff-root",type=Path,required=True)
    args=parser.parse_args();handoff=args.handoff_root.resolve()
    source_paths=["particle_3d/src","particle_3d/tests","particle_3d/scripts","particle_3d/contracts",
                  "particle_3d/PARTICLE2_README.md","particle_3d/pyproject.toml","particle_3d/reports/particle2/literature"]
    if git("status","--porcelain","--",*source_paths):
        raise RuntimeError("Commit implementation first so git_commit identifies the exact tested sources")
    if git("branch","--show-current")!=P2_BRANCH:
        raise RuntimeError("Dedicated Particle-2 branch required")
    if (git("branch","--show-current",cwd=handoff)!=FEM_BRANCH or git("rev-parse","HEAD",cwd=handoff)!=FROZEN_COMMIT
            or git("status","--porcelain",cwd=handoff)):
        raise RuntimeError("Handoff requires the original unchanged frozen checkout; do not repair it here")
    dependencies=check_dependencies(REPO)
    c=load_contract();fem=REPO/"formal_3D_flow_solver/FEM_SimVascular"
    frozen,_,_,_=read_frozen(fem)
    commands=[];LOGS.mkdir(parents=True,exist_ok=True)
    for name,script in [("frozen_validate","validate_frozen.py"),("frozen_manifest","verify_manifest.py")]:
        commands.append(run(name,[sys.executable,"-B","scripts/fem_freeze_sync/"+script],fem))
        if commands[-1]["returncode"]:
            write_json(LOGS/"final_check_commands.json",commands)
            raise RuntimeError("Frozen integrity failure; preserve all original inputs and inspect log")
    frozen_paths=["tests/test_fem_sync_manifest.py"]+[
        str(p.relative_to(handoff)) for p in sorted((handoff/"tests").glob("test_frozen_*.py"))
    ]+["tests/test_particle_handoff_contract.py","tests/test_no_particle_implementation_yet.py","tests/test_main_history_preserved.py"]
    suites=[("frozen",handoff,frozen_paths)]+[(f"particle{i}",REPO,[f"particle_3d/tests/particle{i}"]) for i in range(3)]
    tests={}
    for name,cwd,paths in suites:
        xml=LOGS/f"{name}_junit.xml"
        commands.append(run(name+"_pytest",[sys.executable,"-B","-m","pytest","-q","-p","no:cacheprovider",
                                           f"--junitxml={xml}",*paths],cwd))
        tests[name]=junit(xml)
    write_json(LOGS/"final_check_commands.json",commands)
    check_dependencies(REPO)
    frozen_diff=git("diff",FROZEN_COMMIT,"--","formal_3D_flow_solver/FEM_SimVascular")
    static=read_data("04_static_metrics.json");rotation=read_data("05_rotation_metrics.json")
    shear=[m for d in range(3) for m in read_data(f"06_shear_dt{d}_metrics.json")]
    wide=read_data("07_distribution_wide_metrics.json");real=read_data("08_real_sweep.json")
    meta=read_data("C57BL6_RBC_GEOMETRY_VALIDATION_100000.metadata.json")
    stats=read_data("accepted_statistics.json");repro=read_data("reproducibility.json")
    cdf=read_data("01_distribution_cdf_checks.json");comparison=read_data("10_axis_comparison_summary.json")
    evidence=numerical_evidence(shear,wide,real)
    figures=sorted(str(p.relative_to(REPO)) for p in (REPORT/"figures").glob("*.png"))
    passed=(all(cmd["returncode"]==0 for cmd in commands) and not frozen_diff and
            all(t["failed"]==t["skipped"]==0 for t in tests.values()) and
            tests["particle0"]["passed"]==46 and tests["particle1"]["passed"]==67 and
            tests["frozen"]["passed"]==18 and tests["particle2"]["passed"]>0 and len(figures)==11 and
            static["passed"] and rotation["passed"] and len(shear)==15 and len(wide)==64 and len(real)==15 and
            all(m["passed"] for m in shear+wide+real+cdf) and repro["status"]=="PASS" and
            evidence["max_quaternion_norm_error"]<=NORM_ATOL and evidence["max_p_norm_error"]<=NORM_ATOL and
            evidence["real_max_identity_error_scaled_by_max_one_gradient_frobenius"]<=IDENTITY_RELATIVE_BUDGET)
    acceptance=git("log","-1","--format=%H","--","particle_3d/reports/particle1/PARTICLE1_VALIDATION.json")
    sources=sorted(p for top in [PACKAGE/"src",PACKAGE/"tests",PACKAGE/"scripts"] for p in top.rglob("*.py"))
    sources.extend([PACKAGE/"PARTICLE2_README.md",Path(DEFAULT_CONTRACT),REPORT/"literature/sources.json"])
    v=dict(stage="Particle-2",automated_checks="PASS" if passed else "FAIL",
           stage_result="AUTOMATED_PASS_PENDING_USER_REVIEW" if passed else "AUTOMATED_FAIL",
           git_branch=git("branch","--show-current"),git_commit=git("rev-parse","HEAD"),
           git_commit_semantics="Exact tested committed implementation; later evidence-only commit may save this record.",
           particle0_dependency_commit=P0_COMMIT,particle1_dependency_commit=P1_COMMIT,
           particle1_acceptance_evidence_commit=acceptance,
           particle1_manual_review_evidence=dependencies[1]["manual_review_evidence"],
           particle0_regression="PASS" if tests["particle0"]["failed"]==0 else "FAIL",
           particle1_regression="PASS" if tests["particle1"]["failed"]==0 else "FAIL",tests=tests,
           particle0_source_tests_figures_data_unchanged=True,particle1_source_tests_figures_data_unchanged=True,
           frozen_integrity="PASS" if not frozen_diff and all(cmd["returncode"]==0 for cmd in commands[:2]) else "FAIL",
           frozen_handoff="PASS" if tests["frozen"]["passed"]==18 and tests["frozen"]["failed"]==0 else "FAIL",
           frozen_base_commit=FROZEN_COMMIT,frozen_branch=FEM_BRANCH,frozen_payload_modified=bool(frozen_diff),
           mesh_sha256=frozen["mesh_sha256"],flow_sha256=frozen["flow_sha256"],
           distribution_contract=c["contract_name"],distribution_version=c["version"],
           distribution_contract_path=str(Path(DEFAULT_CONTRACT).relative_to(REPO)),distribution_contract_sha256=sha256(DEFAULT_CONTRACT),
           literature_sources=json.loads((REPORT/"literature/sources.json").read_text()),
           diameter_model="Normal approximation to reported C57BL/6Case single-cell diameter statistics",
           diameter_source_mean_um=c["diameter"]["mean_um"],diameter_source_sd_um=c["diameter"]["sd_um"],
           diameter_source_median_um=c["diameter"]["source_median_um"],diameter_source_iqr_um=c["diameter"]["source_iqr_um"],
           volume_model="Normal V0 approximation; MCV anchor plus user-selected CV",
           volume_mean_fL=c["volume"]["mean_fL"],volume_cv=c["volume"]["cv"],volume_sd_fL=c["volume"]["sd_fL"],
           volume_cv_role=c["volume"]["cv_role"],diameter_guard_um=c["diameter"]["guard_um"],volume_guard_fL=c["volume"]["guard_fL"],
           guard_role=c["diameter"]["guard_role"],D_V_dependence=c["D_V_dependence"],D_V_dependence_role=c["D_V_dependence_role"],
           validation_population_N=meta["N"],validation_population_seed=meta["seed"],validation_population_role=meta["role"],
           formal_hematocrit_population=False,production_particle_population=False,
           candidate_count=meta["candidate_count"],generated_candidate_count=meta["generated_candidate_count"],
           unused_generated_tail_count=meta["unused_generated_tail_count"],candidate_count_semantics=meta["count_semantics"],
           diameter_guard_rejections=meta["diameter_guard_rejections"],volume_guard_rejections=meta["volume_guard_rejections"],
           guard_rejection_count=meta["guard_rejection_count"],shape_rejections=meta["shape_rejection_count"],acceptance_rate=meta["acceptance_rate"],
           accepted_D_statistics=stats["D_um"],accepted_V_statistics=stats["V_fL"],accepted_c_statistics=stats["c_um"],
           accepted_full_thickness_statistics=stats["full_thickness_um"],accepted_r_statistics=stats["r"],
           accepted_lambda_statistics=stats["jeffery_lambda"],accepted_D_V_sample_correlation=stats["accepted_D_V_sample_correlation"],
           latent_and_guard_only_cdf_checks=cdf,reproducibility_status=repro["status"],reproducibility=repro,
           population_sha256=meta["csv_sha256"],population_metadata=meta,
           static_flow_test="PASS" if static["passed"] else "FAIL",rigid_rotation_test="PASS" if rotation["passed"] else "FAIL",
           simple_shear_test="PASS" if all(m["passed"] for m in shear) else "FAIL",
           distribution_wide_orientation_test="PASS" if all(m["passed"] for m in wide) else "FAIL",
           synthetic_validation_budgets=read_data("synthetic_validation_budgets.json"),
           static_metrics=static,rigid_rotation_metrics=rotation,simple_shear_metrics=shear,distribution_wide_metrics=wide,
           real_fem_geometry_count=len({m["rbc_id"] for m in real}),real_fem_replay_count=len(real),
           real_fem_center_finite=all(m["center_and_orientation_finite"] and m["all_centers_inside"] for m in real),
           real_fem_orientation_finite=all(m["center_and_orientation_finite"] for m in real),
           real_fem_outlet_events=[dict(geometry_index=m["geometry_index"],dt_index=m["dt_index"],rbc_id=m["rbc_id"],**m["event"]) for m in real],
           real_fem_cases=real,real_initialization=read_data("08_real_initialization.json"),
           real_row_count=sum(m["row_count"] for m in real),real_checked_segment_count=sum(m["checked_segment_count"] for m in real),
           real_max_velocity_relation_error_m_s=max(m["max_velocity_relation_error_m_s"] for m in real),
           real_center_wall_crossing=any(m["wall_crossing"] for m in real),real_center_inlet_crossing=any(m["inlet_crossing"] for m in real),
           validation_timestep_comparison=comparison,
           real_dt_max_difference_coarse_deg=float(np.degrees(max(m["max_axis_difference_rad"] for m in comparison if m["dt_index"]==0))),
           real_dt_max_difference_middle_deg=float(np.degrees(max(m["max_axis_difference_rad"] for m in comparison if m["dt_index"]==1))),
           real_orientation_timestep_convergence="NOT_ESTABLISHED; significant local sensitivity remains; five max differences decrease, thinnest median slightly increases",
           real_gradient_smoothing=False,finite_size_wall_clearance="NOT_VALIDATED_PARTICLE3",
           production_particle_timestep_frozen=False,production_orientation_distribution_frozen=False,no_rbc_deformation=True,
           manual_visual_review="PENDING_USER_REVIEW",visual_step_jump_review="PENDING_USER_REVIEW",particle3_started=False,
           no_cfd_executed=True,no_gpu_server_started=True,pushed=False,merged_main=False,
           figures=figures,test_files=sorted(str(p.relative_to(REPO)) for p in (PACKAGE/"tests/particle2").glob("test_*.py")),
           source_sha256={str(p.relative_to(REPO)):sha256(p) for p in sources},
           figure_sha256={p:sha256(REPO/p) for p in figures},
           data_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted(DATA.iterdir()) if p.is_file()},
           log_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted(LOGS.iterdir()) if p.is_file() and p.suffix in [".log",".xml",".json"]},
           recorded_at_utc=datetime.now(timezone.utc).isoformat(),
           environment=dict(python_version=platform.python_version(),executable=sys.executable,platform=platform.platform(),
                            packages={n:version(n) for n in ["numpy","vtk","pyvista","matplotlib","pytest","Pillow"]}),**evidence)
    write_json(REPORT/"PARTICLE2_VALIDATION.json",v)
    review(v,stats,static,rotation)
    summaries="\n\n".join(f"{name}: {t['passed']} passed / {t['failed']} failed / {t['skipped']} skipped; [{name}_pytest.log]({name}_pytest.log)" for name,t in tests.items())
    (LOGS/"TEST_SUMMARY.md").write_text(f"# Particle-2 最终检查\n\nAUTOMATED_CHECKS = {v['automated_checks']}\n\n{summaries}\n\n"
        f"Frozen integrity = {v['frozen_integrity']}；P0/P1 原源码、测试、图及数据 SHA 一致。\n\n"
        "开发图像分辨率失败日志保留，最终 suite 为验收依据。真实姿态时间步敏感，未证明 production 收敛。\n\n"
        "MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW；Particle-3 未启动，CFD 未运行。\n")
    changed=set(git("diff","--name-only",acceptance).splitlines())
    changed.update(git("ls-files","--others","--exclude-standard","particle_3d").splitlines())
    changed.add("particle_3d/reports/particle2/FILES_CHANGED.txt")
    (REPORT/"FILES_CHANGED.txt").write_text("\n".join(sorted(changed))+"\n")
    print(json.dumps({k:v[k] for k in ["automated_checks","git_commit","tests","max_quaternion_norm_error",
                                     "max_jeffery_identity_error","max_jeffery_period_relative_error","manual_visual_review"]}),flush=True)
    if not passed:raise SystemExit(1)


if __name__=="__main__":
    main()
