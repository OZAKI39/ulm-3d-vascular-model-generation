# Particle-2：C57BL/6 RBC 分布与单刚性 RBC 姿态

本阶段建立 C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0，然后使用其中的真实抽样几何，
验证单个刚性轴对称扁球的平移和 Jeffery 姿态。100,000 个几何仅作分布验证；
五种形状和 64 种形状的重放都是相互独立的单 RBC 算例，没有多粒子系统或注入。

## 分布层

合同位于 contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json。
文献、来源性质和 V0 模型假设见 reports/particle2/00_distribution_literature_report.md。
D~Normal(6.79,0.93²) µm，V~Normal(47.9,7.0892²) fL，先独立成对抽样，
再依次检查 D guard、V guard、c<a。超出范围或过厚的原始组合直接拒绝，不裁切。
三个标准差范围是模型保护范围，不是实验 min/max。

    from particle_3d.rbc_distribution import sample_rbc_geometries
    from particle_3d.rbc import RBCGeometry
    population = sample_rbc_geometries(n, seed)
    geometry = RBCGeometry.from_population(population, index)

seed 必须显式提供。每个 chunk 固定 4096 行，每行同时生成 D/V。
candidate_count 是处理到第 N 个接受样本的前缀数；最后一个 chunk 未处理尾部另记。
拒绝计数按 D、V、shape 顺序互斥，accepted+rejections=candidate_count。
输出记录环境、种子、contract SHA、原候选编号和 float64 little-endian 数组 SHA。
CSV 是 UTF-8/LF、17 位有效数字，相邻 metadata 绑定 CSV SHA。
只承诺相同实现/N/seed/环境下复现，不承诺跨 NumPy 版本逐位一致。
统计将 latent、guard-only 和最终 accepted 三层区分开。

## 状态和运动接口

    from particle_3d.rbc import RBCState
    from particle_3d.rbc_integrator import advance_single_rbc
    new_state = advance_single_rbc(state, fem_field, dt_s)

核心 SI：轴长和位置 m、体积 m³、速度 m/s、梯度与角速度 s⁻¹。
RBCGeometry 和 RBCState 不可变，状态数组使用只读独立内存。
只有明确标记 SYNTHETIC_ONLY 的数学单元测试可以手填人工轴长。

四元数固定 q=(w,x,y,z)，R(q) 把 body 转到 world，body 短轴 e3=(0,0,1)，
世界短轴 p=R(q)e3。q 与 −q 是同一个旋转；p 与 −p 是同一个轴对称形状姿态。
两种等价关系不同，不能用原始 quaternion 分量相减判定姿态误差。

从 P0 的 G 得到 E=(G+G.T)/2、W=(G−G.T)/2，r=c/a<1，
λ=(r²−1)/(r²+1)<0。世界角速度：
Ω=0.5*vorticity+λ cross(p,E@p)。
独立测试确认 cross(Ω,p)=W@p+λ(E@p−(p.T@E@p)p)。

位置是旧位置速度的显式 Euler；姿态是旧 Ω 的有限旋转增量，
q_new=normalize(delta_q_world ⊗ q_old)，增量在左。
位置和姿态更新后只刷新新状态的 V/Ω，不增加积分阶段。
dt_s 必须由调用者提供，有限且为正；没有默认 production 时间步。

## 验证约定

分布种子 2026092002；64 个等数量 r 分层内的选择种子 2026092064。
五个代表几何取 r 的 5/25/50/75/95% 最近原始秩样本，不制造平均 RBC。
静止梯度和整体旋转使用固定非特殊短轴 (1,2,3)/sqrt(14)。
剪切测试在剪切平面从 p=(1,0,0) 开始，γ=20 s⁻¹ 仅为人工验证参数。
γdt=1/512、1/1024、1/2048，在运行前固定。

简单剪切的形状轴周期 T_axis=π/γ(r+1/r)；有向 p 整圈周期是 2*T_axis。
数值周期通过未折叠方向角第一次跨过 −π 的两个时间点插值得到，
不把解析周期直接赋值为测量结果。每步记录误差最大值；图源时间序列只等距抽取用于显示。
周期误差界预先取 4γdt（一次显式步进的一阶预算），norm 界为 512*float64 eps。

真实 FEM 沿用 P1 预先由几何/速度确定的三个验证步长，五种几何使用同一中心起点和初始短轴。
每条中心线段调用原 P1 边界分类器。首次出口交点处停止，终止姿态使用该线段实际经过的分段 dt。
无反射、推回或 wall response。只验证中心，不验证有限尺寸 RBC 的 wall gap。
原始 G/vorticity/Ω 均不平滑。真实 dt 比较使用 acos(abs(p1·p2))，
包含中心轨迹变化和离散梯度变化，最细结果只是参照，不是真解。

## WSL CPU 复现

从仓库根目录执行：

    P2_PY=/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python
    "$P2_PY" -B particle_3d/scripts/generate_rbc_population.py --n 100000 --seed 2026092002 --output particle_3d/reports/particle2/data
    "$P2_PY" -B particle_3d/scripts/verify_rbc_reproducibility.py --n 100000 --seed 2026092002 --output particle_3d/reports/particle2/data/reproducibility.json
    "$P2_PY" -B particle_3d/scripts/run_particle2_validation.py --stage synthetic
    "$P2_PY" -B particle_3d/scripts/run_particle2_validation.py --stage real
    "$P2_PY" -B particle_3d/scripts/generate_particle2_report.py
    "$P2_PY" -B -m pytest -q -p no:cacheprovider particle_3d/tests/particle0 particle_3d/tests/particle1 particle_3d/tests/particle2

作图支持 --stage 0 到 --stage 10。源码本地提交后再运行最终证据生成：

    "$P2_PY" -B particle_3d/scripts/finalize_particle2.py --handoff-root /home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular

测试覆盖旧 P0/P1、100,000 个分布样本、两独立进程、64 几何剪切重算、
15 次真实轨迹的每行位置/速度/梯度/角速度/四元数更新，以及原始中心边界分类。
原 P0/P1 数值、测试、图、数据和 FEM 均不改动。

## 阶段边界

rigid spheroid 不是双凹膜模型。没有变形、膜节点、弹簧、弯曲能、wall gap、接触、
润滑、投影、RBC-RBC、RBC-MB、多体阻力、旧 RBC drift、Brownian、lift、重力、
LAMMPS、hematocrit injection。毛细血管中真实 RBC 会变形，本 V0 不描述该现象。
production timestep 和入口 orientation distribution 都未冻结。
有限尺寸间隙为 NOT_VALIDATED_PARTICLE3，人工图审为 PENDING_USER_REVIEW。
Particle-2 完成后停止，不进入 Particle-3。
