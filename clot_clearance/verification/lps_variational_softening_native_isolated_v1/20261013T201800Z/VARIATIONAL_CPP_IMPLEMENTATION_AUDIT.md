# 隔离原生实现审查

仅新增 `src/BPM/bond_bpm_peri_variational.{h,cpp}`，以 `bpm/peri/variational` 独立选择。继承服务器锁定 `Bond` API，原最终双修复 `bpm/peri`、`bond_bpm.cpp`、`pair_bpm_peri.cpp` 不修改。官方接口说明仅辅助：[bond style 开发](https://docs.lammps.org/Developer_write_bond.html)、[restart](https://docs.lammps.org/read_restart.html)；实际依据是已冻结本机版本源码。

数学基线为上一轮正式 model.py（非共享内部伸长初稿）。端点 q=(i,b) 的 D=15GV²h/(r m_i)、L=3Vh/m_i、A=V(K−5G/3)、k=9KV²h/(r m_i)；den=Dd+k(1−d)，w=k(1−d)/den，h'=d/den，p_i=ΣLwe/(A⁻¹+ΣL²h')，v_q=we−h'Lp_i，t_q=w(De+Lp_i)。两端相加构成物理键力，沿当前真实方向成对组装。此为同一局部正定二次自由能的精确对角加秩一消元。6 个端点另用完整稠密线性系统核对，非近似优化。d=0、1 均直接进入有限公式，未加 epsilon 或伪逆。

完整 U=Σ_q[D_q v_q²/2+k_q(1−d)d((D_q e+Lp)/den)²/2]+Σ_i Aθ_el,i²/2，θ_el=ΣLv。实际节点力和该 U 同时计算并通过固定历史中心差分。N 与 J 是导出字段单位。原生 BPM 的 atom f 保留历史 N/m³ 约定，乘固定 V 得 N；本研究 style 的完整全局 bond energy 为 U/V，乘 V 得 J，独立接口直接给 U。原旧式 PE tally 不作为完整势能参考。所有路径为规定坐标，不涉及自由动力学质量/加速度认证。

历史：compute 只读 κ/d。显式 C API `lpsvar_action(op=2,accepted_token)` 在 run 段之间执行：汇集 owned 当前坐标→以旧 d 求 U_old→κ=max(κ,e/r,0)→冻结正式双线性 d(κ)→U_new→有符号 ΔD=U_old−U_new→累计 D→一次性失效拓扑移除。重复 token 且坐标相同为无操作；不同坐标或倒退 token 拒绝。run0、差分正负扰动、restart 首次 force 均不提交损伤。API op=1 仅允许在首次物理提交前导入明确标注的固定诊断历史，不能用于重置已发生的物理损伤。

身份策略：本轮不借助容易重排的本地 bond loop 下标，也不把 deleted bond history 丢弃。稳定图保存在 style 的全局原始 ID 表中，每 rank 保有完整小图，使用 owned 坐标集合与唯一 owner 断言，不读取 ghost 的 owned bond/history 数组；施力只加 owned atom，原生 pair 仍遵循原通信。支持原始固定非周期、同体积、单类型试样；重复端点及周期镜像明确拒绝。O(Nb) 本地状态代数，但每 rank 复制图与全局通信，绝非生产并行伸缩性方案。节点局部数组扩容无需搬迁此独立图，但未实际触发额外 nmax 扩容，不虚报该子项。

断键采用双修复版的 owned-only 原则：仅访问 owned atom 的 num_bond/bond_atom，移除对应 1-2 排斥屏蔽；其他 special 权重固定为 1。保留原始图用于弹性内部状态消元，完全断键实际 t=0。部分损伤保持原生 1-2 排斥屏蔽，防止重复接触力；完全失效后下一次正常 neighbor rebuild 启用现有 1.0dx 工程 pair 排斥。49 对真实闭合接触逐对对照原 pair 罚力，未连接壁面或 V3。

重要本构差异：旧二值删键法将失效伸长从 θ 求和中删除；新能量则释放相应内部约束并局部松弛。因此新模型 d∈{0,1} 一般不等于旧断键后力。已在真实历史五键断裂坐标上实测、单列差异；原生新实现准确匹配冻结候选，未为追求旧力而改变自由能。

一次代码安全修正：删除未使用的 lpsvar_sort 辅助入口，因为它直接调用 Atom::sort，绕过源码要求的 exchange/borders 生命周期，且缺少 sortfreq 前置条件。没有科学路径调用此入口；所有排序使用标准 atom_modify+run0。数学核心未变，修正前代码、二进制、数值结果和 diff 全保留。最终库重做必要门槛与原两条路径。首次 MPI 分区未触发迁移，已保存 NOT_EXERCISED；按原坐标选取 -40 µm 分区后实际迁移49个节点，未调材料或扫描。

范围：341节点原生运动学一致性、有限MPI/排序/restart/接触语义；未认证大规模MPI、通用周期图、多体积、多类型、动力学积分、宏观强度、自然裂纹Gc或实验血栓。下一轮仅可进入隔离的小型非均匀裂纹研究。
