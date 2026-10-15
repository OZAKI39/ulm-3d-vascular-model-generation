# 锁定源码及积分时间层审查

本轮依据上轮 source_owner_fix 的实际源码与 SHA，未修改原生断裂算法。

`bond_bpm_peri.cpp:267–276` 先按当前存活键计算 theta；`342–345` 在随后的力循环中执行 stretch > s00 的断键并跳过该键。其余键仍使用该次 force evaluation 之前的 theta。alpha=0，临界伸长不依赖前一步 smin。

`verlet.cpp:255–259` 执行初始半次速度更新及位移更新；`306` 清空 force；`315–325` 先 pair 后 bond；`347` 执行 post_force；`348` 执行 final_integrate；`349` 执行 end_of_step。因此 post_force 时自由 DOF 的 v 是半步速度，end_of_step 才是接受整步速度。run 0 的 setup 会重新评价力，不能在连续动力学事件之间额外插入 run 0 来消除 lagged theta。

`fix_nve.cpp` 使用 dtf=0.5*dt*ftm2v。SI 单位 ftm2v=1。对自由 DOF：

```
v_half[n] = v[n-1] + dt/(2 mass) F[n-1]
x[n] = x[n-1] + dt v_half[n]
v[n] = v_half[n] + dt/(2 mass) F[n]
mass (v[n]-v[n-1]) = dt/2 (F[n-1]+F[n])
```

一次断键的 DeltaF_timelayer[n] 在 n 步接受态首先贡献 dt/2 的代数权重；若继续推进，下一步 initial_integrate 再使用一次旧 F[n]，另贡献 dt/2。到接受态 N 的累计代理必须使用梯形权重：初始和最后力各半权重，内部力全权重。代理评价于实际 lagged 轨迹，不能等同于另一个算法的反事实轨迹差。

原生 BPM/PERI atom->f 是体积力密度。独立 nodal force 为 N，必须且只能乘本点 V 一次。历史静态诊断只在读出时乘 V，因为没有积分器；任何经授权的动力学适配则必须在实际半次积分前完成既有 SI 力转换。此处不是新增材料模型。

本体内力、外部控制反力、阻尼、接触必须单列。历史小试样 pair zero、无阻尼、无背景拖曳；不能从全 V3 自动继承那些项。规定 DOF 不计入自由 DOF 动能；边界功按传入可变形系统的功定义，外部夹具惯性单独说明。

独立势能沿用 macro_lab_v1.Energy 及历史 Python 的展开式：

```
U = 0.5 V (K-5G/3) sum(theta_i^2)
  + 7.5 G V^2 sum_alive_unique(h_ij e_ij^2/r0_ij (1/m_i+1/m_j))
```

其中 m 固定为完整参考邻域，theta 随实际拓扑计算。native PE*V 具有已知偏能统计缺失，不用于能量余额。

带符号拓扑变化定义为 DeltaU_topology=U(x[n],B[n])-U(x[n],B[n-1])。若采用传入系统的边界功 W，则能量残差可定义为 `R=K_free+U-(K_free[0]+U[0])-W-DampingWork-sum(DeltaU_topology)`，其中 DampingWork 带功的符号（阻尼通常为负）。该残差不是 Gc。历史报告有用 `U_pre-U_post` 的相反号约定，本轮输出必须明确转换。

历史五键事件全部是直接写坐标的评价，没有可重建的动态速度、夹具功或真实物理首断时刻。其 nominal dt=0.5 us 不能赋予静态操作真实加载时间。
