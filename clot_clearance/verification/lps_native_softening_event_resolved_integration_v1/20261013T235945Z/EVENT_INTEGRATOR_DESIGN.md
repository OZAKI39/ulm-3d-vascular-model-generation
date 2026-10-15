# 事件分辨原生积分：设计、冻结与实际范围

冻结文件：EVENT_INTEGRATOR_FREEZE.json。只改变研究库中 fix_move::reset_dt；新增 event_state_audit.cpp 只读输出完整 special 邻接表。原生 bond_bpm_peri_variational、动态协调 fix、pair、fix_nve、Verlet 和自由能源码不变。编译器和锁定上游均沿用，不替换生产库。

每个基准区间长度 0.25 µs。真实原生候选步先完成第一半踢、位置漂移、解析夹持运动、POST_INTEGRATE κ/d 与拓扑提交、新力计算及 POST_FORCE 单位转换、第二半踢。Python读取实际候选Δd；超过0.10时整步拒绝。关闭试算实例，以接受态的完整原生restart重建实例和同名fix，再从恢复的x/v/F重新积分两段h/2。后半段使用前半段真正接受的新力和新历史。绝不剪裁d或将一次提交人为拆成多次。

采用只读的VV漂移与冻结损伤公式预测值≥0.075，或当前已处于子步时，才保存完整检查点。预测只控制I/O，不决定正式接受；实际Δd始终是唯一损伤门槛。预测坐标与候选坐标差须≤2e−18m；意外未保存检查点的越界只能停止。本次预测最大坐标误差为0。只有35个完整运行检查点，未在2755个接受步逐步写restart。

拒绝后核对所有点的x/v/F、材料力、κ/d、U/D/lastdrop/事件数/提交数/token、K/Uc、步号、物理时间、实际存活键端点及完整1–2/1–3/1–4排除表。4次完整运行回滚和4次局部回滚全部逐值一致。force-evaluation/cache诊断计数在重建实例时重新开始，它们不属于材料历史；全局trial编号和拒绝计数由驱动单独持久记录。Python账本只在候选接受之后累积，拒绝候选不进入W/D/Jc/R。

二分下限预注册为0.03125 µs；实际仅用到0.125 µs，没有测试/宣称更小步长。每个基准区间通过队列完整覆盖后才恢复0.25 µs；覆盖误差容差1e−15 s，仅处理浮点累计舍入。实际4个相邻基准区间二分，640.25–641.25 µs有8个真实半步；随后恢复基准步长。

锁定版本 Input::timestep 先update_time再改变dt并调用各fix的reset_dt；Update保存atime/atimestep，restart序列化ATIME/ATIMESTEP/TIMESTEP/NTIMESTEP。新fix_move允许六个完全规定的equal变量使用可变dt，其他样式仍拒绝。边界变量读取真实time；没有step×0.25µs换算。原生fix_nve的reset_dt同步两次半踢系数。

物理账本仍为 W=ΔU+ΔK+D+ΔUc−Jc+R。W逐接受步用真实夹持Δx的梯形反力功，加夹持动能变化；动量/VV等式使用实际h；D用同一新坐标的U_old−U_new。负D仅显示，不重复累计。力切换误差独立列出。归一化累计尺度与原轮完全相同，未借调分母掩盖误差。局部窗口仅为对比继承旧0–640.25µs标量账本；正式新路径从t=0重新积分，没有轨迹拼接。

局部真实restart来自最小2561步确定性重放；与旧640.25µs完整FP64逐值一致。局部候选也精确复现旧640.5µs现场。新完整路径在687.75µs之前通过逐步门槛；下一候选发生固定盒子丢点，停止而未扩大域、改材料、放宽门槛或继续计算。

本算法仅在观测到的4次损伤拒绝上认证完整回滚；不承诺任意原生异常可恢复。原生丢点异常暴露驱动finally仍调用完整gather/restart的缺陷，导致部分内存统计未落盘。未用NPZ伪造最终restart。此限制及丢点属于下一轮优先问题。

参考说明以本地锁定源码为准。官方文档用于核对一般接口：[NVE/Velocity-Verlet](https://docs.lammps.org/fix_nve.html)、[时间步命令](https://docs.lammps.org/timestep.html)、[规定运动](https://docs.lammps.org/fix_move.html)、[restart重建及同名fix](https://docs.lammps.org/read_restart.html)、[Fix接口](https://docs.lammps.org/Developer_write_fix.html)。标准[dt/reset](https://docs.lammps.org/fix_dt_reset.html)并不提供本项目的损伤增量事务准入，本轮没有采用它。
