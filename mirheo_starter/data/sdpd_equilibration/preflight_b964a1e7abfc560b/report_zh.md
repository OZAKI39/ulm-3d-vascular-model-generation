# 无驱动 SDPD 稳定性实验

状态：PREPARATION_CPU_VALIDATION_REQUIRED。执行：NOT_RUN_AWAITING_AUTHORIZATION。
程序完成、液体稳定、稳定温度正确、全部物性合格分别记录；selection=null，人工验收 PENDING。

固定 400000 步，dt*=1e-06，连续 12.662309354 µs；正式窗口 t*∈(0.30,0.40]。
0.30* 是事前设计的观察起点，不是已证明的平衡时间。启动峰值完整保留，不按温度偏离 2% 停机。
原参数 μ*=3312.4601269222444，盒 [8,8,8]，N=4096，n*=8，m*=1，kBT*=1，rc*=1，Linear EOS，cs*=120，rho0*=0。
原单位 L0=5e-07 m，M0=1.6499999999999997e-17 kg，t0=3.1655773384195756e-05 s；目标 298.15 K。
初始化只做一次 Uniform、高斯速度及 COM 扣除；保留原生压力、随机和耗散作用，无驱动力、无追加恒温器、无持续速度缩放。

## 预算与实测成本来源

旧无驱动任务：172000 步，实收 210.933464 s。
等比例成本 490.543 s；加 15% 余量及 18 s 初始化/输出/退出余量，预计 582.124 s。
建议单任务上限 600 s，包含 12 s 合作退出预留；预计正常运行必须在 588 s 前完成。
原授权剩余 107.080671112 s；本任务可用 107.080671112 s；需要追加 492.919328888 s，建议申请整数 493 s。
新观测增加少量应力读取、粒子 ID 与固定快照，尚无新实测速度；余量不是耗时或获得稳定平台的保证。预算不足则拒绝启动，不缩短主实验。

## 事前判据

先看整个固定窗口的两半与四个固定块趋势，再估计相关时间。明显持续漂移时不报告可信稳态 ACF。
块长至少 max(5×新窗口及两半的 ACF，0.00783050120411*)；至少 8 个完整块，两半各至少 4 块。
均值与 CI 使用同一批完整块，尾部单列。每块均值用点显示，整条轨迹均值的 CI 用带显示，不伪造单块独立 CI。
温度两半均值漂移加误差上界不超过目标热能的 0.5%；固定四分窗范围不超过 1%。温度初筛：|均值/目标−1|+CI半宽/目标≤2%。
压力漂移尺度为 n*kBT=8.0*（0.263449920 Pa），不以背景压力作分母。
核密度均值以原 n*、方差以 n*²、分位数以 n* 归一化；COM 漂移以 sqrt(kBT/m) 归一化。精确阈值见 run_plan.json。
这些是本次事前筛选设计，不是普适物理定律。压力稳定不等于出口压力覆盖，N*m/V 恒定不证明局部结构稳定。

## 本次证据

NOT_RUN_AWAITING_AUTHORIZATION
液体稳定：None；稳定温度匹配：None。
实际步数：None；温度均值*：None；温度 CI95*：None。
未运行时所有新测量字段为 null；raw_statistics.csv 只有表头。历史曲线单独标注，不能填成新结果。
未新增黏度或 EOS 验证，不能据此宣称全部液体物性合格。

## 唯一下一步

明确追加授权后执行这一项完整无驱动实验。

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration.yaml --execute
```

## 可追溯与保护

诊断来源由交付记录及 manifest 哈希定位；本次计划哈希 `961a4fda2dc853d62b9d52360806b0f85749deb1dfd28e6fe1b06bbaf19bccef`。
复用原 runner 的单 GPU 锁、monotonic 看门狗与本进程组清理。旧用量、新用量、未完成预约统一核算；原授权记录和历史账本不覆盖。
只在用户明确给出秒数和范围后追加授权记录；配置 extra_authorized_gpu_seconds=0 不能自动放行。
普通快照不含完整原生随机/积分状态，不称为 checkpoint，也不拼接旧短轨迹。
压力定义：sum(stresses.xx + stresses.yy + stresses.zz)/(3V) + m*sum(|v-COM|^2)/(3V); compression positive, no second EOS term. Native Stress layout [xx,xy,xz,yy,yz,zz]. Force/density are pre-kick (step-1)*dt; momenta are post-kick step*dt. This is the existing discrete virial definition, with the phase difference explicitly retained.
数据包：/home/lzy/projects/mirheo_starter/data/sdpd_equilibration/preflight_b964a1e7abfc560b
