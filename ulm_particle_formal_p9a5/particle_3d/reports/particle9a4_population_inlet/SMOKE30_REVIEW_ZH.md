# 首 30 个 accepted births smoke：PASS WITH SUPPORTED STATIONARY

在 inlet gate、115 portable、35 new 和15 legacy tests 全部通过后，取100k群体的前30个 accepted births；不按出口筛选。cohort SHA `aebee7de17212e461b50b8a47f20ee79452161d4d3253268bfc3328316af06d2`。服务器独立目录 `/workspace/particle9a4_population_inlet_20260924T233708Z`，Python `/root/particle8_2_runs/env/bin/python`，6 CPU workers、OMP/BLAS/MKL/NUMEXPR=1。没有 GPU 轨迹内核。

结果 completed=28、stationary=2（ID 4、8）、solver fail=0、time-limit=0。MB 出口 O1/O2/O3=0/6/22。两条 stationary 的最后32个已接受求解均有 contact rank=3、平移速度在原 KKT budget 内、非负 multipliers，因此按原审核规则是受约束 stationary；不把其解释成生理 trapping，也没有更改接触模型以消除 stationary。

保存样本 13479，连续 handoff 证书叶节点 14592。原诊断及独立真实 WALL 最近距离复算均给 penetration=0、handoff violation=0。逐线段入口向外相交复查给 escape=0。NaN/Inf=0，出口分类复判全部一致。

最小真实 wall gap=1.9999999999811979e-09 m；最小 g_nf=-1.8802228913636596e-20 m。后者微小负值约 −1.88e−20 m，小于冻结 roundoff 2.0094114631e−17 m，报告原始值而非裁为0；因此无超容差违规。连续证书验证贯穿原 accepted path；额外独立最近距离检查覆盖保存状态。

point tracer 在同一已接受的30个初始中心事后运行：O1=1、O2=12、O3=14、noexit=3（ID18、25、28达到30 s point horizon）。这不是100k accepted群体的精确 basin split，更不能以 fluid O1/O2/O3 比例作 quota。finite-size MB 和 point 的路径模型不同。

MB运行 140.198 s；保存数据审核 12.362 s。部署204个输入/源码的 SHA，回传154个输出逐一 SHA 验证。P9-A.1、积分入口、初始速度规则、dt=0.00025 s、MB horizon=1.5 s 均复用原 Network runner；新输入只有入口群体/身份。原 integrate_admitted 对既定 birth 还有启动一致性验证，不重抽、不重记 source acceptance、不改变 ledger。

[逐轨迹指标](data/smoke30_metrics.json) · [接触与连续证书证据](data/smoke30_contact_audit.json) · [汇总](data/smoke30_summary.json) · [图08](figures/Figure_08_smoke30_trajectories.png) · [服务器日志](logs/server_smoke30.txt)

结论为可供人工 large-sample review；本轮没有运行正式500。
