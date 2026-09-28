# BraVa 最新结果与停止状态

本次 GitHub 同步只读取、保存已有结果，没有运行 CFD 或微泡积分。

## 已完成的流场

几何是用户选择的 BG001 RMCA BALANCED 四端口血管芯，含人工接管、无盒体。candidate 0 与15是同一几何的刚性打印姿态；本模型不含重力，复用同一物理解并旋转位置与矢量。入口为18,000 μL/min（18 mL/min），CFD隐式推进dt=0.01 s。

采用svMultiPhysics P1/P1稳定化牛顿流体模型，ρ=1056 kg/m³、μ=0.00345312 Pa·s，刚性无滑移壁面。网格135,852节点、692,313四面体，未设专门边界层。先进行等流量校准，再以三个压力出口重新求解；表内是正式压力出口求解的实际分流。

| 边界 | 流量 μL/min | 分流 % | 设定牵引压力 Pa |
|---|---:|---:|---:|
| OUTLET_01 | 6000.016253 | 33.333424 | 556.970656 |
| OUTLET_02 | 5998.129842 | 33.322944 | 0.000000 |
| OUTLET_03 | 6001.853905 | 33.343633 | 17.137427 |

正式第81步达到已设定稳态标准，独立检查通过。流场SHA256：`73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`。

- 实际算例：[balanced_pressure_final](cases/balanced_pressure_final/)。
- [candidate 0 流场动画](visualization/candidate_0/OPEN_RESULTS.html) / [candidate 15 流场动画](visualization/candidate_15/OPEN_RESULTS.html)。每姿态各5段，1920×1080、24 fps、432帧，含流线、局部矢量、压力、显示WSS和原始面片WSS；同时保留4K PNG/PDF。
- 原始WSS均值3.058876 Pa，P5/P50/P95=1.538445/2.822987/5.680428 Pa。
- [CPU/GPU等价核验](gpu_solver_fix/backend_equivalence.json)：独立副本修复设备解到主机视图同步后，同网格两步场及WSS相对差约1e-13量级。GPU负责CUDA稀疏Krylov，CPU8负责装配与局部LU；这是后端一致性检查，不是网格独立性验证。

## 已确认的限制

内部人工接管截面与各自端口流量存在约2%–3%的差异；独立体积散度积分与截面通量差一致。边界整体守恒不能代替局部精度证明。详见 [截面CSV](reports/internal_sections/section_flux.csv)、[独立散度CSV](reports/internal_sections/divergence_theorem.csv) 和 [对比图](reports/internal_sections/internal_flux_and_divergence.png)。本轮只有一档人脑网格，不能宣称网格无关；近壁停留和局部WSS不应视为已充分验证的物理预测。

## 用户已停止微泡批次

计划1500条、dt=0.5 ms，保留原SonoVue条件粒径分布。入口1500个样本、CUDA梯度复核，以及原/编译几何核的单条12 s轨迹和审计一致性检查已经完成。该单条为长停留时间截断，不能当作1500条批次结果。

正式批次完成条数 **0**，保留 **8** 条未完成轨迹的记录。任务及独立缓存验证均按用户要求停止：`CANCELLED_BY_USER`；自动重启关闭，主流程检查`USER_STOP.json`后退出。没有本批微泡动画，没有完整批次汇总，没有最终完整交付报告。`write_results_report.py`是待满足交付条件的脚本，不代表报告已生成。

新增`query_cache.py`是未完成验证的草稿；停止时批量进程使用的是缓存修改前已经载入的入口代码。不能把此草稿写成已通过数值验证或已用于正式结果。

停止证据：[USER_STOP.json](microbubble/data/USER_STOP.json)、[说明](microbubble/data/CANCELLED_EXECUTION_NOTES.json)、[原始目录](microbubble/)。所有中间记录原样保存，不补造完成标记。

服务器：`vast4090:/workspace/brava_flow_roi_18mlmin_20260928/`。复现说明见[REPRODUCE.md](REPRODUCE.md)；其中计算命令仅供以后明确授权恢复时参考，本次同步不得执行。GPU输入缓存因体积排除，恢复现有文件的方法见仓库同步[排除说明](../sync_metadata/current_20260929/EXCLUSIONS.md)。
