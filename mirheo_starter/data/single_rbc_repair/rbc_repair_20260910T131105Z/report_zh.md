# 单红细胞修复：实际运行与关口结果

**尚未完成同等质量比较，qualified_speedup=null。** 授权已使用；本次停止原因是物理质量关口失败，不是等待授权或显存不足。

原生库已隔离编译，未替换原安装。新增 local→halo 顺序、空 halo 早退和原生出错状态记录；另外发现并修正项目分段观测：Mirheo 每次 run() 都重新分类并清除对象力，因此改为每阶段一次连续调用、由原生插件保存中间状态。旧代码与错误的 classifier_corrections=0 原始字段保留，解释见 [run_boundary_audit.json](run_boundary_audit.json)。一个 CPU 行为回归先失败、修复后通过。

| 任务 | 原生执行状态 | 最后常规存帧步数 | 准备已计算 t*（存帧） | GPU 计费秒数 |
|---|---|---:|---:|---:|
| A0_shared_installed | STOPPED_OR_FAILED | 4100 | 4.1000 | 27.427140 |
| A1_independent_installed | STOPPED_OR_FAILED | 4800 | 4.8000 | 27.659322 |
| A2_shared_ordered_isolated | COMPLETED | 5000 | 5.0000 | 27.831135 |
| A3_ordered_instrumented | COMPLETED | 5000 | 5.0000 | 26.914817 |
| A4_continuous_native_observation | COMPLETED | 5000 | 5.0000 | 21.895217 |
| A5_continuous_preparation_30 | STOPPED_OR_FAILED | 6000 | 6.0000 | 26.722155 |
| A6_continuous_preparation_half_dt | COMPLETED | 60000 | 30.0000 | 231.720664 |

以上全是 Γ=0 的诊断；常规存帧不是精确崩溃步。A0/A1 分别以 8659/7458 个候选超过 6400 失败。A2/A3 的短完成不能证明全程修复。A4 连续 5000 步完成，三个全粒子存帧未见可靠成员不符。

A5 在零基第 6407 步、t*=6.407 的 local outer 候选生成后报错，13,675 > 6,400；常规最后存帧为 6000 步。出错前旧位置已经有 **182 对自交、WLC 最大伸长比 7.486**；当前反弹前为 **68 对自交、最大伸长比 88.425**，顶点 189 单步移动 101.623。旧/新位置与 ID 同时保存，独立 CPU 法复核，未周期折回这些失效坐标。[出错一步的原始证据与几何](A5_failure_step_geometry.json)。因此本次不是已证实可安全扩容的正常碰撞负载；最初失稳发生在更早的未保存步骤，尚不能断言是弯曲应力、热噪声或反冲中的哪一个首先触发。

A6 把 dt 减半为 0.0005，同步保持物理准备时长，完成 **60,000 步、t*=30、Γ=0**，没有原生溢出。60 个完整液体存帧（t*=0 到 29.5）中，110,592 个粒子及内侧 642 个 ID 保持不变；确认成员不符 0 个点帧，近膜不确定 17 个点帧。未保存的时间、最终 t*=30 的全液体状态不在这项全量结论内，终点另有 192 个探针。[全粒子核查](A6_full_population_membership.json)

但 A6 在 t*=28.5 的面积偏差 **2.294581% > 2%**；最后三组相邻准备窗口的形状残余率为 **0.024283、0.029368、0.028834 > 0.002**。末帧相对共享参考形状的中心去除 RMS/a 为 0.146696，这是参考形状诊断，不冒充双方新配对误差。延长至原定最大准备时间仍没有通过质量要求。[准备判据及逐窗口结果](A6_preparation_quality.json)

两边细胞尚未证明足够相同。本轮未开始新的材料拟合；原始伸长响应差 17.82%、常参考角造成的初始应力差、约 0.989 的膜/排液质量比和半步黏度差 6.47% 仍待处理。半步稳定性不能替代新的黏度和膜响应匹配。材料训练/留出标准仍见 [material_matching.json](material_matching.json)，没有改阈值，也没有声称已穷尽两组候选。

两边本轮到 Γ=4 的合格冷启动次数均为 0，端到端耗时均为空。历史 HemoCell 的 58.171937、60.782516 秒仅保留为旧配置结果，不能与本轮准备诊断相除。A6 的 231.720664 秒仅是这项零剪切诊断的完整进程费用。

本轮实际 GPU **390.170450/8500 秒**（其中失败 81.808618 秒），CPU 求解 **0/600 秒**，隔离编译 **510.347624/1800 秒**。剩余 GPU 8109.829550 秒、CPU 求解 600 秒、编译 1289.652376 秒；A/B/C 的剩余额度分别为 1209.829550/1900/5000 秒，未转移额度、未重置历史账本、未自动重试。单独 CPU 分析、图表与浏览器费用另列，编译控制器时间和墙粒子准备子项没有重复相加。

下一项最小行动是定位 WLC 越界之前的每步膜力与双侧反冲，分离初始弯曲应力、热扰动和积分响应，再决定最多两组原生材料候选。当前 A 的质量关口未通过，B 标定和 C 的 Γ=4 计时按原停止规则不启动。禁止靠容量扩展、删粒子、强制重分类或放宽阈值取得速度比。

[当前结果 JSON](comparison_results.json) · [根因状态](root_cause_evidence.json) · [候选准入](candidate_comparability.json) · [实际 GPU 账本](../../../runs/single_rbc_repair/rbc_repair_20260910T131105Z/gpu/budget_ledger.json) · [离线中文 HTML](../../../test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_review.html)

查看现有页面（不启动求解，不重新渲染）：

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.single_rbc_repair.authorized_result_report --open
```

CPU 预检命令仍是 `.venv/bin/python -B -m py_scripts.repair_single_rbc_benchmark --config py_scripts/single_rbc_benchmark_repaired.yaml --preflight-only`。原冻结 --execute 队列已有失败记录，不能作为自动重试入口。新页面渲染命令为 `.venv/bin/python -B -m py_scripts.single_rbc_repair.authorized_result_report --render`，仅读取已经保存的结果。浏览器验收绑定当前 HTML 哈希，单独见 delivery_receipt.json；人工验收 PENDING。本轮未提交或推送。
