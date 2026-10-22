# 复现本轮诊断

主报告：[CAUSE_REPORT_ZH.md](CAUSE_REPORT_ZH.md)。机器入口：[CAUSE_RESULTS.json](CAUSE_RESULTS.json)。

此目录中的脚本只写此目录。原正式运行、CFD、壁面、生产模块保持只读。不调用 `population_v1.run`，不重跑完整群体。`FROZEN_PLAN.json` 固定四个事件、时间窗及60分钟累计／单窗口10分钟上限。`supervise.py` 复用项目的进程监督器，计入失败尝试，拒绝超出累计预算的新任务；不会自动提交下一项。

服务器上使用已有 `population_v1/env/bin/python`，生产源码来自同一算例的 `population_v1/source/particle_3d/src`。服务器路径以运行 `IDENTITY.json` 为准。这里不保存服务器连接配置。源码SHA对应历史提交 `98d705f13babdfa0cba6c12331623292a588cb77`；195项输入的路径以服务器算例根为相对根，GitHub共享源码对应 `ulm_particle_formal_p9a5/particle_3d/src/particle_3d/`。

首次在该诊断目录执行时，依次运行：

```bash
python supervise.py prepare
python supervise.py diagnose
python supervise.py replay
python supervise.py increments
```

`prepare.py` 校验身份，从原始轨迹取四窗，按原MANIFEST仅读重叠状态块；逐段查实际邻居并保留原始出生状态。`diagnose.py` 复用GPU场、有限三角形壁面、阻力和原CPU回退，重算接受位置RHS，保存矩阵及界面探针。`increments.py` 用原求解器拆分每个相邻状态的平动、curl、剪切、阻力、运动学约束与残差项。

`replay.py` 仅从822的实际出生状态开始闭合重积分，保留其它历史粒子的事件时钟，并连续检查与这些粒子的原路径无近场相遇。原步长已保存时读取该结果，不重复积分。其它三个事件缺少中途安全监视器checkpoint，只做记录位置的RHS重算，不能将其命名为闭合回放。

本次原/半步的差异触发了预留的单次四分之一步控制；实际决策保存在 `QUARTER_DECISION.json`。该控制不是自动队列：

```bash
python supervise.py quarter
python supervise.py finalize_evidence
```

`finalize_evidence.py` 补齐最后接受线段的穿越记录，核对显示均值/掩码并重新计算输入SHA。它跟在新生成的 `result_*.json` 之后执行；不要单独对已补齐的结果重复追加。不同阶段的 `monitoring/*/EXECUTION.json`、日志和进程预算均保留在实际工作目录；提交的 `EXECUTION_RECEIPT.json` 汇总真实执行记录及工程错误，原始监控文件不作为正式科学输出。

在本地用现有科学Python运行新增测试（将 `PYTHONPATH` 设置为现有共享粒子 `src` 目录），再生成交付：

```bash
python -m pytest -q -p no:cacheprovider test_diagnostics.py
python build_delivery.py
```

本次实际测试结果保存在 `TEST_RESULTS.json`。脚本生成的两张PNG是独立分析图，不改旧动画；速度来自真实接受时刻，未滤波。原始短窗是 `original_*.npz`，逐状态记录是 `measured_*.npz`，方程分量是 `increments_*.npz`。回放只保存822的三个短窗，不复制CFD或完整历史轨迹。

`CAUSE_RESULTS.json` 保存完整事件记录及限制；`result_*.json` 是便于逐阶段复核的事件记录；`ledger_*.json` 区分旧日志和本轮邻居测量。所有JSON禁止非有限数；未知项以null和原因保存。过去的失败尝试不作为科学通过，最终有效结果与其执行顺序见receipt。

本轮的针对性测试检验诊断代码及保存结果，不改变原场FAIL、未完成的生产资格验证或原1565条轨迹的终态分类。三个近壁窗口的时间准确性仍未评价。
