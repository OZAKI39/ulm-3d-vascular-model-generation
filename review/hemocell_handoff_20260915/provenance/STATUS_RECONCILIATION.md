# 历史报告与当前交接状态

来源顺序：本次用户明确给出的当前科研状态 + 已密封的最新实际报告；历史报告保持原字节。

| 项目 | 历史原件 | 当前解释 |
|---|---|---|
| Human ParaView | Step3 / formal Step3C 报告为生成时的 PENDING；Step2 原始记录已 PASS | 本次用户要求的当前 HUMAN_PARAVIEW_REVIEW=PASS。该声明未附新的独立人工审查文件；明确归因于用户，不伪造审阅时间或截图。RBC review 仍 PENDING |
| PBS/BSA numerics contract | `STATIC_CONVERSION_PASS_RUNTIME_PENDING` / runtime PENDING | 合同在运行前冻结，不能回写；500/5000 smoke 与 LOCAL final audit 提供运行 PASS |
| Restore summary | PBS_BSA_RUNTIME_SMOKE=NOT_RUN，RBC_STAGE1=NOT_RUN | 表示恢复阶段结束时的状态；后续两个任务的报告是最新证据 |
| RBC Stage1 | FAIL，NEXT_STEP=USER_REVIEW_BEFORE_RBC_STAGE2 | FAIL_GEOMETRY_GATE 是失败类型的准确归纳；用户本次指定下一阶段为 RBC_GEOMETRY_AND_WALL_COMPATIBILITY_STAGE，不授权跳过 Stage1 |
| RBC target volume | nominal 50 µm³ | actual 45.046330078320075 µm³；两者均保留，不声称体积目标已复现 |
| Build PASS | HemoCell CPU library + native mesh probe | 没有耦合 timestep production binary，也没有 runtime monitors 的实现证据 |
| Stage4 vs smoke speed | 109.7730426 vs 163.8231326 steps/s | 介质和监测 cadence 不同，不能作为同口径性能比较 |

本次所复制原件、报告与合同没有改状态或数字。当前汇总是新的文档层。源码中的变量 `token` 是 JSON key 的解析字符串，不是认证 token；它在 driver、生成器和对应 diff 中出现，已按上下文人工检查。
