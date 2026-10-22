# 现有诊断资料

这里保存以前的研究证据；本次没有恢复研究队列，也没有把失败改成成功。可执行的当前诊断数值函数统一在 `src/human_vessel_flow/`，不从本目录动态加载脚本。

2026-10-23 补齐了原始 Gmsh V2 生成器 `reference_pipe/gmsh_reference/reference_v2.py` 及其少量冻结前置输入，并保留早期同步工作区中缺失的 P1 检查点、P2 边界输入和执行记录。该生成器保持原 SHA，属于来源证据，不作为当前运行入口；本轮没有安装新 Gmsh、重新生成网格或恢复试算。新增保留位置和 SHA 见 `docs/cleanup_20261023/preservation_receipt.json`。

| 需要查看的问题 | 原报告 | 应保留的含义 |
|---|---|---|
| 旧失败和异常速度 | [原因分析](field_events/failure_speed_review_20260929/CAUSE_ANALYSIS_ZH.md) | 失败与观察截止单独记录 |
| 原速度急变事件 | [事件原因定位](field_events/speed_event_cause_review/CAUSE_REPORT_ZH.md) | 原事件、账本和模型作用分开理解 |
| 事件区背景流场 | [流场可靠性报告](field_events/field_reliability_event_review/FIELD_RELIABILITY_REPORT_ZH.md) | 全局守恒不能证明局部场准确 |
| 稳定化与局部守恒 | [审计报告](field_events/stabilization_local_conservation_review_20261006/STABILIZATION_AND_LOCAL_CONSERVATION_REVIEW_ZH.md) | 原生实现与误差证据 |
| V2 参考管场误差 | [误差分解](reference_pipe/V2_FIELD_ERROR_DECOMPOSITION.json) | MIXED_FIELD_ERROR；已从保存数组只读复算 |
| P1 对轨迹的影响 | [敏感性记录](reference_pipe/P1_FIELD_PARTICLE_SENSITIVITY.json) | 诊断性结果，不能称正式收敛 |
| P2 曲面压力问题 | [根因报告](p2p1/common_curved_geometry_review_20261006/MIXED_CURVED_GEOMETRY_ROOT_CAUSE_ZH.md) | 保留未通过检查 |
| P2 资源停止 | [最新重执行报告](p2p1/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md) | STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT |
| V3 网格 | [提案文件](reference_pipe/PROPOSED_V3_MESH_POLICY.json) | 仅提案，没有新场 |

原报告中路径、命令、源码身份和当时的 HTML 导航保留原文，可能指向未迁入的旧中间产物。它们是历史证据，不是当前执行说明。当前入口、环境和测试命令请使用根目录 README；新的图表入口为 [迁移核查页](../results/reviews/migration_final/OPEN_RESULTS.html)。
