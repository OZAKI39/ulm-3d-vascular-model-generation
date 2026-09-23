# Stage 1.8 cleanup notice

用户在 Stage 3 请求中明确要求删除 MMG/TetGen 活动实现、专用测试、独立安装环境及原始实验结果。
本阶段的最终科学结论仍为 **KEEP_STAGE017**；用户已确认 Stage 1.8 人工审核 PASS。
原 REPORT.md、assessment_summary.json、final_status.json 和最终 PNG 保持原字节，报告中的历史状态未追写。

原实验代码和独立 mesh-tools 环境已从活动工作区移除。后续 production pipeline 不依赖 MMG/TetGen。
此前归档中的实验代码和原始结果也已删除，仅保留已有最终报告、结论 JSON 和图像，没有建立新的代码备份。
历史源码可由已有 Git 提交 `96db3f4` 和 `28c9413` 追溯。

根据用户要求，旧报告中指向已删除策略、日志、工具环境和原始结果的链接不再有效；最终结论和图像仍可阅读。
通用边界身份检查、minSICN 评估、成本计数和邻域几何模块保留，供 Stage 3 复用。
具体删除范围、依赖审计与执行结果见 [Stage 3 cleanup plan](../stage03/stage018_cleanup_plan.json)
及 [Stage 3 cleanup audit](../stage03/cleanup_audit.json)。
