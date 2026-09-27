# Partial FGW 语义迁移退休记录

曾评估 Partial FGW，以允许 TopBrain 与 BraVa 之间不完整的几何对应，并通过多个 donor 提供血管段语义。

真实 TopBrain Pilot 包含 MRA001–005 双侧。原三分类 NN macro F1 为 0.774811，最佳 FGW（α=0.7，运输质量=0.90）为 0.771080，因此三分类准入状态为 `TOPBRAIN_LABEL_TRANSFER_NOT_VALIDATED`。

研究目标是联合 M2+M3 的 MeVO ROI。将 M2/M3 合并且保留原 UNKNOWN 后，NN MeVO F1 为 **0.977445**，FGW 为 **0.945324**；NN 的 M1 错纳率也更低（0.088139 对 0.095358）。十侧中 NN 胜九侧，FGW 胜一侧。FGW 边界误差稍低，但未形成预定主要终点上的生产优势。最终状态为 `MEVO_BINARY_NN_PREFERRED`。

因此正式退休 Partial FGW 实现、运输集成、参数搜索、矩阵调试、专用测试及 POT 依赖声明。正式路线保留原 Open3D similarity registration 与 multi-donor NN，并增加简单的支持质量拒绝。

原始 comparison JSON/CSV/Markdown、Pilot summary、TopBrain semantic cache、配准矩阵和 donor 排名保留在 `outputs/topbrain_brava_transfer/`。NN 补存预测保持原样，新生产冻结文件另存。历史环境中已安装的 POT 不执行卸载；生产代码和环境检查均不再依赖它。
