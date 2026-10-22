# 实时流水线

1. Frame acquisition：独立源线程，按 PTS 释放或持续消费原桥接 latest。
2. LatestSlot：一个未处理帧；覆盖计数。不追赶历史积压。
3. Observation worker：全通道扫描、原 PIV、连续趋势，GPU 上下文和固定几何复用。
4. State publisher：独立 20 Hz TTL 检查，localhost ZMQ PUB，HWM=1，非阻塞；端口保存在 session/publisher.json。
5. Qt UI：默认 20 Hz 请求快照；默认最新原图持续刷新；较早但未过期的分析图层淡化并标记来源帧与年龄。勾选 Measured frame 可查看逐像素对应的已分析原帧；过期时只保留原图诊断与醒目标志。
6. Recorder：有界 32 项异步队列；默认保存紧凑状态与区域指标，避免逐帧重复固定网格、掩码、多千数组。完整科研分析另存全部候选与 PIV。
7. Resource monitor：独立采样 CPU/RSS/GPU，nvidia-smi 不阻塞 TTL publisher。

检测有效性与运动有效性分开：新鲜检测仍可显示，PIV 不可信时不显示有效速度。TTL 默认 250 ms，是展示上限；报告同时给出年龄分布，不将它解释为曝光到显示延迟。历史 F10X 100 ms gate 对照须使用相同门槛重新统计。

运行内存中的历史曲线最多 600 点，单元时间窗最多 90 点，profile 每阶段最多 6000 值，资源最多 3600 点，源和 recorder 队列有界。完整长期记录在磁盘。长稳态 summary 的分位数是最近 6000 样本，完整统计可从紧凑日志重算。

相机错误无录像 fallback；停止先回收源，再等待在途分析完成、关闭发布者、flush recorder。只有本任务独立 Supervisor 标签被操作，其他服务保持原状。
