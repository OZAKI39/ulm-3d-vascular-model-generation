# 多速率与定速回放

| Run | Frames | p50 / p95 / max ms | Mean valid windows | Fresh group-frame Hz | Drops |
|---|---|---|---|---|---|
| original_offline_regular | 937 | 58.02 / 77.67 / 254.12 | 490.8 | offline | 0 |
| original_offline_sparse1 | 937 | 58.79 / 79.16 / 162.85 | 491.2 | offline | 0 |
| original_offline_sparse3 | 937 | 50.86 / 74.53 / 96.86 | 365.8 | offline | 0 |
| original_offline_sparse6 | 937 | 51.20 / 68.77 / 119.93 | 334.1 | offline | 0 |
| verified_offline_regular | 937 | 51.80 / 70.76 / 209.11 | 490.8 | offline | 0 |
| verified_replay_regular | 778 | 28.98 / 71.43 / 94.30 | 218.8 | 10.42 | 159 |
| verified_replay_sparse3 | 791 | 26.53 / 65.82 / 128.92 | 179.9 | 11.47 | 146 |
| verified_replay_sparse6 | 803 | 30.38 / 63.50 / 114.05 | 174.2 | 12.69 | 134 |


多速率减少的是当前精细观测，也可能降低valid空间覆盖。original四策略使用同一原拟合；verified使用已全量核对、在鲁棒活跃输入上保留原求解器的混合适配。不能把v2/v3未准入加速的快时间当正式性能。

窗口每次输出last_update_pts、age_ms、fresh、stale与reason_code；未更新窗口vx/vy均null，不拿旧场充新鲜值。全通道快速扫描始终在每个被处理的源帧执行。每窗口真实更新间隔可以从完整record序列重建，源PTS是真实时间基准。

因果回放producer依原PTS释放，队列最大1；超过负载会丢旧帧，gap后历史清空并再次预热。低处理p50可能包含无PIV预热帧，不能单独称有效观测帧率。ready_pipeline_ms、所有原始行、drops与consumer事件保留。一次性回放调度有CPU抖动，3/6策略不得仅凭单次Hz宣称普遍更优。
