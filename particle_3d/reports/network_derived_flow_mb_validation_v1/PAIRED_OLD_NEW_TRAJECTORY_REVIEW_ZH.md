# 逐泡配对审核

同一颗泡只替换背景流；N=30，非 population statistics。所有旧结果 SHA 和完整事件元数据匹配才复用。

| ID | 直径 µm | OLD | NEW | OLD transit ms | NEW transit ms | OLD min gap nm | NEW min gap nm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 1.2854 | O2 | O3 | 31.614 | 88.268 | 2.0000 | 2.0000 |
| 5 | 0.9527 | O2 | O3 | 32.346 | 109.946 | 2.0000 | 2.0000 |
| 7 | 1.6312 | O2 | O3 | 27.880 | 82.575 | 2.0000 | 2.0000 |
| 10 | 1.3869 | O2 | O3 | 23.296 | 78.140 | 2.0000 | 2.0000 |
| 11 | 1.0794 | O2 | O3 | 25.804 | 94.794 | 2.0000 | 2.0000 |
| 12 | 1.2066 | O2 | O3 | 29.467 | 91.282 | 2.0000 | 2.0000 |
| 14 | 1.5355 | O2 | O3 | 22.652 | 84.097 | 2.0000 | 2.0000 |
| 15 | 2.4900 | STATIONARY | STATIONARY | — | — | 2.0000 | 2.0000 |
| 17 | 1.0626 | O2 | O3 | 29.223 | 89.409 | 2.0000 | 2.0000 |
| 18 | 0.8075 | O3 | O3 | 367.163 | 86.910 | 2.0000 | 2.0000 |
| 19 | 1.1544 | O2 | O3 | 20.346 | 78.793 | 2.0000 | 2.0000 |
| 22 | 1.0079 | O2 | O3 | 38.072 | 88.474 | 2.0000 | 2.0000 |
| 24 | 1.2906 | O2 | O3 | 21.956 | 83.229 | 2.0000 | 2.0000 |
| 26 | 1.8570 | O2 | O3 | 27.672 | 97.676 | 2.0000 | 2.0000 |
| 27 | 1.2300 | O2 | O3 | 30.910 | 98.389 | 2.0000 | 2.0000 |
| 28 | 1.3119 | O2 | O3 | 25.883 | 79.397 | 2.0000 | 2.0000 |
| 29 | 1.1563 | O2 | O2 | 21.482 | 25.828 | 2.0000 | 2.0000 |
| 32 | 1.1586 | O2 | O2 | 28.179 | 37.258 | 2.0000 | 2.0000 |
| 33 | 1.3339 | O2 | O3 | 21.243 | 78.045 | 2.0000 | 2.0000 |
| 35 | 1.1414 | O2 | O2 | 22.623 | 30.238 | 2.0000 | 2.0000 |
| 36 | 0.9545 | O2 | O2 | 23.618 | 30.674 | 2.0000 | 2.0000 |
| 37 | 1.1857 | O2 | O2 | 25.688 | 30.971 | 2.0000 | 2.0000 |
| 42 | 1.1482 | O2 | O3 | 28.585 | 91.007 | 2.0000 | 2.0000 |
| 44 | 0.9759 | O2 | O2 | 21.409 | 24.655 | 12.0334 | 11.6102 |
| 46 | 1.3290 | O2 | O2 | 20.216 | 26.670 | 40.5099 | 2.0000 |
| 49 | 1.1471 | O2 | O2 | 25.477 | 29.859 | 2.0000 | 2.0000 |
| 50 | 1.7532 | O2 | O3 | 24.395 | 84.033 | 2.0000 | 2.0000 |
| 53 | 1.3014 | O2 | O3 | 28.170 | 87.090 | 2.0000 | 2.0000 |
| 54 | 1.9936 | O2 | O3 | 28.239 | 97.239 | 2.0000 | 2.0000 |
| 57 | 1.5717 | O2 | O2 | 27.293 | 34.279 | 2.0000 | 2.0000 |

## 指标语义

COMPLETED 只用于原出口三角形首次中心穿越；stationary 保留为 STATIONARY，transit time 留空，物理已模拟时间单独记录。contact_events 统计已接受诊断序列中“无接触→有接触”的连续接触段起点；contact_accepted_steps 与最大同时约束数另外保存。handoff_events 只计已接受事件。solver iteration statistics 指 contact active-set 迭代，直接线性求解无迭代数时不编造该指标。planar wall correction 是已存壁面权重乘以目标切向速度与 bulk 切向速度之差的模，另存全部 unconstrained hydro−free 平移修正。所有这些指标只读取原诊断。

原始 20 列轨迹（含速度、角速度、姿态、gap、g_nf、h_lower、状态、dt）和完整 accepted/rejected 诊断保存在 `outputs/OLD`、`outputs/NEW`。所有连续 union 证书的分段范围必须无缝覆盖 [0,1]，逐个叶证书审核。

## point tracer 差异

| ID | point | MB | 出口改变（两者都完成） |
| --- | --- | --- | --- |
| 15 | O3 | STATIONARY | False |
| 26 | O2 | O3 | True |
| 27 | O2 | O3 | True |
| 46 | O1 | O2 | True |

## 图片

[01_new_flow_30_tracks PNG](figures/01_new_flow_30_tracks.png) · [PDF](figures/01_new_flow_30_tracks.pdf)

[02_paired_old_new_trajectories PNG](figures/02_paired_old_new_trajectories.png) · [PDF](figures/02_paired_old_new_trajectories.pdf)

[03_paired_outcome_transition PNG](figures/03_paired_outcome_transition.png) · [PDF](figures/03_paired_outcome_transition.pdf)

[04_transit_time_paired PNG](figures/04_transit_time_paired.png) · [PDF](figures/04_transit_time_paired.pdf)

[05_minimum_gap_and_nearwall_exposure PNG](figures/05_minimum_gap_and_nearwall_exposure.png) · [PDF](figures/05_minimum_gap_and_nearwall_exposure.pdf)

[06_point_tracer_vs_microbubble PNG](figures/06_point_tracer_vs_microbubble.png) · [PDF](figures/06_point_tracer_vs_microbubble.pdf)

[07_fluid_split_and_paired_route_changes PNG](figures/07_fluid_split_and_paired_route_changes.png) · [PDF](figures/07_fluid_split_and_paired_route_changes.pdf)

