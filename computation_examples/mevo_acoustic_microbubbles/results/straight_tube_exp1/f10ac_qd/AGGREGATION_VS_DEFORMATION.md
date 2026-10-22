# 聚集、分裂和形变分别评价

| 候选 | 频率 MHz | Δlog Rg | Δlog aspect | 新初态同号复现 | 全部保持连通 | 绝对主轴伸展超过噪声 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 2.6 | -0.065757 | 0.095171 | True | False | False |
| 1 | 2.4 | -0.001762 | 0.004319 | True | True | False |
| 2 | 2.5 | -0.076256 | 0.079426 | True | False | False |
| OD | 2.5 | -0.286524 | -0.197027 | False | True | False |


几何邻接阈值为 2.1R；硬接触另外统计。长宽比增大必须结合绝对主轴长度判读，不能自动称轴向伸长。Rigid/similarity 残差剔除整体位移、旋转和均匀缩放。局部组的固定 ID、Rg、长宽比、组间距和 2.05/2.1/2.15R 阈值诊断见 [LOCAL_COMPONENT_SHAPES.md](LOCAL_COMPONENT_SHAPES.md) 与 GROUP_DISTRIBUTION_DIAGNOSTICS.json；原 2.1R 准入规则不变。分裂后的空间重排作为 Task A 结果保留，但不能称完整材料团的形变、融合或断裂。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`


| 搜索初态轨迹 | Rg um | 长宽比 | Kabsch 旋转 deg | 剪切范数 | 分量 | 最大分量泡数 | 去相似变换残差 nm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| frequency_2600000_init_3 | 9.0029 | 2.8691 | 0.03287 | NOT_IDENTIFIABLE | 1 | 64 | 1010.033 |
| frequency_2400000_15 | 9.5978 | 2.6199 | 0.00237 | 0.000323 | 1 | 64 | 26.173 |
| bop_v2_10 | 8.9089 | 2.8242 | 0.04164 | 0.035688 | 1 | 64 | 820.284 |
| sobol_v2_20 | 7.2194 | 2.1421 | 0.63044 | NOT_IDENTIFIABLE | 1 | 64 | 1727.274 |

这些是所选轨迹的终点量；非仿射残差过大时保留原仿射有效性拒判，剪切写 NOT_IDENTIFIABLE，不填零。不能把很小旋转和剪切自动列为新可控模式。OD 跨初态验收未通过。