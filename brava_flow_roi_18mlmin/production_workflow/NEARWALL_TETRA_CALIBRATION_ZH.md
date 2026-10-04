# 受控 near-wall tetra 标定

**STOPPED_NW1_MESH_OR_RESOLUTION_GATE**。adapter 与回归 PASS；NW1 surface/volume QA PASS、P0-like bulk PASS，但实际壁面法向高中位 **50.164 µm**，未达到20–30 µm目标，也未达到网格前预声明的≤0.5×P0检查（实测0.66345）。按停止原则，不运行CFD、不生成NW2。没有field winner；未运行oracle、轨迹、dt或BraVa/population。

## 旧场与接口

旧P0/P1原始场零CFD补算完成；所有particle bands 96/96严格包含，u/G/s/curl及绝对误差见[基线](LEGACY_P0_P1_PARTICLE_INPUT_BASELINE.json)与[短表](PARTICLE_CFD_INPUT_AUDIT_ZH.md)。P0近壁u=5.8810%–7.1799%，G=7.1008%–8.1262%；P1的u=7.7551%–8.9494%，G=7.9216%–9.0747%，均不能作为合格field。旧u、du/dn、WSS及FAIL判定未改。

（历史图件已退休；结论保留于本报告/JSON，tracked文件可从清理前算例Git提交91ebc021恢复。）
PDF（历史引用 `figures/legacy_particle_input_baseline.pdf`；结论见本报告和最终摘要，旧路径已清理）。每组96点，G使用完整Frobenius范数；连线仅引导读图，未插值或平滑原始场。

已安装TetGen接受内部点metric；固定box同点云ablation和确定性测试、旧P0表面回归、manufactured boundary QA与4项particle/population imports通过。adapter核查（历史引用 `VOLUMETRIC_METRIC_ADAPTER_REVIEW_ZH.md`；结论见本报告和最终摘要，旧路径已清理）与完整测试（历史引用 `VOLUMETRIC_METRIC_ADAPTER_TEST.json`；结论见本报告和最终摘要，旧路径已清理）。新增薄层只链接原库，不改算法或surface mesher。上一轮“Python公共接口未暴露”的结论保留在Git及策略历史字段；它不再是本轮停止原因。

（历史图件已退休；结论保留于本报告/JSON，tracked文件可从清理前算例Git提交91ebc021恢复。）
PDF（历史引用 `figures/volumetric_metric_fixture.pdf`；结论见本报告和最终摘要，旧路径已清理）。两图是同一固定中截面的实际网格线，无平滑。

## 网格前冻结的NW1设置

R=596.779229982 µm，L=7161.350759789 µm；原512周向/30轴向PLC SHA `dbd071479f3fe3679e94f00d481e5e94d82f01198d0454459fbdeaf1ff906d1b`。原23554个表面点、47104个三角形不动。

设d=R−sqrt(x²+y²)，单位m：`h(d)=30.6186e-6 + (238.711691993e-6−30.6186e-6)*clip((d−25e-6)/50e-6,0,1)`。H_shell=75 µm，无扫参。

内部点仅位于d=25、75 µm两环；每环64个等角点，半格方位偏移；z=kL/30，k=1..29。共3712点，无新增bulk支撑点，环间弦内缩上界0.6887 µm。这是离散metric支撑，不是物理interface，也不是承诺每个新节点都精确满足连续h(d)。[预声明计划](manufactured_pipe/NW1_PLAN.json)。

## 实际QA与分辨率

NW1 **30301 nodes / 121542 tetra / 121204 velocity-pressure DOF**。表面坐标逐字节保留，连接与ID一致；free/non-manifold/zero-volume/negative/all-WALL均0，单一fluid volume。实际finite-facet最大偏差0.011234 µm≤0.05；入口trace L2=0.126824%，Q误差=0.134493%，无flux normalization。[完整网格结果](manufactured_pipe/R1_NW1/mesh/NW1_RESULT.json)、[geometry QA](manufactured_pipe/R1_NW1/mesh/FIDELITY_QA.json)。

中央z/L=0.4–0.6，6144个wall facets分别对应6144个相邻tetra；法向高是相对于实际wall facet的tetra高度：

| 网格 | min / P5 / median / P95 / max (µm) |
|---|---|
| P0 | 0.2995 / 13.4885 / 75.6101 / 126.2785 / 145.9742 |
| NW1 | 14.6587 / 34.8272 / 50.1639 / 54.7719 / 103.0061 |

下表按中央tetra质心解析壁距分组。mean/median edge对组内全部六条边统计；内部normal span为顶点沿质心解析径向法线的投影范围，**不是**相对于壁面的tetra高度，也不是棱柱层厚度。

| 壁距 µm | cells | mean edge µm | median edge µm | normal span min/P5/median/P95/max µm |
|---|---:|---:|---:|---|
| 0–10 | 835 | 93.180 | 37.267 | 14.483 / 14.622 / 31.295 / 37.780 / 39.902 |
| 10–25 | 12967 | 105.009 | 127.864 | 6.453 / 28.101 / 49.860 / 52.122 / 84.501 |
| 25–50 | 2691 | 93.560 | 89.659 | 13.263 / 49.341 / 50.702 / 79.971 / 140.883 |
| 50–75 | 1177 | 100.254 | 108.666 | 33.171 / 49.905 / 50.632 / 113.213 / 140.502 |
| 75–150 | 2112 | 137.084 | 141.527 | 69.131 / 95.078 / 131.584 / 176.853 / 218.205 |

远处bulk（同中央区且质心d≥150 µm）平均边长的单元中位数：P0=183.512 µm，NW1=180.850 µm，比值0.98549，处于预声明[0.75,1.25]。因此**bulk保持；近壁高度有所改善，但未实现要求的约25 µm**。总tetra数量减少不能被解释成全局更粗或更准。

（历史图件已退休；结论保留于本报告/JSON，tracked文件可从清理前算例Git提交91ebc021恢复。）
PDF（历史引用 `figures/controlled_nearwall_nw1.pdf`；结论见本报告和最终摘要，旧路径已清理） · [图文入口](OPEN_RESULTS.html)。PNG 300 dpi；CDF、分组统计和实际切片均来自原网格。

## 为什么输入30.6 µm却没有得到25 µm法向高度

已证实的限制：冻结WALL三角形边长范围7.324–238.824 µm，长短边比中位32.61，全部含超过近壁目标的长边。锁定surface使这些边无法满足各向同性30.6 µm边长要求；输入scalar metric本身也不直接约束tetra法向高度。NW1壁面相邻tetra的对顶点来源：{'PLC': 0, 'explicit_support': 8, 'backend_Steiner': 6136}，说明内部控制点不是每个wall facet的层结构约束。

据此可解释“metric接口确实有效，仍未达到圆管normal-height目标”。**不能仅凭这一个mesh把具体贡献唯一归因于支撑点密度、质量优化或-Y限制，也不能宣称所有锁定PLC的TetGen路径都不可能。** 本轮没有扫描点云/改变优化参数来检验这些替代解释。[只读约束诊断](manufactured_pipe/R1_NW1/mesh/CONSTRAINT_DIAGNOSTIC.json)。

## 资源与停止顺序

网格+QA监督作业约6.015秒，进程树RSS采样峰值0.340 GiB，无OOM。CFD预测按网格前固定规则分别缩放已测P0/P1：831–1228秒、5.60–6.33 GiB，低于30分钟/20 GiB。这是乐观工程外推，未测NW1求解器条件数；**资源不是本次停止原因**。按顺序先命中resolution FAIL，所以新CFD次数0；不能称NW1 field FAIL，也不满足NW2“合法且达到目标后field FAIL”的前提。

## 六问结论与待验证内容

Q1：旧真实particle-input误差已完整补算，两旧场均无资格作winner。Q2：最小adapter可暴露现有volumetric metric，已实测PASS。Q3：本次固定PLC/P0 bulk保持，normal-height只降到50.16 µm，未达到约25 µm。Q4：NW1场精度未运行，不能判断u/G/s/curl是否过门槛。Q5：无winner，未执行oracle或CFD轨迹比较。Q6：未执行时间步收敛，0.5 ms既未获支持也未被否定，默认值不改。

本轮停止于已指定门槛，没有自动后续任务。若以后另行授权继续，需先解决受约束近壁mesh目标；其后仍需field gate→oracle sanity→10组CFD对照→dt收敛。当前证据不足以估计这些未授权修订的总工期。本轮没有对旧1565轨迹作成功重分类，也没有证明原速度急变/数值失败已消除。

## 清理后的保留范围

[最终摘要](FINAL_NW1_SUMMARY.json)。保留合法NW1 mesh及P0/P1 scientific baseline；volumetric metric接口已证明可用，但后续Gmsh不引用其glue、fixtures或support-cloud生成器，故专用实现退休。原生产SimVascular/TetGen不动。
