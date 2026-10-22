# Production workflow 当前结论

2026-10-04。当前任务完成并停止：**MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。原正式field状态继续是 **REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC**。

| 本轮问题 | 结论 |
|---|---|
| 0.377 µm 的8.52%误差来自哪里？ | P1 space=9.759%、CFD−P1=4.688%，cos=−0.488，存在抵消；P1空间项中切向几何项比径向项大。 |
| 29.8–149.2 µm bulk是否同样受P1限制？ | 29.8/59.7 µm主要是离散解差异；149.2 µm支持P1表示限制，不能统一归因。 |
| 为什么G/s/curl通过而部分u失败？ | 近壁u_exact很小且受切向插值亏损影响，导数有有限尺度、采用8%门槛；原粒子带导数实测约4.4–5.0%。 |
| P1表示是否影响生产粒子轨迹？ | A局部transit差5.149%触发LARGE，B/C为MODERATE，D为SMALL；无径向分离或contact类别变化。仅0.2 s单方位诊断。 |
| 若以后只做一次V3，应改多大？ | 条件估计：wall altitude≤16 µm，周向弦长≤7.075 µm（532段时7.048 µm），相关中央core名义bulk target=120–130 µm；全部PROPOSAL_ONLY。 |

[完整诊断报告](V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md) · [机器分解](V2_FIELD_ERROR_DECOMPOSITION.json) · [粒子敏感性](P1_FIELD_PARTICLE_SENSITIVITY.json) · [V3提案](PROPOSED_V3_MESH_POLICY.json) · [两张新图与HTML入口](OPEN_RESULTS.html)。

已完成：冻结输入SHA审计、960点三项向量分解/交叉项/绝对量/分量分析、G/s/curl简表、3840点radial profile、两个discrete-limited bulk带的零新CFD定位、4×2条production动力学诊断、5项针对性测试。粒子服务器计算约12.03 s，4CPU进程、峰值RSS1.269 GiB；无新增CFD或meshing。

原GPU CFD已在88.91 min内达到连续五区间稳态并保存第71步，CUDA/RA身份保持。原field的失败没有被覆写；（历史图件已退休；结论保留于本报告/JSON，tracked文件可从清理前算例Git提交91ebc021恢复。）继续有效。正式CFD-vs-oracle轨迹和三个dt的canonical convergence仍未运行，不能写READY，也不能宣称dt=0.5 ms已获支持。

服务器job已退出、TASK已移除，无自动下一任务。V3建议只是本地设计模型，不是已验证h²收敛率或PASS保证。等待用户新任务指令，不自动改gate、不生成V3、不运行CFD、不进入BraVa。

正式BraVa `cases/balanced_pressure_final` 与旧1565条population保持原身份；141观察截止、2数值失败不改称成功，新固定灌注协议未运行新正式population。未来新BraVa边界设计仍按“中心线阻力估算→出口压力”。
