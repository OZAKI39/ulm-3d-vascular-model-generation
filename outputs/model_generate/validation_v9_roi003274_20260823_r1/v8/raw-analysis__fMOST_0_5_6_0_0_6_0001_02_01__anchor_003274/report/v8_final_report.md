# v8 local smooth-union 最终诊断

最终结论：**KEEP_V7_HARD_MIN**

1. 锯齿最早出现于 S0_raw_flying_edges。
2. 是；缺陷与 winner-branch switching boundary 的近邻覆盖率为 defect→switch 26.808%，switch→defect 93.939%。
3. hard-min gradient jump P99 = 92.9025°。
4. 否。
5. 否；三档参数均未通过 silhouette 消除条件。
6. 无可采纳参数。
7. 所选结果最大 junction local volume 增量为 0.000%。
8. radius P95 = 0.679%，可接受。
9. 最大 branch-local hydraulic resistance 相对变化为 0.000%，未显著改变。
10. KEEP_V7_HARD_MIN

所有 S0–S3 与 k/r 对照图均使用同一相机、flat shading、wireframe；
smooth candidate 的 FlyingEdges、第一次 Newton 和第二次 Newton 始终使用同一 Phi_v8。
端口中心线延伸、端口平面和端口法向未修改。
