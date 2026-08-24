# v7 Unified PolyBall 最终核验

最终决策：**KEEP_V6_HYBRID**

- 全 CFD_DOMAIN centerline/radius 由一个连续 PolyBallLine 场重建。
- 显式管面、局部隐式面拼接、surface Boolean、collar 和 hybrid interface 均为 0。
- 完整 wall 在一个等值面提取后做 Newton 投影、各向同性重网格及再次投影。
- 最终 CFD 端口先由 tail 穿过平面，再精确裁剪并生成平盖。

- radius P95: v6=0.003220318, v7=0.0067853064
- former merge-ring worst P99: v6=22.876316°, v7=13.882031°
- v7 hybrid_interface_edge_count=0
- topology: boundary=0, nonmanifold=0, self-intersection=0, internal faces=0, internal caps=0, components=1

决策同时依据 real-ROI 定量比较与 figures/ 中相同相机 flat/smooth/wireframe/silhouette 图；运行时间不单独否决 v7。
