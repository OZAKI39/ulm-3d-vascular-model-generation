# 下一轮建议

在当前 2.5 MHz 均匀水代理及已测试初态/时间窗内，存在经过数值验证、动作可区分的有限形状响应自由度；局部两参数分析仅支持一个可辨识形状方向，尚未证明多个独立形变方向、跨初态稳健的完整连通团目标变换或可逆控制。

优先研究保持邻近网络的低幅/分段动作和动态连通约束，将“固定集合改变分布”与“完整团变形”分别设目标；保留真实非黏附物理，不靠人工连杆锁团。扩展相位/幅值动作前，先寻找超过噪声的第二形状方向，继续双幅中心差分和新seed验证。再以更长但有限窗口研究恢复/不可逆重排。

在真实实验前校准声压空间梯度与电声映射，尤其Y近场节点；当前数据只适合接口/小规模探索，不启动大型世界模型、RL或MPPI。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
