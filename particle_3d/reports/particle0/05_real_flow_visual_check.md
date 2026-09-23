# P0.5 检查说明

固定 seed=20260925，按四面体体积加权选取单元，再在单元内均匀取点，共 400 个位置。速度范围 1.731e-07–1.063e-03 m/s，有效样本 NaN=0、Inf=0。

vorticity 是速度场在局部的旋转趋势；strain rate 是流动在局部拉伸和剪切的快慢。图中应检查箭头与血管方向是否协调、颜色是否出现孤立异常；箭头长度只作显示缩放，原值在 CSV。

这些图仍等待用户人工审核，自动检查没有代替视觉判断。

![检查图](figures/05_real_flow_velocity_vectors.png)
![检查图](figures/06_real_flow_scalar_diagnostics.png)
