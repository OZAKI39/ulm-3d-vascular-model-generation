![01_pair_gap_geometry_matrix](figures/01_pair_gap_geometry_matrix.png)

应该看什么：看六类配对在分离、接触、穿透时是否返回正确符号。

实际看到什么：六类均有固定方向和旋转几何检查，输入交换后间隙不变、法向反向。

有没有异常：负间隙只用于几何诊断，不作为轨迹初始状态。
