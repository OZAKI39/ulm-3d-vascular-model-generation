# Particle-1 小步骤 3 检查说明

vorticity 是流场的局部旋转趋势，在这个刚体旋转人工场里等于整体旋转速度的两倍。因此球的旋转速度必须乘 1/2；49 个位置的最大误差为 0.000e+00 1/s。专门的永久测试用三个非零分量检查，误写成 Omega=vorticity 会失败。

![检查图](figures/03_pure_rotation_validation.png)
