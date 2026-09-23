应该看什么：逐字段比较三种形状在 checkpoint 前后的几何和状态。

实际看到什么：每个粒子的全部字段差为 0，椭球原旋转矩阵、胶囊轴和保留四元数都恢复。

有没有异常：RETAINED_QUATERNION_NOT_CAPSULE_ORIENTATION；额外保存旋转矩阵是为保留原 P4 状态的最后一位。
