![04_hard_contact_sliding_validation.png](figures/04_hard_contact_sliding_validation.png)

应该看什么：入墙、平行、离墙、无接触和旋转椭球的两种速度箭头。

实际看到什么：只在接触点向墙运动时补偿法向平移，其他情形不改速度。

有没有异常：法向残差和切向变化均在浮点预算内。
