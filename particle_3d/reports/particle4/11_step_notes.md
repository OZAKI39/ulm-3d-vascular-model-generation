![11_wall_plus_particle_multicontact](figures/11_wall_plus_particle_multicontact.png)

应该看什么：看粒子被 WALL 与另一粒子夹住时是否两边都满足约束。

实际看到什么：WALL 和粒子接触进入同一个求解，球和 RBC 混合案例均通过。

有没有异常：固定 WALL 的作用不要求所有粒子的平移修正总和为零。
