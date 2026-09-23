![15_particle4_timestep_comparison](figures/15_particle4_timestep_comparison.png)

应该看什么：看步长减半后接触时间、间隙、末位置、速度和短轴差。

实际看到什么：三个步长均只用于人工混合场景验证，比较的是短轴而非 quaternion 分量。

有没有异常：NOT PRODUCTION TIMESTEP SELECTION；真实 RBC 通行限制仍保留。
