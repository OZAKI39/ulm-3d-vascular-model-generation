### 图 11：11_particle5_timestep_comparison

![图 11](figures/11_particle5_timestep_comparison.png)

应该看什么：看球墙和双球在步长减半后的间隙、速度和耗散差异。

实际看到什么：两类案例完整覆盖物理时间，保持正间隙；末间隙差与独立隐式解析关系的误差随细化减小。

有没有异常：图中所有步长仅用于验证，不能据此选定生产步长。

数据：[11_timestep.json](data/11_timestep.json)、[11_timestep.csv](data/11_timestep.csv)、[11_wall_dt1_states.json](data/11_wall_dt1_states.json)、[11_wall_dt1_states.csv](data/11_wall_dt1_states.csv)、[11_wall_dt2_states.json](data/11_wall_dt2_states.json)、[11_wall_dt2_states.csv](data/11_wall_dt2_states.csv)、[11_wall_dt4_states.json](data/11_wall_dt4_states.json)、[11_wall_dt4_states.csv](data/11_wall_dt4_states.csv)、[11_pair_dt1_states.json](data/11_pair_dt1_states.json)、[11_pair_dt1_states.csv](data/11_pair_dt1_states.csv)、[11_pair_dt2_states.json](data/11_pair_dt2_states.json)、[11_pair_dt2_states.csv](data/11_pair_dt2_states.csv)、[11_pair_dt4_states.json](data/11_pair_dt4_states.json)、[11_pair_dt4_states.csv](data/11_pair_dt4_states.csv)。
