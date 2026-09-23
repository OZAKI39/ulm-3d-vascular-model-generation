### 图 09：09_real_fem_near_wall_mb_lubrication

![图 09](figures/09_real_fem_near_wall_mb_lubrication.png)

应该看什么：分别看原最小间隙样本和原轨迹中的靠墙样本，注意法向速度正负。

实际看到什么：原最小间隙样本正在离墙；另一个原始靠墙样本的法向速度显著减小，切向变化仅为浮点误差。

有没有异常：没有移动中心或缩放球；这是局部法向修正，尚未包含切向、转动耦合、曲率和多墙效应。

数据：[09_real_near_wall.json](data/09_real_near_wall.json)、[09_real_near_wall.csv](data/09_real_near_wall.csv)。
