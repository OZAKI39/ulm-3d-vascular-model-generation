# 按频率重新计算声场

采用真实四圆片尺寸/位置/朝向的均匀水 FNM 声场，Catmull–Rom 缓存的梯度取自同一个压力插值多项式。频率只在完整轨迹之间变化。PDMS 反射、弹性波、声流、自由液面和硬件传递函数均未模拟。APECSS 为原 RPAR + 可逆 Marmottant，当前半径相关壳参数、0.036 N/m 初始表面张力假设及 0.002 Pa s 径向有效黏度不作实验标定解释。

2.4/2.5/2.6 MHz 均重新建立 3 µm 与 1.5 µm 缓存；48/96 阶积分、独立 Rayleigh 面积分、解析导数与有限差分一致性均有记录。2.5 MHz 新缓存与历史节点并非逐位相同；新缓存独立积分及插值误差通过门槛。正式历史轨迹继续使用原缓存，轨迹回归逐位一致。Y 近场零点及 ±5 µm 初态扰动单独保存。

| 频率 MHz | 声场/响应 | 完整动力学准入动作 | 所有测试动作通过 |
| --- | --- | --- | --- |
| 2.4 | True | X180, Y90, FOUR | True |
| 2.5 | True | X180, Y90, FOUR | True |
| 2.6 | True | X180, Y90, FOUR | True |

证据在 ../f10ad_multifrequency/*/FIELD_GATE.json；48 个独立空间点的检验是有限采样数值证据，不是全域误差定理。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
