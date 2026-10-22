# 适用范围与未完成事项

采用真实四圆片尺寸/位置/朝向的均匀水 FNM 声场，Catmull–Rom 缓存的梯度取自同一个压力插值多项式。频率只在完整轨迹之间变化。PDMS 反射、弹性波、声流、自由液面和硬件传递函数均未模拟。APECSS 为原 RPAR + 可逆 Marmottant，当前半径相关壳参数、0.036 N/m 初始表面张力假设及 0.002 Pa s 径向有效黏度不作实验标定解释。

径向振动来自离线稳态孤立泡响应；F2 是原有成对单极子反馈假设，不是完整多泡径向耦合或近场多极理论。

Four-source G3、分段 G4、频率连续优化、变频瞬态、长期可逆导航、千泡生产、硬件控制、AI 训练均 NOT_RUN。群体有限响应张成空间不能当作局部可控性秩。局部 Y 弱响应与频率变化强烈受均匀水近场结构影响。任何场级/数值 PASS 都不是 PDMS 实验标定。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
