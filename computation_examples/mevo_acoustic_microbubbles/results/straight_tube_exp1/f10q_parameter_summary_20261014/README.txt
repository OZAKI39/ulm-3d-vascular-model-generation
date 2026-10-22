F10Q 参数与受力整理（2026-10-14）

仅对已有轨迹做只读后处理，未启动新轨迹、未改动F10Q历史交付。
正式主表：X = NAV_X_DT1；Y = NAV_Y_DT1；两者完整2 ms，dt=1 us，声场缓存3 um。
另列OFF与旧强控制的合法0.883 ms前缀；不包含强控制拒绝提案。

坐标原点为PDMS块/通道中心；X沿32.5 mm长轴，Y沿左右方向，Z向上。
PZT发声面中心、朝向和直径来自原CONFIG_SNAPSHOT.json；四片顺序Y-、Y+、X-、X+。
source_amplitude_Pa = rho*c*v_normal，为源强等效Pa，不是电压或已标定实际PZT表面声压。
相位按Re[P exp(+i omega t)]约定，CSV同时提供有符号角度及0–360度等效角度。
频率2.5 MHz；32泡，半径1.78 um；背景流严格为零。

每个合法时间步的步初坐标均有保存。此处读取全部这些坐标，以原场缓存、原APECSS和原F2公式重新评估力，没有积分新状态。
X/Y各有2000个实际施力状态：0至1.999 ms，每步1 us；不是只统计稀疏保存的20个力采样点。
压力最小/最大值逐步核对原extra_metrics；逐泡主力、次力及总力在原保存时刻与GPU记录核对通过。

平均推力：主声场Primary Bjerknes力的周期平均向量；primary_mean_magnitude_pN为先取每步力模再对时间平均。
聊天中的单泡平均推力再对32个ID平均，不代表每个泡的力完全相同。
最大单对吸引力：原F2模型中指向另一泡的正向成对力峰值；不是硬球约束反力。
secondary_net_max_pN：来自其余31个泡的声学相互作用向量和的模的最大值，不等于单对吸引力之和或最大值。
压力均为峰值复声压幅值|P|，不是瞬时压力或RMS。primary pressure是四路PZT的叠加背景场。
SUMMARY.json另列F2成对反馈有效压力范围；该近场闭合不是完整多泡散射场。
1 pN = 1e-12 N。

pzt_controls.csv：逐片几何、启停、源幅值和相位。
per_bubble_forces.csv：每组32个固定ID的平均主推力、主力峰值、最大单对吸引力及对应伙伴/时间、邻泡声学合力峰值和声压。
*_force_evaluation_data.npz：逐时间步逐泡的复声压、主推力、邻泡合力及最大吸引对。
SUMMARY.json：均值/范围/峰值、原数据和脚本SHA、原GPU一致性核验。

X完整形态/ID数值准入通过；Y的形态/ID未收敛，因此Y吸引峰值只是该dt/cache轨迹的结果，不是已收敛或实验标定的物理常数。
EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / training_ready=false。
