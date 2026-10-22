# 动作—状态数据集

32 个固定 ID，R=1.78 µm，2.5 MHz；背景流严格为零；继承平移阻力黏度 0.001 Pa·s、径向响应有效黏度 0.002 Pa·s、密度 1000 kg/m³、声速 1500 m/s。这是原模型参数，未宣称为本实验纯水校准。
 
每条 raw trajectory.h5、INPUT/RESULT/LAST_STEP/SWITCHES 均保留；dataset/ 是确定性派生，analysis/RUN_ANALYSIS.json 关联两者 SHA。全部 102 条不是 102 个独立随机试验：2 ms 标准步长具有 37 个初态—动作组，细步长、细缓存、短时前缀、延长窗口属于同组。
 
HDF5 保存实际时间、32 个固定 ID、三维位置/去质心位置、Q/Rg/主轴可辨识、两两/最近邻距离、邻接和接触、逐区间速度、径向均值/中位数、动作段 ID/幅值/相位/启停。力采样有自身时间，阻力是原阻力矩阵作用于实际步平均速度的推导量；不伪装成全时间原始采样。逐步接触反力与数值诊断可回溯原 H5 列名。
 
动作来自原 Segment/Schedule，四路顺序 Y-,Y+,X-,X+；合法源压 0..1500 Pa-equivalent。非零声压严格在 APECSS 0.1..40000 Pa 内，OFF 精确零单列处理，无静默裁剪。PAIR_OFF 在 1 ms 精确切 OFF，保持绝对时钟；零流过阻尼模型在关断后静止，不是实验中一定无残余运动。
 
不足时间窗不会填补；数据只标 SIMULATION_3D_TRUTH 与 SIMULATION_XY_PROJECTION。全部单位 SI，图表明确转换 µm/nm/ms。声学幅值不是 AFG 电压。
