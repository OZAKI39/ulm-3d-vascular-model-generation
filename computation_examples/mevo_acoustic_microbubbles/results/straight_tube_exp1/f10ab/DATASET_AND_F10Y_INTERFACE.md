# 自包含动作—状态HDF5

自包含 `portable_dataset/*.h5` 共147条，包括失败前缀。time_s与positions_m/形状指标逐帧对应；force_time_s独立对应实际步前受力；solver_gates保留每个尝试步，gate_time_s给出其时钟。每10接受步+终点导出原始状态，没有插值；全量每步原始HDF5仍不可变保存在服务器。

服务器根：`/workspace/mevo_acoustic_microbubbles/f10ab_shape_20261018/results/straight_tube_exp1/f10ab`。全量清单与SHA：[RAW_DATA_MANIFEST.json](RAW_DATA_MANIFEST.json)。便携数据清单：[portable_dataset/INDEX.json](portable_dataset/INDEX.json)。模式：[DATASET_SCHEMA.json](DATASET_SCHEMA.json)。便携数据无外部链接，可独立读取。另有 `review_dataset/*.h5`：50µs间隔的实际状态加终点，147条全部覆盖，保留每个尝试步的完整Gate。该副本用于本地快速查阅，不能替代全量事件时间分析。10µs便携数据与每步原始数据完整保留在服务器，本地回传范围以evidence/LOCAL_RECOVERY.json为准。

3D真值和XY投影与F10Y可比较COM/Rg/Q等定义，但不是光学图像、不是视觉固定ID精度证明。所有模型/seed/调度/源码/半径/单位/失效原因保留；training_ready=false。覆盖有限、非独立样本不能当大规模训练集。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
