# F10AA后端和历史保护

本地/服务器可信交付基准为 `6028598d0012912cc19aa766fa786c4afb0a91a5`；不是服务器旧 repo HEAD。185 条历史轨迹 SHA 与原源码已逐一核对。

直接继承 `f10aa.trajectory.Campaign.run`；`AuditEngine` 仅附加原始布尔 Gate 输出。力学采用 `WaterAccelerated`/SciPy CPU，声学保持 Warp GPU；没有重新实现粒子动力学。冻结原 src、协议与响应表/缓存 SHA。

证据：[初始审计](BASELINE_AUDIT.json)、[服务器审计](evidence/REMOTE_BASELINE_AUDIT.json)、[保护复核](evidence/FINAL_PROTECTION_AUDIT.json)。原始 INPUT 的父协议哈希属于 F10AA；F10AB 初态、计划、源身份另由阶段 START/COMPLETE、PREREGISTRATION 和 portable HDF5 属性明确记录。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
