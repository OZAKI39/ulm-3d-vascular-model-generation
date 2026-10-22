# 中心差分、SVD与误差门槛

| 初态 | 窗口/ms | 差分检验 | 形状奇异值 | 保守阈值 | 有效形状秩 |
| --- | --- | --- | --- | --- | --- |
| S2_r0 | 0.5 | PASS | 3.756, 0.1437 | 0.8267 | 1 |
| S2_r0 | 2.0 | PASS | 10.94, 0.1795 | 1.437 | 1 |
| S4_r0 | 0.5 | PASS | 3.51, 0.04178 | 0.7147 | 1 |
| S4_r0 | 2.0 | PASS | 10.67, 1.36 | 2.132 | 1 |

输入归一化：1单位=公共X双源幅值5%或X+相位6°；采用±全幅和±半幅中心差分，全部半幅端点另做dt/2与细缓存。输出COM/20nm、log(λ1/λ3)/0.01、log(λ2/λ3)/0.01、Rg变化比例/0.01、similarity RMS/Rg0/0.01。形状子矩阵仅取两个谱比例和去统一尺度后的残差，不包含COM和纯Rg尺度。

阈值=max(两差分模板的谱范数差、5×数值噪声范数×2、1e-6)。四种初态/窗口组合的第二奇异值都未超过保守门槛，因此仅支持1个可辨识局部形状方向。不是全局可控性证明，也不证明其余所有合法动作都无效。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
