# 有限频率与壳敏感性

采用真实四圆片尺寸/位置/朝向的均匀水 FNM 声场，Catmull–Rom 缓存的梯度取自同一个压力插值多项式。频率只在完整轨迹之间变化。PDMS 反射、弹性波、声流、自由液面和硬件传递函数均未模拟。APECSS 为原 RPAR + 可逆 Marmottant，当前半径相关壳参数、0.036 N/m 初始表面张力假设及 0.002 Pa s 径向有效黏度不作实验标定解释。

仅扰动壳弹性系数 ±20%，其余参数保持不变；三个频率均实际求解响应、固定状态受力与选择性完整动力学。三个离散频率不足以定位唯一共振峰；不宣称 2.95 MHz 是共振，也不假定过共振必然翻转整体运动。

实际数值见 ../f10ad_multifrequency/*/RESPONSE_*.json、FORCE_AUDIT.json 和 ROBUSTNESS.json。壳参数未经实验标定，任何排序变化只能作为敏感性证据。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
