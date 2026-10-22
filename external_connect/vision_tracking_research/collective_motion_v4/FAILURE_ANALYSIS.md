# 失败证据

1. OpenPIV wheel未带Rust；已使用官方源码编译修复，原失败日志保留。
2. API local_median_val不支持backend参数，TorchPIV强制GUI；按真实API薄适配修复，没有伪造GPU结果。
3. CuPy s2n_ratio cache返回None；改读已经计算的原始GPU相关比。
4. 稀疏CuPy bool-mask cumsum默认int64造成kernel非法访问。修复dtype_i，旧pilot低耗时是假运行，不用于收益。
5. 初版central-ridge mask排除了人工泡中心；基于已审核标签修订，旧mask和失败数据保留。
6. 解析Jacobian试验v2/v3部分six差2.745。修复copy后仍有稀疏Huber目标近乎平坦的优化路径差异；正式对鲁棒活跃输入回退原求解器，重跑v5；失败v2/v3不进入正式图表。新的异常点回归测试覆盖该问题。
7. 低对比selector只覆盖10/120；正常人工部分ROI仍漏2/505活动中心。没有修改final配置来隐藏失败。
8. final单遍4/44错接（2对称旋转超差、2局部相反运动）。多遍仍1/44非仿射错接。
9. 定速30fps未达到完整33.3ms p95，latest-slot丢帧并重新预热，显著降低有效输出Hz。没有插值或沿用旧速度冒充新观测。
10. 开发/验证多遍各1次ALL_COARSE_VECTORS_INVALID；全部低纹理困难输入保留。
11. 第一版比较视频渲染遇到JSON null到numpy object，改为显式float NaN，未改任何数值结果；编码日志保留。

未连接Chronos、未控制PZT/AFG、未训练模型、未取得独立真实速度真值。以上均不能用软件测试通过替代。
