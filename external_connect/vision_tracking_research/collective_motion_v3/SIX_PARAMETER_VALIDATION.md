# 六参数数学与时间定义

所有方法复用collective_motion_v1/src/affine.py。有限端点变换H的生成元为logm(H)/dt；v0在固定ROI几何中心，图像x右、y下。div=A11+A22，vort=A21−A12，s1=(A11−A22)/2，s2=(A12+A21)/2。平移px/s，其余1/s。没有用协方差变化替代速度场。

共同误差限[5,5,0.2,0.2,0.1,0.1]未变。密集流先在真实检测点取样，经前后向一致性≤1.5 px后拟合。CoTracker用官方online predictor的可见性/置信门槛获取短对应。OT采用低歧义软对应。局部流无可信仿射/空间支持不足则逐字段null。

8项本地与服务器合同测试通过：非方形ROI有限logm解析真值、流场xy顺序和越界NaN、M4 v2与v1适配、无效值必须null、断帧窗口重置、TTL、真实本机ZMQ消费者以及OT已知平移/数值失败。OpenCV多通道border NaN行为已显式修复；旧版本任一通道NaN同样被finite过滤，该修复不改变已冻结有效支持。
