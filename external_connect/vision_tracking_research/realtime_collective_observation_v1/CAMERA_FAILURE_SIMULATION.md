# 故障模拟与恢复

results/software_validation/fault_events.json 保存无新帧、短/长断流、重新连接、重复ID、PTS跳变、错尺寸、空帧、抖动/源丢帧、事件跨断流、在途停止、坏二进制的实际事件与观察结果。unit tests另含慢处理拥塞、source切换、灰度/BGR、会话重置和Qt意外关闭。

results/bridge_mock/summary.json 用真实OS管道测试继承的WindowsLatestBridge循环、多个有限子进程会话与EOF回收；results/windows_ipc_mock/summary.json 则是实际Windows Python→WSL的1280×1024合成BGR传输。后者没有相机，也不是公网带宽测量。

故障测试不表示真实Chronos发生过相同错误。重复帧/坏尺寸拒绝，长无新帧TTL失效；恢复需要合法新帧，跨session不衔接运动。stop/UI close都回收线程和子进程。真实相机端到端/重连/曝光时钟仍NOT_MEASURED。

已保留开发失败：Qt event命名冲突、录屏初版编码PTS冲突、SciPy整通道速度不足；相应修复有单独提交与测试。v4/v5任务在等待旧稳定性任务时被本任务主动停止，用最终v6快照替代，不算正式科研运行。
