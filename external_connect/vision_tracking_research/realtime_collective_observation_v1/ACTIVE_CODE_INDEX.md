# 活跃代码索引

| 文件 | 职责 |
|---|---|
| src/protocol.py | FramePacket、尺寸/时钟检查、连续性、单槽、CV2FRAME 错包测试 |
| src/sources.py | PyAV 源 PTS 调度、0.5/1/2×、暂停恢复、循环新会话 |
| src/live_adapter.py | 继承 WindowsLatestBridge._loop、可配置启动、有限轮换、进程回收 |
| src/observer.py | 原 F10X Channel/Backend、因果背景、有效窗口、同一运动核心 |
| src/aggregation.py | cKDTree、DBSCAN、固定区域时间窗、收缩/汇聚与拒判 |
| src/pipeline.py | 独立分析、TTL 发布、localhost ZMQ、资源采样、异步记录 |
| src/ui.py | 一个 PySide6 界面、图层、曲线、事件、几何审核、实际 Qt 抓屏 |
| scripts/full_study.py | 937 帧逐帧完整分析、1×/2×/0.5×实际节奏 |
| scripts/synthetic_validation.py | 固定种子的数学与真实图块图像两级测试 |
| scripts/software_validation.py | 真实 PTS 核查、故障事件、节奏和会话测试 |
| scripts/windows_ipc_validation.py | 真实 Windows 进程—WSL 原始 BGR 管道；不接触相机 |
| scripts/render_deliverables.py | 六组图、完整真实录像与合成验证视频、解码验收 |
| scripts/remote_job.py | Git 不可变归档、独立 Supervisor 标签、总时限 |

历史只读复用：collective_motion_v4/src/channel.py、backends.py、xcommon.py；collective_motion_v1/src/affine.py；collective_motion_v2 的 splat 与 Writer；smart_online_tracking_v2/live_bridge.py 和 windows_frame_server.py。未重写 FFT、PIV、通用聚类或 RTSP 解码器。
