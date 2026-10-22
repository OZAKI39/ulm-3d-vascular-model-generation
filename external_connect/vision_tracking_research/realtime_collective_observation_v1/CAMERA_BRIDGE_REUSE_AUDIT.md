# 接收与桥接审计

读取并核对两处相机工程：`/home/lzy/projects/external_connect/CHRONOS_control` 与 `/mnt/f/CHRONOS_control`。身份见 evidence/camera_source_identity.json。检查了 chronos_common.py、配置、test_chronos.py、launcher、RTSP relay、check_connection、calibrate_scale 及历史接收/延迟记录。

LatestFrameStream 自己持有接收线程和最新帧，使用原配置，断线清除 latest，并有 reconnect/read/open timeout。其时间是 Windows 解码完成时间，非曝光。历史 1280×1024 / 约 60 Hz 是 H.264 监看预览，不是高速传感器所有帧；历史日志也记录过断流和自动恢复，不能称永久无故障。

ConfiguredWindowsBridge 继承原 WindowsLatestBridge._loop、snapshot 与 read_exact，保留 CV2FRAME、最大 header/payload、单请求在途、最新槽和 sender_age+GET_RTT。只改启动参数和回收封装。原 windows_frame_server.py 未改，仍有限 1–180 s，含紧急 watchdog。默认每 120 s 启动一个有限会话，提前结束并回收子进程再重建；轮换可能产生短暂停顿，必须失效，而不是掩盖。

本轮使用真实 Windows Python 跑合成原始 BGR 测试服务器，验证 1280×1024 Windows→WSL 管道和原 inherited loop。另有实际 OS 子进程多会话 EOF/轮换模拟。都没有导入或启动真实相机接收实例。没有修改原相机配置、RTSP 服务、IP、防火墙、录制控制或功放。

拓扑：实验室 Windows 接收 → 同机 WSL CuPy GPU（RTX4060）或 Rust 后备分析/Qt；远程服务器处理其本地录制文件、做 CuPy 研究验证。**没有部署未经验证的实验室→公网 GPU 原始图像通道。** 本机 BGR 测试结果不能推广到 WAN。当前传输无图像压缩，无有损压缩精度损失；远程 LIVE 加速仍需单独测量和接入。
