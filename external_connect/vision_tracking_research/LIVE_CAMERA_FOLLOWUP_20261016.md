# Chronos 实时连接复核 — 2026-10-16

**已从 WSL 调用 `F:\CHRONOS_control` 的 Windows Python 环境，成功接收并解码真实相机视频。** 原相机工程能够工作；此前 WSL 端口超时只代表此前那一次探测，不能据此判断 Windows 工程不可用。

相机地址 `192.168.0.10`，型号读回 `CR14-1.0`，API 状态 `idle` / `videoState=live`。本次 Windows 和 WSL 均能访问 HTTP 80、RTSP 554。相机 SDP 确认为 H.264，实际解码图像为 1280 × 1024、BGR uint8。

## 实测结果

| 路径 | 测试窗口 | 实际解码 | 稳定阶段接收率 | 结果 |
|---|---:|---:|---:|---|
| Windows 原 `LatestFrameStream`，正确继承启动环境 | 30 s | 1760 帧 | 最后 5 s：60.41 fps | 1 次连接，0 次接收器失败 |
| Windows 原 `test_chronos.py --headless` 独立复测 | 10 s | 579 帧 | 最后 5 s：60.02 fps | 1 次连接，0 次接收器失败，正常退出 |
| WSL 原 HTTP JPEG 适配器 | 3 次顺序请求 | 3 张实际图像 | 请求分别约 466 / 391 / 370 ms | 成功；不能当作 60 fps 视频 |
| WSL 研究环境 OpenCV 4.11.0 | 15 s，另有 3 s 协议追踪 | 0 帧 | 不适用 | UDP 参数被忽略，实际发出 TCP SETUP |
| 独立 WSL OpenCV 4.10.0 环境 | 30 s | 0 帧 | 不适用 | UDP 协商成功，视频读取仍然超时 |

30 秒测试首次出帧约 1.384 s；原启动脚本复测首次出帧约 0.888 s。起始 H.264 解码日志存在 PPS/header 警告，随后取得持续画面，因此“接收器失败为零”不能改写为“完全没有解码警告或网络丢包”。起始探测缓冲释放也会影响全窗口接收率，所以表中使用最后 5 秒统计。

30 秒测试消费到 1760 个不同序号，1759 对相邻抽样图像像素有变化。原启动脚本复测消费 576 帧、记录跳过 2 帧；其单槽取最新帧和停止时刻会影响消费计数，不能据此计算相机或网络丢帧率。

## 两个具体问题

1. 新增探测入口初版没有复用 Windows 原启动脚本的子进程环境处理，30 秒得到零帧和 `Nonmatching transport in server reply`。现已按 `test_chronos.py` 的方式，在子进程启动前传入 `OPENCV_FFMPEG_CAPTURE_OPTIONS`。失败日志、初版脚本和输入身份仍保留在 `windows/`；修正结果另存 `windows_process_env/`，没有覆盖失败证据。
2. WSL 的 OpenCV 4.11.0 存在不同的问题：配置明明为 UDP，实际 RTSP SETUP 却为 `RTP/AVP/TCP`。已核对官方 [4.11.0 源码](https://github.com/opencv/opencv/blob/4.11.0/modules/videoio/src/cap_ffmpeg_impl.hpp#L1063)：非空选项误入默认 TCP 分支。[4.10.0](https://github.com/opencv/opencv/blob/4.10.0/modules/videoio/src/cap_ffmpeg_impl.hpp) 和 [4.12.0](https://github.com/opencv/opencv/blob/4.12.0/modules/videoio/src/cap_ffmpeg_impl.hpp) 的对应判断不同。独立安装 4.10.0 后，实际请求改为 `RTP/AVP/UDP;unicast;client_port=5000-5001`，相机成功回应，但 WSL 仍未收到可解码视频。**WSL 原生 RTSP 直连尚未解决，具体 UDP 阻断位置未确定。** 本次没有改防火墙、路由或相机固件。

可用路径是 **WSL 发起 Windows 接收进程 → 原 Windows 接收器获取相机实时画面**。这次未实现连续图像向 WSL 视觉进程转发；通过 Windows 路径取得并保存的样帧已经能从 WSL 读取。WSL HTTP 取图同样已经实测成功。

## 复用入口和环境

新增 `scripts/probe_live_camera.py` 只导入原项目接收器，支持有限时长、逐帧时序记录、样帧、源码哈希和可选 RTSP 协议追踪。没有重新实现相机采集器。

从 WSL 重复 Windows 检查，输出目录每次须不同：

```bash
cd /home/lzy/projects/external_connect/vision_tracking_research
/mnt/f/CHRONOS_control/.venv/Scripts/python.exe -B \
  "$(wslpath -w "$PWD/scripts/probe_live_camera.py")" \
  --source F:/CHRONOS_control \
  --output "$(wslpath -w "$PWD/results/live_camera_repeat_001")" \
  --seconds 30
```

接收器使用 UDP 5000/5001；探测入口会拒绝已被占用的端口。本轮检查前未发现已有接收进程占用它们，测试结束后已停止自己的接收进程。

独立 WSL 排查环境：`/home/lzy/projects/external_connect/_envs/chronos_live_opencv410`，Python 3.12、NumPy 1.26.4、OpenCV-headless 4.10.0.84，`pip check` 通过。原研究环境仍保留 OpenCV 4.11.0，过去的回放基准不重命名、不冒充当前实测。原 Windows 环境及两套相机源码/配置保持原样。安装时首次误用系统 Python 3.13 触发 NumPy 源码构建，已中止并重新使用 Python 3.12 的二进制包，安装日志保留。

## 证据与适用范围

本次完整记录：`results/live_connection_20261016_1104/`。目录标签只是标识，实际时间以各运行 `identity.json` / `diagnostic.json` 的 UTC 为准。

- `windows_process_env/summary.json`、`observed_frames.jsonl`、三张 PNG：30 秒 Windows 实测。
- `original_launcher/capture_*/summary.json`、`stream_metrics.csv`：原启动脚本独立复测。
- `http_jpeg/`：三张真实 HTTP 图像和请求耗时。
- `wsl_trace/rtsp_control.jsonl`、`wsl_opencv410/rtsp_control.jsonl`：修正前后实际协商；失败结果完整保留。
- `camera_sdp.txt`、`opencv_source_evidence.json` 和对应官方源码：编码与版本诊断依据。
- `camera_state_comparison.json`：型号、序列号、state、videoState、分辨率、framePeriod、exposurePeriod 前后相同。
- `protected_sources_check.json`：Windows 和 WSL 各 5 个核心源码/配置哈希均未变化。
- `visual_inspection.json`：已查看实际样帧，显微场景可见，未进行生成或增强。

本次仅验证相机接流。约 60 fps 是预览视频接收率，不能当作高速传感器采集率。接收后消费帧龄 p95 约 2.29 ms 不包含曝光、相机编码和网络缓冲，**不是端到端延迟**。没有验证真实微泡身份、团形态或物理闭环，没有操作 Tektronix/PZT。测得 live 接流可用不会改变先前关于检测误差、ID 不稳定和物理控制未准入的结论。
