# 实验室接入：从回放切换到真实 Chronos

本轮 `LIVE_CAMERA_PASS=false`。下面是已实现接口的操作步骤，实机质量必须到实验室验收。

1. **分工**：Windows 保持已验证的 F:\CHRONOS_control 接收环境；同机 WSL 运行本程序、CuPy 本机 GPU 分析/Qt（Rust 后备）。服务器用于 CuPy 录像研究，不要求它访问 F 盘。
2. **环境**：本机 `.venv` 已准备好。使用 start 脚本，勿使用系统 Python。本机CuPy 13.6与经过F10X修正的原开源PIV源码已放在独立环境；无需修改原 Windows Python 或重装相机服务。新机器的最小依赖和 Rust 来源见 requirements_runtime.txt 与 environment locks。
3. **先验证回放**：运行 `bash /home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1/start_replay.sh`。确认原图、计数、状态和图层正常；到 EOF 后点 Reconnect 重播。
4. **相机准备**：在实验室按原 CHRONOS_control 使用说明连接和开机，检查已有配置。只用原来可用的 IP/端口/预览设置，不更改 RTSP 服务或防火墙。
5. **避免双重接收**：可先用原预览确认画面，随后正常关闭原预览窗口，再启动本系统；两者不应同时占用 UDP5000/5001。无需另开一个永久旧接收器，本程序会自动调用它。
6. **启动 LIVE**：运行 `bash /home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1/start_live.sh`。确认 Source=LIVE_CAMERA，再点 Start。连接失败应出现错误或 DISCONNECTED，绝不会自动显示录像为 LIVE。
7. **确认真直播**：在光路里做一个安全、可识别的现场变化，核对画面随之变化，并检查原相机状态。禁止仅凭文件名或帧号认定来自实时场景。没有新鲜场景证据前不开始科研记录。
8. **按当前画面绘制区域**：系统会检查区域与实际图像分辨率是否一致。历史预览可能是 1280×1024，不能把 1476×1120 录像 mask 缩放后当作验证通过；不确定时在当前画面重新绘制。
9. **画完立即生效**：点击 **Draw LIVE mask**，沿通道内部边缘按顺序单击，至少三个点，避开管壁和背景。用 **双击最后一点、按 Enter 或点击起点** 完成闭合，程序立即切换分析区域、自动保存，不必再点击保存或修改 JSON。`Backspace` 撤回最后一点；`Esc` 取消本次绘制、保留原有区域。画错可重新 Draw 替换。交叉、退化、越界的区域会显示原因，不会覆盖原区域。**Clear current mask** 清除当前 mask、候选和运动结果；不会删除旧记录。更新区域不会断开或重连相机。
10. **直接观察，不必设置参考**：新画的区域默认使用 `IMAGE_SPACE_UNREGISTERED`，即图像坐标下的运动；画面提示 `Image-space motion | camera drift not corrected`。相机晃动可能混入运动，不能将其解释为已校正的实际流速。程序仍检查帧时效、PIV 有效窗口和候选支持，不会因“画完区域”而制造有效速度。原有配置若有固定壁面 `registration_roi`，仍执行参考检查，不通过时不显示有效速度。**Load saved mask (optional)** 可以直接加载旧 JSON；缺少参考的旧文件会明确转换成图像坐标观测，并另存新快照，不修改旧文件。不同分辨率必须重新绘制，不自动缩放。
11. **查看主要结果**：英文状态区显示 Status、Candidates、Velocity 和 Trend；底部显示录像状态，有漏录或写盘错误时提示。CPU/GPU/内存、内部耗时、后端名称和开发用计数不再显示，默认也不记录资源与性能采样。逐帧时间、识别对应关系和数据质量仍保存在实验文件中。显示流畅不意味着每一帧都完成了识别，录像也不等于相机内部高速原片。
12. **看观测层**：Candidates、PIV、Mask、Regions、Shape、Curves 可独立开关；**Regions 与 Shape 默认不勾选**。默认显示最新接收画面，识别叠加在有效期内保持固定亮度，`Overlay delay` 固定显示时间差。已取消周期性渐隐，到期仍撤下。勾选 Measured frame 可查看对应的已分析帧。PIV 是细线开放式箭头，带 25 px/s 标尺。**Smooth PIV 默认开启**，在固定网格上平滑显示，图例标记 `smoothed`；它有短时响应延迟，右侧数值与原始记录不做平滑。需要逐次检查原始结果时取消勾选，图例变为 `raw`。过时、无效或缺少当前支持的运动不会通过平滑持续显示。当前速度为 px/s，物理尺度需同光路测微尺标定。
13. **自动录像与实验事件**：点击 **Start** 自动开始接收并录制实验图像，**Stop** 自动完成保存，窗口底部显示录像状态。**Open records** 在 Windows 文件资源管理器中直接打开本次 `session_*`：`experiment_raw.mp4` 是无标记的接收画面，`experiment_overlay.mp4` 是带当前识别图层和时间/帧号标注的画面，均不包含软件窗口。逐帧时间与识别对应关系保存在 `video_frames.jsonl`，录制完整性见 `video_manifest.json`。没有收到画面时不会产生视频；录制失败或漏录会显示，不会伪称完整。录制的是相机预览流，不等于相机内部高速原片，UTC 也是接收时间而非曝光时间。手工操作 PZT 时把事件源选 MANUAL_ANNOTATION，选通道、填写条件并点 Mark ON/OFF，按钮只记事件。新 mask 快照在会话上级的 `geometry/live_mask_*.json`，变更日志在会话内的 `geometry_changes.json`；下次启动按当前视野重新绘制或主动加载。保留真实设备日志和焦点/注射条件。
14. **断流恢复与安全停止**：断流会使状态无效，原接收器尝试重连；每约120s会有有限桥接会话轮换，短暂停顿正常记录。需要时点 Reconnect；先 Stop 再关闭窗口，等待线程与进程回收。异常关闭也执行清理。不要删除超时以保持“永久运行”。
15. **无画面排查与实机验收**：检查相机电源/既有网线/Windows 原预览是否正常、UDP端口是否被旧窗口占用、WSL是否能运行 Windows Python、路径是否匹配、windows_stderr.log 与 bridge_status.json。必要时使用原 check_connection.py 只读诊断，但本轮未主动连接相机。实机必须补做现场新鲜度、分辨率/遮罩、真实输入速率、现场到显示光学延迟、长时重连、静止标定及人工 ON/OFF 对照；未通过前不开放闭环控制。

路径变化可用 `--windows-project`、`--windows-python`、`--live-geometry`，或独立 JSON 配置，毋须编辑算法源文件。默认服务只在 localhost 发布，不开公网无认证图像端口。
