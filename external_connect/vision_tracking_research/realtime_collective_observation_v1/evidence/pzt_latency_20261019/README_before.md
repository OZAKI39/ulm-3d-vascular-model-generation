# F10Y · Chronos 实时群体观测

这是一个可启动的 PySide6 实验观测程序。默认输入为真实 Chronos 录像，按原始 PTS 1× 接收；LIVE 接口复用原 Windows 接收器。本轮没有连接真实相机，也没有控制 PZT/AFG。

## 在家一条命令启动

```bash
bash /home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1/start_replay.sh
```

界面显示完整画面、通道遮罩、候选泡、有效 PIV、候选区域与连续趋势。播放结束会显示 EOF，旧运动立即失效。点击 **Reconnect** 从头开始；**Stop → Start** 也会新建会话。Pause/Resume 仅作用于录像。

一键脚本在本机已安装的 RTX 4060/CuPy 环境使用 B4；未配置 GPU runtime 时使用 Rust OpenPIV（B2）。远程 RTX 4090 的正式验证也使用同一 CuPy（B4）。两者共享检测和趋势代码，但 PIV 实现不数值等价，不能把服务器吞吐套到笔记本。显式 `--backend B2` 可使用 Rust 后备；普通 SciPy（B1）只作为慢速诊断备用。首次 CUDA 编译会暂时降低观测更新，原图仍持续显示；本机实测性能单列于报告。

实验室启动：`bash start_live.sh`，选择 LIVE_CAMERA 后点击 Start。点击 **Draw LIVE mask** 沿通道内部描点，**双击、Enter 或点击起点完成**，区域立即生效并自动保存，无需保存按钮或编辑 JSON。**Clear current mask** 清除当前区域和旧分析结果，历史记录保留；`Backspace` 撤回一点，`Esc` 取消本次绘制。可选 **Load saved mask** 复用同分辨率的旧区域。更改区域不会重连相机。详细步骤见 [LAB_CAMERA_INTEGRATION_GUIDE.md](LAB_CAMERA_INTEGRATION_GUIDE.md)。

新绘制的区域默认显示**图像坐标下的运动**（`IMAGE_SPACE_UNREGISTERED`），不要求固定壁面参考，也不会声称已校正相机晃动；画面会明确提示。已有固定壁面参考的配置仍执行参考质量检查。两种模式都不是已校准的真实速度或闭环控制依据。WSL 的 **Open records** 直接调用 Windows 文件资源管理器；无需 Linux 桌面 launcher。

## Start → Stop 自动实验录像

无需另开录像开关。点击 **Start** 后，收到第一帧就自动录制；点击 **Stop** 或正常关闭窗口会完成视频收尾。底部状态栏显示录制/保存状态，有漏录或错误时明确提示。**Open records** 直接打开本次会话目录，包含：

- `experiment_raw.mp4`：原始接收画面，无软件界面或识别标记（H.264 压缩，不是无损传感器数据）。
- `experiment_overlay.mp4`：同一批画面，加上当前勾选的候选泡、PIV、mask、群体范围/形状，以及 UTC、实验时间、候选数和必要的叠加时间差提示。完整帧对应关系在逐帧记录中；画面没有侧栏、按钮或曲线面板。
- `video_frames.jsonl`：每个录像帧对应的接收 UTC、单调时间、源帧号、源时间戳（如有）、MP4 PTS、识别帧号、对应方式和 mask 版本。
- `video_manifest.json`：视频文件、编码信息、帧数、漏录、录制状态与时间基准说明。
- `states.jsonl.gz`：开启自动视频时保存科研观测数据，含逐帧候选位置和有效速度向量；另有事件、mask 变更和运行汇总。

每次 Start 建立新的 `session_*`，不会覆盖上次视频。无图像时不制造空视频。分辨率变化会自动分段为 `_part002.mp4` 等，原始像素不缩放，奇数边长仅在右/下补一像素。再次按 Stop 不会重复写坏文件。

录像编码队列最多缓存 48 帧，用于吸收启动时 GPU 初始化造成的短暂积压；这不会延长检测有效期，也不阻塞相机接收。持续编码不足仍会明确记录漏帧并提示，不能把缓冲当作帧率保证。

录像使用本机实际接收时间间隔，不强行按固定帧率加速或补帧。识别结果优先匹配同一帧；当分析较慢时，仅使用仍在时效内的过去结果，固定位置标注 `Overlay delay` 和时间差（数据记录为 `RECENT_PAST_FRAME`），不把它称为精确同帧识别。过时、未来帧、其他会话或其他 mask 的结果不叠加。图层开关会影响叠加视频；`Measured frame` 仅影响监看画面，录像仍保留接收到的图像序列。原始视频始终不加标记。

时间戳是**预览帧接收时间**，不是未经提供的相机曝光时间；录像记录的是收到的预览流，不等于相机内部的高速原片。暂停回放/断流产生的时间间隔保留，不填造运动。实验录像与旧 `--capture-ui` 软件界面录屏是独立功能，普通 Start 不会录制软件窗口。

## English 界面与 PIV 箭头

界面按钮、说明、错误提示和录像状态统一为英文。主状态区仅显示观测状态、候选数、平均速度向量（px/s）与聚集趋势；底部三条曲线为候选数、群体回转半径和径向运动。LIVE 模式隐藏回放速度、循环与暂停控件，PZT 事件默认采用人工记录。

**Regions 和 Shape 默认关闭**，需要时可分别勾选。候选和箭头在有效期内保持固定亮度，提示条保持固定位置。已撤销旧版最后 80 ms 渐隐：实际更新间隔接近 TTL 时，这种处理反而会形成周期性的明暗脉冲。过期、断流和无效运动仍撤下，不延长原有有效时间。`Measured frame` 的图像和检测状态从同一个锁定快照读取，避免更新时错配。

检测图层按新结果批量绘制并缓存为透明图像，后续画面刷新直接叠加，减少逐泡重画对界面和录像的干扰。缓存随检测、图层、窗口尺寸和有效状态更新，不积累旧位置，也不会延续过期结果。

PIV 改为亮黄色细线和开放式箭头端，使用窄暗色衬边，不再画粗黑边和实心三角。界面线宽约 1.35 px；界面与视频共用箭头几何。长度采用 `0.8 s × 显示速度 px/s`，带 `25 px/s` 标尺，超过 40 px 限长并提示。小于 3 px 的箭头不画，避免近零速度被画成杂点。

**Smooth PIV 默认开启**：在固定 72 px 网格内取当前有效速度的局部均值，再使用 0.12 s 时间常数的因果平滑，减少箭头起点换位和长度/方向噪声。图例标记 `smoothed`；这只是显示处理，有响应延迟，不是额外测量或逐泡轨迹。无当前支持的网格立即撤下；停流、反向、无效结果、超时、重新连接和区域改变不会沿用旧方向。关闭 **Smooth PIV** 可检查 `raw` 箭头，每个 48 px 网格选取一个原始实测向量。右侧速度、科学统计和 `states.jsonl.gz` 中的原始速度向量均不做此平滑；录像逐帧记录显示开关，manifest 保存显示方法，复现时可从完整原始状态顺序重建显示场。

默认不采集或保存 CPU/GPU/RSS 资源采样，不写入内部耗时表、UI 绘制耗时、自动窗口截图或重复的无效窗口网格。视频、时间戳、候选位置、有效速度、mask、事件、质量状态、漏录与错误记录继续保存。`diagnostics: true` 仅供显式工程排查；此次调整不删除历史记录。

## 文件与复现

- `run_replay.py`、`run_live.py`：同一个 Qt 界面和核心。
- `run_headless.py`：无界面有限或循环运行，支持 `--seconds`、`--speed`、`--loop`、`--backend`、`--video`、`--output`。
- `.venv/bin/python run_tests.py`：软件与 Qt 操作测试，不连接相机。
- `configs/default.json`：路径、TTL、线程、回放配置。
- `results/`：实算数据；`figures/`：PNG/PDF；`videos/`：实录与科研视频。
- `runtime/<时间>/session_*`：用户每次操作的独立记录，停止后写 summary。源视频、原标注、F10U–F10X 均未改写。
- [F10Y_DEVELOPMENT_REPORT.md](F10Y_DEVELOPMENT_REPORT.md)：最终结果、全部 35 项问题。

所有速度均为图像像素单位。二维候选靠近、密度增长、区域连通都不是物理融合或三维聚集证据。`TRAINING_READY=false`，`SCIENTIFIC_CONTROL_ADMITTED=false`。

## PZT 状态与实验模型（2026-10-19）

右侧 **PZT status / layout** 使用当前 F10M-R3 的固定装置模型，显示 PDMS、中央通道、两个顶部孔及 X±、Y± 四枚 7 mm PZT。点击 **PZT connections**，只需确认 A/CH1、A/CH2、B/CH1、B/CH2 分别接到哪枚实体 PZT。实际连接与读回参数自动按对应关系显示，不再手填几何坐标；重复分配被拒绝，未确认的关系不猜测。操作说明见 [PZT_MONITOR_GUIDE.md](PZT_MONITOR_GUIDE.md)。

**Read PZT status** 可直接发起四路只读查询并等待完成回执；背景已统一为黑色。每次启动自动恢复最后保存的 PZT 对应关系。详细说明见上述指南。
