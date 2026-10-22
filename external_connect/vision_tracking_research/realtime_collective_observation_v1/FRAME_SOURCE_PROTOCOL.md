# FramePacket 与时钟

两种源都进入 Pipeline.receive → 最新槽 → Continuity → 同一 Observer。必备字段见 src/protocol.py。GRAY8/BGR8、uint8、原始分辨率、ROI 与 image.shape 必须一致；帧号非负；PTS 有限非负；空图像/错包拒绝。

REALTIME_REPLAY：source_pts 是容器里的真实录制时间轴；host_decode_timestamp 是本机解码完成；arrival_monotonic 是按 PTS 调度后的释放/接收时间。回放 age 从释放时刻算，并另报 decode_ms，不能将等待播放时刻称为计算耗时或曝光延迟。1× 为原始节奏；2× 只是压力输入，0.5× 为调试。

LIVE：capture_timestamp=null，EXPOSURE_TIMESTAMP_UNKNOWN；Windows 解码 monotonic 只作为带时钟标签的原始证据。到 WSL 的年龄上界 = Windows sender age + WSL GET RTT + WSL 本地等待；从不直接相减跨主机 monotonic。LIVE PIV dt 使用本地接收间隔近似，明确标注，不冒充硬件曝光间隔。

同会话重复/倒序帧拒绝；会话变更、PTS 回退或超过 250 ms 间断会清空运动历史。循环录像产生新 UUID，从 PTS 0 重新开始，不能跨循环算速度。

解码数、发送数、桥接接收数、桥接覆盖、分析跳过、过期拒绝、处理数、有效数、UI 跳过及异步记录丢弃分别记录。它们不是网络丢包统计。
