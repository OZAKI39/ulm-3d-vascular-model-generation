# v2 录像的 PIV 群体运动可视化

正式结果采用两个 `final/` 目录中的视频。青色小圆圈沿用已通过视频的 OpenCV 候选，黄色箭头表示 PIV 图像运动估计，原有时间戳、Arial 字体和 PZT 面板保留。

| 输入 | 正式视频 | 原始帧数 | 时长 | 代表图 |
|---|---|---:|---:|---|
| experiment_overlay_v1 | [PIV 视频 1](experiment_overlay_v1/final/experiment_overlay_v1_opencv_piv.mp4) | 293 | 7.6463643 s | [PNG](experiment_overlay_v1/final/preview_00146.png) |
| experiment_overlay_v2 | [PIV 视频 2](experiment_overlay_v2/final/experiment_overlay_v2_opencv_piv.mp4) | 3588 | 86.3961843 s | [PNG](experiment_overlay_v2/final/preview_01794.png) |

PIV 在原录像的相机像素上计算，排除顶部文字和右侧 PZT，不在新增的圆圈或箭头上做相关。保存的检测中心用于判断窗口是否有候选支持，不承担位移测量。未重新运行或更改 OpenCV 检测。

方法为已有 F10X B4 OpenPIV GPU FFT 互相关、Gaussian 亚像素峰拟合、96×96 px 窗口、48 px 重叠。使用当前帧与前三帧的源 PTS 实际间隔，得到该时间段内的平均表观速度。预处理为 sigma=4 px 的局部暗特征增强。质量检查包括 peak-to-peak SNR≥1.2、纹理标准差≥2、2.5 px 局部中值容差、位移≤16 px，并要求前后窗口均含至少一个候选。无效矢量不填补、不外推。

速度单位是 px/s，图像 x 向右、y 向下。没有这批图像的人工确认管腔 mask、显微像素物理标尺或固定管壁配准，因此可能包含相机/载台运动；PZT 示意图的 5 mm 标尺不能用来换算显微画面。当前结果是群体图像运动估计，不能当作校准的真实血流或单泡三维速度。

为便于查看慢运动，最终箭头采用 24 px/(px/s) 的显示比例，参考箭头为 1 px/s，最大长度 48 px（对应 2 px/s）；超过该值只限制绘图长度，原始速度不裁剪。低于 0.25 px/s 的显示矢量不画箭头，避免很短的箭头变成噪点。72 px 固定显示网格和 120 ms 因果平滑只影响显示；当前缺少有效结果的网格不会延用旧箭头。短视频有 3 帧预热；长视频有 5 帧预热和 1 帧 PTS 间隔过大，未估计这些帧的 PIV，视频本身仍保留全部帧。

数值检查：整数位移(3,-2) px 的误差中位数 0.0328 px、95 分位 0.2039 px；半像素位移(-2.5,1.5) px 分别为 0.1249 px、0.4561 px。静止测试位移为零，空白图不输出有效矢量；一对真实图像的 CPU/GPU 共同有效窗口差异中位数 0.000785 px、95 分位 0.00310 px。这些是数值核查，不是实验运动真值。附加 0.125 px 位移测试显示中位估计只有约 0.053–0.063 px，存在小位移低估，不能据此宣称极低速度的精确定量。详见 `preflight.json`、`low_displacement_sensitivity.json`。

窗口比较和失败记录保留：64 px 窗口未通过最初半像素误差门槛，改用满足同一门槛的最小测试窗口 96 px；GPU field-mask 路径发生非法内存访问，正式计算使用通过测试的全网格路径，再筛选有候选支持的结果。未改写上游历史代码或放宽质量阈值。

每个 `experiment_overlay_v*/` 目录保存：

- `piv_fields.npz`：原始位移、有效性、支持窗口、SNR、时间戳、配对帧和实际 dt；原始速度可由位移/dt 得到，只使用 valid=true 的位置。
- `frame_statistics.csv`：逐帧有效窗口数、速度统计、运行耗时及无结果原因。
- `manifest.json`、`validation.json`、`source_snapshot/`：输入/代码/结果哈希、Git 身份、参数、时序检查与源代码快照。
- `final/display_vectors.jsonl.gz`：最终显示矢量，与原始数值分开保存。
- `final/manifest.json`：正式 MP4 的 SHA-256、逐帧时间戳误差、PZT 绘图保护检查和数据未修改检查。

正式计算脚本：`smart_online_tracking_v2/scripts/render_piv_video.py`；最终显示脚本：`smart_online_tracking_v2/scripts/render_piv_display.py`。复用 `realtime_collective_observation_v1/.venv/bin/python`。本次没有修改实时界面默认设置、仪器输出、人工标签或训练模型。
