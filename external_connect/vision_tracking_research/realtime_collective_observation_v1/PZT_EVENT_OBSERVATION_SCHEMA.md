# PZT 事件：只记录

schema 字段：event_id、channel（ALL/X±/Y±）、state（ON/OFF）、source（MANUAL_ANNOTATION/DEVICE_TELEMETRY/SIMULATED_EVENT）、command_time、actual_execution_time、camera_time、source_pts、session_id、host_event_monotonic、time_sync_quality、conditions、note、hardware_control=false。

当前 UI 提供模拟与人工标注两类按钮，Mark ON/OFF 只写 JSON 并在曲线中标记。执行时间和曝光相机时间保持 null。人工点击不证明声场实际输出；导入设备证据应保存只读遥测原件，不能把人工时间改称实际执行时间。

`--simulate-events` 仅用于演示：5.1 s 与 18.1 s 定时标记 SIMULATED_EVENT。录像没有真实 PZT 日志，事件前后变化不能被归因为声学。断流时事件仍可记录，但观测窗口无效，需要排除。

将来比较 ON/OFF 应按相同 session 与明确同步误差截取有效窗口，保留条件、焦点、注射、流动、曝光等协变量，先做预登记对照。当前不运行 AFG SCPI、不设置电压、不调用设备输出。

实际比较入口：`.venv/bin/python scripts/event_window_analysis.py <session目录> --output results/event_windows.json`。每个事件按相同session、源PTS划分前后2s，仅汇总有有效运动的样本；没有可信时间映射或不足8个样本标为NOT_OBSERVABLE。模拟事件始终输出SOFTWARE_SEGMENTATION_ONLY，不给声学因果结论。
