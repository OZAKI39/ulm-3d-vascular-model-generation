# 单泡与团级准入

> **最新：第49、50张留出测试已完成。** [本次结果](results/heldout_49_50_20261016_v1/REPORT.md)；此前续研记录保留。

> **最新续研：2026-10-16，已完成离线人工ROI评价。** 以下较早记录按历史保留；当前状态以文末续研和[完整报告](results/supervised_20261016_v1/REPORT.md)为准。

| 层次 | 当前状态 | 依据与限制 |
|---|---|---|
| 实时相机入口 | LIVE_CAMERA_PASS | 原 Windows 接收器与 WSL 有界桥，实际两条20秒试验 |
| 人工工具 | ANNOTATION_TOOL_PASS | 修复后实际 Qt 多对象、保存恢复、换帧、SAM 后台与取消 |
| 专家标签 | EXPERT_GROUND_TRUTH_PENDING | 详见FINAL_GATE.json的带时间快照；验证/测试标签尚缺 |
| 单泡识别 | NOT_EVALUATED | score 未校准；面积分类不是语义真值 |
| 短轨迹 | NOT_EVALUATED | Norfair软件关联可用；无人工连续ID真值，不声称改善断裂 |
| 团轮廓 / 形状 | GROUP_SHAPE_NOT_ADMITTED | group_shape=null；原始图像区域形态单独保存，不进正式接口 |
| 团附近单泡 | NOT_ADMITTED | 有像素邻域几何辅助函数；没有人工准入团，因此正式输出null |
| 并入 / 离开事件 | NOT_EVALUATED | 当前不产生已验证事件；未来只能称APPARENT_INCORPORATION_CANDIDATE，不能推断物理融合 |
| 软件实时 | OpenCV ROI throughput PASS | 并非所有帧低于16.7ms，也非身份准确率PASS |
| 曝光端到端延迟 | NOT_MEASURED | 无LED测量；曝光时间未知 |
| 硬件控制 | DISABLED | 仅本地ZeroMQ mock；无设备控制API |

单泡位置采用当前观测；丢失不输出预测旧坐标。速度为最近观测的像素/源时间秒，直播源时间为host decode时间，不能解释为标定物理速度。LK为局部图像运动而非液体流速。跨session、逆向时间、过大帧间隔重置轨迹。原始mask的质心、面积、Rg、Q和长短比通过独立矩形软件检查，但形状公式正确不证明其输入是真实泡团。

当前训练和科学控制状态：training_ready=false。后续需冻结完整人工标签、用未调参测试组验证单泡和团分割、用短序列检查真实关联，才考虑解除相应语义门控。本轮不会自动将模型预测转成专家标签。


## 2026-10-16 已有标注的离线续研（supervised_20261016_v1）

EXPERT_LABELS_READY仅适用于2 train+2 validation源图的局部单泡掩码。BUBBLE_DETECTION_VALIDATED、SHORT_TRACKING_VALIDATED、GROUP_SEGMENTATION_VALIDATED仍false；没有test、团或短ID真值。group_shape=null，control_permission=false，training_ready=false。

完整参数、失败证据、28项回答和输出路径见[本轮报告](results/supervised_20261016_v1/REPORT.md)。


## 第49、50张留出测试（2026-10-16）

已有train/validation/test各2张源图的局部单泡真值，HELDOUT_ROI_EVALUATION_COMPLETE=true；科学检测、短轨迹、团形态准入仍false。掩码重合度低、场景有限，不能提升为任意实时识别或控制许可。
