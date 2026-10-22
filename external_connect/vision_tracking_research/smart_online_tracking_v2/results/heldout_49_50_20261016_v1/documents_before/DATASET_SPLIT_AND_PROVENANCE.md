# 数据划分与来源

> **最新续研：2026-10-16，已完成离线人工ROI评价。** 以下较早记录按历史保留；当前状态以文末续研和[完整报告](results/supervised_20261016_v1/REPORT.md)为准。

原始录像：`/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4`。SHA-256：`3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c`。原始图像 1476×1120，录像实际 PTS 约 30.03 fps；不是 60 fps 实时流，更不代表传感器采集帧率。没有物镜倍率、适配器、焦平面、像素到 µm 标定，保持 null。

72 张录像候选：train 帧 0–270、validation 360–570、test 660–936，各 24 张；每组含一个连续 8 帧短片段（120–127、450–457、780–787）。组间保留时间间隙，不随机拆分相邻帧。完整候选 manifest 保存原始 PTS / time_base、图像和视频 SHA、尺寸、ROI、类别审核状态。

另 12 张原始 1280×1024 实时帧约每秒采样一次，共 84 张。直播接收器没有暴露传感器 PTS，因此 source_pts=null；另存 Windows host decode time / UTC、相机配置、session、原始 PNG SHA。source_sha 在直播帧中指 PNG 文件，不冒称原始流文件哈希。直播数据全部 QUARANTINE_UNASSIGNED，未经新分组和审核不得训练或评价。现有样本是否包含真实团、突然转向、不同焦平面等需要专家确认；不能把 12 张相近直播画面当作 12 种独立条件。

当前 1 张正式确认、1 张草稿。原 20/21 正例匹配没有混入真值。保存整个源图，任何 ROI 推理另记 xywh，未改变检测分辨率；演示视频只在显示时放大 2 倍。

同视频内部初步划分不证明跨实验泛化：CROSS_SESSION_GENERALIZATION_NOT_VERIFIED。下一次采集应补独立注射轮次 / 实验日、清晰单泡与团、焦平面变化、亮度变化、边缘进出、急转弯和聚集 / 分离片段，并记录真实光学配置。

## 人工完成后继续

先运行 `.venv/bin/python -B scripts/review_labels.py` 核查实际文件、SHA、完整审核和类别冲突。用户完成后可用 `--freeze-name <唯一名字> --export-train` 固定所确认修订并导出 train / validation 的标准 Cellpose PNG + masks TIFF。测试集不会导出到训练目录；含 UNCERTAIN / OUT_OF_FOCUS 的帧默认排除训练，防止标准 Cellpose loss 把未知区域误当背景。未确认帧一律不导出。

评估规则已在 protocol.json 预登记：单泡 6 px 中心匹配、团实例 IoU 0.5，Hungarian 一对一；不确定 / 失焦区域显式忽略，预测与忽略区重叠超过 0.5 时拒判。阈值仅可在训练 / 验证组选择。当前只有部分训练帧确认，验证/测试真值仍不足，尚未运行检测或短轨迹准确率评价，也未据测试结果调参。


## 2026-10-16 已有标注的离线续研（supervised_20261016_v1）

新冻结集vision2_supervised_20261016_v1含4张源图：train90/108（210泡）、validation360/374（133泡），均为区域审核。test未使用。5个验证ROI只算2张源图；模型辅助后人工审核，非盲法双人真值。原第1–5张仅部分标签，不作为全图真值。

完整参数、失败证据、28项回答和输出路径见[本轮报告](results/supervised_20261016_v1/REPORT.md)。
