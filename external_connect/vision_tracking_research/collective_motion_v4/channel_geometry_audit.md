# 通道遮罩审查记录

原始首帧、中帧、末帧在 evidence/channel_overlay_0000.png、0450.png、0936.png，最终 mask 为 evidence/channel_mask.png。人工部分标签只读保存，未将未审核区域当背景。中央条带曾排除真实泡，已经修正；原版仍保留。

当前多边形是基于真实影像的可编辑初始分割，不是物理边界真值。whole-window 100% 内含、配准拒判及遮罩哈希均进入状态；20 px 边界不确定性及固定纹理残留均限制其适用性。视觉检查详细清单见 evidence/visual_inspection.json。
