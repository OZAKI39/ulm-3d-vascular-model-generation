# 旧数据为何看起来先停后走

原生末态到输运初态的组件 COM 位置、速度误差均为 0；TRUE_VELOCITY_RESET_BUG=NOT_OBSERVED。
207 个 1–2 节点组件的初始速度范围 0.0952942–1.5685 m/s，中位数 0.858503 m/s；对应无栓 CFD 速度范围 0.000319625–0.00987695 m/s，中位数 0.00458474 m/s。

旧模型在 2.33 ms 前背景拖曳为零，交接后突然启动 tau=100 µs 的工程松弛；冻结局部流速时滑移速度按 exp(-t/tau) 衰减，500 µs 后仅剩约 0.00674。这是快速减速，不能误称速度清零。旧接受步数据可复现到 2.255e-16 m/s。

近壁无滑移使局部背景速度低；19 个组件初始背景流速小于 1 mm/s；208 个组件在 100 µs 内按背景速度移动都小于旧近景一像素。旧播放率从每视频秒 0.0004 s 物理时间跳到 0.8 s，相差 2000 倍，且漏过快速松弛过程；更早视频还有 0.05 s 重复端点（上一轮已修复）。

另一个真实限制是附着主体在 2.33 ms 后仅显示冻结参考快照，其原生末态 COM 速度约 0.0547157 m/s；这不意味着主体物理静止。

来源：/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance/verification/y_clot_continuous_handoff_and_honest_render_v1/HANDOFF_CONTINUITY_AUDIT.json、handoff_velocity_audit.csv、VIDEO_TIMEWARP_AUDIT.json。本轮只读复用旧审计，不重算 baseline。
