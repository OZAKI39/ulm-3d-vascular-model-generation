# 活动代码索引

| 文件 | 职责 |
|---|---|
| src/backends.py | B0 旧实现、B1 SciPy、B2 Rust、B3 TorchPIV、B4 CuPy；现成内核薄适配 |
| src/channel.py | 可编辑几何、全窗口覆盖、每帧全通道候选/纹理/活动及配准 |
| src/observer.py | 因果背景、调度、窗口 null/年龄、局部候选群六参数 |
| src/fit_adapter.py | 保留 F10U 拟合主体，仅提供每次复制的解析雅可比 |
| src/channel_state.py | 新接口检查、TTL、显式局部群 v2/v1 兼容 |
| src/studies.py | 官方多遍调用、已知位移与六参数误差 |
| scripts/setup_remote.sh, build_rust.sh | 独立环境、官方 Rust 编译 |
| scripts/cupy_mask_compat.py | int64/int32 窗口压缩索引兼容修复与 diff |
| scripts/backend_benchmark.py | 相同帧对三种尺度的同步性能比较 |
| scripts/frozen_baseline.py | 原 F10W 代码及原环境的额外基准 |
| scripts/generate_data.py, accuracy_study.py, freeze_and_final.py | 新程序种子、开发/验证/一次性 final |
| scripts/audit_channel.py | 只读已审核人工区域、全帧扫描及位移拒判检查 |
| scripts/full_offline.py, causal_replay.py | 全 937 帧离线和独立 PTS 定速回放 |
| scripts/diagnostics_and_clip.py | 60 帧五后端对照、筛选/墙边/突变/传输诊断 |
| scripts/analyze_delivery.py | 正式结果汇总、逐帧数值等价、面积与不确定性 |
| scripts/publication_figures.py, render_videos.py, render_backend_video.py | 保存数据驱动的 PNG/PDF/MP4 |
| tests/test_contracts.py | 数学、几何、null、TTL、兼容及鲁棒雅可比回归 |

历史 F10U 数学、F10V PyAV 编码、F10W 旧环境以只读方式复用。third_party 是被 Git 忽略的官方检出；许可证与 commit 在 evidence 中。
