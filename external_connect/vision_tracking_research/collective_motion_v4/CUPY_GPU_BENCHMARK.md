# CuPy GPU-PIV 实测

已执行官方 test0、48/24 单遍、64→32 多遍、整通道与稀疏 mask。设备 RTX 4090，CuPy 13.6.0，实际可运行环境版本锁在 evidence。48-px 窗口内部 FFT padded 到 64，与 SciPy/Rust circular48 并非数值等价算法。

| Backend | ROI p50/p95 ms | Medium p50/p95 ms | Full bbox p50/p95 ms |
|---|---|---|---|
| B0 | 23.04 / 26.84 | 79.39 / 85.72 | 194.37 / 215.73 |
| B1 | 14.88 / 18.27 | 71.29 / 78.06 | 174.16 / 187.73 |
| B2 | 6.54 / 7.65 | 22.88 / 25.27 | 56.15 / 63.18 |
| B3 | 5.98 / 7.82 | 12.27 / 13.00 | 18.73 / 21.95 |
| B4 | 5.89 / 7.23 | 9.02 / 11.61 | 14.79 / 16.39 |


mask 的 int64/int32 kernel ABI 问题修复前导致 illegal address；原 pilot_v1 的极短时间伴随后端失败，不计入性能收益。修复后共同选窗位移完全一致。缓存以 shape/window/overlap 为键，复用官方 solver 和 GPU buffer。

独立传输微基准：两幅1120×912 float32共 8171520 bytes，上传 p95 1.669 ms，下载 p95 1.325 ms。它是同步额外 probe，不是从总时间扣减或对融合官方 API 逐项归因。实际显存/GPU利用率见 PERFORMANCE_AND_RESOURCES。
