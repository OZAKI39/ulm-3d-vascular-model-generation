# TorchPIV GPU 实测

RTX 4090 实际调用官方 PIVbackend.extended_search_area_piv 及 CWS 64→32，多遍和单遍均运行。Headless namespace 仅绕过 GUI 初始化，原相关代码不改。TorchPIV 本地公开接口只给原生 SNR flags，未虚构原始 SNR 数值；最终再做共同局部中值/纹理与六参数检查。

| Backend | ROI p50/p95 ms | Medium p50/p95 ms | Full bbox p50/p95 ms |
|---|---|---|---|
| B0 | 23.04 / 26.84 | 79.39 / 85.72 | 194.37 / 215.73 |
| B1 | 14.88 / 18.27 | 71.29 / 78.06 | 174.16 / 187.73 |
| B2 | 6.54 / 7.65 | 22.88 / 25.27 | 56.15 / 63.18 |
| B3 | 5.98 / 7.82 | 12.27 / 13.00 | 18.73 / 21.95 |
| B4 | 5.89 / 7.23 | 9.02 / 11.61 | 14.79 / 16.39 |


TorchPIV 适合规则大批量；稀疏窗拼装在试验中反而较慢。当前 dtype 经上游归一化为 float64，不能声称这是 float32 对等 benchmark。上传和同步由适配器计时，下载融合在官方调用内部，单独下载归因 NOT_SEPARATELY_MEASURED。全流程计时包含此开销。
