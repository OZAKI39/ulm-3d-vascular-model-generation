# 真实录像回放

源SHA `3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c`；1476×1120，937帧，PTS 0–31.1688s，步长0.0333s，呈现时长31.2021s，容器约31.232s。完整核查见 results/software_validation/source_audit.json。没有修改原视频。

- replay_1x: {'bridge_received': 937, 'processed': 693, 'state_published': 619, 'fresh_valid': 689}；跳过 244。
- stress_2x: {'bridge_received': 937, 'state_published': 310, 'processed': 405, 'fresh_valid': 401}；跳过 532。
- debug_half: {'bridge_received': 937, 'state_published': 1240, 'processed': 937, 'fresh_valid': 933}；跳过 0。

1×真实节奏；2×只是约60.06Hz软件加速输入；0.5×约15.02Hz调试。循环会话重建、暂停/恢复、EOF与停止已有实际测试，不复制帧冒充新运动。
