# F10S 历史基准核查

从最新 F10Y 交付 647eccd83d439f4dc9600097f9b19444f3cbc70c 建立隔离 worktree；继承的科学 F10S 交付为 31a44d3f342e5d0bca55f7f9b83f4d32b8189696。先检查服务器源代码与 /etc/vast-agents-guide.md，未用公开旧版本覆盖服务器。
11 条旧轨迹逐一核验 SHA-256、时间、初态和 ID；详情 evidence/F10S_BASELINE_AUDIT.json。原 3 ms 双源独立重算 X 位移 0.197999446649 µm；dt/cache 通过。F10Z 结果没有改名为旧结果。

32 个固定 ID，R=1.78 µm，2.5 MHz；背景流严格为零；继承平移阻力黏度 0.001 Pa·s、径向响应有效黏度 0.002 Pa·s、密度 1000 kg/m³、声速 1500 m/s。这是原模型参数，未宣称为本实验纯水校准。

保留 Sequential QP 原 active_set 求解；OSQP 1.0.5 可用但未切换后端。硬球几何限 1 nm，原 overlap budget 17.8 nm，接触修正 178 nm，mobility 1e-13 m，实际互补 1e-11 m；墙修正必须为零。连续量比较保留 COM/Rg 20 nm、固定 ID 50 nm、Q/长宽比/XY 1%、分布 RMS/R 0.02、可辨识主轴 2°。相位/压力非法时拒绝，没有裁剪。

12 个受保护源码/配置与冻结缓存、APECSS 表哈希在执行前和最终审查再次核验。服务器执行快照、独立 Supervisor 标签与环境均保存在 logs/ 和 evidence/。
