# 历史基线与冻结输入

F10S/F10Z/scaling 共 126 条历史轨迹逐文件对照保存 SHA；详见 `../BASELINE_AUDIT.json`。Git 原工作树清洁，新工作树从实际 F10Z 最终提交建立，未从旧公开 F10S 分支覆盖。

冻结 16 个初态、32 个 ON/OFF 力学状态：32/64/128 的分散初态、历史终态、合法单对/立方/FCC，以及原 F10S 32 泡近接触群。立方与 FCC 的压缩力 F=ζA(-20x) 是制造解的力学检查，**不是声学形成轨迹**。生成规则、种子、输入哈希先于优化保存。原 REFERENCE.h5 保存稠密 A/B、法向/间隙、声学力、速度、完整接触反力和迭代记录。

128 泡旧 1 ms 仍无活跃硬球接触，不能用于单独验收复杂接触；本轮额外覆盖 532 个同时 FCC 活跃约束、64 泡真实接触延续、32 原近接触轨迹和延长的 128 泡真实声学轨迹。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
