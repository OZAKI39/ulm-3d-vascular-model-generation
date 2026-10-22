# 固定时间 population 实验

主报告：[FIXED_WINDOW_PROTOCOL_REPORT_ZH.md](FIXED_WINDOW_PROTOCOL_REPORT_ZH.md)。

正式合同：[FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json](FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json)。

T_infusion=1.240815660 s；停灌后观察10.0 s；T_end=11.240815660 s；出生区间[0,T_infusion)。停止规则不读取出口覆盖结果。0.5 ms是正式名义步长；少量合成测试使用0.37 s以验证非整除时间窗。

唯一正式入口：`python -m particle_3d.population_v1.run_fixed_window_population`，参数及服务器复现命令见主报告。当前兼容入口`population_v1.run`同样强制固定合同。没有运行正式大批量实验，没有CFD计算。

源提议ledger与准入ledger分开。`RIGHT_CENSORED_OBSERVATION_END`不是失败或出口；accepted births及全部终态始终保留。

迁移来源：旧production无显式固定灌注时长，因此按用户要求固化前次1.240815660 s。旧history的精确停止时间1.2408156603855842和停灌后12s观察保留原样，不能改称新10s协议。旧结果hash不变性见validation/history_immutability.json。

CPU/GPU只做已有模型执行一致性检查；不替代近壁CFD审计或物理精度验证。
