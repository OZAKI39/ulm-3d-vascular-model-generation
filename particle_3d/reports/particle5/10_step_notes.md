### 图 10：10_real_fem_two_mb_resistance_smoke

![图 10](figures/10_real_fem_two_mb_resistance_smoke.png)

应该看什么：看原 P4 双球初值、尺寸和时长下的新轨迹，以及每一步的近场资格。

实际看到什么：完整烟雾回放状态有限且无穿透；墙面近场启用，双球间隙比仍远大于 0.01，pair 润滑没有启用。

有没有异常：新轨迹产生更小墙隙是积分结果，不是改初值；验证窗口没有出口事件，也不代表完整悬浮液。

数据：[10_initialization.json](data/10_initialization.json)、[10_initialization.csv](data/10_initialization.csv)、[10_real_summary.json](data/10_real_summary.json)、[10_real_summary.csv](data/10_real_summary.csv)、[10_real_two_mb_states.json](data/10_real_two_mb_states.json)、[10_real_two_mb_states.csv](data/10_real_two_mb_states.csv)、[10_pair_eligibility.json](data/10_pair_eligibility.json)、[10_pair_eligibility.csv](data/10_pair_eligibility.csv)、[09_real_near_wall.json](data/09_real_near_wall.json)、[09_real_near_wall.csv](data/09_real_near_wall.csv)。
