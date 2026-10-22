# 冻结动作与公平性

| 动作 | 模式 | 幅值 Y-/Y+/X-/X+ / Pa-eq | 相位/deg |
| --- | --- | --- | --- |
| A0_OFF | ONE_HOT_PZT | 0, 0, 0, 0 | 0, 0, 0, 0 |
| A1_XM | ONE_HOT_PZT | 0, 0, 1036.2, 0 | 0, 0, 0, 0 |
| A2_XP | ONE_HOT_PZT | 0, 0, 0, 1036.2 | 0, 0, 0, 0 |
| A3_X0 | AXIS_PAIR_CONTROL | 0, 0, 1036.2, 1036.2 | 0, 0, 0, 0 |
| A4_X90 | AXIS_PAIR_CONTROL | 0, 0, 1036.2, 1036.2 | 0, 0, 0, 90 |
| A5_X180 | AXIS_PAIR_CONTROL | 0, 0, 1036.2, 1036.2 | 0, 0, 0, 180 |
| A6_X270 | AXIS_PAIR_CONTROL | 0, 0, 1036.2, 1036.2 | 0, 0, 0, 270 |
| A7_X210 | AXIS_PAIR_CONTROL | 0, 0, 1036.2, 1036.2 | 0, 0, 0, 210 |
| A8_YM | ONE_HOT_PZT | 1036.2, 0, 0, 0 | 0, 0, 0, 0 |
| A9_YP | ONE_HOT_PZT | 0, 1036.2, 0, 0 | 0, 0, 0, 0 |
| A10_Y90 | AXIS_PAIR_CONTROL | 1036.2, 1036.2, 0, 0 | 0, 90, 0, 0 |
| A11_FOUR | FOUR_CHANNEL_INDEPENDENT | 1036.2, 1036.2, 1036.2, 1036.2 | 0, 90, 0, 210 |
| B3_X0_L2 | AXIS_PAIR_CONTROL | 0, 0, 732.71, 732.71 | 0, 0, 0, 0 |
| B4_X90_L2 | AXIS_PAIR_CONTROL | 0, 0, 732.71, 732.71 | 0, 0, 0, 90 |
| B7_X210_L2 | AXIS_PAIR_CONTROL | 0, 0, 732.71, 732.71 | 0, 0, 0, 210 |

频率2.5MHz，背景流严格为零，最大允许1500 Pa-equivalent。基础X每片1036.201658Pa-eq。B系列每片缩小√2，使Σ|u|²与一个基础单源相等；这只是源幅值预算代理，不是校准声功率。A5/A7幅值预算相同，只改变相位。四片复幅与局部压力分别保存，不将Pa-eq称为泡处压力或AFG电压。

15动作先做105组静态诊断；完整积分有限筛选，没有全网格扫参。Y与四源仅静态探索，动态结论 `NOT_RUN`。所有阶段计划在执行前冻结于 plans/。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
