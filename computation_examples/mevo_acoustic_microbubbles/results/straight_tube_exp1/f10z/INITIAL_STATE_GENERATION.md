# 初态生成与合法性

32 个固定 ID，R=1.78 µm，2.5 MHz；背景流严格为零；继承平移阻力黏度 0.001 Pa·s、径向响应有效黏度 0.002 Pa·s、密度 1000 kg/m³、声速 1500 m/s。这是原模型参数，未宣称为本实验纯水校准。

| 初态 | 种子 | split | 最小距 µm | Rg µm | 分量 | 近接触数 | 合法 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| original | None | historical_baseline | 3.560178 | 5.903887 | 1 | 64 | True |
| low_train | 18001 | train | 12.562421 | 40.852796 | 32 | 0 | True |
| medium_train | 18002 | train | 5.043357 | 18.491836 | 32 | 0 | True |
| medium_validation | 28002 | validation | 5.106813 | 19.808404 | 32 | 0 | True |
| medium_test | 38002 | test | 5.387992 | 18.665312 | 32 | 0 | True |
| local_dense | [18003, 18004] | train | 3.858324 | 20.967787 | 32 | 0 | True |
| multi_cluster | [18010, 18011, 18012, 18013] | train | 4.265083 | 22.353144 | 32 | 0 | True |
| near_contact | 18005 | train | 3.714992 | 6.226581 | 20 | 0 | True |
 
拒绝采样只检查几何，不根据未来结果筛选。low 使用 ±40 µm、最低间距 12 µm；medium 使用 ±18 µm、最低间距 5 µm，去质心后检查缓存边界。local_dense 是 16+16 密度对比，multi_cluster 为四个预置局部空间簇，各自尚无 d≤2.1R 连边。near_contact 是历史布局放大 1.055 倍并小幅扰动，初始已有 12 条几何邻边，必须单列。original 是人为给定的连通对照，不能当作声学形成的团。
 
无硬球重叠且墙间隙为正。CONTACT_FREE_INITIALIZATION_PASS 仅指新增分散初态无实际硬接触；original 本来就在 1 nm 近接触统计阈值内。几何近邻 d≤2.1R，不表示微泡融合。固定输入和种子见 PREREGISTRATION.json、INITIAL_STATES.h5、INITIAL_STATE_AUDIT.json。
