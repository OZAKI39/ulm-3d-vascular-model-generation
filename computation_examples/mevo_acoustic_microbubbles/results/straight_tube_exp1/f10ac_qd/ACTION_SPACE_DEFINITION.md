# 动作空间与停止条件

搜索初态为 F10AB 的 `S2_r0`，64 个半径 1.78 µm 的固定 ID 微泡，几何合法、初始邻接连通；它是人工指定的演示初态，不是本轮声学自发聚集得到的初态。三个档案候选的新增验证种子为 301801/301802/301803；档案外候选另用 401801/401802/401803。搜索及各自候选选择均未使用对应验证初态。N32 用于历史回归和新频率短轨迹准入。背景速度严格为零。

G2: [common X amplitude 200..1100 Pa-eq, X+ relative phase 150..240 deg, common Y amplitude 0..1100 Pa-eq, Y+ phase 0..360 deg]; X- and Y- phases zero。固定相位规约已通过共同相位平移的力不变性测试；最大力差约 1.0339757656912846e-24 N。

G1 完成 OFF/X/Y/四源短测试；G2 为正式低维搜索。G3 七变量完整四路和 G4 分段动作均 NOT_RUN：本轮优先完成数值与频率准入，没有建立额外收益证据。完整一条轨迹不变频。幅值是 Pa-equivalent，不是电压。预算和时限见 PREREGISTRATION.json。

`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `training_ready=false`
