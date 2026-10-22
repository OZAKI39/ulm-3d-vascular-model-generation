# 有效代码、执行与复现入口

源代码在项目 `scripts/straight_tube_exp1/f10ab/`：prepare.py冻结初态/输入；metrics.py独立形状定义；run.py继承原逐泡积分；plan.py有限筛选；analyze.py独立KKT复核；evaluate.py种子/数值准入；sensitivity.py中心差分；target.py预登记目标；figures.py/videos.py实际数据可视化；package.py自包含数据；report.py本报告。

配置：`configs/straight_tube_exp1/f10ab/protocol.json`；测试：`tests/f10ab`与`tests/f10aa`。

执行快照由logs/*EXECUTION.json逐次给出commit、目录、Supervisor标签、超时，logs/*supervisor.conf保存精确命令；复用快照只新增外部wrapper，不修改既有源码。所有计划在plans/，每个阶段START/COMPLETE保存真实身份与耗时。物理环境 `.venv`，绘图 `.postvenv`；均设置PYTHONPATH=src:.和BLAS线程1。

复现到新的独立输出根，勿覆盖已存在的stage目录。示例：`python -m scripts.straight_tube_exp1.f10ab.run pilot`；分析：`python -m scripts.straight_tube_exp1.f10ab.analyze pilot`；绘图：`python -m scripts.straight_tube_exp1.f10ab.post media`。先恢复PREREGISTRATION列出的原响应表和缓存，SHA必须一致。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
