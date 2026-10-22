# 活动代码与结果身份

`src/mevo_mb/fast_sim/accelerated_contact.py`：FP64 全对索引、几何、阻力、PairMatrix、CPU/CuPy 因子、惰性工作集。

`accelerated_sequential.py`：原 Sequential QP 的逐字变换副本，仅替换矩阵/线性代数表示；测试逐字核查其全部 Gate 和重线性化循环。

`accelerated_water.py`：原字段与双壁面检查，显式 opt-in。`accelerated_fields.py`：可选同流打包读回，物理 kernel 不变。

`scripts/straight_tube_exp1/f10aa/prepare.py` 冻结状态与源文件；`profile.py` 原实现 cProfile/line_profiler；`benchmark.py` 固定输入基准；`trajectory.py` 复用事件与拒判逻辑、每步全量力输出；`full.py` 正式轨迹；`parallel.py` spawn 独立任务；`extension.py` 条件延长；`profile_extra.py` 传输与隔离四线程测试；`analyze.py` 全轨迹独立 KKT；`figures.py` / `videos.py` 从保存数据绘制；`tests/f10aa` 原门槛、线性相关约束及非法输入失败闭锁。

输出均在 `results/straight_tube_exp1/f10aa`。不存在新默认生产后端切换；旧 F10S/F10Z 入口继续可用。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
