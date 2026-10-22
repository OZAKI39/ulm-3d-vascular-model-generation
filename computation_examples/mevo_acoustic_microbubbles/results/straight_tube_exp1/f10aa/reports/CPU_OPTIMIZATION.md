# CPU 数学等价优化

粒子对保持 np.triu_indices 的原始顺序；缓存仅不可变索引。距离与法向向量化，润滑对仍采用原 activation=0.1、floor=0.001 和原 φ 表达式。对角块用 add.at 累加，非对角块唯一赋值；不创建原路径未使用的 coupling/groups。

浮点相加顺序可能改变最后几位；A 最大差约 2.3e-13（无量纲），由原 Gate、独立原 A/B 力平衡和全轨迹比较验收。原 A 每步重建；一次步内固定 A 复用因子，绝不跨不同位置复用旧 Cholesky。

单线程为推荐稳健配置。双线程全部冻结状态通过，部分中等 Cholesky 小幅改善，但小矩阵没有稳定整体优势。运行时从1线程切换到4线程曾发生native崩溃；进程启动前固定OPENBLAS_NUM_THREADS=4的独立测试通过。进一步2/4线程完整轨迹和全部冻结状态见threaded；失败切换路线不作性能样本。未额外引入 Numba/MKL，因为已测结构性收益足够，未无界扩工具栈。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
