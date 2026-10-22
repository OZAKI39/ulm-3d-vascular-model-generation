# 官方实现与审查依据

- [SciPy cho_solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html)：因子复用及多 RHS。
- [CuPy Cholesky](https://docs.cupy.dev/en/stable/reference/generated/cupy.linalg.cholesky.html)：GPU SPD 求解；实际版本13.6 API以本机安装源为准。
- [threadpoolctl](https://github.com/joblib/threadpoolctl)：隔离线程数量。
- [NVIDIA Warp](https://github.com/NVIDIA/warp)：本机1.18.0 ScopedTimer CUDA events及数组所有权/同步接口。
- [line_profiler](https://github.com/pyutils/line_profiler)、[Python cProfile](https://docs.python.org/3/library/profile.html)：原实现瓶颈实测。

性能及科学结论来自本轮保存数据，上述文档不作为本轮加速数字的外部证据。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
