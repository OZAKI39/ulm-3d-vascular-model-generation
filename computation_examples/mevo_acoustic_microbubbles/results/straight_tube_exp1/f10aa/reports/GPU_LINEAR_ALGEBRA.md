# 实际 CuPy FP64 负结果

| N | 物理时间 | 原CPU完整墙钟/s | 优化CPU/s | Packed/s | CuPy/s | CPU加速 |
|---|---:|---:|---:|---:|---:|---:|
|32|2 ms|7.7872|2.6455|2.4808|4.6503|2.944×|
|64|2 ms|34.6319|5.7603|5.4442|7.7861|6.012×|
|128|1 ms|44.3742|6.5612|6.5523|7.8433|6.763×|

CuPy 13.6 实际执行 GPU pair distances、FP64 向量化阻力、Cholesky、三角多 RHS 求解。冻结状态和完整轨迹均通过原数值门槛；当前规模上比优化 CPU 慢。GPU 路线为明确的 CPU 接触/GPU 因子混合实验，包含上传 A、每次 RHS 及解下载，不能称全流程常驻 GPU。

Warp/CuPy 共享指针与值一致测试通过，借用前同步、保持原所有者引用；未将这项只读互操作测试误称部署了全流程零拷贝。源码继续使用安全的原流路径。未使用 FP32/TF32；PyTorch 和 OSQP CUDA 为 NOT_RUN，未因没有性能收益强行更换另一种近似 QP。

初版 cp.errstate API 不存在，改用已安装版 cupyx.errstate；修复前失败输入和日志保存，修复后完整通过。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
