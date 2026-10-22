# 传输与同步审计

原 Warp 每步 H2D 3 次、72N bytes；D2H 8 次、112N+12 bytes。包括声场、F2 读回和两个各自独立的原壁面检查。Packed 只把压力/一次力/总力/错误位打包：D2H 5 次、120N+8 bytes，少 3 次同步读回，但略增字节数；不减少力刷新和墙面验收。

所有 Warp kernel、设备内主力复制和 pack 同一当前 stream，numpy() 完成同步。primary_device 保留 F2 累加前的主力，原压力与力计算完全未改。64 泡完整计算有小幅附加收益，128 泡实测差异很小；不声称全规模显著传输加速。

CuPy 额外上传/下载（包含 A 和求解列）在 profile_extra/TRANSFER_ALLOCATION.json 精确按调用记录。benchmark FULL_STEP 的 transfer_bytes 仅是 Warp，不能把它当成 CuPy 总流量。隐式 CUDA API 的每个内部同步没有 Nsight 全轨迹追踪，报告为未单独拆分；已测 wrapper 的等待+复制和 CUDA kernel 时间。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
