# 本次未上传内容

本地原文件和服务器文件均不删除。逐项记录见 `excluded_local.json` 及 `server_inventory.json` 的 omitted 字段。

- `.git`、虚拟环境、Python/测试/渲染缓存、已安装 SimVascular 分发包：可重新安装或本就不是科学源码与结果。
- `input_data/streamlines/data/candidates.npz`：可再生成的候选种子中间缓存；保留已接受流线、实际输入场、计算脚本和验证结果。
- 用户已取消细档真实血管算例的 `SV_MESH/mesh-complete.mesh.vtu`（约92.5MiB）、`SV_MESH/mesh_arrays.npz`（约61.4MiB）、`raw_mesh/volume.vtu`（约62.4MiB）：不属于已完成 CFD 结果。配置、取消状态、日志和其余证据保留，完整网格继续存在于本地原目录。
- 第三方求解器源码中的950个 Git LFS 占位文件：本地只有指针，没有上游大文件内容；完整指针保存在 `upstream_lfs_pointers.json`。不将占位指针伪装成已下载的测试数据，也避免触发错误的 LFS 上传。
- 服务器编译产物、环境、socket、PID 和重复归档按 `tools/collect_server.py` 排除；正式二进制路径、大小和哈希保留在服务器清单。
- 账户密钥和环境认证文件不在科学快照范围内。

通常不收录单个≥50MiB的非关键文件，但两批微泡的 `gpu_mesh_input.npz` 是正式核验依赖，约55.6MiB，**明确保留**。没有因体积较大而删减最新1500条原始轨迹、逐轨迹审计日志、0.5ms导出数据或正式动画。每个 Git 文件均小于100MiB。

既有 `sync_metadata/current_20260927/` 记录上一次同步排除范围，日期及语义不改写。共享几何/第三方代码范围沿用该次已整理集合；本轮读取当前同路径文件刷新，未纳入与该血管工作流无关的大型影像或打印结果。
