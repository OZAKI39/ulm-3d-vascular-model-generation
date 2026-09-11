云端部署、小测试与回传已完成：CLOUD_BUILD_PASS、CLOUD_SMOKE_PASS、RESULTS_RETURN_VERIFIED。

记录时间（UTC）：2026-09-10T22:25:46.385378+00:00。BUILD_ID：`mirheo-sm120-20731713865ae510`；JOB_ID：`cloud-smoke-20260910T221827Z`。

实际来源为上传清单指定的修复树：`/workspace/bloodflow/uploads/20260910T212439Z/source/mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source`。对应原提交 `8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf`；源码清单聚合 SHA-256 为 `2b30d14be357b1bcd7a63805cc164f4907d291b38fd8ae4d4db44980b3f0ac90`。原本地路径为 `/home/lzy/projects/mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source`。

构建使用独立可写副本 `/workspace/bloodflow/work/mirheo-sm120-20731713865ae510/source`。仅适配归档版本信息、预编译库安装路径和旧 CMake 的 sm_120 架构分支；原修复补丁没有重复应用，碰撞、力学和积分源码保持不变。原 vendor 未被选为运行来源。子模块实际文件齐全；原归档注明的两个 Windows NuGet 包装文件差异不影响 Linux 构建。

Python 模块实际加载自 `/workspace/bloodflow/.venv/lib/python3.12/site-packages/mirheo/__init__.py`。计算库实际加载自 `/workspace/bloodflow/.venv/lib/python3.12/site-packages/libmirheo.cpython-312-x86_64-linux-gnu.so`，SHA-256：`d45b4fd1b4498b6f365a012e758eec35bc58a7844dabc66f22e617a7b8cf84ee`。这是云端新二进制身份，不替换历史 WSL 库哈希。

CUDA 12.8.93、GCC 13.3.0、CMake 3.28.3、OpenMPI 4.1.6、并行 HDF5 1.10.10，RTX 5090。96 个 CUDA 编译命令与 cuobjdump 的 cubin 清单均只指向 sm_120；ldd 无 not found。保留 Release、32 位主精度、膜/杆双精度 OFF、fast-math，编译并行 2。

项目 .venv 新增 NumPy 1.26.4 和本轮 Mirheo 1.6.2；pip、setuptools、wheel、packaging 未升级。完整前后版本与安装报告在构建证据目录。/venv/main 的文件身份指纹、受保护工具二进制哈希均保持一致。

实际小测试（两 rank：1 计算 + 1 后处理；单 GPU、逐项执行）：

| 项目 | 实际执行 | 退出码 | 进程墙钟时间 |
|---|---|---:|---:|
| A | 真实库加载、MPI 归约、CUDA 同步 | 0 | 1.001 s |
| B | 512 个 DPD 粒子，100 步，dt*=0.0005，实际 t*=0.0500000024；数值有限、粒子移动、ID 唯一 | 0 | 0.801 s |
| C | 1 个可变形 RBC、642 顶点、110592 个液体粒子；250 步连续零剪切准备，另有 4000 步冻结壁准备；保留内外液体和 bounce-back | 0 | 2.103 s |

没有 Γ=4 剪切、checkpoint/restart 或物理参数搜索。C 保留本地 A6 参数，仅缩短准备时长；只记录膜的原生 HDF5 帧、最终流体统计及小探针，未输出逐步全液体轨迹。正常退出及有限值不能证明材料匹配、长时稳定性或严格不可渗透性。原物理验收仍为 NOT_MATCHED，原共同终点比较未完成的结论保持不变。

结果目录：本地 `/home/lzy/projects/cloud_results/cloud-smoke-20260910T221827Z`；云端 `/workspace/bloodflow/cloud_runs/cloud-smoke-20260910T221827Z`。55 个清单内文件，共 635494 字节；另有清单和 RESULTS_READY 控制文件。自动 rsync 到 .partial 后核对路径边界、数量、大小、SHA-256、执行位并解析 JSON/CSV，再原子归档。CSV 实际有 512 行液体粒子、642 行膜顶点、24 行流场剖面及 1 行最终矩统计。已完成小测试结果未回传量为 0。

构建证据（包括失败）也已归档：`/home/lzy/projects/cloud_results/mirheo-sm120-20731713865ae510`、`/home/lzy/projects/cloud_results/mirheo-sm120-4aca28ba929befae`、`/home/lzy/projects/cloud_results/mirheo-sm120-4aca28ba929befae-a2`。四个额外 MPI 诊断文件经单独白名单下载并逐文件哈希校验，见 `records/mpi_diagnostics_return.json`。编译库/中间对象保留云端，不作为小测试结果重复下载。

首次 MPI hostname 启动无输出并超时；详细日志显示 OpenCL 已禁用而 GL 拓扑组件仍启用。仅在本任务子进程使用 HWLOC_COMPONENTS=-opencl,-gl 后，真实 Allreduce/Bcast/Barrier 通过。随后发现 CMake 3.28 的旧 FindCUDA 只接受一位主架构号；新增受 CUDA >=12.8 检查保护的 sm_120 分支后构建成功。失败耗时和日志保留，没有自动循环重试。参考 [hwloc 环境变量说明](https://www.open-mpi.org/projects/hwloc/doc/v2.14.0/envvar.html)。

原生编译 256.031 s；所有构建入口累计 304.464 s，低于 1800 s。累计小测试/启动诊断进程时间 40.426 s（含失败），低于 180 s；实际 A/B/C 合计约 3.905 s。记录可核查工作窗口 2026-09-10T21:59:57.250650+00:00 至 2026-09-10T22:18:36.507156+00:00，共 1119.257 s，包含检查、排错、编排和等待；不是纯编译时间。租赁开始时间和账单未读取，租赁时长及费用记为未核实，不将上述秒数当作 Vast 费用。

/workspace 已用空间较首次核查增加 1606897664 字节（1.497 GiB），最终可用 18811752448 字节（17.520 GiB）。目录实际分配空间及最大的 20 个新增文件详见 `records/final_environment.json`：

| 目录 | 实际分配字节 |
|---|---:|
| `/workspace/bloodflow/uploads/20260910T212439Z` | 499834880 |
| `/workspace/bloodflow/build` | 1031221248 |
| `/workspace/bloodflow/work` | 282742784 |
| `/workspace/bloodflow/.venv` | 316280832 |
| `/workspace/bloodflow/cloud_runs` | 774144 |

最终核查确认：本地冻结副本、其对应的 8694 个本地正式文件、云端 uploads 的内容和执行位保持一致；本任务构建/测试 PID 均已退出，GPU 无计算进程。

清理预览确认本地归档和远端哈希一致，但 4 个现有 Syncthing 相关进程（1371、1394、1505、1518）的文件句柄不可读，不能完整排除占用。因此可直接放行的清理候选为 0；55 个已校验文件列为等待进程占用核查，详见 `records/cloud-smoke-20260910T221827Z-cleanup-preview.json`。没有删除云端文件，也没有停止这些既有服务。

工具测试 13 项通过；还实际核查了重复 deploy 的二进制哈希复用，以及重复 fetch、verify、cleanup --dry-run。失败任务的封口日志已真实回传并验证。

后续命令（本地 WSL）：

```bash
cd /home/lzy/projects/cloud_compute
python3 -B cloud_run.py status
python3 -B cloud_run.py smoke --new
python3 -B cloud_run.py status cloud-smoke-20260910T221827Z
python3 -B cloud_run.py fetch cloud-smoke-20260910T221827Z
python3 -B cloud_run.py verify cloud-smoke-20260910T221827Z
python3 -B cloud_run.py cleanup cloud-smoke-20260910T221827Z --dry-run
```

`smoke --new` 输出新的 JOB_ID，并在该任务结束后自动 fetch；手动命令应使用它实际输出的新编号。本轮检查账本继续受 180 秒累计限制，不会借用旧 RBC 预算。完整说明见 README_zh.md。

本轮完成的是云端部署和文件回传验证，没有自动运行完整红细胞长实验，没有改变原有科学结论，没有自动清理未获确认的数据或销毁服务器。
