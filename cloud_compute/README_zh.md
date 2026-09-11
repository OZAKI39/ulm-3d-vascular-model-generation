# Vast.ai Mirheo 部署与结果回传

本地 WSL 工具目录：`/home/lzy/projects/cloud_compute`。通过现有 `vast-mirheo` SSH 别名连接。主工程、上传冻结副本、旧预算和科学结论不由这些工具修改。

```bash
cd /home/lzy/projects/cloud_compute
python3 -B cloud_run.py preflight
python3 -B cloud_run.py deploy
python3 -B cloud_run.py status
python3 -B cloud_run.py smoke --new
```

`deploy` 启动专属 tmux 构建并立即返回；`status` 查看编译状态。成功后再次 `deploy` 校验已安装二进制哈希并复用，避免重复编译。

`smoke --new` 必须显式指定新任务，每次产生新的 UTC JOB_ID。它依次运行加载检查、100 步 DPD 液体和可变形 RBC 的 250 步短准备；RBC 另包含两面冻结壁各 2,000 步准备。材料、质量、碰撞、时间步和液体参数来自上传的 A6 配置。`--skip-repair` 可只做 A/B。不会执行 Γ=4 或参数搜索。

任务启动后，远程 tmux 负责限时运行；本地轮询任务完成后自动下载、校验并归档。如果本地断线，服务器保持独立运行，结果保留云端。在本地恢复时使用同一 JOB_ID：

```bash
python3 -B cloud_run.py status JOB_ID
python3 -B cloud_run.py fetch JOB_ID
python3 -B cloud_run.py verify JOB_ID
python3 -B cloud_run.py cleanup JOB_ID --dry-run
```

将 JOB_ID 替换为启动输出或 `cloud_runtime.json` 中的值。`fetch` 只下载已封口的结果，不启动求解器；已有归档会核查清单身份。中断时保留 `cloud_results/JOB_ID.partial`，再次 fetch 同一编号继续。失败任务同样回传关闭的日志和实际输出。

构建日志也能回传；BUILD_ID 取自 `build_plan.json`：

```bash
python3 -B cloud_run.py fetch BUILD_ID --build
python3 -B cloud_run.py verify BUILD_ID --build
```

显式审阅失败日志后才使用 `deploy --resume`，它保留失败回执并复用同一 work/build，受次数和累计构建时间限制。后续尝试的证据编号为 `BUILD_ID-a2` 等。调整构建兼容性后使用 `preflight --revise-build-plan` 生成绑定新补丁的构建身份；旧记录保留。

限制：编译并行 2，构建保护 1,800 秒；小测试单项最多 60 秒，本轮累计最多 180 秒，并发 1。MPI 启动诊断和失败也计入本轮检查耗时。这些是操作保护，不是 Vast 实际计费，也不使用旧 RBC 余额。低于 3 GiB 空闲不启动或继续任务；构建前还需预计增长后保留至少 5 GiB。

本实例在 `HWLOC_COMPONENTS=-opencl` 时仍启用了 GL 拓扑插件并挂住。经过有超时的诊断后，本工具仅在子进程设置 `HWLOC_COMPONENTS=-opencl,-gl`；两 rank 的实际 Allreduce/Bcast/Barrier 随后通过。没有更改全局环境或 MPI 安装。参考 [hwloc 环境变量文档](https://www.open-mpi.org/projects/hwloc/doc/v2.14.0/envvar.html)。

构建适配仅在 work 中：归档版本信息、预编译库安装路径、CMake 旧架构解析器的 sm_120 分支。CUDA 最低版本检查保留，未改物理公式，未替换第三方库。包仅安装到 `/workspace/bloodflow/.venv`；`/venv/main`、驱动和 Toolkit 保持不变。

`cleanup --dry-run` 仅列出本地已验证、远端哈希仍相同且没有进程持有的结果文件。工具没有真实删除入口，不执行 rsync 删除选项，不停止或销毁实例。

```bash
python3 -B -m unittest discover -s tests -v
```

测试覆盖映射、误选旧 vendor、上传副本变化检测、路径越界、链接/特殊文件、未完成结果、内容和执行位损坏、重复 fetch、失败日志、清理限制、空壳输出及真实子进程超时/磁盘停止。具体执行回执在 `records/tool_tests.json`。
