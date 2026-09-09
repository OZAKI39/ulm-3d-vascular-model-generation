# Mirheo / HemoCell 纯流体对照

固定配置：`py_scripts/solver_benchmark.yaml`。所有默认/查看/帮助路径不启动求解器。

最简单的查看命令：

```bash
cd /home/lzy/projects/mirheo_starter && .venv/bin/python -B -m test_code.review_solver_benchmark --open
```

页面是内嵌 Plotly JavaScript 的单文件 HTML，可从 Windows Edge/Chrome 直接打开，
无需 Web 服务。最近一次页面路径在 `data/solver_benchmark/LATEST.json`。

```bash
# 固定版本部署；已有安装核对后复用，绝不删除 Palabos 或重编译 Mirheo。
bash /home/lzy/projects/hemocell_starter/scripts/setup_hemocell.sh

cd /home/lzy/projects/mirheo_starter
# 单位、冻结任务、当前预算和已有证据；无新求解。
.venv/bin/python -B -m py_scripts.benchmark_mirheo_hemocell --config py_scripts/solver_benchmark.yaml --preflight-only
# 执行尚未运行且已合法授权的冻结任务。已有成功任务校验哈希后复用；失败不自动重试。
.venv/bin/python -B -m py_scripts.benchmark_mirheo_hemocell --config py_scripts/solver_benchmark.yaml --execute
# 原始数据回读、共同分析与页面，不编译、不重跑。
.venv/bin/python -B -m py_scripts.benchmark_mirheo_hemocell --config py_scripts/solver_benchmark.yaml --analyze-only
# CPU 测试，不启动 GPU。
.venv/bin/python -B -m unittest test_code.test_solver_benchmark -v
```

运行目录为 `runs/solver_benchmark/solver_benchmark_20260909/{cpu,gpu}/`。
每个任务独立目录含精确参数、命令、XML（LBM）、控制器墙钟、主机内存采样记录、
真实剖面 CSV、分块计时及完成记录。共享原单位、冻结模型与预算计算模块；未更改其行为。
原始文件通过 `output_sha256.json` 封存，分析前校验。新分析用新指纹目录，不覆盖旧报告。

HemoCell 是原生 C++ 可执行文件：
`/home/lzy/projects/hemocell_starter/build/benchmark/pure_fluid_benchmark`。
调度器自动生成独立 XML，并在正确工作目录启动，不需要手工复制配置或 RBC/PLT 文件。
本例根本不创建 cellfield、不注册细胞类型、不加载位置文件，不经过无效细胞回调。

共同定义：4 µm 周期立方盒，ρ=1056 kg/m³，ν=3.27e-6 m²/s，μ=ρν；
下半盒 −x、上半盒 +x，a=50000 m/s²，F粒子=m a。驱动来自小盒设计，不是旧血管泵流量。
0–1 µs 无驱动准备，1–3 µs 建立流场，3–8 µs 统计。每次从新状态开始，
SDPD 保留随机初速度、一次 COM 扣除与原生随机/耗散作用，不读 checkpoint。
在同一个 coordinator 中通过原生 deregisterIntegrator/setIntegrator 切换体力，
粒子状态、时间和随机发生器不重置。

LBM：D3Q19 Guo BGK，τ=1，dx=0.25/0.125 µm。dt 从同一 ν 计算，
不是匹配 SDPD 数值 dt。半格偏移节点构成 N 个周期唯一点，避免 N+1 或变号平面偏置。
使用已含半步外力修正的 computeVelocity，同时输出独立原始动量以核查。
只有 unique bulk 参与 MPI 求和，halo 不重复计数。

SDPD：读取原冻结映射完整精度，L0=5e-7 m、m*=1、n*=8、kBT*=1、
Wendland C2 rc*=1、Linear EOS cs*=120、rho0*=0、dt*=1e-6。
μ* 读取原单位结果，不手工舍入。当前 libmirheo 为单精度；不重编译、不动 .venv。
本次不输出应力/压力，stress=False，保留密度诊断和粒子检查。此输出设置在结果中公开，
不把历史带压力输出的速度伪装成新工作流成本。

LBM 背景温度是物性对应的 298.15 K，不虚构热噪声测量；SDPD 温度单独验收。
CI 使用原修正的完整块均值函数，并在独立块不足或非稳态时输出 null；
所有空间 bin 保留。LBM 统计误差不适用、离散误差实测、物理模型误差未量化分别记录。

CPU 数值测试累计上限 1200 s；单任务≤600 s；一项求解器运行。
编译也使用同一排他锁。MPI 绑定真实可用核心、OMP_NUM_THREADS=1，不 oversubscribe。
对所有 rank 取最大完成时间；GPU 在计时结束前 cudaDeviceSynchronize。
主机 RSS 与 GPU 全设备显存采样分开。编译与下载永不加入 ms/步。

内存采样发现并保留了一项缺测：初版按进程组采样，只捕获 mpirun，未覆盖另建进程组的
OpenMPI ranks。原始记录没有改写；分析将其标记为 `LAUNCHER_ONLY_RANK_RSS_NOT_MEASURED`，
求解器进程树 RSS 为 null。修正为递归追踪子进程后，另外实跑了 4 项 CPU 内存审计，
分别覆盖 16³ 的 1/2/4 rank 和 32³ 的 1 rank，计入同一 CPU 预算，与原计时重复分列。
这四项的进程树 RSS 采样峰值依次约为 37.25/61.40/109.49/43.41 MiB。
Mirheo 主机 rank RSS 仍为未测；GPU 全设备显存采样峰值为 1725 MiB，含桌面及其他程序，
不能称为本任务独占显存。没有为补这项缺测再次运行 GPU。

用户新增 600 GPU 秒和一次显式重试批准均以原文、哈希、范围记录。
第一次启动因缺少 mpi4py 在导入 Mirheo 之前失败，耗时 2.458043783 s，实际流体步数为零。
已替换为现有 OpenMPI 的 C ABI，未安装依赖；双 rank CPU 通信测试是真实传输，
但不作为流体测量。显式重试最多用余下 597.541956217 s。两次尝试都保留并计费。
不得再用同一批准自动重试第三次或更改 campaign 重置用量。

网格核查与局部平均流场通过，不等于真实压力出口或完整物性验收。
圆管标记 `PIPE_BENCHMARK_DEFERRED`：SDPD 近壁核截断、冻结壁面与滑移仍缺验证。
不做 RBC；后续还须匹配细胞尺寸、膜模型、浓度、内外黏度及耦合精度。

来源和许可：HemoCell/Palabos 与派生 C++ 案例 AGPL-3.0-or-later，Mirheo MIT。
原生来源头部保留；安装与补丁证据在 `../hemocell_starter/metadata/`。
两个旧项目只读；Mirheo 源码、二进制、兼容补丁和 Python 环境哈希在保护清单中前后核查。
人工验收始终 `PENDING`；自动 Chrome 检查另存，不能代替人工验收。

本轮交付：19 次 CPU 数值运行（1 冒烟、8 主基准、6 低 I/O、4 内存审计），
CPU 累计 119.365107860 s。SDPD 一次启动失败加一次成功重试，共计费 320.162673519 GPU s，
新增的 600 s 额度还余 279.837326481 s；此余额不授权继续扫描或第三次重试。
主基准全部推进约 8 µs。16³ / 32³ LBM 剖面误差分别为 1.436% / 0.359%，
粗细网格差为 1.077%，通过本轮 `PROPOSED` 局部筛选。
SDPD 完整块描述均值的剖面误差为 20.172%（全窗口描述均值为 17.848%），
仅有 1 个完整剖面块，温度约 325.16 K 且仍有漂移，无有效 CI，因此合格解耗时与速度比为 null。

最终页面目录中的 `结果说明.md`、`delivery_record.json` 和 `solver_benchmark_results.zip`
提供摘要与结果包；ZIP 内含单文件 HTML、本轮原始数据、新增代码、部署证据和哈希清单。
结果包不含编译二进制、虚拟环境或全部旧项目；重新运行依赖本机已有正式模块及固定安装。
解压后可直接用 Windows 浏览器打开根目录 `benchmark_review.html`。
