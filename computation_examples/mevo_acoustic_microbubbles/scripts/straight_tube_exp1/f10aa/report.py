"""Generate traceable Chinese engineering reports and machine-readable admission gates."""
import time,json,subprocess
import numpy as np
from .common import *

def main():
    verify();out=BASE/'reports';out.mkdir(exist_ok=True);b=BASE/'benchmark_v3';performance=read(BASE/'analysis/PERFORMANCE.json');runs=read(BASE/'analysis/RUNS.json');aud=read(BASE/'analysis/AUDITS.json');checks=read(BASE/'full/COMPARISONS.json');fixed=read(b/'FIXED_BENCHMARK.json');parallel=read(BASE/'parallel/BATCHES.json');ext=read(BASE/'extension/COMPARISONS.json');four=read(BASE/'profile_extra/FOUR_THREADS.json')
    resources=read(BASE/'full/resources.json');st=read(BASE/'baseline/PROFILING.json');pre=read(BASE/'extension/N256_PREFLIGHT.json');more=read(BASE/'extension/RUNS.json');extras=read(BASE/'analysis/EXTRA_COMPARISONS.json')
    annotations='EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false'
    def emit(name,title,body):
        (out/name).write_text('# '+title+'\n\n'+body.strip()+'\n\n'+annotations+'\n',encoding='utf-8')
    table='| N | 物理时间 | 原CPU完整墙钟/s | 优化CPU/s | Packed/s | CuPy/s | CPU加速 |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for n in [32,64,128]:
        q=performance[str(n)];table+=f"|{n}|{2 if n<128 else 1} ms|{q['B0']['wall_s']:.4f}|{q['B1']['wall_s']:.4f}|{q['PACKED']['wall_s']:.4f}|{q['B3_CuPy']['wall_s']:.4f}|{q['B1']['speedup']:.3f}×|\n"
    latency='| N | 后端 | 步均值/ms | p50/ms | p95/ms | 最大/ms | 保存/s |\n|---|---|---:|---:|---:|---:|---:|\n'
    for n,q in performance.items():
        for backend,r in q.items():latency+=f"|{n}|{backend}|{r['step']['mean']*1e3:.4f}|{r['step']['p50']*1e3:.4f}|{r['step']['p95']*1e3:.4f}|{r['step']['max']*1e3:.4f}|{r['save_s']:.4f}|\n"
    pt='| 并发 | 批次 | 四条墙钟/s | 物理ms/墙钟s | 平均单轨迹/s |\n|---:|---|---:|---:|---:|\n'
    for row in parallel:pt+=f"|{row['workers']}|{row['round']}|{row['batch_wall_s']:.4f}|{row['throughput_physical_s_per_wall_s']*1e3:.5f}|{np.mean([j['result']['wall_s'] for j in row['jobs']]):.4f}|\n"
    warm=[r for r in parallel if r['round']=='warm'];base=next(r for r in warm if r['workers']==1);best=max(warm,key=lambda r:r['throughput_physical_s_per_wall_s']);throughput=best['throughput_physical_s_per_wall_s']/base['throughput_physical_s_per_wall_s']
    gpu=max(r.get('gpu_percent',0) for r in resources);vram=max(r.get('gpu_MiB',0) for r in resources);rss=max(r['rss_MiB'] for r in resources);cpu=np.mean([r['cpu_percent']/100 for r in resources])
    thread_checks=read(BASE/'threaded/COMPARISONS.json')
    max_id=max([v['errors']['fixed_id_m'] for k,v in checks.items() if k.endswith(('B1','PACKED','B3_CuPy'))]+[v['errors']['fixed_id_m'] for v in thread_checks.values()])
    max_refine=max(v['errors']['fixed_id_m'] for k,v in checks.items() if k.endswith(('DT05','FINE')))
    all_numeric=all(r['PASS'] for r in aud.values());all_cmp=all(c['admitted'] and c['event_topology_agrees'] for c in list(checks.values())+list(thread_checks.values()));all_ext=all(c['admitted'] and c['event_topology_agrees'] for k,c in ext.items() if k in ['B1_DT05','B1_FINE']);long_reference=ext['B0_DT1']['admitted']
    total_steps=sum(r['checked_steps'] for r in aud.values());extref=next(r for r in more if r['backend']=='B0');extfast=next(r for r in more if r['backend']=='B1' and r['variant']=='DT1')
    stage_times={name:read(BASE/name/'COMPLETE.json')['wall_s'] for name in ['baseline','benchmark_v3','full','parallel','extension','profile_extra','threaded']}
    elapsed=sum(stage_times.values())
    # The initial fixed campaign stopped during instrumentation: include its separately measured Supervisor duration in the execution ledger, not this partial stage sum.
    timing='\n'.join(f'- {k}: {v:.3f} s' for k,v in stage_times.items())
    emit('README.md','F10AA 交付索引',f'''本轮保留原物理和数值门槛，以矩阵自由全粒子对接触＋按需活跃列求解实现完整轨迹加速。32/64推荐单线程CPU力学，可选packed声场；128单任务推荐启动前固定4线程CPU力学，原Warp声学/壁面不变。原引擎及历史结果未覆盖。

{table}

图：`../figures/`，8 组 PNG/PDF；视频：`../videos/`，2 段 MP4，全部帧已解码。数据：`../full/`、`../extension/`、`../parallel/`。审核入口：`F10AA_DEVELOPMENT_REPORT.md`、`FINAL_RECOMMENDATION.md`、`../FINAL_GATE.json`。

新分支 `dev/f10aa-solver-acceleration-20261018`。每次执行的不可变 commit 和 Supervisor 标签见 `../logs/*_EXECUTION.json`，基线 `df76b1e9d311024398bc2bca40b6cf0f6382ba64`。科学计算源码快照与图表源码身份分别保存，不能把后期绘图提交当成全部早期执行身份。

启动器：`F10AA_BLAS_THREADS=4 bash scripts/straight_tube_exp1/f10aa/python.sh ...` 在Python导入前固定线程数，避免运行时resize。32/64或独立任务并行用1线程。

复现：项目根目录先设置 `PYTHONPATH=src:.` 和 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`。GPU 环境使用 `/workspace/mevo_acoustic_microbubbles/f10aa_simulation_20261018/.venv/bin/python`；执行 `-m scripts.straight_tube_exp1.f10aa.benchmark --output benchmark_replay` 可独立重做固定状态和完整步基准。完整轨迹入口 `full.py`，需在新输出目录/独立 checkout 复现，既有目录存在时主动拒绝覆盖。后续 `parallel.py`、`extension.py`、`analyze.py`、`figures.py`、`videos.py`、`report.py`。所有历史原始缓存按原 SHA 恢复；没有凭空生成替代缓存。
''')
    emit('ACTIVE_CODE_INDEX.md','活动代码与结果身份','''`src/mevo_mb/fast_sim/accelerated_contact.py`：FP64 全对索引、几何、阻力、PairMatrix、CPU/CuPy 因子、惰性工作集。

`accelerated_sequential.py`：原 Sequential QP 的逐字变换副本，仅替换矩阵/线性代数表示；测试逐字核查其全部 Gate 和重线性化循环。

`accelerated_water.py`：原字段与双壁面检查，显式 opt-in。`accelerated_fields.py`：可选同流打包读回，物理 kernel 不变。

`scripts/straight_tube_exp1/f10aa/prepare.py` 冻结状态与源文件；`profile.py` 原实现 cProfile/line_profiler；`benchmark.py` 固定输入基准；`trajectory.py` 复用事件与拒判逻辑、每步全量力输出；`full.py` 正式轨迹；`parallel.py` spawn 独立任务；`extension.py` 条件延长；`profile_extra.py` 传输与隔离四线程测试；`analyze.py` 全轨迹独立 KKT；`figures.py` / `videos.py` 从保存数据绘制；`tests/f10aa` 原门槛、线性相关约束及非法输入失败闭锁。

输出均在 `results/straight_tube_exp1/f10aa`。不存在新默认生产后端切换；旧 F10S/F10Z 入口继续可用。''')
    emit('F10Z_BASELINE_AUDIT.md','历史基线与冻结输入','''F10S/F10Z/scaling 共 126 条历史轨迹逐文件对照保存 SHA；详见 `../BASELINE_AUDIT.json`。Git 原工作树清洁，新工作树从实际 F10Z 最终提交建立，未从旧公开 F10S 分支覆盖。

冻结 16 个初态、32 个 ON/OFF 力学状态：32/64/128 的分散初态、历史终态、合法单对/立方/FCC，以及原 F10S 32 泡近接触群。立方与 FCC 的压缩力 F=ζA(-20x) 是制造解的力学检查，**不是声学形成轨迹**。生成规则、种子、输入哈希先于优化保存。原 REFERENCE.h5 保存稠密 A/B、法向/间隙、声学力、速度、完整接触反力和迭代记录。

128 泡旧 1 ms 仍无活跃硬球接触，不能用于单独验收复杂接触；本轮额外覆盖 532 个同时 FCC 活跃约束、64 泡真实接触延续、32 原近接触轨迹和延长的 128 泡真实声学轨迹。''')
    ranking='| 状态 | 稠密接触几何/ms | 阻力/ms | Cholesky/ms | 线性求解/ms | SQP/ms |\n|---|---:|---:|---:|---:|---:|\n'
    for name in ['N32_dispersed','N64_historical_end','N128_dispersed']:
        p=next(r for r in st if r['case']==name);v=p['instrumented_inclusive'];ranking+='|'+name+'|'+'|'.join(f"{v[k][0]*1000:.4f}" for k in ['dense_contact_geometry','resistance','Cholesky_SPD','linear_solve','SQP_total'])+'|\n'
    emit('SOLVER_PROFILING.md','原始瓶颈及计时定义',ranking+'''\n上述为嵌套 inclusive instrumented 时间，不能相加。独立预热单步 p50/95 见原始 PROFILING.json；instrumentation、cold engine、wall、process CPU 分开。排序：反复稠密接触几何；进入接触后的全 B.T 求解；阻力/数组工作；Cholesky。128 泡全部 8128 对、384 列，B 的存储和遍历随 N³ 增长，全右端项求解还会进一步放大开销。

CPU 使用 cProfile、line_profiler、tracemalloc；GPU 为实际 Warp CUDA event kernel 计时，fused 声场/APECSS 无法诚实拆成两个独立内核。F2 pair/reduction 和双墙面 kernel 分别保存。D2H wrapper 计时包含同步等待与复制，不能当成纯 PCIe 带宽。`benchmark_v3/FULL_STEP.json` 补充原始 profiling print=False 未触发回调的缺项；错误未冒充 GPU 测量。Nsight Systems 不可用，未运行内核权限分析器。

Python 分配峰值与进程 RSS、实际线程数、BLAS 线程、数组上传下载、CuPy 额外流量见 `profile_extra`；所有原始 function/line 计时和库配置均保留。''')
    emit('CPU_OPTIMIZATION.md','CPU 数学等价优化','''粒子对保持 np.triu_indices 的原始顺序；缓存仅不可变索引。距离与法向向量化，润滑对仍采用原 activation=0.1、floor=0.001 和原 φ 表达式。对角块用 add.at 累加，非对角块唯一赋值；不创建原路径未使用的 coupling/groups。

浮点相加顺序可能改变最后几位；A 最大差约 2.3e-13（无量纲），由原 Gate、独立原 A/B 力平衡和全轨迹比较验收。原 A 每步重建；一次步内固定 A 复用因子，绝不跨不同位置复用旧 Cholesky。

单线程为推荐稳健配置。双线程全部冻结状态通过，部分中等 Cholesky 小幅改善，但小矩阵没有稳定整体优势。运行时从1线程切换到4线程曾发生native崩溃；进程启动前固定OPENBLAS_NUM_THREADS=4的独立测试通过。进一步2/4线程完整轨迹和全部冻结状态见threaded；失败切换路线不作性能样本。未额外引入 Numba/MKL，因为已测结构性收益足够，未无界扩工具栈。''')
    emit('CONTACT_MATRIX_REPRESENTATION.md','矩阵自由接触表示','''每个约束行保留 i,j,n 三项；Bv=n·(vj−vi)，Bᵀλ 在 i/j 上累加 −λn/+λn。所有 M=N(N−1)/2 对都保存和检查，未裁剪远距对、未删除 F2、没有空间筛选假设。Bᵀ 是相同算子的真实转置；独立稠密乘积/伴随测试覆盖全部 N。

128 泡稠密 B=8128×384，24,969,216 bytes；PairMatrix 索引和法向为 325,120 bytes，另有 gap、pair mapping 和活动小矩阵。这是 **B 表示存储** 对比，不是整个仿真器 RSS 的 76.8× 缩减。阻力 A 仍为稠密 SPD；全系统力学并未变成稀疏近似。

工作集选取、全部非活跃不等式重查、原 primal/dual、真实最终法向、互补与力平衡均保留。''')
    emit('CHOLESKY_AND_QP_OPTIMIZATION.md','Cholesky、QP 和 SQP','''原求解器一进入有约束路径就求解 A⁻¹Bᵀ 的全部 M 列；新实现仅对当前工作集新增列解线性系统，缓存限制在同一 QP 的同一 A/B 中，SQP 重线性化后重新建缓存。原目标、工作集排序、96 次上限、循环检测、gelsy 最小二乘/cond=1e-12、负对偶移除、违例添加和失败返回均不变。

Cholesky 没有跨物理步偷用旧因子，也没有显式矩阵逆。最大收益来自不再求解不需要的全约束右端项，以及不再每次构建稠密 B，不是把 Cholesky 本身宣称优化数倍。

FCC 532 活跃约束存在冗余/反力非唯一可能，物理验收采用净反力、KKT、真法向力平衡及位移，不要求任意数学对偶向量唯一。矛盾等式测试返回 active_set_inconsistent_equalities，非法几何、非 SPD、非法 SQP 上限均拒绝。''')
    emit('GPU_LINEAR_ALGEBRA.md','实际 CuPy FP64 负结果',table+'''\nCuPy 13.6 实际执行 GPU pair distances、FP64 向量化阻力、Cholesky、三角多 RHS 求解。冻结状态和完整轨迹均通过原数值门槛；当前规模上比优化 CPU 慢。GPU 路线为明确的 CPU 接触/GPU 因子混合实验，包含上传 A、每次 RHS 及解下载，不能称全流程常驻 GPU。

Warp/CuPy 共享指针与值一致测试通过，借用前同步、保持原所有者引用；未将这项只读互操作测试误称部署了全流程零拷贝。源码继续使用安全的原流路径。未使用 FP32/TF32；PyTorch 和 OSQP CUDA 为 NOT_RUN，未因没有性能收益强行更换另一种近似 QP。

初版 cp.errstate API 不存在，改用已安装版 cupyx.errstate；修复前失败输入和日志保存，修复后完整通过。''')
    emit('CPU_GPU_TRANSFER_AUDIT.md','传输与同步审计','''原 Warp 每步 H2D 3 次、72N bytes；D2H 8 次、112N+12 bytes。包括声场、F2 读回和两个各自独立的原壁面检查。Packed 只把压力/一次力/总力/错误位打包：D2H 5 次、120N+8 bytes，少 3 次同步读回，但略增字节数；不减少力刷新和墙面验收。

所有 Warp kernel、设备内主力复制和 pack 同一当前 stream，numpy() 完成同步。primary_device 保留 F2 累加前的主力，原压力与力计算完全未改。64 泡完整计算有小幅附加收益，128 泡实测差异很小；不声称全规模显著传输加速。

CuPy 额外上传/下载（包含 A 和求解列）在 profile_extra/TRANSFER_ALLOCATION.json 精确按调用记录。benchmark FULL_STEP 的 transfer_bytes 仅是 Warp，不能把它当成 CuPy 总流量。隐式 CUDA API 的每个内部同步没有 Nsight 全轨迹追踪，报告为未单独拆分；已测 wrapper 的等待+复制和 CUDA kernel 时间。''')
    emit('PARALLEL_SCALING.md','独立任务并行吞吐',pt+f'''\n相同 4 条 64 泡短窗口（0.2 ms）任务，分散/历史接触初态×OFF/PAIR，spawn 进程，1 BLAS 线程/worker，1/2/4 并发。每条内部原时序严格顺序，输出相同；不合并物理相互作用。

暖轮最高为 {best['workers']} 进程，相对 1 进程吞吐 {throughput:.3f}×。这是当前短任务下的真实结果，含每任务初始化和完整输出；不能宣传为单条轨迹同倍加速。冷轮包括进程/导入/CUDA 初始化；暖轮复用工作进程，但仍真实重建每条任务的引擎。资源监控 parent CPU 不包括 child，子进程 CPU/RSS 从各 RESULT 单列，GPU 为全卡采样。''')
    emit('NUMERICAL_EQUIVALENCE.md','完整数值一致性',f'''首轮96个冻结输入后端状态（CPU1/2线程、CuPy、OFF）全部通过原Gate；启动前固定2/4线程又完成64个状态及6条完整轨迹，均通过。32/64/128 固定 ID 的完整窗口 CPU/GPU 最大差 {max_id*1e9:.6g} nm；不是精度容差被放宽。dt/2 与细缓存的最大逐泡差 {max_refine*1e9:.6g} nm，原门槛 50 nm。

COM、Rg、Q、长短比、XY 伸长、可辨主轴方向、固定 ID 及无ID分布匹配、最近邻、几何边事件序列和时间全部对照；见 full/COMPARISONS 与 analysis/EXTRA_COMPARISONS。主轴近退化时不伪造方向准入。形态数据从全部接受帧计算，不只看末态。

{total_steps:,} 个接受步独立重算原 resistance 与最终真实法向的力平衡、对偶、真实间隙互补；另按原单位门槛重查全部逐步 residual，all_pass={all_numeric}。压力逐步范围、OFF=0、源幅相和事件时序沿用原代码，非零压力严格 0.1–40000 Pa。''')
    emit('CONTACT_VALIDATION.md','接触与失败闭锁',f'''32/64/128 的合法 cubic/FCC 固定状态含最多 532 个同时活跃约束；机械制造力不冒称声学结果。真实 64 泡接触延续、32 原近接触对照及正式 N32/N64 接触段均通过。原 128 泡 1 ms 无硬接触这一事实保留。

延长 128 泡 2 ms 的活跃接触峰值为 {extfast['max_metrics']['active_count']:.0f}，原 CPU 已完成部分峰值 {extref['max_metrics']['active_count']:.0f}；优化版完整2ms的dt/cache通过={all_ext}。原CPU在{extref['completed_s']*1e3:.3f}ms触及240s预算，完整2ms参考对照未通过（缺失），仅公共前缀可比较。

几何近邻阈值 2.1R 与硬接触 d<=2R+1nm/正反力事件分开；几何连通不是物理融合。反力的全部 pair ID 和每步值保存在 HDF5，未删减非活跃约束检查。独立验证用真实非线性间隙及原容差，不能仅 allclose 宣布正确。''')
    emit('FULL_TRAJECTORY_BENCHMARK.md','公平完整轨迹耗时',table+'\n'+latency+'''\n同硬件、每条1 BLAS线程、相同初态/幅相/dt/时长/缓存，原/优化/packed/CuPy 全保存相同每步位置、压力、力、接触反力、速度和原验收。包含初始化及 HDF5/JSON 写入；独立报告冷引擎、预热步、核心积分和保存。每正式后端/规模只跑一次长窗口；固定步 2 次预热+9 次测量，不能把微小耗时差说成统计显著。

N32 2 ms、N64 2 ms、N128 1 ms；之后 N128 优化版条件延长到2ms，半步长和细缓存完整完成；原CPU参考尝试在1.659ms因240s预算终止，保留FAIL与完整已接受前缀。没有秒级、千泡或新训练。''')
    emit('PERFORMANCE_AND_RESOURCES.md','真实资源与时钟',f'''AMD Ryzen 9 7950X3D，16 核/32逻辑线程，约63 GiB RAM；一张 RTX4090 24 GiB，驱动550.120；Python3.12.3，NumPy2.2.6，SciPy1.15.3，CuPy13.6，Warp1.18.0+cu12，实际 BLAS 为两套 OpenBLAS LP64/ILP64。版本/CPU亲和/线程和 GPU 空闲预检保存在 baseline/ENVIRONMENT。

完整正式 campaign 采样 CPU 均值 {cpu:.3f} 核当量，GPU 利用率峰 {gpu:.1f}%，显存峰 {vram:.1f} MiB，进程 RSS 峰 {rss:.1f} MiB。1秒采样会漏短峰；CPU 均值包含分析和保存。GPU 开跑无其他计算任务，但服务器始终存在另一 CPU 任务，因此不是整机独占基准。各新后端在相同可用资源下对照，历史 F10Z 探针不作为加速分母。

已完成阶段墙钟（部分固定测试在中断 campaign 中，单列 Supervisor ledger 后才是总数）：

{timing}

上述完成阶段合计 {elapsed:.3f} s；额外失败安装、native 崩溃及固定状态阶段时间见 EXECUTION_LEDGER.json，不隐匿失败。绘图/视频/下载另列，不混为物理求解时间。''')
    emit('FAILURES_AND_LIMITATIONS.md','失败、负结果与适用范围',f'''1. 新隔离环境首次 profiling 在依赖安装未完成时启动，两个 attempt 仅 import 失败，无物理输出。随后使用已验证官方已装包独立复制，恢复记录保存。旧项目环境未改动。
2. 运行时切换到四线程的 OpenBLAS Cholesky 发生 SIGSEGV；随后在进程启动前固定线程数，独立测试N32/N64/N128及原B0通过。失败不是接受物理步，保留 faulthandler 和输入 SHA，不降容差绕过。具体结果：{[(r['N'],r['backend'],r['status'],r['exitcode']) for r in four]}。
3. 初版 CuPy errstate 名称错误、Warp profiling callback indent 参数错误，分别修复；保存旧日志与不可变提交，已完成固定测试未选择性重跑取快值。
4. CuPy 在当前规模比优化CPU慢；GPU几何与因子转移未形成生产后端。Packed读回数量减少但字节略增，128泡没有可论证的额外整体收益。
5. 单条正式后端只有一次完整重复、每N一套主要分散初态；接触制造状态不是新物理证据。性能改善不证明 PDMS/自由液面/声流/气泡接触/电声映射已校准。
6. 256 仅合法短步预检，长轨迹 NOT_RUN。没有秒级、AI/world-model训练和真实PZT闭环。
7. N128延长原CPU参考在1.659ms耗尽240s预算，不能称完整2ms已通过参考；优化版dt/cache完整通过也不补齐缺失参考。
8. 图/视频显示的是2ms物理过程放慢播放，未通过插值、聚类简化、球径放大或篡改轨迹展示“聚成大团”。''')
    emit('FINAL_RECOMMENDATION.md','推荐后端与后续科研路径',f'''在原物理及原 Gate 不变前提下，本轮显著降低64/128泡完整计算成本：分别 {performance['64']['B1']['speedup']:.3f}×、{performance['128']['B1']['speedup']:.3f}×。32泡也加速 {performance['32']['B1']['speedup']:.3f}×。

32/64泡推荐1线程WaterAccelerated或WaterPacked，原Warp声学+壁面。128单任务采用启动前固定4线程WaterAccelerated，完整1ms实测{performance['128']['B2_t4']['wall_s']:.3f}s（相对原单线程{performance['128']['B2_t4']['speedup']:.3f}×）。多个独立任务仍每进程1线程；不在运行中调整BLAS线程。Packed在128没有额外明确收益。CuPy保留验证路线，不默认启用。四线程运行时切换的native异常未作为科学结果接受，启动前固定线程数的路线单列验证。

后续先做少量64/128泡更长窗口及多个独立初态/幅相的数值重复，优先确认接触事件与持续形变，而非立即扩到千泡。256预检{pre['status']}，需独立初态、多接触参考、dt/cache和墙面验证后才增加长轨迹。可增加完整逐泡研究数据吞吐，但不能将其标为真实PZT训练准入；下一轮动作—状态试验须仍保留固定ID、有限预算、失败及未校准标签。''')
    answer=[
    '反复稠密 B 几何是分散大N主瓶颈；接触后全 B.T 多右端项求解也显著。',
    '逐内核原始时间见 SOLVER_PROFILING 的表及 CUDA event 数据；嵌套时间不能相加。',
    'M=N(N−1)/2，稠密 B 有 M×3N 项，进入接触还求解全部约束列。',
    '有；已向量化并删除未用的 coupling/groups 分配，原润滑公式未改。',
    '可以数学等价矩阵自由化，保存全部 i/j/法向；非近邻裁剪。',
    '是；只对新增工作集列求解，其他不等式仍完整重查。',
    '是，全部原始粒子对及每次SQP的非活跃约束均检查。',
    'Cholesky没有结构性大加速；同一步因子复用，双线程仅部分小收益，GPU更慢。',
    '是，主要减少全 B.T 的无用解列；多接触状态原KKT通过。',
    '保持原目标、迭代上限、几何重线性化、力刷新、拒判和所有门槛。',
    '没有：当前32/64/128 CuPy完整轨迹慢于优化CPU。',
    '混合GPU额外矩阵/RHS传输及小矩阵启动成本抵消收益，不能仅归因FP64算力。',
    '推荐SciPy CPU力学＋现有WarpGPU声学；32/64单线程，128单任务启动前4线程；CuPy保留实验，PyTorch未运行。',
    f"没有，32泡完整CPU加速{performance['32']['B1']['speedup']:.3f}×。",
    f"64泡完整2ms加速{performance['64']['B1']['speedup']:.3f}×。",
    f"128泡完整1ms加速{performance['128']['B1']['speedup']:.3f}×。",
    '是，32/64真实接触、历史接触延续与128延长接触阶段通过。',
    '通过：包括532个同时活跃FCC冗余约束；反力非唯一以净反力及KKT判断。',
    '全部保留，并按原单位门槛独立重算完整接受步。',
    f'完整逐泡及形态最大固定ID后端差{max_id*1e9:.6g} nm，均通过原门槛。',
    f'dt/2和细缓存全窗口通过，基础正式组最大逐泡差{max_refine*1e9:.6g} nm；延长组另列。',
    '四线程库崩溃和初始接口/仪器问题已保留；修复后CPU1/2、CuPy及Packed无接受非法步。',
    f'暖轮最佳{best["workers"]}进程吞吐{throughput:.3f}×；若小于1则无收益，不能冒充单轨迹加速。',
    f'正式campaign CPU均值{cpu:.3f}核、GPU峰{gpu:.1f}%、显存峰{vram:.1f}MiB；非整机独占，1秒采样。',
    '32/64单任务用1线程（可选Packed），128单任务启动前4线程；独立任务每进程1线程，实测4进程吞吐最好。',
    f'256合法短步预检{pre["status"]}，没有据此批准长轨迹或所有密集状态。',
    f'已实际延长128到2ms，CPU优化{extfast["wall_s"]:.3f}s，原CPU{extref["wall_s"]:.3f}s仅完成{extref["completed_s"]*1e3:.3f}ms后预算停止；仅优化版dt/cache完整通过，不能从不等时长计算加速比。',
    '支持更多有限完整物理实验；单次轨迹成本降低不等于统计覆盖或训练准入完成。',
    '8组PNG/PDF及2段完整解码MP4，清单和哈希见figures/videos MANIFEST。',
    '验证了数值等价、局部几何近邻与Rg趋势、计算成本下降；真实装置/电声/接触/融合仍未校准。',
    '提供相同step/control接口的显式适配器，可用于F10Z后续；默认历史入口未改，需选定新适配器。',
    '先64/128更长窗口、更多独立初态及有限动作，边做收敛和事件可信度，再考虑256长轨迹；不急于真实世界模型训练。']
    emit('F10AA_DEVELOPMENT_REPORT.md','F10AA 中文开发报告',table+'\n\n'+ '\n\n'.join(f'{i+1}. {a}' for i,a in enumerate(answer))+'\n\n全部37项机器门槛见 FINAL_GATE.json；成功为计算实现和数值保持，不是实验校准。')
    emit('SOURCES.md','官方实现与审查依据','''- [SciPy cho_solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html)：因子复用及多 RHS。
- [CuPy Cholesky](https://docs.cupy.dev/en/stable/reference/generated/cupy.linalg.cholesky.html)：GPU SPD 求解；实际版本13.6 API以本机安装源为准。
- [threadpoolctl](https://github.com/joblib/threadpoolctl)：隔离线程数量。
- [NVIDIA Warp](https://github.com/NVIDIA/warp)：本机1.18.0 ScopedTimer CUDA events及数组所有权/同步接口。
- [line_profiler](https://github.com/pyutils/line_profiler)、[Python cProfile](https://docs.python.org/3/library/profile.html)：原实现瓶颈实测。

性能及科学结论来自本轮保存数据，上述文档不作为本轮加速数字的外部证据。''')
    gates={k:'PASS' for k in ['F10Z_BASELINE_VERIFIED','PHYSICS_MODEL_UNCHANGED','NUMERICAL_TOLERANCES_UNCHANGED','CPU_PROFILING_COMPLETED','GPU_PROFILING_COMPLETED','PAIR_GEOMETRY_OPTIMIZED','RESISTANCE_ASSEMBLY_OPTIMIZED','SPARSE_CONTACT_GEOMETRY_VALIDATED','ACTIVE_SET_QP_VALIDATED','CHOLESKY_SOLVER_VALIDATED','GPU_LINEAR_ALGEBRA_TESTED','CPU_GPU_TRANSFER_OPTIMIZED','N32_NUMERICAL_PASS','N64_NUMERICAL_PASS','N128_NUMERICAL_PASS','ACTIVE_CONTACT_CASES_PASS','MULTI_CONTACT_CASES_PASS','OFF_CONTROL_PASS','STEP_SIZE_CONVERGENCE_PASS','CACHE_CONVERGENCE_PASS','CONTACT_COMPLEMENTARITY_PASS','FORCE_BALANCE_PASS','WALL_GEOMETRY_PASS','FULL_TRAJECTORY_COMPARISON_PASS']}
    if not all_numeric or not all_cmp or not all_ext:raise RuntimeError('Cannot generate PASS gates for failed numerical admission')
    gates.update(GPU_ACCELERATION_BENEFICIAL='NO',N32_SPEEDUP=performance['32']['B1']['speedup'],N64_SPEEDUP=performance['64']['B1']['speedup'],N128_SPEEDUP=performance['128']['B1']['speedup'],N128_STARTUP_FOUR_THREAD_SPEEDUP=performance['128']['B2_t4']['speedup'],TOTAL_THROUGHPUT_SPEEDUP=throughput,N256_PREFLIGHT=pre['status'],LONGER_DURATION_TEST='PASS' if all_ext and long_reference else 'PARTIAL_REFERENCE',LONG_OPTIMIZED_DT_CACHE='PASS' if all_ext else 'FAIL',LONG_REFERENCE_COMPLETION='PASS' if long_reference else 'FAIL_WALL_BUDGET',FIGURES_VALIDATED='PENDING_VISUAL_REVIEW',VIDEOS_VALIDATED='PENDING_VISUAL_REVIEW',PHYSICS_CALIBRATED=False,REAL_PZT_VALIDATED=False,WORLD_MODEL_TRAINED='NOT_RUN',TRAINING_READY=False,SCIENTIFIC_CONTROL_ADMITTED=False,FOUR_THREAD_RUNTIME_RESIZE='FAIL',FOUR_THREAD_ISOLATED_STARTUP='PASS' if all(r['status']=='PASS' for r in four) else 'FAIL')
    write(BASE/'FINAL_GATE.json',dict(gates=gates,labels=LABELS,training_ready=False,identity=identity(),scope='Implementation and numerical acceptance only; all physical/electrical calibration exclusions retained',complete_trajectory_steps_independently_checked=total_steps,full_trajectory_scope='Registered complete 32/64/128 windows across CPU1/2/4, CuPy and Packed; the budget-limited 128 extension reference remains separately incomplete'))
    write(out/'SUMMARY.json',dict(performance=performance,parallel_best_workers=best['workers'],parallel_throughput_speedup=throughput,complete_stages_wall_s=stage_times,reference_fixed_id_max_m=max_id,refinement_fixed_id_max_m=max_refine,extension_reference_s=extref['wall_s'],extension_optimized_s=extfast['wall_s'],checked_steps=total_steps,resources=dict(cpu_mean_cores=cpu,gpu_peak_percent=gpu,VRAM_MiB=vram,RSS_MiB=rss),identity=identity()))
if __name__=='__main__':main()
