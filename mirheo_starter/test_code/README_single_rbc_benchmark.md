# 单红细胞简单剪切：本轮实际结果

**没有完成双方从零到共同终点的完整成本比较。** 两方均实际推进了液体和一个原生可变形细胞，HemoCell 两次主运行完成 Γ=4；Mirheo 两次主运行、半步检查和额外受限队列诊断均以原生粗碰撞候选溢出退出。工作流 `PARTIAL`，模型可比性 `PARTIAL`，`BENCHMARK_SCREEN=FAILED`，`qualified_speedup=null`。这些数据不能证明同等物理结果下哪套软件更快；真实研究适用性 `NOT_VALIDATED`，人工验收 `PENDING`。

本轮 campaign：`rbc_shear_20260910`。新增 3600 求解墙钟秒已获用户“批准”，见独立 `authorization.json`。所有求解已经停止，剩余额度不触发自动寻参或重试。

## 实际终点与成本

| 任务 | 状态 | 完成终点或最后存帧 Γ | 进程墙钟 / 秒 | 完整任务含必要分析 / 秒 |
|---|---|---:|---:|---:|
| HemoCell 冷启动 1 | 完成 | 4.0 | 52.437 | 58.172 |
| HemoCell 冷启动 2 | 完成 | 4.0 | 55.055 | 60.783 |
| Mirheo 冷启动 1 | 原生错误退出 | 0.6 | 147.767 | 未获得 |
| Mirheo 冷启动 2 | 原生错误退出 | 2.9 | 620.421 | 未获得 |
| HemoCell 半步 | 完成 | 4.0 | 106.572 | 112.107 |
| Mirheo 半步 | 原生错误退出 | 0.7 | 355.336 | 未获得 |
| Mirheo 受限队列长诊断 | 原生错误退出 | 3.4 | 786.973 | 未获得 |

HemoCell 两次端到端均值 **59.477 秒**，样本标准差 **1.846 秒**，平均 **14.869 秒/累计剪切应变**。确定性轨迹相同，只有两次独立进程重复，不能给出精确 95% 区间。HemoCell 每次含新建输入、进程启动、液体/细胞初始化、零剪切松弛、流体+膜+IBM 演化、输出和共同数值分析。分析在保存结果后独立计时，端到端是这些互不重复阶段的实测和；并非一次连续计时器测得的总等待时间。

Mirheo 的数值检查和失败进程成本全部保留在 `results.json` 的各任务 `costs`，但 `total_s` 仅表示**该失败尝试加核查的成本**，没有冒充到 Γ=4 的完整任务时间。最后存帧是已完成推进的下界，不是精确崩溃时刻；失败时尚未返回的计时块保留在未分解开销内，不补造步数。

首次新案例编译/配置共 **69.474 秒 CPU 进程墙钟**（8 条记录，链接既有 HemoCell 库）。下表另外列出非主运行检查，避免把额外诊断隐藏在最优成绩里：

| 检查类别 | CPU 求解墙钟 / 秒 | GPU 任务墙钟 / 秒 | 本版必要分析及新建输入 / 秒 |
|---|---:|---:|---:|
| 主运行前的预检、材料/膜匹配及空通道 | 82.950 | 555.656 | 7.554 |
| 两方自身半步检查 | 106.572 | 355.336 | 8.602 |
| 额外受限队列诊断（短测和失败长测） | 0.000 | 835.508 | 10.389 |

非主运行求解合计 **1936.023 秒**；加上两方四次主运行，全部求解计费 **2811.703/3600 秒**，剩余 **788.297 秒**。CPU 秒不乘 MPI rank 数；GPU 任务秒包括配套 CPU 和准备阶段。所有求解与编译按时间戳顺序执行，无重叠，最长单任务 786.973 秒。旧 SDPD 的授权、余额和科学结论保持原样。

最终 HTML 的生成耗时单列于 `data/.../review_latest.json` 的 `html_render_s`，证据收集/缓存加载为 `collection_wall_s`。历次真实分析计时保存在 `analysis_events.jsonl`，重新分析不增加冷启动重复数。软件已安装后的本轮测量没有重装成本；没有可靠记录的历史安装人工时间和网格转换零散人工时间记为未知，不估造 Codex 思考时间。

## 共同问题及对应路径

| 项目 | HemoCell | Mirheo DPD |
|---|---|---|
| 液体 | 原生 LBM，2 MPI rank，double | 原生 DPD，1 GPU compute + 1 postprocess rank，single |
| 细胞 | High Order 原生材料力 | WLC + Kantor 原生膜力、Velocity Verlet |
| 双向作用 | 实际 `hemocell.iterate()` 的膜力及 IBM 反馈 | 膜两侧原生 DPD 力及 `bounce_back`，两边液体参数相同 |
| 移动边界 | X/Y 周期；Z 两侧速度边界 | MovingPlane、同速壁粒子及原生壁碰撞 |
| 新增适配 | 独立 C++ 案例、公共参数与实际网格/流场输出 | Python worker、同一 OFF、静止到移动壁状态交接、统一输出 |

X 流向、Z 梯度，X/Y 周期长度 24，有效壁距 H=24，下壁 −0.24、上壁 +0.24，γ=0.02。单细胞中心 (12,12,12)，使用同一导出几何与取向，642 顶点、1280 三角面。参考 A0=129.211626987、V0=81.116318378，a=sqrt(A0/4π)=3.206607954，减缩体积 0.587329785；H/a=L/a=7.484544523。它是有限盒问题，没有称为无限剪切流。

HemoCell 的有效墙面为 z=0、24，格点 48×48×49、dx*=0.5。Mirheo 总盒子 24×24×26，墙面 z=1、25，输出减去 z 偏移 1。实际 helper 只在 X/Y 周期，本轮显式转换其相反壁速符号。HemoCell 以角度 (90°,0°,0°) 生成真实初态，导出的同向 OFF 在 Mirheo 使用标量在前的单位四元数。

每次从无液体、无壁粒子、无细胞状态的新进程开始。先零剪切耦合松弛 t*=5，再共同瞬时启动壁面；主目标 t*=200、Γ=4，每 ΔΓ=0.1 输出。HemoCell 每主运行实际 5535 次 iterate，其中 135 次准备；dt*=1/27，半步为 1/54。Mirheo 主计划 205000 演化步另含 2000 壁准备步，dt*=0.001；半步 0.0005。半步计划终点 Γ=2、共同检查窗 Γ=0.5～2，原配置实际只存到 Γ=0.7。适配器只在计时块边界额外同步，不新增每步强制同步；保留 Mirheo 原生调度。

Mirheo DPD 参数来自固定版本原生案例：n=8、m=1、kBT=1、rc=1、a=10、gamma=10、power=0.5。唯一液体点集按膜内/外分类，内内、外外、内外使用同一 DPD 参数，λ=1；没有叠加生成两份流体、额外恒温器或对流体/RBC 指定剪切轨迹。壁粒子生成、静止准备及同进程状态交接均计时。交接保存实际液体位置/速度和膜状态，保持初始应力参考；阶段内自动分类修正为 0。阶段交接重新分类和探针 ID 变化单列，不能用重新分类声称无泄漏。

统一数值嵌入为 L0=1e−6、E0=4.100531391e−21、M0=1.25e−16、T0=sqrt(M0 L0²/E0)，用于 HemoCell 的 SI 接口。结果按共同 t* 和 Γ 报告，不称为已标定的血液、生理秒或 PDMS 实验。

## 实际质量结论

- DPD 独立体力驱动双向周期 Poiseuille 得 ν*=1.120197342，HemoCell 使用这一冻结目标；半步 ν*=1.052089394，差 6.4736%，未过事前 5% 门槛。短采样有相关热涨落，不能把差异全归为精确的 dt 偏差。Couette 直线没有被用来反推黏度。
- 两边空通道 Γ=2～4 平均速度、剪切率、质量和近壁检查通过；DPD 平均速度 L2 误差约 3.21%，有效剪切率误差约 1.31%。液体材料的半步检查仍未通过。
- 独立原生仿射力探针得到有效剪切响应 20.06819 与 21.30977，剪切约差 6.2%，等体积伸长最大差 17.82%，未过冻结 10% 膜响应门槛。这个总增量力功代理不是精确连续膜 Gs。没有用最终剪切曲线反复拟合参数。
- Re≈0.18358；代理 Ca≈0.02864/0.02697。HO 的局部参考曲率与 Kantor 常角弯曲不同，HemoCell 连续弯曲 B 未核实。Mirheo 膜质量/排开流体质量≈0.9893，不能认为惯性可忽略。LBM 无分子热涨落，DPD 保留噪声；整体热能比约 0.0046～0.0048 不证明局部噪声无关。
- 零剪切准备后的对应顶点 RMS/a 约 7.7%，超过冻结 2% 门槛。即使参考网格相同，两边准备后的细胞也没有达到共同形状质量要求。
- HemoCell 主运行保存帧无非相邻自交、无穿墙；最大 A/V 漂移约 0.1893%/0.0298%。其半步在完整 Γ=0.5～2 窗口平均 D 差 0.0001831，通过当前时间初筛；不证明空间收敛。最终 Γ=4 的 D≈0.46606、θ≈−69.44°。
- Mirheo 原配置两次主运行分别有 40/3456 与 376/7872 探针点帧成员不符，半步为 72/3648。第 2 次 A/V 最大漂移约 2.20%/2.33%。不能宣称无跨膜迹象。受限队列短测为 0/3264；长诊断最终仍有 13/8832 不符、存帧累计 21 对非相邻三角面相交及约 2.35% 面积漂移，未解决稳定性问题。
- 两侧均保存实际膜与局部液体变化，支持经源码核实的双向反馈路径。HemoCell 有/无细胞差使用 Γ=2～4；Mirheo 原配置第 2 次使用 Γ=2～2.9。每张差值图内部窗口相同，两软件窗口不同，DPD 差值还含热噪声，不把它作为无噪声因果力测量或跨模型吻合证明。

形变 D 统一从固定 X-Z 平面顶点协方差特征值平方根计算；θ 取轴向角，主轴比 <1.05 记无效。只做周期整数展开，不作旋转或时间拟合。严格不渗透、空间收敛、完整释放耗散等效及运动周期均未建立。现有数据是实际执行与失败证据，不能标为科学筛选 PASS。

## 原生失败与受限诊断

原配置两次主运行分别报粗碰撞候选 9185>6400、8244>6400；半步报 13057>6400。额外诊断只在子进程设置 `CUDA_DEVICE_MAX_CONNECTIONS=1` 与 `CUDA_DEVICE_MAX_COPY_CONNECTIONS=1`，不修改原生库，长测仍报 7787>6400。任务均以原生错误退出，未耗尽单任务或总预算。

共享碰撞缓冲及多流调度是源码支持的待验证假设，尚未证明根因；几何与力学异常也必须考虑。短测成功不是修复，长诊断不算原配置的第二次成功重复。长诊断失败后，条件计划中的后续半步没有启动。`native_failure_diagnostic_plan.json`、`runtime_candidate_frozen.json` 和 `native_failure_final.json` 保留计划与停止理由。

最值得做的下一项工作是：在后续独立任务中，用最小可复现的膜两侧碰撞案例定位候选溢出与跨膜/自交的根因，再决定原生路径是否可用；目前不应直接增加长测时间或给出速度排名。

## 实际命令与核查文件

```bash
cd /home/lzy/projects/mirheo_starter
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

.venv/bin/python -B -m py_scripts.benchmark_single_rbc --help
.venv/bin/python -B -m py_scripts.benchmark_single_rbc \
  --config py_scripts/single_rbc_benchmark.yaml --preflight-only
.venv/bin/python -B -m py_scripts.benchmark_single_rbc \
  --config py_scripts/single_rbc_benchmark.yaml --execute
.venv/bin/python -B -m test_code.review_single_rbc_benchmark \
  --config py_scripts/single_rbc_benchmark.yaml --open
.venv/bin/python -B -m unittest test_code.test_single_rbc_benchmark -v
```

帮助、预检和 review 不启动 GPU 求解。当前 `--execute` 已实际核验并回读全部 8 个冻结正式任务的封存结果（包括失败），不会自动重跑；两类账本 SHA 和 26 条求解尝试数不变，见 `readonly_cli_final.json`。review 缓存不是新的冷启动重复。`--prepare-cpu` 是已完成的首次新案例编译/短测入口，不是普通查看选项。

- [中文离线 HTML](outputs/single_rbc_benchmark/rbc_shear_20260910/single_rbc_review.html)：内嵌 Plotly 与实际数据，15 张图、同步真实网格、墙面/方向、播放/暂停、重复及单列诊断；缺失终点留空。
- [结果 JSON](outputs/single_rbc_benchmark/rbc_shear_20260910/results.json)：每次运行、全部几何/液体/守恒、分项计时、实际环境及判断。
- [真实浏览器核查](outputs/single_rbc_benchmark/rbc_shear_20260910/browser_delivery_verified/browser_check.json)：以该报告状态与当前 HTML SHA-256 同时吻合为准，人工始终 PENDING；同目录含 overview、cells、dynamics、costs 截图。
- `/home/lzy/projects/mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910/`：共同 OFF、单位/参数、冻结规则、授权、数值证据、26 项 CPU 测试、保护与最终交付记录。
- `/home/lzy/projects/mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910/`：全部 34 条执行记录（26 求解 + 8 构建），配置/源码快照、原始 vertices/profiles/local_flow/timings、探针、控制台、资源和输出 SHA。
- `/home/lzy/projects/hemocell_starter/cases/single_rbc_shear_benchmark/`：新增真实 HemoCell C++ 案例及许可证。

旧环境/数据保护覆盖 20157 个文件，实际校验 PASS。26 项 CPU 测试通过，涵盖单位、边界、几何/周期/ID、原生细胞路径、成本、资格门控、缓存、只读查看和早期失败空 CSV；最终真实 Chrome 浏览器 19 项检查全部 PASS（含实际播放、第二次重复、诊断存帧、无外部请求及图例排版），与当前 HTML 哈希吻合；浏览器功能核查不代表物理筛选通过。

## 固定版本及许可

HemoCell `5a410848bd5c57d5ae1c171112e78eab4a82e650`，Mirheo `8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf`。Palabos 是固定 `05712164d940a42e06afdd705249912fa0c49f14` 的解压归档与现有 HemoCell 补丁，**不是独立 Git checkout**。此前对该目录做 `git -C rev-parse HEAD` 会继承父 HemoCell 的 HEAD；原错误记录保留，明确更正在 `environment_authorized_start_corrected.json`，不是修改真实源码版本。

核实的 Palabos 归档 SHA-256：`5a9c4f22c169ad16de259a72ef2f128b233afaf81c251e6c3baa8eb34d14838f`；现有补丁 SHA-256：`f01533c608ce82bf75e771cd542b3b14f326b5fda7ee60233e1943f58ef64c3e`。HemoCell 库 SHA-256：`27ddaa6e67f1fbe169345f02df17107a8840eed46827f8602f4e7eaa9f9cb0e2`；Mirheo 库 SHA-256：`1bab2b922b43a33dd0ffc84211234f3cf75f8edc852e30adfcd17f68c8c78c87`。

CPU 可见 24 逻辑核，授权启动时 WSL 可用内存约 3.44 GiB，RTX 4060 Laptop 总显存 8188 MiB、空闲 6204 MiB；逐任务 RSS/设备显存采样另存。设备显存是整卡采样，不冒充精确进程峰值。HemoCell 派生案例附 AGPL，Mirheo 派生 worker 附 MIT。没有重新安装、修改原生 Mirheo C++/CUDA 或旧兼容补丁、重编译 Mirheo、覆盖官方示例或提交/推送本轮文件。
