# Mirheo 云端完整单红细胞验证：实际结果

没有完成 Γ=4。准备阶段发生原生碰撞候选溢出并退出，未进入剪切；本轮没有重试。

CLOUD_RBC_RUN_COMPLETE = NOT_COMPLETE
CLOUD_RBC_NUMERICAL_SCREEN = FAILED
RESULTS_RETURN_VERIFIED = PASS

JOB_ID：cloud-rbc-full-20260911T000225Z
本地原始归档：/home/lzy/projects/cloud_results/cloud-rbc-full-20260911T000225Z
云端任务：/workspace/bloodflow/cloud_runs/cloud-rbc-full-20260911T000225Z/

构建 mirheo-sm120-20731713865ae510；实际计算库 /workspace/bloodflow/.venv/lib/python3.12/site-packages/libmirheo.cpython-312-x86_64-linux-gnu.so；SHA-256 d45b4fd1b4498b6f365a012e758eec35bc58a7844dabc66f22e617a7b8cf84ee。两个 MPI rank 的实际加载记录一致，单 GPU、sm_120、原单精度。

修复 local_before_halo_v2；运行方式每阶段一次连续 u.run；本次顶层 spec 的 bouncer_policy=shared，同一个 bounce_back bouncer 绑定两侧液体。嵌套 config.repair 的 independent_per_pv_control 是保留的历史配置字段，worker 实际读取顶层值。原生候选容量固定，未扩容。

从零初始化：液体 110592（outer 109950、inner 642），粒子质量 1；膜 642 顶点 / 1280 三角面。有效域 24³、X/Y 周期、Z 壁面。dt=0.0005，实际 float32 dt 见碰撞原始记录；准备上限 60000 步、正式剪切计划 400000 步（γ̇=0.02，Γ=4）。WLC/Kantor、ks=3.0、kb=8.0，DPD kBT=1；完整参数在 actual_parameters.json 与 full_spec.json。

准备未通过：u.run(60000) 未正常返回，成功返回的准备步数未知（null）。定期完整帧 46，最后准备步标记 45000、名义 t*=22.5，只是进度下界。正式剪切 0 步、Γ=0；4000 步壁粒子准备不计入 RBC 正式应变。无成功 completion，也没有准备末态或剪切末态。

定期保存膜帧：面积最大相对漂移 1.962563%，体积 0.791998%；确证非相邻自交帧 0，近接触帧 0，未见定期帧穿墙、退化或 NaN/Inf。该结论仅覆盖这些帧，不包括报错瞬间快照。

原生失败：准备阶段零基步标签 45040、原生时间标签 22.5200010696；outer/local 粗候选 18074 > 容量 6400。已保存 6400 个候选对，其内部重复 0；不能据此判断被截断部分。fine_count=-1 表示该失败调用尚未执行细筛，不能写成零碰撞。

故障快照已出现严重膜损坏：7 个顶点的实际 Z 越过壁面，最小有效 Z=-126.7130126953125。膜最大速度 918972，outer 液体最大速度 60336.7。周期展开失败（INCONSISTENT_UNWRAPPING_OR_CELL_SPANS_BOX）；至少 121 条边的最短周期像 WLC 延伸值已≥1，最大下界 197.392。

故障快照位于该步积分之后、此 bouncer 碰撞处理之前，未当作正常完成帧。Z 为非周期轴，越墙证据不依赖 X/Y 展开。由于无法构造一致闭合周期曲面，该瞬间自交数、物理 A/V 和成员分类均标为无法可靠判定，不用定期帧的零自交替代。

成员检查仅覆盖真正初态：内外各 96 个按排序 ID 等间距抽取的探针，共 192 点次；确证不符 0，近膜不确定 0。没有故障前完整探针轨迹，不能证明全程或所有粒子不可渗透。strict_impermeability=NOT_VERIFIED。未做事后成员修正，静壁到动壁交接 NOT_REACHED。

速度证据：最后定期帧（beforeForces，45000）膜最大速度约 43.31，之前定期帧最大约 5.99；原生 Stats 在其 afterIntegration 样本 t*=22.5005 记录 inner 最大速度约 94.26、kBT≈5.658，膜最大速度≈59.15。原生 kBT=m·mean(|v|²)/3 包含整体流动动能；不是扣除流速的温度，也不是开尔文。

原生 Stats 各保存 46 行，记录时液体总数/质量不变。局部流场保存 45 个可读、有限的 8³ 原生分箱场；尚无正式剪切响应。因调用未返回，vertices.csv、moments.csv、profiles.csv、local_flow.csv、timings.csv 只有表头；它们不代表零物理量。实际膜轨迹来自原生 HDF5/XMF，派生 geometry.csv 有 46 行，原生统计与局部场保留在 simulation 下。

原因范围：已确认候选容量溢出，且失败调用的碰撞处理之前已存在速度与几何爆发。现有数据不足以确定最初触发来自积分、膜力、共享 bouncer 状态或其他机制；不能仅把增大容量视为修复，也不把该结果归因为 GPU 硬件。

求解进程退出码 255；墙钟 36.903925 秒（1800 秒授权的一次尝试，已消耗 1/1 次）。包括初始化与必要准备，不等于 GPU 活跃时间或 Vast 租赁费用；剩余额度不自动授权第二次运行。setup 已完成 1.171085 秒，其中壁粒子准备 0.809405 秒，不能重复相加。phase_timings 的 relaxation_s=0 是未返回调用的占位值，实际准备阶段耗时未单独闭合计量，不能解释为零耗时。云端事后 CPU 分析 7.606645 秒另计；本地事后分析时间记录于 postmortem.json。

资源覆盖不完整：59 次进程采样只包含 runner 7455 与 mpirun 7458，漏掉实际 MPI rank 7495/7496。因此所见 RSS 峰值 83374080 字节仅为这两个管理进程，不能当作整个求解任务峰值；rank 的 CPU 时间与 RSS 未测。实际 rank 身份由各自 rank_0.json/rank_1.json 证实。GPU 为全设备采样，观测显存峰值 571 MiB；cgroup 为全容器。事后只读复核四个登记 PID 均已不存在。该采样遗漏仅记录为限制，本轮没有为补测重跑。

回传清单通过：970 文件、65021041 字节（不含清单/校验收据自身）；SHA-256 8266cee0a8052247a1ee3eb3d5e4c2ee3fea192d306615287f7b3bb6791be1b9。RESULTS_READY、清单、JSON/CSV 与原始 HDF5 读取核查通过。远端 report_zh.md / numerical_checks.json 中回传 PENDING 是封包前状态；本地 LOCAL_ARCHIVE_VERIFIED.json 和本派生报告给出实际 PASS。

原始云端封包及本地归档保持不变；本地补充分析、HTML 和浏览器证据位于同级 -review 目录。浏览器实际检查见 browser_final/browser_check.json（仅该文件及对应 HTML 哈希一致且状态 PASS 时才算通过），人工验收 PENDING。

PHYSICAL_MODEL_VALIDATION：既有材料 NOT_MATCHED；本轮材料标定、空间/时间收敛及 HemoCell 比较 NOT_TESTED。qualified_speedup=null。

没有重新编译、安装依赖、修改物理参数/原生修复、运行 HemoCell、推送 Git、删除云端数据或停止 Vast 实例。

查看/下载不会启动求解：
cd /home/lzy/projects/cloud_compute
python3 -B cloud_run.py status cloud-rbc-full-20260911T000225Z
python3 -B cloud_run.py fetch cloud-rbc-full-20260911T000225Z
python3 -B cloud_run.py verify cloud-rbc-full-20260911T000225Z
python3 -B cloud_run.py rbc-review cloud-rbc-full-20260911T000225Z
原 rbc-full --preflight-only 和 rbc-full --new --execute 已实现；本轮授权与启动收据已消费，不能用执行命令重复本轮任务。

本地回传校验：RESULTS_RETURN_VERIFIED，970 文件，65021041 字节。
