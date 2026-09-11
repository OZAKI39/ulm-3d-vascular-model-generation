# Mirheo 完整单 RBC 云端验证：执行前核查

本轮正式任务尚未启动。没有新的 GPU 求解、没有新的 Γ=4 结果、没有为本轮创建远程 JOB_ID。当前状态为等待一次独立的云端长任务授权。

## 修复版对应链

- 本地原生源码：`/home/lzy/projects/mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source`。
- 由上传 manifest 解析的云端原生源码：`/workspace/bloodflow/uploads/20260910T212439Z/source/mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source`。
- 原生源文件核对：2324 个文件的本地原件、冻结副本和云端上传副本一致。原有零字节 `googletest/googlemock/build-aux/.keep` 由冻结上传的 `excluded_files.json` 明确排除；没有原生代码差异。
- 原生源码集合哈希：`2b30d14be357b1bcd7a63805cc164f4907d291b38fd8ae4d4db44980b3f0ac90`。
- 构建：`mirheo-sm120-20731713865ae510`；sm_120、96 个 CUDA 编译单元、单精度、Release。云端兼容补丁只涉及归档版本信息、预编译库安装路径和 CMake 对 sm_120 的识别。
- 实际云端库：`/workspace/bloodflow/.venv/lib/python3.12/site-packages/libmirheo.cpython-312-x86_64-linux-gnu.so`。
- 库 SHA-256：`d45b4fd1b4498b6f365a012e758eec35bc58a7844dabc66f22e617a7b8cf84ee`。
- 参数来源：A6_continuous_preparation_half_dt/spec.json；对应交付记录为 `GATE_A_UNRESOLVED_HALF_DT_PREPARATION_COMPLETED_QUALITY_FAILED`。新完整工况由这一组参数和本轮 Γ=4 目标组成；目前不存在已经复验通过的完整 Γ=4 配置。

## 五项版本结论

1. 原生修复把 local 碰撞放在 halo 碰撞之前，避免共享 bouncer 的临时状态被交错调用覆盖；修复版还保留候选计数/容量和故障快照。没有扩容补丁，也没有证据能据此宣布长期稳定。
2. Python 必须配套使用 continuous_protocol.evolve_stage：每阶段一次连续 `u.run`，由原生插件输出。该固定版本在每次 `u.run` 时重新建立任务、重新分离粒子向量和清理对象力，不能用旧分块 worker 做连续状态观测。
3. A6 实际 spec 的 bouncer_policy 是 shared，内外两侧都保留 bounce_back。独立 bouncer 是旧对照分支，并非当前已确认的必要条件。本轮不擅自改为独立 bouncer。
4. 采用 A6 的 dt、膜/液体质量、作用参数、共同网格和准备上限；只把本轮目标设为一次 Γ=4。完整 worker 的 CPU 路径已检查，新增原生输出组合与完整两阶段云端路径仍未经过本轮 GPU 实跑。
5. 已复验：库来源与 sm_120、双 rank MPI/CUDA smoke、短程连续准备和修复后的调用顺序。仍未解决：长准备质量及其失稳初始触发因素、材料匹配、完整端点和收敛。A6 的面积漂移曾为约 2.2946%，超过原有 2% 筛选线，准备残差也未通过。

## 冻结运行参数与门槛

| 项目 | 本轮值 |
|---|---:|
| 有效流体域 | 24 × 24 × 24；X/Y 周期、Z 方向两壁 |
| dt | 0.0005 |
| 原生 float32 dt | 0.0005000000237487257 |
| 壁粒子准备 | 2 × 2000 步；不计入 RBC 应变 |
| RBC 准备上限 | 60000 步，名义 t*=30 |
| 正式剪切 | 400000 步，名义 t*=200、Γ=4 |
| 目标步原生时间/应变 | 200.00000949949026 / 4.000000189989805；浮点表示差异，实际值另记 |
| 剪切率与壁速 | 0.02；±0.24 |
| 液体 | 110592 粒子；质量 1；数密度 8；DPD kBT=1、rc=1、a=10、gamma=10、power=0.5 |
| 膜 | 642 顶点、1280 三角面；质量 1；WLC + Kantor；ks=3、kb=8、x0=0.457 |
| 膜其他参数 | mpow=2、ka=279.693130111、ka_tot=0、kv_tot=1188.746、theta=6.97、gammaC=0、kBT=1 |
| 耦合 | 双向真实液体—可变形膜；膜耦合 a=0、gamma=10；两侧 bounce_back |
| 成员修正 | 每阶段原生初次分类，阶段内 correct_every=0 |
| 准备膜输出 | 每 1000 步 |
| 正式膜输出 | 每 10000 步，即 ΔΓ=0.1 |
| 最多膜帧 | 102，含两个阶段实际返回末态 |

准备判据沿用 material_matching.json：相邻非重叠 2 单位时间窗，去中心平移，不拟合旋转；最后连续 3 个窗间残差≤0.002/a/时间、倾角变化≤0.5°/时间。A/V 变化≤2%，几何、粒子数/质量和末态温度同时检查。参考形状差异按原 preparation_quality.py 作为诊断，不冒充新的 HemoCell 匹配。

A/V 轻微超限会记录筛选失败，准备最多运行到固定上限；准备窗、几何或末态检查未通过就不进入剪切。正式剪切期间轻微筛选超限只记录，不修改参数。原生 CUDA/MPI/碰撞错误、非有限值、明确非邻面自交、穿墙、退化三角形、WLC 非法伸长、液体数量变化以及时间/磁盘/内存保护线触发停止。没有自动重试、半步或额外准备。

## 观测、交接与输出范围

每阶段的原生初态保存实际液体/膜状态；膜按固定间隔输出。完整液体仅保存两个阶段的初态，不输出全粒子轨迹。成员检查每侧最多 96 个确定性 ID，覆盖阶段初态和正常返回末态；报告近膜不确定带和第二方法歧义，不能把零探针不符解释成所有粒子严格不渗透。

静止壁面到移动壁面沿用原有两个顺序协调器。保存/恢复液体与膜坐标、速度、壁粒子位置和原参考网格；粒子 ID、初始分类、膜 oldPositions、力/任务结构和原生随机状态重新生成。保留排序 ID 映射、恢复误差和分类改变前后核查。这是部分状态交接，不能称为完整状态无缝连续；没有加载 smoke 末态或 checkpoint。

失败时保留原始日志与数据；失败调用只报告最后完整帧的进度下界，正常返回的阶段另记真实步数。原生 Stats 的 m·mean(|v|²)/3 含流动动能；末态温度另外扣除 24 个 Z 分箱的平均速度，以 DPD kBT=1 判断，不使用 SDPD 298.15 K。

最多预留/保护输出 536870912 字节（512 MiB），包括初态液体、膜轨迹、原生碰撞诊断、固定局部场、统计、配置和报告；达到输出保护值时停止。当前云端空闲 18811748352 字节，扣除预留仍大于 5 GiB；运行到 3 GiB 空闲保护线时停止。

短 smoke 的 250 步准备耗时外推约 358.2 秒。考虑新观测与 CPU 检查，计划范围约 360–1500 秒；这不是长任务速度实测或完成保证。单次远程求解进程硬上限为 1800 秒。setup 已包含壁准备；原生输出已包含在连续阶段耗时中；交接、末态输出、CPU 分析另记，不重复相加。GPU 显存是全设备采样，任务墙钟不是 Vast 租赁账单。

复用通过验证的 `HWLOC_COMPONENTS=-opencl,-gl` 子进程环境；只设 `-opencl` 在先前部署中会卡在 GL 拓扑发现。没有新增 CUDA 并发限制，没有修改已部署 .venv。

## 预检和测试

- 本地/云端源码和库对应：PASS；所需 Python/输入核对 12 个。
- 云端 GPU：NVIDIA GeForce RTX 5090, GPU-384914d6-f85a-7a2a-6c75-fdf21b3918b3, 595.84, 32607, 15, 32096。
- 自有未结束任务：0；GPU compute apps：空。
- CPU 回归：26 项 PASS，包含源码/库身份、旧协议阻断、步数与 wall 区分、单次授权/启动、未完成进度、几何邻面、HDF5 实读、超时/低磁盘、哈希及导出恢复。
- 新 HDF5 只读适配已在云端现有 Python/系统 HDF5 库读取旧 smoke 原生帧成功；没有安装依赖或运行 GPU。
- 导出测试明确使用旧 smoke 数据，192 点次成员探针；不会当成本轮结果。
- 旧数据页面浏览器检查 9 项 PASS，首末坐标、真实帧数、固定坐标轴、离线请求和 JS 异常检查通过。第一次自动播放冲突保留在测试证据中，修复纯页面导出后通过；没有为此重新仿真。本轮实际结果页面尚未生成，人工验收 PENDING。
- 原冻结包 8726 个条目保护核查 PASS，哈希集合未变化。

## 授权与真实命令

现有云端授权账本仅覆盖部署和累计 180 秒 smoke（实际原账本已用约 40.426 秒）。本轮第二节要求用途变更或额外额度集中确认；因此请求另行批准一次最多 1800 秒的单 GPU 完整任务。该额度按求解进程墙钟计，包含准备和进程内 CPU 工作，不重用本地 RBC 剩余额度，不修改原 smoke 账本。

冻结计划：`/home/lzy/projects/cloud_compute/rbc_full/frozen_plan.json`。
计划 SHA-256（规范化 JSON 内容，不是文件字节哈希）：`b5f63e4a6feccbf33162c5ec5e821f63cf5ca50f53e1d152e72ccdade243627d`。
worker 及输入哈希已纳入计划。完整任务尚无 JOB_ID，未占用新的运行次数。`--new --execute` 在没有独立批准记录时已实际验证会阻断，不把执行选项当作授权。

```bash
cd /home/lzy/projects/cloud_compute
python3 -B cloud_run.py rbc-full --preflight-only
# 下条命令已实现；独立批准记录生成后才允许启动一次：
python3 -B cloud_run.py rbc-full --new --execute
```

任务启动后复用 `status JOB_ID`、`fetch JOB_ID`、`verify JOB_ID`。`rbc-review JOB_ID` 只从已经校验的原始归档导出页面，可用于断线或导出错误恢复，不启动求解。原始结果将存放 cloud_results/JOB_ID，派生 HTML/报告存放相邻的 JOB_ID-review，保证原始 manifest 不被导出修改。云端文件不自动删除，实例不调整，Git 不提交/推送。

本轮当前结论：CLOUD_RBC_RUN_COMPLETE=NOT_RUN；CLOUD_RBC_NUMERICAL_SCREEN=INCONCLUSIVE（未运行）；RESULTS_RETURN_VERIFIED=NOT_RUN（没有本轮结果）；PHYSICAL_MODEL_VALIDATION 保持 NOT_MATCHED/NOT_TESTED；qualified_speedup=null。
