# 固定参数的无驱动 SDPD 稳定性实验

本轮只延长一个周期纯液体任务。诊断来源由
`py_scripts/sdpd_equilibration.yaml` 指定的 delivery_record 和 package manifest
哈希确认，不按文件修改时间选择，不拉取远程代码覆盖本地。

当前没有新增 GPU 授权。原池累计用量为 3492.9193288880488 s，原余额
107.08067111195123 s。完整任务预计 582.1243803313985 s，建议单任务
上限 600 s，需要追加 492.91932888804877 s；建议明确批准 **493 s**，
只用于本计划。实际以预检读取的全部账本为准。

预算不足时状态为 `READY_TO_RUN_AWAITING_BUDGET`，执行入口返回 3，
不初始化 CUDA，不创建缩短的实验，不改变原预算。

## 复用与新增

- 复用 `fluid_physics.analysis.block_statistics` 的完整块均值/CI 一致性逻辑，
  `native_pressure` 的压力定义、现有 CUDA 通道布局和相位核查。
- `gpu_worker.py` 把原分块控制循环抽成 `continuous_chunks`，只加可选观察回调。
  一个协调器、一次初始化、两 rank 匹配 run 序列；每块没有重置时间、速度或原生随机状态。
- `equilibration_worker.py` 加 t=0 哈希与 IDs、稀疏同相位快照、应力及局部密度摘要。
  原生 `stresses` 的顺序为 `[xx,xy,xz,yy,yz,zz]`，不是前三列都为对角元。
  它的 persistence 为 None，故使用 `saved_stresses` 的积分前持久副本，
  不依赖 run 末尾 cell-list 重排后临时应力仍然有效。
- 复用原 runner 的锁、单任务进程组、monotonic 看门狗、资源预留和退出清理；
  `strict_allocation=True` 防止余额变化时自动缩短任务。
- 新预算扩展只读取显式追加的授权记录；原 `limit_s=3600` 和全部旧账本保留。
  追加额先用于指定 campaign/task，其余任务不能消费该范围的未用追加额。
- 页面复用现有 Plotly、表格和 Windows Chrome/CDP 检查流程。

## 冻结设计

从真实 `sdpd_equilibrium/actual_parameters.json` 复制精确参数：
`mu*=3312.4601269222444`，`nu_prior*=414.0575158652805`，
`dt*=1e-6`，`[8,8,8]` 周期盒，4096 粒子，`n*=8`、`m*=1`、
`kBT*=1`、`rc*=1`、WendlandC2 Density/SDPD 核、Linear EOS、
`sound_speed*=120`、`rho_0*=0`。单位原值和初始速度种子 20260908 不变。

原生 `VelocityVerlet` 无外加驱动力。SDPD 自带的随机、耗散和压力作用保留。
高斯速度只抽取一次、COM 只在初始化时扣除一次，不额外恒温或持续缩放速度。

连续 400000 步推进到 `t*=0.40`，物理时间约 12.662309354 µs。
每 200 步采样；正式统计采用固定右端采样区间 `(0.30,0.40]` 的 500 个样本。
初态及完整启动曲线保留。顶层步数、desired_time 和当前 task 字段相互一致；
原任务模板仅存于明确的 `historical_template`。

稀疏快照的 post-step 标记为 200、100000、200000、300000、350000、400000。
每帧积分前位置、速度、力、密度对应 `t*=(step−1)dt`；积分后状态对应 `step*dt`，
两相位的 packed ID 顺序及完整粒子 ID 集合均核查。普通诊断快照不称为 checkpoint。

压力为配对应力迹/(3V) 加一次 COM 扣除后的动能动量通量。
应力来自积分前、速度来自积分后，明确保留 dt 的相位差；不额外再加一次 EOS 压力。
仅按冻结采样间隔下载通道，不每步复制所有粒子到 CPU；重型快照 CPU 复核在 GPU 退出后进行。
新增应力 saver 每步在设备内复制一次，额外成本尚未 GPU 实测，保留于预算余量假设中。
已先用 CPU 假设备复现非持久应力失效，再最小修复观测代码；未启动过本次 GPU 实验。

## 判定边界

先查看固定两半、四分窗、局部核密度分布、动量及同相位快照，再估计 ACF。
明显持续漂移时不把 ACF 标为可靠稳态相关时间。旧 82.44 样本仅为历史规划参考。
只有固定四段呈同向变化、变化超过冻结限值，且两半变化超过段内 RMS 的 3 倍时，
才据趋势标记持续过渡；噪声较大而证据不足时保留“不确定”。该筛选不把样本视为独立。

块长至少 `max(5×本次全窗/两半 ACF, 2×盒黏性衰减先验时间)`，物理下限约
0.007830501204*（向上取到 40 个采样间隔）。全窗至少 8 块、固定两半各至少 4 块。
正式均值和 CI 使用相同完整块；全窗均值、尾部与两半各自使用的样本单独报告。
块数不足时有效 CI 为 null，不为凑块数缩短块长。

温度两半完整块均值差，加两半均值 CI 半宽之和，除以目标热能，要求不超过 0.5%。
该误差带是保守构造，不声称多个指标的联合 95% 覆盖。
稳定且采样充分后才检查 `abs(mean/target−1)+CI_halfwidth/target<=0.02`。
百分比以 298.15 K/目标热能计算。

压力漂移以 `n*kBT=8*≈0.263449920 Pa` 为绝对尺度，禁止用巨大背景压力作分母。
核密度均值/分位数以原 n*、方差以 n*²、COM 速度以热速度归一化。
最近邻分位数和配位数还用三个固定末段快照检查。
精确阈值见 YAML 和 run_plan.json；都是本次事前设计，不是普适定律。

输出状态包括：`TRANSIENT_PERSISTS`、`STATIONARITY_OR_SAMPLING_INCONCLUSIVE`、
`STATIONARY_BUT_THERMALLY_BIASED`、`EQUILIBRIUM_SCREEN_PASS`、
`MEASUREMENT_OR_ANALYSIS_ERROR`、`NUMERICAL_OR_RUNTIME_FAILURE`。
合作停止即使退出码为 0 也标 `PARTIAL`，不作为完成的正式窗口。
不自动换参数、重试、移动窗口、续跑、Poiseuille 拟合或 EOS 扫描。
所有情况下生产参数 `selection=null`；温度初筛不代表全部物性通过。

## 入口

```bash
cd /home/lzy/projects/mirheo_starter

# CPU 准备、核对哈希与全部预算账本；生成不可覆盖的数据包与 HTML。
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration.yaml --preflight-only

# 仅在明确授权已登记、余额足够且 CPU 记录与代码一致后执行。
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration.yaml --execute

# 重新分析/查看均不启动 GPU，也不自动重试失败任务。
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration.yaml --analyze-only
.venv/bin/python -B -m test_code.review_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration.yaml --open

# 真实离线浏览器核查；只有该过程实际成功才写 browser PASS。
.venv/bin/python -B -m test_code.review_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration.yaml --browser-test

# 全部 CPU 回归；不运行旧 GPU 任务。
.venv/bin/python -B -m unittest discover -s test_code -p 'test_*.py' -v
```

没有明确模式时只显示参数错误，`--help` 不导入 Mirheo/CUDA。
重复预检/查看复用经过哈希检查的产物。浏览器自动检查与人工验收分开，人工仍为 PENDING。

## 追加授权的登记

**本交付未创建任何真实追加授权记录。** `extra_authorized_gpu_seconds: 0` 固定不变。
不得通过修改 YAML、写 `approved:true` 或伪造用户文本放行。

只有在用户真实消息明确同意秒数和本任务范围后，才把该原文保存为独立文本证据，
再准备 `USER_CONFIRMATION_JSON`，包含以下字段：

- `authority`: `explicit_user_message`；
- `user_message_file`: 保存真实用户确认原文的绝对路径；
- `additional_seconds`: 用户明确批准的秒数；
- `scope`: 最新 `budget_request.json` 中 `authorization_scope` 原值；
- `budget_request_sha256`: 该预算请求文件的 SHA256。

此登记入口不是批准本身，调用者必须已有真实用户授权。程序不会生成用户同意文本。
记录绑定原文哈希、精确计划、请求哈希及单任务 600 s / 并发 1，历史用量保留。

```bash
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration \
  --config py_scripts/sdpd_equilibration.yaml \
  --register-authorization /absolute/path/to/USER_CONFIRMATION_JSON
```

之后再预检和执行。同一授权不可重复登记；改动证据、任务范围或已批准计划会拒绝执行。
审批登记只追加记录并保存池的前态，不篡改旧成员账本。

## 输出与证据

`data/sdpd_equilibration/<run_id>/` 包含计划、预算请求、独立结果/采样/稳定性 JSON、
raw_statistics.csv、provenance.json、中文报告及完整文件 manifest。
等待授权时 CSV 只有表头，实际步数、温度、压力和新统计均为 null；
`actual_parameters.json` 仅在真正运行后发布。

实际原始运行写入 `runs/sdpd_equilibration/<campaign>/<case>/`，含冻结执行代码、
真实初态哈希、同相位状态、各采样摘要、日志、monotonic 用量与进程组退出记录。
查看页面位于 `test_code/outputs/sdpd_equilibration/<run_id>/<viewer_hash>/`。

CPU 假设备/合成数据只保存在自动清理的测试临时目录，不能作为液体测量证据。
GPU 尚未执行，原生新观察路径目前由源码契约和 CPU 生命周期测试覆盖，物理结果仍未知。
