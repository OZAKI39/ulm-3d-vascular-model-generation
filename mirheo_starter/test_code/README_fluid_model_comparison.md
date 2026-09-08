# 原生 DPD / SDPD 纯液体对照

本模块复用 `fluid_physics` 的单位、来源、文件保护、进程调度、黏度拟合和分块统计。
仅使用现有 Mirheo 库；主配置为 `py_scripts/fluid_model_comparison.yaml`。
源物理工况、精确单位和旧预算由指定结果包、配置与 SHA-256 验证。

## 实际入口

```bash
cd /home/lzy/projects/mirheo_starter

# CPU 来源、单位、接口源码和预算核查；生成不可覆盖的准备包
.venv/bin/python -B -m py_scripts.compare_dpd_sdpd \
  --config py_scripts/fluid_model_comparison.yaml --preflight-only

# 执行已冻结计划，只动用原 3600 秒授权的剩余额度
.venv/bin/python -B -m py_scripts.compare_dpd_sdpd \
  --config py_scripts/fluid_model_comparison.yaml --execute

# 只分析已经完成的记录并打开离线页面，不启动 GPU
.venv/bin/python -B -m test_code.review_fluid_model_comparison \
  --config py_scripts/fluid_model_comparison.yaml --open
```

`--help` 在导入 Mirheo 之前退出。比较入口必须显式指定 `--preflight-only` 或
`--execute`。`--only` 只能选择已冻结任务。成功缓存核对完整原始输出哈希；失败或
输入不同的任务需要显式 `--retry-failed`，新 attempt 继续计费。

## 预算与文件保护

旧来源账本：`runs/fluid_calibration/dpd_round1_20260908/budget_ledger.json`。
共享引用：`runs/fluid_calibration/shared_budget_pool.json`。
新任务目录：`runs/fluid_model_comparison/dpd_sdpd_remaining_20260908/`。
每个 campaign 保留自己的账本，共享引用记录成员、原授权上限、初始来源哈希和
同一 GPU 锁。旧调度命令也会自动识别其已注册的共享预算成员身份。

累计值包括旧任务、本轮完成/失败任务，以及所有未结束任务的全额预约。注册后缺失
任何必需账本都会阻断启动；更换 campaign 名称不会取得新的 3600 秒。
每任务上限 600 秒、GPU 并发 1。合作停止及必要的信号仅作用于调度器自己的进程组。
时间包括 MPI/CUDA 初始化、计算、采样、退出和进程清理。

原始结果不覆盖。新准备包与分析结果保存在 `data/fluid_model_comparison/`，
页面及浏览器检查保存在 `test_code/outputs/fluid_model_comparison/`。
输入、代码或原始运行改变时生成新结果版本；已有包回读校验通过后才能复用。
查看页面会先验证历史与新任务文件哈希，同一结果/源码指纹命中现有结果。

## 本轮模型与观测

主 SDPD 候选使用精确冻结的 `mu_star=3312.4601269222444`、`n_star=8`、
`m_star=1`、`kBT_star=1`、`rc_star=1`。目标 `nu_star=414.0575158652805`。
Density 与 SDPD 分别注册，使用相同 Wendland C2 核和支持半径。
SDPD 粒子对不同时应用 DPD 的整套力或恒温器。

Linear EOS 的源码定义为 `p_i = cs_star**2 * (m_star*d_i - rho_0_star)`，
其中 `d_i=sum_j W(r_ij)` 包含自贡献。本轮 `cs_star=120` 根据目标压差与
不超过 5% 的密度变化初筛范围设计；`rho_0_star=0` 给出正背景配对压力。
背景压力对离散结构、初始化和机械压力的影响由实际测试保留。
它没有被宣称为已知热力学绝对压力。

全盒 `N/V`、`N*m/V`、局部核数密度和 `m*d_i` 分开记录。
原生 ParticleChannelSaver 在 `beforeIntegration` 保存密度和同相位位置；
使用现有 CUDA runtime 的同步内存复制读取原生通道，没有增加第三方 GPU Python 包。
少量真实快照由 CPU 独立重算 Wendland C2 密度，保留自贡献、周期距离、误差和 CPU 耗时。

压力使用原生 SDPD/DPD 的配对维里应力，加一次动能动量通量，除以真实体积；
独立机械压力与输入 EOS、局部密度代入 EOS 的预测分开。
EOS 只使用平衡任务，流动时全局压力迹中包含的宏观流动贡献不冒充平衡 EOS。
Stats/virial 的时间差逐运行核对；六位有效数字输出的量化误差按两个时间标签分别计算。
短探针也保留实际压力字段，不以没有足够平衡样本为由隐藏原始记录。

流动态温度使用实际瞬时局部均速扣除，扣除所拟合的自由度；比较半宽分箱与原分箱。
不减去目标解析速度场来凑温度。剖面仍采用相同 0.25 模拟长度分箱，
使用原有每粒子力/粒子质量/运动黏度公式和二次曲线分箱平均修正。
目标解析线、实测均值和拟合实测黏度的线分别输出。

证据按 `method + candidate_id + test_kind + comparison_group` 选择。
不同候选不能借用 EOS、温度或减半检查。新分析还修正了原 EOS 的固定密度 8、
候选名分支、将宽误差条误当成压力覆盖区间，以及只检验流动而遗漏温度漂移的问题。
旧报告和 CSV 保留，新算法只产生新版本结果。

继承上一轮 PROPOSED 判据：温度 2%、目标黏度/拟合残差/CI/敏感性 10%、
至少 8 个统计块等。新增密度、核密度和温度分箱检查在执行前写入配置。
统计块覆盖相关时间和黏性松弛下限；初始化、排除段和统计段各自记录。
旧血管短/长窗口从冻结文件读取，仅报告时长覆盖，不替代局部材料统计判据。

## 成本比较与适用范围

短成本控制保持粒子数、区域、精度、实际时间步、采样间隔和主要输出相同。
DPD 加被动 Density 的控制用于估计密度路径额外开销；这包括相关 halo、
通道保存及调度，不称为单个 CUDA kernel 的计时。
同步计算计时在 `u.run` 后明确等待 CUDA 完成，另外记录 CPU 传输和输出时间。
数千粒子下的 MPI/调度开销也包含在成本中，不能外推到整个血管。

正式材料测试不强制相同 dt；报告每步成本、推进单位物理时间成本、实际演化时长、
有效统计时长和误差。未取得合格材料时，“达到合格解的时间”为 null，
不按失败的低黏度材料对目标问题作速度排名。

## 自动测试与浏览器

```bash
.venv/bin/python -B -m unittest \
  test_code.test_fluid_model_comparison test_code.test_fluid_physics \
  test_code.test_vessel_geometry test_code.test_vessel_geometry_migration -v

# 需本机已有 Windows Node/Chrome；另存实际交互结果和截图
.venv/bin/python -B -m test_code.review_fluid_model_comparison \
  --config py_scripts/fluid_model_comparison.yaml --browser-test
```

CPU 合成样本仅测试代码，不进入实验对照表。页面为离线 Plotly，不启动 Web 服务。
浏览器记录绑定实际 HTML 哈希，只有真正通过渲染、剖面切换、鼠标缩放、图例及恢复
检查才写 PASS；人工验收始终为 PENDING。

本轮没有进行真实血管、SDF、壁面、RBC、微泡、黏附或开放出口计算。
下一阶段仍需检查 SDPD 近壁核截断、冻结壁面、膜内外液体/膜节点耦合、
局部储液区及压力观测。纯液体结果不自动验证这些功能。
