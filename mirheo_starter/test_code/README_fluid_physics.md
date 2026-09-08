# 第二阶段：旧纯流体对照工况与 DPD 标定

本模块实际读取已验收 Musubi 工况，运行小周期纯液体测试，保存原始数据与预算，再生成离线核查页。科学结果允许为 **PARTIAL**；程序退出成功不等于材料参数合格。`selected_dpd_parameters.json` 无合格候选时为 `selection: null`，人工验收始终为 `PENDING`。

## 最简单的重新查看

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m test_code.review_fluid_physics --config py_scripts/fluid_physics.yaml --open
```

页面查看只读已完成 GPU 任务，可以重新分析和写入新的不可覆盖报告；绝不启动 CUDA 或标定。原始数据未变时复用同一报告。`--help` 也不启动 CUDA。`--open` 使用第一阶段已验证的 `wslpath` 和 Windows Explorer；如果 Windows 打开失败，终端仍会打印完整 HTML 路径。

## 入口和文件职责

```bash
# 只读来源、单位、资源预检；仅在目标项目写入不可覆盖的准备包
.venv/bin/python -B -m py_scripts.prepare_fluid_physics --config py_scripts/fluid_physics.yaml --preflight-only

# 已授权首轮的剩余任务；同一配置下成功结果命中缓存，不重新计费或执行
.venv/bin/python -B -m py_scripts.calibrate_dpd_fluid --config py_scripts/fluid_physics.yaml --execute

# CPU 自动测试（合成数据不会进入标定结果）
.venv/bin/python -B -m unittest test_code.test_fluid_physics -v
```

正式逻辑在 `py_scripts/fluid_physics/`：`legacy_case.py` 核对旧生成 Lua、最终报告、runtime contract、物理测量截面和刚体变换；`units.py` 提供独立 SI/模拟单位正反换算；`runner.py` 负责锁、预算、进程与输出；`gpu_worker.py` 是实际两 MPI rank 的小盒测试；`analysis.py` 做分块统计和黏度拟合；`boundary_plan.py` 定义每个端口的有限控制区域及原生接口缺口；`reporting.py` 导出可复核结果。

数据位置：

- `data/fluid_physics/prepared_<指纹>/`：来源与单位冻结包，保留历史版本。
- `data/fluid_physics/result_<指纹>/`：正式结果包，包含用户要求的 JSON、配置、中文报告和哈希。
- `runs/fluid_calibration/dpd_round1_20260908/<case_id>/`：实际参数、执行脚本及版本、原始 CSV、Mirheo 日志、控制台、退出与资源记录。
- `runs/fluid_calibration/dpd_round1_20260908/budget_ledger.json`：首轮累计预算；这是调度状态，不能删除来重置预算。
- `test_code/outputs/fluid_physics/result_<指纹>/physics_review.html`：离线真实数据页面。

## 来源与物理含义

权威来源是 `ulm_3D_vascular/outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/` 的最终报告及实际 `fresh_initial_segment.lua`。当前配置的物理常数相同，但 Musubi 二进制版本已不同；本模块保留这个差别，以被验收生成快照为准，没有重跑旧求解器。

ρ=1056 kg/m³；ν=3.27e-6 m²/s；μ=ρν=0.00345312 Pa·s；运动体积黏度=2.18e-6 m²/s（普通 DPD 未独立匹配）。入口 Q=2.7369132390905703e-15 m³/s，质量流率=2.890180380479642e-12 kg/s。出口表压依次为 14.544978101274268、132.20454922317552、−13.700626673311461 Pa。

以上边界数值是用户确认的**测试值**。298.15 K 是用户同意的**温度假设**，没有因此调整密度或黏度。本轮对照液体不称为真实小鼠血浆。小鼠 RBC 尺寸、膜参数未知且未生成；2–4 µm 微泡只用于将来分辨率约束。

旧 3,387,510.72 Pa 是 ρ cs² 形式的 LBM **数值偏置**。历史 gauge zero 来自 1D 模型结构末端零表压，未给出测量绝对压力。本轮将零表压对应到实测 DPD 平衡参考，保留全部出口表压差。不得把未知参考擅自设为大气压，不得把负表压改成零。

第一阶段血管保持解剖坐标。旧测量平面位于刚体旋转后的 CFD 坐标；逐面验证旋转壁面后，仅逆变换测量中心、基向量和法向，保留有限轮廓及端口身份。原始物理切口、实际远端端盖、测量截面分别记录。已有约 11.7–15.7 µm 延伸段直接继承，没有再次添加。

## 单位与材料测量

首轮空间设计 L0=0.5 µm、rc*=1、n*=8、m*=1、kBT*=1。密度固定 M0=1.65e-17 kg；物理热能固定 E0=kB×298.15 K；因而 t0≈3.1656e-5 s。目标 ν*≈414；实测是否达到它是独立问题。所有驱动力和时间步比较使用这一套单位，不能重新调 t0 让结果自动命中目标。

平均粒子间距约 0.25 µm；2–4 µm 对象跨越 8–16 个间距和 4–8 个 rc。这不保证狭缝、黏附或 RBC 分辨率。旧 0.20 µm 是几何/LBM 参考，未来 SDF 采样间距仍未定。

Mirheo DPD 的随机力幅值内部已经使用 `sqrt(2 gamma kBT / dt)`，权函数 `wr=(1-r/rc)^power`，耗散用 `wr²`。不叠加另一套随机力或重复时间步因子。基线 a=25、gamma=10；高耗散量级探针 gamma=3700；它们都是候选设计，非生产默认。

周期双向 Poiseuille 的原生输入是**每粒子力**。方向 x，变号沿 y 的 Ly/2：下半区 −F，上半区 +F。体力密度 nF，正确剖面下半区 `u=-F*y*(Ly/2-y)/(2*m*nu)`，上半区相反。拟合包含非单位 m、常量速度偏置和有限分箱的二次曲线精确平均修正。这是无壁材料测量，没有验证实际管壁无滑移。

## 压力、温度与置信度

原生 virial 输出是 `sum(trace(stress))/3`，配对因子 0.5 已包含，没有动能项，也没有除体积。分析补入 COM 扣除后的动能压力，再除真实盒体积并转换 Pa。`Stats.vx/vy/vz` 实际是每粒子平均动量，除以 m 才得到均速。

两个插件在同一个 `afterIntegration` 回调采样，但 Stats 在之后发送时写当前时间，而 virial 保存采样时的时间。因此原始 Stats CSV 时间比 virial 大一个 dt；分析按该源码事实核对对应关系，并保留原始文件。源码和实测偏移均有证据。

平衡温度扣除质心运动并使用 3(N−1) 自由度；流动温度扣除瞬时 y 分箱均速，使用 3(N−非空箱数) 自由度。有限箱宽留下的剪切速度差会使温度偏高，生产前仍需分箱敏感性测试。

先按 Ly²/(4π²ν_prior) 估计黏性松弛，再在允许排除范围内用实际数据判断平衡段。时间块至少覆盖 2 个松弛时间和 5 个积分相关时间；至少 8 块。黏度 CI 从斜率分块 CI 反演，斜率区间跨零时不给虚假的有限黏度误差条。统计不足、持续漂移、CI 过宽分别保留，不用最后两点接近宣布收敛。

PROPOSED 门槛在 GPU 之前冻结：温度 2%；目标黏度、拟合相对 RMS、黏度 CI 相对半宽、时间步/驱动力差异 10%；Mach≤0.15。这些是首轮筛选门槛，不是用户已验收的最终精度。旧短/长物理窗口为 0.244140673 / 0.488281346 ms，窗口够长仅表示有该时间跨度，不表示周期盒执行了真实血管验收。

旧压力 residual 是两个时间检查点的变化，分母 `max(abs(p_current),1 Pa)`，不是与旧压力的 0.5% 误差。旧出口分流比例分母是 `sum(Qout)`，密度流量一致性分母是 `abs(Q_rho_u/rho0)`，均已与接受报告的实际采样重新核对。显著平均回流是 `mean Qout<0` 且其幅度超过 `0.05*abs(mean Qin)`，短长两窗均检查；不限制每颗粒子的瞬时反向速度。

## 边界下一阶段实现

`boundary_plan.json` 逐端口给出实际端盖、原始切口、三个有限测量截面、现有延伸段和局部控制层。计划在现有端盖附近、沿法向 `[−rc,0]` 的有限局部腔体内控制；不是已生成的新缓冲几何。

入口可以复用 `createVelocityInlet` 的静态供液原型。Q/A 只给平均法向速度；当前绑定没有已核实的运行时调供液 setter。下一阶段需要任意有限截面通量测量、可调供液/储液反馈、质量库存记录。

出口可以复用 `createDensityControl` 的平滑密度反馈力与 `createDensityOutlet` 的局部超额密度移除组件，但后者只会移除，不是完整压力储液边界。需要双向局部储液交换、EOS 有效范围检查、局部总应力测量与反馈。无限 `PlaneOutlet` 不用于多分支压力控制。`VelocityInlet` 使用绝对局部通量，必须验证速度向内、其他表面速度为零；源码还有每三角面每步新增小于 5 个粒子的限制。没有虚构缺失 API，没有修改 C++/CUDA。

## 预算与安全停止

首轮 campaign 固定为 `dpd_round1_20260908`；每任务≤600 s，总额≤3600 s，一个 GPU 任务。首次 2000 步实际预检决定后续步数。至少保留 1536 MiB 显存及 1024 MiB 主机可用内存；测试显存估计 512 MiB，不关闭其他程序。

调度器使用 monotonic 计时，先写合作停止文件，两 rank 通过同一分块决策正常结束；必要时仅对自己启动的进程组发送 SIGTERM/SIGKILL，检查无遗留工作进程。正常结束、失败、重试、MPI/CUDA 初始化、输出与退出都计入预算。监测的显存是**设备整体采样峰值**，含其他应用，不能称为真实任务峰值。

成功缓存核对配置、单位、参数、库、精度、脚本及原始数据哈希。不同参数或失败结果不自动重跑；`--retry-failed` 是明确重试操作，会建立新 attempt 并继续累计旧账本。被强制杀掉的调度器留下全额预约，后续拒绝启动，需先核实其记录进程并显式恢复，不能删除账本。追加 campaign 或预算必须重新明确授权，不会由查看页面或重复执行自动发生。

源码、共享库、现有补丁、两个只读来源项目和第一阶段数据均保留。源文件完整性核查与实际浏览器检查写在交付目录。原生 Mirheo 示例的 MIT 许可证随每个 GPU 任务快照保存。
