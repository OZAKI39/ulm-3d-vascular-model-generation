# Shear-Lift Magnitude Audit

本目录完全独立于正式轨迹求解器。结果入口为 `SHEAR_LIFT_MAGNITUDE_REVIEW.md` 和 `OPEN_RESULTS.html`。

## 内容

- `code/`：公式、只读全状态审计、数据清点、统计和可视化；`shared_rotation.py` 逐字复制既有固定 Z 轴旋转显示代码，仅复用显示类。
- `tests/`：公式 convention、标度、零值、有效性、只读、几何和图像输出检查。
- `data/trajectory_inventory.json`：每条原数组/元数据的来源、SHA256、去重关系、终止情况。
- `data/readonly_baseline.json`：正式源码、冻结输入、原轨迹的审计前哈希。
- `data/remote_results/`：远程全量逐状态 NPZ、CSV 统计、JSON、两版 Eulerian 数据和代表轨迹 CSV。
- `figures/`：英语学术图，统计图同时提供矢量 PDF。
- `animations/`：相机绕固定 Z 轴旋转；使用已有保存位置，动画不代表新轨迹积分。
- `literature/`：来源表、阅读/下载状态和可取得的原始资料。
- `logs/`：远程运行和测试日志，包括首次打包遗漏依赖而被完整性检查拦截的记录。

## 复现

使用项目现有 `.venv/bin/python`，设置 `PYTHONDONTWRITEBYTECODE=1`、`OPENBLAS_NUM_THREADS=1`、`OMP_NUM_THREADS=1`。`prepare.py` 清点并打包，正式计算入口为 `run_audit.py data/remote_config.json`；配置中的路径是原服务器隔离目录，迁移机器时仅调整配置路径。正式计算后依次运行 `final_statistics.py data/remote_config.json`、`post_verify_remote.py data/remote_config.json` 生成补充统计和输出哈希。`refine_ratio_guards.py` 仅用于复现本轮初始审计的数值保护修订，当前 `run_audit.py` 已含相同预算，不必重复修订。全量后处理在服务器执行；本地 `local_config.json` 只读取 300 状态小样本。

本地绘图命令：`python code/plot_audit.py`、`python code/render_audit.py`。检查：`python code/validate_inputs.py` 和 `python -m pytest -q tests`。依赖沿用现有 NumPy、SciPy、PyVista/VTK、Matplotlib、Pillow、imageio-ffmpeg、pytest；不安装新 GPU 计算栈。

## 逐状态字段

所有量使用 SI。`endpoint_*` 是保存的当前位置诊断；`evaluation_*`、`local_u_m_s` 和梯度来自产生该保存速度的前一接受位置。初始行为自由流诊断，`initial_diagnostic=true`；dt=0。`E_s_inv`、`W_s_inv` 和 `curl_s_inv` 保留完整张量/向量。

`CANDIDATE_SAFFMAN_MAGNITUDE` 是三维局部应变标量下的刚性球量级代理。`candidate_direction_only` 是归一化的 slip×curl 示意，并非已验证的任意三维升力向量。没有把候选力写回正式动力学。

比值状态码：0=可计算；1=`BOTH_FORCE_MAGNITUDES_NEAR_ZERO`；2=`DENOMINATOR_NEAR_ZERO_NOT_EVALUABLE`；3=`NO_ACTIVE_LUBRICATION`（仅润滑比）。无定义比值在 NPZ/CSV 中为 NaN；JSON 统计忽略 NaN 并记录有效 count，空组用 null，禁止用分母下限制造比值。

数值零值预算从现有线性求解器的 `512×6×machine_epsilon×condition×max(|u|,|v|)` 速度尺度出发。condition=1+ζ_n/(6πμa) 是该独立球/单法向壁面块的 self-scaled 矩阵条件数，用于把后向误差预算传播为保守的前向速度预算；它不是对多接触约束灵敏度的完整认证。drag 分母乘 self resistance；润滑分母取 self 与 lubrication resistance 放大预算的较大值；候选力分子使用自身前因子传播该速度预算。这样 handoff 的近零法向速度不会因大润滑系数而被误认为可分辨的受力。该预算不表示物理忽略阈值。绝对力仍原样保存。`descriptive_screen_*` 仅是渐近分离参数的灵敏度筛查，不能称为已验证有效区域。

`region_far_field` 仅指现有 V1 近场修正关闭；不等于相对于 Saffman 扰动长度远离壁面。分叉、分叉前、入口/出口邻域和主干中心均为明示的描述性 ROI，彼此允许重叠，不能相加成分区概率。全部统计提供状态等权和按 dt 加权两种口径；后者仍是独立轨迹时间加权，并非浓度或生理分布。

Eulerian 文件是每个原四面体的重心；P1 梯度在单元内常数，统计另提供体积权重。表面图显示外表面相邻单元的局部剪切标量，不是 WSS。新流场中只有剪切、Re_G 和力/滑移前因子，未设定实际微泡滑移。敏感性表明示 `HYPOTHETICAL_ONLY`。

最终检查和报告：结果下载后核对 `remote_results/OUTPUT_SHA256.json` 的每项哈希，记录在 `data/delivery_verification.json`；保存测试结果与逐图检查记录，然后运行 `python code/build_report.py`。`supplementary_statistics.json` 包含初始行/接受区间的独立统计、全部区域的两种比值状态码和纯数值零值预算的 0.1/1/10 倍敏感性。
