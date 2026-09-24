# P9-A.2 复现

正式入口：`METHOD_C_SIZE_FIRST_FLUX_CONDITIONAL_V1`。新配置：`particle_3d/contracts/PARTICLE9A2_PRODUCTION_V1.json`。

## 输入与环境

- 当前 2.0 mm/s 冻结 FEM SHA：`129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d`。
- 未修改的 SonoVue 根目录：`/home/lzy/projects/sonovue_size_distribution_v0`。
- 它的 bin 内均匀定义来自原 contract；新的 4 µm contract 只做条件截断。
- 服务器：`/workspace/particle9a2_inlet_20260924T081657Z`。
- Python：`/root/particle8_2_runs/env/bin/python`；实际版本 Python 3.13.15 / NumPy 2.5.3。
- WSL 使用项目 `.venv/bin/python`。跨 NumPy 版本的逐字节一致性没有预先保证；同一服务器环境已做全量、逆序、单 worker 重放。

## 必要顺序

以下命令从项目根目录执行。重新计算必须复制代码与冻结流场到一个新的 `/workspace/particle9a2_inlet_<UTC_TIMESTAMP>`，不得覆盖此次输出。准入数据不能复用 Method B。允许复用只读 FEM 和几何。

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=particle_3d/src
P9A2_PYTHON=/root/particle8_2_runs/env/bin/python
"$P9A2_PYTHON" -m pytest particle_3d/tests/particle9a2_inlet -q
"$P9A2_PYTHON" particle_3d/scripts/run_particle9a2.py audit --workers 6
```

生成前必须读取并保留 `data/protected_before.json`；历史保护测试按其记录比对旧 Method B 源文件。`data/local_gate.json` 必须来自真实测试成功记录，不可手填绕过测试。2000 个事件及其 point tracer 完成后，才可继续：

```bash
"$P9A2_PYTHON" particle_3d/scripts/run_particle9a2.py smoke --workers 6
"$P9A2_PYTHON" particle_3d/scripts/audit_particle9a2.py smoke30 --workers 6
"$P9A2_PYTHON" particle_3d/scripts/check_particle9a2_flux_reference.py
```

检查 `data/smoke30_gate.json` 的 `no_new_systematic_failure`，以及独立 flux reference 的 `all_position_agreement`。只在真实 Gate 通过后运行：

```bash
"$P9A2_PYTHON" particle_3d/scripts/run_particle9a2.py formal --workers 6
"$P9A2_PYTHON" particle_3d/scripts/audit_particle9a2.py formal500 --workers 6
"$P9A2_PYTHON" particle_3d/scripts/check_particle9a2_stationary_equilibrium.py
```

`formal` 会重新生成 seed=2026092494 的 500 个事件，先保存同中心 point 去向，再进行 P9-A.1 积分；已有正式 ledger 会触发拒绝覆盖。audit/smoke seeds 分别是 2026092492/2026092493。单个 birth 由 seed、永久 ID、角色、attempt 确定，与 worker 排序无关。

## 图表与报告

将服务器输出同步回 WSL 后执行：

```bash
export PYTHONPATH=particle_3d/src PYTHONDONTWRITEBYTECODE=1
.venv/bin/python particle_3d/scripts/render_particle9a2_inlet.py
.venv/bin/python particle_3d/scripts/finalize_particle9a2_inlet.py
```

`reference/` 保存旧 B ledger、旧审核 CSV 与流量的只读副本。最终图包含真实 inlet 三角形，不使用理想化管道代替真实入口。全部图提供 PNG/PDF；图源 CSV/JSON 可以独立作图。

## 复核重点

1. `inlet_size_capacity.py` 用三角形细分和真实 WALL 距离上界证明哪些位置可排除；无经验几何余量。
2. `injection_method_c.py` 的直径流和位置流完全分开，尺寸冻结后的 API 没有尺寸抽样调用。
3. 全局尺寸拒绝仅在 radius 大于几何容量上界时发生；括区内继续求精而非猜测。
4. 加速保留区域是可行域的超集，仍由原 flux sampler 提案、原准入 checker 判定。
5. `audit_particle9a2.py` 不推进轨迹，对保存的每个运动段重做连续 handoff 证书；首先保守检查终态接触与顺流几何可行性，再由 `check_particle9a2_stationary_equilibrium.py` 独立求解原机械平衡问题。两种 stationary 证据分别保存，不能把几何可行性和水动力平衡混为一谈。
6. 时间限制、solver failure、无法解释的停止均不伪装成出口；任何新增数值问题导致 FAIL。
7. `delivery_manifest.json` 是所有交付报告、记录、代码和新输出的哈希清单；清单本身及传输确认文件为自引用例外。

最终 stationary Gate 使用原 6D SPD 阻力问题的独立 KKT 平衡证书。几何上存在 FEM 顺流方向，不能单独判定当前仿射水动力模型的静止是数值失败。两阶段审核均保留原始记录；第二阶段不修改或重新积分轨迹。
