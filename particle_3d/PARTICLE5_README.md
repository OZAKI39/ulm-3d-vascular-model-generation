# Particle-5：球形粒子阻力与法向润滑

P5 在新模块中实现完整 6N 的无惯性阻力平衡，复用并锁定 P0–P4、Frozen FEM 和原 WALL。只在本地验证，不执行 CFD，不进入 P6。

`hydrodynamic_resistance.py` 从冻结 summary 和 solver.xml 交叉核对黏度与 rho×nu；真实回放使用核对后的 0.00345312 Pa·s。`resistance_assembly.py` 组装 CSR 矩阵和 `b = R_self U_free`，静止墙面与相对运动 pair 项右端为零。

球的 self 对角为 `[6πμa]*3 + [8πμa³]*3`。墙面和双球法向块分别为 `6πμa²/h · nnᵀ`、`6πμR_eff²/h · nnᵀ`；双球四块是 `[K,-K;-K,K]`，`R_eff=ai aj/(ai+aj)`。合成扫描使用显式 `VALIDATION_GAPS_ONLY`；真实回放只允许 `h/relevant_radius <= 0.01`，这个数是验证资格，不是生产 cutoff。`h <= geometry_roundoff` 转入原硬接触，没有间隙下限。

`resistance_solver.py` 使用 self 对角代数缩放与 SciPy 稀疏 LU，不求逆矩阵。记录缩放系统的谱条件数、原符号最小特征值和相对后向残差。浮点误差估计 `eps * 6N * condition >= 1/512` 时报告 `RESISTANCE_SYSTEM_ILL_CONDITIONED`；这是数值保护，不是物理阈值。条件诊断目前对本轮小规模系统采用稠密谱分解，因此尚未验证大规模求解性能。

硬接触复用 P4 接触记录和确定性对偶 active set，在完整球形 Jacobian 上通过 Cholesky 三角求解构建 `J R^-1 Jᵀ`，求最小 R 度量修正。球的法向接触角速度修正为零，与 P4 一致。所有墙面和 pair 约束一起处理；乘子标为 `KINEMATIC_CONSTRAINT_MULTIPLIER`，不是物理接触力。

`particle5_motion.py` 每个真实时间子区间重新采样场、计算几何、组装和求解，随后使用 P3/P4 连续墙面与 pair 安全证书。失败请求真实左右时间细分；细分仍失败则停止，不冻结粒子后继续计时。原间隙和几何舍入界保持不变。

`particle5_cases.py` 中的 RBC–MB 系数 1 仅是内部归一测试单位，且只产生验证矩阵。正式接口只接受几何规格并自行计算球形系数，拒绝外部验证矩阵；椭球与胶囊明确报告 `NONSPHERICAL_LUBRICATION_NOT_FROZEN`。没有等效球半径替代。

真实近墙静态回放保留两种原 P4 样本：绝对最小间隙样本正在离墙；用于靠近减速验证的是原轨迹第 125 行、ID 203。真实双球动态烟雾完全沿用 P4 的初值、原 SonoVue 尺寸、dt 和 horizon。原轨迹索引、SHA、资格诊断和所有接受子区间均保存。

在仓库根目录使用已有 WSL Python 环境：

```bash
PY=/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python
$PY -B particle_3d/scripts/run_particle5_validation.py
$PY -B particle_3d/scripts/generate_particle5_report.py
$PY -B particle_3d/scripts/finalize_particle5.py --check
```

环境需要 NumPy、SciPy、VTK、PyVista、Matplotlib、Pillow 和 pytest；完整版本记录见报告 logs。运行 `finalize_particle5.py --bind-commit` 只在自动 PASS、源码已提交且 SHA 未变后绑定被测提交；它不改变科学数据。

审核入口：[中文报告](reports/particle5/PARTICLE5_REVIEW.md) 与 [机器记录](reports/particle5/PARTICLE5_VALIDATION.json)。图源是逐项 CSV/JSON，绘图脚本可独立重生成全部 12 张 PNG；永久测试核对字节 SHA。人工审核状态保持 PENDING_USER_REVIEW。
