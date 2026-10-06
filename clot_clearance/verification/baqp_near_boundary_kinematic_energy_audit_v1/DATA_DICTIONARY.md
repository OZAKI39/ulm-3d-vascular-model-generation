# Data dictionary and reproduction

本目录是固定态诊断。科学入口 `scripts/near_boundary_kinematic_energy_audit.py`，独立诊断模块 `pd_clot/pfpd_candidate/nonlocal_boundary/kinematic_energy_diagnostics.py`；绘图和报告入口分别为 `scripts/plot_near_boundary_kinematic_energy.py` 与 `scripts/report_near_boundary_kinematic_energy.py`，后两者只读取已存科学数组。

运行环境：`/home/lzy/projects/computation_examples/_shared/python_environment/bin/python`；项目目录 `PYTHONPATH=.`，`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`，`python -B`。重跑科学入口会覆盖本轮输出，需先使用新输出目录保存独立复现；本轮没有重跑旧 equilibrium。

## Index and units

节点索引沿用保存的 X（6,859 行），键索引与 BOND_GEOMETRY.pairs（690,020 行）完全对齐。数组均未为画图重排。MOMENT_TENSORS 的 K/invK/M3 单位 m⁵/m⁻⁵/m⁶，G 单位 m，X/depth distances 单位 m，volume/omega0/bond_weight 单位 m³。F 及其偏差无量纲，Ψ 单位 Pa，逐键能量 U 单位 J；方向余弦、condition(K)、比值无量纲。

MOMENT_TENSORS 包含物理 R3、uniform R2、lattice R1 及 continuous 参考；节点 invK 直接保存，uniform invK 也保存。SOURCE_CONTRIBUTIONS 是单个非零输出行的三个输入方向分量，support+volume=total；用于剪切，轴向保持其历史一 ulp L 差异。geometry_class 数轴方向近面数量，不用于力学模板闭合判定。

NODE_DIAGNOSTICS 保存完整实际/预测/残差 3×3 张量；标量 F_error 为 Frobenius 范数。BOND_DIAGNOSTICS 的 predicted_F_error 来自实际节点 D 投影；moment_predicted_F_error 来自独立 M3 预测，是 Ψ 因果重建的输入。projection_residual / moment_chain_residual 由完整张量之差计算再取范数。为节省空间未保存重复的全键 Ftilde 张量。

Psi_actual_bias=Psi_actual−Psi_exact；Psi_predicted_bias=Psi_predicted−Psi_exact，均保留符号。Psi_causal_residual=Psi_actual−Psi_predicted。CSV 的 actual_max/predicted_max/residual_max 取绝对值最大；节点/键 RMS 分别用 V_i / 2*wbar*Vi*Vj 权重。

所有 F 的 scaled 字段除以 0.004；Ψ scaled 字段除以各场 Psi_ref_Pa（解析物理节点最大，节点和键共享）。不使用局部趋零分母。K modifier 的安全条件为分母/0.004>1e-10；projection reduction ratio 的条件为 ||D||/0.004>1e-10。未定义值保留 NaN，并有有效计数/标记。

FBAR_MOMENT_CLOSURE.csv 每节点每场一行；FTILDE_PROJECTION_CLOSURE.csv 按深度聚合，完整逐键量在相应 NPZ。PSI_BIAS_RECONSTRUCTION.csv 按场、节点/键、深度聚合。ENERGY_DEPTH_PROFILES.csv 的每层连续 U2 来自内外两个中心立方体积分差；壳界=半 h 中点±h/4，裁剪到 [0,L/2]。全部壳层加总核对全局能量。最大深度半 h 壳可能没有键，需保留完整连续体积分（本次诊断已核对总和）。

## Provenance and figure policy

PREDECLARED_CONFIG 保存全部固定门槛、资源估算及几何/体积规则。PROTECTED_INPUT_HASHES 保护 425 个旧文件，HISTORICAL_INTEGRITY 和最终 QA 验证不变；旧 manifest 未改。SAVED_TARGET_REPRODUCTION 记录 bitwise 目标态复现和 Ftilde SHA。TASK_REQUEST 是本轮原请求。LITERATURE_PROVENANCE 指向原文缓存，不复制大型 PDF。

六幅图见 results 同名目录。FIGURE_SAMPLING 对每个散点面板记录固定 stride，不随机、不修改值；所有 gate 和 profiles 用全量数据。FIGURE_LAYOUT_QA 检查可见文字边界；最终 QA 检查 600 dpi PNG 与 PDF 无栅格 image 对象。资源统计仅覆盖科学诊断，非完整对话或绘图用时。

## Saved array inventory

下表是实际 NPZ schema；完整 dtype/shape 另存 ARRAY_SCHEMA.json。CSV 列名存 CSV_SCHEMA.json。

| File | Arrays |
| --- | --- |
| AXIAL_BOND_DIAGNOSTICS.npz | actual_F_error, predicted_F_error, moment_predicted_F_error, projection_residual, moment_chain_residual, D_norm, projection_reduction_ratio, Psi_actual, Psi_exact, Psi_predicted, Psi_actual_bias, Psi_predicted_bias, Psi_causal_residual, U0_bond_J, U1_bond_J, Uplus0_bond_J, Uplus1_bond_J |
| AXIAL_NODE_DIAGNOSTICS.npz | actual_delta_F, predicted_delta_F, prediction_residual, actual_F_error, predicted_F_error, F_prediction_residual, Psi_actual, Psi_exact, Psi_predicted, Psi_actual_bias, Psi_predicted_bias, Psi_causal_residual, X, volume, depth_over_delta, geometry_class |
| BOND_GEOMETRY.npz | pairs, abs_n, r_over_delta, midpoint_depth_over_delta, bond_weight, orientation_bin |
| FULL_LATTICE_RECONSTRUCTION.npz | X, shear, axial |
| MOMENT_TENSORS.npz | X, volume, K, invK, M3, K_uniform, invK_uniform, M3_uniform, K_lattice, M3_lattice, K_continuous, M3_continuous, G_physical, G_uniform, G_lattice, family_size, omega0, condition_K, distance_face_m, distance_edge_m, distance_corner_m, geometry_class, K_modifier_continuous, K_modifier_lattice, ratio_valid |
| SHEAR_BOND_DIAGNOSTICS.npz | actual_F_error, predicted_F_error, moment_predicted_F_error, projection_residual, moment_chain_residual, D_norm, projection_reduction_ratio, Psi_actual, Psi_exact, Psi_predicted, Psi_actual_bias, Psi_predicted_bias, Psi_causal_residual, U0_bond_J, U1_bond_J, Uplus0_bond_J, Uplus1_bond_J |
| SHEAR_NODE_DIAGNOSTICS.npz | actual_delta_F, predicted_delta_F, prediction_residual, actual_F_error, predicted_F_error, F_prediction_residual, Psi_actual, Psi_exact, Psi_predicted, Psi_actual_bias, Psi_predicted_bias, Psi_causal_residual, X, volume, depth_over_delta, geometry_class |
| SOURCE_CONTRIBUTIONS.npz | support, volume, total, physical_K_error, continuous_K_error, lattice_K_error |
