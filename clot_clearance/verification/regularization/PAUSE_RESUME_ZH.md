# 暂停与恢复记录

用户于本轮开发中明确要求在下一个可停止点暂停，以便电脑休眠。**已暂停，不得自行恢复或启动队列；等待用户新的恢复指令。**

暂停时本任务没有活动计算、测试、渲染或自动续跑进程。627 个受保护旧文件重新核对全部一致。此结论仅针对当前 clot_clearance 开发任务，不代表其他项目的独立任务状态。

## 当前任务和授权

- 继续现有 `/home/lzy/projects/clot_clearance`，不能重建替代工程，不能改写旧文件。
- 最新完整提示词已保存：`provenance/regularization_20260930/USER_PROMPT.txt`（原附件 `cafe2a16-c223-4279-9291-804431804cd1/Pasted text.txt`）。
- 用户补充：局部载荷中心应在 clot 左侧来流方向；完成新模型开发后，按现有黑底、Arial 常规字重、单主图、损伤配色、固定视角及真实位移比例生成新动画，**不显示箭头**，不恢复 N 顶部文字或右侧 Force 图。
- 用户已回复“使用演示值”：本轮选择 `Gc_demo = 0.01 J/m²`，明确未实验标定；需要敏感性比较，不再重复询问授权。
- 新中心配置为 `[-0.000825, 0, -0.000725] m`，即 clot 左面 x=-0.625 mm 上游 0.2 mm，y/z 与 clot 几何中点齐平。这是新计算输入，不能只移动旧结果图上的图标。
- 原模型材料、物理 clot 尺寸、管流参数、0.01 m/s 声流式速度标度和 25,000 代表周期保留；D_break=1。不以产生碎裂为目标调参。

## 已完成及证据

项目没有 Git 仓库，源码身份由文件快照和哈希确定。未发现适用 AGENTS.md。

1. Phase 0：`provenance/regularization_20260930/` 保存旧源码/配置快照、627 个旧源码/配置/结果/可视化文件哈希。旧 26 项测试通过，见 `verification/regularization/BASELINE_TESTS.log`。
2. Phase 1：实际 NOSB 弹性/稳定化能量、动能、阻尼耗散、外功与损伤能量释放账本。固定坐标损伤更新如使储能增加则报错，不能隐藏。3 项测试通过。
3. Phase 2：`pd_fracture_energy_coupon`，固定 0.5 mm 立方体、两半块位移控制循环张开、对称弱平面。跨平面键仅由损伤律失效，没有切键命令。粗档自动标定完成：目标 0.01，实测 **0.00999993791219862 J/m²**，相对误差 **6.208780138128311e-6**；标定临界能量密度 **0.012520992043877192 Pa**。独立力积分外功相对残差 **2.1506829457554293e-15**。`calibration_coarse/CALIBRATION.json` 保存全曲线。Phase 1+2 共 5 项测试通过。中、细档标定尚未做。
4. Phase 3：新 `energy_regularized_cyclic` 损伤类，校验标定的 spacing/horizon/material/C_E/m_E/Gc 身份；D=1 才断裂。旧律包装类遇 D_break<1 发警告并标识 forced_failure_verification；冻结的旧入口仍保留原样。2 项测试通过。
5. Phase 4：已实现 resolved/under-resolved/singleton 诊断分类、P_singleton/P_lowrank、支撑秩分布、局域性、按首次穿越类别单计清除量。无脱落时分辨率标为 NOT_ASSESSED_NO_DETACHMENT，不能当作已证明充分分辨。1 项测试通过。
6. Phase 5：`StreamingFlowProvider`/`ManufacturedStreamingProvider` 及原表面加载几何的适配类。同一 u/p/grad 产生牵引和输运。高斯向量势旋度给出局部涡流，解析梯度/散度/时间导数/衰减/应力符号 3 项测试通过。
7. Phase 6：**刚写出 `pd_clot/regularization/runner.py`，尚未编译/烟雾测试，也未运行任何新的正式 regularized clot 算例。** 不能宣称新模型已跑通或动画已生成。

新功能测试共 11 项（分阶段执行，不含重复运行的 Phase 1 测试）。每阶段日志为 `PHASE1_TESTS.log` 至 `PHASE5_TESTS.log`。

## 新文件

源码：`pd_clot/regularization/{__init__,energy,coupon,damage,diagnostics,flow,runner}.py`。

配置：`configs/streaming_regularized_demo_coarse.json`。

测试：`tests/test_regularization_{energy,coupon,damage,diagnostics,flow}.py`。

新源码身份、旧文件保护核对、暂停时间：`verification/regularization/PAUSED_STATE.json`。旧结果目录仍为：

- `results/straight_pipe_damage_demo` → `results/straight_pipe_demo`。
- `results/streaming_fragmentation_demo` → `runs/fragmentation_pilot_003`。
- 最新已有动画仍是 `visualization/streaming_fragmentation_demo/streaming_traction_v4/`，它是暂停前上一轮成果，仍含牵引箭头，**不是本轮要求的新无箭头动画**。

## 恢复后的顺序

1. 用户明确恢复后，读取本文件及完整提示词，核对 PAUSED_STATE 与旧文件哈希。
2. 先审查并烟雾测试新 runner：能量积分、首次损伤前后账本、无 NaN、合法 Jacobian、连通图及分类守恒。它尚未测试；不要跳过验证直接开大批。
3. 运行 Phase 6 粗档候选，保留所有失败诊断，不覆盖旧目录。
4. Phase 7 周期跳跃 1000/500，实际可行时 250；不在这些对照中重标材料。
5. Phase 8 相同物理尺寸 coarse/medium/fine：h=0.125/0.083333333333/0.0625 mm，分别 10×6×6 / 15×9×9 / 20×12×12 粒子；每档独立 coupon。保留原物理固定基底厚度，注意中档粒子中心恰在固定层边界的舍入/归属，应明确统一规则并披露边界离散误差。
6. 新失效算例 dt 与 dt/2 敏感性。若实际曝光内没有失效，不虚构失效窗口或首脱落事件。
7. Phase 9 optional fragment_drag（当前 runner 有延迟导入占位，`transport.py` 尚未实现；默认 relaxation 可运行）。记录 Re、阻力机制和保守力分配；低分辨碎屑分开处理。
8. Phase 10 ResolvedStreamingFieldProvider：VTU/VTP/legacy VTK/CSV，显式四种单位元数据、插值/梯度/覆盖/法向诊断及独立测试。**尚未实现**。
9. 参数敏感性扫参基础设施，完整回归和旧强制碎裂可复现核对。
10. Phase 11 生成无箭头新动画、规定的收敛/标定/局域性/能量/碎片尺寸 PNG+PDF、四阶段摘要、ZH/EN 报告及复现命令；原目录不改。

## 科学解释与必须保留的局限

- 新局部场由 `A_y=-U L exp(-r²/(2L²))` 的旋度构造，叠加现有管流。局部压力振幅演示值为 0。它不是实际包膜微泡 CFD，没有双向耦合。
- 新载荷从 μ(grad u+grad uᵀ)-pI 推导，可能明显弱于旧独立指定 100 Pa 牵引。不能偷偷提高 U 或降低 D_break 来产生动画碎裂。
- 能量驱动是 `0.5*E_young*(循环键伸长幅)^2*键积分体积`，是可解释的方向应变能代理，并非 NOSB 全自由能的精确逐键分解。键积分体积从初始族按权重分配，求和等于总材料体积。Wcrit 是每档 coupon 标定密度乘该体积，不跨档复用。
- coupon 的 Gc 是指定循环张开弱平面协议下的实测耗散/单个投影裂面面积；通用加载路径和主流场下的网格无关性仍须验证。不能凭 coupon 通过宣布整个模型已收敛。
- 主算例能量账本记录实际代表周期的外功/动能/储能；不会凭空乘 DeltaN 作为完整 25,000 周期能量。循环跳跃疲劳仍未实验标定。
- 损伤耗散取同一坐标更新前后的实际 NOSB 储能下降；逐键账本分配按非负损伤增量×能量驱动加权，属于估计，不能称精确逐键热力学释放率。
- 若没有发生失效/脱落，完整报告未达到相关验收项；若粉碎指标高或网格差异大，如实标 INADEQUATE / NOT CONVERGED。

## 环境与命令

工作目录：`/home/lzy/projects/clot_clearance`。科学 Python：`/home/lzy/projects/computation_examples/_shared/python_environment/bin/python`。

运行统一设置 `PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=$PWD/build_regularization/numba_cache PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`，避免在冻结源码旁生成缓存。pytest 使用 `-p no:cacheprovider`；临时目录放在 `build_regularization`。

已完成粗档标定的复现命令（须换新输出目录）：

```bash
PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR="$PWD/build_regularization/numba_cache" PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B -m pd_clot.regularization.coupon --config configs/streaming_regularized_demo_coarse.json --spacing .000125 --Gc-demo .01 --output verification/regularization/calibration_coarse_replay
```

本任务没有授权子代理，未启动子代理；当前工具指令也禁止自行委派。不要恢复旧 BraVa 微泡批次或其他已停止项目。
