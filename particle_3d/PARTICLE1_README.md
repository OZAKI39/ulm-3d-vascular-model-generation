# Particle-1：单个球形微泡

本阶段在只读 Particle-0 / frozen FEM 上推进一个微泡。没有质量或惯性积分，
平移速度 V=u，旋转速度 Ω=0.5 curl(u)。保留 Stokes 阻力函数
F=-6πμa(V-u)，但当前没有其它力，不做力积分。核心单位全部 SI。

## 接口

    from particle_3d.microbubble import MicrobubbleState, stokes_drag_force_n
    from particle_3d.integrator import advance_single_microbubble
    new_state = advance_single_microbubble(state, fem_field, dt_s)

调用者必须传入有限、正的 dt_s；接口没有默认时间步。位置严格使用旧位置速度：
x_new=x_old+dt_s*u(x_old)。返回前在新位置重新查场，只刷新 V 和 Ω，
不再次推进位置。这样输出的速度和角速度对应输出位置，仍为显式 Euler。
输入 state 不变，三维数组是不可写的 float64 独立副本。
起点或终点查不到场时抛出 FieldUnavailableError；该 API 不提供壁面响应。

## 一次验证的可重复规则

真实起点来自正式 INLET 的全部 179 个相邻 canonical tetra。
比较每个中心的速度与面积加权入口向内法线的点积，取最大者；精确并列取最小 tetra_id。
本次得到 tetra_id=208001。规则不读取后续轨迹结果。

首次运行前固定 dt0=0.25*所选 tetra 最短边长/起点速度，再取 dt0/2、dt0/4。
这三个值只出现在验证数据和程序中；没有写入正式模拟配置，也没有冻结 production dt。
三条轨迹是同一个微泡的三次独立验证重放，不是三微泡系统。

ValidationBoundaryClassifier 检查每条 Euler 线段与全部五个原始边界面的首次接触，
VTK 仅筛选空间候选，交点由直线平面公式和三角形权重确定。
几何误差预算按 float64 eps、坐标尺度与三角形条件数固定。
边角同时接触时优先 WALL，再 INLET，再出口，避免把壁面接触误报为正常出口。
WALL / INLET 事件记为 FAIL；OUTLET 事件记录真实交点、分段时间和原始 trial endpoint 后停止。
没有反弹、修正轨迹、推回或调 dt 重试。边界交点只是终止事件，不是新的壁面模型。
不穿过 WALL 只验证中心折线，不能证明有限半径球与壁面的间隙为正。

SonoVue adapter 只校验并调用外部冻结 sampler，固定 seed=20260920、N=1，
把 diameter_um 转成 diameter_m，再取一半为 radius_m。
所有外部 manifest 文件保持只读，CSV 附带绑定 SHA 的 metadata。
formal_simulation_population=false。本次 NumPy 2.5.3、Python 3.13.11；
原 sampler 合同的跨版本逐位一致性未被承诺。本次只核对同一环境中原函数的精确输出。
原 2D 代码仅作组织方式参考，没有复制二维坐标、µm 核心或边界行为。

## WSL CPU 复现

在仓库根目录执行；此环境已有 Particle-0 使用的依赖：

    P1_PY=/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python
    "$P1_PY" -B particle_3d/scripts/run_particle1_real_validation.py
    "$P1_PY" -B particle_3d/scripts/generate_particle1_report.py
    "$P1_PY" -B -m pytest -q -p no:cacheprovider particle_3d/tests/particle0 particle_3d/tests/particle1

第二条命令可以加 --stage 0 至 --stage 7 单独重画一步。重跑会更新本阶段生成文件，
不会修改 FEM / Particle-0 / SonoVue 原代码。每个步骤的测试、PNG、CSV/JSON 和中文说明永久保留。
真实轨迹测试重新查询每一个保存的位置，并复查每一条保存线段；不只检查汇总 JSON。

源码本地提交并确认工作区源码干净后，生成最终审核记录：

    "$P1_PY" -B particle_3d/scripts/finalize_particle1.py --handoff-root /home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular

该脚本在当前 Particle-1 worktree 校验 frozen payload，重复原有 P0 46 项测试及 P1 全部测试。
原有 18 项 handoff 测试锁定冻结分支名称，因此在原干净冻结 checkout 执行。
脚本只校验，不调用 CFD。JSON 中 git_commit 指向被测试的实现提交，
后续可单独提交报告和日志作为证据；逐文件 SHA 绑定实际测试代码、图和数据。

## 范围与审核

只有一个球，没有 RBC、其它微泡、壁面力、接触、润滑、黏附、重力、浮力、
Brownian、声力、lift、added mass、Basset history、LAMMPS、连续注入或浓度模型。
Particle-0 数学、容差、原有测试、批量实现和全部冻结输入不变。
本轮不做 P0 增强采样测试，不优化 sampler，不做 FEM 时间步研究。
人工图审及“是否有明显一步跳跃”均保留 PENDING_USER_REVIEW。
完成 Particle-1 后停止，Particle-2 未开始。
