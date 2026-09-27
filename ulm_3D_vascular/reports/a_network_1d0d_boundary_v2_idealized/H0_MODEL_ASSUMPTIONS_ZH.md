# H0 模型假设

本轮采用理想化模型：2410 为主动指定的水力源点，其余 122 个结构叶节点为 p_ref=0 Pa 的参考压力终端；0 Pa 是表压基准。这些定义均不是生理标注，结果不代表真实小鼠脑循环 ground truth。

1. structural root node **2410** 被本模型正式指定为 `IDEALIZED_HYDRAULIC_SOURCE`。机器摘要使用 `source_definition=IDEALIZED_STRUCTURAL_ROOT`。无需再以缺少生理入口标签为由停止。
2. 唯一水力域是已经有 provenance 的 analysis A：**7419 个原节点、7418 条原边、1 个连通分量、0 个 SWC 图环**。精确插入 3 个虚拟切点后求解图为 7422 节点、7421 边。完整发布文件其余 42 个分量仅作溯源背景。
3. H0 = `ALL_STRUCTURAL_LEAVES_REFERENCE_PRESSURE`。除 2410 外全部 degree=1 结构叶节点共 **122 个**，均定义为 `IDEALIZED_REFERENCE_PRESSURE_TERMINAL`，p_ref=0 Pa。统一设零是可重复的参考边界假设，不是测量结论。所有末端同压是一项额外的理想化假设；共同平移压力基准并不能把本来不同的末端压力变成相同。
4. O3 原节点 **4484** 本身就是 H0 终端；p_O3(real)=0 Pa 是模型定义。`geometry-derived R3 = NOT_AVAILABLE`；不存在的 A 下游几何不生成有限 R3。
5. 模型为 geometry-informed、steady、rigid、Newtonian、Poiseuille。复用 v1 的 SI 单位、沿边线性半径和稳定解析阻力积分。μ=0.00345312 Pa·s，ρ=1056.0 kg/m³；1D 稳态阻力不依赖密度，3D 保持该密度。
6. 不加入 Fahraeus–Lindqvist、hematocrit、phase separation、非牛顿黏度或血管顺应性。没有为让分流更平均而调节半径或终端压力。
7. 用途是考察 ROI 截断和等出口压力对分流的影响。不能称为 true physiological mouse cerebral flow。

# 工作点与流向

先取源点 1 Pa、所有终端 0 Pa。单位源压力下 ROI 入口 Q=2.968457224259e-18 m³/s，方向与 FEM 入口一致，因此 λ=5226.146254566908。所有压力、流量乘此 λ，使 ROI 入口达到 1.551359160440e-14 m³/s。源点表压 5226.146254567 Pa 仅是该工作点的模型压力尺度。

全 A 源入流 765.226012087 pL/s 与 ROI 入口 15.513591604 pL/s 分开输出。ROI 流向由真实切点几何、保留的 ROI 侧相邻边和 outward normal 定义；SWC parent-child 只作代数存储方向。入口正号表示入 ROI，出口正号表示出 ROI，绝不对入口流量取绝对值以强行匹配。

原始来源、精确映射和 40 项输入 SHA 沿用 v1 并在运行前重验；[v1 provenance](../a_network_1d0d_boundary_v1/A_NETWORK_PROVENANCE_ZH.md) 保持只读。求解数据见 [graph](data/analysis_A_H0_graph_si.npz)、[方向定义](data/roi_port_flow_sign_convention.json)、[unit solve](data/unit_solve_summary.json)。
