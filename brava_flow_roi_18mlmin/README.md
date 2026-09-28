# BraVa BG001 RMCA 四端口流场与微泡（18 mL/min）

本目录是独立新算例，不修改原打印件、小鼠流场或微泡核心代码。

- 几何：用户确认使用 BALANCED 已验收四端口血管芯，包含人工接管，排除五面盒体。
- candidate 0、15 是同一几何的刚性姿态。无重力、刚壁、牛顿流体模型下复用同一物理解，旋转位置、矢量及轨迹以展示两种姿态。
- 入口：18,000 μL/min = 18 mL/min = 3.0e-7 m³/s。
- 材料：沿用现有流程的 rho=1056 kg/m³、mu=0.00345312 Pa·s；这是血液模型参数，并非已测得的打印灌注液参数。
- 分流目标：各出口 1.0e-7 m³/s，约 33.333%。先测量等流量校准的出口压力，再用全压力出口正式求解验证。校准的强制分流不能当作自然分流验证。
- CFD：稳态目标的隐式瞬态推进，dt=0.01 s；微泡积分另设 dt=0.0005 s。二者不同。
- 微泡：约 1500 条，保留原 SonoVue 条件粒径分布，真实入口流量加权采样，不按出口结果挑选、不强制出口。
- 原始单位 mm，计算单位 m；平移原点和两种打印姿态矩阵见 inputs/geometry.json。

## 文件

- `inputs/`：源几何副本、端口身份、中心线、单位/变换、初步管网设计。
- `mesh/generation/`：官方 SimVascular/TetGen 原始生成文件。
- `mesh/SV_MESH/`：审核后的实际求解体网格、面片身份和数组。
- `cases/`：各次校准/正式求解的独立配置、原始结果和日志；失败记录保留。
- `vendor/`：现有几何审核、求解监视和 WSS 计算代码的独立副本。
- `scripts/`：新几何适配与运行脚本。
- `visualization/`：沿用 rotate_visualization 风格的新图及动画。
- `microbubble/`：用户已要求停止。正式批次已停止，8条未完成轨迹与已完成的单条基准/一致性验证记录保留；没有完成1500条生产轨迹，没有生成本批微泡动画。状态见 `microbubble/data/USER_STOP.json`。
- `reports/`、`logs/`：数值检查与运行日志。

服务器：`vast4090:/workspace/brava_flow_roi_18mlmin_20260928/`。
独立 supervisor：该目录下 `supervisord.conf`；不控制其他项目服务。

当前网格为单档，不宣称网格无关。原制造表面采用多边形环及补偿半径，结果不是未经工程修改的原生脑血管生理预测。

流场和两姿态的各5段流场动画已完成：[candidate 0](visualization/candidate_0/OPEN_RESULTS.html) · [candidate 15](visualization/candidate_15/OPEN_RESULTS.html)。微泡停止不影响这些结果。内部截面检查发现约百分之几的局部速度通量缺陷，已通过散度积分核对，见 `reports/internal_sections/`；边界总流量守恒不代表局部输运精度已充分验证。
