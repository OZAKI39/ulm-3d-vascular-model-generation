# Gmsh R1 V2 reference 验证

**GMSH_REFERENCE_MESH_V2_PASS**。唯一V2已通过冻结mesh门槛，永久冻结。后续GPU svMultiPhysics/PETSc production路线已达到steady并正常保存第71步；局部场gate结论仍为 **REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC**，误差分解为 **MIXED_FIELD_ERROR**。旧CPU的1800s timeout仅为历史阶段，不能再作为当前field状态。见[流场报告](MANUFACTURED_TRAJECTORY_FIELD_VALIDATION_ZH.md)与[分解诊断](V2_P1_REPRESENTATION_DIAGNOSIS_ZH.md)。mesh PASS不代表field PASS。

## 唯一端面修正

复用 V1 原始 BREP，R=596.779229982443 µm，L=7161.35075978932 µm，三个25 µm壳层、272周向与30轴向间隔不变。cap target固定30 µm；入口及出口共8条核心径向曲线由3段确定性二分为6段，每段86.963204997 µm。没有扫描其它target。

4个贯穿全长的内部核心径向面原本是transfinite面，保持它们会把6个径向区间铺满整根管道。因此仅这4个面解除transfinite约束，使用既有Algorithm 6及200 µm背景尺寸，让端部细分过渡到中央core。wall、shell圆柱面及壳层径向面保持原约束；CAD未分裂、未移动，仍为16个共形fluid体。这是端面修正所需的局部拓扑过渡，不新增全长径向shell。

| 端面数据 | 旧trial 3（引用） | 本轮唯一2D trial | 最终3D V2双端面 |
|---|---:|---:|---:|
| cap target | 30 µm | 30 µm | 30 µm |
| core径向intervals | 3 | 6 | 6 |
| 单端面节点 / 三角形 | 1918 / 3562 | 2119 / 3964 | 2119 / 3964 |
| 速度trace相对L2 | 0.943571% | 0.210666% | 0.210666% |
| integrated Q相对误差 | 0.563419% | 0.165549% | 0.165549% |
| rim速度 | 0 | 0 | 0 |
| 中心速度 | 0.2 m/s | 0.2 m/s | 0.2 m/s |
| 联合门槛 | FAIL | PASS | PASS |

本次2D耗时0.114282s，无体单元。两端在浮点误差内一致，解析rim坐标残差约2.0e−16 m/s；未归一化流量、移动节点或修改profile。保留本次实际两端网格[NPZ](gmsh_reference/CAP_TOPOLOGY_CORRECTION_MESH.npz)与[修正记录](CAP_TOPOLOGY_CORRECTION.json)。旧3次校准未重跑、记录未修改。

对Poiseuille二次函数，P1插值亏损由三角形边长及面积共同决定。仅缩小面内target无法拆分固定径向长边；本次在相同30 µm target下将它们二分，Q误差从0.563419%降至0.165549%，直接验证了该约束的影响。并未逐项分离拓扑重排与新增节点的独立贡献。

（历史图件已退休；结论保留于本报告/JSON，tracked文件可从清理前算例Git提交91ebc021恢复。）

## 唯一V2与QA

bulk背景目标固定200 µm。Gmsh4.15.2、HXT10、一阶tetra、7线程、7.68-core cgroup；原有优化与容差不变。实际生成281.814s，监督总耗时282.713s，低于600s目标及1200s硬限。进程CPU时间1011.961s，平均使用3.585核；请求7线程不等于持续满核。峰值进程树RSS0.2618GiB，无OOM。

| 检查 | 实测 | 结论 |
|---|---:|---|
| 节点 / tetra | 36211 / 174785 | 全部一阶四面体 |
| 零体积 / 负体积 / 重复tetra | 0 / 0 / 0 | PASS |
| connected fluid | 1 | PASS |
| exterior free / non-manifold edge | 0 / 0 | PASS |
| WALL / INLET / OUTLET | 完整、方向正确 | PASS |
| 内部CAD三角面 | 50546 | 每面2个fluid邻居，PASS |
| 重复节点 / 重复CAD面 | 0 / 0 | PASS |
| 中央有限facet最大径向偏差 | 0.039805292 µm | ≤0.05，PASS |
| 中央wall facet集合 | 3264个，与V1相同 | 仅比较时容许1e−12 m舍入，未修改mesh |
| 实际wall altitude中位数 | 24.998332495 µm | [20,30]，PASS |
| 中央bulk中位数 | 215.420260362 µm | P0的1.173873291倍，PASS |
| final cap trace | 速度0.210666%、Q 0.165549% | ≤1% / ≤0.5%，PASS |
| strict frozen probes | 960 / 960 | 无移动、删点、外推，PASS |

wall altitude按中央z/L=0.4–0.6真实WALL facet → owning tetra → 垂直高度计算，min/P5/median/P95/max（µm）：24.998332495 / 24.998332495 / 24.998332495 / 25.303404266 / 25.303404266。

bulk按同一中央区域、tet质心d_wall≥150 µm，先每tet六边算术均值，再取中位数，共2104个tet。P0冻结参考183.512362009 µm；本轮比值1.173873处于[0.75,1.25]。全长统计仅诊断，不替代该口径。有限facet几何采用三角形上的真实径向极值，不以节点半径或质心弦高代替。任何输出均无snap/weld/project。

（历史图件已退休；结论保留于本报告/JSON，tracked文件可从清理前算例Git提交91ebc021恢复。）

## 冻结身份与后续

BREP与原生MSH使用mm；CFD的VTU、RAW_MESH/QA坐标使用m，仅作既定单位换算。MSH不能未经单位换算直接替代CFD输入。

- [V2 plan](gmsh_reference/GMSH_REFERENCE_MESH_V2_PLAN.json)：`35e7de9f20ef1e966b84a497133f04175725f50b5fd3ecef6df6f6cd44e1ef55`。
- 原始BREP：`7d27d1af3cefa9c5aa3b6bc442f25fdc070d8db12193ae1684c7c7aaaa9450f7`，与V1逐字节一致。
- MSH：`b67a8e1219d07e0ccdeccda4de18b6e65bddab2241cfc826f6142823a123241c`。
- VTU：`3a16731d70547ab4b657000c1c99d9af6aa4f0581f205cd045c37e410c1bbbd7`。
- Gmsh library：`9db3090d3b720c57b76bcbfa01d13854823ae2698c91343c20bdd4c2b81f6317`。
- [完整QA](GMSH_REFERENCE_MESH_V2_RESULT.json) · [运行记录](monitoring/gmsh_reference_v2_cap_topology/EXECUTION.json)。

7项有限facet/冻结reference接口测试通过（1.91s）；实际3D验收由上述测量给出。V1、旧P0/P1及历史3次cap calibration未被改写。

Q1：端面修正通过。Q2：约25 µm真实wall高度与P0-like bulk同时通过。Q3：唯一CFD运行预算失败，局部u/G/s/curl精度未评估。Q4/Q5：按失败停止规则未执行轨迹与dt验证。当前已停止；不生成V3，不进入BraVa。
