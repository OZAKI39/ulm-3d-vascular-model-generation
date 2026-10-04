# PFPD_AT1_SPATIAL_DISCRETIZATION_AUDIT_V1

**INCOMPLETE — NEED_GEOMETRY_DECISION。** 已完成几何预审、同一粗点云的隔离积分比较、旧临界区间复现和基线测试。根据用户提示词第3节hard gate，尚未生成细网格。该停止发生在物理问题定义阶段，不是数值求解失败。

## 1. 为什么做这一轮
检查原约1%应变失稳是否在固定物理问题、仅细化粒子离散后仍然存在。必须先明确“固定的物理问题”是什么。

## 2. h和delta
h是点间距；delta是非局部物理视界。旧h=0.125 mm，delta=0.25125 mm。后续主研究必须固定delta，不能按h缩放。当前尚无跨网格结果。

## 3. 是否保持相同物理问题
尚未确定。旧代码的体素材料域边长0.375 mm；点包围盒和夹持中心距离0.25 mm。直接固定其中一个并复用全部旧规则，会改变另一个。

## 4. Geometry / volume audit
详见[完整预审](preaudit/GEOMETRY_AND_QUADRATURE_AUDIT.md)。旧体积总和0.052734375 mm³，等于0.375³；点包围盒体积0.015625 mm³，二者相差3.375倍。顶面面积也存在0.140625与0.0625 mm²两种解释。代码约束的是z=±0.125 mm处的点，没有声明细化时将其延拓为精确边界面还是完整刚性体素层。

同27点的守恒包围盒积分规则已另做比较，**未采用**。固定仿射位形下，两种规则的能量密度均为约0.1010025 Pa，但总能量和反力改变。因此它不能被隐藏为普通细化。

## 5. 三档mesh信息
|档位|点数|h mm|delta/h|键数|平均family|
|---|---:|---:|---:|---:|---:|
|旧Level 0|27|0.125|2.01|185|13.7037037|
|Level 1|未生成|未确定|未确定|N/A|N/A|
|Level 2|未生成|未确定|未确定|N/A|N/A|
旧27点全部shape rank=3，全部horizon受边界截断，没有完整球形内部family。

## 6. Intact mechanics
隔离粗积分比较的Fbar/Ftilde仿射误差均小于2e-15；独立SVK能量密度检查通过。该检查没有替代每档完整mechanics gate。Julia/upstream核对尚未执行，因为几何hard gate早于细化及外部reference阶段。

## 7. Stability interval
旧Level 0在新目录按本轮11点加载表重新求驻点与曲率，复现区间[0.0099921875, 0.01]，宽7.8125e-06。185 strong→185 weak仍成立，位形、相场、反力、能量、残差和负曲率方向均通过与旧状态的比较。见[复现记录](level_0/COARSE_REFERENCE_REPRODUCTION.json)。

## 8. 临界区间是否收敛
INCONCLUSIVE。Level 1/2未生成，不能给出空间一致性PASS。旧解析阈值为0.01，尚无跨网格threshold检查。

## 9. Reaction / stress是否收敛
INCONCLUSIVE。旧1%应变反力约4.2822421875 μN；按两种候选面积归一化分别是30.4515 Pa和68.515875 Pa。材料域尚未选定，因此不把其中一个命名为本轮权威nominal stress。提议包围盒守恒积分同一仿射位形反力约1.2688125 μN，对0.0625 mm²为20.301 Pa；这不是新的分支稳定性结果。

## 10. Active-set switch
旧Level 0仍为100%键共同切换，其他档位N/A。

## 11. Negative-cone direction
旧参考负曲率方向已复现，原缩放下上端点Rayleigh约−0.01470214。该缩放与DOF数有关，没有冒充跨网格统一物理尺度。精确比较见复现JSON。

## 12. Mode局域化
跨网格PR、top-10%幅值占比、公共网格相关性尚未计算，不将缺失值写成0。

## 13. Boundary association
尚未进行跨网格判定；旧算例全部family截断不能单独证明boundary artifact。

## 14. 3×3×3分类
**INCONCLUSIVE**。旧失稳仍可复现，但其空间稳健性不能由一个分辨率或积分规则替换来证明。

## 15. Tests与保护
本轮实际重跑旧测试：255 passed / 3 known failed / 0 skipped；新增预审测试：6 passed。旧失败ID、消息及traceback与上一轮逐项比较：一致。新增测试不修改旧失败。
保护核验：19341个旧文件、305个链接；状态PASSED。生产修改0，历史修改0。

## 16. 当前状态
完整标签见[STATUS.json](STATUS.json)。主研究INCOMPLETE、物理域一致性NOT_VERIFIED、区间／反力一致性INCONCLUSIVE；模型仍REMAINS_NOT_VERIFIED。完成两个预审PNG/PDF图对：[几何定义](figures/geometry_domain_ambiguity.png)、[同点云积分比较](figures/coarse_quadrature_comparison.png)。提示词要求的九组跨分辨率图待几何决定和实际计算后生成，未用虚构网格或结果占位。

## 17. 唯一下一步
明确固定的材料域和夹持定义。建议选择0.25 mm边界节点立方体、精确上下夹持面、守恒边界体积；旧粗结果独立保留，修正coarse另名并重新通过全部门槛。另一选择是0.375 mm体素域与原0.125 mm厚夹持区，但需据此另行设计点云。等待用户对几何含义的回答，不能把沉默当作选择。
