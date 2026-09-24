# 固定压力 BC：BLOCKED_DO_NOT_USE

结构根 2410 已精确找到，但现有记录只把它作为假定入口；没有找到将它标注为水力源点的数据集证据。外围端点的水力身份也未确认。按本轮源点 STOP 条件，未执行 A 的 operating-point 求解，也未创建或运行新 3D case。

**本轮没有可用于 FEM 的 O1/O2/O3 压力建议值。** JSON 中 raw/shifted pressure、pressure differences、source solution SHA 都是 null；mapping SHA 已保存。旧 preprocessing 的压力不能冒充本轮解，它使用未经本轮确认的 structural-root/all-leaves 假设和不同 operating point。

若后续源点和终端条件明确，先设源点单位压差，计算有符号真实 ROI inlet cut 流量 `Q_unit`，确认 `Q_unit>0`，再用 `λ=1.551359160440232e-14/Q_unit` 缩放整网压力和流量。目标是 ROI inlet cut，不是 A source total flow；不允许取 abs 掩盖反向流。

全网所得 p 位于真实 cut。人工延长段应按 `p_cap=p_real−R_ext Q_outward` 转到 FEM cap；入口的 Q_outward<0，符号自然反转。这是冻结 operating point 的估计，3D 求解后实际分流变动会影响修正精度。应先按当前网格复核 extension profile，再作 gauge shift。

可选择所有 cap 压力减去同一个最小值；每一对压力差必须保持不变。0 Pa 仅为参考压力。gauge invariance、extension 符号和精确变半径阻力均有永久测试。

| Port | 到当前 FEM cap 轴向长 / μm | 当前 cap 等效半径 / μm | 20 截面全长度阻力估计 / Pa·s·m⁻³ |
|---|---:|---:|---:|
| INLET | 15.708205 | 1.571326 | 2.242282455e+16 |
| O1 | 11.797444 | 1.177505 | 5.281376695e+16 |
| O2 | 11.696299 | 1.161572 | 5.470624510e+16 |
| O3 | 12.692440 | 1.267740 | 4.249596780e+16 |


表内数值仅为阻力估计，不是压力建议值。source/terminal gate 未通过，因此没有复制 mean-2p0-mmps，也没有 config diff、新稳态结果或三组分流比较图。
