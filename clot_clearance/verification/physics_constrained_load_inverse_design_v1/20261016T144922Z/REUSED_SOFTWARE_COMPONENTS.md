# 源码复用表

| 复用资产 | 本轮作用 | 修改范围 |
|---|---|---|
| `analysis/v3_softening_pilot_v1/native.py` Native / Contact | 原生读取、SI协调器、输出字段、接触独立核对 | 原文件零修改；AST加载原类到隔离输出目录 |
| 原冻结 `build/liblammps.so` | endpoint-star变分软化、不可逆历史、原生拓扑 | 零编译、零修改，运行时检查SHA |
| `analysis/softening_feasibility_v1/model.py` | 接受前VV坐标下损伤预测；间隔/事件原生力能核对 | 零修改 |
| `analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py` | 存活本体图及锚点可达性 | 零修改 |
| `benchmark_corrected_v3/shear_x/provenance` | 1792材料点、336锚点及原始预裂/完整参考星 | 三个文件逐字节复用；只改变本轮外部力边界 |
| 已验收CFD表面和稀疏映射 | 3611面片、1927原始材料点的空间支持 | 原文件零修改；新增标量仿射权重与转置功率接口 |
| MEVO F10S执行树 `14731fd...` | APECSS、压力、中心轨迹、Bjerknes数据与当前范围 | 只读；未执行MEVO新扫描 |
| 历史MEVO F1 k-Wave基底 | 真实复基底与相位能力核查 | 只读；不冒充当前F10S几何 |

新增模块见本轮analysis目录README。没有确认可删除且影响本轮的废弃实现，因此没有删除历史代码或做大范围重构。
