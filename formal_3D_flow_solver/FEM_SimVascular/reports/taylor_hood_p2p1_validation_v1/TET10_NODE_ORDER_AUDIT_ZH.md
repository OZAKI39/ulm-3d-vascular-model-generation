# TET10 / TRI6 原生局部节点顺序

| local index（0-based） | 几何角色 | corner pair |
|---|---|---|
| 0 | corner | 0 |
| 1 | corner | 1 |
| 2 | corner | 2 |
| 3 | corner | 3 |
| 4 | edge midpoint | 0–1 |
| 5 | edge midpoint | 1–2 |
| 6 | edge midpoint | 2–0 |
| 7 | edge midpoint | 0–3 |
| 8 | edge midpoint | 1–3 |
| 9 | edge midpoint | 2–3 |

TRI6 的前三点为 corners，后续依次为 (0,1)、(1,2)、(2,0) 中点。

依据固定本地源 `nn.cpp:175–176` 的 solver-to-basis permutation，及 `FE/Basis/NodeOrderingConventions.cpp:283–315` 的 FE lattice：TET10 映射为 `{1,2,3,0,5,9,8,4,6,7}`；TRI6 为 `{1,2,0,4,5,3}`。另以 `vtk_xml_parser.cpp:111` 的四个 quadratic faces 交叉核对。

原求解器参考顶点采用 origin-last；`read_msh.cpp:665–703` 会在标准 determinant 为正时交换 corners 0/1，TET10 同时交换 5/6 和 7/8。本轮在生成新 mesh 时明确执行同一角点排序，再按上述 edge pair 生成 midpoints。corner 全局 ID、位置、单元节点集合和拓扑保持。初次准备的方向审计误用了未转换顺序，已保留该不完整目录及失败日志，重新使用时间戳 case；没有把这一约定问题当成物理网格翻转。

永久测试验证 nodal Kronecker 性质、准确 edge map、共享面中点和 15 个原生积分点正 Jacobian。
