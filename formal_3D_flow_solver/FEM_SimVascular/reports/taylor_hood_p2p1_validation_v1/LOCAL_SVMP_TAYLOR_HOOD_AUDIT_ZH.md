# 当前实际 svMultiPhysics 的 Taylor-Hood 支持

本地审阅源码与 baseline 二进制构建源的 13 个关键文件 SHA 一致。当前固定 commit 为 `c3f0bb892b765b718f61069ecd9726dbc6d177fd`；运行源另有原先 Stage Q 的 GPU 改动，本轮全部只读。二进制 SHA 为 `0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7`。

`Add_equation/Use_taylor_hood_type_basis=true` 对 TET10 设 `nFs=2`，速度空间继承 10 节点二次基函数，压力空间切换为 TET4 角点一阶基函数；`fluid.cpp` 仅在 `nFs==1` 时启用当前等阶 VMS 分支。面、边界条件和原生输出均有相应实现。边中点压力槽位用单位约束保留在矩阵存储中，不代表新增 P2 压力自由度。

源码支持不等于实际算例成功；正式结论由运行和独立数值审核决定。完整源码定位和哈希见 [能力记录](data/local_svmp_taylor_hood_capability.json) 与 [逐行摘录](reference/native_source_evidence.txt)。

原生 VTU 会把坐标写成 float32。所有速度分量仍逐节点写为 double。后续若产生解，必须用原始输入的精确 midpoint 坐标和节点映射审核真正 P2 速度，不能对变成弯曲的舍入坐标直接求导，也不能降成 P1。
