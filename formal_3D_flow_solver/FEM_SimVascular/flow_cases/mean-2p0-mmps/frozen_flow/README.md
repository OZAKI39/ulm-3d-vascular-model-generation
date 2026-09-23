# 新背景流场：入口平均速度 2.0 mm/s

主文件 `steady_flow_mean_2p0_mmps.vtu` 是最终原生求解场的无损副本。请使用同目录 manifest.json 的 SHA-256 固定版本。

`flow_arrays_si.npz`：points_m 为原坐标，tetra 为原四面体连接，velocity_m_s 为三分量节点速度，pressure_pa 为节点压力；boundary_triangles 和 facet_tags 保留边界几何与分类。所有数组使用 SI 单位。边界编号：WALL=1、OUTLET_03=2、OUTLET_01=3、INLET=4、OUTLET_02=5。

VTU 保留原生单元的局部顶点顺序；NPZ 保留原正向四面体顺序，两者节点坐标与每个单元的节点集合已核对一致。stFile_*.bin 含原生积分历史，用于同一求解器与并行布局的重启；不能以仅有 VTU 替代完整重启状态。

本目录为独立的新候选，旧 frozen_reference 未替换。详细验收与适用范围见上一级 FLOW_2MMPS_REVIEW.md。
